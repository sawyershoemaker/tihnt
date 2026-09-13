(() => {
  const MAX_DIMENSION = 512;
  const MAX_CELLS = 16384;
  const CAPTURE_INTERVAL = 16;
  const VISIBILITY_CHECKS = { checkOpacity: true, checkVisibilityCSS: true };
  const state = {
    entries: [], elements: [], byElement: new Map(), ancestors: [], cells: [], w: 0, h: 0,
    cacheDirty: true, scanAll: true, dirty: new Set(),
    last: null, inFlight: false, pending: false, generation: 0, lastFull: 0,
    mines: -1, settingsReady: false, captureMode: false, buffer: ''
  };
  let debug = false;
  let timer = null;
  let lastCapture = -Infinity;
  let clipQuery = null;
  let clipMutation = null;
  let clipPlan = null;
  let clipParent;
  let layoutVersion = 0;
  function log(...args) { if (debug) console.log('[tihnt/content]', ...args); }
  function schedule(fullScan = false, layoutChanged = false) {
    if (layoutChanged) ++layoutVersion;
    state.scanAll ||= fullScan;
    if (!state.settingsReady || document.hidden) return;
    if (state.inFlight) { state.pending = true; return; }
    if (timer !== null) return;
    timer = setTimeout(() => {
      timer = null;
      lastCapture = performance.now();
      tick();
    }, Math.max(0, CAPTURE_INTERVAL - (performance.now() - lastCapture)));
  }
  function resync() {
    ++state.generation;
    clipQuery?.();
    state.last = null;
    state.cacheDirty = true;
    schedule(true);
  }
  function rebuild() {
    clipParent = undefined;
    const entries = [];
    let w = 0, h = 0;
    for (const element of document.querySelectorAll('div[id^="cell_"].cell')) {
      const match = /^cell_(\d+)_(\d+)$/.exec(element.id);
      if (!match) continue;
      const x = Number(match[1]), y = Number(match[2]);
      if (!Number.isSafeInteger(x) || !Number.isSafeInteger(y) || x >= MAX_DIMENSION || y >= MAX_DIMENSION) return false;
      entries.push({ element, x, y });
      w = Math.max(w, x + 1); h = Math.max(h, y + 1);
    }
    if (!w || !h || w * h > MAX_CELLS || entries.length !== w * h) return false;
    const occupied = new Set();
    for (const entry of entries) {
      entry.index = entry.y * w + entry.x;
      if (occupied.has(entry.index)) return false;
      occupied.add(entry.index);
    }
    state.entries = entries;
    state.elements = entries.map(({ element }) => element);
    state.byElement = new Map(entries.map((entry) => [entry.element, entry]));
    const ancestors = new Set();
    for (const { element } of entries) for (let parent = element.parentElement; parent && !ancestors.has(parent); parent = parent.parentElement)
      ancestors.add(parent);
    state.ancestors = [...ancestors];
    state.w = w; state.h = h;
    state.cells = new Array(w * h).fill(0);
    state.cacheDirty = false;
    state.scanAll = true;
    return true;
  }
  function decode(element) {
    const classes = String(element.className || '').split(/\s+/);
    let opened = false, flagged = false, number = null;
    for (const name of classes) {
      if (/^[a-z]+_opened$/.test(name)) opened = true;
      // closed_flag is a background for flag mode, not an actual placed flag.
      if (/^[a-z]+_flag$/.test(name)) flagged = true;
      const match = /^[a-z]+_type(\d+)$/.exec(name);
      if (match) number = Number(match[1]);
    }
    if (opened) {
      if (number === 10 || number === 11) return -1; // Game over.
      // Never invent a zero clue when an opened cell's class is incomplete.
      return number !== null && number >= 0 && number <= 8 ? 10 + number : -2;
    }
    return flagged ? 2 : 0;
  }
  function clipsContent(style) {
    const x = style.overflowX, y = style.overflowY, clip = style.clip, content = style.contentVisibility;
    return (x && x !== 'visible') || (y && y !== 'visible') || (clip && clip !== 'auto') ||
      /\b(paint|strict|content)\b/.test(style.contain) || (content && content !== 'visible');
  }
  function intersectionTolerance(a, b) {
    const ulp = (value) => 2 ** (Math.floor(Math.log2(Math.max(1, Math.abs(value)))) - 23);
    return Math.min(.001, Math.max(1e-7, 2 * Math.max(ulp(a), ulp(b))));
  }
  function sameEdge(a, b) {
    return a === b || Math.abs(a - b) <= intersectionTolerance(a, b);
  }
  function sameRect(a, b, rounded = false) {
    if (rounded) return sameEdge(a.left, b.left) && sameEdge(a.top, b.top) && sameEdge(a.right, b.right) && sameEdge(a.bottom, b.bottom);
    return Math.abs(a.left - b.left) <= 1e-7 && Math.abs(a.top - b.top) <= 1e-7 &&
      Math.abs(a.width - b.width) <= 1e-7 && Math.abs(a.height - b.height) <= 1e-7;
  }
  function unsupportedShape(style) {
    const path = style.clipPath, mask = style.maskImage, radius = style.borderRadius;
    return (path && path !== 'none') || (mask && mask !== 'none') ||
      (radius && radius !== '0px' && radius.split(/[\s/]+/).some((value) => value && parseFloat(value) !== 0) && clipsContent(style));
  }
  function geometry(verifyClip) {
    clipPlan = null;
    const sample = [
      document.getElementById('cell_0_0'),
      document.getElementById('cell_' + (state.w - 1) + '_0'),
      document.getElementById('cell_0_' + (state.h - 1)),
      document.getElementById('cell_' + (state.w - 1) + '_' + (state.h - 1))
    ];
    if (sample.some((element) => !element?.checkVisibility(VISIBILITY_CHECKS))) return null;
    const rects = sample.map((element) => element?.getBoundingClientRect());
    if (rects.some((r) => !r || r.width <= 0 || r.height <= 0)) return null;
    const ancestors = new Map();
    let clipped = false;
    for (const element of new Set([...sample, ...state.ancestors])) {
      const style = window.getComputedStyle(element);
      ancestors.set(element, style);
      const clips = clipsContent(style);
      clipped ||= clips;
      // Intersection rectangles cannot represent rounded or arbitrary masks.
      if (unsupportedShape(style) || style.display === 'none' || style.opacity === '0' || style.contentVisibility === 'hidden') return null;
      if (style.transform !== 'none') {
        const matrix = new DOMMatrixReadOnly(style.transform);
        if (!matrix.is2D || matrix.a <= 0 || matrix.d <= 0 || Math.abs(matrix.b) > 1e-7 || Math.abs(matrix.c) > 1e-7) return null;
      }
      if (style.perspective !== 'none' || style.offsetPath !== 'none' ||
          (style.rotate !== 'none' && Math.abs(parseFloat(style.rotate.split(/\s+/).at(-1)) % 360) > 1e-7) ||
          (style.scale !== 'none' && style.scale.split(/\s+/).some((value) => !(Number(value) > 0)))) return null;
      for (const animation of element.getAnimations()) {
        if ((!animation.pending && (animation.playState === 'finished' || animation.playState === 'idle')) ||
            animation.effect?.target !== element || animation.effect.pseudoElement) continue;
        const properties = animation.effect.getKeyframes().flatMap((frame) => Object.keys(frame))
          .filter((key) => !['offset', 'computedOffset', 'easing', 'composite'].includes(key));
        // Motion can start from an identity frame without another DOM mutation.
        if (!properties.length || properties.some((key) => key.startsWith('--') ||
            /^(transform|transformOrigin|translate|rotate|scale|perspective|perspectiveOrigin|offset.*)$/.test(key))) return null;
      }
    }
    const left = Math.min(...rects.map((r) => r.left)), top = Math.min(...rects.map((r) => r.top));
    const right = Math.max(...rects.map((r) => r.right)), bottom = Math.max(...rects.map((r) => r.bottom));
    const cellWidth = (right - left) / state.w, cellHeight = (bottom - top) / state.h;
    if (rects.some((r, i) => Math.abs(r.left - (left + (i % 2 ? state.w - 1 : 0) * cellWidth)) > .5 ||
        Math.abs(r.top - (top + (i > 1 ? state.h - 1 : 0) * cellHeight)) > .5 ||
        Math.abs(r.width - cellWidth) > .5 || Math.abs(r.height - cellHeight) > .5)) return null;
    const vv = window.visualViewport;
    const result = {
      rect_l: left, rect_t: top, rect_w: right - left, rect_h: bottom - top,
      cell_px: Math.round(Math.max(rects[0].width, rects[0].height)),
      ox: 0, oy: 0, vv_x: vv?.offsetLeft || 0, vv_y: vv?.offsetTop || 0,
      vv_scale: vv?.scale || 1, dpr: window.devicePixelRatio || 1
    };
    if (!Object.values(result).every(Number.isFinite) || result.rect_w > 16384 || result.rect_h > 16384 ||
        Math.abs(left) > 1000000 || Math.abs(top) > 1000000 || Math.abs(result.vv_x) > 1000000 || Math.abs(result.vv_y) > 1000000 ||
        result.vv_scale <= 0 || result.vv_scale > 16 || result.dpr <= 0 || result.dpr > 16) return null;
    if (clipped) {
      let cellRects;
      let exactPartition = false;
      if (verifyClip || clipParent === undefined) {
        clipParent = undefined;
        const parent = state.entries[0].element.parentElement;
        let inFlow = true, direct = true;
        cellRects = new Array(state.entries.length);
        for (const { element } of state.entries) {
          const style = ancestors.get(element) || window.getComputedStyle(element);
          if (!element.checkVisibility(VISIBILITY_CHECKS) || unsupportedShape(style)) return null;
          inFlow &&= style.position === 'static' || style.position === 'relative';
          direct &&= element.parentElement === parent;
        }
        for (const { element, index } of state.entries) cellRects[index] = element.getBoundingClientRect();
        const bounds = { left, top, width: right - left, height: bottom - top };
        clipParent = parent && inFlow && direct && !clipsContent(ancestors.get(parent)) && sameRect(parent.getBoundingClientRect(), bounds) ? parent : null;
        exactPartition = !clipParent;
        for (const { index, x, y } of state.entries) {
          const r = cellRects[index], previous = x ? cellRects[index - 1] : null;
          const above = y ? cellRects[index - state.w] : null;
          if (r.width <= 0 || r.height <= 0 || Math.abs(r.left - (left + x * cellWidth)) > .5 ||
              Math.abs(r.top - (top + y * cellHeight)) > .5 || Math.abs(r.width - cellWidth) > .5 || Math.abs(r.height - cellHeight) > .5 ||
              (previous ? !sameEdge(previous.right, r.left) || !sameEdge(previous.top, r.top) || !sameEdge(previous.bottom, r.bottom) : !sameEdge(r.left, left)) ||
              (above ? !sameEdge(above.bottom, r.top) : !sameEdge(r.top, top)) ||
              (x === state.w - 1 && !sameEdge(r.right, right)) || (y === state.h - 1 && !sameEdge(r.bottom, bottom))) return null;
          exactPartition &&= (previous ? previous.right === r.left && previous.top === r.top && previous.bottom === r.bottom : r.left === left) &&
            (above ? above.bottom === r.top : r.top === top) &&
            (x !== state.w - 1 || r.right === right) && (y !== state.h - 1 || r.bottom === bottom);
        }
      }
      clipPlan = {
        targets: clipParent ? [clipParent] : state.elements,
        cellRects,
        exactPartition: !!cellRects && exactPartition,
        // A ratio observer misses clip movement at an unchanged visible area.
        signature: JSON.stringify([result, window.innerWidth, window.innerHeight, ...[...ancestors].map(([element, style]) => {
          const r = element.getBoundingClientRect();
          return [r.left, r.top, r.width, r.height, element.scrollLeft, element.scrollTop,
            element.clientLeft, element.clientTop, element.clientWidth, element.clientHeight,
            style.borderLeftWidth, style.borderTopWidth, style.borderRightWidth, style.borderBottomWidth,
            style.paddingLeft, style.paddingTop, style.paddingRight, style.paddingBottom,
            style.overflowX, style.overflowY, style.overflowClipMargin, style.clip, style.contain, style.contentVisibility,
            style.position, style.transform, style.translate, style.zoom];
        })])
      };
    }
    return result;
  }
  function observeClip(targets) {
    return new Promise((resolve) => {
      let observer, timeout, finished = false;
      const finish = (entries) => {
        finished = true;
        observer?.disconnect(); clearTimeout(timeout);
        if (clipQuery === cancel) clipQuery = null;
        resolve(entries);
      };
      const cancel = () => finish(null);
      clipQuery = cancel;
      clipMutation = () => {
        // Queued records were acquired before this mutation was delivered.
        if (finished || observer?.takeRecords().length) {
          ++layoutVersion;
          if (!finished) cancel();
        }
      };
      try {
        observer = new IntersectionObserver((entries) => finish(entries));
        for (const target of targets) observer.observe(target);
        timeout = setTimeout(cancel, 250);
      } catch { cancel(); }
    });
  }
  function tiledClip(rects, snap, cellRects, left, top, right, bottom) {
    const first = new Int32Array(snap.h).fill(snap.w), last = new Int32Array(snap.h).fill(-1), counts = new Int32Array(snap.h);
    let firstRow = snap.h, lastRow = -1;
    for (const rect of rects) {
      const cell = cellRects[rect.index];
      if (rect.l !== Math.max(cell.left, left) || rect.t !== Math.max(cell.top, top) ||
          rect.r !== Math.min(cell.right, right) || rect.b !== Math.min(cell.bottom, bottom)) return false;
      const x = rect.index % snap.w, y = Math.floor(rect.index / snap.w);
      first[y] = Math.min(first[y], x); last[y] = Math.max(last[y], x); ++counts[y];
      firstRow = Math.min(firstRow, y); lastRow = Math.max(lastRow, y);
    }
    // Unique, consecutive cells in every row must span the box without a missing cell or row.
    for (let y = firstRow; y <= lastRow; ++y) {
      if (!counts[y] || counts[y] !== last[y] - first[y] + 1 ||
          cellRects[y * snap.w + first[y]].left > left || cellRects[y * snap.w + last[y]].right < right) return false;
    }
    return cellRects[firstRow * snap.w].top <= top && cellRects[lastRow * snap.w].bottom >= bottom;
  }
  function rectangularClip(entries, snap, cellRects) {
    const empty = { clip_l: snap.rect_l, clip_t: snap.rect_t, clip_w: 0, clip_h: 0 };
    const rects = [];
    let l = Infinity, t = Infinity, r = -Infinity, b = -Infinity;
    const boardRight = snap.rect_l + snap.rect_w, boardBottom = snap.rect_t + snap.rect_h;
    for (const entry of entries) {
      const rect = entry.intersectionRect;
      const left = Math.max(snap.rect_l, rect.left), top = Math.max(snap.rect_t, rect.top);
      const right = Math.min(boardRight, rect.right), bottom = Math.min(boardBottom, rect.bottom);
      if (entry.isIntersecting && right > left && bottom > top) {
        rects.push({ l: left, t: top, r: right, b: bottom, index: cellRects ? state.byElement.get(entry.target).index : -1 });
        l = Math.min(l, left); t = Math.min(t, top); r = Math.max(r, right); b = Math.max(b, bottom);
      }
    }
    if (!rects.length) return empty;
    if (rects.length > 1 && !(cellRects && tiledClip(rects, snap, cellRects, l, t, r, b))) {
      function coordinates(values) {
        const result = new Map();
        let anchor = -Infinity;
        for (const value of [...new Set(values)].sort((a, b) => a - b)) {
          const tolerance = intersectionTolerance(anchor, value);
          if (value - anchor > tolerance) anchor = value;
          result.set(value, anchor);
        }
        return result;
      }
      // Browser intersection edges lose float32 precision independently at fractional scales.
      const xs = coordinates(rects.flatMap((rect) => [rect.l, rect.r]));
      const snappedY = coordinates(rects.flatMap((rect) => [rect.t, rect.b]));
      const snapped = rects.map((rect) => ({ l: xs.get(rect.l), r: xs.get(rect.r), t: snappedY.get(rect.t), b: snappedY.get(rect.b) }))
        .filter((rect) => rect.r > rect.l && rect.b > rect.t);
      if (!snapped.length) return empty;
      if (Math.min(...snapped.map((rect) => rect.l)) !== xs.get(l) || Math.max(...snapped.map((rect) => rect.r)) !== xs.get(r) ||
          Math.min(...snapped.map((rect) => rect.t)) !== snappedY.get(t) || Math.max(...snapped.map((rect) => rect.b)) !== snappedY.get(b)) return empty;
      // Sweep the union to reject holes instead of filling them with its bounding box.
      const ys = [...new Set(snapped.flatMap((rect) => [rect.t, rect.b]))].sort((a, b) => a - b);
      const indices = new Map(ys.map((y, i) => [y, i]));
      const events = snapped.flatMap((rect) => [[rect.l, 1, indices.get(rect.t), indices.get(rect.b)], [rect.r, -1, indices.get(rect.t), indices.get(rect.b)]])
        .sort((a, b) => a[0] - b[0]);
      const counts = new Int32Array(ys.length * 4), full = new Uint8Array(ys.length * 4);
      function update(node, lo, hi, start, end, change) {
        if (start <= lo && hi <= end) counts[node] += change;
        else {
          const mid = (lo + hi) >> 1;
          if (start < mid) update(node * 2, lo, mid, start, end, change);
          if (end > mid) update(node * 2 + 1, mid, hi, start, end, change);
        }
        full[node] = counts[node] > 0 || (hi - lo > 1 && full[node * 2] && full[node * 2 + 1]) ? 1 : 0;
      }
      let x = events[0][0];
      for (const [next, change, start, end] of events) {
        if (next > x && !full[1]) return empty;
        update(1, 0, ys.length - 1, start, end, change);
        x = next;
      }
    }
    return { clip_l: l, clip_t: t, clip_w: r - l, clip_h: b - t };
  }
  function snapshot(verifyClip = true) {
    if (state.scanAll && state.entries.some(({ element }) => !element.isConnected || element.ownerDocument !== document)) state.cacheDirty = true;
    if (state.cacheDirty || !state.entries.length) {
      if (!rebuild()) {
        state.cacheDirty = true;
        state.entries = []; state.elements = []; state.byElement.clear(); state.ancestors = []; state.cells = []; state.dirty.clear();
        state.w = state.h = 0;
        return { w: 0, h: 0, cells: [], mines_total: state.mines };
      }
    }
    const entries = state.scanAll ? state.entries : state.dirty;
    for (const entry of entries) state.cells[entry.index] = decode(entry.element);
    state.dirty.clear(); state.scanAll = false;
    const rect = state.cells.some((cell) => cell < 0) ? null : geometry(verifyClip);
    if (!rect) { clipParent = undefined; return { w: 0, h: 0, cells: [], mines_total: state.mines }; }
    return { w: state.w, h: state.h, cells: state.cells.slice(), mines_total: state.mines, ...rect };
  }
  function send(message) {
    return new Promise((resolve) => {
      try {
        chrome.runtime.sendMessage(message, (response) => {
          const error = chrome.runtime.lastError;
          if (error) log('send failed', error.message);
          resolve(!error && response?.ok === true);
        });
      } catch (error) { log('extension unavailable', String(error)); resolve(false); }
    });
  }
  function sameGeometry(a, b) {
    return ['rect_l', 'rect_t', 'rect_w', 'rect_h', 'clip_l', 'clip_t', 'clip_w', 'clip_h', 'vv_x', 'vv_y', 'vv_scale', 'dpr', 'mines_total']
      .every((key) => a[key] === b[key]);
  }
  async function tick() {
    if (!state.settingsReady || document.hidden) return;
    if (state.inFlight) { state.pending = true; return; }
    const generation = state.generation;
    state.inFlight = true;
    try {
      // Reuse only the query target; publication always verifies the current cells and geometry.
      let snap = snapshot(false);
      if (snap.w && clipPlan) {
        const plan = clipPlan;
        const layout = layoutVersion;
        const entries = await observeClip(plan.targets);
        if (document.hidden || generation !== state.generation) return;
        // Decode the latest cells after the rendering step, while retaining its geometry.
        snap = snapshot();
        if (snap.w && clipPlan) {
          function validEntries() {
            if (entries.length !== plan.targets.length) return false;
            const remaining = new Set(plan.targets);
            for (const entry of entries) {
              if (!remaining.delete(entry.target)) return false;
              const cell = state.byElement.get(entry.target);
              const rect = cell ? clipPlan.cellRects[cell.index] : entry.target.getBoundingClientRect();
              if (!sameRect(entry.boundingClientRect, rect, true)) return false;
            }
            return true;
          }
          const stale = layout !== layoutVersion || plan.signature !== clipPlan.signature || plan.targets.length !== clipPlan.targets.length ||
            (plan.targets !== clipPlan.targets && plan.targets.some((target, i) => target !== clipPlan.targets[i])) || (entries && !validEntries());
          Object.assign(snap, entries && !stale ? rectangularClip(entries, snap,
            clipPlan.targets === state.elements && clipPlan.exactPartition ? clipPlan.cellRects : null) :
            { clip_l: snap.rect_l, clip_t: snap.rect_t, clip_w: 0, clip_h: 0 });
          // Hide stale hints during layout churn. Repeated identical failures wait for the fallback.
          if (stale && (!state.last || !sameGeometry(snap, state.last))) state.pending = true;
        }
      }
      clipMutation = null;
      const last = state.last;
      let full = !last || snap.w !== last.w || snap.h !== last.h || Date.now() - state.lastFull > 60000;
      let message;
      if (full) message = { type: 'full', ...snap };
      else {
        const updates = [];
        for (let i = 0; i < snap.cells.length; ++i) {
          if (snap.cells[i] !== last.cells[i]) updates.push({ x: i % snap.w, y: Math.floor(i / snap.w), s: snap.cells[i] });
        }
        if (!updates.length && sameGeometry(snap, last)) return;
        // A restart or flood opening is much smaller as a packed full board.
        full = snap.cells.length >= 64 && updates.length * 8 > snap.cells.length;
        if (full) message = { type: 'full', ...snap };
        else {
          const { cells, w, h, ...metadata } = snap;
          message = { type: 'delta', updates, ...metadata };
        }
      }
      const ok = await send(message);
      if (generation === state.generation) {
        state.last = ok ? snap : null;
        if (ok && full) state.lastFull = Date.now();
      }
    } catch (error) { log('capture failed', String(error)); }
    finally {
      clipMutation = null;
      state.inFlight = false;
      if (state.pending || generation !== state.generation) {
        state.pending = false;
        schedule();
      }
    }
  }

  chrome.storage.local.get({ debug: false, mines_total: -1 }).then((settings) => {
    debug = !!settings.debug;
    const value = Number(settings.mines_total);
    state.mines = Number.isInteger(value) && value >= -1 && value <= MAX_CELLS ? value : -1;
  }).catch(() => {}).finally(() => { state.settingsReady = true; schedule(true); });
  chrome.storage.onChanged.addListener((changes, area) => {
    if (area !== 'local') return;
    if (changes.debug) debug = !!changes.debug.newValue;
    if (changes.mines_total) {
      const value = Number(changes.mines_total.newValue);
      state.mines = Number.isInteger(value) && value >= -1 && value <= MAX_CELLS ? value : -1;
      schedule();
    }
  });

  window.addEventListener('keydown', (event) => {
    if (event.isComposing) return;
    const editable = event.target?.closest?.('input, textarea, select, [contenteditable]:not([contenteditable="false"])');
    const shiftM = event.shiftKey && !event.ctrlKey && !event.altKey && !event.metaKey && event.key.toLowerCase() === 'm';
    if (!state.captureMode && (editable || !shiftM)) return;
    if (shiftM && event.repeat) { event.preventDefault(); return; }
    function commit() {
      const value = state.buffer === '' ? -1 : Number(state.buffer);
      state.mines = Number.isInteger(value) && value >= -1 && value <= MAX_CELLS ? value : -1;
      chrome.storage.local.set({ mines_total: state.mines }).catch(() => {});
      state.captureMode = false; state.buffer = ''; schedule();
    }
    if (shiftM) {
      if (state.captureMode) commit();
      else { state.captureMode = true; state.buffer = ''; }
    } else if (event.key === 'Enter') commit();
    else if (event.key === 'Escape') { state.captureMode = false; state.buffer = ''; }
    else if (/^\d$/.test(event.key)) { if (state.buffer.length < 5) state.buffer += event.key; }
    else if (event.key === 'Backspace') state.buffer = state.buffer.slice(0, -1);
    else return;
    event.preventDefault(); event.stopPropagation();
  }, true);

  function primaryCell(element) {
    return /^cell_\d+_\d+$/.test(element.id || '') && element.matches?.('div[id^="cell_"].cell');
  }
  function cellContent(target) {
    if (state.ancestors.includes(target)) return true;
    for (let node = target; node; node = node.parentNode) {
      if (state.byElement.has(node)) return true;
      if (state.ancestors.includes(node)) return false;
    }
    return false;
  }
  function containsCell(node) {
    if (node.nodeType !== 1) return false;
    if (state.byElement.has(node) || primaryCell(node)) return true;
    for (const element of node.querySelectorAll?.('div[id^="cell_"].cell') || [])
      if (primaryCell(element)) return true;
    return false;
  }
  const observer = new MutationObserver((mutations) => {
    let relevant = false, layoutChanged = false, appearanceChanged = false;
    for (const mutation of mutations) {
      if (mutation.type === 'attributes') {
        const entry = state.byElement.get(mutation.target);
        if (entry) {
          layoutChanged ||= mutation.attributeName !== 'class';
          appearanceChanged ||= mutation.attributeName === 'class';
          if (mutation.attributeName === 'id' || !primaryCell(mutation.target)) state.cacheDirty = true;
          if (mutation.attributeName === 'class') state.dirty.add(entry);
          relevant = true;
        } else if (primaryCell(mutation.target)) {
          state.cacheDirty = true; relevant = layoutChanged = true;
        } else if (state.ancestors.includes(mutation.target) || (state.entries[0] && mutation.target.contains?.(state.entries[0].element))) {
          relevant = layoutChanged = true; // An ancestor style/class can move or resize the board.
        } else if (cellContent(mutation.target)) {
          relevant = appearanceChanged = true;
        }
      } else if (mutation.type === 'childList') {
        const changed = [...mutation.addedNodes, ...mutation.removedNodes];
        if (changed.some(containsCell)) {
          state.cacheDirty = true; relevant = layoutChanged = true;
        }
        if (cellContent(mutation.target) && changed.some((node) => node.nodeType === 1 || node.nodeType === 3))
          relevant = appearanceChanged = true;
      } else if (mutation.type === 'characterData' && mutation.target.nodeType === 3 && cellContent(mutation.target.parentNode)) {
        relevant = appearanceChanged = true;
      }
    }
    if (appearanceChanged) clipMutation?.();
    if (relevant) schedule(false, layoutChanged);
  });
  observer.observe(document, { subtree: true, childList: true, characterData: true, attributes: true });
  chrome.runtime.onMessage.addListener((message, sender, respond) => {
    if (message?.type === 'force_full') { resync(); respond({ ok: true }); }
  });
  document.addEventListener('visibilitychange', () => {
    if (document.hidden) { clearTimeout(timer); timer = null; clipQuery?.(); }
    else resync();
  });
  window.addEventListener('pageshow', resync);
  function animationChanged(event) {
    if (event.pseudoElement || !state.entries.length) return;
    const target = event.target;
    if (state.byElement.has(target) || state.byElement.has(target?.closest?.('div[id^="cell_"].cell')) ||
        target?.contains?.(state.entries[0].element) || state.byElement.has(target?.querySelector?.('div[id^="cell_"].cell'))) schedule(false, true);
  }
  for (const name of ['animationstart', 'animationend', 'animationcancel', 'transitionrun', 'transitionstart', 'transitionend', 'transitioncancel'])
    document.addEventListener(name, animationChanged, true);
  document.addEventListener('scroll', () => schedule(false, true), { passive: true, capture: true });
  window.addEventListener('resize', () => schedule(false, true));
  window.visualViewport?.addEventListener('scroll', () => schedule(false, true), { passive: true });
  window.visualViewport?.addEventListener('resize', () => schedule(false, true));
  // Fallback catches missed class changes and layout shifts without DOM mutations.
  setInterval(() => schedule(true), 1000);
})();
