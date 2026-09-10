(() => {
  const MAX_DIMENSION = 512;
  const MAX_CELLS = 16384;
  const CAPTURE_INTERVAL = 50;
  const state = {
    entries: [], byElement: new Map(), cells: [], w: 0, h: 0,
    cacheDirty: true, scanAll: true, dirty: new Set(),
    last: null, inFlight: false, pending: false, generation: 0, lastFull: 0,
    mines: -1, settingsReady: false, captureMode: false, buffer: ''
  };
  let debug = false;
  let timer = null;
  let lastCapture = -Infinity;
  function log(...args) { if (debug) console.log('[tihnt/content]', ...args); }
  function schedule(fullScan = false) {
    state.scanAll ||= fullScan;
    if (timer !== null) return;
    timer = setTimeout(() => {
      timer = null;
      lastCapture = performance.now();
      tick();
    }, Math.max(0, CAPTURE_INTERVAL - (performance.now() - lastCapture)));
  }
  function resync() {
    ++state.generation;
    state.last = null;
    state.cacheDirty = true;
    schedule(true);
  }
  function rebuild() {
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
    state.byElement = new Map(entries.map((entry) => [entry.element, entry]));
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
  function geometry() {
    const sample = [
      document.getElementById('cell_0_0'),
      document.getElementById('cell_' + (state.w - 1) + '_0'),
      document.getElementById('cell_0_' + (state.h - 1)),
      document.getElementById('cell_' + (state.w - 1) + '_' + (state.h - 1))
    ];
    const rects = sample.map((element) => element?.getBoundingClientRect());
    if (rects.some((r) => !r || r.width <= 0 || r.height <= 0)) return null;
    const left = Math.min(...rects.map((r) => r.left)), top = Math.min(...rects.map((r) => r.top));
    const right = Math.max(...rects.map((r) => r.right)), bottom = Math.max(...rects.map((r) => r.bottom));
    const vv = window.visualViewport;
    return {
      rect_l: left, rect_t: top, rect_w: right - left, rect_h: bottom - top,
      cell_px: Math.round(Math.max(rects[0].width, rects[0].height)),
      ox: 0, oy: 0, vv_x: vv?.offsetLeft || 0, vv_y: vv?.offsetTop || 0,
      vv_scale: vv?.scale || 1, dpr: window.devicePixelRatio || 1
    };
  }
  function snapshot() {
    if (state.cacheDirty || !state.entries.length) {
      if (!rebuild()) {
        state.cacheDirty = true;
        return { w: 0, h: 0, cells: [], mines_total: state.mines };
      }
    }
    const entries = state.scanAll ? state.entries : state.dirty;
    for (const entry of entries) state.cells[entry.index] = decode(entry.element);
    state.dirty.clear(); state.scanAll = false;
    const rect = geometry();
    if (!rect || state.cells.some((cell) => cell < 0)) return { w: 0, h: 0, cells: [], mines_total: state.mines };
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
    return ['rect_l', 'rect_t', 'rect_w', 'rect_h', 'vv_x', 'vv_y', 'vv_scale', 'dpr', 'mines_total']
      .every((key) => a[key] === b[key]);
  }
  async function tick() {
    if (!state.settingsReady || document.hidden) return;
    if (state.inFlight) { state.pending = true; return; }
    let snap;
    try { snap = snapshot(); } catch (error) { log('capture failed', String(error)); return; }
    const last = state.last;
    const full = !last || snap.w !== last.w || snap.h !== last.h || Date.now() - state.lastFull > 60000;
    let message;
    if (full) message = { type: 'full', ...snap };
    else {
      const updates = [];
      for (let i = 0; i < snap.cells.length; ++i) {
        if (snap.cells[i] !== last.cells[i]) updates.push({ x: i % snap.w, y: Math.floor(i / snap.w), s: snap.cells[i] });
      }
      if (!updates.length && sameGeometry(snap, last)) return;
      const { cells, w, h, ...metadata } = snap;
      message = { type: 'delta', updates, ...metadata };
    }
    const generation = state.generation;
    state.inFlight = true;
    const ok = await send(message);
    state.inFlight = false;
    if (generation === state.generation) {
      state.last = ok ? snap : null;
      if (ok && full) state.lastFull = Date.now();
    }
    if (state.pending || generation !== state.generation) {
      state.pending = false;
      schedule();
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

  const observer = new MutationObserver((mutations) => {
    let relevant = false;
    for (const mutation of mutations) {
      if (mutation.type === 'attributes') {
        const entry = state.byElement.get(mutation.target);
        if (entry) {
          if (mutation.attributeName === 'id') state.cacheDirty = true;
          state.dirty.add(entry); relevant = true;
        } else if (state.entries[0] && mutation.target.contains?.(state.entries[0].element)) {
          relevant = true; // An ancestor style/class can move or resize the board.
        }
      } else {
        const containsCell = (node) => node.nodeType === 1 &&
          (/^cell_\d+_\d+$/.test(node.id || '') || node.querySelector?.('div[id^="cell_"].cell'));
        if ([...mutation.addedNodes, ...mutation.removedNodes].some(containsCell)) {
          state.cacheDirty = true; relevant = true;
        }
      }
    }
    if (relevant) schedule();
  });
  observer.observe(document.body, { subtree: true, childList: true, attributes: true, attributeFilter: ['class', 'id', 'style'] });
  chrome.runtime.onMessage.addListener((message, sender, respond) => {
    if (message?.type === 'force_full') { resync(); respond({ ok: true }); }
  });
  document.addEventListener('visibilitychange', () => { if (!document.hidden) resync(); });
  window.addEventListener('pageshow', resync);
  document.addEventListener('scroll', () => schedule(), { passive: true, capture: true });
  window.addEventListener('resize', () => schedule());
  window.visualViewport?.addEventListener('scroll', () => schedule(), { passive: true });
  window.visualViewport?.addEventListener('resize', () => schedule());
  // Fallback catches missed class changes and layout shifts without DOM mutations.
  setInterval(() => schedule(true), 1000);
})();
