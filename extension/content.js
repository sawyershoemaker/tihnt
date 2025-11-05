(() => {
  // script detects minesweeper.online boards from DOM

  let debug = false;
  function log(...args) {
    if (!debug) return;
    try { console.log('[mines-ext/content]', ...args); } catch {}
  }
  try {
    chrome.storage?.local?.get?.({ debug: false }).then((res) => { debug = !!res.debug; });
    chrome.storage?.onChanged?.addListener?.((changes, area) => {
      if (area === 'local' && Object.prototype.hasOwnProperty.call(changes, 'debug')) {
        debug = !!changes.debug.newValue;
        log('debug set to', debug);
      }
    });
  } catch {}

  const STATE = {
    last: null,
    w: 0,
    h: 0,
    origin: { x: 0, y: 0 },
    cellPx: 0,
    lastFullSent: 0,
    minesTotal: -1,
    captureMode: false,
    captureBuffer: '',
    rect: { l: 0, t: 0, w: 0, h: 0 },
    geometryDirty: true
  };
  let lastNoBoardLog = 0;

  const CAPTURE_MIN_MS = 50;
  let capturePending = false;
  let lastCaptureStamp = 0;
  let geomPending = false;

  function scheduleCapture(reason = 'scheduled') {
    try {
      if (capturePending) return;
      capturePending = true;
      const runner = () => {
        capturePending = false;
        const now = typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now();
        if (now - lastCaptureStamp < CAPTURE_MIN_MS) {
          const delay = Math.max(0, CAPTURE_MIN_MS - (now - lastCaptureStamp));
          setTimeout(() => scheduleCapture(reason), delay);
          return;
        }
        lastCaptureStamp = now;
        tick(reason);
      };
      if (typeof requestAnimationFrame === 'function') {
        requestAnimationFrame(runner);
      } else {
        setTimeout(runner, 0);
      }
    } catch {}
  }

  function scheduleGeometry(reason = 'geom') {
    try {
      if (geomPending) return;
      geomPending = true;
      const run = () => {
        geomPending = false;
        sendGeometryDelta(reason);
      };
      if (typeof requestAnimationFrame === 'function') {
        requestAnimationFrame(run);
      } else {
        setTimeout(run, 0);
      }
    } catch {}
  }

  let boardObserver = null;
  function ensureMutationObserver() {
    if (boardObserver || typeof MutationObserver !== 'function') return;
    try {
      boardObserver = new MutationObserver((mutations) => {
        let needCapture = false;
        let geometryTouched = false;
        for (const mutation of mutations) {
          if (mutation.type === 'attributes') {
            const target = mutation.target;
            const id = target && target.id;
            if (id && id.startsWith('cell_')) {
              needCapture = true;
            }
          } else if (mutation.type === 'childList') {
            const inspectList = (nodes) => {
              for (const node of nodes) {
                if (node && node.id && node.id.startsWith && node.id.startsWith('cell_')) {
                  geometryTouched = true;
                  needCapture = true;
                  return;
                }
                if (node && node.querySelector) {
                  const inner = node.querySelector('div[id^="cell_"]');
                  if (inner) {
                    geometryTouched = true;
                    needCapture = true;
                    return;
                  }
                }
              }
            };
            if (mutation.addedNodes && mutation.addedNodes.length) inspectList(mutation.addedNodes);
            if (mutation.removedNodes && mutation.removedNodes.length) inspectList(mutation.removedNodes);
          }
          if (needCapture && geometryTouched) break;
        }
        if (geometryTouched) STATE.geometryDirty = true;
        if (needCapture) {
          scheduleCapture(geometryTouched ? 'mutation:geom' : 'mutation');
        }
      });
      boardObserver.observe(document.body, { attributes: true, attributeFilter: ['class'], subtree: true, childList: true });
    } catch (e) {
      boardObserver = null;
      log('observer error', String(e));
    }
  }

  try {
    chrome.storage?.local?.get?.({ mines_total: -1 }).then((res) => {
      const raw = Number(res?.mines_total);
      let v = Number.isFinite(raw) ? Math.floor(raw) : -1;
      if (v < -1) v = -1;
      STATE.minesTotal = v;
      log('loaded mines_total', STATE.minesTotal);
    });
  } catch {}

  // shift + m to enter mines num
  window.addEventListener('keydown', (e) => {
    try {
      const isShiftM = e.shiftKey && (e.key === 'M' || e.key === 'm');
      const commit = () => {
        const trimmed = STATE.captureBuffer.trim();
        let next = -1;
        if (trimmed !== '') {
          const parsed = Math.floor(Number(trimmed));
          if (Number.isFinite(parsed) && parsed >= 0) {
            next = parsed;
          } else if (Number.isFinite(parsed) && parsed < 0) {
            next = -1;
          }
        }
        STATE.minesTotal = next;
        try { chrome.storage?.local?.set?.({ mines_total: next }); } catch {}
        log('capture mines mode: set', next, 'buffer=', STATE.captureBuffer);
        try {
          const msg = { type: 'delta', updates: [], cell_px: STATE.cellPx, ox: STATE.origin.x, oy: STATE.origin.y, mines_total: n };
          chrome.runtime.sendMessage(msg, () => {});
        } catch {}
        STATE.captureMode = false;
        STATE.captureBuffer = '';
      };

      if (isShiftM) {
        if (!STATE.captureMode) {
          STATE.captureMode = true;
          STATE.captureBuffer = '';
          log('capture mines mode: start');
        } else {
          commit();
          log('capture mines mode: commit via Shift+M');
        }
        e.preventDefault();
        e.stopPropagation();
        return;
      }
      if (!STATE.captureMode) return;
      if (e.key === 'Enter') {
        commit();
        e.preventDefault();
        e.stopPropagation();
        return;
      }
      if (e.key === 'Escape') {
        log('capture mines mode: cancel');
        STATE.captureMode = false;
        STATE.captureBuffer = '';
        e.preventDefault();
        e.stopPropagation();
        return;
      }
      if (e.key >= '0' && e.key <= '9') {
        STATE.captureBuffer += e.key;
        e.preventDefault();
        e.stopPropagation();
        return;
      }
      if (e.key === 'Backspace') {
        STATE.captureBuffer = STATE.captureBuffer.slice(0, -1);
        e.preventDefault();
        e.stopPropagation();
        return;
      }
      // ignore other keys
    } catch {}
  }, true);

  const coordRe = /^cell_(\d+)_(\d+)$/;
  const coordOf = (id) => {
    const m = coordRe.exec(id);
    if (!m) return null;
    return { x: parseInt(m[1], 10), y: parseInt(m[2], 10) };
  };

  function measureGeometryFromSamples(w, h) {
    let minLeft = Infinity, minTop = Infinity, maxRight = -Infinity, maxBottom = -Infinity;
    let cellPx = STATE.cellPx || 0;
    let seen = false;
    const samples = [
      [0, 0],
      [Math.max(0, w - 1), 0],
      [0, Math.max(0, h - 1)],
      [Math.max(0, w - 1), Math.max(0, h - 1)]
    ];
    for (const [sx, sy] of samples) {
      const el = document.getElementById(`cell_${sx}_${sy}`);
      if (!el) continue;
      const r = el.getBoundingClientRect();
      if (!r) continue;
      seen = true;
      if (r.width && r.height) {
        const size = Math.round(Math.max(r.width, r.height));
        if (size > cellPx) cellPx = size;
      }
      if (r.left < minLeft) minLeft = r.left;
      if (r.top < minTop) minTop = r.top;
      if (r.right > maxRight) maxRight = r.right;
      if (r.bottom > maxBottom) maxBottom = r.bottom;
    }
    if (!seen || !Number.isFinite(minLeft) || !Number.isFinite(minTop) || !Number.isFinite(maxRight) || !Number.isFinite(maxBottom)) {
      return null;
    }
    return {
      rect: { l: minLeft, t: minTop, w: Math.max(0, maxRight - minLeft), h: Math.max(0, maxBottom - minTop) },
      cellPx
    };
  }

  function measureGeometryFromCells(cells, defaults) {
    let minLeft = Infinity, minTop = Infinity, maxRight = -Infinity, maxBottom = -Infinity;
    let cellPx = defaults?.cellPx || 0;
    for (const el of cells) {
      const r = el.getBoundingClientRect();
      if (!r) continue;
      if (r.width && r.height) {
        const size = Math.round(Math.max(r.width, r.height));
        if (size > cellPx) cellPx = size;
      }
      if (r.left < minLeft) minLeft = r.left;
      if (r.top < minTop) minTop = r.top;
      if (r.right > maxRight) maxRight = r.right;
      if (r.bottom > maxBottom) maxBottom = r.bottom;
    }
    if (!Number.isFinite(minLeft) || !Number.isFinite(minTop) || !Number.isFinite(maxRight) || !Number.isFinite(maxBottom)) {
      return null;
    }
    return {
      rect: { l: minLeft, t: minTop, w: Math.max(0, maxRight - minLeft), h: Math.max(0, maxBottom - minTop) },
      cellPx
    };
  }

  function capture() {
    const cells = Array.from(document.querySelectorAll('div[id^="cell_"].cell'));
    if (cells.length === 0) return null;

    let maxX = -1;
    let maxY = -1;
    for (const el of cells) {
      const c = coordOf(el.id || '');
      if (!c) continue;
      if (c.x > maxX) maxX = c.x;
      if (c.y > maxY) maxY = c.y;
    }
    const w = maxX + 1;
    const h = maxY + 1;
    if (w <= 0 || h <= 0) return null;

    const states = new Array(w * h).fill(0);
    let bombDetected = false;
    const normType = (cls) => {
      const m = /^[a-z]+_type(\d+)$/.exec(cls);
      return m ? parseInt(m[1], 10) : null;
    };
    for (const el of cells) {
      const c = coordOf(el.id || '');
      if (!c) continue;
      const idx = c.y * w + c.x;
      const classes = (el.className || '').split(/\s+/);
      let opened = false, closed = false, tnum = null, hasFlag = false;
      for (const cl of classes) {
        if (/^[a-z]+_opened$/.test(cl)) opened = true;
        else if (/^[a-z]+_closed$/.test(cl)) closed = true;
        const v = normType(cl); if (v !== null) tnum = v;
        if (/_flag$/.test(cl)) hasFlag = true;
      }
      if (opened) {
        if (tnum === 10) { bombDetected = true; }
        else if (tnum !== null && tnum >= 0 && tnum <= 8) { states[idx] = 10 + tnum; }
        else { states[idx] = 10; }
      } else if (closed) {
        states[idx] = hasFlag ? 2 : 0;
      }
    }

    let rect = STATE.rect;
    let cellPx = STATE.cellPx || 20;
    if (STATE.geometryDirty || !rect || rect.w === 0 || rect.h === 0) {
      const sample = measureGeometryFromSamples(w, h);
      const geom = sample || measureGeometryFromCells(cells, { cellPx });
      if (geom) {
        rect = geom.rect;
        cellPx = geom.cellPx;
        STATE.rect = rect;
        STATE.cellPx = cellPx;
        STATE.geometryDirty = false;
      }
    }

    const vv = window.visualViewport;
    const vv_x = vv ? vv.offsetLeft : window.scrollX;
    const vv_y = vv ? vv.offsetTop : window.scrollY;
    const vv_scale = vv ? vv.scale : 1;
    const origin = { x: 0, y: 0 };

    const snap = { w, h, states, origin, cellPx, bombDetected, rect: rect || { l: 0, t: 0, w: 0, h: 0 }, vv: { x: vv_x, y: vv_y, scale: vv_scale } };
    if (debug) {
      log('capture', { w, h, cellPx, bombDetected, cellCount: cells.length });
    }
    return snap;
  }

  function diffAndSend(snap) {
    const now = Date.now();
    const needFull = !STATE.last || now - STATE.lastFullSent > 10000 || STATE.w !== snap.w || STATE.h !== snap.h || snap.bombDetected;
    if (needFull) {
      // if a bomb was detected, clear board to unknown to reset for next round.
      const cells = snap.bombDetected ? new Array(snap.w * snap.h).fill(0) : snap.states;
      const msg = { type: 'full', w: snap.w, h: snap.h, cells, cell_px: snap.cellPx, ox: snap.origin.x, oy: snap.origin.y, mines_total: STATE.minesTotal,
        rect_l: snap.rect.l, rect_t: snap.rect.t, rect_w: snap.rect.w, rect_h: snap.rect.h,
        vv_x: snap.vv.x, vv_y: snap.vv.y, vv_scale: snap.vv.scale, dpr: window.devicePixelRatio || 1 };
      log('send full', { w: snap.w, h: snap.h, cellsLen: cells.length, cellPx: snap.cellPx, origin: snap.origin });
      chrome.runtime.sendMessage(msg, (resp) => {
        if (resp?.ok !== true) {
          const err = (typeof chrome !== 'undefined' && chrome.runtime && chrome.runtime.lastError) ? chrome.runtime.lastError.message : undefined;
          log('bg send failed', { context: 'full', resp, err });
        }
      });
      STATE.last = cells;
      STATE.w = snap.w; STATE.h = snap.h; STATE.origin = snap.origin; STATE.cellPx = snap.cellPx; STATE.lastFullSent = now;
      return;
    }
    const updates = [];
    for (let y=0; y<snap.h; y++){
      for (let x=0; x<snap.w; x++){
        const i = y*snap.w + x; const a = snap.states[i]; const b = STATE.last[i] ?? 0;
        if (a !== b) updates.push({ x, y, s: a });
      }
    }
    if (updates.length > 0) {
      if (debug) {
        const describe = (s) => {
          if (s === 0) return 'unknown';
          if (s === 2) return 'mine';
          if (s >= 10 && s <= 18) return s === 10 ? 'blank' : String(s - 10);
          return String(s);
        };
        for (const u of updates) {
          log('cell update', { x: u.x, y: u.y, state: u.s, value: describe(u.s) });
        }
      }
      const msg = { type: 'delta', updates, cell_px: snap.cellPx, ox: snap.origin.x, oy: snap.origin.y, mines_total: STATE.minesTotal,
        rect_l: snap.rect.l, rect_t: snap.rect.t, rect_w: snap.rect.w, rect_h: snap.rect.h,
        vv_x: snap.vv.x, vv_y: snap.vv.y, vv_scale: snap.vv.scale, dpr: window.devicePixelRatio || 1 };
      log('send delta', { count: updates.length, cellPx: snap.cellPx, origin: snap.origin });
      chrome.runtime.sendMessage(msg, (resp) => {
        if (resp?.ok !== true) {
          const err = (typeof chrome !== 'undefined' && chrome.runtime && chrome.runtime.lastError) ? chrome.runtime.lastError.message : undefined;
          log('bg send failed', { context: 'delta', resp, err });
        }
      });
      STATE.last = snap.states;
    }
  }

  function tick(trigger = 'timer') {
    try {
      const s = capture();
      if (s) {
        diffAndSend(s);
      } else if (debug) {
        const now = Date.now();
        if (now - lastNoBoardLog > 3000) {
          lastNoBoardLog = now;
          log('no board found on page', trigger);
        }
      }
    } catch (e) { log('tick error', String(e)); }
  }

  function sendGeometryDelta(reason = 'geom') {
    try {
      if (STATE.w <= 0 || STATE.h <= 0) {
        scheduleCapture(`${reason}:warmup`);
        return;
      }
      const measurement = measureGeometryFromSamples(STATE.w, STATE.h);
      if (!measurement) {
        STATE.geometryDirty = true;
        scheduleCapture(`${reason}:fallback`);
        return;
      }
      const rect = measurement.rect;
      const measuredCellPx = measurement.cellPx || STATE.cellPx || 0;
      if (!rect || rect.w === 0 || rect.h === 0 || measuredCellPx <= 0) {
        STATE.geometryDirty = true;
        scheduleCapture(`${reason}:invalid`);
        return;
      }
      STATE.rect = rect;
      STATE.cellPx = measuredCellPx;
      STATE.geometryDirty = false;
      const vv = window.visualViewport;
      const msg = { type: 'delta', updates: [], cell_px: STATE.cellPx, ox: STATE.origin.x, oy: STATE.origin.y, mines_total: STATE.minesTotal,
        rect_l: rect.l, rect_t: rect.t, rect_w: rect.w, rect_h: rect.h,
        vv_x: vv ? vv.offsetLeft : window.scrollX, vv_y: vv ? vv.offsetTop : window.scrollY, vv_scale: vv ? vv.scale : 1, dpr: window.devicePixelRatio || 1 };
      if (debug) log('geom delta', { reason, rect, cellPx: STATE.cellPx });
      chrome.runtime.sendMessage(msg, () => {});
    } catch (e) {
      log('sendGeometryDelta error', String(e));
    }
  }

  // respond to background "force_full" command
  try {
    chrome.runtime.onMessage.addListener((msg, sender, respond) => {
      try {
        if (msg && msg.type === 'force_full') {
          STATE.geometryDirty = true;
          const s = capture();
          if (!s) { respond && respond({ ok: false }); return; }
          // Always send full regardless of timer/size change
          const cells = s.bombDetected ? new Array(s.w * s.h).fill(0) : s.states;
          const full = { type: 'full', w: s.w, h: s.h, cells, cell_px: s.cellPx, ox: s.origin.x, oy: s.origin.y, mines_total: STATE.minesTotal,
            rect_l: s.rect.l, rect_t: s.rect.t, rect_w: s.rect.w, rect_h: s.rect.h,
            vv_x: s.vv.x, vv_y: s.vv.y, vv_scale: s.vv.scale, dpr: window.devicePixelRatio || 1 };
          log('force_full send', { w: s.w, h: s.h, cellsLen: cells.length, cellPx: s.cellPx });
          chrome.runtime.sendMessage(full, () => {});
          STATE.last = cells;
          STATE.w = s.w; STATE.h = s.h; STATE.origin = s.origin; STATE.cellPx = s.cellPx; STATE.lastFullSent = Date.now();
          respond && respond({ ok: true });
          return;
        }
      } catch {}
    });
  } catch {}

  ensureMutationObserver();
  scheduleCapture('init');
  setInterval(() => scheduleCapture('fallback'), 1000);

  const vv = window.visualViewport;
  if (vv) {
    vv.addEventListener('scroll', () => { STATE.geometryDirty = true; scheduleGeometry('vv:scroll'); });
    vv.addEventListener('resize', () => { STATE.geometryDirty = true; scheduleGeometry('vv:resize'); });
    if (typeof vv.addEventListener === 'function') {
      vv.addEventListener('zoom', () => { STATE.geometryDirty = true; scheduleGeometry('vv:zoom'); });
    }
  }

  const handleScroll = (reason) => () => { STATE.geometryDirty = true; scheduleGeometry(reason); };
  window.addEventListener('wheel', () => { STATE.geometryDirty = true; scheduleGeometry('wheel'); }, { passive: true });
  document.addEventListener('scroll', handleScroll('doc-scroll'), { passive: true, capture: true });
  window.addEventListener('scroll', handleScroll('win-scroll'), { passive: true });
  window.addEventListener('resize', handleScroll('resize'));
})();


