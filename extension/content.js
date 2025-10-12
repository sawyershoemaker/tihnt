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
    minesTotal: 0,
    captureMode: false,
    captureBuffer: ''
  };
  let lastNoBoardLog = 0;

  try {
    chrome.storage?.local?.get?.({ mines_total: 0 }).then((res) => {
      const v = Number(res?.mines_total || 0);
      STATE.minesTotal = Number.isFinite(v) && v >= 0 ? Math.floor(v) : 0;
      log('loaded mines_total', STATE.minesTotal);
    });
  } catch {}

  // shift + m to enter mines num
  window.addEventListener('keydown', (e) => {
    try {
      const isShiftM = e.shiftKey && (e.key === 'M' || e.key === 'm');
      const commit = () => {
        const n = STATE.captureBuffer.trim() === '' ? 0 : Math.max(0, Math.floor(Number(STATE.captureBuffer)) || 0);
        STATE.minesTotal = n;
        try { chrome.storage?.local?.set?.({ mines_total: n }); } catch {}
        log('capture mines mode: set', n, 'buffer=', STATE.captureBuffer);
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

  function capture() {
    // dom structure
    const cells = Array.from(document.querySelectorAll('div[id^="cell_"].cell'));
    if (cells.length === 0) return null;
    let maxX = -1, maxY = -1;
    let minLeft = Infinity, minTop = Infinity, maxRight = -Infinity, maxBottom = -Infinity;
    let cellPx = 20;
    const coordOf = (id) => {
      const m = /^cell_(\d+)_(\d+)$/.exec(id);
      if (!m) return null;
      // flip x and y
      return { x: parseInt(m[1], 10), y: parseInt(m[2], 10) };
    };
    // first pass: compute extents and origin from bounding rects
    for (const el of cells) {
      const id = el.id || '';
      const c = coordOf(id);
      if (!c) continue;
      if (c.x > maxX) maxX = c.x;
      if (c.y > maxY) maxY = c.y;
      const r = el.getBoundingClientRect();
      if (r.width && r.height) cellPx = Math.round(Math.max(cellPx, Math.max(r.width, r.height)));
      if (r.left < minLeft) minLeft = r.left;
      if (r.top < minTop) minTop = r.top;
      if (r.right > maxRight) maxRight = r.right;
      if (r.bottom > maxBottom) maxBottom = r.bottom;
    }
    const w = maxX + 1;
    const h = maxY + 1;
    if (w <= 0 || h <= 0) return null;
    // geometry using css and visual viewport
    const vv = window.visualViewport;
    const rect = { l: minLeft, t: minTop, w: Math.max(0, maxRight - minLeft), h: Math.max(0, maxBottom - minTop) };
    const vv_x = vv ? vv.offsetLeft : window.scrollX;
    const vv_y = vv ? vv.offsetTop : window.scrollY;
    const vv_scale = vv ? vv.scale : 1;
    const origin = { x: 0, y: 0 }; // legacy, unused when rect/vv present

    const states = new Array(w*h).fill(0);
    let bombDetected = false;
    const normType = (cls) => {
      const m = /^[a-z]+_type(\d+)$/.exec(cls);
      return m ? parseInt(m[1], 10) : null;
    };
    for (const el of cells) {
      const id = el.id || '';
      const c = coordOf(id);
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
        else { states[idx] = 10; } // opened zero as fallback
      } else if (closed) {
        // represent flagged cells as state=2 (game::CellState::Mine) so backend knows flags already exist
        states[idx] = hasFlag ? 2 : 0;
      }
    }
    const snap = { w, h, states, origin, cellPx, bombDetected, rect, vv: { x: vv_x, y: vv_y, scale: vv_scale } };
    log('capture', { w, h, cellPx, bombDetected, cellCount: cells.length });
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

  function tick() {
    try {
      const s = capture();
      if (s) {
        diffAndSend(s);
      } else if (debug) {
        const now = Date.now();
        if (now - lastNoBoardLog > 3000) {
          lastNoBoardLog = now;
          log('no board found on page');
        }
      }
    } catch (e) { log('tick error', String(e)); }
  }

  // respond to background "force_full" command
  try {
    chrome.runtime.onMessage.addListener((msg, sender, respond) => {
      try {
        if (msg && msg.type === 'force_full') {
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

  setInterval(tick, 200);
  const vv = window.visualViewport;
  if (vv) {
    const onGeom = () => {
      try {
        const s = capture();
        if (!s) return;
        const msg = { type: 'delta', updates: [], cell_px: s.cellPx, ox: s.origin.x, oy: s.origin.y, mines_total: STATE.minesTotal,
          rect_l: s.rect.l, rect_t: s.rect.t, rect_w: s.rect.w, rect_h: s.rect.h,
          vv_x: s.vv.x, vv_y: s.vv.y, vv_scale: s.vv.scale, dpr: window.devicePixelRatio || 1 };
        chrome.runtime.sendMessage(msg, () => {});
      } catch {}
    };
    vv.addEventListener('scroll', onGeom);
    vv.addEventListener('resize', onGeom);
    vv.addEventListener('zoom', onGeom);
    window.addEventListener('wheel', () => {
      try { requestAnimationFrame(onGeom); } catch {}
    }, { passive: true });
    document.addEventListener('scroll', () => { try { requestAnimationFrame(onGeom); } catch {} }, { passive: true, capture: true });
    window.addEventListener('scroll', () => { try { requestAnimationFrame(onGeom); } catch {} }, { passive: true });
    // fallback periodic geometry refresh in case events are throttled
    setInterval(onGeom, 250);
  } else {
    const onGeom = () => {
      try {
        const s = capture();
        if (!s) return;
        const msg = { type: 'delta', updates: [], cell_px: s.cellPx, ox: s.origin.x, oy: s.origin.y, mines_total: STATE.minesTotal,
          rect_l: s.rect.l, rect_t: s.rect.t, rect_w: s.rect.w, rect_h: s.rect.h,
          vv_x: s.vv.x, vv_y: s.vv.y, vv_scale: s.vv.scale, dpr: window.devicePixelRatio || 1 };
        chrome.runtime.sendMessage(msg, () => {});
      } catch {}
    };
    window.addEventListener('scroll', onGeom, { passive: true });
    window.addEventListener('resize', onGeom);
  }
})();


