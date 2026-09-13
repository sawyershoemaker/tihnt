const { test } = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const path = require('node:path');
const source = (name) => fs.readFileSync(path.join(__dirname, '..', 'extension', name), 'utf8');
const event = () => ({ listeners: [], addListener(fn) { this.listeners.push(fn); }, fire(...args) { for (const fn of this.listeners) fn(...args); } });

function timers() {
  let time = 100, id = 0;
  const pending = new Map(), intervals = new Map();
  return {
    api: {
      performance: { now: () => time }, Date,
      setTimeout: (fn, delay = 0) => { pending.set(++id, { fn, at: time + delay }); return id; },
      clearTimeout: (key) => pending.delete(key),
      setInterval: (fn, delay) => { intervals.set(++id, { fn, delay }); return id; },
      clearInterval: (key) => intervals.delete(key)
    },
    async flush() {
      for (let i = 0; i < 30; ++i) {
        await Promise.resolve();
        await Promise.resolve();
        if (pending.size) {
          const [key, next] = [...pending].sort((a, b) => a[1].at - b[1].at)[0];
          pending.delete(key); time = next.at; next.fn();
        }
      }
      assert.equal(pending.size, 0, 'timer loop did not settle');
    },
    async advance(ms) {
      const until = time + ms;
      for (let i = 0; i < 100; ++i) {
        await Promise.resolve();
        await Promise.resolve();
        const next = [...pending].filter(([, timer]) => timer.at <= until).sort((a, b) => a[1].at - b[1].at)[0];
        if (next) { pending.delete(next[0]); time = next[1].at; next[1].fn(); }
      }
      time = until;
    },
    fireInterval(delay) { for (const value of [...intervals.values()]) if (value.delay === delay) value.fn(); },
    intervals, get now() { return time; }, get pending() { return pending.size; }
  };
}
function contentHarness(count = 2) {
  const clock = timers(), events = new Map(), observers = [], intersections = [], messages = [], responses = [], sentAt = [];
  let defer = false, ok = true, reads = 0, scans = 0, cellSize = 20, intersectionFailure = '';
  const makeCells = (values) => values.map((value, x) => ({
    id: 'cell_' + x + '_0', nodeType: 1, value, visible: true, isConnected: true, style: {}, children: [],
    get ownerDocument() { return document; },
    get parentElement() { return this.parentNode?.nodeType === 1 ? this.parentNode : null; },
    get className() { ++reads; return this.value; },
    checkVisibility() {
      if (!this.visible || ['hidden', 'collapse'].includes(this.style.visibility)) return false;
      for (let node = this; node; node = node.parentElement)
        if (node.style.display === 'none' || node.style.opacity === '0' || node.style.contentVisibility === 'hidden') return false;
      return true;
    },
    getAnimations: () => [],
    matches() { return /^cell_/.test(this.id) && /(^|\s)cell(\s|$)/.test(this.value); },
    getBoundingClientRect: () => ({ left: x * cellSize, top: 50, right: (x + 1) * cellSize, bottom: 50 + cellSize, width: cellSize, height: cellSize })
  }));
  const cells = makeCells(new Array(count).fill('cell hd_closed'));
  const add = (name, fn) => { if (!events.has(name)) events.set(name, []); events.get(name).push(fn); };
  const document = {
    nodeType: 9, hidden: false, cells,
    querySelectorAll() { ++scans; return this.cells.filter((cell) => cell.matches()); },
    getElementById(id) { return this.cells.find((c) => c.id === id); },
    addEventListener: add
  };
  const container = (members, parentNode) => ({
    nodeType: 1, parentNode, style: {}, scrollLeft: 0, scrollTop: 0,
    get parentElement() { return this.parentNode?.nodeType === 1 ? this.parentNode : null; },
    getAnimations: () => [],
    getBoundingClientRect: () => ({ left: 0, top: 50, right: document.cells.length * cellSize, bottom: 50 + cellSize,
      width: document.cells.length * cellSize, height: cellSize }),
    querySelector() { return members.find((cell) => cell.matches()); },
    querySelectorAll() { return members.filter((cell) => cell.matches()); },
    contains(node) { for (; node; node = node.parentNode) if (node === this) return true; return false; }
  });
  document.documentElement = container(cells, document);
  document.body = container(cells, document.documentElement);
  for (const cell of cells) cell.parentNode = document.body;
  const notify = (mutations) => {
    for (const observer of observers) {
      const observed = (target) => {
        for (let node = target; node; node = observer.options.subtree ? node.parentNode : null)
          if (node === observer.target) return true;
        return false;
      };
      const removed = mutations.filter((mutation) => mutation.type === 'childList' && observed(mutation.target))
        .flatMap((mutation) => mutation.removedNodes);
      const records = mutations.filter((mutation) => {
        if (!observer.options[mutation.type] || (mutation.type === 'attributes' &&
            observer.options.attributeFilter && !observer.options.attributeFilter.includes(mutation.attributeName))) return false;
        if (observed(mutation.target)) return true;
        // A removed subtree remains observed until its removal record is delivered.
        if (observer.options.subtree) for (let node = mutation.target; node; node = node.parentNode)
          if (removed.includes(node)) return true;
        return false;
      });
      if (records.length) observer.fn(records);
    }
  };
  const runtime = { lastError: null, onMessage: event(), sendMessage(message, respond) {
    messages.push(JSON.parse(JSON.stringify(message)));
    sentAt.push(clock.now);
    if (defer) responses.push(respond); else respond({ ok });
  }};
  const storage = { local: { get: async () => ({ debug: false, mines_total: -1 }), set: async () => {} }, onChanged: event() };
  const window = { devicePixelRatio: 1, innerWidth: 1000, innerHeight: 800,
    visualViewport: { offsetLeft: 0, offsetTop: 0, scale: 1, addEventListener: add }, addEventListener: add,
    getComputedStyle: (element) => ({ transform: 'none', rotate: 'none', scale: 'none', perspective: 'none', offsetPath: 'none',
      overflowX: 'visible', overflowY: 'visible', contain: 'none', contentVisibility: 'visible', position: 'static',
      borderRadius: ['borderTopLeftRadius', 'borderTopRightRadius', 'borderBottomRightRadius', 'borderBottomLeftRadius']
        .map((key) => element.style[key] || '0px').join(' '), ...element.style }) };
  const context = { ...clock.api, console, document, window, chrome: { runtime, storage },
    MutationObserver: class {
      constructor(fn) { this.fn = fn; observers.push(this); }
      observe(target, options) { this.target = target; this.options = options; }
    },
    IntersectionObserver: class {
      constructor(fn) {
        if (intersectionFailure === 'constructor') throw new Error('IntersectionObserver unavailable');
        this.fn = fn; this.targets = []; this.records = []; this.disconnections = 0; intersections.push(this);
      }
      observe(target) {
        this.targets.push(target);
        if (intersectionFailure === 'observe') throw new Error('Cannot observe target');
      }
      disconnect() { ++this.disconnections; }
      acquire(entries = this.targets.map((target) => intersectionEntry(target)), time = clock.now) {
        assert.equal(this.disconnections, 0, 'disconnected observers cannot acquire new records');
        const records = entries.map((entry) => ({ time, ...entry }));
        this.records.push(...records);
        return records;
      }
      takeRecords() { return this.records.splice(0); }
      deliver(entries) {
        if (entries === undefined && !this.records.length && !this.disconnections) this.acquire();
        const queued = this.takeRecords();
        // Callback injection remains available after cancellation.
        this.fn((entries || (queued.length ? queued : this.targets.map((target) => intersectionEntry(target))))
          .map((entry) => ({ time: clock.now, ...entry })));
      }
    }
  };
  vm.runInNewContext(source('content.js'), context, { filename: 'content.js' });
  return {
    clock, document, window, runtime, messages, responses, sentAt, intersections, notify,
    set ok(value) { ok = value; }, set defer(value) { defer = value; }, set size(value) { cellSize = value; },
    set intersectionFailure(value) { intersectionFailure = value; },
    get reads() { return reads; }, get scans() { return scans; },
    fire(name, data = {}) { for (const fn of events.get(name) || []) fn(data); },
    change(index, value) {
      document.cells[index].value = value;
      notify([{ type: 'attributes', attributeName: 'class', target: document.cells[index] }]);
    },
    removeBoard() {
      const old = document.cells; document.cells = [];
      for (const cell of old) { cell.isConnected = false; cell.parentNode = null; }
      notify([{ type: 'childList', target: document.body, addedNodes: [], removedNodes: old }]);
    },
    replaceBody(values, { root = false, observed = true } = {}) {
      const oldCells = document.cells, old = root ? document.documentElement : document.body;
      const target = old.parentNode;
      for (const cell of oldCells) cell.isConnected = false;
      old.parentNode = null;
      document.cells = makeCells(values);
      if (root) document.documentElement = container(document.cells, document);
      document.body = container(document.cells, document.documentElement);
      for (const cell of document.cells) cell.parentNode = document.body;
      if (observed) notify([{ type: 'childList', target, addedNodes: [root ? document.documentElement : document.body], removedNodes: [old] }]);
      return oldCells;
    }
  };
}

function contentNode(parentNode, { nodeType = 1, id = '', value = '', data = '' } = {}) {
  const node = { parentNode, nodeType, id, value, data, children: [], style: {},
    get parentElement() { return this.parentNode?.nodeType === 1 ? this.parentNode : null; },
    getAnimations: () => [],
    matches() { return this.nodeType === 1 && /^cell_/.test(this.id) && /(^|\s)cell(\s|$)/.test(this.value); },
    contains(target) { for (; target; target = target.parentNode) if (target === this) return true; return false; },
    querySelectorAll() { return this.children.flatMap((child) => [...(child.matches?.() ? [child] : []), ...(child.querySelectorAll?.() || [])]); },
    querySelector() { return this.querySelectorAll()[0] || null; },
    remove() {
      const children = this.parentNode?.children;
      if (children) children.splice(children.indexOf(this), 1);
      this.parentNode = null;
    }
  };
  parentNode?.children?.push(node);
  return node;
}

const contentMutations = [
  ['text append', (cell) => () => ({ type: 'childList', target: cell,
    addedNodes: [contentNode(cell, { nodeType: 3, data: '12' })], removedNodes: [] })],
  ['text clear', (cell) => {
    const text = contentNode(cell, { nodeType: 3, data: '12' });
    return () => { text.remove(); return { type: 'childList', target: cell, addedNodes: [], removedNodes: [text] }; };
  }],
  ...[['empty text populated', '', '12'], ['text emptied', '12', '']].map(([name, before, after]) => [name, (cell) => {
    const text = contentNode(contentNode(cell), { nodeType: 3, data: before });
    return () => { text.data = after; return { type: 'characterData', target: text }; };
  }]),
  ...['class', 'style', 'id', 'data-state', 'aria-hidden'].map((attributeName) => ['nested ' + attributeName, (cell) => {
    const child = contentNode(contentNode(cell));
    return () => {
      if (attributeName === 'id') child.id = 'decorative-label';
      else if (attributeName === 'class') child.value = 'active';
      else if (attributeName === 'style') child.style.color = 'red';
      return { type: 'attributes', target: child, attributeName };
    };
  }]),
  ['nested child insertion', (cell) => {
    const child = contentNode(cell);
    return () => ({ type: 'childList', target: child, addedNodes: [contentNode(child)], removedNodes: [] });
  }],
  ['untracked cell lookalike', (cell) => {
    const child = contentNode(contentNode(cell, { id: 'cell_0_0_g2', value: 'cell' }));
    return () => ({ type: 'attributes', target: child, attributeName: 'class' });
  }],
  ['ancestor child insertion', (cell, h) => () => ({ type: 'childList', target: h.document.body,
    addedNodes: [contentNode(h.document.body, { value: 'clear' })], removedNodes: [] })],
  ['ancestor direct text mutation', (cell, h) => {
    const text = contentNode(h.document.body, { nodeType: 3, data: '' });
    return () => { text.data = 'label'; return { type: 'characterData', target: text }; };
  }]
];

function rectangle(left, top, width, height) {
  return { left, top, width, height, right: left + width, bottom: top + height };
}
function intersectionEntry(target, rect = target.getBoundingClientRect()) {
  return { target, boundingClientRect: { ...target.getBoundingClientRect() }, intersectionRect: rect,
    isIntersecting: rect.width > 0 && rect.height > 0 };
}
function clipping(message) {
  return ['clip_l', 'clip_t', 'clip_w', 'clip_h'].map((key) => message[key]);
}
function clippedGridHarness(width, height, { left = 100, top = 100, cellWidth = 20, cellHeight = 20 } = {}) {
  const h = contentHarness(width * height);
  h.logicalCells = h.document.cells.slice();
  h.logicalCells.forEach((cell, index) => {
    const x = index % width, y = Math.floor(index / width);
    cell.id = `cell_${x}_${y}`;
    cell.getBoundingClientRect = () => rectangle(left + x * cellWidth, top + y * cellHeight, cellWidth, cellHeight);
  });
  h.document.body.getBoundingClientRect = () => rectangle(left, top, width * cellWidth, height * cellHeight);
  h.document.body.style.overflowX = 'hidden';
  return h;
}
function clippedEntry(target, visible) {
  const bounds = target.getBoundingClientRect();
  const left = Math.max(bounds.left, visible.left), top = Math.max(bounds.top, visible.top);
  const right = Math.min(bounds.right, visible.right), bottom = Math.min(bounds.bottom, visible.bottom);
  return intersectionEntry(target, right > left && bottom > top ? rectangle(left, top, right - left, bottom - top) : rectangle(0, 0, 0, 0));
}

test('ordinary boards retain synchronous capture and clipping creates fresh bounded observations', async () => {
  const h = contentHarness();
  await h.clock.flush();
  assert.equal(h.intersections.length, 0, 'unclipped captures must not wait for a rendering update');
  assert.deepEqual(clipping(h.messages[0]), [undefined, undefined, undefined, undefined]);
  h.document.documentElement.style.overflowX = 'hidden';
  h.fire('resize');
  await h.clock.advance(16);
  assert.equal(h.messages.length, 1);
  const first = h.intersections[0];
  assert.deepEqual(first.targets, [h.document.body], 'bounds-matching container should replace per-cell observations');
  first.deliver([intersectionEntry(h.document.body, rectangle(5, 50, 30, 20))]);
  await h.clock.advance(0);
  assert.equal(first.disconnections, 1);
  assert.equal(h.clock.pending, 0, 'a completed query retained its timeout');
  assert.equal(h.messages[1].type, 'delta');
  assert.deepEqual(h.messages[1].updates, []);
  assert.deepEqual(clipping(h.messages[1]), [5, 50, 30, 20]);
  h.fire('scroll');
  await h.clock.advance(16);
  assert.equal(h.intersections.length, 2, 'each capture must obtain a fresh observation at unchanged intersection ratios');
  h.intersections[1].deliver([intersectionEntry(h.document.body, rectangle(10, 50, 30, 20))]);
  await h.clock.advance(0);
  assert.equal(h.messages[2].type, 'delta');
  assert.deepEqual(h.messages[2].updates, []);
  assert.deepEqual(clipping(h.messages[2]), [10, 50, 30, 20]);
  h.fire('scroll');
  await h.clock.advance(16);
  h.intersections[2].deliver([intersectionEntry(h.document.body, rectangle(10, 50, 30, 20))]);
  await h.clock.advance(0);
  assert.equal(h.messages.length, 3, 'unchanged clipping must not send a redundant delta');
});

test('clipping observations coalesce mutations and decode current cells after their rendering step', async () => {
  const h = contentHarness(16);
  h.document.documentElement.style.overflowY = 'auto';
  await h.clock.advance(0);
  assert.equal(h.intersections.length, 1);
  const reads = h.reads;
  for (let i = 0; i < 1000; ++i) h.change(3, 'cell hd_opened hd_type2');
  assert.equal(h.clock.pending, 1, 'only the observation timeout should remain scheduled');
  await h.clock.advance(30);
  assert.equal(h.intersections.length, 1, 'mutation bursts started parallel observations');
  assert.equal(h.messages.length, 0);
  h.intersections[0].deliver();
  await h.clock.advance(0);
  assert.deepEqual(h.messages[0].cells, Array.from({ length: 16 }, (_, i) => i === 3 ? 12 : 0));
  assert.equal(h.reads - reads, 1, 'post-observation capture repeatedly decoded dirty cells');
  assert.equal(h.intersections.length, 2, 'pending work did not resume after observation delivery');
  h.intersections[1].deliver();
  await h.clock.advance(0);
  assert.equal(h.messages.length, 1);
  assert.equal(h.clock.pending, 0);
});

test('cell mutations before acquisition retain fresh clipping even when timestamps are equal', async () => {
  for (const elapsed of [0, 1]) {
    const h = clippedGridHarness(3, 3);
    await h.clock.advance(0);
    h.change(4, 'cell hd_flag');
    const changedAt = h.clock.now;
    await h.clock.advance(elapsed);
    const observer = h.intersections[0], entries = observer.acquire();
    assert.equal(entries[0].time, changedAt + elapsed);
    await h.clock.advance(1);
    observer.deliver();
    await h.clock.advance(0);
    assert.equal(h.messages[0].cells[4], 2);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, 60, 60]);
  }
});

test('cell mutations after acquisition hide stale clipping even when timestamps are equal', async () => {
  for (const elapsed of [0, 1]) {
    const h = clippedGridHarness(3, 3);
    await h.clock.advance(0);
    const observer = h.intersections[0], entries = observer.acquire();
    await h.clock.advance(elapsed);
    h.change(4, 'cell hd_flag');
    assert.equal(observer.disconnections, 1, 'acquired records must cancel before their callback or timeout');
    observer.deliver(entries);
    await h.clock.advance(0);
    assert.equal(h.messages[0].cells[4], 2);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, 0, 0]);
    await h.clock.advance(16);
    h.intersections[1].deliver();
    await h.clock.advance(0);
    assert.deepEqual(h.messages[1].updates, []);
    assert.deepEqual(clipping(h.messages[1]), [100, 100, 60, 60]);
  }
});

test('cell class ABA invalidates acquisition even when the decoded board returns to its original value', async () => {
  const h = clippedGridHarness(3, 3);
  await h.clock.advance(0);
  const observer = h.intersections[0], entries = observer.acquire();
  h.change(4, 'cell hd_flag race_clip');
  h.change(4, 'cell hd_closed');
  observer.deliver(entries);
  await h.clock.advance(0);
  assert.deepEqual(h.messages[0].cells, new Array(9).fill(0));
  assert.deepEqual(clipping(h.messages[0]), [100, 100, 0, 0]);
  await h.clock.advance(16);
  h.intersections[1].deliver();
  await h.clock.advance(0);
  assert.deepEqual(h.messages[1].updates, []);
  assert.deepEqual(clipping(h.messages[1]), [100, 100, 60, 60]);
});

test('cell mutations between observation callback and publication invalidate its acquisition', async () => {
  for (const elapsed of [0, 1]) {
    const h = clippedGridHarness(3, 3);
    await h.clock.advance(0);
    const observer = h.intersections[0];
    observer.acquire();
    await h.clock.advance(elapsed);
    observer.deliver();
    assert.equal(observer.disconnections, 1);
    h.change(4, 'cell hd_flag');
    await h.clock.advance(0);
    assert.equal(h.messages[0].cells[4], 2);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, 0, 0]);
    await h.clock.advance(16);
    h.intersections[1].deliver();
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[1]), [100, 100, 60, 60]);
  }
});

test('a partial acquired batch cancels immediately and a late mixed batch cannot affect its replacement', async () => {
  for (const index of [0, 4, 8]) {
    const h = clippedGridHarness(3, 3);
    await h.clock.advance(0);
    const observer = h.intersections[0], old = observer.acquire([intersectionEntry(h.logicalCells[index])]);
    h.change(4, 'cell hd_flag');
    const changedAt = h.clock.now;
    assert.equal(observer.disconnections, 1);
    assert.equal(h.clock.pending, 0, 'partial acquisition retained the 250 ms observation timeout');
    await h.clock.advance(0);
    assert.equal(h.sentAt[0], changedAt, 'hiding must not wait for callback delivery or the observation timeout');
    assert.equal(h.messages[0].cells[4], 2);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, 0, 0]);
    await h.clock.advance(16);
    const current = h.intersections[1], entries = current.acquire();
    const mixed = entries.slice(); mixed[index] = old[0];
    observer.deliver(mixed);
    assert.equal(current.disconnections, 0, 'a retired callback canceled the replacement observation');
    assert.equal(h.messages.length, 1);
    current.deliver();
    await h.clock.advance(0);
    assert.deepEqual(h.messages[1].updates, []);
    assert.deepEqual(clipping(h.messages[1]), [100, 100, 60, 60]);
  }
});

test('resync isolates acquisition state without letting retired callbacks reset a new query', async () => {
  for (const phase of ['unchanged', 'before-retired-callback', 'after-retired-callback']) {
    const h = clippedGridHarness(3, 3);
    await h.clock.advance(0);
    const retired = h.intersections[0], old = retired.acquire();
    await h.clock.advance(16);
    h.change(4, 'cell hd_flag');
    h.runtime.onMessage.fire({ type: 'force_full' }, {}, () => {});
    await h.clock.advance(0);
    assert.equal(h.intersections.length, 2);
    const current = h.intersections[1], entries = current.acquire();
    if (phase === 'before-retired-callback') {
      await h.clock.advance(1);
      h.change(4, 'cell hd_closed');
    }
    retired.deliver(old);
    if (phase === 'after-retired-callback') h.change(4, 'cell hd_closed');
    current.deliver(entries);
    await h.clock.advance(0);
    const mutateNewQuery = phase !== 'unchanged';
    assert.equal(h.messages.length, 1);
    assert.equal(h.messages[0].cells[4], mutateNewQuery ? 0 : 2);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, mutateNewQuery ? 0 : 60, mutateNewQuery ? 0 : 60]);
  }
});

test('each query starts with fresh acquisition state after completion, cancellation, or timeout', async () => {
  for (const prior of ['completed', 'canceled', 'timeout']) {
    const h = clippedGridHarness(3, 3);
    await h.clock.advance(0);
    const retired = h.intersections[0], old = retired.acquire();
    if (prior === 'completed') retired.deliver();
    else if (prior === 'canceled') h.change(4, 'cell hd_flag');
    await h.clock.advance(prior === 'timeout' ? 250 : 0);
    assert.equal(h.messages.length, 1);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, prior === 'completed' ? 60 : 0, prior === 'completed' ? 60 : 0]);
    h.fire('scroll');
    await h.clock.advance(16);
    const current = h.intersections[1];
    assert.ok(current, prior + ' query did not allow another acquisition');
    retired.deliver(old);
    h.change(4, 'cell hd_opened hd_type3');
    const changedAt = h.clock.now, entries = current.acquire();
    assert.equal(entries[0].time, changedAt, 'fresh acquisition must be allowed within the same timestamp bucket');
    assert.equal(current.disconnections, 0, 'the previous query contaminated an empty new observation queue');
    current.deliver();
    await h.clock.advance(0);
    assert.equal(h.messages.length, 2);
    assert.deepEqual(h.messages[1].updates, [{ x: 1, y: 1, s: 13 }]);
    assert.deepEqual(clipping(h.messages[1]), [100, 100, 60, 60]);
  }
});

test('cell content and descendant mutations refresh geometry without decoding or rebuilding cells', async () => {
  for (const [name, prepare] of contentMutations) {
    const h = contentHarness(3), mutate = prepare(h.document.cells[1], h);
    await h.clock.flush();
    const reads = h.reads, scans = h.scans;
    h.size = 10;
    h.notify([mutate()]);
    await h.clock.flush();
    assert.equal(h.messages.length, 2, name + ' did not schedule a capture');
    assert.equal(h.messages[1].rect_w, 30, name + ' left stale geometry');
    assert.deepEqual(h.messages[1].updates, []);
    assert.equal(h.reads, reads, name + ' decoded unchanged classes');
    assert.equal(h.scans, scans, name + ' rebuilt unchanged membership');
  }
});

test('content mutations before acquisition coalesce into a fresh observation without hiding', async () => {
  for (const [name, prepare] of contentMutations) {
    const h = clippedGridHarness(3, 3), mutate = prepare(h.logicalCells[4], h);
    await h.clock.advance(0);
    const observer = h.intersections[0], reads = h.reads, scans = h.scans;
    h.notify([mutate()]);
    assert.equal(observer.disconnections, 0, name + ' canceled an unmeasured query');
    observer.acquire();
    observer.deliver();
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, 60, 60], name + ' hid a fresh measurement');
    await h.clock.advance(16);
    assert.equal(h.intersections.length, 2, name + ' lost pending content work');
    h.intersections[1].deliver();
    await h.clock.advance(0);
    assert.equal(h.messages.length, 1, name + ' sent redundant data');
    assert.equal(h.clock.pending, 0);
    assert.equal(h.reads, reads);
    assert.equal(h.scans, scans);
  }
});

test('content mutations after acquisition cancel queued clipping without decoding or rebuilding', async () => {
  for (const [name, prepare] of contentMutations) {
    const h = clippedGridHarness(3, 3), mutate = prepare(h.logicalCells[4], h);
    await h.clock.advance(0);
    const observer = h.intersections[0], entries = observer.acquire(), reads = h.reads, scans = h.scans;
    h.notify([mutate()]);
    assert.equal(observer.disconnections, 1, name + ' retained acquired clipping');
    assert.equal(h.clock.pending, 0, name + ' retained the observation timeout');
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, 0, 0]);
    await h.clock.advance(16);
    const current = h.intersections[1];
    observer.deliver(entries);
    assert.equal(current.disconnections, 0, 'retired content query canceled its replacement');
    current.deliver();
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[1]), [100, 100, 60, 60]);
    assert.deepEqual(h.messages[1].updates, []);
    assert.equal(h.reads, reads, name + ' decoded unchanged classes');
    assert.equal(h.scans, scans, name + ' rebuilt unchanged membership');
  }
});

test('content mutation bursts keep one pending query and no redundant cell decoding', async () => {
  const h = clippedGridHarness(3, 3), text = contentNode(contentNode(h.logicalCells[4]), { nodeType: 3 });
  await h.clock.advance(0);
  const reads = h.reads, scans = h.scans;
  for (let i = 0; i < 1000; ++i) {
    text.data = String(i);
    h.notify([{ type: 'characterData', target: text }]);
  }
  await h.clock.advance(30);
  assert.equal(h.intersections.length, 1);
  assert.equal(h.clock.pending, 1);
  h.intersections[0].deliver();
  await h.clock.advance(0);
  assert.equal(h.intersections.length, 2);
  h.intersections[1].deliver();
  await h.clock.advance(0);
  assert.equal(h.messages.length, 1);
  assert.deepEqual(clipping(h.messages[0]), [100, 100, 60, 60]);
  assert.equal(h.reads, reads);
  assert.equal(h.scans, scans);
  assert.equal(h.clock.pending, 0);
});

test('content mutations between callback completion and publication reject finished clipping', async () => {
  for (const [name, prepare] of contentMutations) {
    const h = clippedGridHarness(3, 3), mutate = prepare(h.logicalCells[4], h);
    await h.clock.advance(0);
    const observer = h.intersections[0], reads = h.reads, scans = h.scans;
    observer.acquire(); observer.deliver();
    h.notify([mutate()]);
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, 0, 0], name + ' accepted a completed stale query');
    assert.equal(h.reads, reads);
    assert.equal(h.scans, scans);
  }
});

test('detached content records retain their known-cell removal context and invalidate clipping', async () => {
  const h = clippedGridHarness(3, 3), child = contentNode(h.logicalCells[4]);
  const text = contentNode(child, { nodeType: 3, data: '12' });
  await h.clock.advance(0);
  const observer = h.intersections[0], reads = h.reads, scans = h.scans;
  observer.acquire();
  child.remove(); text.data = '';
  h.notify([{ type: 'childList', target: h.logicalCells[4], addedNodes: [], removedNodes: [child] },
    { type: 'characterData', target: text }]);
  assert.equal(text.parentNode, child);
  assert.equal(child.parentNode, null);
  assert.equal(observer.disconnections, 1);
  await h.clock.advance(0);
  assert.deepEqual(clipping(h.messages[0]), [100, 100, 0, 0]);
  assert.equal(h.reads, reads);
  assert.equal(h.scans, scans);
});

test('unknown cell and ancestor attributes refresh only geometry', async () => {
  for (const attributeName of ['style', 'hidden', 'data-state', 'aria-hidden']) for (const ancestor of [false, true]) {
    const h = contentHarness(3);
    await h.clock.flush();
    const reads = h.reads, scans = h.scans;
    h.size = 10;
    h.notify([{ type: 'attributes', attributeName, target: ancestor ? h.document.body : h.document.cells[1] }]);
    await h.clock.flush();
    assert.equal(h.messages[1]?.rect_w, 30, attributeName + ' failed to update geometry');
    assert.deepEqual(h.messages[1].updates, []);
    assert.equal(h.reads, reads, attributeName + ' decoded unchanged cell classes');
    assert.equal(h.scans, scans);
  }
});

test('unrelated widget and comment mutations leave acquired clipping and scheduling untouched', async () => {
  const h = clippedGridHarness(3, 3);
  const widget = contentNode(h.document.body), child = contentNode(widget), text = contentNode(child, { nodeType: 3 });
  const comment = contentNode(h.logicalCells[4], { nodeType: 8, data: 'comment' });
  await h.clock.advance(0);
  const observer = h.intersections[0], reads = h.reads, scans = h.scans;
  observer.acquire();
  const records = [
    { type: 'attributes', target: child, attributeName: 'data-state' },
    { type: 'attributes', target: widget, attributeName: 'class' },
    { type: 'characterData', target: text },
    { type: 'characterData', target: comment },
    { type: 'childList', target: child, addedNodes: [contentNode(child)], removedNodes: [] },
    { type: 'childList', target: h.document.body, addedNodes: [contentNode(h.document.body, { nodeType: 8 })], removedNodes: [] }
  ];
  for (let i = 0; i < 100; ++i) h.notify(records);
  assert.equal(observer.disconnections, 0);
  observer.deliver();
  await h.clock.advance(32);
  assert.equal(h.intersections.length, 1, 'unrelated subtree chatter scheduled another query');
  assert.equal(h.messages.length, 1);
  assert.deepEqual(clipping(h.messages[0]), [100, 100, 60, 60]);
  assert.equal(h.clock.pending, 0);
  assert.equal(h.reads, reads);
  assert.equal(h.scans, scans);
});

test('secondary-board mutations do not rebuild or recapture the primary board', async () => {
  const h = contentHarness(3), widget = contentNode(h.document.body), secondary = contentNode(widget);
  h.document.cells[2].id = 'cell_0_0_g2';
  h.document.cells[2].parentNode = secondary; secondary.children.push(h.document.cells[2]);
  await h.clock.flush();
  assert.equal(h.messages[0].w, 2);
  const reads = h.reads, scans = h.scans;
  h.notify([{ type: 'attributes', target: h.document.cells[2], attributeName: 'class' },
    { type: 'childList', target: widget, addedNodes: [secondary], removedNodes: [] }]);
  assert.equal(h.clock.pending, 0);
  await h.clock.flush();
  assert.equal(h.messages.length, 1);
  assert.equal(h.reads, reads);
  assert.equal(h.scans, scans);
});

test('primary membership discovery precedes descendant relevance and skips secondary lookalikes', async () => {
  for (const nested of [false, true]) {
    const h = contentHarness(3), cell = h.document.cells[2];
    cell.value = 'hd_closed';
    const wrapper = contentNode(nested ? h.document.cells[0] : contentNode(h.document.body));
    contentNode(wrapper, { id: 'cell_0_0_g2', value: 'cell' });
    cell.parentNode = wrapper; wrapper.children.push(cell);
    await h.clock.flush();
    assert.equal(h.messages[0].w, 2);
    const scans = h.scans;
    cell.value = 'cell hd_closed';
    h.notify([nested ? { type: 'attributes', attributeName: 'class', target: cell } :
      { type: 'childList', target: wrapper.parentNode, addedNodes: [wrapper], removedNodes: [] }]);
    await h.clock.flush();
    assert.equal(h.scans, scans + 1, 'new primary membership was treated as decoration');
    assert.equal(h.messages[1].w, 3);
    assert.deepEqual(h.messages[1].cells, [0, 0, 0]);
  }
});

test('replacing observed cells retires old target identities before publishing clipping', async () => {
  for (const root of [false, true]) {
    const h = contentHarness();
    h.document.documentElement.style.overflowX = 'hidden';
    await h.clock.advance(0);
    const stale = h.intersections[0];
    const entries = stale.targets.map((target) => intersectionEntry(target, rectangle(0, 50, 5, 20)));
    h.replaceBody(['cell hd_opened hd_type3', 'cell hd_closed'], { root });
    h.document.documentElement.style.overflowX = 'hidden';
    stale.deliver(entries);
    await h.clock.advance(16);
    assert.equal(h.messages.length, 1, 'the replacement must retire visible hints while waiting for fresh clipping');
    assert.deepEqual(h.messages[0].cells, [13, 0]);
    assert.deepEqual(clipping(h.messages[0]), [0, 50, 0, 0], 'old DOM target clipping was published for its replacement');
    assert.equal(stale.disconnections, 1);
    assert.equal(h.intersections.length, 2);
    assert.deepEqual(h.intersections[1].targets, [h.document.body]);
    h.intersections[1].deliver([intersectionEntry(h.document.body, rectangle(10, 50, 30, 20))]);
    await h.clock.advance(0);
    assert.equal(h.messages[1].type, 'delta');
    assert.deepEqual(h.messages[1].updates, []);
    assert.deepEqual(clipping(h.messages[1]), [10, 50, 30, 20]);
  }
});

test('resync cancels fresh observations and rejects late callbacks from the retired generation', async () => {
  for (const signal of ['force_full', 'pageshow']) {
    const h = contentHarness();
    h.document.documentElement.style.overflowX = 'hidden';
    await h.clock.advance(0);
    const stale = h.intersections[0];
    if (signal === 'force_full') h.runtime.onMessage.fire({ type: signal }, {}, () => {});
    else h.fire(signal);
    assert.equal(stale.disconnections, 1);
    stale.deliver([intersectionEntry(h.document.body, rectangle(0, 50, 1, 20))]);
    await h.clock.advance(16);
    assert.equal(h.messages.length, 0);
    assert.equal(h.intersections.length, 2);
    h.change(0, 'cell hd_closed hd_flag');
    await h.clock.advance(1);
    h.intersections[1].deliver();
    await h.clock.advance(0);
    assert.equal(h.messages[0].type, 'full');
    assert.deepEqual(h.messages[0].cells, [2, 0]);
    assert.deepEqual(clipping(h.messages[0]), [0, 50, 40, 20]);
    await h.clock.advance(16);
    h.intersections.at(-1).deliver();
    await h.clock.advance(0);
    assert.equal(h.clock.pending, 0);
  }
});

test('hiding cancels clipping queries and resuming creates a fresh full snapshot', async () => {
  const h = contentHarness();
  h.document.documentElement.style.overflowX = 'hidden';
  await h.clock.advance(0);
  const stale = h.intersections[0];
  h.document.hidden = true;
  h.fire('visibilitychange');
  assert.equal(stale.disconnections, 1);
  assert.equal(h.clock.pending, 0);
  h.change(0, 'cell hd_opened hd_type1');
  h.clock.fireInterval(1000);
  stale.deliver();
  await h.clock.advance(500);
  assert.equal(h.messages.length, 0);
  assert.equal(h.intersections.length, 1);
  h.document.hidden = false;
  h.fire('visibilitychange');
  await h.clock.advance(0);
  assert.equal(h.intersections.length, 2);
  h.intersections[1].deliver();
  await h.clock.advance(0);
  assert.equal(h.messages[0].type, 'full');
  assert.deepEqual(h.messages[0].cells, [11, 0]);
  assert.deepEqual(clipping(h.messages[0]), [0, 50, 40, 20]);
  assert.equal(h.clock.pending, 0);
});

test('clipping timeouts and observer failures hide stale presentation and permit a fresh retry', async () => {
  for (const failure of ['timeout', 'constructor', 'observe']) {
    const h = contentHarness();
    h.document.documentElement.style.overflowX = 'hidden';
    h.intersectionFailure = failure;
    await h.clock.advance(0);
    if (failure === 'timeout') {
      assert.equal(h.messages.length, 0);
      await h.clock.advance(249);
      assert.equal(h.messages.length, 0);
      await h.clock.advance(1);
    }
    assert.equal(h.messages.length, 1, failure);
    assert.equal(h.messages[0].w, 2, 'observation failure must retain valid solver input');
    assert.deepEqual(h.messages[0].cells, [0, 0]);
    assert.deepEqual(clipping(h.messages[0]), [0, 50, 0, 0]);
    assert.equal(h.clock.pending, 0);
    if (h.intersections.length) assert.equal(h.intersections[0].disconnections, 1);
    h.intersectionFailure = '';
    h.clock.fireInterval(1000);
    await h.clock.advance(16);
    h.intersections.at(-1).deliver();
    await h.clock.advance(0);
    assert.equal(h.messages[1].type, 'delta');
    assert.deepEqual(h.messages[1].updates, []);
    assert.deepEqual(clipping(h.messages[1]), [0, 50, 40, 20]);
    assert.equal(h.clock.pending, 0);
  }
});

test('layout changes during observation retry instead of applying stale intersection geometry', async () => {
  for (const change of [
    (h) => { h.document.documentElement.scrollLeft = 10; },
    (h) => { h.window.innerWidth = 800; },
    (h) => { h.window.visualViewport.offsetLeft = 5; },
    (h) => { h.size = 10; }
  ]) {
    const h = contentHarness();
    h.document.documentElement.style.overflowX = 'hidden';
    await h.clock.advance(0);
    const stale = h.intersections[0];
    const entries = stale.targets.map((target) => intersectionEntry(target));
    change(h);
    stale.deliver(entries);
    await h.clock.advance(16);
    assert.equal(h.messages.length, 1);
    assert.deepEqual(h.messages[0].cells, [0, 0]);
    assert.deepEqual(clipping(h.messages[0]), [0, 50, 0, 0], 'stale layout geometry retained visible hints');
    assert.equal(h.messages[0].rect_w, h.document.cells[1].getBoundingClientRect().right);
    assert.equal(h.messages[0].vv_x, h.window.visualViewport.offsetLeft);
    assert.equal(h.intersections.length, 2, 'changed geometry left capture permanently in flight');
    h.intersections[1].deliver();
    await h.clock.advance(0);
    assert.equal(h.messages.length, 2);
    assert.equal(h.messages[1].type, 'delta');
    assert.deepEqual(h.messages[1].updates, []);
    assert.deepEqual(clipping(h.messages[1]), [0, 50, h.messages[0].rect_w, h.messages[0].rect_h]);
    assert.equal(h.clock.pending, 0);
  }
});

test('partial, duplicate, or stale observation records cannot supply clipping for current targets', async () => {
  for (const invalid of ['partial', 'duplicate', 'foreign', 'bounds']) {
    const h = contentHarness();
    h.document.body.style.overflowX = 'hidden';
    await h.clock.advance(0);
    const stale = h.intersections[0];
    const entries = stale.targets.map((target) => intersectionEntry(target));
    if (invalid === 'partial') entries.pop();
    else if (invalid === 'duplicate') entries[1] = entries[0];
    else if (invalid === 'foreign') entries[0].target = { getBoundingClientRect: () => entries[0].boundingClientRect };
    else entries[0].boundingClientRect = rectangle(1, 50, 20, 20);
    stale.deliver(entries);
    await h.clock.advance(16);
    assert.equal(h.messages.length, 1, invalid);
    assert.deepEqual(h.messages[0].cells, [0, 0]);
    assert.deepEqual(clipping(h.messages[0]), [0, 50, 0, 0], invalid);
    assert.equal(h.intersections.length, 2);
    h.intersections[1].deliver();
    await h.clock.advance(0);
    assert.equal(h.messages[1].type, 'delta');
    assert.deepEqual(h.messages[1].updates, []);
    assert.deepEqual(clipping(h.messages[1]), [0, 50, 40, 20]);
    assert.equal(h.clock.pending, 0);
  }
});

test('scrolling away and back during observation invalidates clipping even when geometry returns unchanged', async () => {
  const h = contentHarness();
  h.document.documentElement.style.overflowX = 'hidden';
  await h.clock.advance(0);
  h.intersections[0].deliver([intersectionEntry(h.document.body, rectangle(5, 50, 30, 20))]);
  await h.clock.advance(0);
  h.fire('scroll');
  await h.clock.advance(16);
  const stale = h.intersections[1];
  const entries = [intersectionEntry(h.document.body, rectangle(5, 50, 30, 20))];
  h.document.documentElement.scrollLeft = 10;
  h.fire('scroll');
  h.document.documentElement.scrollLeft = 0;
  h.fire('scroll');
  stale.deliver(entries);
  await h.clock.advance(16);
  assert.equal(h.messages[1].type, 'delta');
  assert.deepEqual(h.messages[1].updates, []);
  assert.deepEqual(clipping(h.messages[1]), [0, 50, 0, 0], 'an unchanged final layout reused an observation from before intervening scrolling');
  assert.equal(h.intersections.length, 3);
  h.intersections[2].deliver([intersectionEntry(h.document.body, rectangle(5, 50, 30, 20))]);
  await h.clock.advance(0);
  assert.deepEqual(clipping(h.messages[2]), [5, 50, 30, 20]);
  assert.deepEqual(h.messages[2].updates, []);
  assert.equal(h.clock.pending, 0);
});

test('identical invalid clipping records stop immediate retries and recover through the fallback', async () => {
  const h = contentHarness();
  h.document.body.style.overflowX = 'hidden';
  await h.clock.advance(0);
  const invalid = () => [intersectionEntry(h.document.cells[0]), intersectionEntry(h.document.cells[0])];
  h.intersections[0].deliver(invalid());
  await h.clock.advance(16);
  assert.equal(h.messages.length, 1);
  assert.deepEqual(h.messages[0].cells, [0, 0]);
  assert.deepEqual(clipping(h.messages[0]), [0, 50, 0, 0]);
  assert.equal(h.intersections.length, 2, 'the first distinct failure did not request a fresh observation');
  h.intersections[1].deliver(invalid());
  await h.clock.advance(1000);
  assert.equal(h.intersections.length, 2, 'identical bad records caused an endless observation loop');
  assert.equal(h.messages.length, 1);
  assert.equal(h.clock.pending, 0);
  h.clock.fireInterval(1000);
  await h.clock.advance(0);
  assert.equal(h.intersections.length, 3);
  h.intersections[2].deliver(invalid());
  await h.clock.advance(16);
  assert.equal(h.intersections.length, 3, 'unchanged fallback failure scheduled an immediate retry');
  assert.equal(h.clock.pending, 0);
  h.clock.fireInterval(1000);
  await h.clock.advance(0);
  h.intersections[3].deliver();
  await h.clock.advance(0);
  assert.equal(h.messages[1].type, 'delta');
  assert.deepEqual(h.messages[1].updates, []);
  assert.deepEqual(clipping(h.messages[1]), [0, 50, 40, 20]);
  assert.equal(h.clock.pending, 0);
});

test('fractional observer bounds compare their rounded edges rather than a separately derived width', async () => {
  const h = contentHarness();
  h.document.documentElement.style.overflowX = 'hidden';
  h.document.cells.forEach((cell, i) => { cell.getBoundingClientRect = () => rectangle(140.8 + i * 20.1, 50.3, 20.1, 20.1); });
  h.document.body.getBoundingClientRect = () => rectangle(140.8, 50.3, 40.2, 20.1);
  await h.clock.advance(0);
  assert.deepEqual(h.intersections[0].targets, [h.document.body]);
  const entry = intersectionEntry(h.document.body);
  const bounds = entry.boundingClientRect;
  const rounded = Object.fromEntries(['left', 'top', 'right', 'bottom'].map((key) => [key, Math.fround(bounds[key])]));
  rounded.width = rounded.right - rounded.left;
  rounded.height = rounded.bottom - rounded.top;
  assert.ok(Math.abs(rounded.width - bounds.width) > 1e-7, 'fixture did not reproduce independently rounded width drift');
  entry.boundingClientRect = rounded;
  h.intersections[0].deliver([entry]);
  await h.clock.advance(0);
  assert.equal(h.messages[0].w, 2);
  assert.ok(h.messages[0].clip_w > 40 && h.messages[0].clip_h > 20, 'representational rounding hid a valid board');
  assert.equal(h.intersections.length, 1);
  assert.equal(h.clock.pending, 0);
});

test('observer edge tolerances accept two float32 steps and reject larger or uncapped coordinate drift', async () => {
  const cases = [
    { origin: 100.25, drift: 2 ** -16, accepted: true },
    { origin: 100.25, drift: 3 * 2 ** -17, accepted: false },
    { origin: 100000.25, drift: .0009, accepted: true },
    { origin: 100000.25, drift: .0011, accepted: false }
  ];
  for (const { origin, drift, accepted } of cases) for (const edge of ['left', 'top', 'right', 'bottom']) {
    const h = contentHarness();
    h.document.documentElement.style.overflowX = 'hidden';
    h.document.cells.forEach((cell, i) => { cell.getBoundingClientRect = () => rectangle(origin + i * 2, origin, 2, 2); });
    h.document.body.getBoundingClientRect = () => rectangle(origin, origin, 4, 2);
    await h.clock.advance(0);
    const entry = intersectionEntry(h.document.body);
    entry.boundingClientRect[edge] += drift;
    entry.boundingClientRect.width = entry.boundingClientRect.right - entry.boundingClientRect.left;
    entry.boundingClientRect.height = entry.boundingClientRect.bottom - entry.boundingClientRect.top;
    h.intersections[0].deliver([entry]);
    await h.clock.advance(0);
    assert.equal(h.messages[0].clip_w > 0, accepted, `${origin} ${edge} drift ${drift}`);
    if (accepted) assert.equal(h.clock.pending, 0);
    else {
      assert.deepEqual(clipping(h.messages[0]), [origin, origin, 0, 0]);
      await h.clock.advance(16);
      h.intersections[1].deliver();
      await h.clock.advance(0);
      assert.deepEqual(clipping(h.messages[1]), [origin, origin, 4, 2]);
      assert.equal(h.clock.pending, 0);
    }
  }
});

test('late timed-out callbacks cannot cancel a replacement query and removed boards clear once', async () => {
  const h = contentHarness();
  h.document.documentElement.style.overflowX = 'hidden';
  await h.clock.advance(250);
  const expired = h.intersections[0];
  assert.deepEqual(clipping(h.messages[0]), [0, 50, 0, 0]);
  h.fire('resize');
  await h.clock.advance(16);
  assert.equal(h.intersections.length, 2);
  expired.deliver();
  await h.clock.advance(0);
  assert.equal(h.messages.length, 1);
  assert.equal(h.clock.pending, 1, 'retired callback cleared the current observation timeout');
  h.removeBoard();
  h.intersections[1].deliver();
  await h.clock.advance(16);
  assert.equal(h.messages[1].w, 0);
  assert.deepEqual(h.messages[1].cells, []);
  assert.equal(h.messages.length, 2);
  assert.equal(h.clock.pending, 0);
});

test('removing clipping during observation restores ordinary metadata without stale clip fields', async () => {
  const h = contentHarness();
  h.document.documentElement.style.overflowX = 'hidden';
  await h.clock.advance(0);
  h.intersections[0].deliver([intersectionEntry(h.document.body, rectangle(5, 50, 20, 20))]);
  await h.clock.advance(0);
  h.fire('scroll');
  await h.clock.advance(16);
  h.document.documentElement.style.overflowX = 'visible';
  h.intersections[1].deliver([intersectionEntry(h.document.body, rectangle(5, 50, 20, 20))]);
  await h.clock.advance(0);
  assert.equal(h.messages[1].type, 'delta');
  assert.deepEqual(h.messages[1].updates, []);
  assert.deepEqual(clipping(h.messages[1]), [undefined, undefined, undefined, undefined]);
  assert.equal(h.clock.pending, 0);
});

test('unsafe container candidates fall back to cells and reject rectangle-union holes', async () => {
  for (const candidate of ['clipping', 'mismatched', 'absolute', 'fixed']) {
    const h = contentHarness();
    h.document.documentElement.style.overflowX = 'hidden';
    if (candidate === 'clipping') h.document.body.style.overflowX = 'hidden';
    else if (candidate === 'mismatched') h.document.body.getBoundingClientRect = () => rectangle(-5, 45, 50, 30);
    else h.document.cells[0].style.position = candidate;
    await h.clock.advance(0);
    assert.deepEqual(h.intersections[0].targets, h.document.cells, candidate);
    h.intersections[0].deliver(h.document.cells.map((target, i) => intersectionEntry(target, rectangle(i ? 20 : 5, 50, 15, 20))));
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[0]), [5, 50, 30, 20], 'adjacent intersections must retain their visible rectangle');
    h.fire('scroll');
    await h.clock.advance(16);
    h.intersections[1].deliver(h.document.cells.map((target, i) => intersectionEntry(target, rectangle(i * 20, 50, 15, 20))));
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[1]), [0, 50, 0, 0], 'a horizontal hole must not be filled by the bounding rectangle');
    h.fire('scroll');
    await h.clock.advance(16);
    h.intersections[2].deliver(h.document.cells.map((target, i) => intersectionEntry(target, rectangle(i * 20, 50 + i * 10, 20, 10))));
    await h.clock.advance(0);
    assert.equal(h.messages.length, 2, 'a diagonal hole must also keep the clipping empty');
    h.fire('scroll');
    await h.clock.advance(16);
    h.intersections[3].deliver(h.document.cells.map((target, i) => intersectionEntry(target, rectangle(i * 20, 50, i ? 20 : 0, 20))));
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[2]), [20, 50, 20, 20], 'a fully clipped cell must not erase an adjacent visible cell');
    assert.deepEqual(h.messages[2].updates, []);
  }
});

test('cached observation targets never bypass changed cell positioning at publication', async () => {
  const h = contentHarness();
  h.document.documentElement.style.overflowX = 'hidden';
  await h.clock.advance(0);
  h.intersections[0].deliver();
  await h.clock.advance(0);
  h.document.cells[0].style.position = 'fixed';
  h.change(0, 'cell hd_flag');
  await h.clock.advance(16);
  assert.deepEqual(h.intersections[1].targets, [h.document.body], 'the prior query target should remain only a provisional choice');
  h.intersections[1].deliver();
  await h.clock.advance(16);
  assert.deepEqual(clipping(h.messages[1]), [0, 50, 0, 0], 'old parent clipping must not be published for escaped children');
  assert.deepEqual(h.messages[1].updates, [{ x: 0, y: 0, s: 2 }]);
  assert.deepEqual(h.intersections[2].targets, h.document.cells);
  h.intersections[2].deliver();
  await h.clock.advance(0);
  assert.deepEqual(clipping(h.messages[2]), [0, 50, 40, 20]);
  assert.deepEqual(h.messages[2].updates, []);
  assert.equal(h.clock.pending, 0);
});

test('clipped grids validate hidden and displaced interior cells even after a parent target is cached', async () => {
  for (const defect of ['moved', 'gap', 'display', 'visibility', 'opacity', 'rounded']) {
    const h = contentHarness(3);
    h.document.documentElement.style.overflowX = 'hidden';
    await h.clock.advance(0);
    h.intersections[0].deliver();
    await h.clock.advance(0);
    const cell = h.document.cells[1], original = cell.getBoundingClientRect;
    if (defect === 'moved') cell.getBoundingClientRect = () => rectangle(40, 50, 20, 20);
    if (defect === 'gap') cell.getBoundingClientRect = () => rectangle(20, 50, 19.998, 20);
    if (defect === 'display') cell.style.display = 'none';
    if (defect === 'visibility') cell.style.visibility = 'hidden';
    if (defect === 'opacity') cell.style.opacity = '0';
    if (defect === 'rounded') Object.assign(cell.style, { overflowX: 'hidden', borderTopLeftRadius: '4px' });
    h.change(1, 'cell hd_flag');
    await h.clock.advance(16);
    assert.deepEqual(h.intersections[1].targets, [h.document.body]);
    h.intersections[1].deliver();
    await h.clock.advance(16);
    assert.equal(h.messages[1].w, 0, defect + ' left hints over an unsupported interior grid');
    assert.equal(h.clock.pending, 0);
    cell.getBoundingClientRect = original; cell.style = {};
    h.change(1, 'cell hd_opened hd_type1');
    await h.clock.advance(16);
    h.intersections[2].deliver();
    await h.clock.advance(0);
    assert.equal(h.messages[2].w, 3);
    assert.deepEqual(h.messages[2].cells, [0, 11, 0]);
    assert.deepEqual(clipping(h.messages[2]), [0, 50, 60, 20]);
    assert.equal(h.clock.pending, 0);
  }
});

test('clipped grid adjacency follows logical coordinates when DOM order differs', async () => {
  const h = contentHarness(3);
  h.document.documentElement.style.overflowX = 'hidden';
  h.document.cells.reverse();
  await h.clock.advance(0);
  h.intersections[0].deliver();
  await h.clock.advance(0);
  assert.equal(h.messages[0].w, 3);
  assert.deepEqual(clipping(h.messages[0]), [0, 50, 60, 20]);
  assert.equal(h.clock.pending, 0);
});

test('two-dimensional clipping follows cell identities across shuffled DOM and observer order', async () => {
  const h = clippedGridHarness(3, 2);
  h.logicalCells[1].value = 'cell hd_closed hd_flag';
  h.logicalCells[4].value = 'cell hd_opened hd_type1';
  h.document.cells = [5, 1, 3, 0, 4, 2].map((index) => h.logicalCells[index]);
  await h.clock.advance(0);
  assert.deepEqual(h.intersections[0].targets, h.document.cells);
  const visible = rectangle(105, 105, 50, 30);
  h.intersections[0].deliver(h.logicalCells.map((cell) => clippedEntry(cell, visible)).reverse());
  await h.clock.advance(0);
  assert.deepEqual(h.messages[0].cells, [0, 2, 0, 0, 11, 0]);
  assert.deepEqual(clipping(h.messages[0]), [105, 105, 50, 30]);
  h.fire('scroll');
  await h.clock.advance(16);
  const narrower = rectangle(122, 102, 28, 36);
  h.intersections[1].deliver([2, 5, 0, 4, 1, 3].map((index) => clippedEntry(h.logicalCells[index], narrower)));
  await h.clock.advance(0);
  assert.equal(h.messages[1].type, 'delta');
  assert.deepEqual(h.messages[1].updates, []);
  assert.deepEqual(clipping(h.messages[1]), [122, 102, 28, 36], 'fully clipped targets must not change logical grid alignment');
  assert.equal(h.clock.pending, 0);
});

test('individually sub-tolerance observation edges cannot fill a real gap between visible cell bands', async () => {
  for (const axis of ['horizontal', 'vertical']) {
    const h = clippedGridHarness(2, 2);
    await h.clock.advance(0);
    const entries = h.logicalCells.map((cell, index) => {
      const bounds = cell.getBoundingClientRect();
      const second = axis === 'horizontal' ? index % 2 : Math.floor(index / 2);
      // Each edge differs by 0.00001 CSS px; together the hole exceeds the float32 allowance at 120px.
      const rect = axis === 'horizontal' ? rectangle(bounds.left + (second ? .00001 : 0), bounds.top, bounds.width - .00001, bounds.height) :
        rectangle(bounds.left, bounds.top + (second ? .00001 : 0), bounds.width, bounds.height - .00001);
      return intersectionEntry(cell, rect);
    });
    h.intersections[0].deliver(entries.reverse());
    await h.clock.advance(0);
    assert.equal(h.messages[0].w, 2); assert.equal(h.messages[0].h, 2);
    assert.deepEqual(h.messages[0].cells, [0, 0, 0, 0]);
    assert.deepEqual(clipping(h.messages[0]), [100, 100, 0, 0], axis + ' gap was replaced by its bounding rectangle');
    h.fire('scroll');
    await h.clock.advance(16);
    h.intersections[1].deliver();
    await h.clock.advance(0);
    assert.deepEqual(clipping(h.messages[1]), [100, 100, 40, 40]);
    assert.deepEqual(h.messages[1].updates, []);
    assert.equal(h.clock.pending, 0);
  }
});

test('missing interior cells or complete rows cannot be mistaken for a contiguous visible grid', async () => {
  for (const missing of ['center', 'row', 'edge-column']) {
    const h = clippedGridHarness(3, 3);
    await h.clock.advance(0);
    h.intersections[0].deliver(h.logicalCells.map((cell, index) => {
      const hidden = missing === 'center' ? index === 4 : missing === 'row' ? Math.floor(index / 3) === 1 : index % 3 === 0;
      return intersectionEntry(cell, hidden ? rectangle(0, 0, 0, 0) : cell.getBoundingClientRect());
    }).reverse());
    await h.clock.advance(0);
    assert.equal(h.messages[0].w, 3); assert.equal(h.messages[0].h, 3);
    assert.deepEqual(clipping(h.messages[0]), missing === 'edge-column' ? [120, 100, 40, 60] : [100, 100, 0, 0], missing);
    assert.equal(h.clock.pending, 0);
  }
});

test('rounded float32 cell observations preserve a valid fractional two-dimensional intersection', async () => {
  const h = clippedGridHarness(3, 2, { left: 140.8, top: 50.3, cellWidth: 20.1, cellHeight: 20.1 });
  await h.clock.advance(0);
  const visible = rectangle(145.25, 55.25, 50, 30);
  const rounded = (rect) => {
    const result = Object.fromEntries(['left', 'top', 'right', 'bottom'].map((key) => [key, Math.fround(rect[key])]));
    result.width = result.right - result.left; result.height = result.bottom - result.top;
    return result;
  };
  const entries = h.logicalCells.map((cell) => {
    const entry = clippedEntry(cell, visible);
    entry.boundingClientRect = rounded(entry.boundingClientRect);
    entry.intersectionRect = rounded(entry.intersectionRect);
    return entry;
  });
  assert.notEqual(entries[0].intersectionRect.right, h.logicalCells[0].getBoundingClientRect().right,
    'fixture must exercise an observed edge that differs from the fresh layout rectangle');
  h.intersections[0].deliver([entries[4], entries[0], entries[5], entries[2], entries[3], entries[1]]);
  await h.clock.advance(0);
  assert.equal(h.messages[0].w, 3); assert.equal(h.messages[0].h, 2);
  assert.deepEqual(clipping(h.messages[0]), [145.25, 55.25, 50, 30]);
  assert.equal(h.clock.pending, 0);
});

test('a tolerated fresh-grid gap remains a hole when only its narrow neighboring fragments are visible', async () => {
  const h = clippedGridHarness(2, 1);
  h.logicalCells[0].getBoundingClientRect = () => rectangle(100, 100, 19.999995, 20);
  h.logicalCells[1].getBoundingClientRect = () => rectangle(120.000005, 100, 19.999995, 20);
  await h.clock.advance(0);
  h.intersections[0].deliver();
  await h.clock.advance(0);
  assert.equal(h.messages[0].w, 2, 'the existing fresh-layout precision allowance must remain supported');
  assert.deepEqual(clipping(h.messages[0]), [100, 100, 40, 20]);
  h.fire('scroll');
  await h.clock.advance(16);
  const narrow = rectangle(119.999983, 100, .000034, 20);
  h.intersections[1].deliver(h.logicalCells.map((cell) => clippedEntry(cell, narrow)));
  await h.clock.advance(0);
  assert.equal(h.messages[1].type, 'delta');
  assert.deepEqual(h.messages[1].updates, []);
  assert.deepEqual(clipping(h.messages[1]), [100, 100, 0, 0], 'layout tolerance cannot establish an exact partition across a real gap');
  h.fire('scroll');
  await h.clock.advance(16);
  h.intersections[2].deliver();
  await h.clock.advance(0);
  assert.deepEqual(clipping(h.messages[2]), [100, 100, 40, 20]);
  assert.equal(h.clock.pending, 0);
});

test('visible cells may override an ancestor visibility rule but not ancestor opacity', async () => {
  const h = contentHarness(3);
  h.document.documentElement.style.overflowX = 'hidden';
  h.document.body.style.visibility = 'hidden';
  for (const cell of h.document.cells) cell.style.visibility = 'visible';
  await h.clock.advance(0);
  h.intersections[0].deliver();
  await h.clock.advance(0);
  assert.equal(h.messages[0].w, 3);
  assert.deepEqual(clipping(h.messages[0]), [0, 50, 60, 20]);
  h.document.body.style.opacity = '0';
  h.fire('resize'); await h.clock.advance(16);
  assert.equal(h.messages[1].w, 0);
  assert.equal(h.clock.pending, 0);
});

test('failed delivery is retried as a full snapshot', async () => {
  const h = contentHarness();
  h.ok = false;
  await h.clock.flush();
  assert.equal(h.messages[0].type, 'full');
  h.ok = true;
  h.change(0, 'cell hd_opened hd_type1');
  await h.clock.flush();
  assert.equal(h.messages[1].type, 'full');
  assert.deepEqual(h.messages[1].cells, [11, 0]);
  h.ok = false;
  h.change(1, 'cell hd_closed hd_flag');
  await h.clock.flush();
  assert.equal(h.messages.at(-1).type, 'delta');
  h.ok = true;
  h.clock.fireInterval(1000);
  await h.clock.flush();
  assert.equal(h.messages.at(-1).type, 'full');
  assert.deepEqual(h.messages.at(-1).cells, [11, 2]);
});

test('display-only flag backgrounds and incomplete clues are decoded safely', async () => {
  const h = contentHarness(3);
  h.document.cells[0].value = 'cell hd_closed_flag';
  h.document.cells[1].value = 'cell hd_closed_flag hd_flag';
  h.document.cells[2].value = 'cell hd_opened hd_type2';
  await h.clock.flush();
  assert.deepEqual(h.messages[0].cells, [0, 2, 12]);
  h.change(2, 'cell hd_opened');
  await h.clock.flush();
  assert.equal(h.messages.at(-1).w, 0, 'an incomplete clue must hide hints');
});

test('terminal boards and removed boards clear the overlay once', async () => {
  for (const kind of ['bomb10', 'bomb11', 'removed']) {
    const h = contentHarness();
    await h.clock.flush();
    if (kind === 'removed') h.removeBoard();
    else h.change(0, 'cell hd_opened hd_type' + kind.slice(4));
    await h.clock.flush();
    assert.equal(h.messages.at(-1).type, 'full');
    assert.equal(h.messages.at(-1).w, 0);
    assert.deepEqual(h.messages.at(-1).cells, []);
    const count = h.messages.length;
    h.clock.fireInterval(1000); await h.clock.flush();
    assert.equal(h.messages.length, count);
  }
});

test('mutations only decode changed cells and geometry can shrink', async () => {
  const h = contentHarness(16);
  await h.clock.flush();
  const before = h.reads;
  h.change(7, 'cell hd_opened hd_type2');
  await h.clock.flush();
  assert.equal(h.reads - before, 1);
  assert.deepEqual(h.messages.at(-1).updates, [{ x: 7, y: 0, s: 12 }]);
  h.size = 10;
  h.fire('resize'); await h.clock.flush();
  assert.equal(h.messages.at(-1).rect_w, 160);
  assert.equal(h.messages.at(-1).cell_px, 10);
  assert.deepEqual(h.messages.at(-1).updates, []);
});

test('CSS animation and transition lifecycle events refresh geometry without rescanning cells', async () => {
  for (const name of ['animationstart', 'animationend', 'animationcancel', 'transitionrun', 'transitionstart', 'transitionend', 'transitioncancel']) {
    const h = contentHarness(3);
    await h.clock.flush();
    const reads = h.reads;
    h.size = 10;
    h.fire(name, { target: h.document.body, pseudoElement: '' });
    await h.clock.flush();
    assert.equal(h.messages.at(-1).rect_w, 30, name + ' must refresh geometry');
    assert.equal(h.reads, reads, name + ' must not decode unchanged cells');
  }
});

test('animation events follow current cells and their subtrees but ignore unrelated and pseudo events', async () => {
  const targets = [
    (h) => h.document.cells[1],
    (h) => h.document.body,
    (h) => ({ closest: () => h.document.cells[1] }),
    (h) => ({ querySelector: () => h.document.cells[1] })
  ];
  for (const target of targets) {
    const h = contentHarness(3);
    await h.clock.flush();
    h.size = 10;
    h.fire('animationstart', { target: target(h), pseudoElement: '' });
    await h.clock.flush();
    assert.equal(h.messages.at(-1).rect_w, 30);
  }
  const h = contentHarness();
  await h.clock.flush();
  h.size = 10;
  h.fire('animationstart', { target: h.document.body, pseudoElement: '::before' });
  h.fire('transitionend', { target: { contains: () => false, querySelector: () => null } });
  h.fire('animationcancel', { target: { closest: () => ({}) } });
  h.fire('animationiteration', { target: h.document.body });
  assert.equal(h.clock.pending, 0);
  await h.clock.flush();
  assert.equal(h.messages.length, 1);
});

test('animation event bursts retain the frame cap, delivery coalescing, and hidden-tab pause', async () => {
  const h = contentHarness();
  await h.clock.flush();
  await h.clock.advance(4);
  h.size = 10;
  for (let i = 0; i < 1000; ++i) h.fire('transitionrun', { target: h.document.body });
  assert.equal(h.clock.pending, 1);
  await h.clock.advance(11);
  assert.equal(h.messages.length, 1);
  await h.clock.advance(1);
  assert.equal(h.messages.length, 2);
  assert.equal(h.sentAt[1] - h.sentAt[0], 16);
  h.defer = true;
  h.size = 15;
  h.fire('animationstart', { target: h.document.body });
  await h.clock.advance(16);
  assert.equal(h.messages.length, 3);
  h.size = 12;
  for (let i = 0; i < 1000; ++i) h.fire('animationcancel', { target: h.document.body });
  assert.equal(h.clock.pending, 0);
  await h.clock.advance(30);
  const acknowledgedAt = h.clock.now;
  h.responses.shift()({ ok: true });
  await h.clock.advance(0);
  assert.equal(h.messages.length, 4);
  assert.equal(h.sentAt[3], acknowledgedAt);
  assert.equal(h.messages.at(-1).rect_w, 24);
  h.responses.shift()({ ok: true });
  await h.clock.flush();
  h.document.hidden = true;
  h.fire('visibilitychange');
  h.size = 20;
  h.fire('transitionend', { target: h.document.body });
  assert.equal(h.clock.pending, 0);
  assert.equal(h.messages.length, 4);
});

test('replacing the body or document root rebuilds cells and keeps observing the new board', async () => {
  for (const root of [false, true]) {
    const h = contentHarness();
    await h.clock.flush();
    h.replaceBody(['cell hd_opened hd_type2', 'cell hd_closed hd_flag'], { root });
    await h.clock.flush();
    assert.deepEqual(h.messages.at(-1).updates, [{ x: 0, y: 0, s: 12 }, { x: 1, y: 0, s: 2 }]);
    const reads = h.reads;
    h.change(1, 'cell hd_opened hd_type1');
    await h.clock.flush();
    assert.deepEqual(h.messages.at(-1).updates, [{ x: 1, y: 0, s: 11 }]);
    assert.equal(h.reads - reads, 1, 'replacement board lost incremental decoding');
  }
});

test('fallback rebuilds detached or adopted entries after missed body replacement records', async () => {
  for (const adopted of [false, true]) {
    const h = contentHarness();
    await h.clock.flush();
    const old = h.replaceBody(['cell hd_opened hd_type3', 'cell hd_closed'], { observed: false });
    if (adopted) for (const cell of old) {
      cell.isConnected = true;
      Object.defineProperty(cell, 'ownerDocument', { value: {} });
    }
    h.clock.fireInterval(1000);
    await h.clock.flush();
    assert.deepEqual(h.messages.at(-1).updates, [{ x: 0, y: 0, s: 13 }]);
  }
});

test('body replacement during delivery clears removed boards and recovers', async () => {
  const h = contentHarness();
  h.defer = true;
  await h.clock.flush();
  h.replaceBody([]);
  await h.clock.flush();
  assert.equal(h.messages.length, 1, 'replacement bypassed the delivery bound');
  h.responses.shift()({ ok: true });
  await h.clock.flush();
  assert.equal(h.messages.at(-1).w, 0);
  assert.deepEqual(h.messages.at(-1).cells, []);
  h.responses.shift()({ ok: true });
  h.defer = false;
  h.replaceBody(['cell hd_opened hd_type1', 'cell hd_closed', 'cell hd_closed']);
  await h.clock.flush();
  assert.equal(h.messages.at(-1).type, 'full');
  assert.deepEqual(h.messages.at(-1).cells, [11, 0, 0]);
});

test('capture coalesces bursts within a frame and caps sustained event traffic', async () => {
  const h = contentHarness(32);
  await h.clock.flush();
  const start = h.sentAt[0], reads = h.reads;
  await h.clock.advance(4);
  for (let i = 0; i < 500; ++i) h.change(i % 32, 'cell hd_opened hd_type1');
  assert.equal(h.clock.pending, 1, 'a burst must have one capture timer');
  await h.clock.advance(11);
  assert.equal(h.messages.length, 1, 'capture ran more than once per frame');
  await h.clock.advance(1);
  assert.equal(h.messages.length, 2, 'capture retained the old 50ms delay');
  assert.equal(h.sentAt[1] - start, 16);
  assert.equal(h.reads - reads, 32, 'burst repeatedly decoded the same cells');
  for (let i = 0; i < 128; ++i) {
    h.change(0, 'cell hd_opened hd_type' + i % 9);
    await h.clock.advance(1);
  }
  await h.clock.flush();
  for (let i = 1; i < h.sentAt.length; ++i)
    assert.ok(h.sentAt[i] - h.sentAt[i - 1] >= 16, 'continuous mutations exceeded the capture rate limit');
  assert.ok(h.messages.length <= 11, 'continuous mutations flooded messages');
});

test('pending delivery coalesces changes without consuming another capture interval', async () => {
  const h = contentHarness();
  h.defer = true;
  await h.clock.flush();
  await h.clock.advance(30);
  for (let i = 0; i < 100; ++i) h.change(0, 'cell hd_opened hd_type1');
  assert.equal(h.clock.pending, 0, 'an outstanding delivery must not start capture timers');
  await h.clock.advance(1);
  const acknowledgedAt = h.clock.now;
  h.responses.shift()({ ok: true });
  await h.clock.advance(0);
  assert.equal(h.messages.length, 2, 'pending capture did not follow acknowledgement');
  assert.equal(h.sentAt[1], acknowledgedAt, 'a blocked capture imposed an extra cooldown');
  assert.deepEqual(h.messages[1].updates, [{ x: 0, y: 0, s: 11 }]);
  h.responses.shift()({ ok: true });
  await h.clock.flush();
});

test('hiding cancels queued work and returning immediately resynchronizes after idle', async () => {
  const h = contentHarness();
  await h.clock.flush();
  await h.clock.advance(4);
  h.change(0, 'cell hd_opened hd_type1');
  assert.equal(h.clock.pending, 1);
  h.document.hidden = true;
  h.fire('visibilitychange');
  assert.equal(h.clock.pending, 0);
  const reads = h.reads;
  h.clock.fireInterval(1000);
  h.change(1, 'cell hd_closed hd_flag');
  assert.equal(h.clock.pending, 0, 'hidden work started a timer');
  await h.clock.advance(100);
  assert.equal(h.reads, reads);
  h.document.hidden = false;
  h.fire('visibilitychange');
  await h.clock.advance(0);
  assert.equal(h.messages.at(-1).type, 'full');
  assert.deepEqual(h.messages.at(-1).cells, [11, 2]);
});

test('resync during an outstanding delivery cannot acknowledge an obsolete baseline', async () => {
  const h = contentHarness();
  h.defer = true;
  await h.clock.flush();
  h.runtime.onMessage.fire({ type: 'force_full' }, {}, () => {});
  h.change(0, 'cell hd_opened hd_type1');
  await h.clock.flush();
  assert.equal(h.messages.length, 1);
  h.responses.shift()({ ok: true });
  await h.clock.flush();
  assert.equal(h.messages[1].type, 'full');
  assert.deepEqual(h.messages[1].cells, [11, 0]);
  h.responses.shift()({ ok: true });
  await h.clock.flush();
});

test('hidden tabs stop polling and resynchronize on return', async () => {
  const h = contentHarness();
  await h.clock.flush();
  h.document.hidden = true;
  const reads = h.reads;
  h.clock.fireInterval(1000); await h.clock.flush();
  assert.equal(h.reads, reads);
  h.document.hidden = false;
  h.fire('visibilitychange'); await h.clock.flush();
  assert.equal(h.messages.at(-1).type, 'full');
  assert.equal(h.messages.length, 2);
});

test('mine-count shortcut ignores editing fields and key repeat', async () => {
  const h = contentHarness();
  await h.clock.flush();
  const key = (value, extra = {}) => h.fire('keydown', {
    key: value, shiftKey: value === 'M', preventDefault() {}, stopPropagation() {}, ...extra
  });
  key('M', { target: { closest: () => ({}) } }); key('7'); key('Enter');
  await h.clock.flush();
  assert.equal(h.messages.length, 1);
  key('M'); key('M', { repeat: true }); key('2'); key('Enter');
  await h.clock.flush();
  assert.equal(h.messages.at(-1).mines_total, 2);
});

test('cell membership changes rebuild the cached board', async () => {
  const h = contentHarness();
  await h.clock.flush();
  h.change(0, 'hd_opened hd_type1');
  await h.clock.flush();
  assert.equal(h.messages.at(-1).w, 0, 'an incomplete grid must clear hints');
  h.change(0, 'cell hd_opened hd_type1');
  await h.clock.flush();
  assert.equal(h.messages.at(-1).type, 'full');
  assert.deepEqual(h.messages.at(-1).cells, [11, 0]);
});

test('a preexisting element becoming a cell expands the captured grid', async () => {
  const h = contentHarness(3);
  h.document.cells[2].value = 'hd_closed';
  await h.clock.flush();
  assert.equal(h.messages.at(-1).w, 2);
  h.change(2, 'cell hd_closed');
  await h.clock.flush();
  assert.equal(h.messages.at(-1).w, 3);
  assert.deepEqual(h.messages.at(-1).cells, [0, 0, 0]);
});

test('non-rendered or oversized geometry clears hints and recovers', async () => {
  const h = contentHarness();
  await h.clock.flush();
  h.document.cells[0].visible = false;
  h.fire('resize'); await h.clock.flush();
  assert.equal(h.messages.at(-1).w, 0);
  h.document.cells[0].visible = true;
  h.size = 20000;
  h.fire('resize'); await h.clock.flush();
  assert.equal(h.messages.at(-1).w, 0, 'native geometry limits must not be sent as a valid board');
  h.size = 20;
  h.fire('resize'); await h.clock.flush();
  assert.equal(h.messages.at(-1).w, 2);
});

test('large flood openings use compact snapshots and preserve the delta baseline', async () => {
  const h = contentHarness(128);
  await h.clock.flush();
  for (let i = 0; i < 64; ++i) h.change(i, 'cell hd_opened hd_type0');
  await h.clock.flush();
  const full = h.messages.at(-1);
  assert.equal(full.type, 'full');
  assert.equal(full.cells.filter((cell) => cell === 10).length, 64);
  const deltaBytes = JSON.stringify({ type: 'delta', updates: Array.from({ length: 64 }, (_, x) => ({ x, y: 0, s: 10 })) }).length;
  assert.ok(JSON.stringify(full).length < deltaBytes / 2, 'flood payload should be substantially smaller');
  h.change(100, 'cell hd_closed hd_flag');
  await h.clock.flush();
  assert.deepEqual(h.messages.at(-1).updates, [{ x: 100, y: 0, s: 2 }]);
});

function backgroundHarness() {
  const clock = timers(), requests = [], sockets = [], selections = [];
  let selected = [{ id: 11, url: 'https://minesweeper.online/game/1' }];
  let focused = true, deferSelection = false;
  class Socket {
    static OPEN = 1; static CONNECTING = 0;
    readyState = 0; bufferedAmount = 0; messages = [];
    constructor() { sockets.push(this); }
    send(text) { this.messages.push(JSON.parse(text)); }
    open() { this.readyState = 1; this.onopen(); }
    close() { this.readyState = 3; this.onclose(); }
  }
  const chrome = {
    runtime: { id: 'a'.repeat(32), lastError: null, onMessage: event() },
    storage: { local: { get: async () => ({ debug: false }), set: async () => {} }, onChanged: event() },
    commands: { onCommand: event() },
    tabs: {
      sendMessage: (id, message, callback) => { requests.push({ id, message }); callback(); },
      onActivated: event(), onRemoved: event(), onUpdated: event()
    },
    windows: { WINDOW_ID_NONE: -1, onFocusChanged: event(),
      getLastFocused: (_query, callback) => {
        const window = { id: 1, focused, tabs: selected.map((tab) => ({ active: true, status: 'complete', ...tab })) };
        if (deferSelection) selections.push(() => callback(window));
        else callback(window);
      }
    }
  };
  vm.runInNewContext(source('background.js'), { ...clock.api, console, URL, chrome, WebSocket: Socket }, { filename: 'background.js' });
  return { clock, sockets, requests, chrome, selections,
    set focused(value) { focused = value; },
    set deferSelection(value) { deferSelection = value; },
    select(tabs) { selected = tabs; chrome.tabs.onActivated.fire({ windowId: 1 }); },
    update(change) { selected = selected.map((tab) => ({ ...tab, ...change })); chrome.tabs.onUpdated.fire(selected[0].id, change); },
    send(message, overrides = {}) {
      let result;
      chrome.runtime.onMessage.fire(message, {
        id: chrome.runtime.id, frameId: 0, tab: { id: 11 }, url: 'https://minesweeper.online/game/1', ...overrides
      }, (response) => { result = response; });
      return result;
    }
  };
}

test('background requires full snapshots after reconnects and backpressure', async () => {
  const h = backgroundHarness();
  h.sockets[0].open();
  assert.equal(h.send({ type: 'delta', updates: [] }).ok, false);
  assert.equal(h.send({ type: 'full', w: 1, h: 1, cells: [0] }).ok, true);
  const first = h.sockets[0];
  first.close(); await h.clock.flush();
  h.sockets[1].open();
  assert.equal(h.send({ type: 'delta', updates: [] }).ok, false);
  assert.equal(h.send({ type: 'full', w: 1, h: 1, cells: [0] }).ok, true);
  h.sockets[1].bufferedAmount = 2 * 1024 * 1024;
  assert.equal(h.send({ type: 'delta', updates: [] }).ok, false);
  h.sockets[1].bufferedAmount = 0;
  assert.equal(h.send({ type: 'delta', updates: [] }).ok, false);
  assert.equal(h.send({ type: 'full', w: 1, h: 1, cells: [0] }).ok, true);
  first.onclose(); // A late old-socket event must not retire the new connection.
  assert.equal(h.send({ type: 'delta', updates: [] }).ok, true);
  h.clock.fireInterval(20000);
  assert.equal(h.sockets[1].messages.at(-1).type, 'ping');
});

test('only the active game tab and top frame can update the overlay', () => {
  const h = backgroundHarness();
  h.sockets[0].open();
  const full = { type: 'full', w: 1, h: 1, cells: [0] };
  assert.equal(h.send(full, { tab: { id: 12 } }).ok, false);
  assert.equal(h.send(full, { frameId: 1 }).ok, false);
  assert.equal(h.send(full, { url: 'https://minesweeper.online.attacker.test' }).ok, false);
  assert.equal(h.send(full).ok, true);
  h.select([{ id: 12, url: 'https://example.test' }]);
  assert.equal(h.sockets[0].messages.at(-1).w, 0);
  assert.equal(h.send(full).ok, false);
  h.select([{ id: 11, url: 'https://minesweeper.online/game/2' }]);
  assert.equal(h.send({ type: 'delta', updates: [] }).ok, false);
  assert.equal(h.send(full).ok, true);
  h.chrome.windows.onFocusChanged.fire(-1);
  assert.equal(h.sockets[0].messages.at(-1).w, 0);
});

test('unfocused windows cannot be reactivated by reconnects or tab updates', async () => {
  const h = backgroundHarness();
  h.sockets[0].open();
  const full = { type: 'full', w: 1, h: 1, cells: [0] };
  assert.equal(h.send(full).ok, true);
  h.focused = false;
  h.chrome.windows.onFocusChanged.fire(-1);
  h.chrome.tabs.onUpdated.fire(11, { status: 'complete' });
  assert.equal(h.send(full).ok, false);
  h.sockets[0].close(); await h.clock.flush();
  h.sockets[1].open();
  assert.equal(h.send(full).ok, false);
  h.focused = true;
  h.chrome.windows.onFocusChanged.fire(1);
  assert.equal(h.send(full).ok, true);
});

test('selection changes reject the old tab before asynchronous queries settle', () => {
  const h = backgroundHarness();
  h.sockets[0].open();
  const full = { type: 'full', w: 1, h: 1, cells: [0] };
  assert.equal(h.send(full).ok, true);
  h.deferSelection = true;
  h.select([{ id: 12, url: 'https://minesweeper.online/game/2' }]);
  assert.equal(h.send(full).ok, false, 'the previously selected tab is already obsolete');
  h.select([{ id: 13, url: 'https://minesweeper.online/game/3' }]);
  h.selections.pop()();
  assert.equal(h.send(full, { tab: { id: 13 } }).ok, true);
  h.selections.shift()();
  assert.equal(h.send(full, { tab: { id: 12 } }).ok, false, 'an older callback must not restore a previous selection');
  h.focused = false;
  h.chrome.windows.onFocusChanged.fire(-1);
  assert.equal(h.send(full, { tab: { id: 13 } }).ok, false);
});

test('navigating and inactive documents cannot restore an obsolete baseline', () => {
  const h = backgroundHarness();
  h.sockets[0].open();
  const full = { type: 'full', w: 1, h: 1, cells: [0] };
  assert.equal(h.send(full).ok, true);
  h.update({ status: 'loading' });
  assert.equal(h.send(full).ok, false);
  h.update({ status: 'complete', url: 'https://minesweeper.online/game/2' });
  assert.equal(h.send({ type: 'delta', updates: [] }).ok, false);
  for (const documentLifecycle of ['prerender', 'cached', 'pending_deletion'])
    assert.equal(h.send(full, { documentLifecycle }).ok, false);
  assert.equal(h.send(full, { documentLifecycle: 'active' }).ok, true);
});
