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
    fireInterval(delay) { for (const value of [...intervals.values()]) if (value.delay === delay) value.fn(); },
    intervals
  };
}
function contentHarness(count = 2) {
  const clock = timers(), events = new Map(), observers = [], messages = [], responses = [];
  let defer = false, ok = true, reads = 0, cellSize = 20;
  const cells = Array.from({ length: count }, (_, x) => ({
    id: 'cell_' + x + '_0', nodeType: 1, value: 'cell hd_closed',
    get className() { ++reads; return this.value; },
    getBoundingClientRect: () => ({ left: x * cellSize, top: 50, right: (x + 1) * cellSize, bottom: 50 + cellSize, width: cellSize, height: cellSize })
  }));
  const add = (name, fn) => { if (!events.has(name)) events.set(name, []); events.get(name).push(fn); };
  const document = {
    hidden: false, body: {}, cells,
    querySelectorAll() { return this.cells; },
    getElementById(id) { return this.cells.find((c) => c.id === id); },
    addEventListener: add
  };
  const runtime = { lastError: null, onMessage: event(), sendMessage(message, respond) {
    messages.push(JSON.parse(JSON.stringify(message)));
    if (defer) responses.push(respond); else respond({ ok });
  }};
  const storage = { local: { get: async () => ({ debug: false, mines_total: -1 }), set: async () => {} }, onChanged: event() };
  const window = { devicePixelRatio: 1, visualViewport: { offsetLeft: 0, offsetTop: 0, scale: 1, addEventListener: add }, addEventListener: add };
  const context = { ...clock.api, console, document, window, chrome: { runtime, storage },
    MutationObserver: class { constructor(fn) { observers.push(fn); } observe() {} }
  };
  vm.runInNewContext(source('content.js'), context, { filename: 'content.js' });
  return {
    clock, document, window, runtime, messages, responses,
    set ok(value) { ok = value; }, set defer(value) { defer = value; }, set size(value) { cellSize = value; },
    get reads() { return reads; },
    fire(name, data = {}) { for (const fn of events.get(name) || []) fn(data); },
    change(index, value) {
      document.cells[index].value = value;
      observers[0]([{ type: 'attributes', attributeName: 'class', target: document.cells[index] }]);
    },
    removeBoard() {
      const old = document.cells; document.cells = [];
      observers[0]([{ type: 'childList', addedNodes: [], removedNodes: old }]);
    }
  };
}

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

function backgroundHarness() {
  const clock = timers(), requests = [], sockets = [];
  let selected = [{ id: 11, url: 'https://minesweeper.online/game/1' }];
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
      query: (_query, callback) => callback(selected),
      sendMessage: (id, message, callback) => { requests.push({ id, message }); callback(); },
      onActivated: event(), onRemoved: event(), onUpdated: event()
    },
    windows: { WINDOW_ID_NONE: -1, onFocusChanged: event() }
  };
  vm.runInNewContext(source('background.js'), { ...clock.api, console, URL, chrome, WebSocket: Socket }, { filename: 'background.js' });
  return { clock, sockets, requests, chrome,
    select(tabs) { selected = tabs; chrome.tabs.onActivated.fire(); },
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
