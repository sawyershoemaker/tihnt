let socket = null;
let debug = false;
let reconnectTimer = null;
let keepaliveTimer = null;
let activeTabId = null;
let hasSnapshot = false;
let selectionVersion = 0;
const ENDPOINT = 'ws://127.0.0.1:8765';
const MAX_BUFFERED_BYTES = 1024 * 1024;
const emptyBoard = () => ({ type: 'full', w: 0, h: 0, cells: [], mines_total: -1 });

function log(...args) { if (debug) console.log('[tihnt/background]', ...args); }
function isGameUrl(value) {
  try {
    const url = new URL(value);
    return (url.protocol === 'https:' || url.protocol === 'http:') &&
      (url.hostname === 'minesweeper.online' || url.hostname.endsWith('.minesweeper.online'));
  } catch { return false; }
}
function send(message) {
  if (!socket || socket.readyState !== WebSocket.OPEN) return false;
  if (socket.bufferedAmount > MAX_BUFFERED_BYTES) {
    // A delta dropped under backpressure invalidates the receiver's baseline.
    hasSnapshot = false;
    return false;
  }
  try { socket.send(JSON.stringify(message)); return true; }
  catch (error) { log('send failed', String(error)); return false; }
}
function requestFull() {
  if (activeTabId === null) return;
  chrome.tabs.sendMessage(activeTabId, { type: 'force_full' }, () => {
    const error = chrome.runtime.lastError;
    if (error) log('snapshot request', error.message);
  });
}
function selectTab(id) {
  if (id === activeTabId) return;
  activeTabId = id;
  hasSnapshot = false;
  if (!send(emptyBoard()) && socket?.readyState === WebSocket.OPEN) socket.close();
  requestFull();
}
function refreshActiveTab() {
  const version = ++selectionVersion;
  chrome.tabs.query({ active: true, lastFocusedWindow: true }, (tabs) => {
    const error = chrome.runtime.lastError;
    if (version !== selectionVersion) return;
    const tab = !error && tabs?.[0];
    selectTab(tab && isGameUrl(tab.url) ? tab.id : null);
  });
}
function reconnect() {
  if (reconnectTimer !== null) return;
  reconnectTimer = setTimeout(() => { reconnectTimer = null; connect(); }, 1000);
}
function connect() {
  if (socket && (socket.readyState === WebSocket.OPEN || socket.readyState === WebSocket.CONNECTING)) return;
  let current;
  try { current = new WebSocket(ENDPOINT); }
  catch (error) { log('connect failed', String(error)); reconnect(); return; }
  socket = current;
  current.onopen = () => {
    if (socket !== current) return;
    hasSnapshot = false;
    requestFull();
    refreshActiveTab();
    clearInterval(keepaliveTimer);
    // Chrome 116+: traffic within 30 seconds keeps the MV3 worker alive.
    keepaliveTimer = setInterval(() => {
      if (socket === current && current.readyState === WebSocket.OPEN) send({ type: 'ping' });
    }, 20000);
  };
  current.onclose = () => {
    if (socket !== current) return;
    socket = null;
    hasSnapshot = false;
    clearInterval(keepaliveTimer);
    keepaliveTimer = null;
    reconnect();
  };
  current.onerror = () => { if (socket === current) current.close(); };
  current.onmessage = (event) => log('received', String(event.data).slice(0, 200));
}

chrome.storage.local.get({ debug: false }).then((settings) => { debug = !!settings.debug; }).catch(() => {});
chrome.storage.onChanged.addListener((changes, area) => {
  if (area === 'local' && changes.debug) debug = !!changes.debug.newValue;
});
chrome.tabs.onActivated.addListener(refreshActiveTab);
chrome.tabs.onUpdated.addListener((id, change) => {
  if (id === activeTabId && (change.url || change.status === 'loading')) {
    hasSnapshot = false;
    if (!send(emptyBoard()) && socket?.readyState === WebSocket.OPEN) socket.close();
  }
  if (change.url || change.status === 'complete') refreshActiveTab();
});
chrome.tabs.onRemoved.addListener((id) => {
  if (id === activeTabId) selectTab(null);
  refreshActiveTab();
});
chrome.windows.onFocusChanged.addListener((id) => {
  if (id === chrome.windows.WINDOW_ID_NONE) { ++selectionVersion; selectTab(null); }
  else refreshActiveTab();
});
chrome.commands.onCommand.addListener((command) => {
  if (command === 'toggle-logging') {
    chrome.storage.local.get({ debug: false }).then(({ debug: current }) =>
      chrome.storage.local.set({ debug: !current })).catch(() => {});
  } else if (command === 'force-resend-board') {
    hasSnapshot = false;
    refreshActiveTab();
    requestFull();
  }
});
chrome.runtime.onMessage.addListener((message, sender, respond) => {
  if (sender.id !== chrome.runtime.id || sender.frameId !== 0 ||
      sender.tab?.id !== activeTabId || !isGameUrl(sender.url) ||
      (message?.type !== 'full' && message?.type !== 'delta')) {
    respond({ ok: false }); return;
  }
  if (!hasSnapshot && message.type !== 'full') {
    requestFull();
    respond({ ok: false }); return;
  }
  const ok = send(message);
  if (ok && message.type === 'full') hasSnapshot = true;
  if (!ok) { hasSnapshot = false; connect(); }
  respond({ ok });
});

refreshActiveTab();
connect();
