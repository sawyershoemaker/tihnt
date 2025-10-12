let socket = null;
let connected = false;
let debug = false;
const endpoints = ['ws://127.0.0.1:8765', 'ws://localhost:8765'];
let endpointIndex = 0;

function log(...args) {
  if (!debug) return;
  try {
    console.log('[mines-ext/bg]', ...args);
  } catch {}
}

chrome.storage.local.get({ debug: false }).then((res) => { debug = !!res.debug; });
chrome.storage.onChanged.addListener((changes, area) => {
  if (area === 'local' && Object.prototype.hasOwnProperty.call(changes, 'debug')) {
    debug = !!changes.debug.newValue;
    log('debug set to', debug);
  }
});

function connect() {
  if (socket && (socket.readyState === WebSocket.OPEN || socket.readyState === WebSocket.CONNECTING)) return;
  const url = endpoints[endpointIndex % endpoints.length];
  log('ws connect', url);
  try {
    socket = new WebSocket(url);
  } catch (e) {
    log('ws ctor error', url, String(e));
    endpointIndex = (endpointIndex + 1) % endpoints.length;
    setTimeout(connect, 1000);
    return;
  }
  socket.onopen = () => { connected = true; log('ws open', url); };
  socket.onclose = () => {
    connected = false;
    log('ws close', url);
    endpointIndex = (endpointIndex + 1) % endpoints.length;
    setTimeout(connect, 1000);
  };
  socket.onerror = (e) => {
    connected = false;
    log('ws error', url, e);
    try { socket.close(); } catch {}
  };
  socket.onmessage = (ev) => {
    log('ws message', ev.data?.slice?.(0, 200));
  };
}

connect();

chrome.commands.onCommand.addListener((cmd) => {
  if (cmd === 'toggle-logging') {
    chrome.storage.local.get({ debug: false }).then(({ debug: cur }) => {
      const next = !cur;
      chrome.storage.local.set({ debug: next });
      console.log('[mines-ext/bg]', 'toggled debug to', next);
    });
  }
  if (cmd === 'force-resend-board') {
    try {
      chrome.tabs.query({ active: true, currentWindow: true }, (tabs) => {
        const tab = tabs && tabs[0];
        if (!tab) { log('force-resend: no active tab'); return; }
        try {
          log('force-resend: requesting full snapshot on tab', tab.id);
          chrome.tabs.sendMessage(tab.id, { type: 'force_full' }, () => {});
        } catch (e) {
          log('force-resend error', String(e));
        }
      });
    } catch (e) {
      log('force-resend cmd error', String(e));
    }
  }
});

chrome.runtime.onMessage.addListener((msg, sender, respond) => {
  if (!connected || !socket || socket.readyState !== WebSocket.OPEN) {
    log('drop message (ws not open)', { readyState: socket?.readyState, connected, endpoint: endpoints[endpointIndex % endpoints.length] });
    respond({ ok: false });
    return;
  }
  try {
    const payload = JSON.stringify(msg);
    log('send', payload.slice(0, 200));
    socket.send(payload);
    respond({ ok: true });
  } catch (e) {
    log('send error', String(e));
    respond({ ok: false, error: String(e) });
  }
});


