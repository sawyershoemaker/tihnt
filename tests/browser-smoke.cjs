// Optional end-to-end check: npm ci && npx playwright install chromium
// TIHNT_WS_FIXTURE must point to the CMake-built native WebSocket test server.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const net = require('node:net');
const { spawn } = require('node:child_process');
const { once } = require('node:events');
const { setTimeout: delay } = require('node:timers/promises');
const { chromium } = require('playwright');

async function waitFor(fn, label, timeout = 10000) {
  const deadline = Date.now() + timeout;
  while (Date.now() < deadline) { if (await fn()) return; await delay(50); }
  throw Error('Timed out: ' + label);
}
function board(width = 3, height = 2) {
  let html = '<!doctype html><title>TIHNT offline fixture</title><style>body{margin:40px}.board{display:grid;grid-template-columns:repeat(' + width +
    ',20px);width:max-content}.cell{box-sizing:border-box;width:20px;height:20px;border:1px solid #777}</style><div class="board">';
  for (let y = 0; y < height; ++y) for (let x = 0; x < width; ++x)
    html += '<div id="cell_' + x + '_' + y + '" class="cell hd_closed"></div>';
  return html + '</div>';
}
(async () => {
  assert.ok(process.env.TIHNT_WS_FIXTURE, 'Set TIHNT_WS_FIXTURE to the built ws_fixture executable');
  const reserve = net.createServer();
  reserve.listen(0, '127.0.0.1'); await once(reserve, 'listening');
  const port = reserve.address().port;
  await new Promise((resolve) => reserve.close(resolve));
  const root = path.resolve(process.env.TIHNT_TEST_WORK_DIR || path.join(__dirname, '..', 'build', 'browser'));
  fs.mkdirSync(root, { recursive: true });
  const run = fs.mkdtempSync(path.join(root, 'run-'));
  const extension = path.join(run, 'extension');
  fs.cpSync(path.join(__dirname, '..', 'extension'), extension, { recursive: true });
  const background = path.join(extension, 'background.js');
  fs.writeFileSync(background, fs.readFileSync(background, 'utf8').replace('127.0.0.1:8765', '127.0.0.1:' + port));
  const fixture = spawn(process.env.TIHNT_WS_FIXTURE, [String(port)], { windowsHide: true, stdio: ['pipe', 'pipe', 'pipe'] });
  fixture.stderr.on('data', () => {});
  let context;
  try {
    const [ready] = await once(fixture.stdout, 'data');
    assert.match(String(ready), /READY/);
    context = await chromium.launchPersistentContext(path.join(run, 'profile'), {
      channel: 'chromium', headless: true,
      args: ['--disable-extensions-except=' + extension, '--load-extension=' + extension]
    });
    await context.route('https://minesweeper.online/**', (route) => {
      const url = route.request().url();
      route.fulfill({ contentType: 'text/html', body: url.includes('/empty') ? '<!doctype html><title>No board</title>' :
        url.includes('/second') ? board(2, 1) : board() });
    });
    const worker = context.serviceWorkers()[0] || await context.waitForEvent('serviceworker');
    await worker.evaluate(() => {
      globalThis.sentForTest = [];
      const send = WebSocket.prototype.send;
      WebSocket.prototype.send = function (text) {
        globalThis.sentForTest.push(JSON.parse(text));
        return send.call(this, text);
      };
    });
    const sent = () => worker.evaluate(() => globalThis.sentForTest);
    const page = await context.newPage();
    const errors = [];
    page.on('pageerror', (error) => errors.push(error.message));
    await page.goto('https://minesweeper.online/first');
    await page.bringToFront();
    await waitFor(async () => (await sent()).some((m) => m.type === 'full' && m.w === 3), 'initial extension snapshot');
    await page.locator('#cell_0_0').evaluate((cell) => { cell.className = 'cell hd_opened hd_type1'; });
    await waitFor(async () => (await sent()).some((m) => m.type === 'delta' && m.updates.some((u) => u.x === 0 && u.y === 0 && u.s === 11)), 'DOM delta');

    const beforeReconnect = (await sent()).length;
    await worker.evaluate('socket.close()');
    await waitFor(async () => (await sent()).slice(beforeReconnect).some((m) => m.type === 'full' && m.cells[0] === 11), 'full resync after reconnect');

    const second = await context.newPage();
    await second.goto('https://minesweeper.online/second');
    await second.bringToFront();
    await waitFor(async () => (await sent()).at(-1)?.w === 2, 'active second board');
    const afterSwitch = (await sent()).length;
    await page.locator('#cell_1_0').evaluate((cell) => { cell.className = 'cell hd_closed hd_flag'; });
    await delay(1200);
    assert.equal((await sent()).slice(afterSwitch).filter((m) => m.type === 'delta' && m.updates?.length).length, 0, 'hidden board overwrote active tab');
    await second.close();
    await page.bringToFront();
    await waitFor(async () => (await sent()).slice(afterSwitch).some((m) => m.type === 'full' && m.w === 3 && m.cells[1] === 2), 'reactivated tab full snapshot');

    await page.locator('#cell_2_0').evaluate((cell) => { cell.className = 'cell hd_opened hd_type11'; });
    await waitFor(async () => (await sent()).at(-1)?.w === 0, 'terminal board clears hints');
    await page.goto('https://minesweeper.online/empty');
    await waitFor(async () => (await sent()).at(-1)?.w === 0, 'navigation clears hints');
    assert.deepEqual(errors, []);
    console.log('Chromium extension smoke passed: native handshake, capture, deltas, reconnect, tab isolation, game over, navigation.');
  } finally {
    if (context) await context.close();
    if (fixture.exitCode === null) {
      const exit = once(fixture, 'exit');
      fixture.stdin.end('stop\n');
      const timer = setTimeout(() => fixture.kill(), 3000);
      await exit; clearTimeout(timer);
    }
  }
})().catch((error) => { console.error(error); process.exitCode = 1; });
