// Optional end-to-end check: npm ci && npx playwright install chromium
// TIHNT_WS_FIXTURE must point to the CMake-built native WebSocket test server.
// TIHNT_CAPTURE_BENCHMARK=1 records optional latency and capture-time metrics.
// TIHNT_CAPTURE_SOURCE benchmarks an alternate script after common capture checks.
const assert = require('node:assert/strict');
const { createHash } = require('node:crypto');
const fs = require('node:fs');
const path = require('node:path');
const net = require('node:net');
const { spawn } = require('node:child_process');
const { once } = require('node:events');
const { setTimeout: delay } = require('node:timers/promises');
const { chromium } = require('playwright');

async function benchmarkCapture(page, worker) {
  const percentile = (values, p) => [...values].sort((a, b) => a - b)[Math.min(values.length - 1, Math.floor(values.length * p))];
  const stats = (values) => ({ samples: values.length, medianMs: percentile(values, .5), p95Ms: percentile(values, .95), maxMs: Math.max(...values) });
  const waitForValue = async (fn) => {
    let value;
    await waitFor(async () => { value = await fn(); return !!value; }, 'benchmark capture', 10000, 1);
    return value;
  };
  const metrics = await page.context().newCDPSession(page);
  await metrics.send('Performance.enable');
  await page.goto('https://minesweeper.online/benchmark');
  await page.bringToFront();
  await waitForValue(() => worker.evaluate(() => sentForTest.at(-1)?.w === 30));
  const latency = { idle: [], rapid: [] };
  let index = 0;
  for (const mode of ['idle', 'rapid']) {
    for (let sample = 0; sample < 45; ++sample) {
      await delay(mode === 'idle' ? 70 : 4);
      const at = await page.evaluate((i) => {
        const at = Date.now();
        document.getElementById('cell_' + i % 30 + '_' + Math.floor(i / 30)).className = 'cell hd_opened hd_type1';
        return at;
      }, index);
      const arrival = await waitForValue(() => worker.evaluate(({ index, at }) => timedForTest.find((event) =>
        event.at >= at && (event.message.type === 'full' ? event.message.cells[index] === 11 :
          event.message.updates?.some((update) => update.x === index % 30 && update.y === Math.floor(index / 30) && update.s === 11))), { index, at }));
      if (sample >= 5) latency[mode].push(arrival.at - at);
      ++index;
    }
  }
  const before = await worker.evaluate(() => timedForTest.length);
  await page.evaluate(() => { captureTimesForTest = []; });
  const beforeMetrics = Object.fromEntries((await metrics.send('Performance.getMetrics')).metrics.map(({ name, value }) => [name, value]));
  const burst = await page.evaluate(async () => {
    const start = Date.now();
    let events = 0, lastMutationAt = start;
    while (Date.now() - start < 500) {
      lastMutationAt = Date.now();
      document.querySelector('.board').style.transform = 'translateX(' + (++events) + 'px)';
      for (let i = 0; i < 480; ++i)
        document.getElementById('cell_' + i % 30 + '_' + Math.floor(i / 30)).className = 'cell hd_opened hd_type' + events % 9;
      await new Promise((resolve) => setTimeout(resolve, 1));
    }
    return { events, durationMs: Date.now() - start, lastMutationAt };
  });
  const final = await waitForValue(() => worker.evaluate(({ before, events }) => timedForTest.slice(before).find((event) =>
    event.message.rect_l === 40 + events), { before, events: burst.events }));
  const settledCount = await worker.evaluate(() => timedForTest.length);
  await delay(100);
  const sent = await worker.evaluate((before) => timedForTest.slice(before), before);
  const times = await page.evaluate(() => captureTimesForTest);
  const afterMetrics = Object.fromEntries((await metrics.send('Performance.getMetrics')).metrics.map(({ name, value }) => [name, value]));
  const extraMessagesAfterSettled = await worker.evaluate((settledCount) => timedForTest.length - settledCount, settledCount);
  await metrics.detach();
  return {
    board: { width: 30, height: 16 }, warmupSamplesPerMode: 5, idleGapMs: 70, rapidGapMs: 4,
    idle: stats(latency.idle), rapid: stats(latency.rapid),
    burst: {
      events: burst.events, durationMs: burst.durationMs, cellsMutatedPerEvent: 480,
      messages: sent.length, payloadBytes: sent.reduce((sum, event) => sum + JSON.stringify(event.message).length, 0),
      snapshot: { ...stats(times), totalMs: times.reduce((sum, value) => sum + value, 0) },
      rendererScriptMs: (afterMetrics.ScriptDuration - beforeMetrics.ScriptDuration) * 1000,
      rendererTaskMs: (afterMetrics.TaskDuration - beforeMetrics.TaskDuration) * 1000,
      finalMutationToSendMs: final.at - burst.lastMutationAt, extraMessagesAfterSettled
    }
  };
}

async function waitFor(fn, label, timeout = 10000, interval = 50) {
  const deadline = Date.now() + timeout;
  while (Date.now() < deadline) { if (await fn()) return; await delay(interval); }
  throw Error('Timed out: ' + label);
}
async function saveCaptureBenchmark(page, worker, captureSource, run) {
  const result = await benchmarkCapture(page, worker);
  result.browser = page.context().browser().version();
  result.captureIntervalMs = Number(/const CAPTURE_INTERVAL = (\d+);/.exec(captureSource)?.[1]);
  result.source = process.env.TIHNT_CAPTURE_SOURCE || path.join(__dirname, '..', 'extension', 'content.js');
  result.sourceSha256 = createHash('sha256').update(captureSource).digest('hex');
  const output = path.join(run, 'capture-benchmark.json');
  fs.writeFileSync(output, JSON.stringify(result, null, 2) + '\n');
  console.log('Capture benchmark: ' + output + '\n' + JSON.stringify(result));
}
function board(width = 3, height = 2) {
  let html = '<!doctype html><title>TIHNT offline fixture</title><style>body{margin:40px}.board{display:grid;grid-template-columns:repeat(' + width +
    ',20px);width:max-content}.cell{box-sizing:border-box;width:20px;height:20px;border:1px solid #777}' +
    '@keyframes testRotate{from{transform:rotate(0deg)}to{transform:rotate(360deg)}}' +
    '@keyframes testMove{from{transform:translateX(0px)}to{transform:translateX(20px)}}' +
    '@keyframes testFade{from{opacity:.8}to{opacity:1}}</style><div class="board">';
  for (let y = 0; y < height; ++y) for (let x = 0; x < width; ++x)
    html += '<div id="cell_' + x + '_' + y + '" class="cell hd_closed"></div>';
  return html + '</div>';
}
async function checkClipping(page, sent, worker) {
  await page.goto('https://minesweeper.online/clipping');
  await waitFor(async () => (await sent()).at(-1)?.w === 3, 'clipping fixture capture');
  await page.evaluate(() => {
    const board = document.querySelector('.board'), clip = document.createElement('div'), outer = document.createElement('div');
    clip.className = 'clip'; outer.className = 'outer';
    board.before(outer); outer.append(clip); clip.append(board);
  });
  const configure = (changes = {}) => page.evaluate((changes) => {
    document.querySelector('.outer').style.cssText = 'position:absolute;left:100px;top:100px;width:180px;height:120px;' + (changes.outer || '');
    document.querySelector('.clip').style.cssText = 'width:45px;height:25px;border:3px solid black;overflow:hidden;' + (changes.clip || '');
    document.querySelector('.board').style.cssText = 'position:relative;left:-8px;top:-6px;' + (changes.board || '');
  }, changes);
  const clipOf = (message) => [message.clip_l, message.clip_t, message.clip_w, message.clip_h];
  const close = (a, b) => a.length === b.length && a.every((value, i) => Math.abs(value - b[i]) < .001);
  const expect = async (expected, label, after = 0) => {
    await waitFor(async () => (await sent()).slice(after).some((m) => close(clipOf(m), expected)), label);
    return (await sent()).filter((m) => close(clipOf(m), expected)).at(-1);
  };
  const cases = [
    [{}, [103, 103, 45, 25], 'bordered overflow'],
    [{ clip: 'overflow:clip' }, [103, 103, 45, 25], 'overflow clip'],
    [{ clip: 'overflow:auto' }, [103, 103, 45, 25], 'overflow auto'],
    [{ outer: 'transform:scale(1.5,.75);transform-origin:top left' }, [104.5, 102.25, 67.5, 18.75], 'scaled overflow'],
    [{ outer: 'width:30px;height:20px;overflow:hidden', clip: 'margin-left:8px;margin-top:4px' }, [111, 107, 19, 13], 'nested overflow'],
    [{ board: 'left:80px;top:60px' }, [183, 163, 0, 0], 'fully clipped'],
    [{ board: 'position:fixed;left:90px;top:90px' }, [90, 90, 60, 40], 'fixed containing-block escape'],
    [{ board: 'position:absolute;left:80px;top:60px' }, [180, 160, 60, 40], 'absolute containing-block escape']
  ];
  for (const [changes, expected, label] of cases) {
    const before = (await sent()).length;
    await configure(changes);
    await worker.evaluate('requestFull()');
    const message = await expect(expected, label, before);
    assert.equal(message.cells?.length ?? 6, 6, label + ' preserves board data');
  }
  await configure();
  await expect([103, 103, 45, 25], 'base clip restores', (await sent()).length - 1);
  let before = (await sent()).length;
  await page.locator('.clip').evaluate((element) => { element.style.width = '30px'; });
  const shrunk = await expect([103, 103, 30, 25], 'parent-only resize sends presentation delta', before);
  assert.equal(shrunk.type, 'delta'); assert.deepEqual(shrunk.updates, []);
  before = (await sent()).length;
  await page.evaluate(() => {
    document.querySelector('.clip').style.transform = 'translateX(10px)';
    document.querySelector('.board').style.left = '-18px';
  });
  const shifted = await expect([113, 103, 30, 25], 'constant-ratio clip shift', before);
  assert.equal(shifted.rect_l, shrunk.rect_l); assert.equal(shifted.rect_w, shrunk.rect_w);

  before = (await sent()).length;
  await configure({ clip: 'overflow:visible' });
  await waitFor(async () => (await sent()).slice(before).some((m) => m.rect_w === 60 && m.clip_w === undefined), 'removing clip restores omitted full-board geometry');
  before = (await sent()).length;
  await configure({ clip: 'display:contents' });
  await expect([92, 94, 60, 40], 'display contents has no overflow box', before);

  before = (await sent()).length;
  await configure({ clip: 'overflow:visible', board: 'left:0;top:0;box-sizing:border-box;width:60px;height:40px;border:3px solid black;overflow:hidden' });
  await page.locator('.board').evaluate((element) => { element.style.gridTemplateColumns = 'repeat(3,20px)';
    for (const cell of element.children) cell.style.transform = 'translate(-3px,-3px)';
  });
  await expect([106, 106, 54, 34], 'parent own-overflow uses cell intersections', before);
  before = (await sent()).length;
  await page.locator('#cell_1_0').evaluate((element) => { element.style.clipPath = 'inset(2px)'; });
  await waitFor(async () => (await sent()).slice(before).some((m) => m.w === 0), 'nonrectangular cell visibility clears unsupported hints');

  before = (await sent()).length;
  await configure({ clip: 'width:30px;height:20px;border:0', board: 'left:0;top:0;width:60px;height:40px' });
  await page.locator('.board').evaluate((element) => {
    for (const cell of element.children) {
      const [, x, y] = cell.id.split('_').map(Number);
      cell.style.cssText = 'position:fixed;left:' + (100 + x * 20) + 'px;top:' + (100 + y * 20) + 'px';
    }
  });
  await expect([100, 100, 60, 40], 'fixed cells escape a bounds-matching parent clip', before);
  before = (await sent()).length;
  await page.locator('.board').evaluate((element) => {
    for (const cell of element.children) cell.style.cssText = '';
    const middle = element.children[1], wrapper = document.createElement('div');
    wrapper.className = 'inner-clip'; wrapper.style.cssText = 'width:20px;height:20px;overflow:hidden';
    middle.before(wrapper); wrapper.append(middle);
  });
  await expect([100, 100, 30, 20], 'wrapped cells retain contiguous visible intersections', before);
  before = (await sent()).length;
  await page.locator('.inner-clip').evaluate((element) => { element.style.width = '0'; });
  await expect([100, 100, 20, 20], 'interior ancestor style changes immediately refresh clipping', before);
  before = (await sent()).length;
  await page.locator('.clip').evaluate((element) => { element.style.width = '60px'; });
  await expect([100, 100, 0, 0], 'disconnected visible cells cannot fill a hole', before);

  await page.goto('https://minesweeper.online/fractional');
  await waitFor(async () => (await sent()).at(-1)?.w === 30, 'fractional fallback initial capture');
  before = (await sent()).length;
  await page.locator('.board').evaluate((element) => {
    element.style.cssText = 'width:600.3px;grid-template-columns:repeat(30,1fr);transform:scale(1.1);transform-origin:top left;overflow:hidden';
    for (const cell of element.children) cell.style.width = 'auto';
  });
  await waitFor(async () => (await sent()).slice(before).some((m) => m.clip_w > 600 && m.clip_h > 300), 'fractional cell intersections remain contiguous');
}
async function checkInteriorClipping(page, sent, worker) {
  await page.goto('https://minesweeper.online/interior-clipping');
  await waitFor(async () => (await sent()).at(-1)?.w === 3 && (await sent()).at(-1)?.h === 3, 'interior clipping fixture capture');
  await page.evaluate(() => {
    const board = document.querySelector('.board'), clip = document.createElement('div');
    clip.className = 'interior-clip';
    clip.style.cssText = 'position:absolute;left:100px;top:100px;overflow:hidden;width:60px;height:60px';
    board.before(clip); clip.append(board);
    board.style.cssText = 'width:60px;height:60px;grid-template-rows:repeat(3,20px)';
    for (const cell of board.children) {
      const [, x, y] = cell.id.split('_').map(Number);
      cell.style.cssText = 'grid-column:' + (x + 1) + ';grid-row:' + (y + 1);
    }
    document.getElementById('cell_1_0').className = 'cell hd_closed hd_flag';
    document.getElementById('cell_1_1').className = 'cell hd_opened hd_type1';
  });
  const fullSnapshot = async (label) => {
    const before = (await sent()).length;
    await worker.evaluate('requestFull()');
    await waitFor(async () => (await sent()).slice(before).some((message) => message.type === 'full'), label);
    return (await sent()).slice(before).filter((message) => message.type === 'full').at(-1);
  };
  const supported = (message, label) => {
    assert.equal(message.w, 3, label); assert.equal(message.h, 3, label);
    assert.deepEqual(message.cells, [0, 2, 0, 0, 11, 0, 0, 0, 0], label + ' preserves cell-to-index mapping');
    const clip = [message.clip_l, message.clip_t, message.clip_w, message.clip_h];
    assert.ok(clip.every((value, i) => Math.abs(value - [100, 100, 60, 60][i]) < .001), label + ' preserves the full visible board');
  };
  supported(await fullSnapshot('regular interior grid is captured'), 'regular interior grid');
  const cases = [
    ['relative center movement', 'position:relative;left:20px', 'moved'],
    ['non-rendered center', 'display:none', 'hidden'],
    ['hidden center', 'visibility:hidden', 'hidden'],
    ['transparent center', 'opacity:0', 'hidden'],
    ['rounded center overflow', 'overflow:hidden;border-radius:50%', 'rounded'],
    ['permuted interior cells', '', 'swapped']
  ];
  for (const [label, style, kind] of cases) {
    await page.locator('#cell_1_1').evaluate((cell, { style, kind }) => {
      cell.style.cssText += ';' + style;
      if (kind === 'swapped') {
        cell.style.gridRow = '1';
        document.getElementById('cell_1_0').style.gridRow = '2';
      }
    }, { style, kind });
    const fixture = await page.locator('#cell_1_1').evaluate((cell) => ({
      visible: cell.checkVisibility({ checkOpacity: true, checkVisibilityCSS: true }),
      left: cell.getBoundingClientRect().left, top: cell.getBoundingClientRect().top,
      cornerHit: document.elementFromPoint(121, 121)?.id || '',
      adjacentTop: document.getElementById('cell_1_0').getBoundingClientRect().top
    }));
    if (kind === 'hidden') assert.equal(fixture.visible, false, label + ' fixture must hide its center');
    if (kind === 'moved') assert.equal(fixture.left, 140, label + ' fixture must leave a center hole');
    if (kind === 'rounded') assert.notEqual(fixture.cornerHit, 'cell_1_1', label + ' fixture must clip the center corner');
    if (kind === 'swapped') assert.deepEqual([fixture.top, fixture.adjacentTop], [100, 120], label + ' fixture must exchange non-corner cells');
    // Request a fresh stable snapshot so a transient empty clip during mutation cannot satisfy this check.
    const invalid = await fullSnapshot(label + ' fresh capture');
    assert.ok(invalid.w === 0 || invalid.h === 0 || invalid.clip_w === 0 || invalid.clip_h === 0,
      label + ' retained visible hints for unsupported interior geometry');
    await page.locator('.board').evaluate((board) => {
      for (const cell of board.children) {
        const [, x, y] = cell.id.split('_').map(Number);
        cell.style.cssText = 'grid-column:' + (x + 1) + ';grid-row:' + (y + 1);
      }
    });
    supported(await fullSnapshot(label + ' restoration capture'), label + ' restoration');
  }
  await page.evaluate(() => {
    document.querySelector('.interior-clip').style.visibility = 'hidden';
    for (const cell of document.querySelectorAll('.cell')) cell.style.visibility = 'visible';
  });
  assert.equal(await page.locator('#cell_1_1').evaluate((cell) => cell.checkVisibility({ checkOpacity: true, checkVisibilityCSS: true })), true,
    'visible descendants of a hidden-visibility ancestor remain rendered');
  supported(await fullSnapshot('descendant visibility override capture'), 'descendant visibility override');
  await page.evaluate(() => {
    document.querySelector('.interior-clip').style.visibility = '';
    for (const cell of document.querySelectorAll('.cell')) cell.style.visibility = '';
  });
  supported(await fullSnapshot('ancestor visibility restoration capture'), 'ancestor visibility restoration');

  await page.evaluate(() => {
    const style = document.createElement('style');
    style.textContent = '#cell_1_1:not(:empty):not(:has(span)),#cell_1_1:has(.clip-trigger),' +
      '#cell_1_1:has(span:not(:empty)),#cell_1_1[data-cut]{clip:rect(0px,10px,20px,0px)}';
    document.head.append(style);
    document.querySelector('.board').style.position = 'relative';
    for (const cell of document.querySelectorAll('.cell')) {
      const [, x, y] = cell.id.split('_').map(Number);
      cell.style.cssText = 'position:absolute;left:' + x * 20 + 'px;top:' + y * 20 + 'px';
    }
  });
  for (const kind of ['text', 'characterData', 'nestedClass', 'nestedChild', 'attribute']) {
    await page.locator('#cell_1_1').evaluate((cell, kind) => {
      cell.replaceChildren(); cell.removeAttribute('data-cut');
      if (kind === 'characterData') cell.append(document.createTextNode(''));
      if (kind.startsWith('nested')) cell.append(document.createElement('span'));
    }, kind);
    supported(await fullSnapshot(kind + ' initial content capture'), kind + ' initial content');
    const beforeContent = (await sent()).length;
    const rect = await page.locator('#cell_1_1').evaluate((cell, kind) => {
      const before = cell.getBoundingClientRect().toJSON();
      if (kind === 'text') cell.append(document.createTextNode('1'));
      else if (kind === 'characterData') cell.firstChild.data = '1';
      else if (kind === 'nestedClass') cell.firstChild.className = 'clip-trigger';
      else if (kind === 'nestedChild') cell.firstChild.append(document.createTextNode('1'));
      else cell.setAttribute('data-cut', '');
      return { before, after: cell.getBoundingClientRect().toJSON(), clip: getComputedStyle(cell).clip };
    }, kind);
    assert.deepEqual(rect.after, rect.before, kind + ' must preserve outer cell bounds');
    assert.notEqual(rect.clip, 'auto', kind + ' must change native clipping');
    await waitFor(async () => (await sent()).slice(beforeContent).some((message) => message.clip_w === 0), kind + ' content refresh');
    const invalid = await fullSnapshot(kind + ' stable content capture');
    assert.equal(invalid.clip_w, 0, kind + ' interior hole retained visible hints');
    await page.locator('#cell_1_1').evaluate((cell) => { cell.replaceChildren(); cell.removeAttribute('data-cut'); });
    supported(await fullSnapshot(kind + ' content restoration'), kind + ' content restoration');
  }
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
  const benchmark = process.env.TIHNT_CAPTURE_BENCHMARK === '1';
  let captureSource;
  if (benchmark) {
    const content = path.join(extension, 'content.js');
    captureSource = fs.readFileSync(process.env.TIHNT_CAPTURE_SOURCE || content, 'utf8');
    const instrument = `
  const snapshotForBenchmark = snapshot;
  snapshot = (...args) => {
    const start = performance.now();
    try { return snapshotForBenchmark(...args); }
    finally { window.postMessage({ tihntSnapshotMsForTest: performance.now() - start }, '*'); }
  };
`;
    assert.ok(captureSource.includes('  chrome.storage.local.get('), 'Cannot instrument snapshot capture');
    fs.writeFileSync(content, captureSource.replace('  chrome.storage.local.get(', instrument + '  chrome.storage.local.get('));
  }
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
    if (benchmark) await context.addInitScript(() => {
      window.captureTimesForTest = [];
      window.addEventListener('message', (event) => {
        if (event.source === window && Number.isFinite(event.data?.tihntSnapshotMsForTest)) captureTimesForTest.push(event.data.tihntSnapshotMsForTest);
      });
    });
    await context.route('https://minesweeper.online/**', (route) => {
      const url = route.request().url();
      route.fulfill({ contentType: 'text/html', body: url.includes('/empty') ? '<!doctype html><title>No board</title>' :
        url.includes('/second') ? board(2, 1) : url.includes('/interior-clipping') ? board(3, 3) :
          url.includes('/benchmark') || url.includes('/fractional') ? board(30, 16) : board() });
    });
    const worker = context.serviceWorkers()[0] || await context.waitForEvent('serviceworker');
    await worker.evaluate((benchmark) => {
      globalThis.sentForTest = [];
      globalThis.timedForTest = [];
      const send = WebSocket.prototype.send;
      WebSocket.prototype.send = function (text) {
        const message = JSON.parse(text);
        globalThis.sentForTest.push(message);
        if (benchmark) globalThis.timedForTest.push({ at: Date.now(), message });
        return send.call(this, text);
      };
    }, benchmark);
    const sent = () => worker.evaluate(() => globalThis.sentForTest);
    const page = await context.newPage();
    const errors = [];
    page.on('pageerror', (error) => errors.push(error.message));
    await page.goto('https://minesweeper.online/first');
    await page.bringToFront();
    await waitFor(async () => (await sent()).some((m) => m.type === 'full' && m.w === 3), 'initial extension snapshot');
    await page.locator('#cell_0_0').evaluate((cell) => { cell.className = 'cell hd_opened hd_type1'; });
    await waitFor(async () => (await sent()).some((m) => m.type === 'delta' && m.updates.some((u) => u.x === 0 && u.y === 0 && u.s === 11)), 'DOM delta');

    if (benchmark && process.env.TIHNT_CAPTURE_SOURCE) {
      await saveCaptureBenchmark(page, worker, captureSource, run);
      assert.deepEqual(errors, []);
      console.log('Alternate-source capture benchmark passed: native handshake, capture, deltas, benchmark. Current-source smoke regressions were skipped for TIHNT_CAPTURE_SOURCE.');
      return;
    }

    for (const hiddenStyle of ['visibility:hidden', 'opacity:0', 'content-visibility:hidden', 'display:none']) {
      const start = (await sent()).length;
      await page.locator('.board').evaluate((element, style) => { element.style.cssText = style; }, hiddenStyle);
      await waitFor(async () => (await sent()).slice(start).some((m) => m.type === 'full' && m.w === 0), 'hidden board clears: ' + hiddenStyle);
      await page.locator('.board').evaluate((element) => { element.style.cssText = ''; });
      await waitFor(async () => (await sent()).slice(start).some((m) => m.type === 'full' && m.w === 3), 'visible board restores: ' + hiddenStyle);
    }

    const beforeMembership = (await sent()).length;
    await page.locator('#cell_1_0').evaluate((cell) => { cell.className = 'hd_closed'; });
    await waitFor(async () => (await sent()).slice(beforeMembership).some((m) => m.type === 'full' && m.w === 0), 'incomplete DOM grid clears hints');
    await page.locator('#cell_1_0').evaluate((cell) => { cell.className = 'cell hd_closed'; });
    await waitFor(async () => (await sent()).slice(beforeMembership).some((m) => m.type === 'full' && m.w === 3 && m.cells[0] === 11), 'cell membership restores snapshot');

    const beforeResize = (await sent()).length;
    await page.locator('.board').evaluate((element) => { element.style.cssText = 'transform:scale(.5);transform-origin:top left'; });
    await waitFor(async () => (await sent()).slice(beforeResize).some((m) => m.rect_w === 30 && m.rect_h === 20 && m.cell_px === 10), 'CSS scaling updates geometry');
    await page.locator('.board').evaluate((element) => { element.style.cssText = ''; });
    await waitFor(async () => (await sent()).slice(beforeResize).some((m) => m.rect_w === 60 && m.rect_h === 40), 'geometry restores');

    const geometryChange = async (action, predicate, label) => {
      const start = (await sent()).length;
      await action();
      await waitFor(async () => (await sent()).slice(start).some(predicate), label);
    };
    const boardStyle = (style) => page.locator('.board').evaluate((element, value) => { element.style.cssText = value; }, style);
    const restoreGeometry = () => geometryChange(() => boardStyle(''),
      (m) => m.rect_l === 40 && m.rect_t === 40 && m.rect_w === 60 && m.rect_h === 40, 'supported geometry restores');
    for (const style of ['transform:rotate(180deg)', 'transform:rotate(90deg)', 'transform:skewX(15deg)',
      'transform:scaleX(-1)', 'rotate:0.5turn', 'rotate:0 0 1 180deg', 'scale:-1 1', 'direction:rtl']) {
      await geometryChange(() => boardStyle(style), (m) => m.type === 'full' && m.w === 0, 'unsupported geometry clears: ' + style);
      await restoreGeometry();
    }
    await geometryChange(() => boardStyle('transform:translate(7.5px,-3px) scale(1.5,.75);transform-origin:top left'),
      (m) => m.rect_l === 47.5 && m.rect_t === 37 && m.rect_w === 90 && m.rect_h === 30, 'positive axis-aligned transform');
    await restoreGeometry();
    await geometryChange(() => boardStyle('animation:testRotate 10s 60s infinite linear'),
      (m) => m.type === 'full' && m.w === 0, 'delayed rotation hides hints at identity');
    assert.equal(await page.locator('.board').evaluate((element) => getComputedStyle(element).transform), 'none');
    await restoreGeometry();
    await geometryChange(() => boardStyle('animation:testRotate 10s infinite linear paused'),
      (m) => m.type === 'full' && m.w === 0, 'paused rotation hides hints at identity');
    await restoreGeometry();
    await geometryChange(() => boardStyle('animation:testMove 100s linear forwards'),
      (m) => m.type === 'full' && m.w === 0, 'moving board hides hints');
    await geometryChange(() => page.locator('.board').evaluate((element) => { element.getAnimations()[0].finish(); }),
      (m) => m.w === 3 && m.rect_l === 60, 'finished animation restores stable geometry');
    await restoreGeometry();
    await geometryChange(() => boardStyle('animation:testRotate 100s linear'),
      (m) => m.type === 'full' && m.w === 0, 'rotation begins');
    await geometryChange(() => page.locator('.board').evaluate((element) => { element.getAnimations()[0].cancel(); }),
      (m) => m.w === 3 && m.rect_l === 40, 'cancelled animation restores geometry');
    await boardStyle('animation:testFade 100s infinite linear');
    const beforeFade = (await sent()).length;
    await page.locator('#cell_2_1').evaluate((cell) => { cell.className = 'cell hd_closed hd_flag'; });
    await waitFor(async () => (await sent()).slice(beforeFade).some((m) => m.type === 'delta' && m.updates.some((u) => u.x === 2 && u.y === 1 && u.s === 2)), 'paint-only animation preserves capture');
    await boardStyle('');

    for (const root of ['body', 'documentElement']) {
      const beforeReplacement = (await sent()).length;
      await page.evaluate((root) => {
        const old = document[root], replacement = old.cloneNode(true);
        replacement.querySelector('#cell_1_0').className = 'cell hd_opened hd_type3';
        old.replaceWith(replacement);
      }, root);
      await waitFor(async () => (await sent()).slice(beforeReplacement).some((m) => m.type === 'delta' && m.updates.some((u) => u.x === 1 && u.y === 0 && u.s === 13)), root + ' replacement updates capture');
      const afterReplacement = (await sent()).length;
      await page.locator('#cell_1_0').evaluate((cell) => { cell.className = 'cell hd_closed'; });
      await waitFor(async () => (await sent()).slice(afterReplacement).some((m) => m.type === 'delta' && m.updates.some((u) => u.x === 1 && u.y === 0 && u.s === 0)), root + ' replacement remains observed');
    }

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
    await page.goto('https://minesweeper.online/fractional');
    await waitFor(async () => (await sent()).at(-1)?.w === 30, 'fractional board initial capture');
    for (const scale of [8, 16]) {
      const beforeFractional = (await sent()).length;
      await page.locator('.board').evaluate((element, scale) => {
        element.style.cssText = 'width:503.9px;grid-template-columns:repeat(30,1fr);transform:scale(' + scale + ');transform-origin:top left';
        for (const cell of element.children) cell.style.width = 'auto';
      }, scale);
      const expected = await page.locator('.board').evaluate((element) => element.getBoundingClientRect().width);
      await waitFor(async () => (await sent()).slice(beforeFractional).some((m) => Math.abs(m.rect_w - expected) < .01), 'fractional grid remains supported at scale ' + scale);
    }
    await checkClipping(page, sent, worker);
    await checkInteriorClipping(page, sent, worker);
    assert.deepEqual(errors, []);
    console.log('Chromium extension smoke passed: native handshake, capture, deltas, CSS visibility, membership, body/root replacement, geometry rejection/restoration, animations, fractional scaling, overflow and interior clipping, reconnect, tab isolation, game over, navigation.');
    if (benchmark) await saveCaptureBenchmark(page, worker, captureSource, run);
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
