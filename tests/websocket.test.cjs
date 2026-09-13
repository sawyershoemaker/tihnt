const { test, before, after } = require('node:test');
const assert = require('node:assert/strict');
const net = require('node:net');
const { spawn } = require('node:child_process');
const { once } = require('node:events');
const { createInterface } = require('node:readline');
const { setTimeout: delay } = require('node:timers/promises');
let fixture, fixtureLines, port;
const sockets = new Set();

function outputLine(expected, timeout = 3000) {
  return new Promise((resolve, reject) => {
    const cleanup = () => {
      clearTimeout(timer);
      fixtureLines.off('line', onLine);
      fixture.off('exit', onExit);
      fixture.off('error', onError);
    };
    const onLine = (line) => { if (line === expected) { cleanup(); resolve(); } };
    const onExit = (code) => { cleanup(); reject(Error('fixture exited early: ' + code)); };
    const onError = (error) => { cleanup(); reject(error); };
    const timer = setTimeout(() => { cleanup(); reject(Error('fixture output timeout: ' + expected)); }, timeout);
    fixtureLines.on('line', onLine);
    fixture.once('exit', onExit);
    fixture.once('error', onError);
  });
}

async function start() {
  const reservation = net.createServer();
  reservation.listen(0, '127.0.0.1');
  await once(reservation, 'listening');
  port = reservation.address().port;
  await new Promise((resolve) => reservation.close(resolve));
  fixture = spawn(process.env.TIHNT_WS_FIXTURE, [String(port)], { windowsHide: true, stdio: ['pipe', 'pipe', 'pipe'] });
  fixtureLines = createInterface({ input: fixture.stdout });
  fixture.stderr.on('data', () => {});
  await outputLine('READY');
}
before(start);
after(async () => {
  for (const socket of sockets) socket.destroy();
  if (fixture && fixture.exitCode === null) {
    const exit = once(fixture, 'exit');
    fixture.stdin.end('stop\n');
    const timer = setTimeout(() => fixture.kill(), 3000);
    await exit;
    clearTimeout(timer);
  }
});

class Reader {
  buffer = Buffer.alloc(0); closed = false; wake = null;
  constructor(socket) {
    this.socket = socket;
    socket.on('data', (data) => { this.buffer = Buffer.concat([this.buffer, data]); this.wake?.(); });
    for (const event of ['close', 'end', 'error']) socket.on(event, () => { this.closed = true; this.wake?.(); });
  }
  async wait() {
    if (this.closed) throw Error('connection closed');
    await new Promise((resolve, reject) => {
      const timeout = setTimeout(() => { this.wake = null; reject(Error('read timeout')); }, 3000);
      this.wake = () => { clearTimeout(timeout); this.wake = null; resolve(); };
    });
  }
  async read(n) {
    while (this.buffer.length < n) await this.wait();
    const result = this.buffer.subarray(0, n);
    this.buffer = this.buffer.subarray(n);
    return result;
  }
  async headers() {
    while (!this.buffer.includes('\r\n\r\n')) await this.wait();
    const size = this.buffer.indexOf('\r\n\r\n') + 4;
    return String(await this.read(size));
  }
  async frame() {
    const header = await this.read(2);
    assert.equal(header[1] & 128, 0, 'server frames must be unmasked');
    let length = header[1] & 127;
    if (length === 126) length = (await this.read(2)).readUInt16BE();
    if (length === 127) length = Number((await this.read(8)).readBigUInt64BE());
    return { opcode: header[0] & 15, payload: await this.read(length) };
  }
}
function handshake(origin = 'chrome-extension://' + 'a'.repeat(32), extra = '') {
  return 'GET / HTTP/1.1\r\nHost: 127.0.0.1:' + port +
    '\r\nUpgrade: websocket\r\nConnection: keep-alive, Upgrade\r\nSec-WebSocket-Version: 13\r\n' +
    'Sec-WebSocket-Key: dGhlIHNhbXBsZSBub25jZQ==\r\nOrigin: ' + origin + '\r\n' + extra + '\r\n';
}
function frame(opcode, data, fin = true) {
  const payload = Buffer.isBuffer(data) ? data : Buffer.from(data);
  const header = Buffer.alloc(payload.length < 126 ? 2 : payload.length <= 65535 ? 4 : 10);
  header[0] = (fin ? 128 : 0) | opcode;
  if (header.length === 2) header[1] = 128 | payload.length;
  else if (header.length === 4) { header[1] = 128 | 126; header.writeUInt16BE(payload.length, 2); }
  else { header[1] = 128 | 127; header.writeBigUInt64BE(BigInt(payload.length), 2); }
  const mask = Buffer.from([17, 31, 47, 61]);
  const masked = Buffer.from(payload);
  for (let i = 0; i < masked.length; ++i) masked[i] ^= mask[i % 4];
  return Buffer.concat([header, mask, masked]);
}
async function tcp() {
  const socket = net.connect(port, '127.0.0.1');
  sockets.add(socket); socket.once('close', () => sockets.delete(socket));
  const reader = new Reader(socket);
  await once(socket, 'connect');
  return reader;
}
async function connect(extra = Buffer.alloc(0)) {
  const reader = await tcp();
  reader.socket.write(Buffer.concat([Buffer.from(handshake()), extra]));
  assert.match(await reader.headers(), /HTTP\/1.1 101/);
  return reader;
}
async function close(reader) {
  const event = reader.closed ? Promise.resolve() : once(reader.socket, 'close');
  reader.socket.destroy();
  await event;
}

test('handshake uses RFC accept key and preserves a coalesced first frame', async () => {
  const reader = await tcp();
  reader.socket.write(Buffer.concat([Buffer.from(handshake()), frame(1, 'first')]));
  assert.match(await reader.headers(), /Sec-WebSocket-Accept: s3pPLMBiTxaQ9kYGzzhZRbK\+xOo=/);
  assert.equal(String((await reader.frame()).payload), 'first');
  await close(reader);
});

test('TCP may split every header, mask, and extended-length byte', async () => {
  const reader = await connect();
  const payload = 'x'.repeat(160);
  const bytes = frame(1, payload);
  for (let i = 0; i < 8; ++i) { reader.socket.write(bytes.subarray(i, i + 1)); await delay(3); }
  reader.socket.write(bytes.subarray(8));
  const result = await reader.frame();
  assert.equal(result.opcode, 1);
  assert.equal(String(result.payload), payload);
  await close(reader);
});

test('fragmented text supports intervening ping and close handshakes', async () => {
  const reader = await connect();
  reader.socket.write(Buffer.concat([frame(1, 'hello ', false), frame(9, 'alive'), frame(0, 'world')]));
  const pong = await reader.frame();
  assert.equal(pong.opcode, 10); assert.equal(String(pong.payload), 'alive');
  assert.equal(String((await reader.frame()).payload), 'hello world');
  reader.socket.write(frame(8, Buffer.from([3, 232])));
  assert.equal((await reader.frame()).opcode, 8);
  await close(reader);
});

test('large valid text frames round-trip through partial sends', async () => {
  const reader = await connect();
  const text = 'hello \u{1f600}'.repeat(20000);
  reader.socket.write(frame(1, text));
  assert.equal(String((await reader.frame()).payload), text);
  await close(reader);
});

test('coalesced frames retain ordering across buffered reads and messages', async () => {
  const reader = await connect();
  const messages = Array.from({ length: 1000 }, (_, i) => 'message ' + i);
  reader.socket.write(Buffer.concat(messages.map((text) => frame(1, text))));
  for (const text of messages) {
    const result = await reader.frame();
    assert.equal(result.opcode, 1);
    assert.equal(String(result.payload), text);
  }
  await close(reader);
});

test('maximum-size messages and Unicode split across fragments remain valid', async () => {
  const reader = await connect();
  const text = 'a'.repeat(1024 * 1024 - 4) + '\u{1f600}';
  const bytes = Buffer.from(text);
  reader.socket.write(Buffer.concat([
    frame(1, bytes.subarray(0, bytes.length - 2), false),
    frame(9, 'boundary'), frame(0, bytes.subarray(bytes.length - 2))
  ]));
  assert.equal((await reader.frame()).opcode, 10);
  assert.equal(String((await reader.frame()).payload), text);
  await close(reader);
});

test('ordinary web origins and malformed upgrade headers are rejected', async () => {
  const requests = [
    handshake('https://minesweeper.online'), handshake('https://attacker.test'),
    handshake('chrome-extension://' + 'z'.repeat(32)),
    handshake().replace('Version: 13', 'Version: 12'),
    handshake().replace('Connection: keep-alive, Upgrade', 'Connection: close'),
    handshake().replace('dGhlIHNhbXBsZSBub25jZQ==', 'bad-key'),
    handshake().replace('127.0.0.1:' + port, 'attacker.test:' + port),
    handshake('chrome-extension://' + 'a'.repeat(32), 'Origin: https://attacker.test\r\n')
  ];
  for (const request of requests) {
    const reader = await tcp();
    const closed = once(reader.socket, 'close');
    reader.socket.write(request);
    await closed;
    assert.equal(reader.buffer.length, 0);
  }
});

test('protocol violations and allocation bombs close without delivering data', async () => {
  const malformed = [
    Buffer.from([0x81, 0]), // Unmasked client frame.
    Buffer.from([0xc1, 0x80]), // Reserved bit.
    frame(0, 'orphan'), frame(9, 'fragmented ping', false),
    frame(2, 'binary'), frame(1, Buffer.from([0xc0, 0xaf])),
    frame(8, Buffer.from([1])), frame(8, Buffer.from([3, 237])), // Reserved close code 1005.
    Buffer.from([0x81, 0xff, 0, 0, 0, 1, 0, 0, 0, 0]), // 4 GiB advertised payload.
    Buffer.from([0x81, 0xfe, 0, 1]) // Noncanonical short length.
  ];
  for (const bytes of malformed) {
    const reader = await connect();
    reader.socket.write(bytes);
    // A TCP reset can overtake the Close frame when rejected payload bytes
    // remain unread. Either outcome must close without delivering text.
    let response;
    try { response = await reader.frame(); }
    catch (error) { assert.ok(reader.closed, error.message); }
    if (response) {
      assert.equal(response.opcode, 8);
      assert.ok([1002, 1007, 1009].includes(response.payload.readUInt16BE()));
    }
    await close(reader);
  }
  const reader = await connect();
  reader.socket.write(Buffer.concat([frame(1, 'a', false), frame(1, 'b')]));
  let finalResponse;
  try { finalResponse = await reader.frame(); }
  catch (error) { assert.ok(reader.closed, error.message); }
  if (finalResponse) assert.equal(finalResponse.opcode, 8);
  await close(reader);
});

test('disconnects, replacements, and callback failures leave the server usable', async () => {
  for (let i = 0; i < 20; ++i) {
    const reader = await connect();
    reader.socket.write(frame(1, 'message ' + i));
    assert.equal(String((await reader.frame()).payload), 'message ' + i);
    await close(reader);
  }
  const previous = await connect();
  const replacement = await connect();
  replacement.socket.write(frame(1, 'new client'));
  assert.equal(String((await replacement.frame()).payload), 'new client');
  await close(previous);
  replacement.socket.write(frame(1, 'throw'));
  if (!replacement.closed) await once(replacement.socket, 'close');
  const healthy = await connect();
  healthy.socket.write(frame(1, 'recovered'));
  assert.equal(String((await healthy.frame()).payload), 'recovered');
  await close(healthy);
});

test('nonstandard message and connection callback exceptions preserve server lifecycle', async () => {
  const reader = await connect();
  reader.socket.write(frame(1, 'throw nonstandard'));
  if (!reader.closed) await once(reader.socket, 'close');
  const arm = await connect();
  arm.socket.write(frame(1, 'throw connections'));
  assert.equal(String((await arm.frame()).payload), 'throw connections');
  await close(arm);
  for (let i = 0; i < 3; ++i) {
    const healthy = await connect();
    healthy.socket.write(frame(1, 'connection callback recovered ' + i));
    assert.equal(String((await healthy.frame()).payload), 'connection callback recovered ' + i);
    await close(healthy);
  }
});

test('a replacement interrupts a blocked sender without delaying its first message', async () => {
  const previous = await connect();
  previous.socket.pause();
  const sending = outputLine('SENDING');
  const done = outputLine('FLOOD_DONE');
  previous.socket.write(frame(1, 'flood'));
  await sending;
  await delay(75);
  const start = performance.now();
  const replacement = await connect();
  replacement.socket.write(frame(1, 'replacement is ready'));
  assert.equal(String((await replacement.frame()).payload), 'replacement is ready');
  assert.ok(performance.now() - start < 750, 'replacement blocked behind a socket write');
  await done;
  await close(previous); await close(replacement);
});

test('partial headers and fragmented messages expire even while pings arrive', async () => {
  for (const fragmented of [false, true]) {
    const reader = await connect();
    const closed = new Promise((resolve) => reader.socket.once('close', resolve));
    const start = performance.now();
    reader.socket.write(fragmented ? frame(1, 'unfinished', false) : Buffer.from([0x81]));
    const keepalive = fragmented ? setInterval(() => {
      if (!reader.closed) reader.socket.write(frame(9, 'still waiting'));
    }, 100) : null;
    let forcedClose = false;
    const deadline = setTimeout(() => { forcedClose = true; reader.socket.destroy(); }, 6500);
    await closed;
    clearTimeout(deadline); clearInterval(keepalive);
    assert.ok(!forcedClose, 'partial frame ignored its receive deadline');
    assert.ok(performance.now() - start >= 4500, 'a valid partial frame was rejected before its deadline');
    while (reader.buffer.length) assert.equal((await reader.frame()).opcode, 10, 'incomplete text reached the callback');
  }
  const healthy = await connect();
  healthy.socket.write(frame(1, 'receive deadline recovered'));
  assert.equal(String((await healthy.frame()).payload), 'receive deadline recovered');
  await close(healthy);
});

test('an expired write closes the stream before it can send another frame', async () => {
  const reader = await connect();
  reader.socket.pause();
  const done = outputLine('FLOOD_DONE', 8000);
  reader.socket.write(frame(1, 'flood'));
  await done;
  // Discard the buffered portion of the intentionally interrupted large frame.
  reader.socket.removeAllListeners('data');
  const closed = reader.closed ? Promise.resolve() : once(reader.socket, 'close');
  reader.socket.resume();
  let forcedClose = false;
  const deadline = setTimeout(() => { forcedClose = true; reader.socket.destroy(); }, 2000);
  await closed;
  clearTimeout(deadline);
  assert.ok(!forcedClose, 'a failed send left a broken stream open');
  const healthy = await connect();
  healthy.socket.write(frame(1, 'write deadline recovered'));
  assert.equal(String((await healthy.frame()).payload), 'write deadline recovered');
  await close(healthy);
});

test('shutdown interrupts a blocked sender and an incomplete HTTP handshake', async () => {
  const idle = await connect();
  idle.socket.pause();
  const sending = outputLine('SENDING');
  idle.socket.write(frame(1, 'flood'));
  await sending;
  await delay(75);
  const slow = await tcp();
  slow.socket.write('GET / HTTP/1.1\r\n');
  await delay(25);
  const exit = once(fixture, 'exit');
  const start = performance.now();
  fixture.stdin.end('stop\n');
  const [code] = await exit;
  assert.equal(code, 0);
  assert.ok(performance.now() - start < 750, 'shutdown blocked behind a socket write');
  await close(idle); await close(slow);
});
