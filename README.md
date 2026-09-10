<p align="center">
  <img src="https://repository-images.githubusercontent.com/1067421994/aa3a03f7-aaf1-4f5c-aab6-59676f1ec9ff" alt="TIHNT repository banner" width="640">
</p>

#

**TIHNT** or the *Tile Information & Hidden Node Toolkit* (holy backronym) is a project that serves as an overlay/cheat for https://minesweeper.online.

## 2-part system

- desktop overlay written in cpp with a small Windows UI layer and websocket server for game data
- browser extension (Manifest V3) that injects content scripts to get data for the overlay

## FEATURES!! woo

- on screen overlay with capture exclusion where Windows and the capture app support it
- discrete keyboard shortcuts for toggles
- advanced chording (toggleable), and a safe mode that blocks left clicks outside suggested cells
- bounded solver searches, background calculations, and incremental board updates

## keybinds

- **Shift + M** — Enter mine amount (helps with guesses); Enter or Shift + M saves, Escape cancels, and an empty entry resets to unknown
- **Ctrl + Alt + X** — Close overlay
- **Ctrl + P** — Toggle chording
- **Ctrl + Alt + Shift + W** — Toggle capture exclusion
- **Ctrl + Alt + S** — Toggle safety
- **Ctrl + Q** — Resend full board
- **Ctrl + Shift + L** — Toggle extension logging

The mine amount is remembered, so update or clear it when changing difficulty. Flags are assumed correct. Green is safe, red is a mine, and orange is a guess. Blue marks show a chord and mines to flag first; yellow marks show a chord ready to use. **Guesses can still hit mines, including with safety enabled.** Right clicks and scrolling pass through.

## compiling..

Windows, Visual Studio 2022 Build Tools with the C++ desktop workload, and CMake 3.20+:

```powershell
cmake -S . -B build -G "Visual Studio 17 2022" -A x64
cmake --build build --config Release --parallel
.\build\cpp\Release\tihnt.exe
```

Release enables compiler optimizations without requiring AVX2. The exe needs the Microsoft Visual C++ x64 runtime.

In Chrome/Chromium 116+, open the extensions page, enable Developer mode, click **Load unpacked**, and select `extension`. Open a game and run the exe. The extension connects to `127.0.0.1:8765` and reconnects automatically. Only the active game tab supplies the board; switching away or disconnecting clears its hints.

Capture exclusion uses [Windows display affinity](https://learn.microsoft.com/en-us/windows/win32/api/winuser/nf-winuser-setwindowdisplayaffinity), so support depends on Windows and the capture app. The overlay follows a visible Chromium rendering window and hides when another app is foreground.

## checks

Node.js 20+ enables the extension and websocket tests alongside the native solver, parser, and rendering tests:

```powershell
ctest --test-dir build -C Release --output-on-failure
.\build\cpp\Release\core_tests.exe --benchmark
```

For the real extension/native server browser check:

```powershell
npm ci --ignore-scripts
npx playwright install chromium
$env:TIHNT_WS_FIXTURE = (Resolve-Path .\build\cpp\Release\ws_fixture.exe).Path
npm run test:browser
```

GitHub Actions checks Windows and Linux, including the browser test and sanitizers. Linux builds the portable solver/parser tests; the overlay is Windows-only.
