# TIHNT

TIHNT or the Tile Information & Hidden Node Toolkit (holy backronym) is a project that serves as an overlay/cheat for https://minesweeper.online. 

## 2-part system
- desktop overlay written in cpp with a small Windows UI layer and websocket server for game data
- browser extension (Manifest V3) that injects content scripts to get data for the overlay

## FEATURES!! woo
- on screen overlay (STREAMPROOOF, literally invisible on obs/discord/whatever respects the flag)
- discrete keyboard shortcuts for toggles
- advanced chording (toggleable), and a safe mode (prevents you from clicking mines)
- optimized asf! minimal cpu overhead compared to calculations and such an advanced solver

## keybinds
- **Shift + M** — Enter mine amount (helps with guesses)
- **Ctrl + Alt + X** — Close overlay
- **Ctrl + P** — Toggle chording
- **Ctrl + W**— Toggle stream-proof mode
- **Ctrl + Alt + S** — Toggle safety
- **Ctrl + Q**— Resend full board

## compiling..
```cmake -S . -B build -G "Ninja" -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_STANDARD_REQUIRED=ON -DCMAKE_CXX_FLAGS_RELEASE="/O2 /Ot /Ob3 /DNDEBUG /arch:AVX2 /fp:precise"``` 

and make sure you're compiling for x64, also if u cant use AVX2, use `/arch:SSE2`

---

>

