# Emacs + CMake Workflow (Quick Setup)

This note gives a practical Emacs workflow for this repository so you can run
builds/tests with short commands (`b`, `t`) and one-key compile reruns.

## 1) Bash functions (recommended)

Put these in your `~/.bashrc` (or equivalent shell startup file):

```bash
# CMake helpers (repo-root aware, with PWD fallback)

function c() {  # configure
    local r
    r="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
    cmake \
        -S "$r" \
        -B "$r/build" \
        -G Ninja \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        "$@"
}

function b() {  # rebuild from existing build folder
    local r
    r="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
    cmake --build "$r/build" -j "$@"
}

function f() {  # fresh rebuild
    local r
    r="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
    cmake --build "$r/build" -j --clean-first "$@"
}

function t() {  # run tests using ctest
    local r
    r="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
    ctest --test-dir "$r/build" --output-on-failure "$@"
}
```

Then reload shell config:

```bash
source ~/.bashrc
```

If you previously used aliases `b/c/f/t`, remove them so functions are used:

```bash
unalias b c f t 2>/dev/null
source ~/.bashrc
type -a b
```

Expected: `b is a function`.

## 2) Emacs keybinding

In your Emacs init file (`~/.emacs`, `~/.emacs.d/init.el`, etc.), add:

```elisp
(global-set-key (kbd "<f5>") #'compile)
(setq compilation-scroll-output t)
```

Optional test hotkey:

```elisp
(global-set-key
 (kbd "<f6>")
 (lambda ()
   (interactive)
   (compile "t")))
```

## 3) Per-project compile command

At repo root, create `.dir-locals.el` with:

```elisp
((nil . ((compile-command . "b"))))
```

Now inside this repo:

- `M-x compile` (or `<f5>`) runs `b`
- in `*compilation*`, press `g` to rerun the same command quickly

## 4) If Emacs cannot find `b`

Some Emacs setups start shell commands without loading interactive bash
functions. If that happens, set `compile-command` to the full command instead:

```elisp
((nil . ((compile-command . "cmake --build build -j"))))
```

You can do the same for tests:

```elisp
(compile "ctest --test-dir build --output-on-failure")
```

## 5) Typical daily flow

From anywhere inside repo:

1. `c` once (or when CMake config changes)
2. `b` after edits
3. `t` to run deterministic testsuite (`ctest`)

From Emacs:

1. `<f5>` to build
2. `g` in compilation buffer to rebuild quickly
3. `<f6>` (if configured) to run tests
