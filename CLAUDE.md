# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Harp practice, hymnal, and music analysis tools.

## Project Structure

```
harp/
  harp-hymnal/        Rust eframe/egui native hymnal app
  harp-trainer/       Rust eframe/egui native drill app (chord ID, sight-reading)
  harp-pitch-wasm/    Rust WASM polyphonic pitch detector (FFT, pedal-constrained)
  harp-voicing-wasm/  Rust WASM voicing engine
  harp-midi/          Rust MIDI generator for harp voicings
  html/               Browser apps (harp_trainer.html, harp_voicing.html, HarpHymnal.html)
  scripts/            Python tools (analysis, PDF generation, voicing pipeline)
  data/               ABC source files (OpenHymnal.abc, harp_leadsheets.abc), voicings JSON
  docs/               Specs (SRS), notes, reference PDFs, desktop launcher
```

## Build Commands

```bash
# Harp Hymnal (native app — double-click or cargo run)
cd harp-hymnal && cargo run

# Harp Trainer (native app)
cd harp-trainer && cargo run

# WASM pitch detector
cd harp-pitch-wasm && wasm-pack build --target web

# HTML apps — serve from repo root
python3 -m http.server 8000
# then open localhost:8000/html/HarpHymnal.html

# PDF drill sheet
python3 scripts/drills.py
```

## Key Architecture

- **harp-hymnal**: Parses SATB voices from `data/OpenHymnal.abc`, expands 4-note harmony to 6-8 note harp voicings (R1-R10 scoring rules), displays leadsheet with string numbers per finger. Modules: `music.rs`, `chord.rs`, `voicing.rs`, `abc.rs`.
- **harp_trainer.html**: Self-contained browser app. VexFlow notation, pedal state drives key/pitches, optional WASM mic input. Three drill modes.
- **HarpHymnal.html**: Browser port of the hymnal — same voicing engine in JS. Auto-loads `../data/OpenHymnal.abc`.
- **Voicing engine**: Melody on RH thumb, LH carries the load (3-4 notes), 10-string max hand span, voice leading optimization, chromatic-to-diatonic snapping.
- **Chord notation**: 7-bit ASCII Roman numeral system (e.g., `V7`, `ii%7`, `bVII+`). Supports inversions (`^1`, `^2`), suspensions (`~2`, `~4`).
- **`data/OpenHymnal.abc`**: 284 SATB hymns transposed to C major. Original is `OpenHymnal2014.06.abc`.
