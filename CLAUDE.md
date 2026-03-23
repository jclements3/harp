# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

Harp practice and music analysis tools. Three main components:

1. **`harp_trainer.html`** — Single-file browser app (2300+ lines) for harp chord/sight-reading drills. Uses VexFlow 4.2.3 for notation, optional WASM polyphonic pitch detection. Three drill modes: Chord ID, Interval→Chord, and Sight Read (with shape notes, metronome, mic input). Pedal diagram UI controls key/tuning. Mobile-optimized (PWA-capable).

2. **`harp-pitch-wasm/`** — Rust→WASM polyphonic pitch detector. FFT-based, constrained to harp pedal state (7 available pitch classes). Built with `wasm-pack build --target web`. Used by the HTML trainer's mic input mode.

3. **`harp-trainer/`** — Native desktop port of the trainer using `eframe`/`egui` + `cpal` for audio. Mirrors the HTML version's drill logic. Modules: `audio.rs`, `drill.rs`, `music.rs`, `staff.rs`, `metronome.rs`.

## Build Commands

```bash
# WASM pitch detector (requires wasm-pack)
cd harp-pitch-wasm && wasm-pack build --target web

# Native desktop trainer
cd harp-trainer && cargo run

# HTML trainer — just open harp_trainer.html in a browser (or serve locally for WASM)
# python3 -m http.server   # from repo root, then open localhost:8000/harp_trainer.html

# PDF drill sheet
python3 drills.py   # outputs harp_practice_drill.pdf (requires reportlab)
```

## Music Analysis Scripts (Python)

- **`analyze_hymnal.py`** — Derives chord progressions from SATB voices in `OpenHymnal-C.abc`. Outputs Roman numeral analysis. Detects phrase structure patterns.
- **`analyze_octave_corruption.py`** — Compares original vs transposed ABC files to find octave comma redistribution bugs from `abctool.py` transposition.
- **`render_mutopia_shapes.py`** — Batch-renders Mutopia SATB hymns with Aiken 7-shape noteheads via LilyPond (`~/local/lilypond-2.24.4/`).
- **`abctool.py`** — Third-party ABC notation Swiss army knife (transposition, format conversion, etc.). Wraps `abcm2ps`, `abc2midi`, `gs`.

## Key Architecture Details

- **`harp_trainer.html`** is entirely self-contained (CSS + JS in one file). The JS is organized into numbered sections: VexFlow imports → config → music constants → shape notes → pedal state → drill logic → sight-read engine → rendering.
- **Chord ambiguity filtering**: `unambiguousChordTypes` auto-excludes chord type pairs that share a pitch class set (e.g., sus2/sus4). Used in Chord ID mode; Interval→Chord mode shows all types with multi-select.
- **Pedal state** drives everything: key signature, available pitches for drills, and pitch class constraints for the WASM detector.
- The WASM detector (`harp-pitch-wasm/src/detector.rs`) takes the 7 pedal pitch classes as input to constrain FFT peak matching to harp-playable notes.
- `OpenHymnal-C.abc` is the entire OpenHymnal transposed to C major via `abctool.py -t`. `OpenHymnal2014.06.abc` is the original.
