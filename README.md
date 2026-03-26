# Harp

Practice, hymnal, and music analysis tools for pedal harp.

## Apps

| App | Type | Description |
|-----|------|-------------|
| **harp-hymnal** | Rust native | SATB hymnal with voicing engine — double-click to run |
| **harp-trainer** | Rust native | Chord ID, interval, and sight-reading drills |
| **HarpHymnal.html** | Browser | Hymnal with leadsheet display and chord table |
| **harp_trainer.html** | Browser | Drill app with VexFlow notation and optional WASM mic input |
| **harp_voicing.html** | Browser | Interactive voicing helper |

## Quick Start

```bash
# Native hymnal app
cd harp-hymnal && cargo run

# Native trainer app
cd harp-trainer && cargo run

# Browser apps
python3 -m http.server 8000
# open http://localhost:8000
```

## Project Structure

```
harp-hymnal/        Rust eframe/egui native hymnal app
harp-trainer/       Rust eframe/egui native drill app
harp-pitch-wasm/    Rust WASM polyphonic pitch detector
harp-voicing-wasm/  Rust WASM voicing engine
harp-midi/          Rust MIDI generator
html/               Browser apps
scripts/            Python tools (analysis, PDF generation)
data/               ABC source files, voicings JSON
docs/               Specs, notes, reference PDFs
```

## Voicing Engine

Expands 4-note SATB hymn harmony to 6-8 note harp voicings:

- Melody on RH thumb, LH carries the load (3-4 notes)
- 10-string max hand span per hand
- Voice leading optimization between chords
- Chromatic-to-diatonic snapping for non-harp tones
- Scoring rules R1-R10 rank voicing quality

Input: `data/OpenHymnal.abc` (284 SATB hymns in C major)

## Chord Notation

7-bit ASCII Roman numeral system for diatonic harp:

| Symbol | Meaning |
|--------|---------|
| `V7` | dominant seventh |
| `ii` | minor (lowercase) |
| `ii%7` | half-diminished seventh |
| `bVII+` | flat-seven augmented |
| `I^1` | first inversion |
| `V~4` | sus4 |

## Build Requirements

- Rust stable (for native apps)
- `wasm-pack` (for WASM crates)
- Python 3 + `reportlab` (for PDF scripts)
