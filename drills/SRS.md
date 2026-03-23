# Harp Voicing Engine — Software Requirements Specification

## Purpose

Given a hymn melody and chord progression, generate optimal 6-8 note harp voicings for both hands that maximize fullness and musical quality within physical hand span constraints.

## Inputs

- Melody line (soprano) as ABC notation with chord symbols
- Source: `harp_leadsheets.abc` (293 hymns, generated from OpenHymnal SATB)
- Mode: Ionian through Locrian (offsets string numbers by 0-6)

## Output Format

Per chord change, the engine produces:

    'I0|IM7^1       terse chord notation (LH | RH)
     12              string numbers stacked top-to-bottom
     10              (relative to chord root)
      8
      |              hand separator
      5
      3
      1

### Notation Grammar

    [']LH_chord GAP| RH_chord

- `'` prefix = LH uses 3 fingers instead of 4
- GAP = 0-3 (strings between LH thumb and RH ring)
- `|` separates hands
- Each chord from the 391 diatonic vocabulary

### Chord Symbol Components

| Symbol | Meaning |
|--------|---------|
| Case | uppercase = major 3rd, lowercase = minor 3rd |
| `+` | augmented (raised 5th) |
| `-` | diminished (lowered 5th) |
| `M7` | major 7th |
| `7` | dominant 7th |
| `%7` | half-diminished 7th |
| `-7` | diminished 7th |
| `~2` / `~4` | suspended 2nd / 4th |
| `^1` `^2` `^3` | 1st, 2nd, 3rd inversion |
| `2,4,6,8,9,A` | add note (diatonic degree, A=10 hex) |

Sus + add: first digit after `~` is sus type, next is add (`I~42` = sus4 + add2).

## Physical Constraints

### Hand Rules
- 4 fingers per hand: thumb (1), index (2), middle (3), ring (4). No pinky.
- RH thumb = melody note (highest RH note, fixed from lead sheet).
- Ring finger = lowest note in each hand.
- Fingers can skip strings within the span.

### Span Limits
- Max hand span: 10 diatonic strings per hand.
- Gap between hands: 0-3 strings (LH thumb to RH ring).
- Max total voicing span: 10 + 3 + 10 = 23 strings.
- LH and RH must not overlap.

### Harp Size
- 29 strings minimum: 23 (max voicing) + 6 (mode offset for all 7 church modes).
- All voicings must fit within 29 strings regardless of mode.

## Voicing Rules

### R1: Maximize finger usage
Always use all 8 fingers (4 per hand) unless the chord has fewer than 4 unique tones AND no beneficial add tone exists. A 3-note chord should be filled to 4 fingers with an add tone (2, 4, 6, 8, 9, or A).

### R2: Maximize span
Prefer open voicings that spread chord tones across the full available span. A voicing using 5 strings when 10 are available must score lower than one using 8-10.

### R3: Melody avoidance in accompaniment
Do not double the melody note's pitch class in the RH accompaniment fingers. The melody should stand out, not be buried. LH doubling of the melody pitch class at a lower octave is acceptable.

### R4: Add tone selection
When filling a 3-note chord to 4 fingers, choose add tones that complement the melody:
- If melody is the root: prefer add 6, 9, or A.
- If melody is the 3rd: prefer add 2, 6, or 9.
- If melody is the 5th: prefer add 2, 4, or 9.
- Avoid add tones that create a semitone clash with the melody.

### R5: Bass note
- Root in bass (LH ring finger) by default.
- If inversion is specified (`^N`), the Nth chord tone must be in the bass.
- Root doubling in the bass register adds fullness.

### R6: Voice leading (context-aware)
The engine must consider the previous and next chord when choosing a voicing:
- Minimize total finger movement between consecutive chords.
- Prefer smooth bass motion (stepwise or small intervals over large leaps).
- Common tones between consecutive chords should stay in the same hand position when possible.
- Penalize bass leaps larger than a 5th unless harmonically motivated (e.g., V-I).

### R7: Gap usage
- Gap 0: dense, powerful sound. Use for climactic moments, final chords.
- Gap 1-2: balanced, standard voicing. Default preference.
- Gap 3: open, ethereal. Use for quiet passages, transitions.
- The engine should vary the gap based on musical context (dynamics, phrase position).

### R8: Even distribution
Chord tones should be distributed evenly across the register. Penalize voicings where all notes cluster in one octave with a gap elsewhere. The ideal voicing fills the space between LH bass and RH melody without holes.

### R9: Chromatic chord substitution
The diatonic harp cannot play accidentals without a pedal/lever change. When the SATB analysis produces a chromatic chord (root has `b` or `#` prefix), the engine must substitute the nearest diatonic chord by:
1. Finding diatonic chords that share the most pitch classes with the chromatic chord.
2. Preferring substitutes that preserve the harmonic function (dominant, subdominant, tonic).
3. Choosing the substitute whose root is closest (by diatonic step) to the original root.

Common substitutions (key of C):

| Chromatic | Function | Substitute | Reason |
|-----------|----------|------------|--------|
| bVII (Bb) | dominant approach | V or vi | shares B-D with V, shares common tones with vi |
| #iv- (F#dim) | dominant prep | V or vii- | leading tone function |
| bII (Db) | Neapolitan | ii or I | subdominant function |
| bIII (Eb) | mediant | iii or IV | closest diatonic root |
| bVI (Ab) | submediant | vi or IV | shared tones |
| #IV+ | chromatic passing | IV or V | closest root |
| bVI+ | augmented 6th | V or vi | resolves to dominant |

### R10: Chromatic melody substitution
When the soprano melody contains an accidental (a note not in the current diatonic scale):
1. Replace with the nearest diatonic note (prefer the note the chromatic tone resolves to).
2. The chromatic melody note is "implied" — the listener's ear fills it in from context.
3. Prefer the diatonic note in the direction of melodic motion (e.g., C# resolving up to D → play D; Bb resolving down to A → play A).
4. If the chromatic note is a passing tone between two diatonic notes, either neighbor is acceptable.

## Scoring Function

Each candidate voicing receives a score. Higher = better.

| Factor | Points | Description |
|--------|--------|-------------|
| Finger usage | +5 per finger used | 8 fingers = 40 points |
| Span width | +2 per string of total span | 23 strings = 46 points |
| Chord tone coverage | +15 per unique chord tone present | All 4 tones = 60 points |
| Bass correctness | +20 | Correct bass note per inversion |
| Root in bass register | +10 per root below MIDI 60 | Bass reinforcement |
| Both hands used | +15 | LH and RH both have notes |
| Voice leading | +20 max | Smooth motion from previous chord |
| Melody not doubled in RH | +10 | Melody stands out |
| Gap > octave between adjacent notes | -3 per excess string | Penalize holes |
| Unused fingers | -10 per unused finger | Penalize thin voicings |
| Bass leap > 5th from previous | -5 per excess string | Penalize clunky bass |

## Mode Support

The 7 church modes are supported via string number offset:

| Mode | Offset | Root degree |
|------|--------|-------------|
| Ionian | +0 | 1 |
| Dorian | +1 | 2 |
| Phrygian | +2 | 3 |
| Lydian | +3 | 4 |
| Mixolydian | +4 | 5 |
| Aeolian | +5 | 6 |
| Locrian | +6 | 7 |

When mode changes:
- All string numbers shift by the mode offset.
- Chord names are re-analyzed relative to the new root.
- Shapes and finger patterns remain identical.

## Architecture

### Rust/WASM Crate: `harp-voicing-wasm/`

    src/
      lib.rs            WASM API (HarpVoicer)
      chord_symbol.rs   chord symbol parser (391 chord vocabulary)
      voicing.rs        voicing enumerator + scorer
      harp.rs           physical constraints, types
      music.rs          pitch class, scale, string mapping

### WASM API

    HarpVoicer.new(key_root: u8) -> Self
    HarpVoicer.set_key_root(key_root: u8)
    HarpVoicer.suggest_voicings(melody_midi: u8, chord: &str, max: usize) -> JSON
    HarpVoicer.suggest_voicings_in_context(
        melody_midi: u8,
        chord: &str,
        prev_voicing: &str,   // JSON of previous voicing for voice leading
        max: usize
    ) -> JSON

### HTML App: `harp_voicing.html`

- Hymn selector dropdown
- Mode selector (7 buttons: Ion Dor Phr Lyd Mix Aeo Loc)
- Melody flow: notes with chord symbols above, barlines between measures
- Stacked string numbers below each chord
- Terse LH|RH chord notation above string numbers
- Reference table: ASCII notation to traditional notation with pronunciation

### Data Flow

    OpenHymnal.abc
      → harp_leadsheet.py (Python, SATB → melody + chords)
      → harp_leadsheets.abc
      → harp_voicing.html loads ABC
      → WASM engine picks LH|RH voicings per chord
      → display: chord names + stacked string numbers + mode offset
