# Harp Voicing Engine — Software Requirements Specification

## Purpose

Given SATB hymn voices from the OpenHymnal, expand the 4-note harmony to 6-8 note harp voicings that maximize fullness and musical quality within physical hand constraints. The soprano melody is placed on the RH thumb; the remaining notes are distributed across both hands with the LH carrying the load.

## Inputs

- SATB voices from `OpenHymnal.abc` (284 hymns, 4 voices: S1V1, S1V2, S2V1, S2V2)
- Fallback: `harp_leadsheets.abc` (melody + Roman numeral chord symbols)
- Mode: Ionian through Locrian (offsets string numbers by 0-6)

## Output Format

Per chord change, the engine produces:

    RH chord label     (e.g., II~2)
    LH chord label     (e.g., V)
    RH string numbers  thumb (melody), index, middle, ring — high to low
    LH string numbers  thumb, index, middle, ring — high to low

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
- RH thumb = soprano melody note (highest RH note).
- Ring finger = lowest note in each hand.
- Fingers can skip strings within the span.

### Span Limits
- Max hand span: 10 diatonic strings per hand.
- Gap between hands: 0-3 strings (RH ring to LH thumb).
- Max total voicing span: 10 + 3 + 10 = 23 strings.
- LH and RH must not overlap.

### Harp Size
- 29 strings minimum: 23 (max voicing) + 6 (mode offset for all 7 church modes).
- All voicings must fit within 29 strings regardless of mode.
- **The engine must verify: highest string number + 6 (max mode offset) ≤ 29.**

## Voicing Rules

### R1: LH carries the load
LH plays 3-4 notes, RH plays 2-3 notes (melody on thumb + 1-2 support notes). The melody hand stays light so the harpist can focus on the tune. A 6-note chord should split LH=4 + RH=2. An 8-note chord: LH=4 + RH=4.

### R2: Maximize span — use inversions
Prefer open voicings that spread chord tones across the full 10-string span per hand. Use 1st, 2nd, and 3rd inversions to avoid:
- Duplicate pitch classes within the same hand
- Simple octave doublings when an inversion would add a new chord tone
- Clustered voicings when spread voicings are available

A voicing using 5 strings when 10 are available must score lower than one using 8-10.

### R3: Melody stands alone in RH
Do not double the melody note's pitch class in the RH accompaniment fingers. The melody should stand out, not be buried. LH doubling of the melody pitch class at a lower octave is acceptable and adds fullness.

### R4: Add tone selection
When filling a hand to use all fingers, choose add tones that complement the harmony:
- If melody is the root: prefer add 6, 9, or A.
- If melody is the 3rd: prefer add 2, 6, or 9.
- If melody is the 5th: prefer add 2, 4, or 9.
- Avoid add tones that create a semitone clash with the melody.
- Prefer tones already present in the SATB harmony.

### R5: Bass note
- Root in bass (LH ring finger) by default.
- If inversion is specified (`^N`), the Nth chord tone must be in the bass.
- Root doubling in the bass register adds fullness.

### R6: Voice leading (context-aware)
The engine must consider the previous chord when choosing a voicing:
- Minimize total finger movement between consecutive chords.
- Prefer smooth bass motion (stepwise or small intervals over large leaps).
- Common tones between consecutive chords should stay in the same hand position when possible.
- Penalize bass leaps larger than a 5th unless harmonically motivated (e.g., V→I).

### R7: Gap usage
- Gap 0: dense, powerful sound. Use for climactic moments, final chords.
- Gap 1-2: balanced, standard voicing. Default preference.
- Gap 3: open, ethereal. Use for quiet passages, transitions.

### R8: Even distribution
Chord tones should be distributed evenly across the register. Penalize voicings where all notes cluster in one octave with a gap elsewhere. The ideal voicing fills the space between LH bass and RH melody without holes larger than an octave.

### R9: Chromatic chord substitution
The diatonic harp cannot play accidentals without a pedal/lever change. When the SATB produces a chromatic chord, substitute the nearest diatonic chord by:
1. Finding diatonic chords that share the most pitch classes.
2. Preferring substitutes that preserve harmonic function (dominant, subdominant, tonic).
3. Choosing the substitute whose root is closest by diatonic step.

### R10: Chromatic melody substitution
When the soprano melody contains an accidental:
1. Replace with the nearest diatonic note (prefer the resolution direction).
2. Prefer the diatonic note in the direction of melodic motion.

## Scoring Function

Each candidate voicing receives a score. Higher = better.

| Factor | Points | Description |
|--------|--------|-------------|
| LH has 4 notes | +20 | Heavy lifting |
| LH has 3 notes | +10 | Acceptable |
| RH has 2-3 notes | +10 | Melody stays light |
| RH has 4 notes + LH < 3 | -15 | Wrong hand balance |
| Span width | +3 per string of hand span | Reward open voicings |
| Chord tone coverage | +15 per unique pitch class | All tones covered |
| Bass correctness | +20 | Correct bass note per inversion |
| Root in bass register | +10 per root below MIDI 60 | Bass reinforcement |
| Both hands used | +15 | Full sound |
| Melody not doubled in RH | +10 | Melody stands out |
| Melody doubled in RH | -15 | Buried melody |
| Gap > octave between adjacent | -5 per excess string | Penalize holes |
| Total notes 6-8 | +20 | Sweet spot |
| Total notes < 5 | -20 | Too thin |
| Voice leading from prev chord | +20 max | Smooth motion |
| Bass leap > 5th from prev | -5 per excess string | Penalize clunky bass |
| Voicing exceeds 29 strings with mode | -100 | Unplayable |

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

### Pure JavaScript: `HarpHymnal.html`

Single-file browser app. No build step, no WASM, no server-side processing.

    Music Theory    — pitch class, scale, string mapping, all octaves
    Chord Parser    — Roman numeral chord symbol parser (391 chord vocabulary)
    Chord Identifier — identifies chord name from a set of MIDI notes
    SATB Parser     — extracts 4 voices from OpenHymnal ABC format
    Voicing Engine  — expands SATB to 6-8 notes, splits for harp hands
    Lead Sheet Parser — fallback for melody + chord symbol ABC files
    Renderer        — measure-fit line wrapping, chord fraction display

### Data Flow

    OpenHymnal.abc (SATB, 284 hymns)
      → JS SATB parser extracts 4 voices per hymn
      → Voicing engine: SATB notes → expand to 6-8 → split RH/LH
      → Chord identifier: label each hand's chord shape
      → Renderer: chord fractions + string number stacks + barlines

### Implementation Status

| Rule | Status | Notes |
|------|--------|-------|
| R1: LH heavy lifting | ✅ Implemented | Scoring prefers LH 3-4 |
| R2: Maximize span | ⚠️ Partial | Scored but not aggressive enough |
| R3: Melody avoidance | ❌ Not implemented | Need to penalize RH melody doubling |
| R4: Add tone selection | ❌ Not implemented | Currently adds any chord tone |
| R5: Bass note | ✅ Implemented | Prefers root in bass |
| R6: Voice leading | ❌ Not implemented | No previous-chord context |
| R7: Gap usage | ⚠️ Partial | Fixed preference for gap 1 |
| R8: Even distribution | ❌ Not implemented | No gap penalty between notes |
| R9: Chromatic substitution | ❌ Not implemented | Chromatic notes dropped |
| R10: Chromatic melody | ❌ Not implemented | Non-diatonic melody skipped |
| 29-string check | ❌ Not implemented | No total range validation |
