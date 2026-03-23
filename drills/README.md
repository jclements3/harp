# Harp Chord Notation

7-bit ASCII Roman numeral chord notation for diatonic harp.
391 chords, all within 10-string hand span, same shapes in any key.

## Notation

| Symbol | Meaning |
|--------|---------|
| Case | uppercase = major 3rd, lowercase = minor 3rd |
| `+` | augmented (raised 5th) |
| `-` | diminished (lowered 5th) |
| `M7` | major 7th |
| `7` | dominant 7th |
| `%7` | half-diminished 7th |
| `-7` | diminished 7th |
| `~2` | suspended 2nd (replaces 3rd) |
| `~4` | suspended 4th (replaces 3rd) |
| `^1` `^2` `^3` | 1st, 2nd, 3rd inversion |
| digit `2,4,6,8,9,A` | add note (scale degree, A=10) |

## Parsing Rules

- `~` is always followed by `2` or `4` (sus type)
- Additional digit after sus type = add note: `I~42` = sus4 + add 2
- Single digit/hex on a triad = add note: `I6` = major triad + add 6th
- `7` on a triad = 7th chord (not add): `V7` = dominant 7th
- `^N` at the end = Nth inversion

## Diatonic Chords (Key of C shown, shapes identical in all keys)

### Triads (7)

| Chord | Notes | Strings |
|-------|-------|---------|
| I | C E G | 5 |
| ii | D F A | 5 |
| iii | E G B | 5 |
| IV | F A C | 5 |
| V | G B D | 5 |
| vi | A C E | 5 |
| vii- | B D F | 5 |

### Seventh Chords (7)

| Chord | Notes | Strings |
|-------|-------|---------|
| IM7 | C E G B | 7 |
| ii7 | D F A C | 7 |
| iii7 | E G B D | 7 |
| IVM7 | F A C E | 7 |
| V7 | G B D F | 7 |
| vi7 | A C E G | 7 |
| vii%7 | B D F A | 7 |

### Suspended Chords (12)

| Chord | Notes |
|-------|-------|
| I~2 | C D G |
| I~4 | C F G |
| ii~2 | D E A |
| ii~4 | D G A |
| iii~4 | E A B |
| IV~2 | F G C |
| IV~4 | F B C |
| V~2 | G A D |
| V~4 | G C D |
| vi~2 | A B E |
| vi~4 | A D E |
| vii~2 | B C F |

### Add Chords (102)

3-note chord + 1 added scale degree. Uses 4 fingers.

Triads add 2, 4, 6, 8, 9, A (6 each = 42):

    I2 I4 I6 I8 I9 IA
    ii2 ii4 ii6 ii8 ii9 iiA
    iii2 iii4 iii6 iii8 iii9 iiiA
    IV2 IV4 IV6 IV8 IV9 IVA
    V2 V4 V6 V8 V9 VA
    vi2 vi4 vi6 vi8 vi9 viA
    vii-2 vii-4 vii-6 vii-8 vii-9 vii-A

Sus2 add 4, 6, 8, 9, A (5 each = 30):

    I~24 I~26 I~28 I~29 I~2A
    ii~24 ii~26 ii~28 ii~29 ii~2A
    IV~24 IV~26 IV~28 IV~29 IV~2A
    V~24 V~26 V~28 V~29 V~2A
    vi~24 vi~26 vi~28 vi~29 vi~2A
    vii~24 vii~26 vii~28 vii~29 vii~2A

Sus4 add 2, 6, 8, 9, A (5 each = 30):

    I~42 I~46 I~48 I~49 I~4A
    ii~42 ii~46 ii~48 ii~49 ii~4A
    iii~42 iii~46 iii~48 iii~49 iii~4A
    IV~42 IV~46 IV~48 IV~49 IV~4A
    V~42 V~46 V~48 V~49 V~4A
    vi~42 vi~46 vi~48 vi~49 vi~4A

## Inversions

Every chord can be played in multiple positions:

- 3-note chords: root, ^1, ^2 (3 positions)
- 4-note chords: root, ^1, ^2, ^3 (4 positions)

All inversions fit within 10-string hand span.

## Chord Count

| Type | Base | x Positions | Total |
|------|------|-------------|-------|
| Triads | 7 | x 3 | 21 |
| Sevenths | 7 | x 4 | 28 |
| Sus | 12 | x 3 | 36 |
| Triad adds | 42 | x 3 | 126 |
| Sus adds | 60 | x 3 | 180 |
| **Total** | | | **391** |

## Harp Hand Rules

- 4 fingers per hand: thumb (1), index (2), middle (3), ring (4). No pinky.
- RH thumb = melody note (highest, fixed).
- Ring finger = lowest note in each hand.
- Max hand span: 10 diatonic strings.
- Fingers can skip strings within the span.
- LH thumb sits at or below RH ring finger (no overlap).
- Changing key = move pedals, same shapes.
