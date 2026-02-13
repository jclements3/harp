# Transposition Corruption - Session Notes

## Project Goal
Create a printed shaped-note hymnal book with Aiken 7-shape noteheads on all
SATB (4-part harmony) music. Maximize song count. Include Roman numeral chord
analysis and form labels (letter-tagged phrases with Nx repeat notation).

## Current Song Sources

### 1. OpenHymnal (293 SATB hymns)
- File: `OpenHymnal-C.abc` — all transposed to K:C
- Rendered: `OpenHymnal-C.pdf` (514 pages)
- Shape notes via PostScript overrides in `%%beginps`/`%%endps` block
- Chord analysis + form labels already added via `analyze_hymnal.py`
- Original (pre-transposition): `OpenHymnal2014.06.abc` (various keys)

### 2. Mutopia Project (93 hymn tunes + ~8 Bach SATB chorales)
- Cloned to `~/projects/harp/mutopia/` from GitHub
- Rendered with `render_mutopia_shapes.py` → `~/projects/harp/mutopia_shapes/pdf/`
- 454 total PDFs rendered, 93 are hymn tunes usable for the book
- Uses LilyPond's native `\aikenHeads` command (installed at `~/local/lilypond-2.24.4/`)
- Some overlap with OpenHymnal (e.g., Amazing Grace, Old 100th) — needs dedup

## Two Systematic Bugs in abctool.py Transposer

### Bug 1: Rhythm Corruption
- **177 songs affected** (165 rhythm-only + 22 with both bugs)
- **~1,199 corrupted patterns**
- Pattern: `X/Y/` (two eighth notes) became `XY//` (quarter + sixteenth)
- Examples: `A/G/` → `AG//`, `C/B,/` → `CB,//`
- Song #43 "Savior Of The Nations Come" was manually fixed (commit 8505622)

### Bug 2: Octave Comma Corruption in Chords
- **34 songs affected** (12 octave-only + 22 with both bugs)
- **350 corrupted chords**
- Pattern: commas inside `[...]` chords get redistributed to the last note
- Example: `[F,,F,]` (F2+F3, normal bass) → `[FF,,,]` (F4+F0, absurd)
- Happens even for songs already in K:C (zero transposition!)
- "Away In A Manger" (X:51) is a clear example of bass clef abuse
- Analysis script: `analyze_octave_corruption.py` (commit 9dc5a00)

### Combined Impact
| Bug type | Songs |
|---|---|
| Only rhythm (`XY//`) | 155 |
| Only octave (chord commas) | 12 |
| Both | 22 |
| **Total corrupted** | **177 out of 293** |

## Recommended Fix Approach
Since 177 songs (60%) are affected, the cleanest fix is to **fix `abctool.py`'s
transpose() function and re-transpose from the original**. Both bugs are in the
same code path (around line 1909+ in abctool.py). Alternative: write a repair
script that compares original vs transposed and patches corrupted patterns.

## Key Files
| File | Purpose |
|---|---|
| `OpenHymnal-C.abc` | Main hymnal, 293 SATB in K:C (has corruption) |
| `OpenHymnal2014.06.abc` | Original hymnal, various keys (correct) |
| `OpenHymnal-C.pdf` | Rendered hymnal with shapes, 514 pages |
| `abctool.py` | Transposition tool (source of both bugs) |
| `analyze_hymnal.py` | Chord analysis, form labels, cycle detection |
| `analyze_octave_corruption.py` | Octave comma corruption analysis |
| `render_mutopia_shapes.py` | Batch Mutopia rendering with Aiken shapes |
| `TODO_rhythm_fix.md` | This file |
| `mutopia/` | Cloned Mutopia Project (17,137 files) |
| `mutopia_shapes/` | Rendered Mutopia PDFs (454 files) |
| `~/local/lilypond-2.24.4/` | LilyPond binary install |

## Shape Note Implementation Details
- **System**: Aiken 7-shape (modern standard, NOT Funk)
- **Shapes**: do=triangle-up, re=crescent, mi=diamond, fa=flag, sol=oval, la=square, ti=inverted-triangle
- **ABC pipeline**: PostScript overrides for `/hd`, `/Hd`, `/HD` notehead functions in `%%beginps`/`%%endps`
- **LilyPond pipeline**: Native `\aikenHeads` command, no hacking needed
- Degree formula: `((staffPos + offset) mod 7 + 7) mod 7 + 1` where offset=2 for treble, 4 for bass

## Next Steps (Priority Order)
1. **Fix transposer and re-transpose** — Fix both bugs in abctool.py, re-run on OpenHymnal2014.06.abc
2. **Deduplicate** — Identify overlap between OpenHymnal and Mutopia hymns
3. **Combine into unified book** — Merge all PDFs with table of contents
4. **Re-render OpenHymnal-C.pdf** — After fixing corruption, regenerate with shapes + analysis
