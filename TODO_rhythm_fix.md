# Rhythm Corruption Fix - TODO

## What was fixed
- Song #43 "Savior Of The Nations Come" - manually corrected all 9 corrupted patterns

## What still needs fixing
The transposition from OpenHymnal2014.06.abc (K:Bb etc) to OpenHymnal-C.abc (K:C) 
corrupted eighth-note pairs throughout the file:

- **165 songs affected** out of 293
- **~1,199 corrupted patterns** across 660 voice lines
- Pattern: `X/Y/` (two eighth notes) became `XY//` (quarter + sixteenth)
- Examples: `A/G/` → `AG//`, `C/B,/` → `CB,//`, `G/A/` → `GA//`

## How to fix
Compare each `[V:]` line in OpenHymnal-C.abc against the corresponding line in
OpenHymnal2014.06.abc. Wherever the original has `X/Y/`, restore that pattern
in the transposed file. The note NAMES may differ (due to transposition) but the
rhythm pattern (`/` placement) should match the original.

A Python script could automate this by:
1. Parsing both files song-by-song 
2. Matching voice lines by song number + voice name + phrase index
3. Copying the `/` placement from original to transposed

## Files
- Original: OpenHymnal2014.06.abc (various keys, correct rhythms)
- Transposed: OpenHymnal-C.abc (K:C, corrupted rhythms)
- Analysis script: analyze_hymnal.py
