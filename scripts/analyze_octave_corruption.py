#!/usr/bin/env python3
"""
Analyze octave comma corruption in transposed OpenHymnal ABC file.

Compares OpenHymnal2014.06.abc (original, various keys) with
OpenHymnal-C.abc (transposed to K:C) to find chord corruption
where octave commas get incorrectly redistributed among notes.

Example corruption:
  Original:  [F,,F,]   (two notes: F,, and F,)
  Transposed: [FF,,,]   (two notes: F and F,,, -- commas piled on last note)
"""

import re
import sys
from collections import defaultdict

ORIGINAL = "/home/james.clements/projects/harp/OpenHymnal2014.06.abc"
TRANSPOSED = "/home/james.clements/projects/harp/OpenHymnal-C.abc"


def parse_songs(filepath):
    """Parse an ABC file into a dict of songs keyed by X number.
    Each song is a dict with 'title', 'key', 'voices' (dict of voice_name -> list of music lines),
    and 'bass_voices' (set of voice names with clef=bass).
    """
    songs = {}
    current_song = None
    current_voice = None
    header_done = False

    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            stripped = line.strip()

            # New song starts with X:
            if stripped.startswith('X:'):
                x_num = stripped[2:].strip()
                current_song = {
                    'x': x_num,
                    'title': '',
                    'key': '',
                    'voices': defaultdict(list),
                    'bass_voices': set(),
                    'all_voice_defs': {},
                }
                songs[x_num] = current_song
                header_done = False
                current_voice = None

            elif current_song is None:
                continue

            # Title
            elif stripped.startswith('T:') and not current_song['title']:
                current_song['title'] = stripped[2:].strip()

            # Voice definition (in header area)
            elif stripped.startswith('V:') and not header_done:
                # Parse voice name and check for clef=bass
                rest = stripped[2:].strip()
                voice_name = re.split(r'\s+', rest)[0]
                current_song['all_voice_defs'][voice_name] = stripped
                if 'clef=bass' in stripped.lower():
                    current_song['bass_voices'].add(voice_name)

            # Key signature marks end of header
            elif stripped.startswith('K:') and not header_done:
                current_song['key'] = stripped[2:].strip().split('%')[0].strip()
                header_done = True

            # Voice music line: [V: name] music...
            elif stripped.startswith('[V:') and header_done:
                m = re.match(r'\[V:\s*(\S+)\s*\]', stripped)
                if m:
                    current_voice = m.group(1)
                    # Get the music content after the voice tag
                    music = stripped[m.end():].strip()
                    if music:
                        current_song['voices'][current_voice].append(music)

            # Empty line resets song context (song boundary)
            elif stripped == '' and current_song:
                # Don't reset -- songs have blank lines between sections sometimes
                pass

            # %%newpage signals song boundary
            elif stripped == '%%newpage':
                current_song = None
                current_voice = None
                header_done = False

    return songs


def find_abc_chords(music_line):
    """Find all bracketed chords [...] in a music line.
    Returns list of (chord_string, start_pos, end_pos).
    A chord is [...] where the content contains notes (not [V:...] or [K:...]).
    """
    chords = []
    i = 0
    while i < len(music_line):
        if music_line[i] == '[':
            # Check if this is a special bracket like [V:, [K:, [M:, [I:, [P:, [Q:
            rest = music_line[i+1:]
            if rest and rest[0] in 'VKMIPQLr':
                # Skip to closing bracket
                j = music_line.find(']', i+1)
                if j >= 0:
                    i = j + 1
                else:
                    i += 1
                continue
            # Find matching ]
            j = music_line.find(']', i+1)
            if j >= 0:
                chord_content = music_line[i:j+1]
                # Check this looks like a note chord (contains note letters)
                if re.search(r'[a-gA-G]', chord_content):
                    chords.append((chord_content, i, j+1))
                i = j + 1
            else:
                i += 1
        else:
            i += 1
    return chords


def parse_chord_notes(chord_str):
    """Parse an ABC chord like [F,,F,] into individual notes.
    Returns list of note strings.
    """
    # Strip brackets
    content = chord_str[1:-1]
    notes = []
    # ABC note pattern: optional accidentals, note letter, optional octave markers, optional duration
    note_pattern = re.compile(r'([_^=]*[a-gA-G][,\']*)([\d]*[/]?[\d]*)')
    pos = 0
    while pos < len(content):
        m = note_pattern.match(content, pos)
        if m:
            note_with_dur = m.group(0)
            notes.append(note_with_dur)
            pos = m.end()
        else:
            pos += 1  # skip non-note chars
    return notes


def count_commas(note_str):
    """Count the net octave lowering in a note string (commas minus apostrophes)."""
    return note_str.count(',') - note_str.count("'")


def note_pitch_class(note_str):
    """Extract just the pitch letter (uppercase) from a note."""
    m = re.search(r'[a-gA-G]', note_str)
    if m:
        return m.group(0).upper()
    return '?'


def note_is_lowercase(note_str):
    """Check if the note letter is lowercase (higher octave)."""
    m = re.search(r'[a-gA-G]', note_str)
    if m:
        return m.group(0).islower()
    return False


def note_to_octave(note_str):
    """Convert ABC note to approximate octave number.
    C, = C3 (cello range), C = C4 (middle C), c = C5, c' = C6
    """
    m = re.search(r'[a-gA-G]', note_str)
    if not m:
        return -1
    letter = m.group(0)
    if letter.isupper():
        base_octave = 4  # C-B = octave 4
    else:
        base_octave = 5  # c-b = octave 5

    commas = note_str.count(',')
    apostrophes = note_str.count("'")
    return base_octave - commas + apostrophes


def has_extreme_low_note(chord_str, min_commas=3):
    """Check if a chord has any note with 3+ commas (at or below C1/B0)."""
    notes = parse_chord_notes(chord_str)
    for note in notes:
        total_commas = note.count(',')
        total_apos = note.count("'")
        net_low = total_commas - total_apos
        # Check for uppercase letter with many commas, or lowercase with even more
        m = re.search(r'[a-gA-G]', note)
        if m:
            letter = m.group(0)
            if letter.isupper() and net_low >= min_commas:
                return True
            if letter.islower() and net_low >= (min_commas + 1):
                return True
    return False


def total_comma_count(chord_str):
    """Count total commas in the chord (across all notes)."""
    content = chord_str[1:-1]  # strip brackets
    return content.count(',')


def comma_distribution(chord_str):
    """Return the distribution of commas per note in a chord."""
    notes = parse_chord_notes(chord_str)
    return [n.count(',') for n in notes]


def analyze_chord_corruption(orig_chord, trans_chord):
    """Analyze whether a transposed chord shows comma redistribution corruption.

    The corruption pattern: in original, commas are distributed among notes correctly,
    e.g. [F,,F,] (2 commas on first note, 1 on second).
    In transposed, all commas pile up on one note:
    e.g. [FF,,,] (0 commas on first note, 3 on second).

    The total comma count is preserved but the distribution changes.
    """
    orig_notes = parse_chord_notes(orig_chord)
    trans_notes = parse_chord_notes(trans_chord)

    if len(orig_notes) != len(trans_notes):
        return None  # Can't compare if note count differs

    if len(orig_notes) < 2:
        return None  # Single-note "chords" can't have redistribution

    orig_commas = [n.count(',') for n in orig_notes]
    trans_commas = [n.count(',') for n in trans_notes]

    orig_total = sum(orig_commas)
    trans_total = sum(trans_commas)

    # Check if comma distribution changed
    if orig_commas != trans_commas and (orig_total > 0 or trans_total > 0):
        # Check if it looks like redistribution (total may differ due to transposition
        # changing octaves, but the key signal is commas piling on one note)
        max_orig = max(orig_commas) if orig_commas else 0
        max_trans = max(trans_commas) if trans_commas else 0

        return {
            'orig_chord': orig_chord,
            'trans_chord': trans_chord,
            'orig_notes': orig_notes,
            'trans_notes': trans_notes,
            'orig_commas': orig_commas,
            'trans_commas': trans_commas,
            'orig_total': orig_total,
            'trans_total': trans_total,
            'max_orig': max_orig,
            'max_trans': max_trans,
            'extreme_low': has_extreme_low_note(trans_chord),
        }

    return None


def extract_chords_from_voice_lines(voice_lines):
    """Extract all chords from a list of voice music lines, preserving order."""
    all_chords = []
    for line in voice_lines:
        chords = find_abc_chords(line)
        for chord_str, start, end in chords:
            all_chords.append(chord_str)
    return all_chords


def check_rhythm_corruption(orig_line, trans_line):
    """Check if a voice line pair has the X/Y/ -> XY// rhythm corruption."""
    # Pattern: in original, X/Y/ (two eighth notes), in transposed XY//
    orig_eighth_pairs = re.findall(r'[_^=]*[a-gA-G][,\']*/[_^=]*[a-gA-G][,\']*/(?!\d)', orig_line)
    trans_corrupted = re.findall(r'[_^=]*[a-gA-G][,\']*[_^=]*[a-gA-G][,\']*//', trans_line)
    return len(orig_eighth_pairs) > 0 or len(trans_corrupted) > 0


def main():
    print("=" * 80)
    print("OCTAVE COMMA CORRUPTION ANALYSIS")
    print("Comparing original OpenHymnal2014.06.abc vs transposed OpenHymnal-C.abc")
    print("=" * 80)
    print()

    # Parse both files
    print("Parsing original file...")
    orig_songs = parse_songs(ORIGINAL)
    print(f"  Found {len(orig_songs)} songs")

    print("Parsing transposed file...")
    trans_songs = parse_songs(TRANSPOSED)
    print(f"  Found {len(trans_songs)} songs")
    print()

    # Find matching songs
    common_x = set(orig_songs.keys()) & set(trans_songs.keys())
    print(f"Songs in common (by X: number): {len(common_x)}")
    print()

    # Analyze each song
    songs_with_octave_corruption = []
    songs_with_rhythm_corruption = set()
    total_corrupted_chords = 0
    total_extreme_low_chords = 0
    all_corruptions = []

    # Track co-occurrence
    songs_octave_only = set()
    songs_rhythm_only = set()
    songs_both = set()

    for x_num in sorted(common_x, key=lambda x: int(x) if x.isdigit() else 0):
        orig = orig_songs[x_num]
        trans = trans_songs[x_num]

        song_has_octave_corruption = False
        song_chord_corruptions = []
        song_has_rhythm_corruption = False

        # Check ALL voices (not just bass), since chord corruption can happen anywhere
        all_voices = set(orig['voices'].keys()) & set(trans['voices'].keys())

        for voice in sorted(all_voices):
            orig_lines = orig['voices'][voice]
            trans_lines = trans['voices'][voice]

            # Check rhythm corruption for this voice
            for ol, tl in zip(orig_lines, trans_lines):
                if check_rhythm_corruption(ol, tl):
                    song_has_rhythm_corruption = True

            # Extract chords and compare
            orig_chords = extract_chords_from_voice_lines(orig_lines)
            trans_chords = extract_chords_from_voice_lines(trans_lines)

            # Compare chord by chord (same position)
            for idx, (oc, tc) in enumerate(zip(orig_chords, trans_chords)):
                corruption = analyze_chord_corruption(oc, tc)
                if corruption:
                    corruption['song_x'] = x_num
                    corruption['song_title'] = orig['title']
                    corruption['voice'] = voice
                    corruption['orig_key'] = orig['key']
                    corruption['is_bass'] = voice in orig['bass_voices']
                    song_has_octave_corruption = True
                    song_chord_corruptions.append(corruption)
                    total_corrupted_chords += 1
                    if corruption['extreme_low']:
                        total_extreme_low_chords += 1
                    all_corruptions.append(corruption)

        if song_has_octave_corruption:
            songs_with_octave_corruption.append({
                'x': x_num,
                'title': orig['title'],
                'orig_key': orig['key'],
                'corruptions': song_chord_corruptions,
            })

        if song_has_rhythm_corruption:
            songs_with_rhythm_corruption.add(x_num)

        # Track co-occurrence
        if song_has_octave_corruption and song_has_rhythm_corruption:
            songs_both.add(x_num)
        elif song_has_octave_corruption:
            songs_octave_only.add(x_num)
        elif song_has_rhythm_corruption:
            songs_rhythm_only.add(x_num)

    # =========================================================================
    # REPORT
    # =========================================================================

    print("=" * 80)
    print("RESULTS")
    print("=" * 80)
    print()

    print(f"Total songs analyzed:                {len(common_x)}")
    print(f"Songs with octave comma corruption:  {len(songs_with_octave_corruption)}")
    print(f"Total corrupted chords:              {total_corrupted_chords}")
    print(f"  - with extreme low notes (>=3 commas): {total_extreme_low_chords}")
    print()

    # =========================================================================
    # Breakdown by voice type
    # =========================================================================
    bass_corruptions = [c for c in all_corruptions if c['is_bass']]
    treble_corruptions = [c for c in all_corruptions if not c['is_bass']]
    print(f"Corrupted chords in bass voices:     {len(bass_corruptions)}")
    print(f"Corrupted chords in treble voices:   {len(treble_corruptions)}")
    print()

    # =========================================================================
    # Breakdown by corruption severity
    # =========================================================================
    print("--- Corruption severity breakdown ---")
    severity_counts = defaultdict(int)
    for c in all_corruptions:
        max_trans = c['max_trans']
        max_orig = c['max_orig']
        diff = max_trans - max_orig
        severity_counts[diff] += 1

    for diff in sorted(severity_counts.keys()):
        print(f"  Max comma increase of {diff:+d} on a single note: {severity_counts[diff]} chords")
    print()

    # =========================================================================
    # Example comparisons
    # =========================================================================
    print("=" * 80)
    print("EXAMPLE CORRUPTIONS (first 30)")
    print("=" * 80)
    print()

    shown = 0
    for c in all_corruptions:
        if shown >= 30:
            break
        print(f"  Song X:{c['song_x']} \"{c['song_title']}\" [V:{c['voice']}] "
              f"(orig key: {c['orig_key']}, bass={c['is_bass']})")
        print(f"    Original:   {c['orig_chord']}")
        print(f"      notes: {c['orig_notes']}  commas per note: {c['orig_commas']}")
        print(f"    Transposed: {c['trans_chord']}")
        print(f"      notes: {c['trans_notes']}  commas per note: {c['trans_commas']}")
        print(f"    Extreme low (>=3 commas on one note): {c['extreme_low']}")
        print()
        shown += 1

    # =========================================================================
    # Specific: chords with 3+ commas on a single note (almost certainly wrong)
    # =========================================================================
    print("=" * 80)
    print("CHORDS WITH 3+ COMMAS ON A SINGLE NOTE (definite artifacts)")
    print("=" * 80)
    print()

    extreme_examples = [c for c in all_corruptions if c['extreme_low']]
    print(f"Total: {len(extreme_examples)}")
    print()

    shown = 0
    for c in extreme_examples:
        if shown >= 20:
            print(f"  ... and {len(extreme_examples) - shown} more")
            break
        print(f"  X:{c['song_x']} \"{c['song_title']}\" [V:{c['voice']}]: "
              f"{c['orig_chord']} -> {c['trans_chord']}  "
              f"(commas: {c['orig_commas']} -> {c['trans_commas']})")
        shown += 1
    print()

    # =========================================================================
    # The specific corruption pattern: commas pushed to one note
    # =========================================================================
    print("=" * 80)
    print("COMMA REDISTRIBUTION PATTERN ANALYSIS")
    print("=" * 80)
    print()

    # Categorize the types of redistribution
    pattern_counts = defaultdict(int)
    pattern_examples = defaultdict(list)

    for c in all_corruptions:
        orig_dist = tuple(c['orig_commas'])
        trans_dist = tuple(c['trans_commas'])
        pattern = f"{orig_dist} -> {trans_dist}"
        pattern_counts[pattern] += 1
        if len(pattern_examples[pattern]) < 3:
            pattern_examples[pattern].append(c)

    print("Most common redistribution patterns:")
    for pattern, count in sorted(pattern_counts.items(), key=lambda x: -x[1])[:25]:
        print(f"  {pattern}: {count} occurrences")
        for ex in pattern_examples[pattern]:
            print(f"    e.g. {ex['orig_chord']} -> {ex['trans_chord']}")
    print()

    # =========================================================================
    # Co-occurrence with rhythm corruption
    # =========================================================================
    print("=" * 80)
    print("CO-OCCURRENCE WITH RHYTHM CORRUPTION (X/Y/ -> XY//)")
    print("=" * 80)
    print()

    print(f"Songs with ONLY octave corruption:   {len(songs_octave_only)}")
    print(f"Songs with ONLY rhythm corruption:   {len(songs_rhythm_only)}")
    print(f"Songs with BOTH corruptions:         {len(songs_both)}")
    print(f"Songs with rhythm corruption (total): {len(songs_with_rhythm_corruption)}")
    print()

    if songs_octave_only:
        print("Songs with octave corruption but NO rhythm corruption:")
        for x in sorted(songs_octave_only, key=lambda x: int(x) if x.isdigit() else 0):
            s = orig_songs[x]
            print(f"  X:{x} \"{s['title']}\" (key: {s['key']})")
        print()

    if songs_both:
        print("Songs with BOTH octave and rhythm corruption:")
        for x in sorted(songs_both, key=lambda x: int(x) if x.isdigit() else 0)[:20]:
            s = orig_songs[x]
            print(f"  X:{x} \"{s['title']}\" (key: {s['key']})")
        if len(songs_both) > 20:
            print(f"  ... and {len(songs_both) - 20} more")
        print()

    # =========================================================================
    # By original key
    # =========================================================================
    print("=" * 80)
    print("CORRUPTION BY ORIGINAL KEY")
    print("=" * 80)
    print()

    key_counts = defaultdict(lambda: {'songs': 0, 'chords': 0})
    for song_info in songs_with_octave_corruption:
        key = song_info['orig_key']
        key_counts[key]['songs'] += 1
        key_counts[key]['chords'] += len(song_info['corruptions'])

    total_by_key = defaultdict(int)
    for x in common_x:
        total_by_key[orig_songs[x]['key']] += 1

    print(f"{'Key':<15} {'Affected/Total Songs':<25} {'Corrupted Chords':<20}")
    print("-" * 60)
    for key in sorted(key_counts.keys()):
        affected = key_counts[key]['songs']
        total = total_by_key[key]
        chords = key_counts[key]['chords']
        print(f"{key:<15} {affected}/{total:<23} {chords:<20}")
    print()

    # =========================================================================
    # Songs in K:C originally (should have zero transposition, zero corruption)
    # =========================================================================
    print("=" * 80)
    print("SANITY CHECK: Songs originally in K:C")
    print("=" * 80)
    print()

    c_songs = [x for x in common_x if orig_songs[x]['key'].startswith('C')]
    c_corrupted = [x for x in c_songs if any(
        s['x'] == x for s in songs_with_octave_corruption)]
    print(f"Songs originally in C: {len(c_songs)}")
    print(f"  of those, with octave corruption: {len(c_corrupted)}")
    if c_corrupted:
        for x in c_corrupted:
            print(f"  X:{x} \"{orig_songs[x]['title']}\"")
    print()

    # =========================================================================
    # Summary
    # =========================================================================
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print()
    print(f"The transposition from various keys to K:C introduced octave comma")
    print(f"corruption in {len(songs_with_octave_corruption)} songs "
          f"({total_corrupted_chords} chords affected).")
    print()
    print(f"The corruption pattern is: in a chord like [F,,F,], the commas get")
    print(f"redistributed so they pile up on one note, e.g. [FF,,,].")
    print(f"This is because abctool's transposer processes notes inside chords")
    print(f"individually, and the normalize_octave() function changes the case")
    print(f"of the note letter but doesn't properly track that the commas/apostrophes")
    print(f"need to stay attached to their respective notes inside a chord context.")
    print()

    if len(songs_octave_only) == 0 and len(songs_both) > 0:
        print("FINDING: Octave corruption ALWAYS co-occurs with rhythm corruption.")
        print("         Every song with octave chord corruption also has rhythm corruption.")
    elif len(songs_octave_only) > 0:
        print(f"FINDING: Octave corruption is PARTIALLY INDEPENDENT of rhythm corruption.")
        print(f"         {len(songs_octave_only)} song(s) have octave corruption without rhythm corruption,")
        print(f"         while {len(songs_both)} have both types.")
    else:
        print("FINDING: No octave corruption detected (unexpected).")
    print()


if __name__ == '__main__':
    main()
