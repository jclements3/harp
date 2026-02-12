#!/usr/bin/env python3
"""
Analyze OpenHymnal-C.abc: derive chord progressions from 4-part SATB harmony.

Detects phrase structure (A A B, A B A B, etc.), identifies chords per phrase
using downbeat template matching, and outputs Roman numeral analysis.
"""

import re
import sys
from collections import defaultdict

# Pitch class mapping (C=0)
NOTE_TO_PC = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
PC_TO_NAME = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'F#', 'G', 'Ab', 'A', 'Bb', 'B']

# Roman numeral scale degrees for key of C
PC_TO_ROMAN = {
    0: 'I', 1: 'bII', 2: 'II', 3: 'bIII', 4: 'III', 5: 'IV',
    6: '#IV', 7: 'V', 8: 'bVI', 9: 'VI', 10: 'bVII', 11: 'VII',
}

# Chord quality -> (lowercase for roman?, suffix)
QUALITY_MAP = {
    '':      (False, ''),
    'm':     (True,  ''),
    '7':     (False, '7'),
    'm7':    (True,  '7'),
    'maj7':  (False, 'maj7'),
    'dim':   (True,  '\u00b0'),
    'aug':   (False, '+'),
    'sus4':  (False, 'sus4'),
    'sus2':  (False, 'sus2'),
    'dim7':  (True,  '\u00b07'),
    'm7b5':  (True,  '\u00f87'),
}

# Chord templates: intervals from root -> quality suffix
CHORD_TYPES = [
    (frozenset({0, 4, 7}), ''),
    (frozenset({0, 3, 7}), 'm'),
    (frozenset({0, 4, 7, 10}), '7'),
    (frozenset({0, 3, 7, 10}), 'm7'),
    (frozenset({0, 4, 7, 11}), 'maj7'),
    (frozenset({0, 3, 6}), 'dim'),
    (frozenset({0, 4, 8}), 'aug'),
    (frozenset({0, 5, 7}), 'sus4'),
    (frozenset({0, 2, 7}), 'sus2'),
    (frozenset({0, 3, 6, 9}), 'dim7'),
    (frozenset({0, 3, 6, 10}), 'm7b5'),
]


def chord_to_roman(chord_name):
    """Convert a letter chord name like 'Am' or 'G7' to Roman numeral in C."""
    if not chord_name:
        return chord_name
    m = re.match(r'^([A-G][b#]?)(.*)', chord_name)
    if not m:
        return chord_name
    root_name, quality = m.group(1), m.group(2)
    root_pc = None
    for i, name in enumerate(PC_TO_NAME):
        if name == root_name:
            root_pc = i
            break
    if root_pc is None:
        return chord_name
    roman = PC_TO_ROMAN[root_pc]
    lower, suffix = QUALITY_MAP.get(quality, (False, quality))
    if lower:
        roman = roman.lower()
    return roman + suffix


def abc_note_to_pc(token):
    """Convert an ABC note token to pitch class (0-11), or None."""
    accidental = 0
    i = 0
    while i < len(token) and token[i] in '^_=':
        if token[i] == '^':
            accidental += 1
        elif token[i] == '_':
            accidental -= 1
        i += 1
    if i >= len(token):
        return None
    letter = token[i].upper()
    if letter not in NOTE_TO_PC:
        return None
    return (NOTE_TO_PC[letter] + accidental) % 12


def extract_measures(line, downbeat_only=False):
    """Extract list of pitch-class lists, one per measure, from a voice line."""
    line = re.sub(r'\[V:\s*[^\]]+\]', '', line)
    line = re.sub(r'\[[A-Z]:[^\]]*\]', '', line)
    line = re.sub(r'![a-zA-Z0-9]+!', '', line)

    measures = line.split('|')
    result = []
    for m in measures:
        pcs = []
        for match in re.finditer(r'[_^=]*[A-Ga-g][,\']*', m):
            pc = abc_note_to_pc(match.group())
            if pc is not None:
                pcs.append(pc)
        if downbeat_only and pcs:
            pcs = [pcs[0]]
        result.append(pcs)
    return result


def identify_chord(pitch_classes, bass_pc=None):
    """Identify best-matching chord from pitch classes with optional bass hint."""
    pcs = set(pitch_classes)
    if len(pcs) < 2:
        return PC_TO_NAME[list(pcs)[0]] if pcs else None

    best = None
    best_score = -999

    for root in range(12):
        intervals = frozenset((pc - root) % 12 for pc in pcs)
        for idx, (template, suffix) in enumerate(CHORD_TYPES):
            match = len(intervals & template)
            extra = len(intervals - template)
            missing = len(template - intervals)
            score = match * 3 - extra * 2 - missing
            if bass_pc is not None and root == bass_pc:
                score += 4
            if len(template) == 3:
                score += 2
            score -= idx * 0.1
            if score > best_score:
                best_score = score
                best = PC_TO_NAME[root] + suffix

    return best


def parse_file(filename):
    """Parse ABC file into list of song dicts with phrase-grouped voices."""
    with open(filename) as f:
        lines = f.readlines()

    songs = []
    song = None

    for line in lines:
        line = line.rstrip('\n')

        if line.startswith('X:'):
            if song:
                songs.append(song)
            song = {
                'number': line.split(':', 1)[1].strip(),
                'title': '',
                'key': 'C',
                'voice_lines': [],   # raw [V:] lines in order
                'first_voice': None, # name of the first voice seen
            }
        elif song is None:
            continue
        elif line.startswith('T:') and not song['title']:
            song['title'] = line.split(':', 1)[1].strip()
        elif line.startswith('K:'):
            song['key'] = line.split(':', 1)[1].strip()
        elif line.startswith('[V:'):
            vm = re.match(r'\[V:\s*(\S+)\]', line)
            if vm:
                vname = vm.group(1)
                if song['first_voice'] is None:
                    song['first_voice'] = vname
                song['voice_lines'].append((vname, line))

    if song:
        songs.append(song)
    return songs


def group_into_phrases(song):
    """Split a song's voice lines into phrases.

    A new phrase starts each time the first voice name reappears.
    Returns list of phrases, each phrase is a dict: {voice_name: music_line}.
    """
    first_voice = song['first_voice']
    if not first_voice:
        return []

    phrases = []
    current = {}

    for vname, line in song['voice_lines']:
        if vname == first_voice and current:
            phrases.append(current)
            current = {}
        current[vname] = line

    if current:
        phrases.append(current)

    return phrases


def phrase_chords(phrase_voices, downbeat_only=True):
    """Compute chord progression for a single phrase.

    Returns list of chord names (letter names, not Roman).
    """
    voice_names = list(phrase_voices.keys())
    if not voice_names:
        return []

    bass_voice = voice_names[-1]

    voice_measures = {}
    for vname, line in phrase_voices.items():
        voice_measures[vname] = extract_measures(line, downbeat_only=downbeat_only)

    max_m = max((len(m) for m in voice_measures.values()), default=0)

    chords = []
    for i in range(max_m):
        all_pcs = []
        bass_pcs = []
        for vname, measures in voice_measures.items():
            if i < len(measures):
                all_pcs.extend(measures[i])
                if vname == bass_voice:
                    bass_pcs.extend(measures[i])

        if not all_pcs:
            continue
        bass_pc = bass_pcs[0] if bass_pcs else None
        chord = identify_chord(all_pcs, bass_pc)
        if chord:
            chords.append(chord)

    return chords


def analyze_form(song):
    """Analyze a song's phrase structure and return form + per-section chords.

    Returns:
        form_str: e.g. "A A B" or "A B A B"
        sections: dict of label -> Roman numeral chord string
                  e.g. {'A': 'I - IV - V - I', 'B': 'IV - V - vi - I'}
    """
    phrases = group_into_phrases(song)
    if not phrases:
        return '', {}

    # Get chord tuple for each phrase
    phrase_chords_list = []
    for p in phrases:
        chords = phrase_chords(p)
        phrase_chords_list.append(tuple(chords))

    # Assign letters: identical chord tuples get the same letter
    label_map = {}  # chord_tuple -> letter
    labels = []
    next_letter = 0

    for ct in phrase_chords_list:
        if ct not in label_map:
            label_map[ct] = chr(ord('A') + next_letter)
            next_letter += 1
        labels.append(label_map[ct])

    form_str = ' '.join(labels)

    # Build sections dict: label -> chord progression as Roman numerals
    sections = {}
    for ct, label in label_map.items():
        roman = [chord_to_roman(c) for c in ct]
        # Show with bar lines between measures
        sections[label] = ' | '.join(roman)

    return form_str, sections


def detect_cycle(all_chords, min_len=3, threshold=0.60):
    """Detect a repeating cycle in a chord sequence.

    Tries every cycle length and starting offset.
    Returns (cycle_list, approx_repeats, match_rate) or None.
    """
    n = len(all_chords)
    if n < min_len * 2:
        return None

    best = None
    best_score = 0

    for k in range(min_len, n // 2 + 1):
        # Try each starting offset within the first period
        for offset in range(k):
            if offset + k > n:
                break
            cycle = all_chords[offset:offset + k]
            matches = sum(
                1 for i in range(n)
                if all_chords[i] == cycle[(i - offset) % k]
            )
            rate = matches / n

            if rate >= threshold:
                repeats = n / k
                # Prefer shorter cycles (more repeats) with high match rate
                score = rate * (repeats ** 0.5)
                if score > best_score:
                    best_score = score
                    best = (list(cycle), repeats, rate)

    return best


def format_song_analysis(song):
    """Return chord progression as Roman numerals with letter-labeled phrases.

    Each unique phrase gets a letter (A-Z). First occurrence shows chords:
        A(V V IV I) 2x B(IV V vi I)
    Subsequent identical phrases reuse the letter without chords.
    Cycle detection used when no exact phrase repeats exist (3+ phrases).
    """
    form_str, sections = analyze_form(song)
    if not form_str:
        return None

    labels = form_str.split()

    # Run-length encode consecutive repeated phrases
    runs = []
    i = 0
    while i < len(labels):
        label = labels[i]
        count = 1
        while i + count < len(labels) and labels[i + count] == label:
            count += 1
        runs.append((label, count))
        i += count

    # Check if any runs have repeats or any label appears more than once
    has_repeats = any(count > 1 for _, count in runs)
    has_reuse = len(set(labels)) < len(labels)

    # If no exact phrase repeats and 3+ unique phrases, try cycle detection
    if not has_repeats and not has_reuse and len(labels) >= 3:
        all_chords = []
        for label in labels:
            roman = sections[label].replace(' | ', ' ').split()
            all_chords.extend(roman)

        cycle = detect_cycle(all_chords)
        if cycle:
            cycle_chords, repeats, rate = cycle
            cycle_str = ' '.join(cycle_chords)
            repeat_count = round(repeats)
            if rate >= 0.90:
                return f"{cycle_str} {repeat_count}x"
            else:
                return f"{cycle_str} ~{repeat_count}x"

    # Build letter-labeled output
    defined = set()  # labels whose chords have been shown
    parts = []
    for label, count in runs:
        chords = sections[label].replace(' | ', ' ')
        if label not in defined:
            part = f"{label}({chords})"
            defined.add(label)
        else:
            part = label
        if count > 1:
            part += f" {count}x"
        parts.append(part)

    # If only one unique phrase, skip the letter label
    if len(set(labels)) == 1:
        chords = sections[labels[0]].replace(' | ', ' ')
        total = len(labels)
        if total > 1:
            return f"{chords} {total}x"
        else:
            return chords

    return ' '.join(parts)


def inject_chords(filename, output_filename):
    """Insert chord analysis as T: subtitle into the ABC file."""
    songs = parse_file(filename)

    chord_map = {}
    for song in songs:
        analysis = format_song_analysis(song)
        if analysis:
            chord_map[song['number']] = analysis

    with open(filename) as f:
        lines = f.readlines()

    out_lines = []
    current_song_num = None
    inserted = set()

    for line in lines:
        out_lines.append(line)
        stripped = line.strip()

        if stripped.startswith('X:'):
            current_song_num = stripped.split(':', 1)[1].strip()
        elif (stripped.startswith('T:')
              and current_song_num
              and current_song_num not in inserted
              and current_song_num in chord_map):
            out_lines.append(f'T: ({chord_map[current_song_num]})\n')
            inserted.add(current_song_num)

    with open(output_filename, 'w') as f:
        f.writelines(out_lines)

    print(f"Wrote {len(inserted)} chord analyses to {output_filename}")


def main():
    filename = sys.argv[1] if len(sys.argv) > 1 else 'OpenHymnal-C.abc'
    output = sys.argv[2] if len(sys.argv) > 2 else filename

    if '--print' in sys.argv:
        songs = parse_file(filename)
        for song in songs:
            analysis = format_song_analysis(song)
            if not analysis:
                print(f"#{song['number']:>3}  {song['title']:<50} (no data)")
                continue
            print(f"#{song['number']:>3}  {song['title']}")
            print(f"       {analysis}")
            print()
    else:
        inject_chords(filename, output)


if __name__ == '__main__':
    main()
