#!/usr/bin/env python3
"""
Precompute 6-8 note harp voicings for all hymn lead sheets.

For each chord change, generates optimal LH (3-4 notes) + RH (2-4 notes)
voicing within 10-string hand span per hand, 0-3 string gap between hands.
Outputs JSON for the HTML viewer.
"""

import json
import re
import sys
from collections import defaultdict

# ── Music constants ────────────────────────────────────────

MAJOR_SCALE = [0, 2, 4, 5, 7, 9, 11]
NOTE_TO_PC = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
PC_NAMES = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'F#', 'G', 'Ab', 'A', 'Bb', 'B']
ABC_BASE = {'C': 48, 'D': 50, 'E': 52, 'F': 53, 'G': 55, 'A': 57, 'B': 59}

CHORD_TEMPLATES = [
    (frozenset({0, 4, 7}),      '', False, ''),
    (frozenset({0, 3, 7}),      'm', True, ''),
    (frozenset({0, 4, 7, 10}),  '7', False, '7'),
    (frozenset({0, 3, 7, 10}),  'm7', True, '7'),
    (frozenset({0, 4, 7, 11}),  'maj7', False, 'M7'),
    (frozenset({0, 3, 6}),      'dim', True, '-'),
    (frozenset({0, 4, 8}),      'aug', False, '+'),
    (frozenset({0, 5, 7}),      'sus4', False, '~4'),
    (frozenset({0, 2, 7}),      'sus2', False, '~2'),
    (frozenset({0, 3, 6, 9}),   'dim7', True, '-7'),
    (frozenset({0, 3, 6, 10}),  'm7b5', True, '%7'),
]

INVERSION_ORDER = {
    '': [0, 4, 7], 'm': [0, 3, 7], '7': [0, 4, 7, 10],
    'm7': [0, 3, 7, 10], 'maj7': [0, 4, 7, 11], 'dim': [0, 3, 6],
    'aug': [0, 4, 8], 'sus4': [0, 5, 7], 'sus2': [0, 2, 7],
    'dim7': [0, 3, 6, 9], 'm7b5': [0, 3, 6, 10],
}

KEY_TO_PC = {'C':0,'Db':1,'D':2,'Eb':3,'E':4,'F':5,'F#':6,'Gb':6,
             'G':7,'Ab':8,'A':9,'Bb':10,'B':11,'Cb':11,'C#':1}

MAX_SPAN = 10
HARP_LOW = 36   # C2
HARP_HIGH = 96  # C7


def pc(midi):
    return midi % 12


def midi_to_string(midi, key_root):
    """Convert MIDI to diatonic string number, or None."""
    p = pc(midi)
    octave = midi // 12 - 1
    for deg, semi in enumerate(MAJOR_SCALE):
        if (key_root + semi) % 12 == p:
            return octave * 7 + deg
    return None


def string_to_midi(string_num, key_root):
    octave = string_num // 7
    deg = string_num % 7
    return (octave + 1) * 12 + (key_root + MAJOR_SCALE[deg]) % 12


def all_octaves(pitch_class, low, high):
    notes = []
    midi = pitch_class
    if midi < low:
        midi += ((low - midi + 11) // 12) * 12
    while midi <= high:
        notes.append(midi)
        midi += 12
    return notes


def abc_note_to_midi(token):
    acc = 0
    i = 0
    while i < len(token) and token[i] in '^_=':
        if token[i] == '^': acc += 1
        elif token[i] == '_': acc -= 1
        i += 1
    if i >= len(token):
        return None
    ch = token[i]
    is_lower = ch.islower()
    letter = ch.upper()
    if letter not in ABC_BASE:
        return None
    midi = ABC_BASE[letter] + acc
    if is_lower:
        midi += 12
    i += 1
    while i < len(token):
        if token[i] == ',': midi -= 12
        elif token[i] == "'": midi += 12
        else: break
        i += 1
    return midi


# ── Chord parsing ──────────────────────────────────────────

def parse_chord_symbol(s):
    """Parse chord symbol, return (root_degree, root_accidental, intervals, inversion)."""
    pos = 0
    root_acc = 0
    if s[pos:pos+1] == 'b': root_acc = -1; pos += 1
    elif s[pos:pos+1] == '#': root_acc = 1; pos += 1

    # Roman numeral
    upper = s[pos:].upper()
    degree = None
    for numeral, deg in [('VII',6),('VI',5),('IV',3),('V',4),('III',2),('II',1),('I',0)]:
        if upper.startswith(numeral):
            chunk = s[pos:pos+len(numeral)]
            is_minor = chunk == chunk.lower()
            degree = deg
            pos += len(numeral)
            break
    if degree is None:
        return None

    intervals = [0, 3, 7] if is_minor else [0, 4, 7]
    rest = s[pos:]

    if rest.startswith('+'):
        intervals = [0, 4, 8] if not is_minor else [0, 3, 8]
        pos += 1; rest = s[pos:]
    elif rest.startswith('-') and not rest.startswith('-7'):
        intervals = [0, 3, 6]
        pos += 1; rest = s[pos:]

    if rest.startswith('M7'):
        intervals.append(11); pos += 2
    elif rest.startswith('%7'):
        intervals = [0, 3, 6, 10]; pos += 2
    elif rest.startswith('-7'):
        intervals = [0, 3, 6, 9]; pos += 2
    elif rest.startswith('7'):
        intervals.append(10); pos += 1

    rest = s[pos:]
    if rest.startswith('~2'):
        intervals = [x for x in intervals if x not in (3,4)]
        intervals.append(2); intervals.sort(); pos += 2
    elif rest.startswith('~4'):
        intervals = [x for x in intervals if x not in (3,4)]
        intervals.append(5); intervals.sort(); pos += 2

    inversion = 0
    rest = s[pos:]
    if rest.startswith('^') and len(rest) > 1 and rest[1].isdigit():
        inversion = int(rest[1])

    return (degree, root_acc, intervals, inversion)


# ── Voicing engine ─────────────────────────────────────────

def generate_voicing(melody_midi, chord_str, key_root):
    """Generate optimal 6-8 note voicing for melody + chord.

    Returns dict with rh (list of midi), lh (list of midi), or None.
    """
    parsed = parse_chord_symbol(chord_str)
    if parsed is None:
        return None

    degree, root_acc, intervals, inversion = parsed
    root_pc = (key_root + MAJOR_SCALE[degree] + root_acc) % 12
    chord_pcs = [(root_pc + iv) % 12 for iv in intervals]

    # Required bass pc
    if inversion > 0 and inversion < len(intervals):
        bass_pc = (root_pc + intervals[inversion]) % 12
    else:
        bass_pc = root_pc

    melody_str = midi_to_string(melody_midi, key_root)
    if melody_str is None:
        return None

    # Generate all candidate notes (chord tones + root doublings)
    candidates = set()
    for cpc in chord_pcs:
        candidates.update(all_octaves(cpc, HARP_LOW, HARP_HIGH))
    candidates.update(all_octaves(root_pc, HARP_LOW, HARP_HIGH))
    # Filter to diatonic strings
    candidates = [m for m in sorted(candidates) if midi_to_string(m, key_root) is not None]

    # RH: melody thumb + up to 3 below within span
    rh_eligible = [m for m in candidates
                   if m < melody_midi
                   and midi_to_string(m, key_root) is not None
                   and melody_str - midi_to_string(m, key_root) < MAX_SPAN]
    rh_eligible.sort(reverse=True)  # closest to melody first
    rh_eligible = rh_eligible[:7]   # limit combos

    best_voicing = None
    best_score = -9999

    # Try RH with 1-3 additional notes
    rh_configs = [[melody_midi]]
    for i, a in enumerate(rh_eligible):
        rh_configs.append([melody_midi, a])
        for j in range(i+1, len(rh_eligible)):
            rh_configs.append([melody_midi, a, rh_eligible[j]])
            for k in range(j+1, len(rh_eligible)):
                lo_s = midi_to_string(rh_eligible[k], key_root)
                if lo_s is not None and melody_str - lo_s < MAX_SPAN:
                    rh_configs.append([melody_midi, a, rh_eligible[j], rh_eligible[k]])

    for rh in rh_configs:
        rh_lowest = min(rh)
        rh_lowest_str = midi_to_string(rh_lowest, key_root)
        if rh_lowest_str is None:
            continue

        # LH must be below RH lowest
        lh_upper_str = rh_lowest_str - 1
        if lh_upper_str < 0:
            continue
        lh_upper_midi = string_to_midi(lh_upper_str, key_root)

        lh_eligible = [m for m in candidates if m <= lh_upper_midi and m >= HARP_LOW]
        lh_eligible.sort(reverse=True)
        lh_eligible = lh_eligible[:7]

        # Try LH with 1-4 notes
        lh_configs = []
        for ti, thumb in enumerate(lh_eligible):
            ts = midi_to_string(thumb, key_root)
            if ts is None: continue
            below = [m for m in lh_eligible[ti+1:]
                     if midi_to_string(m, key_root) is not None
                     and ts - midi_to_string(m, key_root) < MAX_SPAN]

            # Check bass constraint
            def bass_ok(hand):
                if not hand: return False
                return pc(min(hand)) == bass_pc

            configs_for_thumb = [[thumb]]
            for i, a in enumerate(below):
                configs_for_thumb.append([thumb, a])
                for j in range(i+1, len(below)):
                    configs_for_thumb.append([thumb, a, below[j]])
                    for k in range(j+1, len(below)):
                        configs_for_thumb.append([thumb, a, below[j], below[k]])

            for lh in configs_for_thumb:
                if not bass_ok(lh):
                    continue
                lh_configs.append(lh)

        if not lh_configs:
            # Try without bass constraint
            for ti, thumb in enumerate(lh_eligible):
                ts = midi_to_string(thumb, key_root)
                if ts is None: continue
                below = [m for m in lh_eligible[ti+1:]
                         if midi_to_string(m, key_root) is not None
                         and ts - midi_to_string(m, key_root) < MAX_SPAN]
                lh_configs.append([thumb])
                for i, a in enumerate(below):
                    lh_configs.append([thumb, a])
                    for j in range(i+1, len(below)):
                        lh_configs.append([thumb, a, below[j]])
                        for k in range(j+1, len(below)):
                            lh_configs.append([thumb, a, below[j], below[k]])

        for lh in lh_configs:
            total = len(rh) + len(lh)
            if total < 5:
                continue  # want 6-8 notes

            score = score_voicing(rh, lh, root_pc, chord_pcs, bass_pc)
            if score > best_score:
                best_score = score
                best_voicing = {'rh': sorted(rh, reverse=True),
                                'lh': sorted(lh, reverse=True)}

    return best_voicing


def score_voicing(rh, lh, root_pc, chord_pcs, bass_pc):
    all_midis = sorted(rh + lh)
    if not all_midis:
        return -9999

    score = 0.0

    # Total span
    score += (all_midis[-1] - all_midis[0]) * 2.0

    # Bass is correct
    if pc(all_midis[0]) == bass_pc:
        score += 20

    # Root in bass register
    for m in all_midis:
        if m < 60 and pc(m) == root_pc:
            score += 10

    # Chord tone coverage
    covered = set(pc(m) for m in all_midis)
    score += sum(15 for cpc in chord_pcs if cpc in covered)

    # Total notes (prefer 6-8)
    total = len(all_midis)
    score += total * 8
    if 6 <= total <= 8:
        score += 20

    # Gaps
    for i in range(len(all_midis) - 1):
        gap = all_midis[i+1] - all_midis[i]
        if gap > 12:
            score -= (gap - 12) * 3

    # Both hands
    if lh and rh:
        score += 15

    # LH bias: prefer LH 3-4, RH 2-3
    if len(lh) >= 3: score += 20
    if len(lh) == 4: score += 10
    if len(rh) == 2: score += 10
    if len(rh) == 3: score += 5
    if len(rh) >= 4 and len(lh) < 3: score -= 15

    return score


# ── Chord identification ──────────────────────────────────

def identify_chord_from_notes(midis, key_root):
    """Identify chord symbol from a list of MIDI notes."""
    pcs = set(pc(m) for m in midis)
    if len(pcs) < 2:
        return ''

    bass = pc(min(midis))
    best_root = 0
    best_quality = ''
    best_score = -999

    for root in range(12):
        ivs = frozenset((p - root) % 12 for p in pcs)
        for tmpl, quality, is_lower, ascii_suffix in CHORD_TEMPLATES:
            matched = len(ivs & tmpl)
            extra = len(ivs - tmpl)
            missing = len(tmpl - ivs)
            s = matched * 3 - extra * 2 - missing
            if root == bass: s += 4
            if len(tmpl) == 3: s += 2
            if s > best_score:
                best_score = s
                best_root = root
                best_quality = quality

    # Find degree
    deg = None
    for d, semi in enumerate(MAJOR_SCALE):
        if (key_root + semi) % 12 == best_root:
            deg = d
            break
    if deg is None:
        return ''

    for tmpl, quality, is_lower, ascii_suffix in CHORD_TEMPLATES:
        if quality == best_quality:
            roman = ['I','II','III','IV','V','VI','VII'][deg]
            if is_lower:
                roman = roman.lower()
            # Inversion
            inv_order = INVERSION_ORDER.get(quality, [])
            bass_iv = (bass - best_root) % 12
            inv = 0
            if bass_iv in inv_order:
                inv = inv_order.index(bass_iv)
            result = roman + ascii_suffix
            if inv > 0:
                result += f'^{inv}'
            return result

    return ''


# ── Lead sheet parsing ─────────────────────────────────────

def parse_leadsheets(filename):
    with open(filename, encoding='latin-1') as f:
        text = f.read()

    hymns = []
    for tune in text.split('\nX:'):
        if not tune.strip():
            continue
        if not tune.startswith('X:'):
            tune = 'X:' + tune

        lines = tune.split('\n')
        number = title = meter = key = ''
        def_len = '1/4'
        tempo = 120
        music_lines = []

        for line in lines:
            if line.startswith('X:'): number = line[2:].strip()
            elif line.startswith('T:'): title = title or line[2:].strip()
            elif line.startswith('M:'): meter = line[2:].strip().split('%')[0].strip()
            elif line.startswith('L:'): def_len = line[2:].strip().split('%')[0].strip()
            elif line.startswith('K:'): key = line[2:].strip().split('%')[0].strip()
            elif line.strip() and not line.startswith('%'):
                music_lines.append(line)

        music_str = ' '.join(music_lines)

        # Extract tempo
        tm = re.search(r'\[Q:\s*\d+/\d+=(\d+)\]', music_str)
        if tm:
            tempo = int(tm.group(1))

        # Parse beats
        beats = []
        last_chord = None
        for m in re.finditer(r'("([^"]*)")?\s*([_^=]*[A-Ga-g][,\']*[\d/]*|z[\d/]*|\[x\d*\]|\|[\|\]\[:]*)', music_str):
            chord_sym = m.group(2) or None
            token = m.group(3)
            if chord_sym:
                last_chord = chord_sym
            if token.startswith('|'):
                beats.append({'type': 'bar'})
                continue
            if token.startswith('z') or token.startswith('[x'):
                continue
            pm = re.match(r'^[_^=]*[A-Ga-g][,\']*', token)
            if not pm:
                continue
            midi = abc_note_to_midi(pm.group())
            if midi is None:
                continue
            beats.append({
                'type': 'note',
                'midi': midi,
                'chord': chord_sym,
                'active_chord': last_chord,
            })

        if beats:
            hymns.append({
                'number': number, 'title': title, 'key': key,
                'meter': meter, 'tempo': tempo, 'beats': beats,
            })

    return hymns


# ── Main: generate voicings ───────────────────────────────

def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'harp_leadsheets.abc'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'harp_voicings.json'

    hymns = parse_leadsheets(input_file)
    print(f"Parsed {len(hymns)} hymns")

    results = []
    total_chords = 0
    total_voiced = 0

    for hymn in hymns:
        key_root = KEY_TO_PC.get(hymn['key'], 0)
        key_str = midi_to_string(key_root + 48, key_root)  # reference C for string numbering

        voiced_beats = []
        for beat in hymn['beats']:
            if beat['type'] == 'bar':
                voiced_beats.append({'type': 'bar'})
                continue

            entry = {
                'type': 'note',
                'midi': beat['midi'],
            }

            if beat.get('chord'):
                total_chords += 1
                voicing = generate_voicing(beat['midi'], beat['chord'], key_root)
                if voicing:
                    total_voiced += 1
                    # Compute string numbers relative to key root
                    ref_str = midi_to_string(key_root + 36, key_root) or 0  # lowest key root
                    entry['rh'] = [midi_to_string(m, key_root) - ref_str + 1
                                   for m in voicing['rh']
                                   if midi_to_string(m, key_root) is not None]
                    entry['lh'] = [midi_to_string(m, key_root) - ref_str + 1
                                   for m in voicing['lh']
                                   if midi_to_string(m, key_root) is not None]
                    entry['rh_midi'] = voicing['rh']
                    entry['lh_midi'] = voicing['lh']
                    # Identify LH and RH chord shapes
                    entry['rh_chord'] = identify_chord_from_notes(voicing['rh'], key_root)
                    entry['lh_chord'] = identify_chord_from_notes(voicing['lh'], key_root)
                    entry['chord'] = beat['chord']
                    entry['total_notes'] = len(voicing['rh']) + len(voicing['lh'])

            voiced_beats.append(entry)

        results.append({
            'number': hymn['number'],
            'title': hymn['title'],
            'key': hymn['key'],
            'meter': hymn['meter'],
            'tempo': hymn['tempo'],
            'beats': voiced_beats,
        })

    with open(output_file, 'w') as f:
        json.dump(results, f)

    print(f"Voiced {total_voiced}/{total_chords} chord changes")
    print(f"Wrote {output_file}")


if __name__ == '__main__':
    main()
