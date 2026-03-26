#!/usr/bin/env python3
"""Convert harp_voicings.json to multi-voice ABC for PDF rendering.

Generates a 3-staff ABC file: Melody (treble), RH chords (treble), LH chords (bass).
All notes are diatonic (no accidentals) — uses harp string numbers to determine
note letters, guaranteeing correct key-signature spelling.

Render with:
  python3 voicings_to_abc.py
  abcm2ps harp_voicings_score.abc -O - | ps2pdf - harp_voicings_score.pdf
"""

import json
import sys
from fractions import Fraction

MAJOR_SCALE = [0, 2, 4, 5, 7, 9, 11]
NOTE_NAMES = ['C', 'D', 'E', 'F', 'G', 'A', 'B']
NATURAL_PCS = [0, 2, 4, 5, 7, 9, 11]

KEY_ROOT_PC = {
    'C': 0, 'Db': 1, 'D': 2, 'Eb': 3, 'E': 4, 'F': 5,
    'F#': 6, 'Gb': 6, 'G': 7, 'Ab': 8, 'A': 9, 'Bb': 10, 'B': 11,
}

# Key root PC -> note letter index (0=C..6=B)
KEY_ROOT_LETTER = {
    0: 0, 1: 1, 2: 1, 3: 2, 4: 2, 5: 3,
    6: 3, 7: 4, 8: 5, 9: 5, 10: 6, 11: 6,
}

HARP_LOW = 36   # C2
HARP_HIGH = 96  # C7


def midi_to_raw_string(midi, key_root_pc):
    """Convert MIDI to raw (absolute) harp string number, or None if not diatonic."""
    p = midi % 12
    octave = midi // 12 - 1
    for deg, semi in enumerate(MAJOR_SCALE):
        if (key_root_pc + semi) % 12 == p:
            return octave * 7 + deg
    return None


def raw_string_to_midi(string_num, key_root_pc):
    """Convert raw string number to MIDI note."""
    octave = string_num // 7
    deg = string_num % 7
    return (octave + 1) * 12 + (key_root_pc + MAJOR_SCALE[deg]) % 12


def compute_ref_str(key_root_pc):
    """Find the raw string number of the lowest root on the harp."""
    root_midis = [m for m in range(HARP_LOW, HARP_HIGH + 1) if m % 12 == key_root_pc]
    if root_midis:
        return midi_to_raw_string(root_midis[0], key_root_pc)
    return 0


def relative_to_raw(rel_str, ref_str):
    """Convert 1-indexed relative string to raw string."""
    return rel_str + ref_str - 1


def raw_string_to_abc(raw_str, key_root_pc):
    """Convert a raw harp string number to ABC note. Never produces accidentals."""
    deg = raw_str % 7
    root_letter = KEY_ROOT_LETTER[key_root_pc]
    letter_idx = (root_letter + deg) % 7
    letter = NOTE_NAMES[letter_idx]

    midi = raw_string_to_midi(raw_str, key_root_pc)
    midi_octave = midi // 12 - 1

    if midi_octave <= 4:
        abc_letter = letter.upper()
        suffix = ',' * (4 - midi_octave)
    else:
        abc_letter = letter.lower()
        suffix = "'" * (midi_octave - 5)

    return abc_letter + suffix


def midi_to_abc_diatonic(midi, key_root_pc):
    """Convert MIDI to ABC, snapping to nearest diatonic note if needed."""
    raw = midi_to_raw_string(midi, key_root_pc)
    if raw is not None:
        return raw_string_to_abc(raw, key_root_pc)

    # Snap to nearest diatonic note
    pc = midi % 12
    best_dist = 99
    best_deg = 0
    for deg, semi in enumerate(MAJOR_SCALE):
        scale_pc = (key_root_pc + semi) % 12
        dist = min(abs(pc - scale_pc), 12 - abs(pc - scale_pc))
        if dist < best_dist:
            best_dist = dist
            best_deg = deg

    octave = midi // 12 - 1
    raw = octave * 7 + best_deg
    return raw_string_to_abc(raw, key_root_pc)


def duration_to_abc(dur_tuple, default_len=(1, 4)):
    """Convert duration (num, den) in quarter-note beats to ABC duration string."""
    if dur_tuple is None:
        return ''

    dur = Fraction(dur_tuple[0], dur_tuple[1])
    dl = Fraction(default_len[0], default_len[1]) / Fraction(1, 4)
    ratio = dur / dl

    if ratio == 1:
        return ''
    elif ratio.denominator == 1:
        return str(ratio.numerator)
    elif ratio.numerator == 1:
        return f'/{ratio.denominator}'
    else:
        return f'{ratio.numerator}/{ratio.denominator}'


def beats_to_abc_voice(beats, key_root_pc, ref_str, default_len, voice='melody'):
    """Convert beats to an ABC voice string."""
    parts = []
    measure_notes = []

    for beat in beats:
        bt = beat.get('type', '')

        if bt == 'bar':
            if measure_notes:
                parts.append(' '.join(measure_notes))
                measure_notes = []
            continue

        if bt == 'rest':
            dur_str = duration_to_abc(beat.get('duration'), default_len)
            measure_notes.append(f'z{dur_str}')
            continue

        if bt != 'note':
            continue

        dur = beat.get('duration')
        dur_str = duration_to_abc(dur, default_len)

        if voice == 'rh':
            # RH includes melody (thumb = top note) + chord notes
            rel_strings = beat.get('rh')
            chord_sym = beat.get('chord', '')
            prefix = f'"{chord_sym}"' if chord_sym else ''
            if not rel_strings:
                # No voicing change — just the melody note
                midi = beat.get('midi')
                if midi is None:
                    measure_notes.append(f'{prefix}z{dur_str}')
                else:
                    note = midi_to_abc_diatonic(midi, key_root_pc)
                    measure_notes.append(f'{prefix}{note}{dur_str}')
            else:
                raw_strings = [relative_to_raw(s, ref_str) for s in rel_strings]
                notes = [raw_string_to_abc(s, key_root_pc) for s in sorted(raw_strings)]
                if len(notes) == 1:
                    measure_notes.append(f'{prefix}{notes[0]}{dur_str}')
                else:
                    measure_notes.append(f'{prefix}[{"".join(notes)}]{dur_str}')

        elif voice == 'lh':
            rel_strings = beat.get('lh')
            if not rel_strings:
                measure_notes.append(f'z{dur_str}')
            else:
                raw_strings = [relative_to_raw(s, ref_str) for s in rel_strings]
                notes = [raw_string_to_abc(s, key_root_pc) for s in sorted(raw_strings)]
                if len(notes) == 1:
                    measure_notes.append(f'{notes[0]}{dur_str}')
                else:
                    measure_notes.append(f'[{"".join(notes)}]{dur_str}')

    if measure_notes:
        parts.append(' '.join(measure_notes))

    return ' | '.join(parts) + ' |]'


def hymn_to_abc(hymn, idx):
    """Convert a hymn dict to a multi-voice ABC tune."""
    key = hymn['key']
    key_root_pc = KEY_ROOT_PC.get(key, 0)
    ref_str = compute_ref_str(key_root_pc)
    meter = hymn['meter']
    tempo = hymn.get('tempo', 100)
    default_len = (1, 4)

    rh = beats_to_abc_voice(hymn['beats'], key_root_pc, ref_str, default_len, 'rh')
    lh = beats_to_abc_voice(hymn['beats'], key_root_pc, ref_str, default_len, 'lh')

    lines = [
        f'X:{idx}',
        f'T:{hymn["title"]}',
        f'M:{meter}',
        f'L:1/4',
        f'Q:1/4={tempo}',
        f'%%staves {{1 2}}',
        f'V:1 clef=treble name="RH"',
        f'V:2 clef=bass name="LH"',
        f'K:{key}',
        f'[V:1] {rh}',
        f'[V:2] {lh}',
    ]
    return '\n'.join(lines)


def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'harp_voicings.json'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'harp_voicings_score.abc'

    with open(input_file) as f:
        hymns = json.load(f)

    print(f"Converting {len(hymns)} hymns...")

    abc_tunes = []
    for i, hymn in enumerate(hymns, 1):
        abc_tunes.append(hymn_to_abc(hymn, i))

    with open(output_file, 'w') as f:
        f.write('\n\n'.join(abc_tunes) + '\n')

    print(f"Wrote {len(abc_tunes)} tunes to {output_file}")


if __name__ == '__main__':
    main()
