#!/usr/bin/env python3
"""
Generate harp lead sheets from SATB hymns in OpenHymnal-C.abc.

Parses 4-part SATB harmony, derives per-beat chords using Roman numeral
analysis with 7-bit ASCII notation, and generates LaTeX + ABC output.

ASCII chord notation:
  Case    = major/minor 3rd (uppercase=major, lowercase=minor)
  +       = augmented (raised 5th)
  -       = diminished (lowered 5th)
  M7      = major 7th
  7       = dominant 7th
  %7      = half-diminished 7th
  -7      = diminished 7th
  ~2/~4   = suspended 2nd/4th
  ^1/^2/^3 = inversions
"""

import re
import sys
from fractions import Fraction
from pathlib import Path

# ── Pitch class data ──────────────────────────────────────────────

NOTE_TO_PC = {'C': 0, 'D': 2, 'E': 4, 'F': 5, 'G': 7, 'A': 9, 'B': 11}
PC_TO_NAME = ['C', 'Db', 'D', 'Eb', 'E', 'F', 'F#', 'G', 'Ab', 'A', 'Bb', 'B']

PC_TO_ROMAN = {
    0: 'I', 1: 'bII', 2: 'II', 3: 'bIII', 4: 'III', 5: 'IV',
    6: '#IV', 7: 'V', 8: 'bVI', 9: 'VI', 10: 'bVII', 11: 'VII',
}

# Chord templates: intervals from root -> quality key
CHORD_TYPES = [
    (frozenset({0, 4, 7}),      ''),
    (frozenset({0, 3, 7}),      'm'),
    (frozenset({0, 4, 7, 10}),  '7'),
    (frozenset({0, 3, 7, 10}),  'm7'),
    (frozenset({0, 4, 7, 11}),  'maj7'),
    (frozenset({0, 3, 6}),      'dim'),
    (frozenset({0, 4, 8}),      'aug'),
    (frozenset({0, 5, 7}),      'sus4'),
    (frozenset({0, 2, 7}),      'sus2'),
    (frozenset({0, 3, 6, 9}),   'dim7'),
    (frozenset({0, 3, 6, 10}),  'm7b5'),
]

# quality key -> (lowercase roman?, ASCII suffix)
QUALITY_TO_ASCII = {
    '':      (False, ''),
    'm':     (True,  ''),
    '7':     (False, '7'),
    'm7':    (True,  '7'),
    'maj7':  (False, 'M7'),
    'dim':   (True,  '-'),
    'aug':   (False, '+'),
    'sus4':  (False, '~4'),
    'sus2':  (False, '~2'),
    'dim7':  (True,  '-7'),
    'm7b5':  (True,  '%7'),
}

# Interval index within a chord template for computing inversions
# e.g., for major triad {0,4,7}: 0=root, 4=3rd, 7=5th
INVERSION_ORDER = {
    '':      [0, 4, 7],
    'm':     [0, 3, 7],
    '7':     [0, 4, 7, 10],
    'm7':    [0, 3, 7, 10],
    'maj7':  [0, 4, 7, 11],
    'dim':   [0, 3, 6],
    'aug':   [0, 4, 8],
    'sus4':  [0, 5, 7],
    'sus2':  [0, 2, 7],
    'dim7':  [0, 3, 6, 9],
    'm7b5':  [0, 3, 6, 10],
}


# ── ABC note parsing ──────────────────────────────────────────────

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


def parse_abc_duration(suffix, default_len=Fraction(1, 4)):
    """Parse ABC duration suffix into beats (as Fraction).

    default_len is the L: value (e.g., Fraction(1,4) for L:1/4).
    Duration is returned in quarter-note beats.
    """
    beats_per_default = default_len / Fraction(1, 4)

    if not suffix or suffix.isspace():
        return beats_per_default

    # Count leading slashes for shorthand: / = /2, // = /4
    if suffix.startswith('/'):
        slash_only = suffix.rstrip('0123456789')
        num_part = suffix[len(slash_only):]
        if not num_part:
            divisor = 2 ** len(slash_only)
        else:
            divisor = int(num_part)
        return beats_per_default * Fraction(1, divisor)

    # n/m or just n
    if '/' in suffix:
        parts = suffix.split('/')
        return beats_per_default * Fraction(int(parts[0]), int(parts[1]))

    return beats_per_default * int(suffix)


# Note token regex: optional accidentals, letter, optional octave marks, optional duration
NOTE_RE = re.compile(
    r'([_^=]*[A-Ga-g][,\']*)'   # pitch part (group 1)
    r'(\d+(?:/\d+)?|/+\d*)?'    # duration part (group 2)
)

REST_RE = re.compile(
    r'(z)'                       # rest (group 1)
    r'(\d+(?:/\d+)?|/+\d*)?'    # duration part (group 2)
)

CHORD_BRACKET_RE = re.compile(
    r'\[([_^=]*[A-Ga-g](?:[_^=]*[A-Ga-g,\']*)*)\]'  # chord bracket
    r'(\d+(?:/\d+)?|/+\d*)?'    # duration part
)


def extract_timed_notes(music_str, default_len=Fraction(1, 4)):
    """Parse a voice's music string into timed events per measure.

    Returns list of measures, each measure is a list of
    (beat_offset, pitch_class_or_None, duration_in_quarter_beats).
    """
    # Strip inline fields, decorations, grace notes
    s = re.sub(r'\[V:\s*[^\]]+\]', '', music_str)
    s = re.sub(r'\[[A-Z]:[^\]]*\]', '', s)
    s = re.sub(r'![a-zA-Z0-9_]+!', '', s)
    s = re.sub(r'\{[^}]*\}', '', s)  # grace notes
    s = re.sub(r'"[^"]*"', '', s)     # chord annotations

    measures_str = s.split('|')
    all_measures = []

    for mstr in measures_str:
        events = []
        beat_offset = Fraction(0)
        pos = 0
        while pos < len(mstr):
            # Skip whitespace, ties, slurs, dots (staccato)
            if mstr[pos] in ' \t\n-().~':
                pos += 1
                continue

            # Chord bracket [CEG]
            cm = CHORD_BRACKET_RE.match(mstr, pos)
            if cm:
                inner = cm.group(1)
                dur_str = cm.group(2) or ''
                dur = parse_abc_duration(dur_str, default_len)
                # Extract first note's pitch class (for soprano this is what matters)
                first_note = re.match(r'[_^=]*[A-Ga-g][,\']*', inner)
                if first_note:
                    pc = abc_note_to_pc(first_note.group())
                    events.append((beat_offset, pc, dur))
                beat_offset += dur
                pos = cm.end()
                continue

            # Rest
            rm = REST_RE.match(mstr, pos)
            if rm:
                dur_str = rm.group(2) or ''
                dur = parse_abc_duration(dur_str, default_len)
                events.append((beat_offset, None, dur))
                beat_offset += dur
                pos = rm.end()
                continue

            # Note
            nm = NOTE_RE.match(mstr, pos)
            if nm:
                pitch_str = nm.group(1)
                dur_str = nm.group(2) or ''
                dur = parse_abc_duration(dur_str, default_len)
                pc = abc_note_to_pc(pitch_str)
                events.append((beat_offset, pc, dur))
                beat_offset += dur
                pos = nm.end()
                continue

            # Skip unrecognized character
            pos += 1

        all_measures.append(events)

    return all_measures


# ── Chord identification ─────────────────────────────────────────

def identify_chord_ext(pitch_classes, bass_pc=None):
    """Identify chord, returning (root_pc, quality, inversion).

    Returns (None, None, 0) if no chord can be identified.
    """
    pcs = set(pitch_classes)
    if len(pcs) < 2:
        if pcs:
            return (list(pcs)[0], '', 0)
        return (None, None, 0)

    best_root = None
    best_quality = None
    best_score = -999

    for root in range(12):
        intervals = frozenset((pc - root) % 12 for pc in pcs)
        for idx, (template, quality) in enumerate(CHORD_TYPES):
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
                best_root = root
                best_quality = quality

    if best_root is None:
        return (None, None, 0)

    # Compute inversion from bass note
    inversion = 0
    if bass_pc is not None and best_quality in INVERSION_ORDER:
        intervals = INVERSION_ORDER[best_quality]
        bass_interval = (bass_pc - best_root) % 12
        if bass_interval in intervals:
            inv_idx = intervals.index(bass_interval)
            if inv_idx > 0:
                inversion = inv_idx

    return (best_root, best_quality, inversion)


def chord_to_ascii_roman(root_pc, quality, inversion=0):
    """Convert chord to 7-bit ASCII Roman numeral notation."""
    if root_pc is None:
        return '?'

    roman = PC_TO_ROMAN.get(root_pc, '?')
    lower, suffix = QUALITY_TO_ASCII.get(quality, (False, quality or ''))
    if lower:
        roman = roman.lower()

    result = roman + suffix
    if inversion > 0:
        result += f'^{inversion}'
    return result


# ── Per-beat chord analysis ──────────────────────────────────────

def parse_time_sig(m_field):
    """Parse M: field into (numerator, denominator) or None for free meter."""
    m_field = m_field.strip()
    # Strip comments
    if '%' in m_field:
        m_field = m_field[:m_field.index('%')].strip()
    if m_field == 'none' or m_field == 'C|':
        if m_field == 'C|':
            return (2, 2)
        return None
    if m_field == 'C':
        return (4, 4)
    m = re.match(r'(\d+)/(\d+)', m_field)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    return None


def parse_default_len(l_field):
    """Parse L: field into a Fraction."""
    l_field = l_field.strip()
    if '%' in l_field:
        l_field = l_field[:l_field.index('%')].strip()
    if '\t' in l_field:
        l_field = l_field[:l_field.index('\t')].strip()
    m = re.match(r'(\d+)/(\d+)', l_field)
    if m:
        return Fraction(int(m.group(1)), int(m.group(2)))
    return Fraction(1, 4)


def beat_chords_for_phrase(phrase_voices, default_len=Fraction(1, 4),
                           beat_unit=Fraction(1)):
    """Compute per-beat chords from all voices in a phrase.

    beat_unit: duration of one beat in quarter-note units (1 for quarter, 1.5 for dotted quarter).
    Returns list of (measure_idx, beat_offset, ascii_chord_str).
    Only emits when chord changes from previous beat.
    """
    voice_names = list(phrase_voices.keys())
    if not voice_names:
        return []

    bass_voice = voice_names[-1]

    # Parse all voices into timed events
    voice_measures = {}
    for vname, line in phrase_voices.items():
        voice_measures[vname] = extract_timed_notes(line, default_len)

    max_m = max((len(m) for m in voice_measures.values()), default=0)

    results = []
    prev_chord = None

    for mi in range(max_m):
        # Find all beat positions in this measure across all voices
        beat_positions = set()
        for vname in voice_names:
            measures = voice_measures[vname]
            if mi < len(measures):
                for offset, pc, dur in measures[mi]:
                    # Quantize to beat grid
                    quantized = (offset / beat_unit).limit_denominator(8) * beat_unit
                    beat_positions.add(quantized)

        for beat_pos in sorted(beat_positions):
            all_pcs = []
            bass_pc = None

            for vname in voice_names:
                measures = voice_measures[vname]
                if mi >= len(measures):
                    continue
                # Find the note sounding at this beat position
                sounding_pc = None
                for offset, pc, dur in measures[mi]:
                    if offset <= beat_pos < offset + dur:
                        sounding_pc = pc
                        break
                if sounding_pc is not None:
                    all_pcs.append(sounding_pc)
                    if vname == bass_voice:
                        bass_pc = sounding_pc

            if len(all_pcs) < 2:
                continue

            root_pc, quality, inversion = identify_chord_ext(all_pcs, bass_pc)
            if root_pc is None:
                continue

            chord_str = chord_to_ascii_roman(root_pc, quality, inversion)

            if chord_str != prev_chord:
                results.append((mi, beat_pos, chord_str))
                prev_chord = chord_str

    return results


# ── Soprano melody annotation ────────────────────────────────────

def annotate_soprano(soprano_line, beat_chords, default_len=Fraction(1, 4)):
    """Inject chord annotations into soprano melody ABC string.

    beat_chords: list of (measure_idx, beat_offset, chord_str)
    Returns annotated ABC string.
    """
    # Strip [V:...] prefix
    line = re.sub(r'\[V:\s*[^\]]+\]', '', soprano_line)
    # Strip [Q:...] tempo
    line = re.sub(r'\[Q:[^\]]*\]', '', line)
    # Strip lyrics-related stuff and decorations
    line = re.sub(r'![a-zA-Z0-9_]+!', '', line)

    # Build a lookup: (measure_idx, beat_offset) -> chord_str
    chord_at = {}
    for mi, beat_off, chord_str in beat_chords:
        chord_at[(mi, beat_off)] = chord_str

    # Walk through the melody, splitting by measures
    measures = line.split('|')
    annotated_measures = []

    for mi, mstr in enumerate(measures):
        beat_offset = Fraction(0)
        out = ''
        pos = 0

        while pos < len(mstr):
            # Skip whitespace
            if mstr[pos] in ' \t':
                out += mstr[pos]
                pos += 1
                continue

            # Skip ties, slurs, staccato
            if mstr[pos] in '-().~':
                out += mstr[pos]
                pos += 1
                continue

            # Grace notes - pass through
            gm = re.match(r'\{[^}]*\}', mstr[pos:])
            if gm:
                out += gm.group()
                pos += gm.end()
                continue

            # Inline fields - pass through
            fm = re.match(r'\[[A-Z]:[^\]]*\]', mstr[pos:])
            if fm:
                out += fm.group()
                pos += fm.end()
                continue

            # Chord annotation already present - strip it
            am = re.match(r'"[^"]*"', mstr[pos:])
            if am:
                pos += am.end()
                continue

            # Chord bracket
            cm = CHORD_BRACKET_RE.match(mstr, pos)
            if cm:
                dur_str = cm.group(2) or ''
                dur = parse_abc_duration(dur_str, default_len)
                chord = chord_at.get((mi, beat_offset))
                if chord:
                    out += f'"{chord}"'
                out += cm.group()
                beat_offset += dur
                pos = cm.end()
                continue

            # Rest
            rm = REST_RE.match(mstr, pos)
            if rm:
                dur_str = rm.group(2) or ''
                dur = parse_abc_duration(dur_str, default_len)
                chord = chord_at.get((mi, beat_offset))
                if chord:
                    out += f'"{chord}"'
                out += rm.group()
                beat_offset += dur
                pos = rm.end()
                continue

            # Note
            nm = NOTE_RE.match(mstr, pos)
            if nm:
                pitch_str = nm.group(1)
                dur_str = nm.group(2) or ''
                dur = parse_abc_duration(dur_str, default_len)
                chord = chord_at.get((mi, beat_offset))
                if chord:
                    out += f'"{chord}"'
                out += nm.group()
                beat_offset += dur
                pos = nm.end()
                continue

            # Pass through unrecognized
            out += mstr[pos]
            pos += 1

        annotated_measures.append(out)

    return '|'.join(annotated_measures)


# ── ABC file parsing ─────────────────────────────────────────────

def parse_file(filename):
    """Parse ABC file into list of song dicts."""
    with open(filename, encoding='latin-1') as f:
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
                'titles': [],
                'key': 'C',
                'meter': '4/4',
                'default_len': '1/4',
                'voice_lines': [],
                'first_voice': None,
            }
        elif song is None:
            continue
        elif line.startswith('T:'):
            t = line.split(':', 1)[1].strip()
            song['titles'].append(t)
            if not song['title']:
                song['title'] = t
        elif line.startswith('K:'):
            song['key'] = line.split(':', 1)[1].strip()
        elif line.startswith('M:'):
            song['meter'] = line.split(':', 1)[1].strip()
        elif line.startswith('L:'):
            song['default_len'] = line.split(':', 1)[1].strip()
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
    """Split voice lines into phrases (new phrase when first voice reappears)."""
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


# ── Lead sheet generation ────────────────────────────────────────

def generate_leadsheet_abc(song):
    """Generate a single-voice ABC lead sheet with chord annotations for a song."""
    phrases = group_into_phrases(song)
    if not phrases:
        return None

    default_len = parse_default_len(song['default_len'])
    time_sig = parse_time_sig(song['meter'])

    # Determine beat unit in quarter-note beats
    if time_sig is None:
        beat_unit = Fraction(1)  # free meter: use quarter note
    elif time_sig[0] in (6, 9, 12) and time_sig[1] == 8:
        beat_unit = Fraction(3, 2)  # compound time: dotted quarter
    else:
        beat_unit = Fraction(1)  # simple time: quarter note

    soprano_voice = song['first_voice']
    melody_lines = []

    for phrase in phrases:
        # Get per-beat chords
        chords = beat_chords_for_phrase(phrase, default_len, beat_unit)

        # Get soprano line and annotate
        soprano_line = phrase.get(soprano_voice, '')
        if not soprano_line:
            continue

        annotated = annotate_soprano(soprano_line, chords, default_len)
        # Clean up: remove [V:...] prefix if still present
        annotated = re.sub(r'\[V:\s*[^\]]+\]', '', annotated).strip()
        melody_lines.append(annotated)

    if not melody_lines:
        return None

    # Build clean meter string
    meter = song['meter']
    if '%' in meter:
        meter = meter[:meter.index('%')].strip()
    if '\t' in meter:
        meter = meter[:meter.index('\t')].strip()

    dl = song['default_len']
    if '%' in dl:
        dl = dl[:dl.index('%')].strip()
    if '\t' in dl:
        dl = dl[:dl.index('\t')].strip()

    # Build ABC block
    abc_lines = [
        f'X:{song["number"]}',
        f'T:{song["title"]}',
        f'M:{meter}',
        f'L:{dl}',
        f'K:{song["key"].split("%")[0].split("\t")[0].strip()}',
    ]
    for ml in melody_lines:
        abc_lines.append(ml)

    return '\n'.join(abc_lines)


def generate_latex(songs, output_path, max_songs=None):
    """Generate LaTeX document with ABC lead sheets."""
    preamble = r"""\documentclass[11pt,letterpaper]{article}
\usepackage[generate,ps2eps]{abc}
\usepackage[margin=0.75in]{geometry}
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\renewcommand{\headrulewidth}{0pt}
\cfoot{\thepage}
\renewcommand{\abcwidth}{0.95\linewidth}

\begin{document}
\begin{center}
{\Large\bfseries Harp Lead Sheets}\\[0.5em]
{\normalsize OpenHymnal -- Key of C}
\end{center}
\vspace{1em}
"""

    body_parts = []
    count = 0

    for song in songs:
        abc_block = generate_leadsheet_abc(song)
        if abc_block is None:
            continue

        count += 1
        if max_songs and count > max_songs:
            break

        safe_name = re.sub(r'[^a-zA-Z0-9]', '', f"h{song['number']}")

        body_parts.append(f'\\begin{{abc}}[name={safe_name}]')
        body_parts.append(abc_block)
        body_parts.append('\\end{abc}')
        body_parts.append('\\newpage')
        body_parts.append('')

    postamble = r"\end{document}"

    full_doc = preamble + '\n'.join(body_parts) + '\n' + postamble

    with open(output_path, 'w') as f:
        f.write(full_doc)

    print(f"Wrote {count} lead sheets to {output_path}")


# ── Main ─────────────────────────────────────────────────────────

def generate_abc_file(songs, output_path, max_songs=None):
    """Generate a single ABC file with all lead sheets."""
    parts = []
    count = 0

    for song in songs:
        abc_block = generate_leadsheet_abc(song)
        if abc_block is None:
            continue
        count += 1
        if max_songs and count > max_songs:
            break
        parts.append(abc_block)

    with open(output_path, 'w') as f:
        f.write('\n\n'.join(parts) + '\n')

    print(f"Wrote {count} lead sheets to {output_path}")


def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else 'OpenHymnal-C.abc'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'harp_leadsheets.abc'

    max_songs = None
    for arg in sys.argv[1:]:
        if arg.startswith('--max='):
            max_songs = int(arg.split('=')[1])

    if '--preview' in sys.argv:
        songs = parse_file(input_file)
        n = min(5, len(songs))
        for song in songs[:n]:
            abc = generate_leadsheet_abc(song)
            if abc:
                print(abc)
                print()
        return

    if '--latex' in sys.argv:
        tex_out = output_file.replace('.abc', '.tex')
        songs = parse_file(input_file)
        generate_latex(songs, tex_out, max_songs)
        print(f"Compile with: pdflatex -shell-escape {tex_out}")
        return

    songs = parse_file(input_file)
    generate_abc_file(songs, output_file, max_songs)


if __name__ == '__main__':
    main()
