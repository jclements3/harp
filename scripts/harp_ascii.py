#!/usr/bin/env python3
"""Generate ASCII harp hymnal — 100 columns wide, measure-aligned wrapping."""

import json
import sys


PAGE_WIDTH = 100
PAGE_HEIGHT = 66  # lines per page (standard terminal/print page)
MARGIN_TOP = 2
MARGIN_BOTTOM = 2
USABLE_HEIGHT = PAGE_HEIGHT - MARGIN_TOP - MARGIN_BOTTOM
STACK_ROWS = 10  # chord_rh, chord_lh, rh[0..3], lh[0..3]
MIN_TITLE_SPACE = STACK_ROWS + 3  # title + key + one stack + gap


def format_hymn(hymn):
    """Format one hymn into lines of text. Returns list of strings."""
    beats = hymn['beats']

    # Build columns with measure boundaries
    cols = []
    for beat in beats:
        if beat['type'] == 'bar':
            cols.append({'type': 'bar'})
        else:
            rh = beat.get('rh', [])
            lh = beat.get('lh', [])
            rc = beat.get('rh_chord', '') if beat.get('chord') else ''
            lc = beat.get('lh_chord', '') if beat.get('chord') else ''
            tick = "'" if lc and len(lh) > 0 and len(lh) < 4 else ''

            cols.append({
                'type': 'beat',
                'rh_chord': rc,
                'lh_chord': f"{tick}{lc}",
                'rh': rh,
                'lh': lh,
                'is_chord': bool(beat.get('chord')),
            })

    # Split into measures (groups between barlines), skip empty measures
    measures = []
    current = []
    for c in cols:
        if c['type'] == 'bar':
            if current:
                measures.append(current)
                current = []
        else:
            current.append(c)
    if current:
        measures.append(current)
    measures = [m for m in measures if m]  # drop empties

    # Compute column width for each beat: max width of any cell in that column
    def measure_width(measure):
        """Width of a measure including trailing barline."""
        w = 0
        for beat in measure:
            col_w = 1  # minimum
            for val in [beat['rh_chord'], beat['lh_chord']]:
                col_w = max(col_w, len(val))
            for s in beat['rh'] + beat['lh']:
                col_w = max(col_w, len(str(s)))
            w += col_w + 1  # +1 for space between columns
        return w + 2  # +2 for "| " barline

    # Group measures into lines that fit within PAGE_WIDTH
    lines_of_measures = []
    current_line = []
    current_width = 0

    for m in measures:
        mw = measure_width(m)
        if current_width + mw > PAGE_WIDTH and current_line:
            lines_of_measures.append(current_line)
            current_line = []
            current_width = 0
        current_line.append(m)
        current_width += mw

    if current_line:
        lines_of_measures.append(current_line)

    # Render each line of measures into STACK_ROWS of text
    output_lines = []

    # Title
    key = hymn['key']
    tempo = hymn['tempo']
    title = f"{hymn['number']}. {hymn['title']} ({key}, {tempo}bpm)"
    output_lines.append(title)

    for line_measures in lines_of_measures:
        # Build the 10 rows for this line
        rows = [[] for _ in range(STACK_ROWS)]

        for mi, measure in enumerate(line_measures):
            if mi > 0:
                # Barline separator
                for r in rows:
                    r.append('|')

            for beat in measure:
                # Compute this column's width
                col_w = 1
                for val in [beat['rh_chord'], beat['lh_chord']]:
                    col_w = max(col_w, len(val))
                for s in beat['rh'] + beat['lh']:
                    col_w = max(col_w, len(str(s)))

                # Row 0: RH chord
                rows[0].append(beat['rh_chord'].rjust(col_w))
                # Row 1: LH chord
                rows[1].append(beat['lh_chord'].rjust(col_w))
                # Rows 2-5: RH fingers
                for i in range(4):
                    if beat['is_chord'] and i < len(beat['rh']):
                        rows[2 + i].append(str(beat['rh'][i]).rjust(col_w))
                    elif i == 0 and beat['rh']:
                        rows[2 + i].append(str(beat['rh'][0]).rjust(col_w))
                    else:
                        rows[2 + i].append(' ' * col_w)
                # Rows 6-9: LH fingers
                for i in range(4):
                    if beat['is_chord'] and i < len(beat['lh']):
                        rows[6 + i].append(str(beat['lh'][i]).rjust(col_w))
                    else:
                        rows[6 + i].append(' ' * col_w)

        # Join each row
        for r in rows:
            line = ' '.join(r).rstrip()
            output_lines.append(line)

        # Blank line between stacks
        output_lines.append('')

    return output_lines


def generate_ascii(voicings_file, output_file, max_hymns=None):
    hymns = json.load(open(voicings_file))
    if max_hymns:
        hymns = hymns[:max_hymns]

    pages = []
    current_page = []
    current_line_count = 0
    header = 'OpenHymnal -- Harp Lead Sheets'

    def flush_page():
        nonlocal current_page, current_line_count
        pages.append(current_page)
        current_page = []
        current_line_count = 0

    for hymn in hymns:
        hymn_lines = format_hymn(hymn)

        # Only break if title + at least one stack won't fit
        if current_line_count > 0 and current_line_count + MIN_TITLE_SPACE > USABLE_HEIGHT:
            flush_page()

        # Split hymn into blocks (title + stacks separated by blank lines)
        blocks = []
        current_block = []
        for line in hymn_lines:
            if line == '' and current_block:
                blocks.append(current_block)
                current_block = []
            else:
                current_block.append(line)
        if current_block:
            blocks.append(current_block)

        # Add blocks, breaking pages between blocks only
        for block in blocks:
            # If block won't fit, break page first
            if current_line_count > 0 and current_line_count + len(block) + 1 > USABLE_HEIGHT:
                flush_page()
            for line in block:
                current_page.append(line)
                current_line_count += 1
            # Blank line after block
            if current_line_count < USABLE_HEIGHT:
                current_page.append('')
                current_line_count += 1

    if current_page:
        flush_page()

    # Write output with headers and footers
    with open(output_file, 'w') as f:
        for pi, page in enumerate(pages):
            if pi > 0:
                f.write('\f')
            # Header
            page_num = f'Page {pi + 1}'
            pad = PAGE_WIDTH - len(header) - len(page_num)
            f.write(header + ' ' * max(pad, 1) + page_num + '\n')
            f.write('-' * PAGE_WIDTH + '\n')
            # Content
            for line in page:
                f.write(line + '\n')

    print(f"Wrote {len(hymns)} hymns, {len(pages)} pages to {output_file}")


if __name__ == '__main__':
    voicings = sys.argv[1] if len(sys.argv) > 1 else 'harp_voicings.json'
    output = sys.argv[2] if len(sys.argv) > 2 else 'OpenHymnal-harp.txt'
    max_hymns = None
    for arg in sys.argv[1:]:
        if arg.startswith('--max='):
            max_hymns = int(arg.split('=')[1])
    generate_ascii(voicings, output, max_hymns)
