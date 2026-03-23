#!/usr/bin/env python3
"""Generate harp lead sheet PDF from precomputed voicings."""

import json
import sys
from reportlab.lib.pagesizes import letter, landscape
from reportlab.lib import colors
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Spacer, Paragraph, PageBreak
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT

PAGE = landscape(letter)
MARGIN = 0.5 * inch


def build_hymn_table(hymn, page_width=None):
    """Build a reportlab Table for one hymn's lead sheet."""
    if page_width is None:
        page_width = PAGE[0] - 2 * MARGIN
    beats = hymn['beats']

    # Collect columns
    cols = []
    for beat in beats:
        if beat['type'] == 'bar':
            cols.append({'type': 'bar'})
        else:
            rh = beat.get('rh', [])
            lh = beat.get('lh', [])
            cols.append({
                'type': 'beat',
                'chord': beat.get('chord', ''),
                'rh_chord': beat.get('rh_chord', ''),
                'lh_chord': beat.get('lh_chord', ''),
                'rh': rh,
                'lh': lh,
                'lh_tick': len(lh) > 0 and len(lh) < 4,
            })

    # Build 10 rows: rh_chord, lh_chord, rh[0..3], lh[0..3]
    row_rh_chord = []
    row_lh_chord = []
    rh_rows = [[] for _ in range(4)]
    lh_rows = [[] for _ in range(4)]

    for c in cols:
        if c['type'] == 'bar':
            row_rh_chord.append('|')
            row_lh_chord.append('|')
            for r in rh_rows: r.append('|')
            for r in lh_rows: r.append('|')
        else:
            row_rh_chord.append(c['rh_chord'] if c['chord'] else '')
            tick = "'" if c['lh_tick'] else ''
            row_lh_chord.append(f"{tick}{c['lh_chord']}" if c['chord'] else '')

            for i in range(4):
                rh_rows[i].append(str(c['rh'][i]) if i < len(c['rh']) and c['chord'] else
                                  (str(c['rh'][0]) if i == 0 and c['rh'] else ''))
                lh_rows[i].append(str(c['lh'][i]) if i < len(c['lh']) and c['chord'] else '')

    # Combine all rows
    all_rows = [row_rh_chord, row_lh_chord] + rh_rows + lh_rows
    num_cols = len(all_rows[0]) if all_rows else 0

    if num_cols == 0:
        return None

    # Column widths
    col_w = 0.28 * inch
    max_cols_per_line = int(page_width / col_w)

    # Split into lines that wrap
    tables = []
    for start in range(0, num_cols, max_cols_per_line):
        end = min(start + max_cols_per_line, num_cols)
        chunk_data = []
        for row in all_rows:
            chunk_data.append(row[start:end])

        chunk_cols = end - start
        t = Table(chunk_data, colWidths=[col_w] * chunk_cols, hAlign='LEFT')

        style_cmds = [
            ('FONTNAME', (0, 0), (-1, -1), 'Courier'),
            ('FONTSIZE', (0, 0), (-1, -1), 6),
            ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('TOPPADDING', (0, 0), (-1, -1), 0),
            ('BOTTOMPADDING', (0, 0), (-1, -1), 0),
            ('LEFTPADDING', (0, 0), (-1, -1), 1),
            ('RIGHTPADDING', (0, 0), (-1, -1), 1),
            # RH chord row - blue
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
            ('FONTNAME', (0, 0), (-1, 0), 'Courier-Bold'),
            # LH chord row - blue
            ('TEXTCOLOR', (0, 1), (-1, 1), colors.black),
            ('FONTNAME', (0, 1), (-1, 1), 'Courier-Bold'),
            # LH thumb row highlight
            ('TEXTCOLOR', (0, 6), (-1, 6), colors.black),
            ('FONTNAME', (0, 6), (-1, 6), 'Courier-Bold'),
            # Light line between RH and LH
            ('LINEABOVE', (0, 6), (-1, 6), 0.3, colors.HexColor('#999999')),
        ]

        # Color barline columns
        for ci in range(chunk_cols):
            if chunk_data[0][ci] == '|':
                for ri in range(len(all_rows)):
                    style_cmds.append(('TEXTCOLOR', (ci, ri), (ci, ri), colors.HexColor('#999999')))

        t.setStyle(TableStyle(style_cmds))
        tables.append(t)

    return tables


def generate_pdf(voicings_file, output_file, max_hymns=None, use_landscape=True):
    hymns = json.load(open(voicings_file))
    if max_hymns:
        hymns = hymns[:max_hymns]

    page = landscape(letter) if use_landscape else letter

    doc = SimpleDocTemplate(
        output_file,
        pagesize=page,
        leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=MARGIN, bottomMargin=MARGIN,
    )

    title_style = ParagraphStyle('Title',
        fontName='Helvetica-Bold', fontSize=9,
        textColor=colors.black,
        spaceAfter=2, alignment=TA_LEFT)

    key_style = ParagraphStyle('Key',
        fontName='Courier-Bold', fontSize=7,
        textColor=colors.black,
        spaceAfter=1, alignment=TA_LEFT)

    elements = []

    for i, hymn in enumerate(hymns):
        # Title
        elements.append(Paragraph(
            f"{hymn['number']}. {hymn['title']} ({hymn['key']}, {hymn['tempo']}bpm)", title_style))

        # Lead sheet table(s)
        tables = build_hymn_table(hymn, page_width=page[0] - 2 * MARGIN)
        if tables:
            for t in tables:
                elements.append(t)
                elements.append(Spacer(1, 1))

        # Small gap between hymns (no page break — let reportlab flow)
        elements.append(Spacer(1, 8))

    doc.build(elements)
    print(f"Wrote {len(hymns)} hymns to {output_file}")


if __name__ == '__main__':
    voicings = sys.argv[1] if len(sys.argv) > 1 else 'harp_voicings.json'
    output = sys.argv[2] if len(sys.argv) > 2 else 'OpenHymnal-harp.pdf'
    max_hymns = None
    use_landscape = '--landscape' in sys.argv
    for arg in sys.argv[1:]:
        if arg.startswith('--max='):
            max_hymns = int(arg.split('=')[1])
    generate_pdf(voicings, output, max_hymns, use_landscape)
    if not use_landscape:
        print("(portrait mode)")
