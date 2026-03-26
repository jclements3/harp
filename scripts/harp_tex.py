#!/usr/bin/env python3
"""Generate harp lead sheet LaTeX from precomputed voicings."""

import json
import sys


def escape_tex(s):
    return s.replace('#', r'\#').replace('~', r'\textasciitilde{}').replace('%', r'\%').replace('^', r'\^{}').replace('_', r'\_')


def build_hymn_latex(hymn, max_cols=35):
    """Build LaTeX tabular blocks for one hymn."""
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

    # Build 10 rows
    row_rh_chord = []
    row_lh_chord = []
    rh_rows = [[] for _ in range(4)]
    lh_rows = [[] for _ in range(4)]

    for c in cols:
        if c['type'] == 'bar':
            row_rh_chord.append('\\textbar')
            row_lh_chord.append('\\textbar')
            for r in rh_rows: r.append('\\textbar')
            for r in lh_rows: r.append('\\textbar')
        else:
            rc = escape_tex(c['rh_chord']) if c['chord'] else ''
            tick = "'" if c['lh_tick'] else ''
            lc = escape_tex(f"{tick}{c['lh_chord']}") if c['chord'] else ''
            row_rh_chord.append(rc)
            row_lh_chord.append(lc)

            for i in range(4):
                if c['chord'] and i < len(c['rh']):
                    rh_rows[i].append(str(c['rh'][i]))
                elif i == 0 and c['rh']:
                    rh_rows[i].append(str(c['rh'][0]))
                else:
                    rh_rows[i].append('')

                if c['chord'] and i < len(c['lh']):
                    lh_rows[i].append(str(c['lh'][i]))
                else:
                    lh_rows[i].append('')

    all_rows = [row_rh_chord, row_lh_chord] + rh_rows + lh_rows
    num_cols = len(all_rows[0]) if all_rows else 0
    if num_cols == 0:
        return ''

    # Split into chunks that fit page width
    chunks = []
    for start in range(0, num_cols, max_cols):
        end = min(start + max_cols, num_cols)
        n = end - start
        col_spec = '@{}' + 'c@{\\,}' * n + '@{}'

        lines = []
        lines.append(f'\\begin{{tabular}}{{{col_spec}}}')

        for ri, row in enumerate(all_rows):
            cells = row[start:end]
            # Bold chord rows
            if ri < 2:
                cells = [f'\\textbf{{{c}}}' if c and c != '\\textbar' else c for c in cells]
            # Bold LH thumb row (ri=6)
            if ri == 6:
                cells = [f'\\textbf{{{c}}}' if c and c != '\\textbar' else c for c in cells]

            line = ' & '.join(cells) + ' \\\\'
            # Add thin rule between RH and LH
            if ri == 5:
                line += '[1pt]'
            lines.append(line)

        lines.append('\\end{tabular}')
        chunks.append('\n'.join(lines))

    return '\n\\par\\noindent\n'.join(chunks)


def generate_tex(voicings_file, output_file, max_hymns=None):
    hymns = json.load(open(voicings_file))
    if max_hymns:
        hymns = hymns[:max_hymns]

    parts = []
    parts.append(r"""\documentclass[12pt]{article}
\usepackage[letterpaper,margin=0.4in]{geometry}
\usepackage[T1]{fontenc}
\usepackage{inconsolata}
\renewcommand{\familydefault}{\ttdefault}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0pt}
\pagestyle{plain}
\begin{document}
\begin{center}
{\Large\bfseries OpenHymnal -- Harp Lead Sheets}\\[2pt]
{\small 292 hymns, 6--8 note voicings, LH|RH notation}
\end{center}
\vspace{4pt}
""")

    for i, hymn in enumerate(hymns):
        title = escape_tex(hymn['title'])
        key = escape_tex(hymn['key'])
        tempo = hymn['tempo']

        parts.append(f'\\noindent\\textbf{{{hymn["number"]}. {title}}} ({key}, {tempo}bpm)\\par')
        parts.append('{\\scriptsize')

        table_tex = build_hymn_latex(hymn)
        if table_tex:
            parts.append('\\noindent')
            parts.append(table_tex)

        parts.append('}')
        parts.append('\\vspace{6pt}')
        parts.append('')

    parts.append(r'\end{document}')

    with open(output_file, 'w') as f:
        f.write('\n'.join(parts))

    print(f"Wrote {len(hymns)} hymns to {output_file}")
    print(f"Compile: pdflatex {output_file}")


if __name__ == '__main__':
    voicings = sys.argv[1] if len(sys.argv) > 1 else 'harp_voicings.json'
    output = sys.argv[2] if len(sys.argv) > 2 else 'OpenHymnal-harp.tex'
    max_hymns = None
    for arg in sys.argv[1:]:
        if arg.startswith('--max='):
            max_hymns = int(arg.split('=')[1])
    generate_tex(voicings, output, max_hymns)
