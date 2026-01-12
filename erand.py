#!/usr/bin/env python3
"""
Erard Concert Harp - CadQuery Model

Reads erand.json and generates multi-view SVG using cq_plate().
Colors: C strings red, F strings blue, others gray.
String thickness varies by diameter.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import cadquery as cq
from cq_plate import plate


def load_data():
    with open(Path(__file__).parent / "erand.json") as f:
        return json.load(f)


def note_for_string(num):
    """Return note letter for string number (1-47)."""
    return ['C', 'D', 'E', 'F', 'G', 'A', 'B'][(num - 1) % 7]


def color_for_note(note):
    """Return color for note: C=red, F=blue, others=gray."""
    if note == 'C':
        return '#cc0000'
    elif note == 'F':
        return '#0000cc'
    return '#666666'


def diameter_to_stroke(diameter, min_dia=0.64, max_dia=2.64):
    """Scale diameter to visible stroke width.

    At SVG scale ~0.09, a 14.5mm gap becomes ~1.3mm.
    Use stroke range 0.1 to 0.5 to show variation without overlap.
    """
    t = (diameter - min_dia) / (max_dia - min_dia)
    return 0.1 + t * 0.4


def build_groups(data):
    """Build geometry groups for cq_plate().

    Returns list of (geometry, color, stroke_width) tuples.
    Each string gets its own group for individual stroke width.
    """
    groups = []

    # Soundboard line
    sb = data['soundboard']
    sb_edge = cq.Edge.makeLine(
        cq.Vector(sb['x1'], sb['y1'], 0),
        cq.Vector(sb['x2'], sb['y2'], 0)
    )
    groups.append((
        cq.Workplane("XY").newObject([cq.Compound.makeCompound([sb_edge])]),
        '#666666',
        0.3
    ))

    # Collect markers by color (fixed thin stroke)
    red_markers = []
    blue_markers = []
    gray_markers = []

    for s in data['strings']:
        note = note_for_string(s['num'])
        color = color_for_note(note)
        b, t = s['b'], s['t']
        stroke = diameter_to_stroke(s['diameter'])

        # String edges (all same diameter and color)
        # Main string: b to t
        # Peg segment: part of string path
        string_edges = []

        # Main string line (add tiny X offset to avoid CadQuery division by zero)
        string_edges.append(cq.Edge.makeLine(
            cq.Vector(b['x'], b['y'], 0),
            cq.Vector(t['x'] + 0.01, t['y'], 0)
        ))

        # Peg segment is continuation of string
        peg = s.get('peg')
        if peg and 'x2' in peg:
            string_edges.append(cq.Edge.makeLine(
                cq.Vector(peg['x'], peg['y'], 0),
                cq.Vector(peg['x2'], peg['y2'], 0)
            ))

        groups.append((
            cq.Workplane("XY").newObject([cq.Compound.makeCompound(string_edges)]),
            color,
            stroke
        ))

        # peg, pin, ndisc, sdisc tangent markers (1mm vertical lines on right of string)
        marker_edges = []
        marker_half_height = 0.5  # mm (1mm total)
        for key in ['peg', 'pin', 'ndisc', 'sdisc']:
            m = s.get(key)
            if m:
                # Create vertical marker line at tangent point (right side of string)
                marker_edges.append(cq.Edge.makeLine(
                    cq.Vector(m['x'] + 0.01, m['y'] - marker_half_height, 0),
                    cq.Vector(m['x'] + 0.01, m['y'] + marker_half_height, 0)
                ))

        if note == 'C':
            red_markers.extend(marker_edges)
        elif note == 'F':
            blue_markers.extend(marker_edges)
        else:
            gray_markers.extend(marker_edges)

    # Add marker groups with fixed thin stroke
    marker_stroke = 0.15
    if red_markers:
        groups.append((
            cq.Workplane("XY").newObject([cq.Compound.makeCompound(red_markers)]),
            '#cc0000', marker_stroke
        ))
    if blue_markers:
        groups.append((
            cq.Workplane("XY").newObject([cq.Compound.makeCompound(blue_markers)]),
            '#0000cc', marker_stroke
        ))
    if gray_markers:
        groups.append((
            cq.Workplane("XY").newObject([cq.Compound.makeCompound(gray_markers)]),
            '#666666', marker_stroke
        ))

    return groups


def main():
    data = load_data()
    groups = build_groups(data)

    svg = plate(groups, UL='XZ', LL='XY', UR='ISO', LR='YZ',
                grid='light', units='mm')

    output = Path(__file__).parent / "erand.svg"
    output.write_text(svg)

    c_count = sum(1 for s in data['strings'] if note_for_string(s['num']) == 'C')
    f_count = sum(1 for s in data['strings'] if note_for_string(s['num']) == 'F')

    print(f"Generated {output}")
    print(f"  Strings: {len(data['strings'])}")
    print(f"  C strings (red): {c_count}")
    print(f"  F strings (blue): {f_count}")
    print(f"  Other (gray): {len(data['strings']) - c_count - f_count}")

    # Show diameter range
    dias = [s['diameter'] for s in data['strings']]
    print(f"  Diameter range: {min(dias):.2f}mm - {max(dias):.2f}mm")


if __name__ == '__main__':
    main()
