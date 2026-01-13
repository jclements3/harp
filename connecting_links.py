#!/usr/bin/env python3
"""
Connecting Links - Bell Crank to Bell Crank Linkage

Links connect adjacent bell cranks to transfer pedal motion along the neck.
Each note has a chain of bell cranks connected by these links.

Naming convention:
    c1c2ln = link between C1 and C2 natural bell cranks
    c1c2ls = link between C1 and C2 sharp bell cranks

Link geometry:
    Simple bar with pivot holes at each end.
    Length spans string spacing (finger_gap).

    +--( )------------------( )--+
       ^                      ^
    pivot hole            pivot hole
    (to C1 crank)         (to C2 crank)

Counts:
    - 7 notes (C, D, E, F, G, A, B)
    - 6 links per note per action (1-2, 2-3, 3-4, 4-5, 5-6, 6-7)
    - 2 actions (natural, sharp)
    - Total: 7 × 6 × 2 = 84 links
"""

import math
from pathlib import Path

import cadquery as cq
from cq_plate import plate


# Link specifications by register
# Material: 316 Stainless Steel sheet (1mm)
# Link length is approximately finger_gap (29.16mm)
# 7 rows of links (A-G) stacked in Z between front and back plates
LINK_SPECS = {
    'bass': {
        'strings': '1-9',
        'length': 29.0,           # mm - center-to-center of pivot holes
        'width': 5.0,             # mm - bar width
        'thickness': 1.0,         # mm - 1mm stainless sheet
        'hole_diameter': 2.5,     # mm - pivot hole
        'end_radius': 3.5,        # mm - rounded end radius
    },
    'mid': {
        'strings': '10-28',
        'length': 29.0,
        'width': 4.5,
        'thickness': 1.0,
        'hole_diameter': 2.0,
        'end_radius': 3.0,
    },
    'treble': {
        'strings': '29-38',
        'length': 29.0,
        'width': 4.0,
        'thickness': 1.0,
        'hole_diameter': 1.5,
        'end_radius': 2.5,
    },
    'high': {
        'strings': '39-47',
        'length': 29.0,
        'width': 3.5,
        'thickness': 1.0,
        'hole_diameter': 1.5,
        'end_radius': 2.0,
    },
}

# All link names
# Row order matches pedal arrangement: D C B (gap) E F G A from +Z to -Z
# Rods cross over Z axis, so row order front plate to back plate is:
ROW_ORDER = ['d', 'c', 'b', 'e', 'f', 'g', 'a']  # Front plate to back plate
NOTES = ['c', 'd', 'e', 'f', 'g', 'a', 'b']  # Alphabetical for iteration
OCTAVE_PAIRS = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)]
ACTIONS = ['n', 's']  # natural, sharp


def get_link_name(note, oct1, oct2, action):
    """Generate link name.

    Args:
        note: 'c', 'd', 'e', 'f', 'g', 'a', 'b'
        oct1: First octave (1-6)
        oct2: Second octave (2-7)
        action: 'n' (natural) or 's' (sharp)

    Returns:
        String like 'c1c2ln' or 'd3d4ls'
    """
    return f"{note}{oct1}{note}{oct2}l{action}"


def get_all_link_names():
    """Get list of all 84 link names."""
    names = []
    for note in NOTES:
        for oct1, oct2 in OCTAVE_PAIRS:
            for action in ACTIONS:
                names.append(get_link_name(note, oct1, oct2, action))
    return names


def get_register_for_string(string_num):
    """Get register name for a string number."""
    if string_num <= 9:
        return 'bass'
    elif string_num <= 28:
        return 'mid'
    elif string_num <= 38:
        return 'treble'
    else:
        return 'high'


def get_register_for_link(note, octave):
    """Get register for a link based on note and lower octave.

    Links span between two strings, use the lower string's register.
    """
    # String numbering: C1=1, D1=2, E1=3, F1=4, G1=5, A1=6, B1=7, C2=8, ...
    note_index = NOTES.index(note)
    string_num = (octave - 1) * 7 + note_index + 1
    return get_register_for_string(string_num)


def get_row_z_offset(note):
    """Get Z offset for a note's link row relative to front plate back face.

    Returns negative Z value (going toward back plate).
    Row order: D, C, B, E, F, G, A (front to back)
    """
    row_index = ROW_ORDER.index(note.lower())
    thickness = 1.0  # mm - all links are 1mm thick
    spacing = 0.6    # mm - between rows

    # First row starts with small clearance from front plate
    clearance = 0.5
    z_offset = -(clearance + row_index * (thickness + spacing) + thickness / 2)

    return z_offset


def get_row_positions():
    """Get Z positions for all 7 link rows.

    Returns dict mapping note to Z offset from front plate back face.
    """
    return {note: get_row_z_offset(note) for note in ROW_ORDER}


def make_connecting_link(spec):
    """Create a connecting link bar.

    Args:
        spec: Dict with length, width, thickness, hole_diameter, end_radius

    Returns:
        CadQuery solid of the link
    """
    length = spec['length']
    width = spec['width']
    thickness = spec['thickness']
    hole_d = spec['hole_diameter']
    end_r = spec['end_radius']

    # Ensure end radius accommodates the hole
    if end_r < hole_d / 2 + 1.0:
        end_r = hole_d / 2 + 1.5

    # Create stadium-shaped bar (slot shape)
    # Slot2D creates a stadium with given length and width
    bar = (cq.Workplane("XY")
           .slot2D(length, width, angle=0)
           .extrude(thickness))

    # Center at origin (slot2D is already centered)
    bar = bar.translate((0, 0, -thickness / 2))

    # Add pivot holes at each end
    hole_offset = (length - width) / 2  # Distance from center to hole center

    # Left hole
    bar = (bar.faces(">Z")
           .workplane()
           .moveTo(-hole_offset, 0)
           .circle(hole_d / 2)
           .cutThruAll())

    # Right hole
    bar = (bar.faces(">Z")
           .workplane()
           .moveTo(hole_offset, 0)
           .circle(hole_d / 2)
           .cutThruAll())

    return bar


def make_link_by_name(link_name):
    """Create a link by its name (e.g., 'c1c2ln').

    Args:
        link_name: String like 'c1c2ln' or 'd3d4ls'

    Returns:
        CadQuery solid of the link
    """
    # Parse name: note, oct1, note, oct2, 'l', action
    note = link_name[0]
    oct1 = int(link_name[1])
    # link_name[2] is the note again
    oct2 = int(link_name[3])
    # link_name[4] is 'l'
    action = link_name[5]  # 'n' or 's'

    register = get_register_for_link(note, oct1)
    spec = LINK_SPECS[register]

    return make_connecting_link(spec)


def make_dimensioned_views(register):
    """Create multi-view plate for manufacturing."""
    spec = LINK_SPECS[register]
    link = make_connecting_link(spec)

    groups = [(link, '#cc8844', 0.3)]  # Brass color

    svg = plate(
        groups,
        UL='XZ',   # Side view (shows thickness)
        LL='XY',   # Top view (shows length and width)
        UR='ISO',  # Isometric
        LR='YZ',   # End view
        grid='light',
        units='mm'
    )

    return svg


def generate_all_plates():
    """Generate manufacturing plates for all link sizes."""
    output_dir = Path(__file__).parent

    print("=== Connecting Links ===")
    print()
    print("Naming: {note}{oct1}{note}{oct2}l{action}")
    print("  e.g., c1c2ln = C1 to C2 natural link")
    print()

    for register, spec in LINK_SPECS.items():
        svg = make_dimensioned_views(register)

        output_path = output_dir / f"link_{register}.svg"
        output_path.write_text(svg)

        print(f"{register.upper()} (strings {spec['strings']}):")
        print(f"  Length: {spec['length']}mm (center-to-center)")
        print(f"  Width: {spec['width']}mm")
        print(f"  Thickness: {spec['thickness']}mm")
        print(f"  Pivot hole: Ø{spec['hole_diameter']}mm")
        print(f"  Output: {output_path}")
        print()

    # Print link count summary
    all_links = get_all_link_names()
    print(f"Total links: {len(all_links)}")
    print()
    print("Links per note:")
    for note in NOTES:
        note_links = [n for n in all_links if n.startswith(note)]
        print(f"  {note.upper()}: {len(note_links)} ({', '.join(note_links[:3])}...)")


def generate_link_report():
    """Generate detailed report of all links."""
    print("=== Complete Link Inventory ===")
    print()

    # Show row order and Z positions
    print("Row order (front plate to back plate):")
    print("  Pedals +Z to -Z: D C B (gap) E F G A")
    print("  Rods cross over Z axis")
    print()
    positions = get_row_positions()
    print("  Row   Z offset (from front plate)")
    print("  ---   --------------------------")
    for note in ROW_ORDER:
        z = positions[note]
        print(f"  {note.upper()}     {z:.1f}mm")
    print()

    all_links = get_all_link_names()

    # Group by note in row order
    print("Links by row (front to back):")
    for note in ROW_ORDER:
        print(f"{note.upper()} Links (z={positions[note]:.1f}mm):")
        for action in ACTIONS:
            action_name = "Natural" if action == 'n' else "Sharp"
            links = [n for n in all_links if n.startswith(note) and n.endswith(action)]
            print(f"  {action_name}: {', '.join(links)}")
        print()

    # Count by register
    print("Count by register:")
    for register in LINK_SPECS.keys():
        count = 0
        for note in NOTES:
            for oct1, oct2 in OCTAVE_PAIRS:
                reg = get_register_for_link(note, oct1)
                if reg == register:
                    count += 2  # natural + sharp
        print(f"  {register}: {count} links")


def main():
    import sys

    if len(sys.argv) > 1:
        arg = sys.argv[1]
        if arg == '--report':
            generate_link_report()
        elif arg in LINK_SPECS:
            svg = make_dimensioned_views(arg)
            output = Path(__file__).parent / f"link_{arg}.svg"
            output.write_text(svg)
            print(f"Generated {output}")
        else:
            print(f"Unknown argument: {arg}")
            print("Usage: python connecting_links.py [bass|mid|treble|high|--report]")
    else:
        generate_all_plates()


if __name__ == '__main__':
    main()
