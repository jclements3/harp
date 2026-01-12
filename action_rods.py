#!/usr/bin/env python3
"""
Action Rod System - Pedal to Disc Linkage

Routes action rods from 7 pedals through pillar to bell cranks on neck.

System overview:
- 7 pedals (C, D, E, F, G, A, B) at base
- 3 pedal positions: flat (up), natural (middle), sharp (down)
- Action rods run up through hollow pillar
- At neck, rods connect to bell cranks for each disc
- Each note controls ~7 strings across octaves

Pedal positions and disc states:
- Pedal UP (flat): both discs disengaged, string at full length (flat)
- Pedal MIDDLE (natural): natural disc engaged, string shortened 1 semitone
- Pedal DOWN (sharp): sharp disc engaged, string shortened 2 semitones

Rod routing (from video analysis):
1. Pedal connects to vertical rod in pillar base
2. Rod runs up through pillar (follows curve)
3. At neck junction, rod enters channel in neck
4. Rod runs along neck in groove, connecting to bell cranks
5. Each bell crank actuates one disc

For double-action (natural + sharp):
- Two parallel rods per note, or
- Single rod with intermediate linkage

Lyon & Healy / Erard use slightly different mechanisms but same principle.
"""

import json
import math
from pathlib import Path

import cadquery as cq
from cq_plate import plate


# Action rod specifications
ACTION_ROD_SPECS = {
    'rod_diameter_mm': 3.0,        # Main action rod diameter
    'rod_material': '316 Stainless Steel',
    'thread': 'M3',                # End threads for adjustment
    'notes': ['C', 'D', 'E', 'F', 'G', 'A', 'B'],
    'octaves': 7,                  # C1-C7, D1-D7, etc. (47 strings total)

    # Pillar dimensions (approximate)
    'pillar_height_mm': 1600,      # From base to neck junction
    'pillar_curve_radius_mm': 800, # Radius of pillar curve

    # Neck channel dimensions
    'channel_width_mm': 5.0,       # Width of rod channel in neck
    'channel_depth_mm': 8.0,       # Depth of channel
    'channel_spacing_mm': 7.0,     # Spacing between parallel channels
}


def get_strings_for_note(note, num_strings=47):
    """Get string numbers for a given note.

    Harp tuning: C D E F G A B repeating
    String 1 = C1, String 2 = D1, ... String 7 = B1
    String 8 = C2, etc.
    """
    notes = ['C', 'D', 'E', 'F', 'G', 'A', 'B']
    note_index = notes.index(note)

    strings = []
    for octave in range(7):  # Octaves 1-7
        string_num = octave * 7 + note_index + 1
        if string_num <= num_strings:
            strings.append(string_num)

    return strings


def calculate_rod_lengths():
    """Calculate action rod lengths for each note.

    Rod length depends on:
    - Path through pillar (curved)
    - Path along neck to furthest disc
    """
    specs = ACTION_ROD_SPECS

    # Approximate pillar path length (quarter circle arc)
    pillar_arc = math.pi * specs['pillar_curve_radius_mm'] / 2

    rod_lengths = {}
    for note in specs['notes']:
        strings = get_strings_for_note(note)

        # Neck path length (from pillar junction to furthest string)
        # Approximate: furthest string position along neck
        max_string = max(strings)
        neck_path = max_string * 30  # ~30mm per string position (approximate)

        total_length = pillar_arc + neck_path

        rod_lengths[note] = {
            'strings': strings,
            'pillar_length_mm': pillar_arc,
            'neck_length_mm': neck_path,
            'total_length_mm': total_length,
            'num_discs': len(strings) * 2,  # Natural + sharp per string
        }

    return rod_lengths


def make_rod_segment(length, diameter=3.0):
    """Create a straight rod segment."""
    return (cq.Workplane("XY")
            .circle(diameter / 2)
            .extrude(length))


def make_rod_channel(length, width=5.0, depth=8.0):
    """Create a channel/groove for rod routing."""
    return (cq.Workplane("XY")
            .rect(width, depth)
            .extrude(length))


def make_linkage_arm(length, thickness=3.0, width=8.0, hole_d=3.0):
    """Create a simple linkage arm with pivot holes at each end."""
    arm = (cq.Workplane("XY")
           .rect(length, width)
           .extrude(thickness))

    # Holes at each end
    hole_offset = length / 2 - hole_d
    arm = (arm.faces(">Z").workplane()
           .moveTo(hole_offset, 0)
           .circle(hole_d / 2)
           .cutThruAll())
    arm = (arm.faces(">Z").workplane()
           .moveTo(-hole_offset, 0)
           .circle(hole_d / 2)
           .cutThruAll())

    return arm


def make_pedal_box_schematic():
    """Create schematic of pedal box with 7 pedals.

    Standard pedal arrangement (player's view):
    Left foot: D C B
    Right foot: E F G A
    """
    box_width = 300   # mm
    box_depth = 150   # mm
    box_height = 100  # mm

    # Pedal box
    box = (cq.Workplane("XY")
           .rect(box_width, box_depth)
           .extrude(box_height))

    # Pedal slots (7 pedals)
    pedal_width = 25
    pedal_depth = 80
    pedal_spacing = 40

    # Left foot pedals: D, C, B (positions 0, 1, 2 from left)
    # Right foot pedals: E, F, G, A (positions 3, 4, 5, 6)
    pedal_positions = {
        'D': -120, 'C': -80, 'B': -40,
        'E': 20, 'F': 60, 'G': 100, 'A': 140
    }

    for note, x_pos in pedal_positions.items():
        slot = (cq.Workplane("XY")
                .workplane(offset=box_height)
                .moveTo(x_pos, 0)
                .rect(pedal_width, pedal_depth)
                .extrude(-box_height + 10))  # Leave bottom
        box = box.cut(slot)

    return box


def generate_rod_report():
    """Generate report of action rod specifications."""

    specs = ACTION_ROD_SPECS
    rod_lengths = calculate_rod_lengths()

    print("=== Action Rod System Specifications ===")
    print()
    print(f"Rod diameter: {specs['rod_diameter_mm']}mm")
    print(f"Material: {specs['rod_material']}")
    print(f"Thread: {specs['thread']}")
    print()
    print("Pedal arrangement (player's view):")
    print("  Left foot:  D  C  B")
    print("  Right foot: E  F  G  A")
    print()
    print("Rod lengths by note:")
    print("-" * 60)

    total_rod_length = 0
    total_discs = 0

    for note in specs['notes']:
        info = rod_lengths[note]
        print(f"  {note}: {info['total_length_mm']:.0f}mm total")
        print(f"     Strings: {info['strings']}")
        print(f"     Discs to actuate: {info['num_discs']}")
        print()
        total_rod_length += info['total_length_mm']
        total_discs += info['num_discs']

    print("-" * 60)
    print(f"Total rod length: {total_rod_length/1000:.1f}m")
    print(f"Total discs: {total_discs}")
    print()

    # Channel specifications
    print("Neck channel specifications:")
    print(f"  Width: {specs['channel_width_mm']}mm")
    print(f"  Depth: {specs['channel_depth_mm']}mm")
    print(f"  Spacing: {specs['channel_spacing_mm']}mm")
    print(f"  Total channels: 7 (one per note)")


def make_neck_cross_section():
    """Create cross-section of neck showing rod channels."""

    specs = ACTION_ROD_SPECS

    # Neck profile (simplified rectangle)
    neck_width = 80   # mm
    neck_height = 60  # mm

    neck = (cq.Workplane("XY")
            .rect(neck_width, neck_height)
            .extrude(10))  # Thin slice for visualization

    # Cut 7 channels for action rods
    channel_w = specs['channel_width_mm']
    channel_d = specs['channel_depth_mm']
    channel_spacing = specs['channel_spacing_mm']

    total_channel_width = 7 * channel_w + 6 * (channel_spacing - channel_w)
    start_x = -total_channel_width / 2 + channel_w / 2

    for i in range(7):
        x_pos = start_x + i * channel_spacing
        channel = (cq.Workplane("XY")
                   .workplane(offset=10)
                   .moveTo(x_pos, neck_height/2 - channel_d/2)
                   .rect(channel_w, channel_d)
                   .extrude(-10))
        neck = neck.cut(channel)

    return neck


def make_schematic_plate():
    """Generate schematic showing action rod system."""

    # Pedal box
    pedal_box = make_pedal_box_schematic()

    # Neck cross-section
    neck_section = make_neck_cross_section()
    neck_section = neck_section.translate((0, 200, 0))

    groups = [
        (pedal_box, '#666666', 0.3),
        (neck_section, '#996633', 0.4),
    ]

    svg = plate(
        groups,
        UL='XY',   # Top view
        LL='XZ',   # Front view
        UR='ISO',
        LR='YZ',
        grid='light',
        units='mm'
    )

    return svg


def main():
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == '--svg':
        svg = make_schematic_plate()
        output = Path(__file__).parent / 'action_rods_schematic.svg'
        output.write_text(svg)
        print(f"Generated {output}")
    else:
        generate_rod_report()


if __name__ == '__main__':
    main()
