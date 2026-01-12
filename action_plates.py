#!/usr/bin/env python3
"""
Action Plate Assembly - CadQuery Models

Front and back plates that hold the disc axles and bell crank pivots.
Based on analysis of Erard and Lyon & Healy mechanism designs.

Manufacturing: 316 Stainless Steel, waterjet or laser cut + machined holes
"""

import json
import math
from pathlib import Path

import cadquery as cq
from cq_plate import plate


def load_config():
    """Load harp configuration."""
    config_path = Path(__file__).parent / 'erand.json'
    with open(config_path) as f:
        return json.load(f)


def compute_disc_positions(config):
    """Compute disc positions for all strings.

    Returns list of dicts with x, y, z positions for natural and sharp discs.
    """
    sp = config['sp']
    strings = config['strings']
    geometry = config['geometry']
    physics = config['physics']
    hardware = config['hardware']

    c1b = (sp['c1b']['x'], sp['c1b']['y'])
    g7b = (sp['g7b']['x'], sp['g7b']['y'])

    dx = g7b[0] - c1b[0]
    dy = g7b[1] - c1b[1]
    sp_length = math.sqrt(dx*dx + dy*dy)
    ux, uy = dx/sp_length, dy/sp_length

    finger_gap = geometry['finger_gap']
    string_tilt = math.radians(geometry['string_tilt_deg'])
    semitone_ratio = physics['semitone_ratio']

    # String direction
    string_dx = math.sin(string_tilt)
    string_dy = math.cos(string_tilt)

    positions = []
    for i, s in enumerate(strings):
        string_num = i + 1

        # Soundboard position
        bx = c1b[0] + ux * i * finger_gap
        by = c1b[1] + uy * i * finger_gap

        vib_length = s['vib_length']

        # Flat pin position
        fx = bx + string_dx * vib_length
        fy = by + string_dy * vib_length

        # Natural disc position (one semitone shorter)
        vib_n = vib_length / semitone_ratio
        nx = bx + string_dx * vib_n
        ny = by + string_dy * vib_n

        # Sharp disc position (two semitones shorter)
        vib_s = vib_length / (semitone_ratio ** 2)
        sx = bx + string_dx * vib_s
        sy = by + string_dy * vib_s

        # Get disc specs for this string
        if string_num <= 9:
            disc = hardware['discs']['bass']
            register = 'bass'
        elif string_num <= 28:
            disc = hardware['discs']['mid']
            register = 'mid'
        elif string_num <= 38:
            disc = hardware['discs']['treble']
            register = 'treble'
        else:
            disc = hardware['discs']['high']
            register = 'high'

        sphere_r = disc['sphere_radius_mm']
        thickness = disc['thickness_mm']
        axle_d = disc['axle_diameter_mm']

        # Z position for disc center
        disc_z = -(sphere_r + thickness / 2)

        positions.append({
            'string': string_num,
            'note': s['note'],
            'register': register,
            'n': {'x': nx, 'y': ny, 'z': disc_z},
            's': {'x': sx, 'y': sy, 'z': disc_z},
            'axle_diameter': axle_d,
            'sphere_radius': sphere_r,
        })

    return positions


def make_front_plate(config, positions):
    """Create front plate with disc axle holes and bell crank pivot holes.

    The plate follows the curved path of disc positions.
    """
    action = config['hardware']['action']
    thickness = action['front_plate']['thickness_mm']
    clearance = action['front_plate']['hole_clearance_mm']

    # Calculate bounding box from disc positions
    all_x = [p['n']['x'] for p in positions] + [p['s']['x'] for p in positions]
    all_y = [p['n']['y'] for p in positions] + [p['s']['y'] for p in positions]

    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)

    # Add margin
    margin = 20  # mm
    plate_width = max_x - min_x + 2 * margin
    plate_height = max_y - min_y + 2 * margin

    # Create base plate
    front_plate = (cq.Workplane("XY")
                   .rect(plate_width, plate_height)
                   .extrude(thickness)
                   .translate((
                       (min_x + max_x) / 2,
                       (min_y + max_y) / 2,
                       -thickness / 2
                   )))

    # Cut holes for each disc axle (natural and sharp)
    for p in positions:
        axle_r = (p['axle_diameter'] + clearance) / 2

        # Natural disc hole
        front_plate = (front_plate
                       .faces("<Z")
                       .workplane()
                       .moveTo(p['n']['x'], p['n']['y'])
                       .circle(axle_r)
                       .cutThruAll())

        # Sharp disc hole
        front_plate = (front_plate
                       .faces("<Z")
                       .workplane()
                       .moveTo(p['s']['x'], p['s']['y'])
                       .circle(axle_r)
                       .cutThruAll())

    return front_plate


def make_back_plate(config, positions):
    """Create back plate with disc axle engagement holes.

    The back plate is positioned behind the bell crank zone.
    """
    action = config['hardware']['action']
    thickness = action['back_plate']['thickness_mm']
    gap = action['back_plate']['gap_from_front_mm']

    # Calculate bounding box
    all_x = [p['n']['x'] for p in positions] + [p['s']['x'] for p in positions]
    all_y = [p['n']['y'] for p in positions] + [p['s']['y'] for p in positions]

    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)

    margin = 15  # mm
    plate_width = max_x - min_x + 2 * margin
    plate_height = max_y - min_y + 2 * margin

    # Z position: front_plate at z=0, back plate at z = -(front_thickness + gap)
    front_t = config['hardware']['action']['front_plate']['thickness_mm']
    back_z = -(front_t + gap + thickness / 2)

    back_plate = (cq.Workplane("XY")
                  .rect(plate_width, plate_height)
                  .extrude(thickness)
                  .translate((
                      (min_x + max_x) / 2,
                      (min_y + max_y) / 2,
                      back_z
                  )))

    # Cut axle engagement holes (smaller, for bearing surface)
    for p in positions:
        axle_r = p['axle_diameter'] / 2 + 0.05  # Tight fit

        back_plate = (back_plate
                      .faces("<Z")
                      .workplane()
                      .moveTo(p['n']['x'], p['n']['y'])
                      .circle(axle_r)
                      .cutThruAll())

        back_plate = (back_plate
                      .faces("<Z")
                      .workplane()
                      .moveTo(p['s']['x'], p['s']['y'])
                      .circle(axle_r)
                      .cutThruAll())

    return back_plate


def make_plate_section(config, positions, string_range, name):
    """Create a section of the plates for a specific string range.

    Useful for visualizing or manufacturing in sections.
    """
    # Filter positions for this range
    start, end = string_range
    section_pos = [p for p in positions if start <= p['string'] <= end]

    if not section_pos:
        return None

    front = make_front_plate(config, section_pos)
    back = make_back_plate(config, section_pos)

    return front, back


def generate_plate_views():
    """Generate SVG plates showing front and back plate designs."""
    config = load_config()
    positions = compute_disc_positions(config)

    # Create plates for a section (bass register for visualization)
    bass_positions = [p for p in positions if p['string'] <= 9]

    front = make_front_plate(config, bass_positions)
    back = make_back_plate(config, bass_positions)

    # Combine for visualization
    groups = [
        (front, '#4488cc', 0.3),  # Blue - front plate
        (back, '#44cc88', 0.3),   # Green - back plate
    ]

    svg = plate(
        groups,
        UL='XY',   # Top view
        LL='XZ',   # Side view (shows gap between plates)
        UR='ISO',  # Isometric
        LR='YZ',   # End view
        grid='light',
        units='mm'
    )

    output_dir = Path(__file__).parent
    output_path = output_dir / 'action_plates_bass.svg'
    output_path.write_text(svg)

    print(f"Generated {output_path}")
    print(f"Front plate: 6mm thick, {len(bass_positions) * 2} holes")
    print(f"Back plate: 4mm thick, gap 9mm from front")

    return svg


def print_summary():
    """Print summary of plate specifications."""
    config = load_config()
    positions = compute_disc_positions(config)
    action = config['hardware']['action']

    print("=== Action Plate Specifications ===")
    print()
    print(f"FRONT PLATE:")
    print(f"  Thickness: {action['front_plate']['thickness_mm']}mm")
    print(f"  Material: {action['front_plate']['material']}")
    print(f"  Hole clearance: {action['front_plate']['hole_clearance_mm']}mm")
    print(f"  Total holes: {len(positions) * 2} (N + S per string)")
    print()
    print(f"BACK PLATE:")
    print(f"  Thickness: {action['back_plate']['thickness_mm']}mm")
    print(f"  Gap from front: {action['back_plate']['gap_from_front_mm']}mm")
    print(f"  Material: {action['back_plate']['material']}")
    print()
    print(f"BELL CRANK ZONE:")
    print(f"  Clamp length: {action['bell_crank']['clamp_length_mm']}mm")
    print(f"  Arm length: {action['bell_crank']['arm_length_mm']}mm")
    print(f"  Pivot diameter: {action['bell_crank']['pivot_diameter_mm']}mm")
    print()
    print("Hole sizes by register:")
    for reg in ['bass', 'mid', 'treble', 'high']:
        disc = config['hardware']['discs'].get(reg)
        if disc:
            print(f"  {reg.upper()}: Ø{disc['axle_diameter_mm']}mm axle holes")


def main():
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == '--svg':
        generate_plate_views()
    else:
        print_summary()
        print()
        print("Use --svg to generate plate visualization")


if __name__ == '__main__':
    main()
