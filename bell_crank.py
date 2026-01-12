#!/usr/bin/env python3
"""
Bell Crank Assembly - CadQuery Models

Erard-style yoke bell crank that converts action rod linear motion
to disc rotational motion. Based on analysis of Fusion 360 reference models.

Manufacturing: Brass, machined or investment cast
"""

import math
from pathlib import Path

import cadquery as cq
from cq_plate import plate


# Bell crank specifications by register (matching disc sizes)
BELL_CRANK_SPECS = {
    'bass': {
        'strings': '1-9',
        'pivot_diameter': 3.0,      # mm - pivot pin hole
        'arm_length': 18.0,         # mm - action rod connection arm
        'fork_length': 12.0,        # mm - fork arms that wrap disc axle
        'fork_gap': 5.0,            # mm - gap between fork tines (fits Ø4 axle)
        'fork_width': 3.0,          # mm - width of each fork tine
        'body_thickness': 4.0,      # mm - main body thickness
        'body_width': 8.0,          # mm - body width at pivot
        'action_rod_hole': 2.5,     # mm - hole for action rod pin
    },
    'mid': {
        'strings': '10-28',
        'pivot_diameter': 3.0,
        'arm_length': 15.0,
        'fork_length': 10.0,
        'fork_gap': 4.0,            # fits Ø3 axle
        'fork_width': 2.5,
        'body_thickness': 3.5,
        'body_width': 7.0,
        'action_rod_hole': 2.0,
    },
    'treble': {
        'strings': '29-38',
        'pivot_diameter': 2.5,
        'arm_length': 12.0,
        'fork_length': 8.0,
        'fork_gap': 3.0,            # fits Ø2 axle
        'fork_width': 2.0,
        'body_thickness': 3.0,
        'body_width': 6.0,
        'action_rod_hole': 1.5,
    },
    'high': {
        'strings': '39-47',
        'pivot_diameter': 2.0,
        'arm_length': 10.0,
        'fork_length': 6.0,
        'fork_gap': 2.5,            # fits Ø1.5 axle
        'fork_width': 1.5,
        'body_thickness': 2.5,
        'body_width': 5.0,
        'action_rod_hole': 1.5,
    },
}


def make_bell_crank(spec):
    r"""Create Erard-style yoke bell crank.

    Geometry (top view, pivot at origin):

           fork_length
           <-------->
           +--------+   ^
          /          \  | fork_width
         |   (axle)   | v
          \          /
           +--+  +--+
              |  |      <- fork_gap (fits disc axle)
           +--+  +--+
          /          \
         |   (axle)   |
          \          /
           +--------+
              |
              |  <- body extends to pivot
           [PIVOT]  <- pivot_diameter hole
              |
              |  <- arm extends to action rod
              O  <- action_rod_hole

    The bell crank pivots around the center hole. When the action rod
    pushes/pulls, the fork rotates the disc via the axle.
    """
    pivot_d = spec['pivot_diameter']
    arm_len = spec['arm_length']
    fork_len = spec['fork_length']
    fork_gap = spec['fork_gap']
    fork_w = spec['fork_width']
    body_t = spec['body_thickness']
    body_w = spec['body_width']
    rod_hole = spec['action_rod_hole']

    # Build the bell crank as a 2D profile extruded
    # Pivot at origin, fork extends in +Y, arm extends in -Y

    # Main body rectangle from pivot to arm end
    body_half_w = body_w / 2

    # Create the main body profile
    # Start with the arm (extends in -Y direction)
    arm_end_y = -arm_len

    # Fork base starts at pivot and extends in +Y
    fork_base_y = body_w / 2  # Small gap above pivot
    fork_end_y = fork_base_y + fork_len

    # Fork tine positions (centered on Y axis with gap between)
    tine_inner = fork_gap / 2
    tine_outer = tine_inner + fork_w

    # Build as a series of operations
    # Main body: rectangle from arm end to fork base
    body = (cq.Workplane("XY")
            .rect(body_w, arm_len + body_w)
            .extrude(body_t)
            .translate((0, (body_w - arm_len) / 2, 0)))

    # Fork tines: two rectangles extending from body
    # Left tine
    left_tine = (cq.Workplane("XY")
                 .rect(fork_w, fork_len)
                 .extrude(body_t)
                 .translate((-tine_inner - fork_w/2, fork_base_y + fork_len/2, 0)))

    # Right tine
    right_tine = (cq.Workplane("XY")
                  .rect(fork_w, fork_len)
                  .extrude(body_t)
                  .translate((tine_inner + fork_w/2, fork_base_y + fork_len/2, 0)))

    # Fork bridge at the end (connects the two tines)
    bridge_width = fork_gap + 2 * fork_w
    bridge = (cq.Workplane("XY")
              .rect(bridge_width, fork_w)
              .extrude(body_t)
              .translate((0, fork_end_y - fork_w/2, 0)))

    # Combine all parts
    crank = body.union(left_tine).union(right_tine).union(bridge)

    # Cut the pivot hole at origin
    crank = (crank
             .faces(">Z")
             .workplane()
             .moveTo(0, 0)
             .circle(pivot_d / 2)
             .cutThruAll())

    # Cut the action rod connection hole at arm end
    rod_hole_y = arm_end_y + rod_hole + 2  # 2mm from end
    crank = (crank
             .faces(">Z")
             .workplane()
             .moveTo(0, rod_hole_y)
             .circle(rod_hole / 2)
             .cutThruAll())

    # Add rounded corners/fillets for strength
    # (Skip if geometry is too small for fillets)
    try:
        crank = crank.edges("|Z").fillet(min(1.0, fork_w * 0.3))
    except:
        pass  # Fillet failed, keep sharp corners

    return crank


def make_dimensioned_plate(name, spec):
    """Create manufacturing plate with dimensions."""

    crank = make_bell_crank(spec)

    groups = [(crank, '#cc8844', 0.4)]  # Brass color

    svg = plate(
        groups,
        UL='XY',   # Top view (shows fork shape)
        LL='XZ',   # Side view (shows thickness)
        UR='ISO',  # Isometric
        LR='YZ',   # End view
        grid='light',
        units='mm'
    )

    return svg


def generate_all_plates():
    """Generate manufacturing plates for all bell crank sizes."""

    output_dir = Path(__file__).parent

    print("=== Bell Crank Manufacturing Plates ===")
    print()

    for name, spec in BELL_CRANK_SPECS.items():
        svg = make_dimensioned_plate(name, spec)

        output_path = output_dir / f"bell_crank_{name}.svg"
        output_path.write_text(svg)

        print(f"{name.upper()} (strings {spec['strings']}):")
        print(f"  Pivot: Ø{spec['pivot_diameter']}mm")
        print(f"  Arm length: {spec['arm_length']}mm")
        print(f"  Fork length: {spec['fork_length']}mm, gap: {spec['fork_gap']}mm")
        print(f"  Thickness: {spec['body_thickness']}mm")
        print(f"  Output: {output_path}")
        print()


def main():
    import sys

    if len(sys.argv) > 1 and sys.argv[1] in BELL_CRANK_SPECS:
        name = sys.argv[1]
        spec = BELL_CRANK_SPECS[name]
        svg = make_dimensioned_plate(name, spec)
        output = Path(__file__).parent / f"bell_crank_{name}.svg"
        output.write_text(svg)
        print(f"Generated {output}")
    else:
        generate_all_plates()


if __name__ == '__main__':
    main()
