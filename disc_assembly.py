#!/usr/bin/env python3
"""
Disc Assembly - CadQuery Models

Elliptical disc body with prongs at major axis ends.
Lathe-turned from 316SS rod.

Prong features:
- Dome tip (hemisphere) - no sharp edges
- Scoop groove - smooth string seat
- Prongs at ends of ellipse major axis (accounts for prong diameter)

Top view (XY plane, looking down axle):

          ⌒ prong (dome + scoop)
         ╱ ╲
        ╱   ╲
       (     ) ← elliptical disc body
        ╲   ╱
         ╲ ╱
          ⌒ prong
           │
         axle (into page, -Z)

Side view (XZ plane):

    prong ──┬── prong
            │
         ═══════  ← ellipse disc (thin)
            │
          axle
"""

import math
from pathlib import Path

import cadquery as cq
from cq_plate import plate


AXLE_LENGTH = 20.0  # mm

DISC_SPECS = {
    'bass': {
        'strings': '1-9',
        'sphere_radius': 23,        # mm - vibration clearance
        'prong_diameter': 3.0,      # mm
        'prong_length': 23.0,       # mm
        'axle_diameter': 4.0,       # mm
        'axle_length': AXLE_LENGTH,
        'disc_thickness': 3.0,      # mm - ellipse extrusion height
        'minor_radius': 5.0,        # mm - ellipse minor radius (structural)
        'scoop_radius': 1.0,        # mm
    },
    'mid': {
        'strings': '10-28',
        'sphere_radius': 20,
        'prong_diameter': 2.5,
        'prong_length': 20.0,
        'axle_diameter': 3.0,
        'axle_length': AXLE_LENGTH,
        'disc_thickness': 2.5,
        'minor_radius': 4.0,
        'scoop_radius': 0.8,
    },
    'treble': {
        'strings': '29-38',
        'sphere_radius': 6,
        'prong_diameter': 2.0,
        'prong_length': 6.0,
        'axle_diameter': 2.0,
        'axle_length': AXLE_LENGTH,
        'disc_thickness': 2.0,
        'minor_radius': 3.0,
        'scoop_radius': 0.5,
    },
    'high': {
        'strings': '39-47',
        'sphere_radius': 3,
        'prong_diameter': 1.5,
        'prong_length': 3.0,
        'axle_diameter': 1.5,
        'axle_length': AXLE_LENGTH,
        'disc_thickness': 1.5,
        'minor_radius': 2.0,
        'scoop_radius': 0.4,
    },
}


def make_prong(diameter, length, scoop_radius):
    """Create prong with dome tip and scoop groove.

    Prong extends along +Z from origin.
    """
    r = diameter / 2

    # Groove parameters
    groove_z = length - diameter
    groove_width = scoop_radius * 2
    groove_depth = min(scoop_radius * 0.7, r * 0.35)

    if groove_z < r:
        groove_z = r

    # Base cylinder
    base_height = max(0.1, groove_z - groove_width / 2)
    base = (cq.Workplane("XY")
            .circle(r)
            .extrude(base_height))

    # Groove section (reduced diameter)
    groove_r = r - groove_depth
    groove = (cq.Workplane("XY")
              .workplane(offset=base_height)
              .circle(groove_r)
              .extrude(groove_width))

    # Upper section
    upper_start = base_height + groove_width
    upper_height = max(0.1, length - upper_start - r)
    upper = (cq.Workplane("XY")
             .workplane(offset=upper_start)
             .circle(r)
             .extrude(upper_height))

    # Dome tip (hemisphere)
    dome_z = upper_start + upper_height
    dome = (cq.Workplane("XY")
            .workplane(offset=dome_z)
            .sphere(r)
            .cut(
                cq.Workplane("XY")
                .box(r * 3, r * 3, r * 2)
                .translate((0, 0, dome_z - r))
            ))

    prong = base.union(groove).union(upper).union(dome)

    # Fillet groove transitions
    try:
        prong = prong.edges("|Z").edges(
            cq.selectors.BoxSelector(
                (-r * 2, -r * 2, base_height - 0.1),
                (r * 2, r * 2, upper_start + 0.1)
            )
        ).fillet(min(groove_depth * 0.8, scoop_radius * 0.5))
    except:
        pass

    return prong


def make_disc_assembly(spec, include_engagement_pin=True):
    """Create elliptical disc with prongs at major axis ends.

    Geometry:
    - Ellipse in XY plane, centered at origin
    - Major axis along X (through prong centers)
    - Minor axis along Y
    - Prongs extend in +Z direction from ellipse ends
    - Axle extends in -Z direction from center
    - Engagement pin on axle for bell crank fork

    Major radius = sphere_radius + prong_diameter (prong outer edge)
    """
    prong_d = spec['prong_diameter']
    prong_l = spec['prong_length']
    prong_r = prong_d / 2
    axle_d = spec['axle_diameter']
    axle_l = spec['axle_length']
    disc_t = spec['disc_thickness']
    minor_r = spec['minor_radius']
    scoop_r = spec['scoop_radius']
    sphere_r = spec['sphere_radius']

    # Prong positioning:
    # - Prong inner edge at sphere_radius (vibration clearance)
    # - Prong center at sphere_radius + prong_radius
    # - Prong outer edge at sphere_radius + prong_diameter
    prong_center = sphere_r + prong_r
    major_r = sphere_r + prong_d  # Extends to prong outer edge

    # Ensure minor radius is at least enough for axle
    if minor_r < axle_d:
        minor_r = axle_d + 1

    # Create elliptical disc body
    disc = (cq.Workplane("XY")
            .ellipse(major_r, minor_r)
            .extrude(disc_t)
            .translate((0, 0, -disc_t / 2)))

    # Add axle extending from back face (-Z)
    axle = (cq.Workplane("XY")
            .circle(axle_d / 2)
            .extrude(axle_l)
            .translate((0, 0, -disc_t / 2 - axle_l)))

    # Add prongs at major axis ends (+X and -X)
    right_prong = make_prong(prong_d, prong_l, scoop_r)
    right_prong = right_prong.translate((prong_center, 0, disc_t / 2))

    left_prong = make_prong(prong_d, prong_l, scoop_r)
    left_prong = left_prong.translate((-prong_center, 0, disc_t / 2))

    # Combine base assembly
    assembly = disc.union(axle).union(right_prong).union(left_prong)

    # Add engagement pin for bell crank fork
    # Position: in bell crank zone (front_plate + clearance + half bell_crank_zone)
    # Stack: front_plate(6) + clearance(1) + bell_crank_zone(8)/2 = 11mm from disc back
    if include_engagement_pin:
        front_plate_t = 6.0
        clearance = 1.0
        bell_crank_zone = 8.0
        pin_z_offset = front_plate_t + clearance + bell_crank_zone / 2  # 11mm

        # Pin dimensions
        pin_d = axle_d * 0.8  # Slightly smaller than axle
        pin_length = axle_d * 2  # Extends both sides of axle

        # Pin center position (on axle, in bell crank zone)
        pin_z = -disc_t / 2 - pin_z_offset

        # Create pin perpendicular to axle (along Y axis)
        engagement_pin = (cq.Workplane("XZ")
                          .circle(pin_d / 2)
                          .extrude(pin_length)
                          .translate((0, -pin_length / 2, pin_z)))

        assembly = assembly.union(engagement_pin)

    return assembly


def make_dimensioned_views(name, spec):
    """Create multi-view plate for manufacturing.

    Note: Engagement pin is omitted from manufacturing views since it's
    part of the mechanism assembly context, not the turned disc itself.
    """
    assembly = make_disc_assembly(spec, include_engagement_pin=False)

    groups = [(assembly, '#4488cc', 0.3)]

    svg = plate(
        groups,
        UL='XZ',   # Side view (prong length, disc, axle)
        LL='XY',   # Top view (ellipse shape)
        UR='ISO',  # Isometric
        LR='YZ',   # End view
        grid='light',
        units='mm'
    )

    return svg


def generate_all_plates():
    """Generate manufacturing plates for all disc sizes."""

    output_dir = Path(__file__).parent

    print("=== Elliptical Disc Assembly ===")
    print()

    for name, spec in DISC_SPECS.items():
        svg = make_dimensioned_views(name, spec)

        output_path = output_dir / f"disc_{name}.svg"
        output_path.write_text(svg)

        prong_d = spec['prong_diameter']
        major_r = spec['sphere_radius'] + prong_d  # To prong outer edge
        inner_gap = spec['sphere_radius'] * 2

        print(f"{name.upper()} (strings {spec['strings']}):")
        print(f"  Ellipse: {major_r * 2:.1f}mm x {spec['minor_radius'] * 2:.1f}mm x {spec['disc_thickness']}mm")
        print(f"  Prongs: Ø{spec['prong_diameter']}mm x {spec['prong_length']}mm (dome tip)")
        print(f"  Prong inner gap: {inner_gap}mm")
        print(f"  Axle: Ø{spec['axle_diameter']}mm x {spec['axle_length']}mm")
        print(f"  Scoop: R{spec['scoop_radius']}mm")
        print(f"  Output: {output_path}")
        print()


def main():
    import sys

    if len(sys.argv) > 1 and sys.argv[1] in DISC_SPECS:
        name = sys.argv[1]
        spec = DISC_SPECS[name]
        svg = make_dimensioned_views(name, spec)
        output = Path(__file__).parent / f"disc_{name}.svg"
        output.write_text(svg)
        print(f"Generated {output}")
    else:
        generate_all_plates()


if __name__ == '__main__':
    main()
