#!/usr/bin/env python3
"""
Mechanism Assembly - CadQuery Models

Shows disc + bell crank + plates together to verify fit.
Single string unit assembly for visualization.

Assembly stack (looking from string side, +Z toward viewer):

    STRING (at z=0)
        │
    ┌───┴───┐  prongs (dome tips in scoop grooves)
    │ DISC  │  elliptical body
    └───┬───┘
        │ axle passes through:
    ════╪════  FRONT PLATE (6mm)
        │
      ┌─┴─┐    BELL CRANK (pivots, fork engages axle)
      └─┬─┘
        │
    ════╪════  BACK PLATE (4mm)
        │ axle end (engagement)

"""

import math
from pathlib import Path

import cadquery as cq
from cq_plate import plate

# Import component generators
from disc_assembly import make_disc_assembly, DISC_SPECS
from bell_crank import make_bell_crank, BELL_CRANK_SPECS
from connecting_links import make_connecting_link, LINK_SPECS


def make_front_plate_section(width, height, thickness, axle_hole_d):
    """Create a section of front plate with one axle hole."""
    plate_section = (cq.Workplane("XY")
                     .rect(width, height)
                     .extrude(thickness))

    # Axle hole at center
    plate_section = (plate_section
                     .faces(">Z")
                     .workplane()
                     .circle(axle_hole_d / 2 + 0.1)  # Clearance fit
                     .cutThruAll())

    return plate_section


def make_back_plate_section(width, height, thickness, axle_hole_d):
    """Create a section of back plate with axle engagement hole."""
    plate_section = (cq.Workplane("XY")
                     .rect(width, height)
                     .extrude(thickness))

    # Axle hole (tighter fit for bearing)
    plate_section = (plate_section
                     .faces(">Z")
                     .workplane()
                     .circle(axle_hole_d / 2 + 0.05)
                     .cutThruAll())

    return plate_section


def make_torsion_spring(inner_d=5.0, wire_d=0.5, num_coils=3, leg_length=8.0):
    """Create a simple torsion spring model.

    Spring coils around axle, legs extend to anchor points.
    """
    # Simplified: just show the coil as a torus section and two legs
    coil_r = (inner_d + wire_d) / 2  # Center of wire

    # Create coil as a helix-like shape (simplified as stacked rings)
    coil_height = num_coils * wire_d * 1.5
    coil = (cq.Workplane("XY")
            .circle(coil_r + wire_d/2)
            .circle(coil_r - wire_d/2)
            .extrude(coil_height))

    # Leg 1 - extends radially outward
    leg1 = (cq.Workplane("XY")
            .circle(wire_d / 2)
            .extrude(leg_length)
            .rotate((0, 0, 0), (0, 1, 0), 90)
            .translate((coil_r, 0, coil_height)))

    # Leg 2 - extends opposite direction
    leg2 = (cq.Workplane("XY")
            .circle(wire_d / 2)
            .extrude(leg_length)
            .rotate((0, 0, 0), (0, 1, 0), -90)
            .translate((-coil_r, 0, 0)))

    spring = coil.union(leg1).union(leg2)
    return spring


def make_mechanism_assembly(register='bass'):
    """Create complete single-string mechanism assembly.

    Positions all components relative to string at z=0.
    """
    disc_spec = DISC_SPECS[register]
    crank_spec = BELL_CRANK_SPECS[register]

    # Key dimensions
    sphere_r = disc_spec['sphere_radius']
    disc_t = disc_spec['disc_thickness']
    axle_d = disc_spec['axle_diameter']
    axle_l = disc_spec['axle_length']
    prong_l = disc_spec['prong_length']

    # Plate dimensions (from erand.json action specs)
    front_plate_t = 6.0
    back_plate_t = 4.0
    bell_crank_zone = 8.0
    clearance = 1.0

    # Calculate z positions
    # String at z=0
    # Prong tips at z=0 (touching string)
    # Disc front face at z = -prong_length (since prongs extend +Z from disc face)
    # But we built disc with front face at z=disc_t/2, prongs starting there
    # So we need to translate disc so prong tips are at z=0

    # In disc_assembly, prongs start at disc front face (z=disc_t/2) and extend prong_length
    # So prong tips are at z = disc_t/2 + prong_length
    # To put tips at z=0: translate by -(disc_t/2 + prong_length)
    disc_z_offset = -(disc_t / 2 + prong_l)

    # Axle extends from disc back face in -Z
    # Disc back face is at z = disc_z_offset - disc_t/2
    # Axle end is at z = disc_z_offset - disc_t/2 - axle_l

    # Front plate front face should be at disc back face
    front_plate_z = disc_z_offset - disc_t / 2

    # Back plate position
    back_plate_z = front_plate_z - front_plate_t - clearance - bell_crank_zone - clearance

    # Create components
    disc = make_disc_assembly(disc_spec)
    disc = disc.translate((0, 0, disc_z_offset))

    # Front plate section
    plate_size = disc_spec['major_radius_mm'] * 2 + 20 if 'major_radius_mm' in disc_spec else 60
    front_plate = make_front_plate_section(plate_size, 20, front_plate_t, axle_d)
    front_plate = front_plate.translate((0, 0, front_plate_z - front_plate_t))

    # Back plate section
    back_plate = make_back_plate_section(plate_size, 15, back_plate_t, axle_d)
    back_plate = back_plate.translate((0, 0, back_plate_z - back_plate_t))

    # Bell crank - position in the gap between plates
    crank = make_bell_crank(crank_spec)
    # Rotate to correct orientation (fork toward axle)
    crank = crank.rotate((0, 0, 0), (1, 0, 0), 90)  # Rotate to lie flat
    crank_z = front_plate_z - front_plate_t - clearance - bell_crank_zone / 2
    crank = crank.translate((0, crank_spec['arm_length'] / 2, crank_z))

    # String reference - just a horizontal line marker (not 3D cylinder)
    # Note: String runs perpendicular to viewing angle, represented as small sphere
    string = (cq.Workplane("XY")
              .sphere(1.0))  # Small sphere marker at string position

    # Return spring - coils around axle in bell crank zone
    spring = make_torsion_spring(
        inner_d=axle_d + 1,  # Slightly larger than axle
        wire_d=0.5,
        num_coils=3,
        leg_length=8.0
    )
    # Position in bell crank zone
    spring_z = front_plate_z - front_plate_t - clearance - 2  # Just below front plate
    spring = spring.rotate((0, 0, 0), (1, 0, 0), 90)  # Rotate to align with axle
    spring = spring.translate((0, 0, spring_z))

    return {
        'disc': disc,
        'front_plate': front_plate,
        'back_plate': back_plate,
        'bell_crank': crank,
        'spring': spring,
        'string': string,
    }


def make_assembly_plate(register='bass'):
    """Generate SVG plate showing assembly."""

    components = make_mechanism_assembly(register)

    groups = [
        (components['disc'], '#4488cc', 0.4),        # Blue - disc
        (components['front_plate'], '#888888', 0.3), # Gray - front plate
        (components['back_plate'], '#666666', 0.3),  # Dark gray - back plate
        (components['bell_crank'], '#cc8844', 0.4),  # Brass - bell crank
        (components['spring'], '#44aa44', 0.5),      # Green - return spring
        (components['string'], '#222222', 0.8),      # Black - string
    ]

    svg = plate(
        groups,
        UL='XZ',   # Side view (shows stack-up)
        LL='XY',   # Top view
        UR='ISO',  # Isometric
        LR='YZ',   # End view
        grid='light',
        units='mm'
    )

    return svg


def make_bell_crank_chain(register='bass', num_cranks=3):
    """Create a chain of bell cranks connected by links.

    Shows how motion transfers from one crank to the next via connecting links.

    Args:
        register: 'bass', 'mid', 'treble', or 'high'
        num_cranks: Number of bell cranks in chain (2-4)

    Returns:
        Dict with 'cranks' list and 'links' list
    """
    crank_spec = BELL_CRANK_SPECS[register]
    link_spec = LINK_SPECS[register]

    # Finger gap spacing between strings
    finger_gap = 29.16  # mm from erand.json

    # Bell crank dimensions
    arm_length = crank_spec['arm_length']

    cranks = []
    links = []

    for i in range(num_cranks):
        # Position each crank along X axis (string direction)
        x_pos = i * finger_gap

        crank = make_bell_crank(crank_spec)
        # Rotate so fork points up (+Z) and arm extends in -Y
        crank = crank.rotate((0, 0, 0), (1, 0, 0), 90)
        crank = crank.translate((x_pos, 0, 0))
        cranks.append(crank)

        # Add connecting link between this crank and the next
        if i < num_cranks - 1:
            link = make_connecting_link(link_spec)
            # Position link to connect arm ends of adjacent cranks
            # Link lies flat in XY plane, spanning between cranks
            link_x = x_pos + finger_gap / 2
            link_y = -arm_length  # At end of bell crank arm
            link = link.translate((link_x, link_y, 0))
            links.append(link)

    return {
        'cranks': cranks,
        'links': links
    }


def make_chain_plate(register='bass', num_cranks=3):
    """Generate SVG plate showing bell crank chain with links."""

    chain = make_bell_crank_chain(register, num_cranks)

    groups = []

    # Add all cranks (brass color)
    for crank in chain['cranks']:
        groups.append((crank, '#cc8844', 0.4))

    # Add all links (slightly different brass)
    for link in chain['links']:
        groups.append((link, '#ddaa66', 0.4))

    svg = plate(
        groups,
        UL='XZ',   # Side view
        LL='XY',   # Top view (shows chain layout)
        UR='ISO',  # Isometric
        LR='YZ',   # End view
        grid='light',
        units='mm'
    )

    return svg


def generate_all_assemblies():
    """Generate assembly views for all registers."""

    output_dir = Path(__file__).parent

    print("=== Mechanism Assembly Views ===")
    print()

    for register in ['bass', 'mid', 'treble', 'high']:
        svg = make_assembly_plate(register)

        output_path = output_dir / f"mechanism_assembly_{register}.svg"
        output_path.write_text(svg)

        print(f"{register.upper()}:")
        print(f"  Output: {output_path}")

    print()
    print("Assembly shows: disc + front plate + bell crank + back plate + string")

    # Also generate bell crank chain views
    print()
    print("=== Bell Crank Chain Views ===")
    print()

    for register in ['bass', 'mid', 'treble', 'high']:
        svg = make_chain_plate(register, num_cranks=3)

        output_path = output_dir / f"bell_crank_chain_{register}.svg"
        output_path.write_text(svg)

        print(f"{register.upper()} chain:")
        print(f"  Output: {output_path}")

    print()
    print("Chain shows: 3 bell cranks connected by 2 links")


def main():
    import sys

    if len(sys.argv) > 1 and sys.argv[1] in DISC_SPECS:
        register = sys.argv[1]
        svg = make_assembly_plate(register)
        output = Path(__file__).parent / f"mechanism_assembly_{register}.svg"
        output.write_text(svg)
        print(f"Generated {output}")
    else:
        generate_all_assemblies()


if __name__ == '__main__':
    main()
