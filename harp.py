#!/usr/bin/env python3
"""
HARP.PY - Parameterized Harp Design with SVG Output and BOM Generation

A proper parameterized design system for the CLEMENTS47 pedal harp using svgwrite.
Components are modeled as classes with parameters loaded from harp.json.

Usage:
    python3 harp.py                    # Generate SVG and BOM
    python3 harp.py --even-tuning      # Evenly distributed tuning pins
    python3 harp.py --bom-only         # Only print BOM, no SVG

Modules:
    harp_models.py    - Data classes (Disc, String, Harp, etc.)
    harp_physics.py   - String physics and tension calculations
    harp_geometry.py  - Geometry calculations and data loading
    harp_validation.py - Constraint validation
    harp_bom.py       - Bill of materials generation
    harp_renderer.py  - SVG rendering
"""

import math
import argparse

# Import from refactored modules
from harp_models import Harp
from harp_geometry import load_harp_from_json
from harp_validation import validate_constraints, print_constraint_report
from harp_bom import generate_bom, print_bom
from harp_renderer import HarpRenderer
from harp_physics import print_physics_analysis


def main():
    parser = argparse.ArgumentParser(description='Generate harp SVG and BOM')
    parser.add_argument('--json', default='/home/james.clements/projects/harp/harp.json',
                        help='Path to harp JSON file')
    parser.add_argument('--output', default='/home/james.clements/projects/harp/harp.svg',
                        help='Output SVG path')
    parser.add_argument('--even-tuning', action='store_true',
                        help='Distribute tuning pins evenly along neck')
    parser.add_argument('--constant-angle', action='store_true',
                        help='Position tuning pins for constant angle from flat pins')
    parser.add_argument('--target-angle', type=float, default=50.0,
                        help='Target angle in degrees for constant-angle mode (default: 50)')
    parser.add_argument('--bom-only', action='store_true',
                        help='Only print BOM, do not generate SVG')
    parser.add_argument('--scale', type=float, default=0.5,
                        help='SVG scale factor (default: 0.5)')
    parser.add_argument('--neck-offset', type=float, default=80.0,
                        help='Neck X offset in mm for string angle (default: 80mm for ~45-60 deg)')
    parser.add_argument('--neck-thickness', type=float, default=50.0,
                        help='Neck thickness in mm - tuning pins sit on top (default: 50mm)')
    parser.add_argument('--finger-gap', type=float, default=14.0,
                        help='Perpendicular finger gap between strings in mm (default: 14.0)')
    parser.add_argument('--skip-validation', action='store_true',
                        help='Skip constraint validation')
    parser.add_argument('--tolerance-mm', type=float, default=0.5,
                        help='Tolerance for length/gap constraints in mm (default: 0.5)')
    parser.add_argument('--tolerance-deg', type=float, default=1.0,
                        help='Tolerance for angle constraints in degrees (default: 1.0)')
    parser.add_argument('--g7-drop', type=float, default=0.0,
                        help='Additional drop in mm at g7 end of neck (increases downward slope)')
    parser.add_argument('--g7-x-shift', type=float, default=0.0,
                        help='Shift g7 end of neck in X (negative brings g7 tuning pin closer to flat pin)')
    parser.add_argument('--linear-angle', action='store_true',
                        help='Interpolate angles linearly from c1 to g7')
    parser.add_argument('--endpoint-angle', action='store_true',
                        help='Position c1/g7 by angle, distribute others evenly between')
    parser.add_argument('--c1-angle', type=float, default=60.0,
                        help='Target angle for c1 in linear-angle mode (default: 60)')
    parser.add_argument('--g7-angle', type=float, default=45.0,
                        help='Target angle for g7 in linear-angle mode (default: 45)')
    parser.add_argument('--physics', action='store_true',
                        help='Print physics analysis (forces, tensions)')

    args = parser.parse_args()

    # Determine tuning pin mode
    if args.endpoint_angle:
        tuning_mode = "endpoint_angle"
    elif args.linear_angle:
        tuning_mode = "linear_angle"
    elif args.constant_angle:
        tuning_mode = "constant_angle"
    elif args.even_tuning:
        tuning_mode = "even"
    else:
        tuning_mode = "angle_45"

    print(f"Loading harp from {args.json}...")
    print(f"Tuning pin mode: {tuning_mode}")
    if tuning_mode == "constant_angle":
        print(f"Target angle: {args.target_angle}°")
    elif tuning_mode == "linear_angle":
        print(f"Angle range: {args.c1_angle}° (c1) to {args.g7_angle}° (g7)")
    elif tuning_mode == "endpoint_angle":
        print(f"Endpoint angles: c1={args.c1_angle}°, g7={args.g7_angle}°, g7 drop={args.g7_drop}%")
    print(f"Neck X offset: {args.neck_offset} mm")
    print(f"Finger gap: {args.finger_gap} mm (perpendicular)")

    harp = load_harp_from_json(args.json, tuning_pin_mode=tuning_mode,
                               neck_x_offset=args.neck_offset, target_angle=args.target_angle,
                               neck_thickness=args.neck_thickness, finger_gap_mm=args.finger_gap,
                               g7_drop=args.g7_drop, g7_x_shift=args.g7_x_shift,
                               c1_angle=args.c1_angle, g7_angle=args.g7_angle)

    print(f"Loaded {harp.name} with {harp.string_count} strings")

    # Show actual tuning pin positions
    c1_tp = harp.strings[0].tuning_pin
    g7_tp = harp.strings[-1].tuning_pin
    tp_dx = g7_tp.x_mm - c1_tp.x_mm
    tp_dz = g7_tp.z_mm - c1_tp.z_mm
    tp_length = math.sqrt(tp_dx**2 + tp_dz**2)
    print(f"Tuning pins: c1=({c1_tp.x_mm:.1f}, {c1_tp.z_mm:.1f}) to g7=({g7_tp.x_mm:.1f}, {g7_tp.z_mm:.1f})")
    print(f"Tuning pin span: {tp_length:.1f} mm")

    if tuning_mode == "even":
        spacing = harp.neck.length_mm() / (harp.string_count - 1)
        print(f"Tuning pin spacing: {spacing:.2f} mm")

    # Calculate and print string angles (flat pin to tuning pin)
    angles = []
    for s in harp.strings:
        dx = s.tuning_pin.x_mm - s.x_flat_mm
        dz = s.tuning_pin.z_mm - s.z_flat_mm
        angle_deg = math.degrees(math.atan2(dx, dz))  # angle from vertical
        angles.append(angle_deg)

    print(f"String angles (flat to tuning, from vertical):")
    print(f"  c1: {angles[0]:.1f}°, g7: {angles[-1]:.1f}°")
    print(f"  Range: {min(angles):.1f}° to {max(angles):.1f}°")

    # Validate constraints
    if not args.skip_validation:
        violations = validate_constraints(harp,
                                          tolerance_mm=args.tolerance_mm,
                                          tolerance_deg=args.tolerance_deg)
        print_constraint_report(violations, harp)

        if any(v.severity == "error" for v in violations):
            print("WARNING: Design has constraint violations!")
    else:
        print("\n(Constraint validation skipped)")

    # Generate and print BOM
    bom = generate_bom(harp)
    print_bom(bom)

    # Physics analysis (optional)
    if args.physics:
        print_physics_analysis(harp)

    # Generate SVG unless --bom-only
    if not args.bom_only:
        renderer = HarpRenderer(harp, scale=args.scale)

        # Generate three SVGs for different pedal positions
        base_output = args.output.replace('.svg', '')

        # harp0.svg - Flat position (all pedals in neutral)
        renderer.render(f"{base_output}0.svg", pedal_position="flat")

        # harp1.svg - Natural position (natural disc engaged)
        renderer.render(f"{base_output}1.svg", pedal_position="natural")

        # harp2.svg - Sharp position (both discs engaged)
        renderer.render(f"{base_output}2.svg", pedal_position="sharp")

        # harp.svg - Normal view with reaction force vectors
        # Full hardware visible, normal string paths, force vectors overlaid
        renderer.render(args.output, pedal_position="flat",
                       show_force_vectors=True, force_scale=0.05)

        print(f"\nGenerated SVGs:")
        print(f"  {base_output}0.svg - Flat position (strings straight)")
        print(f"  {base_output}1.svg - Natural position (natural disc engaged, string deflected)")
        print(f"  {base_output}2.svg - Sharp position (both discs engaged, string deflected twice)")
        print(f"  {args.output} - Normal view with reaction force vectors")


if __name__ == '__main__':
    main()
