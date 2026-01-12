#!/usr/bin/env python3
"""
Export Models - STEP and STL Export for Manufacturing

Exports disc assemblies and other components to industry-standard formats:
- STEP (.step) - for CNC/CAM and professional CAD interchange
- STL (.stl) - for 3D printing

Usage:
    python export_models.py                    # Export all disc sizes
    python export_models.py bass               # Export single size
    python export_models.py --format step      # STEP only
    python export_models.py --format stl       # STL only
    python export_models.py --output-dir ./out # Custom output directory
"""

import sys
from pathlib import Path

import cadquery as cq

from disc_assembly import make_disc_assembly, DISC_SPECS
from bell_crank import make_bell_crank, BELL_CRANK_SPECS
from connecting_links import make_connecting_link, LINK_SPECS


def export_step(solid, filepath):
    """Export CadQuery solid to STEP format."""
    cq.exporters.export(solid, str(filepath), exportType='STEP')


def export_stl(solid, filepath, tolerance=0.01, angular_tolerance=0.1):
    """Export CadQuery solid to STL format.

    Args:
        solid: CadQuery solid to export
        filepath: Output file path
        tolerance: Linear tolerance for mesh (smaller = finer mesh)
        angular_tolerance: Angular tolerance in radians
    """
    cq.exporters.export(
        solid, str(filepath), exportType='STL',
        tolerance=tolerance, angularTolerance=angular_tolerance
    )


def export_disc_assembly(register, output_dir, formats=('step', 'stl')):
    """Export a disc assembly to specified formats.

    Args:
        register: 'bass', 'mid', 'treble', or 'high'
        output_dir: Directory for output files
        formats: Tuple of formats to export ('step', 'stl')

    Returns:
        List of exported file paths
    """
    spec = DISC_SPECS[register]
    disc = make_disc_assembly(spec, include_engagement_pin=True)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    exported = []

    if 'step' in formats:
        step_path = output_dir / f"disc_{register}.step"
        export_step(disc, step_path)
        exported.append(step_path)
        print(f"  STEP: {step_path}")

    if 'stl' in formats:
        stl_path = output_dir / f"disc_{register}.stl"
        export_stl(disc, stl_path)
        exported.append(stl_path)
        print(f"  STL:  {stl_path}")

    return exported


def export_bell_crank(register, output_dir, formats=('step', 'stl')):
    """Export a bell crank to specified formats."""
    spec = BELL_CRANK_SPECS[register]
    crank = make_bell_crank(spec)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    exported = []

    if 'step' in formats:
        step_path = output_dir / f"bell_crank_{register}.step"
        export_step(crank, step_path)
        exported.append(step_path)
        print(f"  STEP: {step_path}")

    if 'stl' in formats:
        stl_path = output_dir / f"bell_crank_{register}.stl"
        export_stl(crank, stl_path)
        exported.append(stl_path)
        print(f"  STL:  {stl_path}")

    return exported


def export_connecting_link(register, output_dir, formats=('step', 'stl')):
    """Export a connecting link to specified formats."""
    spec = LINK_SPECS[register]
    link = make_connecting_link(spec)

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    exported = []

    if 'step' in formats:
        step_path = output_dir / f"link_{register}.step"
        export_step(link, step_path)
        exported.append(step_path)
        print(f"  STEP: {step_path}")

    if 'stl' in formats:
        stl_path = output_dir / f"link_{register}.stl"
        export_stl(link, stl_path)
        exported.append(stl_path)
        print(f"  STL:  {stl_path}")

    return exported


def export_all(output_dir=None, formats=('step', 'stl')):
    """Export all disc assemblies, bell cranks, and connecting links.

    Args:
        output_dir: Output directory (default: ./exports)
        formats: Tuple of formats to export
    """
    if output_dir is None:
        output_dir = Path(__file__).parent / "exports"

    output_dir = Path(output_dir)

    print("=== Exporting 3D Models ===")
    print(f"Output directory: {output_dir}")
    print(f"Formats: {', '.join(formats).upper()}")
    print()

    all_exported = []

    # Export disc assemblies
    print("Disc Assemblies:")
    for register in DISC_SPECS.keys():
        print(f"  {register.upper()}:")
        exported = export_disc_assembly(register, output_dir / "discs", formats)
        all_exported.extend(exported)
    print()

    # Export bell cranks
    print("Bell Cranks:")
    for register in BELL_CRANK_SPECS.keys():
        print(f"  {register.upper()}:")
        exported = export_bell_crank(register, output_dir / "bell_cranks", formats)
        all_exported.extend(exported)
    print()

    # Export connecting links
    print("Connecting Links:")
    for register in LINK_SPECS.keys():
        print(f"  {register.upper()}:")
        exported = export_connecting_link(register, output_dir / "links", formats)
        all_exported.extend(exported)
    print()

    print(f"Total files exported: {len(all_exported)}")
    return all_exported


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Export disc assemblies and components to STEP/STL formats'
    )
    parser.add_argument(
        'register', nargs='?',
        choices=['bass', 'mid', 'treble', 'high', 'all'],
        default='all',
        help='Register to export (default: all)'
    )
    parser.add_argument(
        '--format', '-f',
        choices=['step', 'stl', 'both'],
        default='both',
        help='Export format (default: both)'
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        default=None,
        help='Output directory (default: ./exports)'
    )
    parser.add_argument(
        '--component', '-c',
        choices=['disc', 'bell_crank', 'all'],
        default='all',
        help='Component type to export (default: all)'
    )

    args = parser.parse_args()

    # Determine formats
    if args.format == 'both':
        formats = ('step', 'stl')
    else:
        formats = (args.format,)

    # Determine output directory
    output_dir = args.output_dir or Path(__file__).parent / "exports"

    if args.register == 'all':
        export_all(output_dir, formats)
    else:
        print(f"=== Exporting {args.register.upper()} ===")
        print(f"Formats: {', '.join(formats).upper()}")
        print()

        if args.component in ('disc', 'all'):
            print(f"Disc Assembly ({args.register}):")
            export_disc_assembly(args.register, output_dir / "discs", formats)
            print()

        if args.component in ('bell_crank', 'all'):
            print(f"Bell Crank ({args.register}):")
            export_bell_crank(args.register, output_dir / "bell_cranks", formats)
            print()


if __name__ == '__main__':
    main()
