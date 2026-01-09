#!/usr/bin/env python3
"""
HARP_BOM.PY - Bill of Materials generation for harp designs

Contains:
- generate_bom: Create BOM from harp configuration
- print_bom: Print formatted BOM to console
"""

from collections import Counter
from typing import Dict

from harp_models import Harp


def generate_bom(harp: Harp) -> Dict:
    """Generate Bill of Materials from harp configuration."""

    bom = {
        "strings": [],
        "flat_pins": Counter(),
        "tuning_pins": Counter(),
        "discs": {
            "natural": Counter(),
            "sharp": Counter(),
        },
        "prongs": Counter(),
        "summary": {},
    }

    for s in harp.strings:
        # String
        string_spec = f"{s.core_material}"
        if s.wrap_material:
            string_spec += f"/{s.wrap_material}"
        string_spec += f" {s.outer_diameter_mm:.3f}mm"
        bom["strings"].append({
            "number": s.number,
            "note": s.note,
            "spec": string_spec,
            "length_mm": s.length_mm,
        })

        # Flat pin
        fp_key = f"{s.flat_pin.thread} x {s.flat_pin.shaft_diameter_mm}mm"
        bom["flat_pins"][fp_key] += 1

        # Tuning pin
        tp_key = f"{s.tuning_pin.thread} x {s.tuning_pin.diameter_mm}mm"
        bom["tuning_pins"][tp_key] += 1

        # Discs (only natural and sharp for standard pedal harp)
        for disc_type, disc in [
            ("natural", s.natural_disc),
            ("sharp", s.sharp_disc),
        ]:
            if disc:
                disc_key = f"{disc.major_radius_mm:.1f}x{disc.minor_radius_mm:.1f}mm"
                bom["discs"][disc_type][disc_key] += 1

                # Prongs (2 per disc)
                prong_key = f"{disc.prong_diameter_mm:.1f}mm"
                bom["prongs"][prong_key] += 2

    # Summary counts
    bom["summary"] = {
        "total_strings": len(harp.strings),
        "total_flat_pins": sum(bom["flat_pins"].values()),
        "total_tuning_pins": sum(bom["tuning_pins"].values()),
        "total_discs": sum(sum(d.values()) for d in bom["discs"].values()),
        "total_prongs": sum(bom["prongs"].values()),
    }

    return bom


def print_bom(bom: Dict):
    """Print formatted BOM to console."""

    print("\n" + "="*60)
    print("BILL OF MATERIALS - CLEMENTS47 HARP")
    print("="*60)

    print("\n--- STRINGS ---")
    print(f"{'#':>3} {'Note':>4} {'Spec':<30} {'Length':>10}")
    print("-" * 50)
    for s in bom["strings"]:
        print(f"{s['number']:>3} {s['note']:>4} {s['spec']:<30} {s['length_mm']:>8.1f}mm")

    print("\n--- FLAT PINS ---")
    for spec, count in sorted(bom["flat_pins"].items()):
        print(f"  {spec}: {count}")

    print("\n--- TUNING PINS ---")
    for spec, count in sorted(bom["tuning_pins"].items()):
        print(f"  {spec}: {count}")

    print("\n--- DISCS ---")
    for disc_type, discs in bom["discs"].items():
        if discs:
            print(f"  {disc_type.upper()}:")
            for spec, count in sorted(discs.items()):
                print(f"    {spec}: {count}")

    print("\n--- PRONGS ---")
    for spec, count in sorted(bom["prongs"].items()):
        print(f"  {spec} diameter: {count}")

    print("\n--- SUMMARY ---")
    for key, value in bom["summary"].items():
        print(f"  {key.replace('_', ' ').title()}: {value}")

    print("="*60 + "\n")
