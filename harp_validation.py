#!/usr/bin/env python3
"""
HARP_VALIDATION.PY - Constraint validation for harp designs

Contains:
- ConstraintViolation: Data class for constraint violations
- validate_constraints: Check all design constraints
- print_constraint_report: Print formatted validation report
"""

import math
from dataclasses import dataclass
from typing import List, Optional

from harp_models import Harp


@dataclass
class ConstraintViolation:
    """A constraint violation found during validation."""
    constraint: str
    message: str
    severity: str = "error"  # "error" or "warning"
    string_number: Optional[int] = None
    actual_value: Optional[float] = None
    expected_value: Optional[float] = None
    tolerance: Optional[float] = None


def validate_constraints(harp: Harp, tolerance_mm: float = 0.5, tolerance_deg: float = 1.0) -> List[ConstraintViolation]:
    """
    Validate harp design constraints:
    1. All strings must remain parallel (same tilt angle)
    2. Finger gap must be maintained (14mm perpendicular gap)
    3. Flat string lengths cannot change from designed values

    Returns list of violations (empty if all constraints pass).
    """
    violations = []

    # Get design parameters from metadata
    design_finger_gap = harp.metadata.get('finger_gap_perpendicular_mm', 14.0)
    design_tilt_deg = harp.metadata.get('string_tilt_degrees', -3.0)

    # ---------------------------------------------------------------------
    # 1. STRING PARALLELISM CHECK
    # All strings should have the same tilt angle (soundboard to flat pin)
    # ---------------------------------------------------------------------
    string_angles = []
    for s in harp.strings:
        # Calculate angle from vertical (soundboard to flat pin)
        dx = s.x_flat_mm - s.x_soundboard_mm
        dz = s.z_flat_mm - s.z_soundboard_mm
        angle_deg = math.degrees(math.atan2(dx, dz))
        string_angles.append((s.number, angle_deg))

    if string_angles:
        # All strings should have same angle (within tolerance)
        angles_only = [a[1] for a in string_angles]
        mean_angle = sum(angles_only) / len(angles_only)

        for string_num, angle in string_angles:
            deviation = abs(angle - mean_angle)
            if deviation > tolerance_deg:
                violations.append(ConstraintViolation(
                    constraint="string_parallelism",
                    message=f"String {string_num} not parallel: {angle:.2f}° (mean: {mean_angle:.2f}°)",
                    string_number=string_num,
                    actual_value=angle,
                    expected_value=mean_angle,
                    tolerance=tolerance_deg
                ))

        # Check against design tilt
        if abs(mean_angle - design_tilt_deg) > tolerance_deg:
            violations.append(ConstraintViolation(
                constraint="string_tilt",
                message=f"Mean string tilt {mean_angle:.2f}° differs from design {design_tilt_deg:.2f}°",
                severity="warning",
                actual_value=mean_angle,
                expected_value=design_tilt_deg,
                tolerance=tolerance_deg
            ))

    # ---------------------------------------------------------------------
    # 2. FINGER GAP CHECK
    # Perpendicular distance between strings on the SAME plate must be 14mm
    # Adjacent strings (1-2, 2-3) are on OPPOSITE plates (+Y/-Y)
    # So we check strings 1-3, 2-4, 3-5, etc. (same plate pairs)
    # ---------------------------------------------------------------------
    for i in range(len(harp.strings) - 2):
        s1 = harp.strings[i]
        s2 = harp.strings[i + 2]  # Skip one string to get same plate

        # Verify they're on the same plate
        if s1.natural_disc.plate != s2.natural_disc.plate:
            violations.append(ConstraintViolation(
                constraint="plate_assignment",
                message=f"Strings {s1.number} and {s2.number} should be on same plate but aren't",
                severity="error",
                string_number=s1.number,
            ))
            continue

        # Calculate perpendicular gap at flat pin height
        dx = s2.x_flat_mm - s1.x_flat_mm
        dz = s2.z_flat_mm - s1.z_flat_mm

        # Get string direction vector (normalized)
        s_dx = s1.x_flat_mm - s1.x_soundboard_mm
        s_dz = s1.z_flat_mm - s1.z_soundboard_mm
        s_len = math.sqrt(s_dx*s_dx + s_dz*s_dz)
        if s_len > 0:
            s_dx /= s_len
            s_dz /= s_len

        # Perpendicular distance = |displacement × string_direction|
        # For 2D: displacement.x * dir.z - displacement.y * dir.x
        perp_gap = abs(dx * s_dz - dz * s_dx)

        # For same-plate pairs (skipping 1 string), the design gap is 2x finger gap
        # because there's one string in between on the other plate
        design_same_plate_gap = design_finger_gap * 2

        gap_error = abs(perp_gap - design_same_plate_gap)
        if gap_error > tolerance_mm:
            violations.append(ConstraintViolation(
                constraint="finger_gap_same_plate",
                message=f"Gap between strings {s1.number}-{s2.number} (plate {s1.natural_disc.plate}): {perp_gap:.2f}mm (design: {design_same_plate_gap}mm)",
                string_number=s1.number,
                actual_value=perp_gap,
                expected_value=design_same_plate_gap,
                tolerance=tolerance_mm
            ))

    # Also check adjacent string gaps (opposite plates) for reference
    for i in range(len(harp.strings) - 1):
        s1 = harp.strings[i]
        s2 = harp.strings[i + 1]

        dx = s2.x_flat_mm - s1.x_flat_mm
        dz = s2.z_flat_mm - s1.z_flat_mm

        s_dx = s1.x_flat_mm - s1.x_soundboard_mm
        s_dz = s1.z_flat_mm - s1.z_soundboard_mm
        s_len = math.sqrt(s_dx*s_dx + s_dz*s_dz)
        if s_len > 0:
            s_dx /= s_len
            s_dz /= s_len

        perp_gap = abs(dx * s_dz - dz * s_dx)

        gap_error = abs(perp_gap - design_finger_gap)
        if gap_error > tolerance_mm:
            violations.append(ConstraintViolation(
                constraint="finger_gap_adjacent",
                message=f"Gap between strings {s1.number}-{s2.number}: {perp_gap:.2f}mm (design: {design_finger_gap}mm)",
                string_number=s1.number,
                actual_value=perp_gap,
                expected_value=design_finger_gap,
                tolerance=tolerance_mm,
                severity="warning"  # Warning since adjacent strings are on opposite plates
            ))

    # ---------------------------------------------------------------------
    # 3. FLAT STRING LENGTH CHECK
    # The designed string length must match the geometry
    # ---------------------------------------------------------------------
    for s in harp.strings:
        # Calculate geometric length (soundboard to flat pin)
        dx = s.x_flat_mm - s.x_soundboard_mm
        dz = s.z_flat_mm - s.z_soundboard_mm
        calc_length = math.sqrt(dx*dx + dz*dz)

        # Compare to designed length
        length_error = abs(calc_length - s.length_mm)
        if length_error > tolerance_mm:
            violations.append(ConstraintViolation(
                constraint="flat_string_length",
                message=f"String {s.number} ({s.note}): calculated {calc_length:.2f}mm vs design {s.length_mm:.2f}mm",
                string_number=s.number,
                actual_value=calc_length,
                expected_value=s.length_mm,
                tolerance=tolerance_mm
            ))

    # ---------------------------------------------------------------------
    # 4. NECK CLEARANCE CHECK
    # Neck bottom line must be 30mm above the highest flat pin
    # ---------------------------------------------------------------------
    design_neck_clearance = harp.neck.clearance_mm  # Default 30mm
    highest_flat_z = max(s.z_flat_mm for s in harp.strings)

    # Neck bottom at the X position of the highest flat pin
    # Find which string has the highest flat pin
    highest_flat_string = max(harp.strings, key=lambda s: s.z_flat_mm)
    highest_flat_x = highest_flat_string.x_flat_mm

    # Interpolate neck Z at that X position
    neck = harp.neck
    if abs(neck.g7_x_mm - neck.c1_x_mm) > 0.001:
        t = (highest_flat_x - neck.c1_x_mm) / (neck.g7_x_mm - neck.c1_x_mm)
        t = max(0, min(1, t))  # Clamp to [0,1]
        neck_z_at_highest_flat = neck.c1_z_mm + t * (neck.g7_z_mm - neck.c1_z_mm)
    else:
        neck_z_at_highest_flat = neck.c1_z_mm

    actual_clearance = neck_z_at_highest_flat - highest_flat_z

    if actual_clearance < design_neck_clearance - tolerance_mm:
        violations.append(ConstraintViolation(
            constraint="neck_clearance",
            message=f"Neck clearance {actual_clearance:.1f}mm < required {design_neck_clearance}mm above highest flat pin (string {highest_flat_string.number})",
            severity="error",
            actual_value=actual_clearance,
            expected_value=design_neck_clearance,
            tolerance=tolerance_mm
        ))

    return violations


def print_constraint_report(violations: List[ConstraintViolation], harp: Harp):
    """Print a formatted constraint validation report."""

    print("\n" + "="*60)
    print("CONSTRAINT VALIDATION REPORT")
    print("="*60)

    # Design parameters
    design_finger_gap = harp.metadata.get('finger_gap_perpendicular_mm', 14.0)
    design_tilt_deg = harp.metadata.get('string_tilt_degrees', -3.0)

    print(f"\nDesign Parameters:")
    print(f"  Finger gap (perpendicular): {design_finger_gap} mm")
    print(f"  String tilt from vertical: {design_tilt_deg}°")
    print(f"  Constraint: {harp.metadata.get('constraint', 'N/A')}")

    if not violations:
        print("\n✓ All constraints PASSED")
        print("="*60 + "\n")
        return

    # Group by constraint type
    by_constraint = {}
    for v in violations:
        if v.constraint not in by_constraint:
            by_constraint[v.constraint] = []
        by_constraint[v.constraint].append(v)

    errors = [v for v in violations if v.severity == "error"]
    warnings = [v for v in violations if v.severity == "warning"]

    print(f"\n✗ Found {len(errors)} errors, {len(warnings)} warnings")

    for constraint, vlist in by_constraint.items():
        print(f"\n--- {constraint.upper().replace('_', ' ')} ---")
        for v in vlist[:5]:  # Limit to first 5 per category
            marker = "✗" if v.severity == "error" else "⚠"
            print(f"  {marker} {v.message}")
        if len(vlist) > 5:
            print(f"  ... and {len(vlist) - 5} more")

    print("="*60 + "\n")
