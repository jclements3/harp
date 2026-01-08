#!/usr/bin/env python3
"""
HARP.PY - Parameterized Harp Design with SVG Output and BOM Generation

A proper parameterized design system for the CLEMENTS47 pedal harp using svgwrite.
Components are modeled as classes with parameters loaded from harp.json.

Usage:
    python3 harp.py                    # Generate SVG and BOM
    python3 harp.py --even-tuning      # Evenly distributed tuning pins
    python3 harp.py --bom-only         # Only print BOM, no SVG
"""

import json
import math
import argparse
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
from collections import Counter

import svgwrite
from svgwrite import mm, px


# =============================================================================
# COMPONENT CLASSES
# =============================================================================

@dataclass
class Prong:
    """A prong on a disc that engages the string."""
    x_mm: float
    z_mm: float
    diameter_mm: float
    material: str = "Stainless Steel 316"

    def radius_mm(self) -> float:
        return self.diameter_mm / 2


@dataclass
class Disc:
    """A semitone disc (natural, sharp, double-flat, or double-sharp)."""
    x_mm: float
    z_mm: float
    y_mm: float  # +10 for +Y plate, -10 for -Y plate
    major_radius_mm: float
    minor_radius_mm: float
    disc_type: str  # "natural", "sharp", "double_flat", "double_sharp"
    prong_diameter_mm: float
    plate: str = "+Y"  # "+Y" or "-Y" neck plate
    material: str = "Stainless Steel 316"
    thickness_mm: float = 4.0
    rotation_degrees: float = 0.0  # 0=flat, 45=natural engaged, 90=sharp engaged
    string_ux: float = 0.0  # String direction unit vector X component
    string_uz: float = 1.0  # String direction unit vector Z component

    def prong_inset_distance(self) -> float:
        """Distance from disc center to prong center."""
        return self.major_radius_mm - self.prong_diameter_mm

    def prong_radius(self) -> float:
        """Radius of each prong."""
        return self.prong_diameter_mm / 2

    def engaged_prong_position(self, pedal_rotation: float) -> Tuple[float, float]:
        """Get position of the engaging (+X) prong when disc is rotated.

        Args:
            pedal_rotation: Pedal rotation in degrees (0=flat, 45=natural, 90=sharp)

        Returns:
            (x, z) position of the engaging prong center
        """
        # The prong starts perpendicular to string, then rotates toward it
        # String direction is (string_ux, string_uz)
        # Perpendicular (90° CCW) is (-string_uz, string_ux)

        # Base angle: perpendicular to string in flat position
        base_angle = math.atan2(-self.string_ux, self.string_uz)

        # Add pedal rotation (clockwise, so subtract in standard math convention)
        total_angle = base_angle - math.radians(pedal_rotation)

        inset = self.prong_inset_distance()
        prong_x = self.x_mm + inset * math.cos(total_angle)
        prong_z = self.z_mm + inset * math.sin(total_angle)

        return prong_x, prong_z

    def prongs(self, rotation_override: float = None) -> List[Prong]:
        """Return the two prongs positions based on disc rotation.

        Args:
            rotation_override: If provided, use this rotation angle instead of self.rotation_degrees
        """
        rotation = rotation_override if rotation_override is not None else self.rotation_degrees

        # Get engaging prong position
        plus_prong_x, plus_prong_z = self.engaged_prong_position(rotation)

        # Opposite prong is 180° away
        inset = self.prong_inset_distance()
        minus_prong_x = 2 * self.x_mm - plus_prong_x
        minus_prong_z = 2 * self.z_mm - plus_prong_z

        return [
            Prong(plus_prong_x, plus_prong_z, self.prong_diameter_mm, self.material),
            Prong(minus_prong_x, minus_prong_z, self.prong_diameter_mm, self.material),
        ]


@dataclass
class FlatPin:
    """A flat pin where the string bends before going to the tuning pin."""
    x_mm: float
    z_mm: float
    shaft_diameter_mm: float
    thread: str  # e.g., "M6", "M8"
    material: str = "Stainless Steel 316"

    def radius_mm(self) -> float:
        return self.shaft_diameter_mm / 2


@dataclass
class TuningPin:
    """A tuning pin in the neck for adjusting string tension."""
    x_mm: float
    z_mm: float
    y_mm: float  # +10 for +Y plate, -10 for -Y plate
    diameter_mm: float
    thread: str  # e.g., "M4", "M5", "M6"
    plate: str = "+Y"  # "+Y" or "-Y" neck plate (same as disc plate for this string)
    material: str = "Stainless Steel 316"

    def radius_mm(self) -> float:
        return self.diameter_mm / 2


@dataclass
class String:
    """A harp string with all its associated hardware."""
    number: int
    note: str
    frequency_hz: float

    # Geometry
    x_soundboard_mm: float
    z_soundboard_mm: float
    x_flat_mm: float
    z_flat_mm: float
    length_mm: float

    # Physical properties
    outer_diameter_mm: float
    core_diameter_mm: float
    core_material: str
    wrap_material: Optional[str]
    tension_n: float

    # Associated hardware
    flat_pin: FlatPin
    tuning_pin: TuningPin
    natural_disc: Disc
    sharp_disc: Disc

    # Optional double accidental discs
    double_flat_disc: Optional[Disc] = None
    double_sharp_disc: Optional[Disc] = None

    def color(self) -> str:
        """Return SVG color based on note name."""
        note_letter = self.note[0].lower()
        colors = {
            'c': '#FF0000',  # red
            'd': '#999999',  # gray
            'e': '#999999',  # gray
            'f': '#0000FF',  # blue
            'g': '#999999',  # gray
            'a': '#999999',  # gray
            'b': '#999999',  # gray
        }
        return colors.get(note_letter, '#999999')

    def stroke_width(self, scale: float = 0.5) -> float:
        """Return stroke width matching actual string diameter at given scale."""
        return self.outer_diameter_mm * scale


@dataclass
class Neck:
    """The neck of the harp where tuning pins are located."""
    c1_x_mm: float
    c1_z_mm: float
    g7_x_mm: float
    g7_z_mm: float
    clearance_mm: float = 30.0
    thickness_mm: float = 50.0  # Neck thickness - tuning pins sit on top

    def length_mm(self) -> float:
        dx = self.g7_x_mm - self.c1_x_mm
        dz = self.g7_z_mm - self.c1_z_mm
        return math.sqrt(dx*dx + dz*dz)

    def angle_deg(self) -> float:
        dx = self.g7_x_mm - self.c1_x_mm
        dz = self.g7_z_mm - self.c1_z_mm
        return math.degrees(math.atan2(dz, dx))

    def top_z_at_x(self, x_mm: float) -> float:
        """Get the Z coordinate of the top of the neck at a given X position."""
        # Interpolate along neck line, then add thickness
        if abs(self.g7_x_mm - self.c1_x_mm) < 0.001:
            base_z = self.c1_z_mm
        else:
            t = (x_mm - self.c1_x_mm) / (self.g7_x_mm - self.c1_x_mm)
            base_z = self.c1_z_mm + t * (self.g7_z_mm - self.c1_z_mm)
        return base_z + self.thickness_mm


@dataclass
class Soundboard:
    """The soundboard with spline control points."""
    control_points: List[Dict[str, float]]
    curve_type: str = "8-point cubic spline"


@dataclass
class Harp:
    """The complete harp assembly."""
    name: str
    string_count: int
    strings: List[String]
    neck: Neck
    soundboard: Soundboard
    metadata: Dict

    # Tuning pin distribution mode
    tuning_pin_mode: str = "angle_45"  # "angle_45" or "even"


# =============================================================================
# DATA LOADING
# =============================================================================

def recalculate_string_positions(strings_data: List[Dict], finger_gap_mm: float, string_tilt_deg: float) -> List[Dict]:
    """
    Recalculate string positions to enforce constraints:

    CRITICAL CONSTRAINTS (in order of priority):
    1. STRING LENGTH - NEVER CHANGES (determines pitch)
    2. Finger gap - 14mm perpendicular distance between adjacent strings
    3. String tilt - all strings parallel at θ degrees from vertical
    4. Soundboard Z - fixed from soundboard curve

    Given fixed string length L and tilt angle θ:
        dz = L * cos(θ)  → z_top = z_bottom + L * cos(θ)
        dx = L * sin(θ)  → x_top = x_bottom + L * sin(θ)

    The finger gap constraint determines x_bottom for each string.
    """
    if not strings_data:
        return strings_data

    tilt_rad = math.radians(string_tilt_deg)
    cos_tilt = math.cos(tilt_rad)
    sin_tilt = math.sin(tilt_rad)

    # Keep first string position as reference
    first_string = strings_data[0]
    base_x_bottom = first_string['x_bottom_mm']
    base_z_bottom = first_string['soundboard_z_mm']

    for i, s in enumerate(strings_data):
        # STRING LENGTH IS SACRED - never change it
        string_length = s['length_mm']

        # Z offset from first string at soundboard (due to curved soundboard)
        dz_from_first = s['soundboard_z_mm'] - base_z_bottom

        # Calculate X position at soundboard that gives exactly i*finger_gap perpendicular distance
        # perp_dist = dx*cos(θ) - dz*sin(θ) = i*finger_gap
        # dx = (i*finger_gap + dz*sin(θ)) / cos(θ)
        if abs(cos_tilt) > 0.001:
            dx_from_first = (i * finger_gap_mm + dz_from_first * sin_tilt) / cos_tilt
        else:
            dx_from_first = i * finger_gap_mm

        s['x_bottom_mm'] = base_x_bottom + dx_from_first

        # Calculate top position from FIXED string length and tilt angle
        # For a string of length L at tilt θ from vertical:
        #   dz = L * cos(θ)
        #   dx = L * sin(θ)
        dz_string = string_length * cos_tilt
        dx_string = string_length * sin_tilt

        s['z_top_mm'] = s['soundboard_z_mm'] + dz_string
        s['x_top_mm'] = s['x_bottom_mm'] + dx_string

        # Also update x_mm (general position reference)
        s['x_mm'] = s['x_bottom_mm']

    return strings_data


def tangent_point_on_circle(px: float, pz: float, cx: float, cz: float, r: float) -> Tuple[float, float]:
    """Calculate tangent point on circle from external point, for LEFT side entry (CW wrap)."""
    dx = cx - px
    dz = cz - pz
    dist = math.sqrt(dx*dx + dz*dz)
    if dist <= r:
        return cx, cz
    theta = math.atan2(dz, dx)
    alpha = math.asin(r / dist)
    tp_angle = theta + math.pi/2 - alpha
    return cx + r * math.cos(tp_angle), cz + r * math.sin(tp_angle)


def calculate_disc_position_from_physics(
    disc_data: Dict,
    string_length_mm: float,
    x_soundboard: float, z_soundboard: float,
    x_tangent: float, z_tangent: float,
    semitones_from_flat: int,
    rotation_degrees: float
) -> Dict:
    """
    Calculate disc position from physics.

    The disc CENTER is positioned ON the actual string line (soundboard to
    tangent point on flat pin) at the distance corresponding to the semitone
    pitch change.

    Args:
        disc_data: Dict with disc properties (major_radius_mm, prong_diameter_mm, etc.)
        string_length_mm: Total string length (vibrating length)
        x_soundboard, z_soundboard: String anchor at soundboard
        x_tangent, z_tangent: Tangent point on flat pin (where string actually goes)
        semitones_from_flat: 1 for natural, 2 for sharp
        rotation_degrees: Disc rotation angle when engaged (45° for natural, 90° for sharp)
    """
    SEMITONE_RATIO = 2 ** (1/12)  # ≈ 1.0595

    # Calculate where disc center should be on the string
    # Each semitone shortens effective length by factor of 1/SEMITONE_RATIO
    engagement_distance = string_length_mm / (SEMITONE_RATIO ** semitones_from_flat)

    # String direction unit vector (from soundboard to tangent point - the actual string path)
    dx = x_tangent - x_soundboard
    dz = z_tangent - z_soundboard
    string_geom_len = math.sqrt(dx*dx + dz*dz)
    ux, uz = dx / string_geom_len, dz / string_geom_len

    # Disc center is ON the string at the engagement distance from soundboard
    disc_x = x_soundboard + engagement_distance * ux
    disc_z = z_soundboard + engagement_distance * uz

    # Store string direction for rendering (unit vector)
    disc_data['string_ux'] = ux
    disc_data['string_uz'] = uz

    # Update disc data
    disc_data['x_mm'] = disc_x
    disc_data['z_mm'] = disc_z
    disc_data['rotation_degrees'] = rotation_degrees

    return disc_data


# Disc rotation angles for each position
NATURAL_ROTATION_DEG = 45.0  # 3 o'clock to ~4:30
SHARP_ROTATION_DEG = 90.0    # 3 o'clock to 6 o'clock
STRING_DEFLECTION_MM = 1.5   # How far prong pushes string sideways


def load_harp_from_json(json_path: str, tuning_pin_mode: str = "angle_45", neck_x_offset: float = 65.0, target_angle: float = 50.0, neck_thickness: float = 50.0, finger_gap_mm: float = 14.0, g7_drop: float = 0.0, g7_x_shift: float = 0.0, c1_angle: float = 60.0, g7_angle: float = 45.0) -> Harp:
    """Load harp configuration from JSON file."""

    with open(json_path, 'r') as f:
        data = json.load(f)

    meta = data['metadata']
    strings_data = data['strings']

    # Get string tilt angle from metadata
    string_tilt_deg = meta.get('string_tilt_degrees', -3.0)

    # Recalculate string positions to enforce finger gap constraint
    strings_data = recalculate_string_positions(strings_data, finger_gap_mm, string_tilt_deg)

    # Calculate disc positions from physics (semitone ratio and rotation geometry)
    for s in strings_data:
        x_soundboard = s['x_bottom_mm']
        z_soundboard = s['soundboard_z_mm']
        x_flat = s['x_top_mm']
        z_flat = s['z_top_mm']
        string_length = s['length_mm']

        # Calculate flat pin radius based on string range
        if s['number'] <= 14:  # Bass
            flat_r = 3.5  # M8 pin
        elif s['number'] <= 28:  # Mid
            flat_r = 2.5  # M6 pin
        else:  # Treble
            flat_r = 2.0  # M5 pin

        # Calculate tangent point on flat pin (where string actually touches)
        x_tangent, z_tangent = tangent_point_on_circle(
            x_soundboard, z_soundboard,
            x_flat, z_flat, flat_r
        )

        # Natural disc: 1 semitone up from flat, 45° rotation
        if 'natural_disc' in s:
            s['natural_disc'] = calculate_disc_position_from_physics(
                s['natural_disc'],
                string_length,
                x_soundboard, z_soundboard,
                x_tangent, z_tangent,
                semitones_from_flat=1,
                rotation_degrees=NATURAL_ROTATION_DEG
            )

        # Sharp disc: 2 semitones up from flat, 90° rotation
        if 'sharp_disc' in s:
            s['sharp_disc'] = calculate_disc_position_from_physics(
                s['sharp_disc'],
                string_length,
                x_soundboard, z_soundboard,
                x_tangent, z_tangent,
                semitones_from_flat=2,
                rotation_degrees=SHARP_ROTATION_DEG
            )

    # Build neck plates - parallel to original neckref line, 30mm above highest flat pin
    # Both +Y and -Y plates are identical (same slope) except mirrored in Y
    neck_clearance_mm = 30.0  # Fixed 30mm clearance above highest flat pin

    # Find the highest flat pin and its position
    highest_flat_z = max(s['z_top_mm'] for s in strings_data)
    highest_flat_string = max(strings_data, key=lambda s: s['z_top_mm'])
    highest_flat_x = highest_flat_string['x_top_mm']

    # Get original neckref slope from JSON
    neck_data = meta['neck_line']
    orig_dx = neck_data['g7_x_mm'] - neck_data['c1_x_mm']
    orig_dz = neck_data['g7_z_mm'] - neck_data['c1_z_mm']
    neck_slope = orig_dz / orig_dx if abs(orig_dx) > 0.001 else 0

    # Neck X range with offset
    # g7_x_shift: negative brings g7 tuning pin closer to its flat pin (-X direction)
    neck_c1_x = neck_data['c1_x_mm'] + neck_x_offset
    neck_g7_x = neck_data['g7_x_mm'] + neck_x_offset + g7_x_shift

    # Position plates parallel to neckref, 30mm above highest flat pin
    # At the X position of the highest flat pin, neck Z = highest_flat_z + 30
    neck_z_at_highest = highest_flat_z + neck_clearance_mm

    # Calculate c1 and g7 Z using the slope from the reference point
    # Apply g7_drop to lower the g7 end (increases downward slope)
    neck_c1_z = neck_z_at_highest + neck_slope * (neck_c1_x - highest_flat_x)
    neck_g7_z = neck_z_at_highest + neck_slope * (neck_g7_x - highest_flat_x) - g7_drop

    neck = Neck(
        c1_x_mm=neck_c1_x,
        c1_z_mm=neck_c1_z,
        g7_x_mm=neck_g7_x,
        g7_z_mm=neck_g7_z,
        clearance_mm=neck_clearance_mm,
        thickness_mm=neck_thickness
    )

    # Build soundboard
    soundboard = Soundboard(
        control_points=meta['soundboard_spline_control_points'],
        curve_type=meta.get('soundboard_curve', '8-point cubic spline')
    )

    # Calculate tuning pin positions on neckref line
    # Tuning pins rest on their neck plates along the neckref line
    def neckref_z_at_x(x):
        """Get Z position on neckref line at given X."""
        if abs(neck.g7_x_mm - neck.c1_x_mm) < 0.001:
            return neck.c1_z_mm
        t = (x - neck.c1_x_mm) / (neck.g7_x_mm - neck.c1_x_mm)
        return neck.c1_z_mm + t * (neck.g7_z_mm - neck.c1_z_mm)

    n_strings = len(strings_data)
    if tuning_pin_mode == "even":
        # Evenly distributed along neckref line
        tuning_positions = []
        for i in range(n_strings):
            t = i / (n_strings - 1)
            tp_x = neck.c1_x_mm + t * (neck.g7_x_mm - neck.c1_x_mm)
            tp_z = neckref_z_at_x(tp_x)
            tuning_positions.append((tp_x, tp_z))
    elif tuning_pin_mode == "constant_angle":
        # Each tuning pin positioned to achieve target angle from flat pin to tuning pin CENTER
        # Account for tuning pin resting on neckref (center is at neckref_z + radius)
        target_angle_rad = math.radians(target_angle)
        tuning_positions = []
        for s in strings_data:
            flat_x = s['x_top_mm']
            flat_z = s['z_top_mm']

            # Determine tuning pin radius based on string number
            if s['number'] <= 14:
                tp_radius = 3.0  # M6
            elif s['number'] <= 28:
                tp_radius = 2.5  # M5
            else:
                tp_radius = 2.0  # M4

            # Iteratively find tp_x that gives target angle
            # angle = atan2(dx, dz) where dz = (neckref_z(tp_x) + tp_radius) - flat_z
            # Start with initial guess
            tp_x = flat_x + 50  # Initial guess
            for _ in range(10):  # Iterate to converge
                neckref_z = neckref_z_at_x(tp_x)
                tp_center_z = neckref_z + tp_radius  # Pin center, not neckref
                dz = tp_center_z - flat_z
                # For target angle: tan(angle) = dx/dz, so dx = dz * tan(angle)
                dx = dz * math.tan(target_angle_rad)
                tp_x = flat_x + dx

            tp_z = neckref_z_at_x(tp_x)  # Store neckref z, radius added later
            tuning_positions.append((tp_x, tp_z))
    elif tuning_pin_mode == "linear_angle":
        # Angles interpolated linearly from c1_angle to g7_angle
        tuning_positions = []
        n = len(strings_data)
        for i, s in enumerate(strings_data):
            # Interpolate target angle
            t = i / (n - 1) if n > 1 else 0
            target_angle_deg = c1_angle + t * (g7_angle - c1_angle)
            target_angle_rad = math.radians(target_angle_deg)

            flat_x = s['x_top_mm']
            flat_z = s['z_top_mm']

            # Determine tuning pin radius based on string number
            if s['number'] <= 14:
                tp_radius = 3.0  # M6
            elif s['number'] <= 28:
                tp_radius = 2.5  # M5
            else:
                tp_radius = 2.0  # M4

            # Iteratively find tp_x that gives target angle
            tp_x = flat_x + 50  # Initial guess
            for _ in range(10):
                neckref_z = neckref_z_at_x(tp_x)
                tp_center_z = neckref_z + tp_radius
                dz = tp_center_z - flat_z
                dx = dz * math.tan(target_angle_rad)
                tp_x = flat_x + dx

            tp_z = neckref_z_at_x(tp_x)
            tuning_positions.append((tp_x, tp_z))
    elif tuning_pin_mode == "endpoint_angle":
        # c1 and g7 positioned by their angles FROM THE NECKREF LINE
        # g7 side lowered by g7_drop percentage, others evenly distributed
        n = len(strings_data)

        # Neckref line angle from horizontal
        neckref_dx = neck.g7_x_mm - neck.c1_x_mm
        neckref_dz = neck.g7_z_mm - neck.c1_z_mm
        neckref_angle = math.atan2(neckref_dz, neckref_dx)  # Angle from horizontal

        # Calculate c1 tuning pin position
        # c1_angle is measured from neckref line toward vertical
        c1_flat_x = strings_data[0]['x_top_mm']
        c1_flat_z = strings_data[0]['z_top_mm']
        c1_tp_radius = 3.0  # Bass string M6

        # String angle from horizontal = neckref_angle + c1_angle
        c1_string_angle = neckref_angle + math.radians(c1_angle)
        c1_tp_x = c1_flat_x + 50  # Initial guess
        for _ in range(10):
            neckref_z = neckref_z_at_x(c1_tp_x)
            c1_tp_center_z = neckref_z + c1_tp_radius
            # Distance from flat pin to tuning pin center
            dist = math.sqrt((c1_tp_x - c1_flat_x)**2 + (c1_tp_center_z - c1_flat_z)**2)
            # Position based on angle from horizontal
            c1_tp_x = c1_flat_x + dist * math.cos(c1_string_angle)
            # Recalculate Z on neckref
            neckref_z = neckref_z_at_x(c1_tp_x)
            c1_tp_center_z = neckref_z + c1_tp_radius
        c1_tp_z = neckref_z_at_x(c1_tp_x)

        # Calculate g7 tuning pin position
        g7_flat_x = strings_data[-1]['x_top_mm']
        g7_flat_z = strings_data[-1]['z_top_mm']
        g7_tp_radius = 2.0  # Treble string M4

        g7_string_angle = neckref_angle + math.radians(g7_angle)
        g7_tp_x = g7_flat_x + 50  # Initial guess
        for _ in range(10):
            neckref_z = neckref_z_at_x(g7_tp_x)
            g7_tp_center_z = neckref_z + g7_tp_radius
            dist = math.sqrt((g7_tp_x - g7_flat_x)**2 + (g7_tp_center_z - g7_flat_z)**2)
            g7_tp_x = g7_flat_x + dist * math.cos(g7_string_angle)
            neckref_z = neckref_z_at_x(g7_tp_x)
            g7_tp_center_z = neckref_z + g7_tp_radius
        g7_tp_z = neckref_z_at_x(g7_tp_x)

        # Apply g7_drop as percentage (lower g7 side)
        # g7_drop is a percentage (e.g., 30 means 30%)
        z_diff = c1_tp_z - g7_tp_z  # How much lower g7 is than c1
        g7_tp_z = c1_tp_z - z_diff * (1 + g7_drop / 100)

        # Distribute other tuning pins evenly between c1 and g7
        tuning_positions = []
        for i in range(n):
            t = i / (n - 1) if n > 1 else 0
            tp_x = c1_tp_x + t * (g7_tp_x - c1_tp_x)
            tp_z = c1_tp_z + t * (g7_tp_z - c1_tp_z)
            tuning_positions.append((tp_x, tp_z))
    else:
        # Default: tuning pin at +X offset from flat pin, on neckref line
        tuning_positions = []
        for s in strings_data:
            flat_x = s['x_top_mm']
            tp_x = flat_x + neck_x_offset
            tp_z = neckref_z_at_x(tp_x)
            tuning_positions.append((tp_x, tp_z))

    # Build strings with all components
    strings = []
    for i, s in enumerate(strings_data):
        # Determine pin sizes based on string range
        if s['number'] <= 14:  # Bass (C1-B2)
            flat_thread = "M8"
            flat_shaft = 7.0
            tuning_thread = "M6"
            tuning_diameter = 6.0
        elif s['number'] <= 28:  # Mid (C3-B4)
            flat_thread = "M6"
            flat_shaft = 5.0
            tuning_thread = "M5"
            tuning_diameter = 5.0
        else:  # Treble (C5-G7)
            flat_thread = "M5"
            flat_shaft = 4.0
            tuning_thread = "M4"
            tuning_diameter = 4.0

        # Prong diameter scales with disc
        nd = s['natural_disc']
        sd = s['sharp_disc']
        prong_d = nd.get('prong_diameter_mm', 5.0)

        # Build discs
        # Plate assignment: odd strings on +Y, even strings on -Y
        expected_plate = "+Y" if s['number'] % 2 == 1 else "-Y"

        natural_disc = Disc(
            x_mm=nd['x_mm'],
            z_mm=nd['z_mm'],
            y_mm=nd.get('y_mm', 10 if expected_plate == "+Y" else -10),
            major_radius_mm=nd['major_radius_mm'],
            minor_radius_mm=nd['minor_radius_mm'],
            disc_type="natural",
            prong_diameter_mm=prong_d,
            plate=nd.get('plate', expected_plate),
            rotation_degrees=nd.get('rotation_degrees', NATURAL_ROTATION_DEG),
            string_ux=nd.get('string_ux', 0),
            string_uz=nd.get('string_uz', 1),
        )

        sharp_disc = Disc(
            x_mm=sd['x_mm'],
            z_mm=sd['z_mm'],
            y_mm=sd.get('y_mm', 10 if expected_plate == "+Y" else -10),
            major_radius_mm=sd['major_radius_mm'],
            minor_radius_mm=sd['minor_radius_mm'],
            disc_type="sharp",
            prong_diameter_mm=prong_d,
            plate=sd.get('plate', expected_plate),
            rotation_degrees=sd.get('rotation_degrees', SHARP_ROTATION_DEG),
            string_ux=sd.get('string_ux', 0),
            string_uz=sd.get('string_uz', 1),
        )

        # Double-flat and double-sharp discs are not used in standard pedal harp
        double_flat_disc = None
        double_sharp_disc = None

        # Build pins
        flat_pin = FlatPin(
            x_mm=s['x_top_mm'],
            z_mm=s['z_top_mm'],
            shaft_diameter_mm=flat_shaft,
            thread=flat_thread,
        )

        tp_x, tp_z = tuning_positions[i]
        # Tuning pin on same plate as disc (odd strings +Y, even strings -Y)
        # Offset Z by pin radius so the pin RESTS on the neckref line (bottom touches line)
        tuning_pin = TuningPin(
            x_mm=tp_x,
            z_mm=tp_z + tuning_diameter / 2,
            y_mm=10 if expected_plate == "+Y" else -10,
            diameter_mm=tuning_diameter,
            thread=tuning_thread,
            plate=expected_plate,
        )

        # Build string
        string = String(
            number=s['number'],
            note=s['note'],
            frequency_hz=s['frequency_hz'],
            x_soundboard_mm=s['x_bottom_mm'],
            z_soundboard_mm=s['soundboard_z_mm'],
            x_flat_mm=s['x_top_mm'],
            z_flat_mm=s['z_top_mm'],
            length_mm=s['length_mm'],
            outer_diameter_mm=s['outer_diameter_mm'],
            core_diameter_mm=s['core_diameter_mm'],
            core_material=s['core_material'],
            wrap_material=s.get('wrap_material'),
            tension_n=s['tension_n'],
            flat_pin=flat_pin,
            tuning_pin=tuning_pin,
            natural_disc=natural_disc,
            sharp_disc=sharp_disc,
            double_flat_disc=double_flat_disc,
            double_sharp_disc=double_sharp_disc,
        )
        strings.append(string)

    return Harp(
        name=meta.get('name', 'CLEMENTS47'),
        string_count=meta['string_count'],
        strings=strings,
        neck=neck,
        soundboard=soundboard,
        metadata=meta,
        tuning_pin_mode=tuning_pin_mode,
    )


# =============================================================================
# CONSTRAINT VALIDATION
# =============================================================================

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


# =============================================================================
# BOM GENERATION
# =============================================================================

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


# =============================================================================
# SVG RENDERING
# =============================================================================

class HarpRenderer:
    """Renders a Harp to SVG using svgwrite."""

    # Disc colors by type
    DISC_COLORS = {
        "natural": ("#B8A8D4", "#55438B"),      # Purple fill, dark purple stroke
        "sharp": ("#D4C4A8", "#8B7355"),        # Tan fill, brown stroke
        "double_flat": ("#E8D4A8", "#8B7355"),  # Light tan fill
        "double_sharp": ("#C4A8E8", "#55438B"), # Light purple fill
    }

    def __init__(self, harp: Harp, scale: float = 0.5, padding: float = 20):
        self.harp = harp
        self.scale = scale
        self.padding = padding

        # Calculate bounds
        self._calculate_bounds()

    def _calculate_bounds(self):
        """Calculate SVG bounds from harp geometry."""
        all_x = []
        all_z = []

        for s in self.harp.strings:
            all_x.extend([s.x_soundboard_mm, s.x_flat_mm, s.tuning_pin.x_mm])
            all_z.extend([s.z_soundboard_mm, s.z_flat_mm, s.tuning_pin.z_mm])
            all_x.extend([s.natural_disc.x_mm, s.sharp_disc.x_mm])
            all_z.extend([s.natural_disc.z_mm, s.sharp_disc.z_mm])

        self.min_x = min(all_x) - 50
        self.max_x = max(all_x) + 50
        self.min_z = min(all_z) - 50
        self.max_z = max(all_z) + 50

        self.width = int((self.max_x - self.min_x) * self.scale + 2 * self.padding)
        self.height = int((self.max_z - self.min_z) * self.scale + 2 * self.padding)

    def _tx(self, x_mm: float) -> float:
        """Transform X coordinate to SVG space."""
        return (x_mm - self.min_x) * self.scale + self.padding

    def _tz(self, z_mm: float) -> float:
        """Transform Z coordinate to SVG space (flip for SVG Y-down)."""
        return self.height - ((z_mm - self.min_z) * self.scale + self.padding)

    def _scale(self, val_mm: float) -> float:
        """Scale a dimension to SVG space."""
        return val_mm * self.scale

    def render(self, output_path: str, pedal_position: str = "flat"):
        """Render the harp to an SVG file.

        Args:
            output_path: Path to save the SVG file
            pedal_position: "flat", "natural", or "sharp"
                - flat: all discs at 0° rotation, strings straight
                - natural: natural disc at 45°, strings deflect at natural engagement
                - sharp: natural at 45°, sharp at 90°, strings deflect at both
        """

        dwg = svgwrite.Drawing(output_path, size=(self.width, self.height))

        # White background
        dwg.add(dwg.rect((0, 0), (self.width, self.height), fill='white'))

        # Reference lines
        self._draw_reference_lines(dwg)

        # Draw discs at appropriate rotation (commented out for now)
        # self._draw_discs(dwg, pedal_position)

        # Draw strings with deflection at engagement points
        self._draw_strings(dwg, pedal_position)

        dwg.save()
        print(f"SVG saved to {output_path} (pedal position: {pedal_position})")

    def _draw_reference_lines(self, dwg):
        """Draw reference lines."""

        # Green dashed line: flat pins c1 to g7
        c1 = self.harp.strings[0]
        g7 = self.harp.strings[-1]

        dwg.add(dwg.line(
            (self._tx(c1.x_flat_mm), self._tz(c1.z_flat_mm)),
            (self._tx(g7.x_flat_mm), self._tz(g7.z_flat_mm)),
            stroke='#90EE90',
            stroke_width=0.5,
            stroke_dasharray='4,2'
        ))

        # Brown solid line: neckref from c1 tuning pin to g7 tuning pin (bottom of pins)
        # Use actual tuning pin positions, offset down by pin radius so line is at pin bottoms
        c1_tp = c1.tuning_pin
        g7_tp = g7.tuning_pin
        dwg.add(dwg.line(
            (self._tx(c1_tp.x_mm), self._tz(c1_tp.z_mm - c1_tp.radius_mm())),
            (self._tx(g7_tp.x_mm), self._tz(g7_tp.z_mm - g7_tp.radius_mm())),
            stroke='#8B4513',
            stroke_width=0.8
        ))

    def _draw_disc(self, dwg, disc: Disc, rotation_override: float = None):
        """Draw a single disc with prongs.

        Args:
            disc: The disc to draw
            rotation_override: If provided, draw disc at this pedal rotation angle.
                             Use 0 for flat position, None for engaged position.
        """
        # Pedal rotation (0 = flat, 45 = natural, 90 = sharp)
        pedal_rotation = rotation_override if rotation_override is not None else disc.rotation_degrees

        # String direction in world coords is stored in disc (ux, uz)
        # Convert to SVG direction: SVG_x = world_x, SVG_y = -world_z (flipped)
        svg_dx = disc.string_ux * self.scale  # X direction same
        svg_dy = -disc.string_uz * self.scale  # Z flipped to Y

        # String angle in SVG from horizontal (positive X axis)
        string_angle_svg = math.degrees(math.atan2(svg_dy, svg_dx))

        # For prongs to be perpendicular to string:
        # String is at string_angle_svg. Perpendicular is string_angle_svg + 90°.
        # Prong axis starts at 0° (horizontal). To get perpendicular to string, rotate by (string_angle_svg + 90°).
        base_rotation = string_angle_svg + 90.0

        # Total rotation = base (perpendicular to string) + pedal rotation
        total_rotation = base_rotation + pedal_rotation

        fill, stroke = self.DISC_COLORS.get(disc.disc_type, ('#888', '#444'))

        # Disc center in SVG coordinates
        cx_svg = self._tx(disc.x_mm)
        cy_svg = self._tz(disc.z_mm)

        # Create a group for the disc and prongs, rotated about disc center
        disc_group = dwg.g(transform=f"rotate({total_rotation:.1f}, {cx_svg:.1f}, {cy_svg:.1f})")

        # Disc (ellipse) - thin stroke to show true boundary
        disc_group.add(dwg.ellipse(
            center=(cx_svg, cy_svg),
            r=(self._scale(disc.major_radius_mm), self._scale(disc.minor_radius_mm)),
            fill=fill,
            stroke=stroke,
            stroke_width=0.15
        ))

        # Prongs at 3 o'clock and 9 o'clock (in disc's local coordinate system)
        # These are drawn at 0° rotation, then the group rotation handles the rest
        prong_r = self._scale(disc.prong_diameter_mm / 2)
        inset_distance = disc.major_radius_mm - disc.prong_diameter_mm

        # +X prong (3 o'clock in local coords, perpendicular to string after base rotation)
        plus_x_prong_cx = cx_svg + self._scale(inset_distance)
        plus_x_prong_cy = cy_svg
        disc_group.add(dwg.circle(
            center=(plus_x_prong_cx, plus_x_prong_cy),
            r=prong_r,
            fill='#888',
            stroke='#444',
            stroke_width=0.1
        ))

        # -X prong (9 o'clock in local coords)
        minus_x_prong_cx = cx_svg - self._scale(inset_distance)
        minus_x_prong_cy = cy_svg
        disc_group.add(dwg.circle(
            center=(minus_x_prong_cx, minus_x_prong_cy),
            r=prong_r,
            fill='#888',
            stroke='#444',
            stroke_width=0.1
        ))

        dwg.add(disc_group)

    def _draw_discs(self, dwg, pedal_position: str = "flat"):
        """Draw all discs for all strings.

        Args:
            pedal_position: "flat" (0°), "natural" (45°), or "sharp" (90°)
        """
        for s in self.harp.strings:
            if pedal_position == "flat":
                # Both discs at 0° pedal rotation (prongs perpendicular to string)
                self._draw_disc(dwg, s.natural_disc, rotation_override=0)
                self._draw_disc(dwg, s.sharp_disc, rotation_override=0)
            elif pedal_position == "natural":
                # Natural disc engaged (45°), sharp disc at 0°
                self._draw_disc(dwg, s.natural_disc, rotation_override=NATURAL_ROTATION_DEG)
                self._draw_disc(dwg, s.sharp_disc, rotation_override=0)
            elif pedal_position == "sharp":
                # Both discs engaged: natural at 45°, sharp at 90°
                self._draw_disc(dwg, s.natural_disc, rotation_override=NATURAL_ROTATION_DEG)
                self._draw_disc(dwg, s.sharp_disc, rotation_override=SHARP_ROTATION_DEG)

    def _tangent_point_left(self, px: float, pz: float, cx: float, cz: float, r: float) -> Tuple[float, float]:
        """Calculate LEFT side tangent point on circle centered at (cx,cz) from external point (px,pz)."""
        dx = cx - px
        dz = cz - pz
        dist = math.sqrt(dx*dx + dz*dz)
        if dist <= r:
            return cx, cz
        theta = math.atan2(dz, dx)
        alpha = math.asin(r / dist)
        tp_angle = theta + math.pi/2 - alpha
        return cx + r * math.cos(tp_angle), cz + r * math.sin(tp_angle)

    def _tangent_point_right(self, px: float, pz: float, cx: float, cz: float, r: float) -> Tuple[float, float]:
        """Calculate RIGHT side tangent point on circle centered at (cx,cz) from external point (px,pz)."""
        dx = cx - px
        dz = cz - pz
        dist = math.sqrt(dx*dx + dz*dz)
        if dist <= r:
            return cx, cz
        theta = math.atan2(dz, dx)
        alpha = math.asin(r / dist)
        tp_angle = theta - math.pi/2 + alpha
        return cx + r * math.cos(tp_angle), cz + r * math.sin(tp_angle)

    def _tangent_point(self, px: float, pz: float, cx: float, cz: float, r: float) -> Tuple[float, float]:
        """Calculate tangent point - defaults to LEFT side (used for flat pin)."""
        return self._tangent_point_left(px, pz, cx, cz, r)

    def _external_tangent_left(
        self, c1x: float, c1z: float, r1: float, c2x: float, c2z: float, r2: float
    ) -> Tuple[float, float, float, float]:
        """
        Calculate external tangent points between two circles, LEFT side approach.

        Returns tangent points on c1 and c2 for string going on LEFT side of both circles.

        Args:
            c1x, c1z, r1: First circle center and radius
            c2x, c2z, r2: Second circle center and radius

        Returns:
            (t1x, t1z, t2x, t2z) - tangent points on circle 1 and circle 2
        """
        dx = c2x - c1x
        dz = c2z - c1z
        dist = math.sqrt(dx*dx + dz*dz)

        if dist < abs(r1 - r2):
            # Circles overlap, return centers
            return c1x, c1z, c2x, c2z

        # Angle from c1 to c2
        theta = math.atan2(dz, dx)

        # For external tangent with same-side approach
        if dist > 0:
            beta = math.asin((r2 - r1) / dist) if abs(r2 - r1) < dist else 0
        else:
            beta = 0

        # Perpendicular angle for LEFT side (CCW wrap direction)
        perp_angle = theta + math.pi/2 + beta

        t1x = c1x + r1 * math.cos(perp_angle)
        t1z = c1z + r1 * math.sin(perp_angle)
        t2x = c2x + r2 * math.cos(perp_angle)
        t2z = c2z + r2 * math.sin(perp_angle)

        return t1x, t1z, t2x, t2z

    def _external_tangent_right(
        self, c1x: float, c1z: float, r1: float, c2x: float, c2z: float, r2: float
    ) -> Tuple[float, float, float, float]:
        """
        Calculate external tangent points between two circles, RIGHT side approach.

        Returns tangent points on c1 and c2 for string going on RIGHT side of both circles.

        Args:
            c1x, c1z, r1: First circle center and radius
            c2x, c2z, r2: Second circle center and radius

        Returns:
            (t1x, t1z, t2x, t2z) - tangent points on circle 1 and circle 2
        """
        dx = c2x - c1x
        dz = c2z - c1z
        dist = math.sqrt(dx*dx + dz*dz)

        if dist < abs(r1 - r2):
            # Circles overlap, return centers
            return c1x, c1z, c2x, c2z

        # Angle from c1 to c2
        theta = math.atan2(dz, dx)

        # For external tangent with same-side approach
        if dist > 0:
            beta = math.asin((r2 - r1) / dist) if abs(r2 - r1) < dist else 0
        else:
            beta = 0

        # Perpendicular angle for RIGHT side
        perp_angle = theta - math.pi/2 - beta

        t1x = c1x + r1 * math.cos(perp_angle)
        t1z = c1z + r1 * math.sin(perp_angle)
        t2x = c2x + r2 * math.cos(perp_angle)
        t2z = c2z + r2 * math.sin(perp_angle)

        return t1x, t1z, t2x, t2z

    def _get_engaged_prong_center(self, disc: Disc, pedal_rotation: float) -> Tuple[float, float]:
        """
        Get the world coordinates of the engaging prong center when disc is rotated.

        Args:
            disc: The disc
            pedal_rotation: Pedal rotation in degrees (0=flat, 45=natural, 90=sharp)

        Returns:
            (x, z) coordinates of the engaging prong center
        """
        # Use the disc's method to calculate prong position
        return disc.engaged_prong_position(pedal_rotation)

    def _calculate_string_wrap_both_prongs(
        self,
        from_x: float, from_z: float,
        disc: Disc,
        pedal_rotation: float,
        to_x: float, to_z: float,
        string_r: float
    ) -> Tuple[Tuple[float, float], Tuple[float, float], Tuple[float, float], Tuple[float, float]]:
        """
        Calculate tangent points for string grabbed by BOTH prongs of a disc.

        The disc rotates CW (from neck side view). String path order:
        - From soundboard: first hits 9 o'clock (upper), then 3 o'clock (lower), then to flat pin
        - The 3 o'clock prong is approached from the TOP (from 9 o'clock direction)

        Args:
            from_x, from_z: Starting point (soundboard or previous prong exit)
            disc: The disc with both prongs
            pedal_rotation: Rotation angle in degrees
            to_x, to_z: Ending point (flat pin or next prong entry)
            string_r: String radius

        Returns:
            ((first_prong_data), (second_prong_data)) - tangent points on both prongs in path order
        """
        prong_r = disc.prong_radius()
        effective_r = prong_r + string_r

        # Get both prong positions
        prongs = disc.prongs(rotation_override=pedal_rotation)
        # prongs[0] is the +X prong (3 o'clock in local coords before rotation)
        # prongs[1] is the -X prong (9 o'clock in local coords before rotation)

        lower_prong = prongs[0]  # 3 o'clock - second contact (closer to flat pin)
        upper_prong = prongs[1]  # 9 o'clock - first contact (closer to soundboard)

        # String path: soundboard → upper (9 o'clock) → lower (3 o'clock) → flat pin

        # Upper prong (9 o'clock): string comes from soundboard, wraps around first
        upper_entry_x, upper_entry_z = self._tangent_point_left(
            from_x, from_z, upper_prong.x_mm, upper_prong.z_mm, effective_r)
        upper_exit_x, upper_exit_z = self._tangent_point_right(
            lower_prong.x_mm, lower_prong.z_mm, upper_prong.x_mm, upper_prong.z_mm, effective_r)

        # Lower prong (3 o'clock): approached from TOP (from 9 o'clock), exits toward flat pin
        lower_entry_x, lower_entry_z = self._tangent_point_left(
            upper_prong.x_mm, upper_prong.z_mm, lower_prong.x_mm, lower_prong.z_mm, effective_r)
        lower_exit_x, lower_exit_z = self._tangent_point_right(
            to_x, to_z, lower_prong.x_mm, lower_prong.z_mm, effective_r)

        # Return in path order: first prong (upper/9 o'clock), then second prong (lower/3 o'clock)
        return (
            (upper_entry_x, upper_entry_z, upper_exit_x, upper_exit_z, upper_prong.x_mm, upper_prong.z_mm),
            (lower_entry_x, lower_entry_z, lower_exit_x, lower_exit_z, lower_prong.x_mm, lower_prong.z_mm)
        )

    def _draw_arc_around_prong(self, dwg, prong_x: float, prong_z: float, prong_r: float,
                                entry_x: float, entry_z: float,
                                exit_x: float, exit_z: float,
                                color: str, stroke_width: float, string_r: float):
        """Draw an arc around a prong from entry point to exit point - goes outside the prong."""
        effective_r = prong_r + string_r

        # SVG arc path - sweep-flag=0 for CCW direction (goes outside the prong)
        path_data = f"M {self._tx(entry_x):.1f},{self._tz(entry_z):.1f} "
        path_data += f"A {self._scale(effective_r):.1f},{self._scale(effective_r):.1f} 0 0,0 "
        path_data += f"{self._tx(exit_x):.1f},{self._tz(exit_z):.1f}"
        dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=stroke_width))

    def _draw_strings(self, dwg, pedal_position: str = "flat"):
        """Draw all strings with their pins.

        Args:
            pedal_position: "flat", "natural", or "sharp"
                - flat: strings go straight from soundboard to flat pin
                - natural: strings wrap around engaged natural disc prong
                - sharp: strings wrap around both engaged prongs (sharp then natural)
        """

        for s in self.harp.strings:
            color = s.color()
            sw = s.stroke_width(self.scale)

            flat_r = s.flat_pin.radius_mm()
            tuning_r = s.tuning_pin.radius_mm()
            string_r = s.outer_diameter_mm / 2  # String radius
            nat_prong_r = s.natural_disc.prong_radius()
            sharp_prong_r = s.sharp_disc.prong_radius()

            # Calculate tangent point on flat pin (where string touches in flat position)
            flat_tp_x, flat_tp_z = self._tangent_point(
                s.x_soundboard_mm, s.z_soundboard_mm,
                s.x_flat_mm, s.z_flat_mm, flat_r + string_r
            )

            # Draw string path based on pedal position
            if pedal_position == "flat":
                # String path: 6 o'clock on tuning pin → CCW arc → exit tangent to flat pin
                # Then tangent line to flat pin → CCW arc around flat pin → line to soundboard

                eff_tp_r = tuning_r + string_r  # effective radius for string centerline
                eff_flat_r = flat_r + string_r

                # Calculate external tangent between tuning pin and flat pin (both at effective radii)
                tp_exit_x, tp_exit_z, flat_entry_x, flat_entry_z = self._external_tangent_right(
                    s.tuning_pin.x_mm, s.tuning_pin.z_mm, eff_tp_r,
                    s.x_flat_mm, s.z_flat_mm, eff_flat_r
                )

                # Calculate tuning pin exit angle (for arc segments)
                tp_exit_angle = math.atan2(tp_exit_z - s.tuning_pin.z_mm, tp_exit_x - s.tuning_pin.x_mm)

                # 6 o'clock angle = -π/2 (pointing down in world coords)
                tp_start_angle = -math.pi / 2

                # Tangent from soundboard to flat pin (left side approach)
                flat_exit_x, flat_exit_z = self._tangent_point_left(
                    s.x_soundboard_mm, s.z_soundboard_mm,
                    s.x_flat_mm, s.z_flat_mm, eff_flat_r
                )

                # Calculate flat pin entry/exit angles
                flat_entry_angle = math.atan2(flat_entry_z - s.z_flat_mm, flat_entry_x - s.x_flat_mm)
                flat_exit_angle = math.atan2(flat_exit_z - s.z_flat_mm, flat_exit_x - s.x_flat_mm)

                # Build tuning pin arc as line segments (guarantees following the circle)
                tp_arc_segments = 24
                # CCW in world coords: angles increase
                # Adjust end angle to be greater than start for CCW
                tp_end_angle = tp_exit_angle
                if tp_end_angle < tp_start_angle:
                    tp_end_angle += 2 * math.pi

                tp_arc_path = f"M {self._tx(s.tuning_pin.x_mm + eff_tp_r * math.cos(tp_start_angle)):.2f},{self._tz(s.tuning_pin.z_mm + eff_tp_r * math.sin(tp_start_angle)):.2f}"
                for i in range(1, tp_arc_segments + 1):
                    angle = tp_start_angle + (tp_end_angle - tp_start_angle) * i / tp_arc_segments
                    x = s.tuning_pin.x_mm + eff_tp_r * math.cos(angle)
                    z = s.tuning_pin.z_mm + eff_tp_r * math.sin(angle)
                    tp_arc_path += f" L {self._tx(x):.2f},{self._tz(z):.2f}"

                # Build flat pin arc as line segments
                flat_arc_segments = 16
                # CCW in world coords
                fp_end_angle = flat_exit_angle
                if fp_end_angle < flat_entry_angle:
                    fp_end_angle += 2 * math.pi

                flat_arc_path = ""
                for i in range(1, flat_arc_segments + 1):
                    angle = flat_entry_angle + (fp_end_angle - flat_entry_angle) * i / flat_arc_segments
                    x = s.x_flat_mm + eff_flat_r * math.cos(angle)
                    z = s.z_flat_mm + eff_flat_r * math.sin(angle)
                    flat_arc_path += f" L {self._tx(x):.2f},{self._tz(z):.2f}"

                # Complete path: tuning pin arc → line to flat pin → flat pin arc → line to soundboard
                path_data = (
                    tp_arc_path +
                    f" L {self._tx(flat_entry_x):.2f},{self._tz(flat_entry_z):.2f}" +
                    flat_arc_path +
                    f" L {self._tx(s.x_soundboard_mm):.2f},{self._tz(s.z_soundboard_mm):.2f}"
                )
                dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=sw))

            elif pedal_position == "natural":
                # String wraps around BOTH prongs of natural disc
                # Path: soundboard → 9 o'clock (first) → 3 o'clock (second) → flat pin
                first, second = self._calculate_string_wrap_both_prongs(
                    s.x_soundboard_mm, s.z_soundboard_mm,
                    s.natural_disc, NATURAL_ROTATION_DEG,
                    flat_tp_x, flat_tp_z,
                    string_r
                )
                first_entry_x, first_entry_z, first_exit_x, first_exit_z, first_px, first_pz = first
                second_entry_x, second_entry_z, second_exit_x, second_exit_z, second_px, second_pz = second

                # 1. Soundboard to first prong (9 o'clock) entry
                dwg.add(dwg.line(
                    (self._tx(s.x_soundboard_mm), self._tz(s.z_soundboard_mm)),
                    (self._tx(first_entry_x), self._tz(first_entry_z)),
                    stroke=color,
                    stroke_width=sw
                ))

                # 2. Arc around first prong (9 o'clock)
                self._draw_arc_around_prong(
                    dwg, first_px, first_pz, nat_prong_r,
                    first_entry_x, first_entry_z,
                    first_exit_x, first_exit_z,
                    color, sw, string_r
                )

                # 3. First prong exit to second prong (3 o'clock) entry - approaches from top
                dwg.add(dwg.line(
                    (self._tx(first_exit_x), self._tz(first_exit_z)),
                    (self._tx(second_entry_x), self._tz(second_entry_z)),
                    stroke=color,
                    stroke_width=sw
                ))

                # 4. Arc around second prong (3 o'clock)
                self._draw_arc_around_prong(
                    dwg, second_px, second_pz, nat_prong_r,
                    second_entry_x, second_entry_z,
                    second_exit_x, second_exit_z,
                    color, sw, string_r
                )

                # 5. Second prong exit to flat pin tangent
                dwg.add(dwg.line(
                    (self._tx(second_exit_x), self._tz(second_exit_z)),
                    (self._tx(flat_tp_x), self._tz(flat_tp_z)),
                    stroke=color,
                    stroke_width=sw
                ))

            elif pedal_position == "sharp":
                # String wraps around BOTH prongs of sharp disc, then BOTH prongs of natural disc
                # Path: soundboard → sharp 9 o'clock → sharp 3 o'clock → nat 9 o'clock → nat 3 o'clock → flat pin

                # Get natural disc first prong position (9 o'clock - needed for sharp disc exit target)
                nat_prongs = s.natural_disc.prongs(rotation_override=NATURAL_ROTATION_DEG)
                nat_upper_prong = nat_prongs[1]  # 9 o'clock is the first prong hit from sharp disc

                # Sharp disc: both prongs (first=9 o'clock, second=3 o'clock)
                sharp_first, sharp_second = self._calculate_string_wrap_both_prongs(
                    s.x_soundboard_mm, s.z_soundboard_mm,
                    s.sharp_disc, SHARP_ROTATION_DEG,
                    nat_upper_prong.x_mm, nat_upper_prong.z_mm,  # Exit toward natural disc 9 o'clock
                    string_r
                )
                s1_entry_x, s1_entry_z, s1_exit_x, s1_exit_z, s1_px, s1_pz = sharp_first
                s2_entry_x, s2_entry_z, s2_exit_x, s2_exit_z, s2_px, s2_pz = sharp_second

                # Natural disc: both prongs (coming from sharp disc second prong exit)
                nat_first, nat_second = self._calculate_string_wrap_both_prongs(
                    s2_exit_x, s2_exit_z,  # Coming from sharp disc 3 o'clock exit
                    s.natural_disc, NATURAL_ROTATION_DEG,
                    flat_tp_x, flat_tp_z,
                    string_r
                )
                n1_entry_x, n1_entry_z, n1_exit_x, n1_exit_z, n1_px, n1_pz = nat_first
                n2_entry_x, n2_entry_z, n2_exit_x, n2_exit_z, n2_px, n2_pz = nat_second

                # Draw the string path: 9 segments total

                # 1. Soundboard to sharp first prong (9 o'clock) entry
                dwg.add(dwg.line(
                    (self._tx(s.x_soundboard_mm), self._tz(s.z_soundboard_mm)),
                    (self._tx(s1_entry_x), self._tz(s1_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                # 2. Arc around sharp first prong (9 o'clock)
                self._draw_arc_around_prong(
                    dwg, s1_px, s1_pz, sharp_prong_r,
                    s1_entry_x, s1_entry_z, s1_exit_x, s1_exit_z,
                    color, sw, string_r
                )

                # 3. Sharp first (9 o'clock) to sharp second (3 o'clock)
                dwg.add(dwg.line(
                    (self._tx(s1_exit_x), self._tz(s1_exit_z)),
                    (self._tx(s2_entry_x), self._tz(s2_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                # 4. Arc around sharp second prong (3 o'clock)
                self._draw_arc_around_prong(
                    dwg, s2_px, s2_pz, sharp_prong_r,
                    s2_entry_x, s2_entry_z, s2_exit_x, s2_exit_z,
                    color, sw, string_r
                )

                # 5. Sharp second (3 o'clock) to natural first (9 o'clock)
                dwg.add(dwg.line(
                    (self._tx(s2_exit_x), self._tz(s2_exit_z)),
                    (self._tx(n1_entry_x), self._tz(n1_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                # 6. Arc around natural first prong (9 o'clock)
                self._draw_arc_around_prong(
                    dwg, n1_px, n1_pz, nat_prong_r,
                    n1_entry_x, n1_entry_z, n1_exit_x, n1_exit_z,
                    color, sw, string_r
                )

                # 7. Natural first (9 o'clock) to natural second (3 o'clock)
                dwg.add(dwg.line(
                    (self._tx(n1_exit_x), self._tz(n1_exit_z)),
                    (self._tx(n2_entry_x), self._tz(n2_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                # 8. Arc around natural second prong (3 o'clock)
                self._draw_arc_around_prong(
                    dwg, n2_px, n2_pz, nat_prong_r,
                    n2_entry_x, n2_entry_z, n2_exit_x, n2_exit_z,
                    color, sw, string_r
                )

                # 9. Natural second (3 o'clock) to flat pin
                dwg.add(dwg.line(
                    (self._tx(n2_exit_x), self._tz(n2_exit_z)),
                    (self._tx(flat_tp_x), self._tz(flat_tp_z)),
                    stroke=color, stroke_width=sw
                ))

            # For natural/sharp positions, draw the remaining path to tuning pin
            # (flat position handles the complete continuous path above)
            if pedal_position != "flat":
                # Arc around flat pin (simplified)
                arc_end_x = flat_tp_x + 0.5
                arc_end_z = flat_tp_z + 1.5
                path_data = f"M {self._tx(flat_tp_x):.1f},{self._tz(flat_tp_z):.1f} "
                path_data += f"A {self._scale(flat_r):.1f},{self._scale(flat_r):.1f} 0 0,1 "
                path_data += f"{self._tx(arc_end_x):.1f},{self._tz(arc_end_z):.1f}"
                dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=sw))

                # Line from flat pin to tuning pin
                dwg.add(dwg.line(
                    (self._tx(arc_end_x), self._tz(arc_end_z)),
                    (self._tx(s.tuning_pin.x_mm - tuning_r * 0.7), self._tz(s.tuning_pin.z_mm + tuning_r * 0.7)),
                    stroke=color,
                    stroke_width=sw
                ))

                # Arc around tuning pin (simplified)
                tp_arc_start_x = s.tuning_pin.x_mm - tuning_r * 0.7
                tp_arc_start_z = s.tuning_pin.z_mm + tuning_r * 0.7
                tp_arc_end_x = tp_arc_start_x + 0.3
                tp_arc_end_z = tp_arc_start_z + 1.0
                path_data = f"M {self._tx(tp_arc_start_x):.1f},{self._tz(tp_arc_start_z):.1f} "
                path_data += f"A {self._scale(tuning_r):.1f},{self._scale(tuning_r):.1f} 0 0,1 "
                path_data += f"{self._tx(tp_arc_end_x):.1f},{self._tz(tp_arc_end_z):.1f}"
                dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=sw))

            # Draw flat pin - thin stroke to show true boundary
            dwg.add(dwg.circle(
                center=(self._tx(s.x_flat_mm), self._tz(s.z_flat_mm)),
                r=self._scale(flat_r),
                fill='#888',
                stroke='#444',
                stroke_width=0.15
            ))

            # 6. Draw tuning pin - thin stroke to show true boundary
            dwg.add(dwg.circle(
                center=(self._tx(s.tuning_pin.x_mm), self._tz(s.tuning_pin.z_mm)),
                r=self._scale(tuning_r),
                fill='#654321',
                stroke='#333',
                stroke_width=0.15
            ))


# =============================================================================
# MAIN
# =============================================================================

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

        # Also save harp.svg as the flat position
        renderer.render(args.output, pedal_position="flat")

        print(f"\nGenerated SVGs:")
        print(f"  {base_output}0.svg - Flat position (strings straight)")
        print(f"  {base_output}1.svg - Natural position (natural disc engaged, string deflected)")
        print(f"  {base_output}2.svg - Sharp position (both discs engaged, string deflected twice)")


if __name__ == '__main__':
    main()
