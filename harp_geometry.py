#!/usr/bin/env python3
"""
HARP_GEOMETRY.PY - Geometry calculations and data loading

Contains:
- recalculate_string_positions: Enforce string position constraints
- tangent_point_on_circle: Calculate tangent points
- calculate_disc_position_from_physics: Position discs based on semitone ratio
- load_harp_from_json: Load complete harp from JSON file
"""

import json
import math
from typing import List, Dict, Tuple

from harp_models import (
    Prong, Disc, FlatPin, TuningPin, String, Neck, Soundboard, Harp,
    NATURAL_ROTATION_DEG, SHARP_ROTATION_DEG
)


# Additional constant used in geometry calculations
STRING_DEFLECTION_MM = 1.5   # How far prong pushes string sideways


# =============================================================================
# GEOMETRY FUNCTIONS
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

    The disc center is positioned so that:
    1. Its projection onto the string line is at the engagement distance
    2. It is offset perpendicular to the string so the prong contacts
       the string when rotated to engagement angle

    Args:
        disc_data: Dict with disc properties (major_radius_mm, prong_diameter_mm, etc.)
        string_length_mm: Total string length (vibrating length)
        x_soundboard, z_soundboard: String anchor at soundboard
        x_tangent, z_tangent: Tangent point on flat pin (where string actually goes)
        semitones_from_flat: 1 for natural, 2 for sharp
        rotation_degrees: Disc rotation angle when engaged (45° for natural, 90° for sharp)
    """
    SEMITONE_RATIO = 2 ** (1/12)  # ≈ 1.0595

    # Calculate where prong should contact string (engagement point)
    # Each semitone shortens effective length by factor of 1/SEMITONE_RATIO
    engagement_distance = string_length_mm / (SEMITONE_RATIO ** semitones_from_flat)

    # String direction unit vector (from soundboard to tangent point)
    dx = x_tangent - x_soundboard
    dz = z_tangent - z_soundboard
    string_geom_len = math.sqrt(dx*dx + dz*dz)
    ux, uz = dx / string_geom_len, dz / string_geom_len

    # Disc center is ON the string line at engagement distance from soundboard
    # The prongs rotate in Y (perpendicular to XZ plane), so in XZ view the disc center IS on the string
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


# =============================================================================
# STRING PHYSICS COMPUTATION
# =============================================================================

# Material properties for string tension/diameter calculation
STRING_MATERIALS = {
    "Steel": {"density_kg_m3": 7800, "min_diameter_mm": 0.3, "max_diameter_mm": 1.2},
    "Nylon": {"density_kg_m3": 1150, "min_diameter_mm": 0.5, "max_diameter_mm": 3.0},
    "Bronze": {"density_kg_m3": 8800},  # For wrap wire
}

def compute_frequency(note: str, base_freq: float = 32.7) -> float:
    """Compute frequency from note name using equal temperament.

    Base is c1 = 32.7 Hz. Notes go c, d, e, f, g, a, b with octave number.
    """
    note_offsets = {'c': 0, 'd': 2, 'e': 4, 'f': 5, 'g': 7, 'a': 9, 'b': 11}
    letter = note[0].lower()
    octave = int(note[1])

    # Semitones from c1
    semitones = note_offsets[letter] + (octave - 1) * 12
    return base_freq * (2 ** (semitones / 12))


def compute_string_diameter(length_mm: float, freq_hz: float, tension_n: float,
                           core_material: str, wrap_material: str = None) -> dict:
    """Compute string diameter to achieve target tension.

    From vibrating string equation: T = 4 * L² * f² * μ
    where μ = mass per unit length = ρ * A = ρ * π * (d/2)²

    Solving for d: d = 2 * sqrt(T / (π * ρ * 4 * L² * f²))
    """
    import math

    # Convert to SI units
    L = length_mm / 1000  # meters
    f = freq_hz
    T = tension_n

    if wrap_material and wrap_material != "None":
        # Wound string: core + wrap
        # Effective density is higher due to wrap
        core_density = STRING_MATERIALS[core_material]["density_kg_m3"]
        wrap_density = STRING_MATERIALS.get(wrap_material, {}).get("density_kg_m3", core_density)

        # For wound strings, assume wrap adds ~60% to effective density
        # This is approximate - real calculation depends on wrap pitch
        effective_density = core_density * 0.4 + wrap_density * 0.6

        # Compute outer diameter
        d_m = 2 * math.sqrt(T / (math.pi * effective_density * 4 * L**2 * f**2))
        d_mm = d_m * 1000

        # Core is typically 40-50% of outer diameter for wound strings
        core_d_mm = d_mm * 0.45
        wrap_d_mm = (d_mm - core_d_mm) / 2  # Wrap wire diameter

        return {
            "outer_diameter_mm": round(d_mm, 3),
            "core_diameter_mm": round(core_d_mm, 3),
            "wrap_diameter_mm": round(wrap_d_mm, 3)
        }
    else:
        # Plain string
        density = STRING_MATERIALS[core_material]["density_kg_m3"]
        d_m = 2 * math.sqrt(T / (math.pi * density * 4 * L**2 * f**2))
        d_mm = d_m * 1000

        return {
            "outer_diameter_mm": round(d_mm, 3),
            "core_diameter_mm": round(d_mm, 3),
            "wrap_diameter_mm": None
        }


def compute_target_tension(string_number: int, length_mm: float) -> float:
    """Compute target tension based on string position and length.

    Bass strings (long): higher tension ~220-240N
    Mid strings: ~150-200N
    Treble strings (short): ~50-100N

    Tension roughly proportional to length for consistent feel.
    """
    # Linear interpolation based on typical harp tensions
    # c1 (1514mm) -> 235N, g7 (60mm) -> 50N
    max_len, max_tension = 1515, 235
    min_len, min_tension = 60, 50

    t = (length_mm - min_len) / (max_len - min_len)
    t = max(0, min(1, t))  # Clamp to [0, 1]

    return min_tension + t * (max_tension - min_tension)


def compute_disc_geometry(string_diameter_mm: float) -> dict:
    """Compute disc major/minor radius based on string diameter.

    Disc needs to be large enough for string to wrap around prongs
    without excessive bending stress.
    """
    # Major radius: 2.5-3x string diameter, min 6mm, max 10mm
    major_r = max(6.0, min(10.0, string_diameter_mm * 3.0))

    # Minor radius: about half of major, sized for prong spacing
    minor_r = max(3.0, min(5.5, major_r * 0.5))

    return {
        "major_radius_mm": round(major_r, 1),
        "minor_radius_mm": round(minor_r, 1)
    }


def expand_simplified_json(data: dict) -> dict:
    """Expand simplified JSON to full format with computed values."""

    meta = data['metadata']
    geometry = data.get('geometry', {})
    materials = data.get('materials', {})
    tuning = data.get('tuning', {})

    # Build material lookup by string number
    def get_material(string_num):
        for mat_name, mat_info in materials.items():
            str_range = mat_info.get('strings', [1, 47])
            if str_range[0] <= string_num <= str_range[1]:
                return mat_info['core'], mat_info.get('wrap')
        return "Nylon", None

    base_freq = tuning.get('base_frequency_hz', 32.7)

    # Expand each string
    expanded_strings = []
    for s in data['strings']:
        num = s['number']
        note = s['note']

        # Get positions (handle both old and new field names)
        flat_x = s.get('flat_x_mm', s.get('x_top_mm', 0))
        flat_z = s.get('flat_z_mm', s.get('z_top_mm', 0))
        sb_x = s.get('soundboard_x_mm', s.get('x_bottom_mm', 0))
        sb_z = s.get('soundboard_z_mm', 0)
        length = s.get('length_mm', 1000)

        # Compute frequency
        freq = compute_frequency(note, base_freq)

        # Get material
        core_mat, wrap_mat = get_material(num)

        # Compute tension
        tension = compute_target_tension(num, length)

        # Compute string diameter
        diameters = compute_string_diameter(length, freq, tension, core_mat, wrap_mat)

        # Compute disc geometry
        disc_geom = compute_disc_geometry(diameters['outer_diameter_mm'])

        # Build expanded string entry
        expanded = {
            "number": num,
            "note": note,
            "frequency_hz": round(freq, 1),
            "x_mm": sb_x,
            "x_bottom_mm": sb_x,
            "x_top_mm": flat_x,
            "soundboard_z_mm": sb_z,
            "z_top_mm": flat_z,
            "length_mm": length,
            "core_material": core_mat,
            "wrap_material": wrap_mat,
            "outer_diameter_mm": diameters['outer_diameter_mm'],
            "core_diameter_mm": diameters['core_diameter_mm'],
            "tension_n": round(tension, 1),
            "natural_disc": {
                "x_mm": 0,  # Will be computed by calculate_disc_position_from_physics
                "z_mm": 0,
                "y_mm": 10 if num % 2 == 1 else -10,
                "major_radius_mm": disc_geom['major_radius_mm'],
                "minor_radius_mm": disc_geom['minor_radius_mm'],
                "plate": "+Y" if num % 2 == 1 else "-Y"
            },
            "sharp_disc": {
                "x_mm": 0,
                "z_mm": 0,
                "y_mm": 10 if num % 2 == 1 else -10,
                "major_radius_mm": disc_geom['major_radius_mm'],
                "minor_radius_mm": disc_geom['minor_radius_mm'],
                "plate": "+Y" if num % 2 == 1 else "-Y"
            }
        }
        expanded_strings.append(expanded)

    # Build expanded metadata (compatible with old format)
    expanded_meta = {
        "name": meta.get('name', 'CLEMENTS47'),
        "string_count": meta.get('string_count', 47),
        "string_tilt_degrees": meta.get('string_tilt_degrees', -3),
        "finger_gap_mm": meta.get('finger_gap_mm', 14.0),
        "finger_gap_perpendicular_mm": meta.get('finger_gap_perpendicular_mm', 14.0),
        "constraint": f"{meta.get('finger_gap_mm', 14)}mm finger gap, parallel strings",
        "soundboard_spline_control_points": geometry.get('soundboard_curve', []),
        "soundboard_curve": "8-point cubic spline",
        "neckref": geometry.get('neckref', {}),
        "neck_line": geometry.get('neck_line', geometry.get('neckref', {})),
    }

    return {
        "metadata": expanded_meta,
        "strings": expanded_strings
    }


# =============================================================================
# DATA LOADING
# =============================================================================

def load_harp_from_json(
    json_path: str,
    tuning_pin_mode: str = "angle_45",
    neck_x_offset: float = 65.0,
    target_angle: float = 50.0,
    neck_thickness: float = 50.0,
    finger_gap_mm: float = 14.0,
    g7_drop: float = 0.0,
    g7_x_shift: float = 0.0,
    c1_angle: float = 60.0,
    g7_angle: float = 45.0
) -> Harp:
    """Load harp configuration from JSON file.

    Supports both simplified format (with computed values) and full format.
    Simplified format has 'geometry' and 'materials' sections instead of
    full per-string data.
    """

    with open(json_path, 'r') as f:
        data = json.load(f)

    # Detect simplified format and expand if needed
    if 'geometry' in data or 'materials' in data:
        data = expand_simplified_json(data)

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

        # Calculate flat pin radius based on string range and force analysis
        # Bass strings have highest tension (~237N), need larger pins for bending strength
        if s['number'] <= 9:  # Bass (steel/bronze wound, highest tension)
            flat_r = 4.0  # M10 pin, 8mm shaft
        elif s['number'] <= 28:  # Mid (nylon wound and plain)
            flat_r = 3.0  # M8 pin, 6mm shaft
        else:  # Treble (plain nylon, lowest tension)
            flat_r = 2.5  # M6 pin, 5mm shaft

        # String radius (outer diameter / 2)
        string_r = s['outer_diameter_mm'] / 2

        # Calculate tangent point on flat pin (where string centerline touches)
        # Use flat_r + string_r for effective radius (string wraps around pin)
        x_tangent, z_tangent = tangent_point_on_circle(
            x_soundboard, z_soundboard,
            x_flat, z_flat, flat_r + string_r
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

            # Determine tuning pin radius based on string number (matches force analysis)
            if s['number'] <= 9:  # Bass
                tp_radius = 3.0  # M6
            elif s['number'] <= 28:  # Mid
                tp_radius = 2.5  # M5
            else:  # Treble
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

            # Determine tuning pin radius based on string number (matches force analysis)
            if s['number'] <= 9:  # Bass
                tp_radius = 3.0  # M6
            elif s['number'] <= 28:  # Mid
                tp_radius = 2.5  # M5
            else:  # Treble
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
        # Constraints:
        # 1. c1 tuning pin: 60° from neckref AND 50mm delta Z
        # 2. g7 tuning pin: 45° from neckref AND 25mm delta Z
        # 3. Other tuning pins evenly distributed between c1 and g7
        #
        # NECKREF = line between flat pins (fixed reference)
        # NECK PLATES = line between tuning pins (what we compute)
        n = len(strings_data)

        # Neckref angle (line between FLAT pins, from JSON)
        neckref_c1_x = neck_data['c1_x_mm']
        neckref_c1_z = neck_data['c1_z_mm']
        neckref_g7_x = neck_data['g7_x_mm']
        neckref_g7_z = neck_data['g7_z_mm']
        neckref_angle = math.atan2(neckref_g7_z - neckref_c1_z, neckref_g7_x - neckref_c1_x)

        # c1: 60° from neckref AND 50mm delta Z
        c1_flat_x = strings_data[0]['x_top_mm']
        c1_flat_z = strings_data[0]['z_top_mm']
        c1_tp_radius = 3.0  # Bass string M6
        c1_delta_z = 50.0

        c1_string_angle = neckref_angle + math.radians(c1_angle)
        c1_delta_x = c1_delta_z / math.tan(c1_string_angle)
        c1_tp_x = c1_flat_x + c1_delta_x
        c1_tp_center_z = c1_flat_z + c1_delta_z
        c1_tp_z = c1_tp_center_z - c1_tp_radius  # Neck position (pin rests on it)

        # g7: 45° from neckref AND 25mm delta Z
        g7_flat_x = strings_data[-1]['x_top_mm']
        g7_flat_z = strings_data[-1]['z_top_mm']
        g7_tp_radius = 2.0  # Treble string M4
        g7_delta_z = 25.0

        g7_string_angle = neckref_angle + math.radians(g7_angle)
        g7_delta_x = g7_delta_z / math.tan(g7_string_angle)
        g7_tp_x = g7_flat_x + g7_delta_x
        g7_tp_center_z = g7_flat_z + g7_delta_z
        g7_tp_z = g7_tp_center_z - g7_tp_radius  # Neck position (pin rests on it)

        # Update neck to be the line between tuning pins
        neck.c1_x_mm = c1_tp_x
        neck.c1_z_mm = c1_tp_z
        neck.g7_x_mm = g7_tp_x
        neck.g7_z_mm = g7_tp_z

        # Redefine neckref_z_at_x for the neck plates line
        def neckref_z_at_x(x):
            if abs(neck.g7_x_mm - neck.c1_x_mm) < 0.001:
                return neck.c1_z_mm
            t = (x - neck.c1_x_mm) / (neck.g7_x_mm - neck.c1_x_mm)
            return neck.c1_z_mm + t * (neck.g7_z_mm - neck.c1_z_mm)

        # Distribute tuning pins with EQUAL GAPS (not equal center spacing)
        # First, get all tuning pin diameters (matches force analysis)
        tuning_diameters = []
        for s in strings_data:
            if s['number'] <= 9:  # Bass
                tuning_diameters.append(6.0)  # M6
            elif s['number'] <= 28:  # Mid
                tuning_diameters.append(5.0)  # M5
            else:  # Treble
                tuning_diameters.append(4.0)  # M4

        # Total length from c1 to g7 (center to center)
        total_length = math.sqrt((g7_tp_x - c1_tp_x)**2 + (g7_tp_z - c1_tp_z)**2)

        # Sum of all diameters
        sum_diameters = sum(tuning_diameters)

        # Gap = (total_length - sum_diameters + first_radius + last_radius) / (n-1)
        # Because first and last pins only use half their diameter at the ends
        first_r = tuning_diameters[0] / 2
        last_r = tuning_diameters[-1] / 2
        gap = (total_length - sum_diameters + first_r + last_r) / (n - 1) if n > 1 else 0

        # Direction vector along neck plates
        dx = g7_tp_x - c1_tp_x
        dz = g7_tp_z - c1_tp_z
        ux = dx / total_length
        uz = dz / total_length

        # Position each pin
        tuning_positions = []
        current_pos = 0  # Distance from c1 center along the line
        for i in range(n):
            if i == 0:
                pos = 0
            else:
                prev_r = tuning_diameters[i-1] / 2
                curr_r = tuning_diameters[i] / 2
                current_pos += prev_r + gap + curr_r
                pos = current_pos

            tp_x = c1_tp_x + pos * ux
            tp_z = c1_tp_z + pos * uz
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
        # Determine pin sizes based on string range and force analysis
        # Flat pins sized for bending strength (cantilever load from string wrap)
        # Tuning pins sized for pull-out resistance and tuning torque
        if s['number'] <= 9:  # Bass (C1-D2) - steel/bronze wound, ~220-237N tension
            flat_thread = "M10"
            flat_shaft = 8.0   # SF=3.1 at 237N
            tuning_thread = "M6"
            tuning_diameter = 6.0
            prong_diameter = 5.0
            disc_thickness = 5.0
        elif s['number'] <= 28:  # Mid (E2-B4) - nylon wound and plain, ~100-220N
            flat_thread = "M8"
            flat_shaft = 6.0   # SF=2.1 at 180N
            tuning_thread = "M5"
            tuning_diameter = 5.0
            prong_diameter = 4.0
            disc_thickness = 4.0
        else:  # Treble (C5-G7) - plain nylon, ~50-100N
            flat_thread = "M6"
            flat_shaft = 5.0   # SF=3+ at 100N
            tuning_thread = "M4"
            tuning_diameter = 4.0
            prong_diameter = 3.0
            disc_thickness = 3.0

        # Get disc geometry from JSON (positions calculated earlier)
        nd = s['natural_disc']
        sd = s['sharp_disc']

        # Build discs with computed prong diameter and thickness based on force analysis
        # Plate assignment: odd strings on +Y, even strings on -Y
        expected_plate = "+Y" if s['number'] % 2 == 1 else "-Y"

        natural_disc = Disc(
            x_mm=nd['x_mm'],
            z_mm=nd['z_mm'],
            y_mm=nd.get('y_mm', 10 if expected_plate == "+Y" else -10),
            major_radius_mm=nd['major_radius_mm'],
            minor_radius_mm=nd['minor_radius_mm'],
            disc_type="natural",
            prong_diameter_mm=prong_diameter,  # From force analysis
            plate=nd.get('plate', expected_plate),
            thickness_mm=disc_thickness,  # From force analysis
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
            prong_diameter_mm=prong_diameter,  # From force analysis
            plate=sd.get('plate', expected_plate),
            thickness_mm=disc_thickness,  # From force analysis
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
