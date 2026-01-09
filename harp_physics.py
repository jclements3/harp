#!/usr/bin/env python3
"""
HARP_PHYSICS.PY - String physics and tension calculations

Contains:
- StringMaterial: Physical properties of string materials
- ContactPoint: Points where string contacts hardware
- StringSegment: Straight segments of string
- StringPath: Complete path through harp
- VibrationAnalysis: String vibration and clearance analysis
- Physics analysis functions
"""

import math
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from harp_models import String, Harp


# =============================================================================
# STRING PHYSICS CLASSES
# =============================================================================

@dataclass
class StringMaterial:
    """Physical properties of a string material.

    A string is an entity of potentially infinite length with these properties.
    The actual path and tension are determined by how it's routed through hardware.
    """
    diameter_mm: float
    core_material: str  # "Steel", "Nylon", etc.
    wrap_material: Optional[str] = None  # "Bronze", "Silver", etc. for wound strings

    # Material properties (defaults for steel)
    density_kg_m3: float = 7850.0  # kg/m³ (steel)
    youngs_modulus_gpa: float = 200.0  # GPa (steel)
    tensile_strength_mpa: float = 2000.0  # MPa (music wire)

    def linear_density(self) -> float:
        """Mass per unit length (kg/m) - μ in physics equations."""
        radius_m = self.diameter_mm / 2000.0
        area_m2 = math.pi * radius_m ** 2
        return self.density_kg_m3 * area_m2

    def breaking_tension(self) -> float:
        """Maximum tension before string breaks (Newtons)."""
        radius_m = self.diameter_mm / 2000.0
        area_m2 = math.pi * radius_m ** 2
        return self.tensile_strength_mpa * 1e6 * area_m2

    def tension_for_frequency(self, frequency_hz: float, vibrating_length_mm: float) -> float:
        """Calculate tension needed for target frequency (Newtons).

        From: f = (1/2L) * sqrt(T/μ)
        Solving for T: T = 4 * L² * f² * μ
        """
        L = vibrating_length_mm / 1000.0  # meters
        mu = self.linear_density()
        return 4 * L**2 * frequency_hz**2 * mu

    def frequency_for_tension(self, tension_n: float, vibrating_length_mm: float) -> float:
        """Calculate frequency for given tension and length (Hz).

        f = (1/2L) * sqrt(T/μ)
        """
        L = vibrating_length_mm / 1000.0
        mu = self.linear_density()
        if mu <= 0 or L <= 0:
            return 0
        return (1 / (2 * L)) * math.sqrt(tension_n / mu)

    def safety_factor(self, tension_n: float) -> float:
        """Ratio of breaking tension to actual tension."""
        return self.breaking_tension() / tension_n if tension_n > 0 else float('inf')


@dataclass
class ContactPoint:
    """A point where the string contacts hardware (pin, prong, etc.).

    The string wraps around this point with a certain wrap angle,
    creating forces on the hardware.
    """
    x_mm: float
    z_mm: float
    radius_mm: float  # Radius of the pin/prong
    contact_type: str  # "tuning_pin", "flat_pin", "prong", "soundboard"

    # Calculated after path is built
    wrap_angle_rad: float = 0.0  # How much string wraps around (radians)
    entry_angle_rad: float = 0.0  # Angle string enters
    exit_angle_rad: float = 0.0  # Angle string exits

    def force_on_contact(self, tension_n: float) -> Tuple[float, float]:
        """Calculate force vector on this contact point (Newtons).

        For a string wrapping around a pin:
        F = 2 * T * sin(θ/2) pointing toward center of wrap

        Returns (Fx, Fz) in Newtons.
        """
        if self.wrap_angle_rad == 0:
            return (0.0, 0.0)

        # Force magnitude
        F = 2 * tension_n * math.sin(abs(self.wrap_angle_rad) / 2)

        # Direction: bisector of entry/exit angles, pointing inward
        bisector_angle = (self.entry_angle_rad + self.exit_angle_rad) / 2
        # Force points opposite to bisector (toward center)
        Fx = -F * math.cos(bisector_angle)
        Fz = -F * math.sin(bisector_angle)

        return (Fx, Fz)

    def tangent_entry_point(self, from_x: float, from_z: float, string_radius: float = 0) -> Tuple[float, float]:
        """Calculate tangent point where string from (from_x, from_z) touches this contact.

        Returns (x, z) of tangent point on the effective radius (contact + string radius).
        """
        eff_r = self.radius_mm + string_radius
        dx = self.x_mm - from_x
        dz = self.z_mm - from_z
        dist = math.sqrt(dx*dx + dz*dz)

        if dist <= eff_r:
            return (self.x_mm, self.z_mm)

        theta = math.atan2(dz, dx)
        alpha = math.asin(eff_r / dist)

        # Left-side tangent (CCW wrap)
        tp_angle = theta + math.pi/2 - alpha
        return (self.x_mm + eff_r * math.cos(tp_angle),
                self.z_mm + eff_r * math.sin(tp_angle))


@dataclass
class StringSegment:
    """A straight segment of string between two points.

    Each segment has uniform tension along its length.
    """
    start_x: float
    start_z: float
    end_x: float
    end_z: float
    tension_n: float = 0.0

    def length_mm(self) -> float:
        """Length of this segment."""
        dx = self.end_x - self.start_x
        dz = self.end_z - self.start_z
        return math.sqrt(dx*dx + dz*dz)

    def angle_rad(self) -> float:
        """Angle of segment from horizontal (radians)."""
        return math.atan2(self.end_z - self.start_z, self.end_x - self.start_x)

    def direction(self) -> Tuple[float, float]:
        """Unit vector in direction of segment."""
        length = self.length_mm()
        if length == 0:
            return (0, 0)
        return ((self.end_x - self.start_x) / length,
                (self.end_z - self.start_z) / length)

    def force_at_start(self) -> Tuple[float, float]:
        """Force vector at start point (pulling toward end)."""
        ux, uz = self.direction()
        return (self.tension_n * ux, self.tension_n * uz)

    def force_at_end(self) -> Tuple[float, float]:
        """Force vector at end point (pulling toward start)."""
        ux, uz = self.direction()
        return (-self.tension_n * ux, -self.tension_n * uz)


@dataclass
class StringPath:
    """Complete path of a string through the harp.

    Built incrementally by adding contact points. Calculates segments,
    wrap angles, tensions, and forces.
    """
    material: StringMaterial
    target_frequency_hz: float

    contacts: List[ContactPoint] = field(default_factory=list)
    segments: List[StringSegment] = field(default_factory=list)

    # Anchor points
    soundboard_x: float = 0.0
    soundboard_z: float = 0.0
    tuning_pin_x: float = 0.0
    tuning_pin_z: float = 0.0

    def add_contact(self, contact: ContactPoint):
        """Add a contact point to the path."""
        self.contacts.append(contact)

    def build_path(self, string_radius: float = 0):
        """Build segments and calculate geometry after all contacts are added.

        Path goes: soundboard → [contacts in order] → tuning_pin
        """
        self.segments = []

        # Build list of waypoints: soundboard, contacts, tuning pin
        waypoints = [(self.soundboard_x, self.soundboard_z, None)]  # (x, z, contact_or_none)
        for c in self.contacts:
            waypoints.append((c.x_mm, c.z_mm, c))
        waypoints.append((self.tuning_pin_x, self.tuning_pin_z, None))

        # Calculate tangent points and segments
        prev_x, prev_z = waypoints[0][0], waypoints[0][1]

        for i in range(1, len(waypoints)):
            curr_x, curr_z, contact = waypoints[i]

            if contact is not None:
                # Calculate tangent entry point
                entry_x, entry_z = contact.tangent_entry_point(prev_x, prev_z, string_radius)
                contact.entry_angle_rad = math.atan2(entry_z - contact.z_mm, entry_x - contact.x_mm)

                # Segment from previous point to tangent entry
                self.segments.append(StringSegment(prev_x, prev_z, entry_x, entry_z))

                # Calculate exit tangent if there's a next point
                if i < len(waypoints) - 1:
                    next_x, next_z, _ = waypoints[i + 1]
                    # For now, simple exit calculation
                    exit_x, exit_z = contact.tangent_entry_point(next_x, next_z, string_radius)
                    contact.exit_angle_rad = math.atan2(exit_z - contact.z_mm, exit_x - contact.x_mm)
                    contact.wrap_angle_rad = contact.exit_angle_rad - contact.entry_angle_rad
                    # Normalize wrap angle
                    while contact.wrap_angle_rad > math.pi:
                        contact.wrap_angle_rad -= 2 * math.pi
                    while contact.wrap_angle_rad < -math.pi:
                        contact.wrap_angle_rad += 2 * math.pi

                    prev_x, prev_z = exit_x, exit_z
                else:
                    prev_x, prev_z = entry_x, entry_z
            else:
                # Direct segment to this point
                self.segments.append(StringSegment(prev_x, prev_z, curr_x, curr_z))
                prev_x, prev_z = curr_x, curr_z

    def total_length_mm(self) -> float:
        """Total length of string path."""
        return sum(seg.length_mm() for seg in self.segments)

    def vibrating_length_mm(self) -> float:
        """Vibrating length (soundboard to first contact that stops vibration)."""
        # For flat position, this is to the flat pin (or soundboard to tuning if no flat pin)
        if not self.segments:
            return 0
        # First segment is soundboard to first contact
        return self.segments[0].length_mm()

    def calculate_tensions(self):
        """Calculate tension based on target frequency and vibrating length."""
        vib_length = self.vibrating_length_mm()
        if vib_length <= 0:
            return

        tension = self.material.tension_for_frequency(self.target_frequency_hz, vib_length)

        # All segments have same tension (assuming frictionless contacts)
        for seg in self.segments:
            seg.tension_n = tension

    def get_tension(self) -> float:
        """Get the string tension (same for all segments)."""
        if self.segments:
            return self.segments[0].tension_n
        return 0

    def force_on_soundboard(self) -> Tuple[float, float]:
        """Force exerted on soundboard by string."""
        if not self.segments:
            return (0, 0)
        return self.segments[0].force_at_start()

    def force_on_tuning_pin(self) -> Tuple[float, float]:
        """Force exerted on tuning pin by string."""
        if not self.segments:
            return (0, 0)
        return self.segments[-1].force_at_end()

    def forces_on_contacts(self) -> List[Tuple[ContactPoint, float, float]]:
        """Get forces on all contact points."""
        tension = self.get_tension()
        result = []
        for c in self.contacts:
            fx, fz = c.force_on_contact(tension)
            result.append((c, fx, fz))
        return result

    def safety_analysis(self) -> Dict:
        """Analyze safety factors and breaking risk."""
        tension = self.get_tension()
        return {
            'tension_n': tension,
            'breaking_tension_n': self.material.breaking_tension(),
            'safety_factor': self.material.safety_factor(tension),
            'percent_of_breaking': 100 * tension / self.material.breaking_tension() if self.material.breaking_tension() > 0 else 0
        }


@dataclass
class VibrationAnalysis:
    """Analyze string vibration and clearances.

    When a string vibrates, it forms a standing wave pattern.
    This class calculates the vibration envelope and clearances to nearby objects.
    """
    string_path: StringPath

    def amplitude_at_position(self, distance_from_soundboard_mm: float,
                              max_amplitude_mm: float = 2.0) -> float:
        """Calculate vibration amplitude at a position along the string.

        For fundamental mode: amplitude = A_max * sin(π * x / L)
        where x is distance from fixed end, L is vibrating length.
        """
        L = self.string_path.vibrating_length_mm()
        if L <= 0 or distance_from_soundboard_mm < 0 or distance_from_soundboard_mm > L:
            return 0

        return max_amplitude_mm * math.sin(math.pi * distance_from_soundboard_mm / L)

    def max_amplitude_position(self) -> float:
        """Position of maximum amplitude (center of vibrating length)."""
        return self.string_path.vibrating_length_mm() / 2

    def clearance_to_point(self, point_x: float, point_z: float,
                           point_radius: float = 0,
                           max_amplitude_mm: float = 2.0) -> Dict:
        """Calculate clearance between vibrating string and a point (e.g., prong).

        Returns dict with:
        - static_clearance: distance when string is at rest
        - min_clearance: minimum distance during vibration
        - touches: whether string touches the point during vibration
        """
        # Find closest point on string path to the given point
        # For now, simplified: assume string is a straight line from soundboard
        if not self.string_path.segments:
            return {'static_clearance': float('inf'), 'min_clearance': float('inf'), 'touches': False}

        seg = self.string_path.segments[0]  # First segment (soundboard toward first contact)

        # Project point onto segment line
        ux, uz = seg.direction()
        dx = point_x - seg.start_x
        dz = point_z - seg.start_z

        # Distance along segment
        t = dx * ux + dz * uz
        t = max(0, min(seg.length_mm(), t))

        # Closest point on string
        closest_x = seg.start_x + t * ux
        closest_z = seg.start_z + t * uz

        # Perpendicular distance (static clearance)
        perp_dist = math.sqrt((point_x - closest_x)**2 + (point_z - closest_z)**2)
        static_clearance = perp_dist - point_radius

        # Amplitude at this position
        amplitude = self.amplitude_at_position(t, max_amplitude_mm)

        # Minimum clearance during vibration
        min_clearance = static_clearance - amplitude

        return {
            'static_clearance': static_clearance,
            'min_clearance': min_clearance,
            'amplitude_at_point': amplitude,
            'touches': min_clearance <= 0,
            'distance_along_string': t
        }

    def prong_engagement_analysis(self, prong_x: float, prong_z: float,
                                   prong_radius: float,
                                   disc_rotation_deg: float) -> Dict:
        """Analyze what happens as a disc rotates toward the string.

        Returns analysis of when prong contacts string, deflection amount, etc.
        """
        clearance = self.clearance_to_point(prong_x, prong_z, prong_radius)

        return {
            'clearance': clearance,
            'disc_rotation_deg': disc_rotation_deg,
            'prong_engaged': clearance['static_clearance'] <= 0,
            'would_touch_during_vibration': clearance['touches']
        }


# =============================================================================
# PHYSICS ANALYSIS FUNCTIONS
# =============================================================================

def create_string_path(string: 'String', include_flat_pin: bool = True) -> StringPath:
    """Create a StringPath from a String object for physics analysis.

    This builds the physical path of the string through the harp:
    - Starts at soundboard
    - Goes through flat pin (if include_flat_pin=True)
    - Ends at tuning pin

    Args:
        string: The String object with hardware positions
        include_flat_pin: If True, include flat pin in path; if False, direct to tuning pin

    Returns:
        StringPath ready for physics analysis
    """
    # Set material properties based on core material type
    core = string.core_material.lower() if string.core_material else "steel"

    if "nylon" in core or "gut" in core:
        # Nylon/gut string properties
        density = 1150.0  # kg/m³ for nylon
        youngs_mod = 3.0  # GPa for nylon
        tensile = 80.0  # MPa for nylon (much lower than steel)
    else:
        # Steel/metal core properties
        density = 7850.0  # kg/m³ for steel
        youngs_mod = 200.0  # GPa for steel
        tensile = 2000.0  # MPa for music wire

    # Create material from string properties
    material = StringMaterial(
        diameter_mm=string.outer_diameter_mm,
        core_material=string.core_material,
        wrap_material=string.wrap_material,
        density_kg_m3=density,
        youngs_modulus_gpa=youngs_mod,
        tensile_strength_mpa=tensile,
    )

    # Create the path
    path = StringPath(
        material=material,
        target_frequency_hz=string.frequency_hz,
        soundboard_x=string.x_soundboard_mm,
        soundboard_z=string.z_soundboard_mm,
        tuning_pin_x=string.tuning_pin.x_mm,
        tuning_pin_z=string.tuning_pin.z_mm,
    )

    # Add flat pin as a contact point
    if include_flat_pin:
        flat_contact = ContactPoint(
            x_mm=string.flat_pin.x_mm,
            z_mm=string.flat_pin.z_mm,
            radius_mm=string.flat_pin.radius_mm(),
            contact_type="flat_pin"
        )
        path.add_contact(flat_contact)

    # Build the path geometry
    path.build_path(string_radius=string.outer_diameter_mm / 2)

    # Calculate tensions
    path.calculate_tensions()

    return path


def analyze_string_physics(string: 'String') -> Dict:
    """Perform complete physics analysis on a string.

    Returns dict with:
    - path: The StringPath object
    - tension: Calculated tension in Newtons
    - forces: Forces on each contact point
    - safety: Safety factor analysis
    - vibration: Vibration analysis at key points
    """
    path = create_string_path(string)

    # Safety analysis using design tension from JSON (more reliable than calculated)
    # The calculated tension from physics may differ from actual string behavior
    design_tension = string.tension_n
    calculated_tension = path.get_tension()

    # Use design tension for safety (it's what the string will actually experience)
    safety = {
        'tension_n': design_tension,
        'calculated_tension_n': calculated_tension,
        'breaking_tension_n': path.material.breaking_tension(),
        'safety_factor': path.material.breaking_tension() / design_tension if design_tension > 0 else float('inf'),
        'percent_of_breaking': 100 * design_tension / path.material.breaking_tension() if path.material.breaking_tension() > 0 else 0,
        'tension_ratio': calculated_tension / design_tension if design_tension > 0 else 0,
    }

    # Vibration analysis
    vib = VibrationAnalysis(path)

    # Analyze clearance to disc prongs (when in flat position)
    prong_clearances = []

    # Natural disc prongs
    for prong in string.natural_disc.prongs(rotation_override=0):  # Flat position
        clearance = vib.clearance_to_point(
            prong.x_mm, prong.z_mm,
            prong.radius_mm(),
            max_amplitude_mm=2.0
        )
        prong_clearances.append({
            'disc': 'natural',
            'prong_x': prong.x_mm,
            'prong_z': prong.z_mm,
            **clearance
        })

    # Sharp disc prongs
    for prong in string.sharp_disc.prongs(rotation_override=0):  # Flat position
        clearance = vib.clearance_to_point(
            prong.x_mm, prong.z_mm,
            prong.radius_mm(),
            max_amplitude_mm=2.0
        )
        prong_clearances.append({
            'disc': 'sharp',
            'prong_x': prong.x_mm,
            'prong_z': prong.z_mm,
            **clearance
        })

    return {
        'path': path,
        'total_length_mm': path.total_length_mm(),
        'vibrating_length_mm': path.vibrating_length_mm(),
        'tension_n': path.get_tension(),
        'force_on_soundboard': path.force_on_soundboard(),
        'force_on_tuning_pin': path.force_on_tuning_pin(),
        'forces_on_contacts': path.forces_on_contacts(),
        'safety': safety,
        'max_vibration_amplitude_position_mm': vib.max_amplitude_position(),
        'prong_clearances': prong_clearances,
    }


def print_physics_analysis(harp: 'Harp'):
    """Print physics analysis report for all strings."""
    print("\n" + "=" * 70)
    print("STRING PHYSICS ANALYSIS")
    print("=" * 70)

    total_soundboard_fx = 0.0
    total_soundboard_fz = 0.0
    total_tuning_fx = 0.0
    total_tuning_fz = 0.0
    total_tension = 0.0

    min_safety_factor = float('inf')
    min_safety_string = None

    clearance_warnings = []

    for s in harp.strings:
        analysis = analyze_string_physics(s)

        # Accumulate forces
        sb_fx, sb_fz = analysis['force_on_soundboard']
        tp_fx, tp_fz = analysis['force_on_tuning_pin']
        total_soundboard_fx += sb_fx
        total_soundboard_fz += sb_fz
        total_tuning_fx += tp_fx
        total_tuning_fz += tp_fz
        total_tension += analysis['tension_n']

        # Track minimum safety factor
        sf = analysis['safety']['safety_factor']
        if sf < min_safety_factor:
            min_safety_factor = sf
            min_safety_string = s

        # Check prong clearances
        for pc in analysis['prong_clearances']:
            if pc['min_clearance'] < 1.0:  # Less than 1mm clearance
                clearance_warnings.append({
                    'string': s.number,
                    'note': s.note,
                    'disc': pc['disc'],
                    'static_clearance': pc['static_clearance'],
                    'min_clearance': pc['min_clearance'],
                })

    # Print summary
    print(f"\nTotal tension across all strings: {total_tension:.1f} N ({total_tension/4.448:.1f} lbf)")

    sb_total = math.sqrt(total_soundboard_fx**2 + total_soundboard_fz**2)
    print(f"\nForce on soundboard:")
    print(f"  X component: {total_soundboard_fx:.1f} N")
    print(f"  Z component: {total_soundboard_fz:.1f} N")
    print(f"  Total:       {sb_total:.1f} N ({sb_total/4.448:.1f} lbf)")

    tp_total = math.sqrt(total_tuning_fx**2 + total_tuning_fz**2)
    print(f"\nForce on tuning pins (total):")
    print(f"  X component: {total_tuning_fx:.1f} N")
    print(f"  Z component: {total_tuning_fz:.1f} N")
    print(f"  Total:       {tp_total:.1f} N ({tp_total/4.448:.1f} lbf)")

    print(f"\nSafety factors:")
    print(f"  Minimum:  {min_safety_factor:.1f}x (string {min_safety_string.number} - {min_safety_string.note})")
    print(f"  Typical recommended: 3-5x for music wire")

    if clearance_warnings:
        print(f"\nClearance warnings (prongs within 1mm of vibrating string):")
        for w in clearance_warnings:
            print(f"  String {w['string']} ({w['note']}): {w['disc']} disc - "
                  f"static {w['static_clearance']:.2f}mm, min {w['min_clearance']:.2f}mm")
    else:
        print(f"\nAll prong clearances OK (>1mm from vibrating string)")

    print("=" * 70)
