#!/usr/bin/env python3
"""
HARP_MODELS.PY - Data classes for harp components

Contains all the data structures representing harp components:
- Prong, Disc, FlatPin, TuningPin
- String, Neck, Soundboard, Harp
"""

import math
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional


# =============================================================================
# DISC ROTATION CONSTANTS
# =============================================================================

NATURAL_ROTATION_DEG = 45.0  # 3 o'clock to ~4:30
SHARP_ROTATION_DEG = 90.0    # 3 o'clock to 6 o'clock


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

    def prongs(self, rotation_override: float = None) -> List['Prong']:
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
