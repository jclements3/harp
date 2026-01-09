#!/usr/bin/env python3
"""
HARP_RENDERER.PY - SVG rendering for harp designs

Contains:
- HarpRenderer: Complete SVG renderer with disc, string, and force vector drawing
"""

import math
from typing import Tuple

import svgwrite

from harp_models import Harp, Disc, NATURAL_ROTATION_DEG, SHARP_ROTATION_DEG


class HarpRenderer:
    """Renders a Harp to SVG using svgwrite."""

    # Plate color palettes: +Y (left/blue) and -Y (right/warm)
    PLATE_COLORS = {
        "+Y": {
            "disc_fill": "#A8C4E8",      # Light blue
            "disc_stroke": "#4A6B8B",    # Dark blue
            "prong_fill": "#6688AA",     # Medium blue
            "prong_stroke": "#334455",   # Dark blue-gray
            "pin_fill": "#7799BB",       # Blue-gray
            "pin_stroke": "#445566",     # Dark blue-gray
        },
        "-Y": {
            "disc_fill": "#E8C4A8",      # Light orange/peach
            "disc_stroke": "#8B6B4A",    # Dark orange/brown
            "prong_fill": "#AA8866",     # Medium orange
            "prong_stroke": "#554433",   # Dark brown
            "pin_fill": "#BB9977",       # Orange-brown
            "pin_stroke": "#665544",     # Dark brown
        },
    }

    # Legacy disc colors (fallback)
    DISC_COLORS = {
        "natural": ("#B8A8D4", "#55438B"),
        "sharp": ("#D4C4A8", "#8B7355"),
        "double_flat": ("#E8D4A8", "#8B7355"),
        "double_sharp": ("#C4A8E8", "#55438B"),
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

    def render(self, output_path: str, pedal_position: str = "flat",
                show_discs: bool = True, show_flat_pins: bool = True,
                show_reference_lines: bool = True, hardware_opacity: float = 1.0,
                string_path_mode: str = "normal", show_force_vectors: bool = False,
                force_scale: float = 0.1):
        """Render the harp to an SVG file.

        Args:
            output_path: Path to save the SVG file
            pedal_position: "flat", "natural", or "sharp"
                - flat: all discs at 0° rotation, strings straight
                - natural: natural disc at 45°, strings deflect at natural engagement
                - sharp: natural at 45°, sharp at 90°, strings deflect at both
            show_discs: If True, draw discs and prongs
            show_flat_pins: If True, draw flat pins
            show_reference_lines: If True, draw reference lines (neckref, flat pin line)
            hardware_opacity: Opacity for non-essential hardware (0.0-1.0)
            string_path_mode: "normal" (wraps flat pin) or "direct" (soundboard to tuning pin only)
            show_force_vectors: If True, draw force vectors at soundboard and tuning pins
            force_scale: Scale factor for force vector length (pixels per Newton)
        """

        # Use viewBox to zoom in on upper part (discs and pegs)
        # Pegs are around y=44, discs around y=100-150, so show y=0 to y=200
        view_height = 200
        aspect_ratio = self.width / view_height
        height_vw = 100 / aspect_ratio  # height in vw to maintain aspect ratio
        dwg = svgwrite.Drawing(output_path,
                               size=('100vw', f'{height_vw:.1f}vw'),
                               viewBox=f'0 0 {self.width} {view_height}',
                               preserveAspectRatio='xMinYMin meet',
                               debug=False)  # Disable validation for vw units

        # White background
        dwg.add(dwg.rect((0, 0), (self.width, self.height), fill='white'))

        # Reference lines (soundboard always drawn, other lines optional)
        self._draw_reference_lines(dwg, show_all=show_reference_lines, opacity=hardware_opacity)

        # Draw discs at appropriate rotation
        if show_discs:
            self._draw_discs(dwg, pedal_position, opacity=hardware_opacity)

        # Draw strings with deflection at engagement points
        self._draw_strings(dwg, pedal_position, show_flat_pins=show_flat_pins,
                          flat_pin_opacity=hardware_opacity, string_path_mode=string_path_mode)

        # Draw force vectors
        if show_force_vectors:
            self._draw_force_vectors(dwg, force_scale, string_path_mode)

        dwg.save()
        print(f"SVG saved to {output_path} (pedal position: {pedal_position})")

    def _draw_reference_lines(self, dwg, show_all: bool = True, opacity: float = 1.0):
        """Draw reference lines."""
        c1 = self.harp.strings[0]
        g7 = self.harp.strings[-1]

        if show_all:
            # Green dashed line: flat pins c1 to g7
            dwg.add(dwg.line(
                (self._tx(c1.x_flat_mm), self._tz(c1.z_flat_mm)),
                (self._tx(g7.x_flat_mm), self._tz(g7.z_flat_mm)),
                stroke='#90EE90',
                stroke_width=0.5,
                stroke_dasharray='4,2',
                opacity=opacity
            ))

            # Brown solid line: neckref from c1 tuning pin to g7 tuning pin (bottom of pins)
            c1_tp = c1.tuning_pin
            g7_tp = g7.tuning_pin
            dwg.add(dwg.line(
                (self._tx(c1_tp.x_mm), self._tz(c1_tp.z_mm - c1_tp.radius_mm())),
                (self._tx(g7_tp.x_mm), self._tz(g7_tp.z_mm - g7_tp.radius_mm())),
                stroke='#8B4513',
                stroke_width=0.8,
                opacity=opacity
            ))

        # Black thin line: soundboard spline (through string bottom positions) - always draw
        soundboard_points = []
        for s in self.harp.strings:
            soundboard_points.append((self._tx(s.x_soundboard_mm), self._tz(s.z_soundboard_mm)))

        if len(soundboard_points) >= 2:
            path_data = f"M {soundboard_points[0][0]:.1f},{soundboard_points[0][1]:.1f}"
            for px, pz in soundboard_points[1:]:
                path_data += f" L {px:.1f},{pz:.1f}"
            dwg.add(dwg.path(d=path_data, fill='none', stroke='black', stroke_width=0.5))

    def _draw_disc(self, dwg, disc: Disc, rotation_override: float = None, opacity: float = 1.0):
        """Draw a single disc with prongs using stadium shape."""
        pedal_rotation = rotation_override if rotation_override is not None else disc.rotation_degrees

        svg_dx = disc.string_ux * self.scale
        svg_dy = -disc.string_uz * self.scale

        string_angle_svg = math.degrees(math.atan2(svg_dy, svg_dx))
        base_rotation = string_angle_svg + 90.0
        total_rotation = base_rotation + pedal_rotation

        plate_colors = self.PLATE_COLORS.get(disc.plate, self.PLATE_COLORS["+Y"])
        fill = plate_colors["disc_fill"]
        stroke = plate_colors["disc_stroke"]
        prong_fill = plate_colors["prong_fill"]
        prong_stroke = plate_colors["prong_stroke"]

        if pedal_rotation > 0:
            opacity = opacity * 0.5

        cx_svg = self._tx(disc.x_mm)
        cy_svg = self._tz(disc.z_mm)

        disc_group = dwg.g(transform=f"rotate({total_rotation:.1f}, {cx_svg:.1f}, {cy_svg:.1f})",
                          opacity=opacity)

        prong_r_mm = disc.prong_diameter_mm / 2
        inset_distance = disc.major_radius_mm - disc.prong_diameter_mm
        arc_r_mm = disc.minor_radius_mm

        prong_r = self._scale(prong_r_mm)
        arc_r = self._scale(arc_r_mm)
        inset = self._scale(inset_distance)

        p1_cx = cx_svg + inset
        p1_cy = cy_svg
        p2_cx = cx_svg - inset
        p2_cy = cy_svg

        path_data = (
            f"M {p1_cx:.2f},{p1_cy - arc_r:.2f} "
            f"A {arc_r:.2f},{arc_r:.2f} 0 0,1 {p1_cx:.2f},{p1_cy + arc_r:.2f} "
            f"L {p2_cx:.2f},{p2_cy + arc_r:.2f} "
            f"A {arc_r:.2f},{arc_r:.2f} 0 0,1 {p2_cx:.2f},{p2_cy - arc_r:.2f} "
            f"L {p1_cx:.2f},{p1_cy - arc_r:.2f} "
            "Z"
        )

        disc_group.add(dwg.path(d=path_data, fill=fill, stroke=stroke, stroke_width=0.15))

        disc_group.add(dwg.circle(center=(p1_cx, p1_cy), r=prong_r,
                                   fill=prong_fill, stroke=prong_stroke, stroke_width=0.1))
        disc_group.add(dwg.circle(center=(p2_cx, p2_cy), r=prong_r,
                                   fill=prong_fill, stroke=prong_stroke, stroke_width=0.1))

        dwg.add(disc_group)

    def _draw_discs(self, dwg, pedal_position: str = "flat", opacity: float = 1.0):
        """Draw all discs for all strings."""
        for s in self.harp.strings:
            if pedal_position == "flat":
                self._draw_disc(dwg, s.natural_disc, rotation_override=0, opacity=opacity)
                self._draw_disc(dwg, s.sharp_disc, rotation_override=0, opacity=opacity)
            elif pedal_position == "natural":
                self._draw_disc(dwg, s.natural_disc, rotation_override=NATURAL_ROTATION_DEG, opacity=opacity)
                self._draw_disc(dwg, s.sharp_disc, rotation_override=0, opacity=opacity)
            elif pedal_position == "sharp":
                self._draw_disc(dwg, s.natural_disc, rotation_override=NATURAL_ROTATION_DEG, opacity=opacity)
                self._draw_disc(dwg, s.sharp_disc, rotation_override=SHARP_ROTATION_DEG, opacity=opacity)

    def _tangent_point_left(self, px: float, pz: float, cx: float, cz: float, r: float) -> Tuple[float, float]:
        """Calculate LEFT side tangent point on circle."""
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
        """Calculate RIGHT side tangent point on circle."""
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
        """Calculate tangent point - defaults to LEFT side."""
        return self._tangent_point_left(px, pz, cx, cz, r)

    def _external_tangent_left(self, c1x: float, c1z: float, r1: float,
                                c2x: float, c2z: float, r2: float) -> Tuple[float, float, float, float]:
        """Calculate external tangent points between two circles, LEFT side approach."""
        dx = c2x - c1x
        dz = c2z - c1z
        dist = math.sqrt(dx*dx + dz*dz)

        if dist < abs(r1 - r2):
            return c1x, c1z, c2x, c2z

        theta = math.atan2(dz, dx)
        if dist > 0:
            beta = math.asin((r2 - r1) / dist) if abs(r2 - r1) < dist else 0
        else:
            beta = 0

        perp_angle = theta + math.pi/2 + beta
        t1x = c1x + r1 * math.cos(perp_angle)
        t1z = c1z + r1 * math.sin(perp_angle)
        t2x = c2x + r2 * math.cos(perp_angle)
        t2z = c2z + r2 * math.sin(perp_angle)

        return t1x, t1z, t2x, t2z

    def _external_tangent_right(self, c1x: float, c1z: float, r1: float,
                                 c2x: float, c2z: float, r2: float) -> Tuple[float, float, float, float]:
        """Calculate external tangent points between two circles, RIGHT side approach."""
        dx = c2x - c1x
        dz = c2z - c1z
        dist = math.sqrt(dx*dx + dz*dz)

        if dist < abs(r1 - r2):
            return c1x, c1z, c2x, c2z

        theta = math.atan2(dz, dx)
        if dist > 0:
            beta = math.asin((r2 - r1) / dist) if abs(r2 - r1) < dist else 0
        else:
            beta = 0

        perp_angle = theta - math.pi/2 - beta
        t1x = c1x + r1 * math.cos(perp_angle)
        t1z = c1z + r1 * math.sin(perp_angle)
        t2x = c2x + r2 * math.cos(perp_angle)
        t2z = c2z + r2 * math.sin(perp_angle)

        return t1x, t1z, t2x, t2z

    def _get_engaged_prong_center(self, disc: Disc, pedal_rotation: float) -> Tuple[float, float]:
        """Get the world coordinates of the engaging prong center when disc is rotated."""
        return disc.engaged_prong_position(pedal_rotation)

    def _calculate_string_wrap_both_prongs(self, from_x: float, from_z: float,
                                            disc: Disc, pedal_rotation: float,
                                            to_x: float, to_z: float, string_r: float):
        """Calculate tangent points for string grabbed by BOTH prongs of a disc."""
        prong_r = disc.prong_radius()
        effective_r = prong_r + string_r

        prongs = disc.prongs(rotation_override=pedal_rotation)
        lower_prong = prongs[0]
        upper_prong = prongs[1]

        upper_entry_x, upper_entry_z = self._tangent_point_left(
            from_x, from_z, upper_prong.x_mm, upper_prong.z_mm, effective_r)
        upper_exit_x, upper_exit_z = self._tangent_point_right(
            lower_prong.x_mm, lower_prong.z_mm, upper_prong.x_mm, upper_prong.z_mm, effective_r)

        lower_entry_x, lower_entry_z = self._tangent_point_left(
            upper_prong.x_mm, upper_prong.z_mm, lower_prong.x_mm, lower_prong.z_mm, effective_r)
        lower_exit_x, lower_exit_z = self._tangent_point_right(
            to_x, to_z, lower_prong.x_mm, lower_prong.z_mm, effective_r)

        return (
            (upper_entry_x, upper_entry_z, upper_exit_x, upper_exit_z, upper_prong.x_mm, upper_prong.z_mm),
            (lower_entry_x, lower_entry_z, lower_exit_x, lower_exit_z, lower_prong.x_mm, lower_prong.z_mm)
        )

    def _draw_arc_around_prong(self, dwg, prong_x: float, prong_z: float, prong_r: float,
                                entry_x: float, entry_z: float, exit_x: float, exit_z: float,
                                color: str, stroke_width: float, string_r: float):
        """Draw an arc around a prong from entry point to exit point."""
        effective_r = prong_r + string_r
        path_data = f"M {self._tx(entry_x):.1f},{self._tz(entry_z):.1f} "
        path_data += f"A {self._scale(effective_r):.1f},{self._scale(effective_r):.1f} 0 0,0 "
        path_data += f"{self._tx(exit_x):.1f},{self._tz(exit_z):.1f}"
        dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=stroke_width))

    def _draw_strings(self, dwg, pedal_position: str = "flat", show_flat_pins: bool = True,
                      flat_pin_opacity: float = 1.0, string_path_mode: str = "normal"):
        """Draw all strings with their pins."""
        for s in self.harp.strings:
            color = s.color()
            sw = s.stroke_width(self.scale)

            flat_r = s.flat_pin.radius_mm()
            tuning_r = s.tuning_pin.radius_mm()
            string_r = s.outer_diameter_mm / 2
            nat_prong_r = s.natural_disc.prong_radius()
            sharp_prong_r = s.sharp_disc.prong_radius()

            flat_tp_x, flat_tp_z = self._tangent_point(
                s.x_soundboard_mm, s.z_soundboard_mm,
                s.x_flat_mm, s.z_flat_mm, flat_r + string_r
            )

            if string_path_mode == "direct":
                eff_tp_r = tuning_r + string_r
                tp_tangent_x, tp_tangent_z = self._tangent_point_left(
                    s.x_soundboard_mm, s.z_soundboard_mm,
                    s.tuning_pin.x_mm, s.tuning_pin.z_mm, eff_tp_r
                )

                tp_tangent_angle = math.atan2(tp_tangent_z - s.tuning_pin.z_mm,
                                               tp_tangent_x - s.tuning_pin.x_mm)
                tp_start_angle = -math.pi / 2

                tp_arc_segments = 24
                tp_end_angle = tp_tangent_angle
                if tp_end_angle < tp_start_angle:
                    tp_end_angle += 2 * math.pi

                tp_arc_path = f"M {self._tx(s.tuning_pin.x_mm + eff_tp_r * math.cos(tp_start_angle)):.2f},{self._tz(s.tuning_pin.z_mm + eff_tp_r * math.sin(tp_start_angle)):.2f}"
                for i in range(1, tp_arc_segments + 1):
                    angle = tp_start_angle + (tp_end_angle - tp_start_angle) * i / tp_arc_segments
                    x = s.tuning_pin.x_mm + eff_tp_r * math.cos(angle)
                    z = s.tuning_pin.z_mm + eff_tp_r * math.sin(angle)
                    tp_arc_path += f" L {self._tx(x):.2f},{self._tz(z):.2f}"

                path_data = tp_arc_path + f" L {self._tx(s.x_soundboard_mm):.2f},{self._tz(s.z_soundboard_mm):.2f}"
                dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=sw))

            elif pedal_position == "flat":
                eff_tp_r = tuning_r + string_r
                eff_flat_r = flat_r + string_r

                tp_exit_x, tp_exit_z, flat_entry_x, flat_entry_z = self._external_tangent_right(
                    s.tuning_pin.x_mm, s.tuning_pin.z_mm, eff_tp_r,
                    s.x_flat_mm, s.z_flat_mm, eff_flat_r
                )

                tp_exit_angle = math.atan2(tp_exit_z - s.tuning_pin.z_mm, tp_exit_x - s.tuning_pin.x_mm)
                tp_start_angle = -math.pi / 2

                flat_exit_x, flat_exit_z = self._tangent_point_left(
                    s.x_soundboard_mm, s.z_soundboard_mm,
                    s.x_flat_mm, s.z_flat_mm, eff_flat_r
                )

                flat_entry_angle = math.atan2(flat_entry_z - s.z_flat_mm, flat_entry_x - s.x_flat_mm)
                flat_exit_angle = math.atan2(flat_exit_z - s.z_flat_mm, flat_exit_x - s.x_flat_mm)

                tp_arc_segments = 24
                tp_end_angle = tp_exit_angle
                if tp_end_angle < tp_start_angle:
                    tp_end_angle += 2 * math.pi

                tp_arc_path = f"M {self._tx(s.tuning_pin.x_mm + eff_tp_r * math.cos(tp_start_angle)):.2f},{self._tz(s.tuning_pin.z_mm + eff_tp_r * math.sin(tp_start_angle)):.2f}"
                for i in range(1, tp_arc_segments + 1):
                    angle = tp_start_angle + (tp_end_angle - tp_start_angle) * i / tp_arc_segments
                    x = s.tuning_pin.x_mm + eff_tp_r * math.cos(angle)
                    z = s.tuning_pin.z_mm + eff_tp_r * math.sin(angle)
                    tp_arc_path += f" L {self._tx(x):.2f},{self._tz(z):.2f}"

                flat_arc_segments = 16
                fp_end_angle = flat_exit_angle
                if fp_end_angle < flat_entry_angle:
                    fp_end_angle += 2 * math.pi

                flat_arc_path = ""
                for i in range(1, flat_arc_segments + 1):
                    angle = flat_entry_angle + (fp_end_angle - flat_entry_angle) * i / flat_arc_segments
                    x = s.x_flat_mm + eff_flat_r * math.cos(angle)
                    z = s.z_flat_mm + eff_flat_r * math.sin(angle)
                    flat_arc_path += f" L {self._tx(x):.2f},{self._tz(z):.2f}"

                path_data = (
                    tp_arc_path +
                    f" L {self._tx(flat_entry_x):.2f},{self._tz(flat_entry_z):.2f}" +
                    flat_arc_path +
                    f" L {self._tx(s.x_soundboard_mm):.2f},{self._tz(s.z_soundboard_mm):.2f}"
                )
                dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=sw))

            elif pedal_position == "natural":
                first, second = self._calculate_string_wrap_both_prongs(
                    s.x_soundboard_mm, s.z_soundboard_mm,
                    s.natural_disc, NATURAL_ROTATION_DEG,
                    flat_tp_x, flat_tp_z, string_r
                )
                first_entry_x, first_entry_z, first_exit_x, first_exit_z, first_px, first_pz = first
                second_entry_x, second_entry_z, second_exit_x, second_exit_z, second_px, second_pz = second

                dwg.add(dwg.line(
                    (self._tx(s.x_soundboard_mm), self._tz(s.z_soundboard_mm)),
                    (self._tx(first_entry_x), self._tz(first_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                self._draw_arc_around_prong(dwg, first_px, first_pz, nat_prong_r,
                    first_entry_x, first_entry_z, first_exit_x, first_exit_z, color, sw, string_r)

                dwg.add(dwg.line(
                    (self._tx(first_exit_x), self._tz(first_exit_z)),
                    (self._tx(second_entry_x), self._tz(second_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                self._draw_arc_around_prong(dwg, second_px, second_pz, nat_prong_r,
                    second_entry_x, second_entry_z, second_exit_x, second_exit_z, color, sw, string_r)

                dwg.add(dwg.line(
                    (self._tx(second_exit_x), self._tz(second_exit_z)),
                    (self._tx(flat_tp_x), self._tz(flat_tp_z)),
                    stroke=color, stroke_width=sw
                ))

            elif pedal_position == "sharp":
                nat_prongs = s.natural_disc.prongs(rotation_override=NATURAL_ROTATION_DEG)
                nat_upper_prong = nat_prongs[1]

                sharp_first, sharp_second = self._calculate_string_wrap_both_prongs(
                    s.x_soundboard_mm, s.z_soundboard_mm,
                    s.sharp_disc, SHARP_ROTATION_DEG,
                    nat_upper_prong.x_mm, nat_upper_prong.z_mm, string_r
                )
                s1_entry_x, s1_entry_z, s1_exit_x, s1_exit_z, s1_px, s1_pz = sharp_first
                s2_entry_x, s2_entry_z, s2_exit_x, s2_exit_z, s2_px, s2_pz = sharp_second

                nat_first, nat_second = self._calculate_string_wrap_both_prongs(
                    s2_exit_x, s2_exit_z,
                    s.natural_disc, NATURAL_ROTATION_DEG,
                    flat_tp_x, flat_tp_z, string_r
                )
                n1_entry_x, n1_entry_z, n1_exit_x, n1_exit_z, n1_px, n1_pz = nat_first
                n2_entry_x, n2_entry_z, n2_exit_x, n2_exit_z, n2_px, n2_pz = nat_second

                dwg.add(dwg.line(
                    (self._tx(s.x_soundboard_mm), self._tz(s.z_soundboard_mm)),
                    (self._tx(s1_entry_x), self._tz(s1_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                self._draw_arc_around_prong(dwg, s1_px, s1_pz, sharp_prong_r,
                    s1_entry_x, s1_entry_z, s1_exit_x, s1_exit_z, color, sw, string_r)

                dwg.add(dwg.line(
                    (self._tx(s1_exit_x), self._tz(s1_exit_z)),
                    (self._tx(s2_entry_x), self._tz(s2_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                self._draw_arc_around_prong(dwg, s2_px, s2_pz, sharp_prong_r,
                    s2_entry_x, s2_entry_z, s2_exit_x, s2_exit_z, color, sw, string_r)

                dwg.add(dwg.line(
                    (self._tx(s2_exit_x), self._tz(s2_exit_z)),
                    (self._tx(n1_entry_x), self._tz(n1_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                self._draw_arc_around_prong(dwg, n1_px, n1_pz, nat_prong_r,
                    n1_entry_x, n1_entry_z, n1_exit_x, n1_exit_z, color, sw, string_r)

                dwg.add(dwg.line(
                    (self._tx(n1_exit_x), self._tz(n1_exit_z)),
                    (self._tx(n2_entry_x), self._tz(n2_entry_z)),
                    stroke=color, stroke_width=sw
                ))

                self._draw_arc_around_prong(dwg, n2_px, n2_pz, nat_prong_r,
                    n2_entry_x, n2_entry_z, n2_exit_x, n2_exit_z, color, sw, string_r)

                dwg.add(dwg.line(
                    (self._tx(n2_exit_x), self._tz(n2_exit_z)),
                    (self._tx(flat_tp_x), self._tz(flat_tp_z)),
                    stroke=color, stroke_width=sw
                ))

            if pedal_position != "flat":
                arc_end_x = flat_tp_x + 0.5
                arc_end_z = flat_tp_z + 1.5
                path_data = f"M {self._tx(flat_tp_x):.1f},{self._tz(flat_tp_z):.1f} "
                path_data += f"A {self._scale(flat_r):.1f},{self._scale(flat_r):.1f} 0 0,1 "
                path_data += f"{self._tx(arc_end_x):.1f},{self._tz(arc_end_z):.1f}"
                dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=sw))

                dwg.add(dwg.line(
                    (self._tx(arc_end_x), self._tz(arc_end_z)),
                    (self._tx(s.tuning_pin.x_mm - tuning_r * 0.7), self._tz(s.tuning_pin.z_mm + tuning_r * 0.7)),
                    stroke=color, stroke_width=sw
                ))

                tp_arc_start_x = s.tuning_pin.x_mm - tuning_r * 0.7
                tp_arc_start_z = s.tuning_pin.z_mm + tuning_r * 0.7
                tp_arc_end_x = tp_arc_start_x + 0.3
                tp_arc_end_z = tp_arc_start_z + 1.0
                path_data = f"M {self._tx(tp_arc_start_x):.1f},{self._tz(tp_arc_start_z):.1f} "
                path_data += f"A {self._scale(tuning_r):.1f},{self._scale(tuning_r):.1f} 0 0,1 "
                path_data += f"{self._tx(tp_arc_end_x):.1f},{self._tz(tp_arc_end_z):.1f}"
                dwg.add(dwg.path(d=path_data, fill='none', stroke=color, stroke_width=sw))

            if show_flat_pins:
                fp_plate = s.tuning_pin.plate
                fp_colors = self.PLATE_COLORS.get(fp_plate, self.PLATE_COLORS["+Y"])
                dwg.add(dwg.circle(
                    center=(self._tx(s.x_flat_mm), self._tz(s.z_flat_mm)),
                    r=self._scale(flat_r),
                    fill=fp_colors["pin_fill"],
                    stroke=fp_colors["pin_stroke"],
                    stroke_width=0.15,
                    opacity=flat_pin_opacity
                ))

            tp_plate = s.tuning_pin.plate
            tp_colors = self.PLATE_COLORS.get(tp_plate, self.PLATE_COLORS["+Y"])
            dwg.add(dwg.circle(
                center=(self._tx(s.tuning_pin.x_mm), self._tz(s.tuning_pin.z_mm)),
                r=self._scale(tuning_r),
                fill=tp_colors["pin_fill"],
                stroke=tp_colors["pin_stroke"],
                stroke_width=0.15
            ))

    def _draw_force_vectors(self, dwg, force_scale: float, string_path_mode: str):
        """Draw reaction force vectors at soundboard, pins, and pegs."""
        PEG_COLOR = '#654321'
        PIN_COLOR = '#888888'
        SB_COLOR = '#000000'

        total_sb_fx, total_sb_fz = 0.0, 0.0
        total_pin_fx, total_pin_fz = 0.0, 0.0
        total_peg_fx, total_peg_fz = 0.0, 0.0

        for s in self.harp.strings:
            T = s.tension_n

            sb_to_pin_dx = s.x_flat_mm - s.x_soundboard_mm
            sb_to_pin_dz = s.z_flat_mm - s.z_soundboard_mm
            sb_to_pin_len = math.sqrt(sb_to_pin_dx**2 + sb_to_pin_dz**2)
            sb_ux = sb_to_pin_dx / sb_to_pin_len
            sb_uz = sb_to_pin_dz / sb_to_pin_len

            sb_fx = -T * sb_ux
            sb_fz = -T * sb_uz
            total_sb_fx += sb_fx
            total_sb_fz += sb_fz

            pin_entry_ux = sb_ux
            pin_entry_uz = sb_uz
            pin_to_peg_dx = s.tuning_pin.x_mm - s.x_flat_mm
            pin_to_peg_dz = s.tuning_pin.z_mm - s.z_flat_mm
            pin_to_peg_len = math.sqrt(pin_to_peg_dx**2 + pin_to_peg_dz**2)
            pin_exit_ux = pin_to_peg_dx / pin_to_peg_len
            pin_exit_uz = pin_to_peg_dz / pin_to_peg_len

            pin_fx = T * (pin_entry_ux - pin_exit_ux)
            pin_fz = T * (pin_entry_uz - pin_exit_uz)
            total_pin_fx += pin_fx
            total_pin_fz += pin_fz

            peg_entry_ux = pin_exit_ux
            peg_entry_uz = pin_exit_uz
            peg_exit_ux = 0.0
            peg_exit_uz = -1.0

            peg_fx = T * (peg_entry_ux - peg_exit_ux)
            peg_fz = T * (peg_entry_uz - peg_exit_uz)
            total_peg_fx += peg_fx
            total_peg_fz += peg_fz

            self._draw_arrow(dwg,
                self._tx(s.x_soundboard_mm), self._tz(s.z_soundboard_mm),
                sb_fx * force_scale, -sb_fz * force_scale,
                color=SB_COLOR, width=0.5, head_size=2, opacity=0.25
            )

            self._draw_arrow(dwg,
                self._tx(s.x_flat_mm), self._tz(s.z_flat_mm),
                pin_fx * force_scale, -pin_fz * force_scale,
                color=PIN_COLOR, width=0.5, head_size=2, opacity=0.25
            )

            self._draw_arrow(dwg,
                self._tx(s.tuning_pin.x_mm), self._tz(s.tuning_pin.z_mm),
                peg_fx * force_scale, -peg_fz * force_scale,
                color=PEG_COLOR, width=0.5, head_size=2, opacity=0.25
            )

        c1 = self.harp.strings[0]
        g7 = self.harp.strings[-1]

        sb_center_x = (c1.x_soundboard_mm + g7.x_soundboard_mm) / 2
        sb_center_z = (c1.z_soundboard_mm + g7.z_soundboard_mm) / 2
        pin_center_x = (c1.x_flat_mm + g7.x_flat_mm) / 2
        pin_center_z = (c1.z_flat_mm + g7.z_flat_mm) / 2
        peg_center_x = (c1.tuning_pin.x_mm + g7.tuning_pin.x_mm) / 2
        peg_center_z = (c1.tuning_pin.z_mm + g7.tuning_pin.z_mm) / 2

        total_scale = force_scale * 0.05

        self._draw_arrow(dwg,
            self._tx(sb_center_x), self._tz(sb_center_z),
            total_sb_fx * total_scale, -total_sb_fz * total_scale,
            color=SB_COLOR, width=1.5, head_size=5
        )

        self._draw_arrow(dwg,
            self._tx(pin_center_x), self._tz(pin_center_z),
            total_pin_fx * total_scale, -total_pin_fz * total_scale,
            color=PIN_COLOR, width=1.5, head_size=5
        )

        self._draw_arrow(dwg,
            self._tx(peg_center_x), self._tz(peg_center_z),
            total_peg_fx * total_scale, -total_peg_fz * total_scale,
            color=PEG_COLOR, width=1.5, head_size=5
        )

        sb_total = math.sqrt(total_sb_fx**2 + total_sb_fz**2)
        pin_total = math.sqrt(total_pin_fx**2 + total_pin_fz**2)
        peg_total = math.sqrt(total_peg_fx**2 + total_peg_fz**2)

        label_x = self.width - 120
        label_y = self.height - 50

        dwg.add(dwg.text(f"sb: {sb_total:.0f}N ({sb_total/4.448:.0f}lbf)",
            insert=(label_x, label_y),
            font_size="9px", fill=SB_COLOR, font_family='sans-serif'
        ))

        dwg.add(dwg.text(f"pins: {pin_total:.0f}N ({pin_total/4.448:.0f}lbf)",
            insert=(label_x, label_y + 12),
            font_size="9px", fill=PIN_COLOR, font_family='sans-serif'
        ))

        dwg.add(dwg.text(f"pegs: {peg_total:.0f}N ({peg_total/4.448:.0f}lbf)",
            insert=(label_x, label_y + 24),
            font_size="9px", fill=PEG_COLOR, font_family='sans-serif'
        ))

    def _draw_arrow(self, dwg, x: float, y: float, dx: float, dy: float,
                    color: str = 'black', width: float = 1.0, head_size: float = 4,
                    opacity: float = 1.0):
        """Draw an arrow from (x,y) with direction (dx,dy)."""
        end_x = x + dx
        end_y = y + dy
        dwg.add(dwg.line((x, y), (end_x, end_y), stroke=color, stroke_width=width, opacity=opacity))

        length = math.sqrt(dx*dx + dy*dy)
        if length > 0:
            ux, uy = dx/length, dy/length
            px, py = -uy, ux

            head_back = head_size * 0.8
            head_width = head_size * 0.4
            p1 = (end_x - head_back*ux + head_width*px, end_y - head_back*uy + head_width*py)
            p2 = (end_x - head_back*ux - head_width*px, end_y - head_back*uy - head_width*py)

            points = [(end_x, end_y), p1, p2]
            dwg.add(dwg.polygon(points, fill=color, stroke=color, stroke_width=0.5, opacity=opacity))
