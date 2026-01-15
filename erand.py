#!/usr/bin/env python3
"""
Erard Concert Harp - Parametric CadQuery Model

Computes string and hardware positions from:
- sp: soundboard path (line from c1b to g7b)
- finger_gap: perpendicular spacing between strings
- per-string: vib_length (to flat pin), diameter

Positions *n and *s are derived from *f using semitone ratio.
"""

import json
import math
from pathlib import Path

import cadquery as cq
from cq_plate import plate
from disc_assembly import (
    make_disc_assembly, DISC_SPECS,
    make_dual_position_disc, DUAL_POSITION_SPECS,
    get_scoop_z_offset, get_dual_scoop_z_offset
)


# =============================================================================
# Parametric Geometry Computation
# =============================================================================

def compute_sp(config):
    """Compute soundboard path from config.

    Returns (start, end, unit_vector, length, extension).
    sp passes through c1b to g7b, extended by 2 finger gaps on each end.
    """
    c1b = config['sp']['c1b']
    g7b = config['sp']['g7b']
    finger_gap = config['geometry']['finger_gap']

    dx = g7b['x'] - c1b['x']
    dy = g7b['y'] - c1b['y']
    length = math.sqrt(dx*dx + dy*dy)
    ux, uy = dx/length, dy/length

    extension = 2 * finger_gap
    start = (c1b['x'] - ux * extension, c1b['y'] - uy * extension)
    end = (g7b['x'] + ux * extension, g7b['y'] + uy * extension)

    return {
        'start': start,
        'end': end,
        'unit': (ux, uy),
        'length': length,
        'extension': extension,
        'c1b': (c1b['x'], c1b['y']),
        'g7b': (g7b['x'], g7b['y'])
    }


def compute_b_positions(sp, num_strings, finger_gap):
    """Compute *b positions along sp with equal finger gap spacing.

    Returns list of (x, y) tuples for each string's soundboard attachment.
    String 1 at c1b, string N at position determined by equal spacing.
    """
    ux, uy = sp['unit']
    c1b = sp['c1b']

    positions = []
    for i in range(num_strings):
        t = i * finger_gap
        x = c1b[0] + ux * t
        y = c1b[1] + uy * t
        positions.append((x, y))

    return positions


def compute_string_geometry(b, vib_length, string_tilt_deg, semitone_ratio):
    """Compute *f, *n, *s positions for a single string.

    Args:
        b: (x, y) soundboard attachment point
        vib_length: distance from *b to *f (flat pin, determines pitch)
        string_tilt_deg: angle from vertical (0 = vertical)
        semitone_ratio: 2^(1/12) for computing *n, *s

    Returns dict with positions for f, n, s (flat pin, natural disc, sharp disc).
    """
    # String direction (tilt from vertical)
    tilt_rad = math.radians(string_tilt_deg)
    # For vertical strings going up: dx=0, dy=1
    # Tilt rotates this: dx=sin(tilt), dy=cos(tilt)
    string_dx = math.sin(tilt_rad)
    string_dy = math.cos(tilt_rad)

    # *f position (flat pin) - defines vibrating length
    f_x = b[0] + string_dx * vib_length
    f_y = b[1] + string_dy * vib_length

    # *n position (natural disc) - one semitone shorter
    vib_n = vib_length / semitone_ratio
    n_x = b[0] + string_dx * vib_n
    n_y = b[1] + string_dy * vib_n

    # *s position (sharp disc) - two semitones shorter
    vib_s = vib_length / (semitone_ratio ** 2)
    s_x = b[0] + string_dx * vib_s
    s_y = b[1] + string_dy * vib_s

    return {
        'b': b,
        'f': (f_x, f_y),
        'n': (n_x, n_y),
        's': (s_x, s_y),
        'vib_f': vib_length,
        'vib_n': vib_n,
        'vib_s': vib_s,
        'string_dir': (string_dx, string_dy)
    }


def get_hardware_spec(string_num, hardware):
    """Get hardware dimensions for a string number."""
    # Flat pins
    if string_num <= 9:
        pin_d = hardware['flat_pins']['bass']['shaft_diameter_mm']
    elif string_num <= 28:
        pin_d = hardware['flat_pins']['mid']['shaft_diameter_mm']
    else:
        pin_d = hardware['flat_pins']['treble']['shaft_diameter_mm']

    # Tuning pegs
    if string_num <= 9:
        peg_d = hardware['tuning_pins']['bass']['diameter_mm']
    elif string_num <= 28:
        peg_d = hardware['tuning_pins']['mid']['diameter_mm']
    else:
        peg_d = hardware['tuning_pins']['treble']['diameter_mm']

    # Discs (4 ranges for monolithic design)
    if string_num <= 9:
        disc = hardware['discs']['bass']
    elif string_num <= 28:
        disc = hardware['discs']['mid']
    elif string_num <= 38:
        disc = hardware['discs']['treble']
    else:
        disc = hardware['discs']['high']

    return {
        'pin_diameter': pin_d,
        'peg_diameter': peg_d,
        'disc_prong': disc['prong_diameter_mm'],
        'disc_thickness': disc['thickness_mm'],
        'disc_major_radius': disc['major_radius_mm'],
        'disc_minor_radius': disc['minor_radius_mm'],
        'disc_sphere_radius': disc['sphere_radius_mm'],
        'disc_prong_length': disc['prong_length_mm']
    }


def compute_hardware_positions(string_geom, string_num, peg_offset, hardware):
    """Compute pin, peg, and disc positions for a string.

    Args:
        string_geom: output from compute_string_geometry
        string_num: 1-47
        peg_offset: {'dx': ..., 'dy': ...} offset from pin to peg
        hardware: hardware specs from config

    Returns dict with pin, peg, ndisc, sdisc positions and dimensions.
    """
    f = string_geom['f']
    n = string_geom['n']
    s = string_geom['s']

    hw = get_hardware_spec(string_num, hardware)

    # Pin at *f position
    pin = {
        'x': f[0],
        'y': f[1],
        'diameter': hw['pin_diameter'],
        'center_x': f[0] + hw['pin_diameter']/2,
        'center_y': f[1]
    }

    # Peg offset from pin
    peg = {
        'x': f[0] + peg_offset['dx'],
        'y': f[1] + peg_offset['dy'],
        'diameter': hw['peg_diameter'],
        'center_x': f[0] + peg_offset['dx'] + hw['peg_diameter']/2,
        'center_y': f[1] + peg_offset['dy'],
        'tangent_x': f[0] + peg_offset['dx'],
        'tangent_y': f[1] + peg_offset['dy']
    }

    # Discs at *n and *s (monolithic design with per-register sizing)
    # Disc z-offset: center at -(sphere_radius + thickness/2)
    # Prongs extend from disc face toward string at z=0
    disc_rotation = -70.3  # degrees
    sphere_r = hw['disc_sphere_radius']
    disc_z = -(sphere_r + hw['disc_thickness'] / 2)

    ndisc = {
        'x': n[0],
        'y': n[1],
        'z': disc_z,
        'major_radius': hw['disc_major_radius'],
        'minor_radius': hw['disc_minor_radius'],
        'prong_diameter': hw['disc_prong'],
        'prong_length': hw['disc_prong_length'],
        'thickness': hw['disc_thickness'],
        'sphere_radius': sphere_r,
        'rotation_deg': disc_rotation
    }

    sdisc = {
        'x': s[0],
        'y': s[1],
        'z': disc_z,
        'major_radius': hw['disc_major_radius'],
        'minor_radius': hw['disc_minor_radius'],
        'prong_diameter': hw['disc_prong'],
        'prong_length': hw['disc_prong_length'],
        'thickness': hw['disc_thickness'],
        'sphere_radius': sphere_r,
        'rotation_deg': disc_rotation
    }

    return {
        'pin': pin,
        'peg': peg,
        'ndisc': ndisc,
        'sdisc': sdisc
    }


def compute_all_strings(config):
    """Compute full geometry for all strings from config.

    This is the main parametric computation function.
    """
    sp = compute_sp(config)
    num_strings = len(config['strings'])
    finger_gap = config['geometry']['finger_gap']
    string_tilt = config['geometry']['string_tilt_deg']
    semitone_ratio = config['physics']['semitone_ratio']
    peg_offset = config['geometry']['peg_offset']
    hardware = config['hardware']

    # Compute *b positions along sp
    b_positions = compute_b_positions(sp, num_strings, finger_gap)

    # Compute full geometry for each string
    strings_computed = []
    for i, string_input in enumerate(config['strings']):
        b = b_positions[i]
        vib_length = string_input['vib_length']

        # Compute string geometry (*b, *f, *n, *s)
        geom = compute_string_geometry(b, vib_length, string_tilt, semitone_ratio)

        # Compute hardware positions
        hw_pos = compute_hardware_positions(geom, i + 1, peg_offset, hardware)

        strings_computed.append({
            'num': i + 1,
            'note': string_input['note'],
            'diameter': string_input['diameter'],
            'b': {'x': b[0], 'y': b[1]},
            'pin': hw_pos['pin'],
            'peg': hw_pos['peg'],
            'ndisc': hw_pos['ndisc'],
            'sdisc': hw_pos['sdisc'],
            # Keep computed values for reference
            '_geom': geom
        })

    return {
        'sp': sp,
        'strings': strings_computed
    }


# =============================================================================
# CadQuery 3D Model Builders
# =============================================================================

def get_disc_spec_for_string(string_num):
    """Get the disc spec dict for a string number.

    Maps string number to bass/mid/treble/high register.
    Optimized boundaries based on collision analysis.
    """
    if string_num <= 9:
        return DISC_SPECS['bass']
    elif string_num <= 28:
        return DISC_SPECS['mid']
    elif string_num <= 38:
        return DISC_SPECS['treble']
    else:
        return DISC_SPECS['high']


def make_positioned_disc(cx, cy, rotation_deg, string_num):
    """Create an elliptical disc assembly positioned for a string.

    The disc is positioned so the prong SCOOP is at z=0 (string level).

    Args:
        cx, cy: Position of disc center on string (at *n or *s position)
        rotation_deg: Rotation angle around Z axis
        string_num: String number (1-47) to determine disc size

    Returns:
        CadQuery solid of positioned disc (no engagement pin for rendering)
    """
    spec = get_disc_spec_for_string(string_num)

    # Create disc assembly (without engagement pin for clean rendering)
    disc = make_disc_assembly(spec, include_engagement_pin=False)

    # Position so scoop groove is at z=0 (where string sits)
    z_offset = get_scoop_z_offset(spec)

    # Apply rotation around Z axis (disc orientation)
    disc = disc.rotate((0, 0, 0), (0, 0, 1), rotation_deg)

    # Position at string location
    disc = disc.translate((cx, cy, z_offset))

    return disc


def make_positioned_dual_disc(cx, cy, rotation_deg, prong_reach, string_num):
    """Create a dual-position disc positioned at the N-S midpoint.

    This disc handles both natural and sharp positions with a single disc.
    The disc is centered between N and S, with prongs reaching to each.
    Scoop grooves are positioned at z=0.

    Args:
        cx, cy: Position of disc center (midpoint between N and S)
        rotation_deg: Rotation angle around Z axis
        prong_reach: Distance from disc center to prong center (N-S distance / 2)
        string_num: String number to determine disc specs

    Returns:
        CadQuery solid of positioned dual-position disc
    """
    spec = DUAL_POSITION_SPECS['high']

    # Create dual-position disc (without engagement pin for clean rendering)
    disc = make_dual_position_disc(prong_reach, spec, include_engagement_pin=False)

    # Position so scoop groove is at z=0
    z_offset = get_dual_scoop_z_offset(spec)

    # Apply rotation around Z axis
    disc = disc.rotate((0, 0, 0), (0, 0, 1), rotation_deg)

    # Position at N-S midpoint
    disc = disc.translate((cx, cy, z_offset))

    return disc


def make_disc_solid_legacy(cx, cy, cz, major_radius, minor_radius, prong_diameter,
                           thickness, rotation_deg, prong_length=None, string_diameter=1.0):
    """Legacy: Create a 3D monolithic disc with integrated prongs.

    DEPRECATED: Use make_positioned_disc() for new elliptical disc design.

    Disc body centered at (cx, cy, cz).
    Prongs extend from disc front face (+Z direction) toward string at z=0.
    Prong groove is sized proportionally to prong diameter.

    Args:
        cx, cy, cz: Disc center position (cz is negative, behind string)
        prong_length: Length of prongs (should equal sphere_radius to reach z=0)
    """
    # Prong inset from disc center
    inset = major_radius - minor_radius

    # Stadium-shaped disc body centered at origin
    slot_length = 2 * inset + 2 * minor_radius
    disc = (cq.Workplane("XY")
            .slot2D(slot_length, 2 * minor_radius, angle=0)
            .extrude(thickness)
            .translate((0, 0, -thickness/2)))

    prong_r = prong_diameter / 2

    # Use provided prong_length or default
    if prong_length is None:
        prong_length = thickness * 2.0

    # Groove sized to prong - max 30% of prong radius
    groove_depth = min(string_diameter * 0.5, prong_r * 0.3)
    groove_width = min(string_diameter * 1.0, prong_r * 0.8)
    groove_r = max(prong_r - groove_depth, prong_r * 0.5)

    # Prongs extend from disc front face (+Z) toward string
    # Groove near the tip where string engages
    groove_start = prong_length - groove_width - prong_r

    for x_offset in [inset, -inset]:
        # Section 1: base to groove
        if groove_start > 0:
            base = (cq.Workplane("XY")
                    .circle(prong_r)
                    .extrude(groove_start)
                    .translate((x_offset, 0, thickness/2)))
            disc = disc.union(base)

        # Section 2: groove (reduced diameter)
        groove = (cq.Workplane("XY")
                  .circle(groove_r)
                  .extrude(groove_width)
                  .translate((x_offset, 0, thickness/2 + max(0, groove_start))))
        disc = disc.union(groove)

        # Section 3: tip above groove
        tip_length = prong_length - max(0, groove_start) - groove_width
        if tip_length > 0:
            tip = (cq.Workplane("XY")
                   .circle(prong_r)
                   .extrude(tip_length)
                   .translate((x_offset, 0, thickness/2 + max(0, groove_start) + groove_width)))
            disc = disc.union(tip)

    disc = disc.rotate((0, 0, 0), (0, 0, 1), rotation_deg)
    disc = disc.translate((cx, cy, cz))
    return disc


def make_tuning_peg(cx, cy, cz, peg_diameter, string_diameter, height):
    """Create a tuning peg with a string hole."""
    peg_r = peg_diameter / 2
    hole_d = string_diameter * 1.5 + 0.5
    hole_r = hole_d / 2

    peg = (cq.Workplane("XY")
           .circle(peg_r)
           .extrude(height)
           .translate((0, 0, -height/2)))

    hole_z = height * 0.2
    hole_length = peg_diameter * 1.5
    hole_cyl = (cq.Workplane("XZ")
                .circle(hole_r)
                .extrude(hole_length)
                .translate((0, -hole_length/2, hole_z)))

    peg = peg.cut(hole_cyl)
    peg = peg.translate((cx, cy, cz))
    return peg


def make_flat_pin(cx, cy, cz, pin_diameter, string_diameter, height):
    """Create a flat pin with a groove for the string."""
    pin_r = pin_diameter / 2
    groove_depth = string_diameter * 0.6
    groove_width = string_diameter * 1.2
    groove_r = pin_r - groove_depth
    section_height = (height - groove_width) / 2

    bottom = (cq.Workplane("XY")
              .circle(pin_r)
              .extrude(section_height))
    groove = (cq.Workplane("XY")
              .workplane(offset=section_height)
              .circle(groove_r)
              .extrude(groove_width))
    top = (cq.Workplane("XY")
           .workplane(offset=section_height + groove_width)
           .circle(pin_r)
           .extrude(section_height))

    pin = bottom.union(groove).union(top)
    pin = pin.translate((cx, cy, cz - height/2))
    return pin


def make_base_2d(config, sp):
    """Create 2D outlines of the pedal base in XY plane.

    Returns list of (wire, color, stroke_width) tuples.
    """
    base_cfg = config['hardware']['base']

    width = base_cfg['width_mm']   # Pedal arrangement (along Z, not visible in XY)
    depth = base_cfg['depth_mm']   # Front-to-back (along X)
    height = base_cfg['height_mm'] # Vertical (along Y)

    # Position: pillar touches sp start, base extends in +X toward player
    sp_start_x = sp['start'][0]
    c1b_y = sp['c1b'][1]

    # Base box - floor level below c1b
    base_y = c1b_y - 150
    x_left = sp_start_x - base_cfg['pillar']['diameter_mm']  # Pillar side
    x_right = x_left + depth  # Player side
    y_bottom = base_y
    y_top = base_y + height

    wires = []

    # Base outline
    base_wire = (cq.Workplane("XY")
                 .moveTo(x_left, y_bottom)
                 .lineTo(x_right, y_bottom)
                 .lineTo(x_right, y_top)
                 .lineTo(x_left, y_top)
                 .close())
    wires.append((base_wire, '#8B4513', 0.8))

    # Pedal - comes through base, extends in +X
    pedal_cfg = base_cfg['pedals']
    pedal_l = pedal_cfg['length_mm']
    pedal_t = pedal_cfg['thickness_mm']

    pivot_x = x_right - 20
    pivot_y = y_bottom + height * 0.6

    # Pedal in natural position (slight upward angle)
    pedal_angle = 10
    pedal_end_x = pivot_x + pedal_l * math.cos(math.radians(pedal_angle))
    pedal_end_y = pivot_y + pedal_l * math.sin(math.radians(pedal_angle))

    pedal_wire = (cq.Workplane("XY")
                  .moveTo(pivot_x, pivot_y - pedal_t/2)
                  .lineTo(pedal_end_x, pedal_end_y - pedal_t/2)
                  .lineTo(pedal_end_x, pedal_end_y + pedal_t/2)
                  .lineTo(pivot_x, pivot_y + pedal_t/2)
                  .close())
    wires.append((pedal_wire, '#2F4F4F', 0.5))

    # Slot in base
    slot_wire = (cq.Workplane("XY")
                 .moveTo(x_right, pivot_y - 6)
                 .lineTo(x_right + 5, pivot_y - 6)
                 .lineTo(x_right + 5, pivot_y + 6)
                 .lineTo(x_right, pivot_y + 6)
                 .close())
    wires.append((slot_wire, '#2F4F4F', 0.3))

    # Three notches for pedal positions
    notch_x = pivot_x + pedal_l * 0.7
    for i, offset in enumerate([-15, 0, 15]):
        notch_y = pivot_y + offset
        notch_wire = (cq.Workplane("XY")
                      .moveTo(notch_x, notch_y - 2)
                      .lineTo(notch_x + 3, notch_y)
                      .lineTo(notch_x, notch_y + 2)
                      .close())
        wires.append((notch_wire, '#2F4F4F', 0.3))

    # Pillar - touches sp start, extends up to c1s (sharp disc of lowest string)
    c1_vib_length = config['strings'][0]['vib_length']
    semitone_ratio = config['physics']['semitone_ratio']
    c1s_y = c1b_y + c1_vib_length / (semitone_ratio ** 2)

    pillar_d = base_cfg['pillar']['diameter_mm']
    pillar_x_right = sp_start_x  # Touch soundboard
    pillar_x_left = pillar_x_right - pillar_d
    pillar_y_bottom = y_top
    pillar_y_top = c1s_y + 20  # Extend slightly above c1s for neck

    pillar_wire = (cq.Workplane("XY")
                   .moveTo(pillar_x_left, pillar_y_bottom)
                   .lineTo(pillar_x_right, pillar_y_bottom)
                   .lineTo(pillar_x_right, pillar_y_top)
                   .lineTo(pillar_x_left, pillar_y_top)
                   .close())
    wires.append((pillar_wire, '#8B4513', 0.8))

    return wires


def make_wire_segment(x1, y1, z1, x2, y2, z2, diameter):
    """Create a cylindrical wire segment (for strings)."""
    dx, dy, dz = x2 - x1, y2 - y1, z2 - z1
    length = math.sqrt(dx*dx + dy*dy + dz*dz)

    if length < 0.001:
        return None

    r = diameter / 2
    wire = (cq.Workplane("XY")
            .circle(r)
            .extrude(length))

    ux, uy, uz = dx/length, dy/length, dz/length

    if abs(uz) > 0.9999:
        if uz < 0:
            wire = wire.rotate((0, 0, 0), (1, 0, 0), 180)
    else:
        angle = math.degrees(math.acos(uz))
        axis = (-uy, ux, 0)
        axis_len = math.sqrt(axis[0]**2 + axis[1]**2)
        if axis_len > 0.0001:
            wire = wire.rotate((0, 0, 0), axis, angle)

    wire = wire.translate((x1, y1, z1))
    return wire


# =============================================================================
# Rendering
# =============================================================================

def note_for_string(num):
    """Return note letter for string number (1-47)."""
    return ['C', 'D', 'E', 'F', 'G', 'A', 'B'][(num - 1) % 7]


def color_for_note(note):
    """Return color for note: C=red, F=blue, others=gray."""
    if note == 'C':
        return '#cc0000'
    elif note == 'F':
        return '#0000cc'
    return '#666666'


def build_groups(computed, config=None):
    """Build geometry groups for cq_plate() from computed positions."""
    groups = []
    sp = computed['sp']
    strings = computed['strings']

    # Soundboard path (sp)
    sp_wire = (cq.Workplane("XY")
               .moveTo(sp['start'][0], sp['start'][1])
               .lineTo(sp['end'][0], sp['end'][1]))
    groups.append((sp_wire, '#8B4513', 1.0))

    # Base assembly (2D outlines)
    if config and 'hardware' in config and 'base' in config['hardware']:
        base_wires = make_base_2d(config, sp)
        groups.extend(base_wires)

    # Collect 3D solid models
    peg_solids = []
    pin_solids = []
    ndisc_solids = []  # Separate N discs for strings 1-41
    sdisc_solids = []  # Separate S discs for strings 1-41
    dual_disc_solids = []  # Combined N+S discs for strings 42-47

    peg_height = 30
    pin_height = 15

    for s in strings:
        string_d = s['diameter']
        string_num = s['num']

        # Tuning pegs
        peg = s.get('peg')
        if peg and 'center_x' in peg:
            peg_solid = make_tuning_peg(
                peg['center_x'], peg['center_y'], 0,
                peg['diameter'], string_d, peg_height
            )
            peg_solids.append(peg_solid)

        # Flat pins
        pin = s.get('pin')
        if pin and 'center_x' in pin:
            pin_solid = make_flat_pin(
                pin['center_x'], pin['center_y'], 0,
                pin['diameter'], string_d, pin_height
            )
            pin_solids.append(pin_solid)

        # Discs: use dual-position for strings 42-47, separate for 1-41
        ndisc = s.get('ndisc')
        sdisc = s.get('sdisc')

        if string_num >= 39:
            # Dual-position disc for high register (strings 39-47)
            if ndisc and sdisc:
                # Calculate midpoint between N and S positions
                mid_x = (ndisc['x'] + sdisc['x']) / 2
                mid_y = (ndisc['y'] + sdisc['y']) / 2

                # Calculate prong reach from N-S distance
                ns_distance = math.sqrt(
                    (ndisc['x'] - sdisc['x'])**2 +
                    (ndisc['y'] - sdisc['y'])**2
                )
                prong_reach = ns_distance / 2

                dual_disc = make_positioned_dual_disc(
                    mid_x, mid_y,
                    rotation_deg=45.0,  # Show in natural position
                    prong_reach=prong_reach,
                    string_num=string_num
                )
                dual_disc_solids.append(dual_disc)
        else:
            # Separate N and S discs for bass/mid/treble (strings 1-38)
            if ndisc:
                ndisc_solid = make_positioned_disc(
                    ndisc['x'], ndisc['y'],
                    ndisc.get('rotation_deg', 45.0),  # Natural position
                    string_num
                )
                ndisc_solids.append(ndisc_solid)

            if sdisc:
                sdisc_solid = make_positioned_disc(
                    sdisc['x'], sdisc['y'],
                    sdisc.get('rotation_deg', 90.0),  # Sharp position
                    string_num
                )
                sdisc_solids.append(sdisc_solid)

    # Union solids by type
    if peg_solids:
        all_pegs = peg_solids[0]
        for p in peg_solids[1:]:
            all_pegs = all_pegs.union(p)
        groups.append((all_pegs, '#009900', 0.3))

    if pin_solids:
        all_pins = pin_solids[0]
        for p in pin_solids[1:]:
            all_pins = all_pins.union(p)
        groups.append((all_pins, '#990099', 0.3))

    if ndisc_solids:
        all_ndiscs = ndisc_solids[0]
        for d in ndisc_solids[1:]:
            all_ndiscs = all_ndiscs.union(d)
        groups.append((all_ndiscs, '#ff6600', 0.2))

    if sdisc_solids:
        all_sdiscs = sdisc_solids[0]
        for d in sdisc_solids[1:]:
            all_sdiscs = all_sdiscs.union(d)
        groups.append((all_sdiscs, '#0066ff', 0.2))

    if dual_disc_solids:
        all_dual_discs = dual_disc_solids[0]
        for d in dual_disc_solids[1:]:
            all_dual_discs = all_dual_discs.union(d)
        groups.append((all_dual_discs, '#00aa88', 0.3))  # Teal for dual-position

    # Build strings as cylindrical wires
    for s in strings:
        b = s['b']
        pin = s.get('pin')
        peg = s.get('peg')

        if not pin or not peg:
            continue

        note = note_for_string(s['num'])
        color = color_for_note(note)
        diameter = s['diameter']

        # String segment 1: soundboard to pin
        seg1 = make_wire_segment(
            b['x'], b['y'], 0,
            pin['x'], pin['y'], 0,
            diameter
        )

        # String segment 2: pin to peg
        peg_tx = peg.get('tangent_x', peg['x'])
        peg_ty = peg.get('tangent_y', peg['y'])
        seg2 = make_wire_segment(
            pin['x'], pin['y'], 0,
            peg_tx, peg_ty, 0,
            diameter
        )

        if seg1:
            groups.append((seg1, color, 0.3))
        if seg2:
            groups.append((seg2, color, 0.3))

    return groups


# =============================================================================
# Main
# =============================================================================

def load_config():
    """Load parametric config from erand.json."""
    with open(Path(__file__).parent / "erand.json") as f:
        return json.load(f)


def main():
    import sys

    config = load_config()
    computed = compute_all_strings(config)
    groups = build_groups(computed, config)

    # Check for --neck flag to zoom into neck area
    if '--neck' in sys.argv:
        svg = plate(groups, UL='XZ', LL='XY', UR='ISO', LR='YZ',
                    grid='light', units='mm',
                    zoom='LL',
                    zoom_region=(200, 1600, 700, 1900))
        output = Path(__file__).parent / "erand_neck.svg"
    else:
        svg = plate(groups, UL='XZ', LL='XY', UR='ISO', LR='YZ',
                    grid='light', units='mm')
        output = Path(__file__).parent / "erand.svg"

    output.write_text(svg)

    strings = computed['strings']
    c_count = sum(1 for s in strings if note_for_string(s['num']) == 'C')
    f_count = sum(1 for s in strings if note_for_string(s['num']) == 'F')

    print(f"Generated {output}")
    print(f"  Strings: {len(strings)}")
    print(f"  C strings (red): {c_count}")
    print(f"  F strings (blue): {f_count}")
    print(f"  Other (gray): {len(strings) - c_count - f_count}")

    dias = [s['diameter'] for s in strings]
    print(f"  Diameter range: {min(dias):.2f}mm - {max(dias):.2f}mm")


if __name__ == '__main__':
    main()
