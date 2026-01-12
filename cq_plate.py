"""
CadQuery plate() extension for multi-view mechanical drawings.

This module vendors the CadQuery SVG export code to ensure stability
across CadQuery version changes.

Usage:
    from cq_plate import plate
    import cadquery as cq

    box = cq.Workplane("XY").box(10, 20, 5)
    svg = plate(box, UL='XY', LL='XZ', UR='ISO', LR='YZ', grid='light', units='mm')
"""

import re
import io
import xml.etree.ElementTree as ET

import cadquery as cq
from cadquery.occ_impl.shapes import Shape, Compound, TOLERANCE
from cadquery.occ_impl.geom import BoundBox

# OpenCascade imports (stable low-level API)
from OCP.gp import gp_Ax2, gp_Pnt, gp_Dir
from OCP.BRepLib import BRepLib
from OCP.HLRBRep import HLRBRep_Algo, HLRBRep_HLRToShape
from OCP.HLRAlgo import HLRAlgo_Projector
from OCP.GCPnts import GCPnts_QuasiUniformDeflection


# ============================================================================
# Vendored CadQuery SVG export code (from cadquery.occ_impl.exporters.svg)
# ============================================================================

DISCRETIZATION_TOLERANCE = 1e-3

SVG_TEMPLATE = """<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg
   xmlns:svg="http://www.w3.org/2000/svg"
   xmlns="http://www.w3.org/2000/svg"
   width="%(width)s"
   height="%(height)s"

>
    <g transform="scale(%(unitScale)s, -%(unitScale)s)   translate(%(xTranslate)s,%(yTranslate)s)" stroke-width="%(strokeWidth)s"  fill="none">
       <!-- hidden lines -->
       <g  stroke="rgb(%(hiddenColor)s)" fill="none" stroke-dasharray="%(strokeWidth)s,%(strokeWidth)s" >
%(hiddenContent)s
       </g>

       <!-- solid lines -->
       <g  stroke="rgb(%(strokeColor)s)" fill="none">
%(visibleContent)s
       </g>
    </g>
    %(axesIndicator)s
</svg>
"""

AXES_TEMPLATE = """<g transform="translate(20,%(textboxY)s)" stroke="rgb(0,0,255)">
        <line x1="30" y1="-30" x2="75" y2="-33" stroke-width="3" stroke="#000000" />
         <text x="80" y="-30" style="stroke:#000000">X </text>

        <line x1="30" y1="-30" x2="30" y2="-75" stroke-width="3" stroke="#000000" />
         <text x="25" y="-85" style="stroke:#000000">Y </text>

        <line x1="30" y1="-30" x2="58" y2="-15" stroke-width="3" stroke="#000000" />
         <text x="65" y="-5" style="stroke:#000000">Z </text>
    </g>"""

PATHTEMPLATE = '\t\t\t<path d="%s" />\n'


class UNITS:
    MM = "mm"
    IN = "in"


def _guessUnitOfMeasure(shape):
    """Guess the unit of measure of a shape."""
    bb = BoundBox._fromTopoDS(shape.wrapped)
    dimList = [bb.xlen, bb.ylen, bb.zlen]
    if max(dimList) > 10:
        return UNITS.MM
    if min(dimList) < 0.1:
        return UNITS.IN
    if sum(dimList) < 10:
        return UNITS.IN
    return UNITS.MM


def _makeSVGedge(e):
    """Creates an SVG edge from an OCCT edge."""
    cs = io.StringIO()
    curve = e._geomAdaptor()
    start = curve.FirstParameter()
    end = curve.LastParameter()

    points = GCPnts_QuasiUniformDeflection(curve, DISCRETIZATION_TOLERANCE, start, end)

    if points.IsDone():
        point_it = (points.Value(i + 1) for i in range(points.NbPoints()))
        p = next(point_it)
        cs.write("M{},{} ".format(p.X(), p.Y()))
        for p in point_it:
            cs.write("L{},{} ".format(p.X(), p.Y()))

    return cs.getvalue()


def _getPaths(visibleShapes, hiddenShapes):
    """Collects the visible and hidden edges from shapes."""
    hiddenPaths = []
    visiblePaths = []

    for s in visibleShapes:
        for e in s.Edges():
            visiblePaths.append(_makeSVGedge(e))

    for s in hiddenShapes:
        for e in s.Edges():
            hiddenPaths.append(_makeSVGedge(e))

    return (hiddenPaths, visiblePaths)


def getSVG(shape, opts=None):
    """
    Export a shape to SVG text.

    Vendored from cadquery.occ_impl.exporters.svg to ensure stability.

    :param shape: A CadQuery shape object to convert to an SVG string.
    :param opts: Dictionary with options:
        width: Width of the resulting image (None to fit based on height).
        height: Height of the resulting image (None to fit based on width).
        marginLeft: Inset margin from the left side of the document.
        marginTop: Inset margin from the top side of the document.
        projectionDir: Direction the camera will view the shape from.
        showAxes: Whether to show the axes indicator.
        strokeWidth: Width of visible edge lines (-1 for auto).
        strokeColor: RGB tuple for visible edges.
        hiddenColor: RGB tuple for hidden edges.
        showHidden: Whether to show hidden lines.
        focus: If specified, creates perspective SVG.
    """
    d = {
        "width": 800,
        "height": 240,
        "marginLeft": 200,
        "marginTop": 20,
        "projectionDir": (-1.75, 1.1, 5),
        "showAxes": True,
        "strokeWidth": -1.0,
        "strokeColor": (0, 0, 0),
        "hiddenColor": (160, 160, 160),
        "showHidden": True,
        "focus": None,
    }

    if opts:
        d.update(opts)

    uom = _guessUnitOfMeasure(shape)

    width = float(d["width"]) if d["width"] is not None else None
    height = float(d["height"]) if d["height"] is not None else None
    marginLeft = float(d["marginLeft"])
    marginTop = float(d["marginTop"])
    projectionDir = tuple(d["projectionDir"])
    showAxes = bool(d["showAxes"])
    strokeWidth = float(d["strokeWidth"])
    strokeColor = tuple(d["strokeColor"])
    hiddenColor = tuple(d["hiddenColor"])
    showHidden = bool(d["showHidden"])
    focus = float(d["focus"]) if d.get("focus") else None

    hlr = HLRBRep_Algo()
    hlr.Add(shape.wrapped)

    coordinate_system = gp_Ax2(gp_Pnt(), gp_Dir(*projectionDir))

    if focus is not None:
        projector = HLRAlgo_Projector(coordinate_system, focus)
    else:
        projector = HLRAlgo_Projector(coordinate_system)

    hlr.Projector(projector)
    hlr.Update()
    hlr.Hide()

    hlr_shapes = HLRBRep_HLRToShape(hlr)

    visible = []
    visible_sharp_edges = hlr_shapes.VCompound()
    if not visible_sharp_edges.IsNull():
        visible.append(visible_sharp_edges)

    visible_smooth_edges = hlr_shapes.Rg1LineVCompound()
    if not visible_smooth_edges.IsNull():
        visible.append(visible_smooth_edges)

    visible_contour_edges = hlr_shapes.OutLineVCompound()
    if not visible_contour_edges.IsNull():
        visible.append(visible_contour_edges)

    hidden = []
    hidden_sharp_edges = hlr_shapes.HCompound()
    if not hidden_sharp_edges.IsNull():
        hidden.append(hidden_sharp_edges)

    hidden_contour_edges = hlr_shapes.OutLineHCompound()
    if not hidden_contour_edges.IsNull():
        hidden.append(hidden_contour_edges)

    # Fix underlying geometry to prevent segfaults
    for el in visible:
        BRepLib.BuildCurves3d_s(el, TOLERANCE)
    for el in hidden:
        BRepLib.BuildCurves3d_s(el, TOLERANCE)

    # Convert to native CQ objects
    visible = list(map(Shape, visible))
    hidden = list(map(Shape, hidden))
    (hiddenPaths, visiblePaths) = _getPaths(visible, hidden)

    # Get bounding box in 2D space
    bb = Compound.makeCompound(hidden + visible).BoundingBox()

    if width is None or height is None:
        if width is None:
            width = (height - (2.0 * marginTop)) * (bb.xlen / bb.ylen) + 2.0 * marginLeft
        else:
            height = (width - 2.0 * marginLeft) * (bb.ylen / bb.xlen) + 2.0 * marginTop
        unitScale = (width - 2.0 * marginLeft) / bb.xlen
    else:
        bb_scale = 0.75
        unitScale = min(width / bb.xlen * bb_scale, height / bb.ylen * bb_scale)

    (xTranslate, yTranslate) = (
        (0 - bb.xmin) + marginLeft / unitScale,
        (0 - bb.ymax) - marginTop / unitScale,
    )

    if strokeWidth == -1.0:
        strokeWidth = 1.0 / unitScale

    hiddenContent = ""
    if showHidden:
        for p in hiddenPaths:
            hiddenContent += PATHTEMPLATE % p

    visibleContent = ""
    for p in visiblePaths:
        visibleContent += PATHTEMPLATE % p

    if showAxes and projectionDir == (-1.75, 1.1, 5):
        axesIndicator = AXES_TEMPLATE % (
            {"unitScale": str(unitScale), "textboxY": str(height - 30), "uom": str(uom)}
        )
    else:
        axesIndicator = ""

    svg = SVG_TEMPLATE % (
        {
            "unitScale": str(unitScale),
            "strokeWidth": str(strokeWidth),
            "strokeColor": ",".join([str(x) for x in strokeColor]),
            "hiddenColor": ",".join([str(x) for x in hiddenColor]),
            "hiddenContent": hiddenContent,
            "visibleContent": visibleContent,
            "xTranslate": str(xTranslate),
            "yTranslate": str(yTranslate),
            "width": str(width),
            "height": str(height),
            "textboxY": str(height - 30),
            "uom": str(uom),
            "axesIndicator": axesIndicator,
        }
    )

    return svg


# ============================================================================
# Plate layout code
# ============================================================================

# Projection directions for each view type
PROJECTIONS = {
    'XY': (0, 0, 1),   # Front view - looking along Z axis
    'XZ': (0, 1, 0),   # Top view - looking down Y axis
    'YZ': (1, 0, 0),   # Side view - looking along X axis
    'ISO': (-1.75, 1.1, 5),  # Isometric view
}

# Grid style colors (green)
GRID_STYLES = {
    'light': {
        'minor': '#b0d8b0',
        'major': '#70b070',
    },
    'medium': {
        'minor': '#90c890',
        'major': '#50a050',
    },
    'heavy': {
        'minor': '#70b070',
        'major': '#309030',
    },
}


def _get_shape(workplane):
    """Extract the underlying shape from a Workplane or Assembly."""
    if hasattr(workplane, 'val'):
        return workplane.val()
    elif hasattr(workplane, 'toCompound'):
        return workplane.toCompound()
    return workplane


def _extract_paths_and_bounds(svg_str):
    """Extract path data and bounding info from SVG output."""
    root = ET.fromstring(svg_str)
    ns = {'svg': 'http://www.w3.org/2000/svg'}

    width = float(root.get('width', 800))
    height = float(root.get('height', 240))

    main_g = root.find('.//svg:g[@transform]', ns)
    if main_g is None:
        main_g = root.find('.//{http://www.w3.org/2000/svg}g[@transform]')

    transform = main_g.get('transform') if main_g is not None else ''

    scale_match = re.search(r'scale\(([-\d.]+),\s*([-\d.]+)\)', transform)
    translate_match = re.search(r'translate\(([-\d.]+),\s*([-\d.]+)\)', transform)

    unit_scale = float(scale_match.group(1)) if scale_match else 1.0
    x_translate = float(translate_match.group(1)) if translate_match else 0.0
    y_translate = float(translate_match.group(2)) if translate_match else 0.0

    hidden_paths = []
    visible_paths = []

    for g in root.iter('{http://www.w3.org/2000/svg}g'):
        stroke = g.get('stroke', '')
        dasharray = g.get('stroke-dasharray', '')

        for path in g.findall('{http://www.w3.org/2000/svg}path'):
            d = path.get('d', '')
            if d:
                if dasharray:
                    hidden_paths.append(d)
                elif 'rgb' in stroke:
                    visible_paths.append(d)

    return {
        'hidden_paths': hidden_paths,
        'visible_paths': visible_paths,
        'unit_scale': unit_scale,
        'x_translate': x_translate,
        'y_translate': y_translate,
        'width': width,
        'height': height,
    }


def _swap_coordinates(path_d):
    """Swap X,Y coordinates in path data for view alignment."""
    def swap_pair(match):
        x, y = match.group(1), match.group(2)
        return f'{y},{x}'
    return re.sub(r'(-?[\d.]+),(-?[\d.]+)', swap_pair, path_d)


def _negate_y(path_d):
    """Negate Y coordinates in path data."""
    def negate(match):
        x = match.group(1)
        y = float(match.group(2))
        return f'{x},{-y}'
    return re.sub(r'(-?[\d.]+),(-?[\d.]+)', negate, path_d)


def _get_path_bounds(paths):
    """Calculate bounding box from path data."""
    all_coords = []
    for path in paths:
        coords = re.findall(r'-?[\d.]+', path)
        for i in range(0, len(coords) - 1, 2):
            x, y = float(coords[i]), float(coords[i+1])
            all_coords.append((x, y))

    if not all_coords:
        return (0, 0, 100, 100)

    xs = [c[0] for c in all_coords]
    ys = [c[1] for c in all_coords]
    return (min(xs), min(ys), max(xs), max(ys))


def _render_grid(width, height, grid_style, units):
    """Render engineering grid as SVG elements."""
    if grid_style not in GRID_STYLES:
        return ''

    colors = GRID_STYLES[grid_style]
    lines = []

    if units == 'mm':
        minor_spacing = 1.0
        major_spacing = 10.0
    else:
        minor_spacing = 2.54
        major_spacing = 25.4

    # Minor grid lines (vertical)
    x = 0
    while x <= width:
        if x % major_spacing != 0:
            lines.append(f'<line x1="{x}" y1="0" x2="{x}" y2="{height}" stroke="{colors["minor"]}" stroke-width="0.0625"/>')
        x += minor_spacing

    # Minor grid lines (horizontal)
    y = 0
    while y <= height:
        if y % major_spacing != 0:
            lines.append(f'<line x1="0" y1="{y}" x2="{width}" y2="{y}" stroke="{colors["minor"]}" stroke-width="0.0625"/>')
        y += minor_spacing

    # Major grid lines (vertical)
    x = 0
    while x <= width:
        lines.append(f'<line x1="{x}" y1="0" x2="{x}" y2="{height}" stroke="{colors["major"]}" stroke-width="0.125"/>')
        x += major_spacing

    # Major grid lines (horizontal)
    y = 0
    while y <= height:
        lines.append(f'<line x1="0" y1="{y}" x2="{width}" y2="{y}" stroke="{colors["major"]}" stroke-width="0.125"/>')
        y += major_spacing

    return '\n    '.join(lines)


def _extract_axes_indicator(svg_str):
    """Extract axes indicator group from SVG."""
    match = re.search(r'(<g transform="translate\(20[^)]*\)"[^>]*>.*?</g>)\s*</svg>', svg_str, re.DOTALL)
    if match:
        return match.group(1)
    return None


def _generate_view_svg(shape, view_type, show_axes=False):
    """Generate SVG for a single view."""
    proj_dir = PROJECTIONS.get(view_type, PROJECTIONS['ISO'])

    opts = {
        'projectionDir': proj_dir,
        'showAxes': show_axes and view_type == 'ISO',
        'showHidden': True,
        'width': None,
        'height': 200,
        'marginLeft': 10,
        'marginTop': 10,
    }

    try:
        svg_str = getSVG(shape, opts)
        return svg_str
    except Exception:
        return '''<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200">
</svg>'''


def _transform_paths(paths, scale, offset_x, offset_y, view_type='XY'):
    """Transform path coordinates with scale and offset."""
    result = []
    for path_d in paths:
        if view_type == 'XZ':
            path_d = _swap_coordinates(path_d)
        elif view_type == 'YZ':
            path_d = _negate_y(path_d)

        def transform_coord(match):
            x = float(match.group(1)) * scale + offset_x
            y = float(match.group(2)) * scale + offset_y
            return f'{x:.4f},{y:.4f}'

        transformed = re.sub(r'(-?[\d.]+),(-?[\d.]+)', transform_coord, path_d)
        result.append(transformed)

    return result


def plate(workplane, UL='XZ', LL='XY', UR='ISO', LR='YZ',
          width=297, height=420, grid=None, units='mm', show_axes=False,
          colors=None):
    """
    Generate a multi-view mechanical drawing plate.

    Layout (3rd angle projection):
        UL: Top (XZ)   | UR: ISO
        ---------------+----------
        LL: Front (XY) | LR: Side (YZ)

    Args:
        workplane: CadQuery Workplane, OR list of (workplane, color, stroke_width) tuples
        UL: View type for upper-left quadrant (default 'XZ' = Top)
        LL: View type for lower-left quadrant (default 'XY' = Front)
        UR: View type for upper-right quadrant (default 'ISO')
        LR: View type for lower-right quadrant (default 'YZ' = Side)
        width: Plate width in mm (default 297 for A3 portrait)
        height: Plate height in mm (default 420 for A3 portrait)
        grid: Engineering grid style (None, 'light', 'medium', 'heavy')
        units: 'mm' (metric) or 'in' (imperial) - affects grid spacing
        show_axes: Show axis indicator on ISO view (default False)

    Returns:
        SVG string with origin at lower-left corner
    """
    # Support list of (geometry, color) or (geometry, color, stroke_width) tuples
    if isinstance(workplane, list):
        geometry_groups = []
        for item in workplane:
            if len(item) == 3:
                geometry_groups.append((item[0], item[1], item[2]))
            else:
                geometry_groups.append((item[0], item[1], 0.5))
    else:
        geometry_groups = [(workplane, '#000000', 0.5)]

    shape = _get_shape(geometry_groups[0][0])

    quadrant_w = width / 2
    quadrant_h = height / 2
    margin = 5

    views = {
        'UL': {'type': UL, 'x': 0, 'y': quadrant_h},
        'LL': {'type': LL, 'x': 0, 'y': 0},
        'UR': {'type': UR, 'x': quadrant_w, 'y': quadrant_h},
        'LR': {'type': LR, 'x': quadrant_w, 'y': 0},
    }

    view_data = {}
    for quad, info in views.items():
        view_data[quad] = {
            'view_type': info['type'],
            'quad_x': info['x'],
            'quad_y': info['y'],
            'groups': [],
            'axes_svg': None
        }

        for geom, color, stroke_width in geometry_groups:
            shape = _get_shape(geom)
            svg_str = _generate_view_svg(shape, info['type'], show_axes=False)
            paths_data = _extract_paths_and_bounds(svg_str)
            view_data[quad]['groups'].append({
                'visible': paths_data['visible_paths'],
                'hidden': paths_data['hidden_paths'],
                'color': color,
                'stroke_width': stroke_width
            })

        if show_axes and info['type'] == 'ISO':
            shape = _get_shape(geometry_groups[0][0])
            svg_str = _generate_view_svg(shape, info['type'], show_axes=True)
            view_data[quad]['axes_svg'] = _extract_axes_indicator(svg_str)

    ortho_views = [q for q in views if views[q]['type'] != 'ISO']

    def correct_paths_for_view(paths, vtype):
        if vtype == 'XZ':
            return [_swap_coordinates(p) for p in paths]
        elif vtype == 'YZ':
            return [_negate_y(p) for p in paths]
        return paths

    max_extent = 0
    for quad in ortho_views:
        data = view_data[quad]
        all_paths = []
        for group in data['groups']:
            all_paths.extend(group['visible'] + group['hidden'])
        if all_paths:
            corrected = correct_paths_for_view(all_paths, data['view_type'])
            bounds = _get_path_bounds(corrected)
            extent = max(bounds[2] - bounds[0], bounds[3] - bounds[1])
            max_extent = max(max_extent, extent)

    available_size = min(quadrant_w, quadrant_h) - 2 * margin
    ortho_scale = available_size / max_extent if max_extent > 0 else 1.0

    svg_parts = []

    svg_parts.append(f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{width}mm" height="{height}mm"
     viewBox="0 0 {width} {height}">
  <!-- Flip Y axis: origin at lower-left -->
  <g transform="translate(0,{height}) scale(1,-1)">
''')

    if grid:
        grid_svg = _render_grid(width, height, grid, units)
        svg_parts.append(f'    <!-- Grid -->\n    {grid_svg}\n')

    svg_parts.append(f'''    <!-- Quadrant dividers -->
    <line x1="{quadrant_w}" y1="0" x2="{quadrant_w}" y2="{height}" stroke="#cccccc" stroke-width="0.5"/>
    <line x1="0" y1="{quadrant_h}" x2="{width}" y2="{quadrant_h}" stroke="#cccccc" stroke-width="0.5"/>
''')

    for quad, data in view_data.items():
        view_type = data['view_type']
        quad_x = data['quad_x']
        quad_y = data['quad_y']

        all_paths = []
        for group in data['groups']:
            all_paths.extend(group['visible'] + group['hidden'])

        if not all_paths:
            continue

        def correct_path(p):
            if view_type == 'XZ':
                return _swap_coordinates(p)
            elif view_type == 'YZ':
                return _negate_y(p)
            return p

        corrected_paths = [correct_path(p) for p in all_paths]
        bounds = _get_path_bounds(corrected_paths)

        path_w = bounds[2] - bounds[0]
        path_h = bounds[3] - bounds[1]

        if view_type == 'ISO':
            scale = min((quadrant_w - 2*margin) / path_w,
                       (quadrant_h - 2*margin) / path_h) if path_w > 0 and path_h > 0 else 1.0
        else:
            scale = ortho_scale

        scaled_w = path_w * scale
        scaled_h = path_h * scale
        offset_x = quad_x + (quadrant_w - scaled_w) / 2 - bounds[0] * scale
        offset_y = quad_y + (quadrant_h - scaled_h) / 2 - bounds[1] * scale

        svg_parts.append(f'    <!-- {quad} view: {view_type} -->')
        svg_parts.append(f'    <g>')

        for group in data['groups']:
            color = group['color']
            stroke_width = group.get('stroke_width', 0.5)
            hidden = _transform_paths(group['hidden'], scale, offset_x, offset_y, view_type=view_type)
            visible = _transform_paths(group['visible'], scale, offset_x, offset_y, view_type=view_type)

            if hidden:
                hidden_sw = stroke_width * 0.6
                svg_parts.append(f'      <g stroke="{color}" stroke-opacity="0.4" stroke-dasharray="2,2" stroke-width="{hidden_sw}" fill="none">')
                for p in hidden:
                    svg_parts.append(f'        <path d="{p}"/>')
                svg_parts.append(f'      </g>')

            if visible:
                svg_parts.append(f'      <g stroke="{color}" stroke-width="{stroke_width}" fill="none">')
                for p in visible:
                    svg_parts.append(f'        <path d="{p}"/>')
                svg_parts.append(f'      </g>')

        if data.get('axes_svg') and view_type == 'ISO':
            axes_x = quad_x + 10
            axes_y = quad_y + 100
            svg_parts.append(f'      <g transform="translate({axes_x},{axes_y}) scale(1,-1)">')
            svg_parts.append(f'        {data["axes_svg"]}')
            svg_parts.append(f'      </g>')

        svg_parts.append(f'    </g>')

    svg_parts.append('  </g>')
    svg_parts.append('</svg>')

    return '\n'.join(svg_parts)


if __name__ == '__main__':
    # Test with L-shape
    L = (cq.Workplane("XY")
         .box(4, 8, 1)
         .faces(">X")
         .workplane()
         .move(0, 2)
         .box(3, 4, 1, centered=False))

    svg = plate(L, UL='XY', LL='XZ', UR='ISO', LR='YZ', grid='light')

    with open('/tmp/test_plate.svg', 'w') as f:
        f.write(svg)

    print("Wrote /tmp/test_plate.svg")
