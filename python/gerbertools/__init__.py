from gerbertools._gerbertools import *
import os
import re
import sys
import argparse
import textwrap

__all__ = ['get_version', 'read', 'CircuitBoard', 'Netlist', 'Shape', 'color', 'main']

def read(prefix, outline_hint=None):
    """Tries to be smart about file detection to load a PCB. If necessary,
    outline can be used to specify the Gerber file extension used for the
    outline; if not specified, GKO or the lowest GM* mechanical layer will be
    used."""
    if os.path.isdir(prefix):
        dirname = prefix
        prefix = ''
    else:
        dirname = os.path.dirname(prefix)
        prefix = os.path.basename(prefix)

    top_copper = ''
    mid_copper = {}
    bot_copper = ''
    top_mask = ''
    bot_mask = ''
    top_silk = ''
    bot_silk = ''
    drill = ''
    drill_nonplated = ''
    outline = {}

    if outline_hint is None:
        outline_hint = 'gko'
    else:
        outline_hint = outline_hint.lower()

    for name in os.listdir(dirname):
        if not name.startswith(prefix):
            continue
        fname = os.path.join(dirname, name)
        if not os.path.isfile(fname):
            continue
        ext = name.split('.')[-1].lower()
        if ext == 'gtl':
            top_copper = fname
        elif ext == 'gbl':
            bot_copper = fname
        elif re.fullmatch(r'g[0-9]+', ext):
            mid_copper[int(ext[1:])] = fname
        elif ext == 'gts':
            top_mask = fname
        elif ext == 'gbs':
            bot_mask = fname
        elif ext == 'gto':
            top_silk = fname
        elif ext == 'gbo':
            bot_silk = fname
        elif ext in ('txt', 'drl'):
            if 'npth' in name.lower():
                drill_nonplated = fname
            else:
                drill = fname
        elif ext == outline_hint:
            outline[ext] = fname
        elif re.fullmatch(r'gm[0-9]+', ext):
            outline[ext] = fname

    if not outline:
        raise ValueError('Could not find any file that looks like an outline!')
    if outline_hint in outline:
        outline = outline[outline_hint]
    else:
        outline = outline[min(outline)]

    substrate_thickness = 1.5 / (1 + len(mid_copper))

    pcb = CircuitBoard('', outline, drill, drill_nonplated);
    if bot_mask:
        pcb.add_mask_layer(bot_mask, bot_silk);
    if bot_copper:
        pcb.add_copper_layer(bot_copper)
    pcb.add_substrate_layer(substrate_thickness)
    for _, fname in reversed(sorted(mid_copper.items())):
        pcb.add_copper_layer(fname)
        pcb.add_substrate_layer(substrate_thickness)
    if top_copper:
        pcb.add_copper_layer(top_copper)
    if top_mask:
        pcb.add_mask_layer(top_mask, top_silk);
    pcb.add_surface_finish()

    return pcb

def main(args=sys.argv[1:]):
    """Main function for the CLI."""

    parser = argparse.ArgumentParser(
        description=textwrap.dedent("""
        GerberTools {version} by Jeroen van Straten
        See https://github.com/jvanstraten/gerbertools

        Tool for generating SVGs and 3D models of circuit boards by means of
        Gerber and NC drill files.

        The first argument is the path prefix or directory in which the board
        files live. If no additional positional arguments are specified,
        GerberTools tries to be smart about finding and loading the files based on
        the usual file extensions (GTL for top copper, TXT for drill, etc.).
        The outline is taken from either GKO or the lowest-indexed GM* file.
        *Usually* this is good enough. If not, --stackup may be used to specify
        the board stackup explicitly. Its value consists of a semicolon-separated
        list of arguments.

        The first of these specifies the suffix for the board outline (such
        that prefix concatenated with suffix gives the complete filename).
        Optionally, a colon may be used to specify two filename suffixes, in
        case milling and outline are on two different layers.

        The second argument specifies the drill file suffix(es); again,
        optionally, a colon may be used to specify two filenames, in which case
        the first specifies plated holes, and the second specifies non-plated.
        If only a single file is specified, plating is derived based on
        "comments" in the drill file.

        The remaining arguments define the board stackup from bottom to top.
        - "m:<mask>[:<silk>]": adds a mask layer with optional silkscreen.
        - "c:<copper>[:<thickness>]": adds a copper layer. Thickness is specified
        in millimeters, and defaults to 0.0348 (= 1 oz).
        - "s:<thickness>": adds a substrate layer with the specified thickness in
        millimeters.

        For example, a typical 2-layer stackup would be:

            --stackup "GKO;TXT;m:GBS:GBO;c:GBL;s:1.5;c:GTL;m:GTS:GTO"

        Without anything else specified, GerberTools will just try to load the
        circuit board and quit. To do something useful, use --bounds, --size,
        --front, --back, --obj, and/or --drc.

        For --drc, the input file format for the netlist must be a CSV file with
        columns X, Y, layer index, and netname. The coordinates are in millimeters.
        Layer indices start at 0, starting from the bottom copper layer. The DRC
        checks that the given locations actually contain copper, and that said
        copper connects to exactly those other points in the file that have the
        same netname. A global minimum clearance between nets can also be specified
        (--clearance), in which case copper for different nets must be at least
        that far apart. Finally, a minimum annular ring size for vias/pads can be
        specified using --annular.
        """).format(version=get_version()),
        formatter_class=argparse.RawTextHelpFormatter
    )

    def parse_stackup(stackup):
        try:
            stackup = stackup.split(';')

            outline = stackup[0].split(':')
            if len(outline) > 2:
                raise ValueError('at most two outline layers may be specified')
            elif len(outline) == 1:
                outline.append('')

            drill = stackup[1].split(':') if len(stackup) > 0 else []
            if len(drill) > 2:
                raise ValueError('at most two drill files may be specified')
            while len(drill) < 2:
                drill.append('')

            constructor_args = (outline[0], drill[0], drill[1], outline[1])

            layers = []
            for index, cmd in enumerate(stackup[2:]):
                args = cmd.split(':')
                if len(args) < 2:
                    raise ValueError('missing m/c/s command for layer {}'.format(index))
                if args[0] == 'm':
                    if len(args) > 3:
                        raise ValueError(
                            'm(ask) command for layer {} must have 1 or 2 filenames'
                            .format(index)
                        )
                    layers.append((CircuitBoard.add_mask_layer, args[1:]))
                elif args[0] == 'c':
                    if len(args) > 3:
                        raise ValueError(
                            'c(opper) command for layer {} must have 1 filename and '
                            'an optional thickness'
                            .format(index)
                        )
                    layer = args[1]
                    if len(args) > 2:
                        try:
                            thickness = float(args[2])
                        except ValueError:
                            raise ValueError(
                                'c(opper) command for layer {} has invalid thickness'
                                .format(index)
                            )
                    else:
                        thickness = 0.0348
                    layers.append((CircuitBoard.add_copper_layer, (layer, thickness)))
                elif args[0] == 's':
                    if len(args) > 2:
                        raise ValueError(
                            's(ubstrate) command for layer {} must have a single'
                            'thickness argument'
                            .format(index)
                        )
                    try:
                        thickness = float(args[1])
                    except ValueError:
                        raise ValueError(
                            's(ubstrate) command for layer {} has invalid thickness'
                            .format(index)
                        )
                    layers.append((CircuitBoard.add_substrate_layer, (thickness,)))
                else:
                    raise ValueError('unknown command {} for layer {}'.format(cmd[0], index))

            return constructor_args, layers
        except Exception as e:
            print(e)
            raise

    parser.add_argument(
        'prefix',
        metavar='prefix',
        help='path prefix or directory containing the board files'
    )
    parser.add_argument(
        '--outline',
        metavar='GKO',
        help='overrides the file extension used for the outline when --stackup is not used'
    )
    parser.add_argument(
        '--stackup',
        metavar='...',
        type=parse_stackup,
        help='specifies an explicit board stackup; see above for syntax'
    )
    parser.add_argument(
        '--bounds',
        action='store_true',
        help='prints the boundbox of the PCB to stdout in '
        'min_x,min_y,max_x,max_y format (millimeters)'
    )
    parser.add_argument(
        '--size',
        action='store_true',
        help='prints the size of the PCB to stdout in width,height format (millimeters)'
    )
    parser.add_argument(
        '--front',
        metavar='out.svg',
        help='renders the front of the PCB to the given SVG file'
    )
    parser.add_argument(
        '--back',
        metavar='out.svg',
        help='renders the back of the PCB to the given SVG file'
    )
    parser.add_argument(
        '--svg-scale',
        metavar='out.svg',
        type=float,
        default=1.0,
        help='SVG scale: the number of SVG millimeters corresponding with a PCB millimeter'
    )
    parser.add_argument(
        '--obj',
        metavar='out.obj',
        help='renders the PCB to the given Wavefront OBJ file'
    )
    parser.add_argument(
        '--drc',
        metavar='in.csv',
        help='performs DRC for, as described above'
    )
    parser.add_argument(
        '--clearance',
        metavar='X',
        type=float,
        default=0,
        help='specifies the global copper-to-copper clearance for --drc'
    )
    parser.add_argument(
        '--annular',
        metavar='X',
        type=float,
        default=0,
        help='specifies the minimum annular ring for --drc'
    )

    args = parser.parse_args(args)

    # Gather results.
    code = 0
    output = []
    def println(msg):
        output.append(msg)
        print(msg)

    # Load the PCB.
    print('loading PCB...', file=sys.stderr)
    if args.stackup is None:
        pcb = read(args.prefix, args.outline)
    else:
        pcb = CircuitBoard(args.prefix, *args.stackup[0])
        for f, a in args.stackup[1]:
            f(pcb, *a)
        pcb.add_surface_finish()

    # Get bounds/size.
    if args.bounds:
        print('computing bounds...', file=sys.stderr)
        println('{:.3f},{:.3f},{:.3f},{:.3f}'.format(*pcb.get_bounds()))
    if args.size:
        print('computing size...', file=sys.stderr)
        b = pcb.get_bounds()
        println('{:.3f},{:.3f}'.format(b[2] - b[0], b[3] - b[1]))

    # Output SVG(s).
    if args.front is not None:
        print('writing front SVG...', file=sys.stderr)
        pcb.write_svg(args.front, False, args.svg_scale)
    if args.back is not None:
        print('writing back SVG...', file=sys.stderr)
        pcb.write_svg(args.back, True, args.svg_scale)

    # Output object file.
    if args.obj is not None:
        print('writing OBJ file...', file=sys.stderr)
        pcb.write_obj(args.obj)

    # Perform DRC.
    if args.drc is not None:
        print('reading netlist...', file=sys.stderr)
        nets = []
        with open(args.drc, 'r') as f:
            for index, line in enumerate(f.readlines()):
                line = line.strip()
                if not line:
                    continue
                try:
                    a = line.split(',', maxsplit=3)
                    if len(a) != 4:
                        raise ValueError('missing argument')
                    x, y, layer, net = a
                    x = float(x.strip())
                    y = float(y.strip())
                    layer = int(layer.strip())
                    nets.append(((x, y), layer, net.strip()))
                except ValueError as e:
                    raise ValueError('on line {} of {}: {}'.format(index+1, args.drc, e))
        print('building physical netlist...', file=sys.stderr)
        physical_netlist = pcb.build_netlist(nets, args.clearance, args.annular)
        print('running DRC...', file=sys.stderr)
        violations = physical_netlist.drc()
        for violation in violations:
            println(violation)
            code = 1

    print('done!', file=sys.stderr)

    return output, code
