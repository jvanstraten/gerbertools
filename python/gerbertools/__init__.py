from gerbertools._gerbertools import *
import os
import re

def read(prefix, outline_hint='GKO'):
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
    outline = {}

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
        elif ext == 'txt':
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

    pcb = CircuitBoard('', outline, drill);
    if bot_mask:
        pcb.add_mask_layer(bot_mask, bot_silk);
    if bot_copper:
        pcb.add_copper_layer(bot_copper)
    pcb.add_substrate_layer(substrate_thickness)
    for _, fname in mid_copper.items():
        pcb.add_copper_layer(fname)
        pcb.add_substrate_layer(substrate_thickness)
    if top_copper:
        pcb.add_copper_layer(top_copper)
    if top_mask:
        pcb.add_mask_layer(top_mask, top_silk);
    pcb.add_surface_finish()

    return pcb
