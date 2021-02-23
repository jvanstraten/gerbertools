# GerberTools

Work in progress.

## Goals

 - `x` Pick a name. I hate naming things. [note: yes, I have no imagination.]

 - `x` Interpret an [RS-274X file](https://www.ucamco.com/files/downloads/file_en/416/the-gerber-file-format-specification-revision-2020-09-update_en.pdf)
   to an independent polygon (using NonZero/OddEven winding rules for holes)
   for each non-connected pour. This is hard, the Gerber file format is
   stupidly esoteric for what it has to do. I suspect this is why there isn't
   much out there that can read Gerber, or at least I couldn't easily find
   anything that does it in the way I want it to. Note that while the result
   is vector, arcs are linearly approximated with some tolerance parameter.

 - `x` As above for NC drill files. At least, the superset of XNC and what I've
   seen Altium output.

 - `x` Ability to export to SVG for manual verification. This is a bit redundent
   as there are plenty of tools out there that can render Gerbers (although
   not usually to other vector formats), but very important for developing
   this actual library.

 - `x` Ability to export to OBJ files or something similar for nice 3D models of
   circuit boards that don't care about poly count; I've yet to find a tool
   that does this in a way that includes the copper traces, soldermask, and
   silkscreen.

 - `x` Ability to perform basic unconnected net, short circuit, and net clearance
   DRC. Clearance will be done by "growing" all polygons by half the clearance
   and then doing a short circuit check, which means that the clearance value
   is global; you can't set it per net.

 - ` ` (Command line) user interface and C++ library.

 - `x` Make a basic Python API for the above, wrap the thing in a Python wheel...

 - ` ` ...and publish on PyPI.

## Why?

 - Because I want nice PCB renders in Blender.

 - Because I have a project in mind that basically requires me to procedurally
   generate PCBs, and I want electrical DRC to check myself. I don't care so
   much about mechanical DRC, since at least the hobby-oriented fabs perform
   such a DRC before you can even make the order anyway.

## Copyright stuff

The polygon geometry processing is done using the
[Clipper](http://www.angusj.com/delphi/clipper.php) library, lazily included
as source files in this repository. It uses the Boost software license. The
same thing goes for [Earcut](https://github.com/mapbox/earcut.hpp/), which uses
the ICS license. The relevant license text is at the top of the copied files.

Stuff I add will be licensed under MIT, which I'm more familiar with. Either
way, it'll be permissively licensed.
