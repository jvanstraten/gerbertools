/**
 * MIT License
 *
 * Copyright (c) 2021 Jeroen van Straten
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/** \file
 * Defines the top-level class for parsing Gerber files.
 */

#pragma once

#include <string>
#include <list>
#include <map>
#include "coord.hpp"
#include "plot.hpp"
#include "aperture.hpp"
#include "aperture_macro.hpp"

/**
 * Namespace for the Gerber file reader.
 */
namespace gerber {

using Path = plot::Path;
using Paths = plot::Paths;

/**
 * Gerber interpolation mode.
 */
enum class InterpolationMode {

    /**
     * Mode has not been set yet.
     */
    UNDEFINED,

    /**
     * Linear interpolation (G01).
     */
    LINEAR,

    /**
     * Clockwise circular interpolation (G02).
     */
    CIRCULAR_CW,

    /**
     * Counterclockwise circular interpolation (G03).
     */
    CIRCULAR_CCW
};

/**
 * Gerber quadrant mode for circular interpolation.
 */
enum class QuadrantMode {

    /**
     * Mode has not been set yet.
     */
    UNDEFINED,

    /**
     * Insane single-quadrant mode (G74).
     */
    SINGLE,

    /**
     * Sort of reasonable multi-quadrant mode (G75).
     */
    MULTI
};

/**
 * Main class for parsing Gerber files. The input is parsed during construction
 * of the object.
 */
class Gerber {
private:

    /**
     * Mapping from aperture index (10..infinity) to Aperture objects. Each
     * aperture is either one of the StandardAperture specializations or a
     * CustomAperture, representing either an aperture macro or a block aperture
     * depending on how it was built.
     */
    std::map<size_t, aperture::Ref> apertures;

    /**
     * Mapping from name to aperture macro data. When an aperture macro is
     * instantiated, a plot is built on-the-fly to replay the macro instructions
     * with the given parameters.
     */
    std::map<std::string, std::shared_ptr<aperture_macro::ApertureMacro>> aperture_macros;

    /**
     * When non-null, an aperture macro is being constructed. Commands are
     * expected to be in the comma-separated aperture macro command syntax, to
     * be variable assignment statements of the form $<var>=<expr>, or to be
     * aperture macro comments. The attribute termination character (%) always
     * resets this to nullptr.
     */
    std::shared_ptr<aperture_macro::ApertureMacro> am_builder;

    /**
     * Stack of Plot images. The first entry is the actual image; this always
     * exists. An additional plot is pushed onto the stack when a block aperture
     * is started. Graphics objects are always pushed onto the last entry.
     */
    std::list<plot::Ref> plot_stack;

    /**
     * Coordinate format information, including expected digit counts and
     * mm/inch switch. Note that only leading zeros may be stripped; stripping
     * trailing zeros is a deprecated Gerber feature and is not supported.
     */
    coord::Format fmt;

    /**
     * Stores linear vs. circular interpolation mode.
     */
    InterpolationMode imode;

    /**
     * Stores circular interpolation quadrant mode.
     */
    QuadrantMode qmode;

    /**
     * The currently selected aperture, or nullptr if none have been selected
     * yet.
     */
    aperture::Ref aperture;

    /**
     * Current coordinate.
     */
    coord::CPt pos;

    /**
     * Current polarity. True for dark, false for clear. Essentially, when
     * polarity is false, drawn objects remove previously drawn features.
     */
    bool polarity;

    /**
     * Whether apertures should be mirrored on the X axis.
     */
    bool ap_mirror_x;

    /**
     * Whether apertures should be mirrored on the Y axis.
     */
    bool ap_mirror_y;

    /**
     * Current rotation for apertures, in radians counter-clockwise from the
     * positive X axis.
     */
    double ap_rotate;

    /**
     * Current scale factor for apertures.
     */
    double ap_scale;

    /**
     * When set, we're inside a G36/G37 region block. Segments drawn via D01
     * commands are added to region_accum instead of being rendered as lines.
     */
    bool region_mode;

    /**
     * Path accumulator for the region command. When we're inside a G36/G37
     * block, segments drawn via D01 are added to the rather than being
     * rendered. When a G37 or D02 is encountered, the path is immediately
     * committed to the current plot and the accumulator is cleared via
     * commit_region().
     */
    Path region_accum;

    /**
     * Accumulates all linear or circular interpolations drawn onto the topmost
     * plot via D01 commands, be they inside a region or not. After the Gerber
     * has been read, the contents of this can be used to look for loops, which
     * may then be used to reconstruct board outline and milling data when the
     * Gerber in question is the board outline and/or milling layer.
     */
    Paths outline;

    /**
     * Whether the outline has been constructed yet. If false, outline contains
     * the accumulated paths. If true, outline represents the board shape.
     */
    bool outline_constructed;

    /**
     * Render the current aperture to the current plot, taking into
     * consideration all configured aperture transformations.
     */
    void draw_aperture();

    /**
     * Render an interpolation of the current aperture to the current plot, or
     * add the interpolation to the current region if we're inside a G36/G37
     * block.
     */
    void interpolate(coord::CPt dest, coord::CPt center);

    /**
     * Commits the region path accumulated in region_accum via a G36/G37 block
     * to the current plot.
     */
    void commit_region();

    /**
     * Handles a Gerber command. Returns true to continue, false if the command
     * marks the end, or throws a runtime error if the command is not
     * recognized.
     */
    bool command(const std::string &cmd, bool is_attrib);

    /**
     * Handles the end of an attribute.
     */
    void end_attrib();

public:

    /**
     * Loads a gerber file from the given stream. This reads until the end
     * command; the stream may not be EOF if there is more data at the end of
     * the stream.
     */
    explicit Gerber(std::istream &s);

    /**
     * Returns the paths representing the Gerber file.
     */
    const Paths &get_paths() const;

    /**
     * Attempts to interpret the Gerber file data as the board outline and/or
     * milling data, returning polygons that follow the center of closed loops
     * of traces/regions in the file rather than the traces themselves. This is
     * a bit sensitive to round-off error and probably not work right if the
     * file isn't a proper outline; your mileage may vary.
     */
    const Paths &get_outline_paths();

};

} // namespace gerber
