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
 * Defines the Plot class, representing a monochrome image by means of polygons.
 * Also includes some utility functions.
 */

#include <cmath>
#include "plot.hpp"

namespace plot {

/**
 * Commits paths in the accumulator to dark/clear using the given fill type.
 * No-op if the accumulator is empty.
 */
void Plot::commit_paths(FillRule fill_rule) const {
    if (accum_paths.empty()) return;
    ClipperLib::SimplifyPolygons(accum_paths, fill_rule);
    ClipperLib::Clipper cld, clc;
    cld.AddPaths(dark, ClipperLib::ptSubject, true);
    clc.AddPaths(clear, ClipperLib::ptSubject, true);
    cld.AddPaths(accum_paths, ClipperLib::ptClip, true);
    clc.AddPaths(accum_paths, ClipperLib::ptClip, true);
    cld.Execute(
        accum_polarity ? ClipperLib::ctUnion : ClipperLib::ctDifference,
        dark, FillRule::pftNonZero, fill_rule
    );
    clc.Execute(
        accum_polarity ? ClipperLib::ctDifference : ClipperLib::ctUnion,
        clear, FillRule::pftNonZero, fill_rule
    );
    simplified = false;
    accum_paths.clear();
}

/**
 * Simplifies the dark/clear paths. No-op if they have already been
 * simplified.
 */
void Plot::simplify() const {
    if (simplified) return;
    ClipperLib::SimplifyPolygons(dark, FillRule::pftNonZero);
    ClipperLib::SimplifyPolygons(clear, FillRule::pftNonZero);
    simplified = true;
}

/**
 * Constructs a plot, optionally with initial conditions for the dark and
 * clear surfaces.
 */
Plot::Plot(const Paths &dark, const Paths &clear) :
    accum_polarity(true),
    dark(dark),
    clear(clear),
    simplified(false)
{
}

/**
 * Adds paths to the plot with the given polarity. It is assumed that the
 * paths are all positively wound, and that nonzero/positive winding rules
 * apply.
 */
void Plot::draw_paths(const Paths &ps, bool polarity) {
    if (ps.empty()) return;

    // If the polarity is not the same as the accumulator, we have to commit
    // the accumulator first.
    if (polarity != accum_polarity) commit_paths();
    accum_polarity = polarity;

    // Simply add to the accumulator.
    accum_paths.insert(accum_paths.end(), ps.begin(), ps.end());
}

/**
 * Advanced method for adding paths, allowing coordinate transformations and
 * the fill rule to be specified. The transformation order is translate,
 * rotate, mirror/scale.
 */
void Plot::draw_paths(
    const Paths &ps, bool polarity,
    coord::CInt translate_x, coord::CInt translate_y,
    bool mirror_x, bool mirror_y,
    double rotate, double scale,
    bool special_fill_type, FillRule fill_rule
) {
    if (ps.empty()) return;

    // If we need to apply a special fill rule, commit first, so previously
    // drawn paths won't interfere.
    if (special_fill_type) commit_paths();

    // Draw the paths normally.
    draw_paths(ps, polarity);

    // Compute transformation matrix.
    double ixx = mirror_x ? -scale : scale;
    double iyy = mirror_y ? -scale : scale;
    double si = std::sin(rotate);
    double co = std::cos(rotate);
    double xx = ixx * co;
    double xy = ixx * si;
    double yx = iyy * -si;
    double yy = iyy * co;

    // Transform the paths we just added.
    for (auto it = accum_paths.end() - ps.size(); it != accum_paths.end(); ++it) {
        for (auto &c : *it) {
            double cx = c.X * xx + c.Y * yx;
            double cy = c.X * xy + c.Y * yy;
            c.X = std::round(cx) + translate_x;
            c.Y = std::round(cy) + translate_y;
        }
    }

    // Maintain winding direction in spite of mirroring.
    if (mirror_x != mirror_y) {
        for (auto it = accum_paths.end() - ps.size(); it != accum_paths.end(); ++it) {
            ClipperLib::ReversePath(*it);
        }
    }

    // If we need to apply a special fill rule, commit immediately with said
    // fill rules.
    if (special_fill_type) commit_paths(fill_rule);

}

/**
 * Add an entire subplot to this plot with the given transformation. The
 * transformation order is translate, rotate, mirror/scale.
 */
void Plot::draw_plot(
    const Plot &plt, bool polarity,
    coord::CInt translate_x, coord::CInt translate_y,
    bool mirror_x, bool mirror_y,
    double rotate, double scale
) {
    draw_paths(
        plt.get_dark(), polarity,
        translate_x, translate_y,
        mirror_x, mirror_y,
        rotate, scale
    );
    draw_paths(
        plt.get_clear(), !polarity,
        translate_x, translate_y,
        mirror_x, mirror_y,
        rotate, scale
    );
}

/**
 * Returns the surface that was made dark as a simplified
 * nonzero/odd-even-filled polygon.
 */
const Paths &Plot::get_dark() const {
    commit_paths();
    simplify();
    return dark;
}

/**
 * Returns the surface that was explicitly made clear as a simplified
 * nonzero/odd-even-filled polygon.
 */
const Paths &Plot::get_clear() const {
    commit_paths();
    simplify();
    return clear;
}

/**
 * Renders an open path to a polygon by applying thickness.
 */
Paths render_path(const Path &path, double thickness, const coord::Format &fmt, bool square) {
    auto co = fmt.build_clipper_offset();
    co.AddPath(
        path,
        square ? ClipperLib::jtMiter : ClipperLib::jtRound,
        square ? ClipperLib::etOpenButt : ClipperLib::etOpenRound
    );
    Paths paths;
    co.Execute(paths, thickness * 0.5);
    return paths;
}

} // namespace plot
