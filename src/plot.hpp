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

#pragma once

#include <memory>
#include "coord.hpp"
#include "clipper.hpp"

namespace plot {

/**
 * Represents a single path, either open or closed (based on context).
 */
using Path = ClipperLib::Path;

/**
 * Represents multiple paths.
 */
using Paths = ClipperLib::Paths;

/**
 * Polygon fill rule.
 */
using FillRule = ClipperLib::PolyFillType;

/**
 * Represents a vector image, tracking which parts were explicitly cleared as
 * well as which parts were made dark.
 */
class Plot {
private:

    /**
     * Accumulator for incoming paths. It is assumed that all paths have
     * positive winding or cancel out previously drawn paths and that the
     * nonzero/positive fill rule shall be applicable unless otherwise
     * specified. When different fill rules are needed, commit_paths() must be
     * called priot to adding the new paths to the accumulator, and then
     * commit_paths() must be called again with the desired fill rule.
     */
    mutable Paths accum_paths;

    /**
     * Whether the paths in the accumulator are to be interpreted as making the
     * plot dark (true) or clear (false).
     */
    mutable bool accum_polarity;

    /**
     * Set of committed paths that make the plot dark. Doesn't intersect clear.
     */
    mutable Paths dark;

    /**
     * Set of committed paths that make the plot clear. Doesn't intersect dark.
     */
    mutable Paths clear;

    /**
     * Whether dark and clear have been simplified yet.
     */
    mutable bool simplified;

    /**
     * Commits paths in the accumulator to dark/clear using the given fill type.
     * No-op if the accumulator is empty.
     */
    void commit_paths(FillRule fill_rule = FillRule::pftNonZero) const;

    /**
     * Simplifies the dark/clear paths. No-op if they have already been
     * simplified.
     */
    void simplify() const;

public:

    /**
     * Constructs a plot, optionally with initial conditions for the dark and
     * clear surfaces.
     */
    explicit Plot(const Paths &dark = {}, const Paths &clear = {});

    /**
     * Adds paths to the plot with the given polarity. It is assumed that the
     * paths are all positively wound, and that nonzero/positive winding rules
     * apply.
     */
    void draw_paths(const Paths &ps, bool polarity = true);

    /**
     * Advanced method for adding paths, allowing coordinate transformations and
     * the fill rule to be specified. The transformation order is translate,
     * rotate, mirror/scale.
     */
    void draw_paths(
        const Paths &ps,
        bool polarity,
        coord::CInt translate_x,
        coord::CInt translate_y = 0,
        bool mirror_x = false,
        bool mirror_y = false,
        double rotate = 0.0,
        double scale = 1.0,
        bool special_fill_type = false,
        FillRule fill_rule = FillRule::pftNonZero
    );

    /**
     * Add an entire subplot to this plot with the given transformation. The
     * transformation order is translate, rotate, mirror/scale.
     */
    void draw_plot(
        const Plot &plt,
        bool polarity = true,
        coord::CInt translate_x = 0,
        coord::CInt translate_y = 0,
        bool mirror_x = false,
        bool mirror_y = false,
        double rotate = 0.0,
        double scale = 1.0
    );

    /**
     * Returns the surface that was made dark as a simplified
     * nonzero/odd-even-filled polygon.
     */
    const Paths &get_dark() const;

    /**
     * Returns the surface that was explicitly made clear as a simplified
     * nonzero/odd-even-filled polygon.
     */
    const Paths &get_clear() const;

};

/**
 * Reference to a plot.
 */
using Ref = std::shared_ptr<Plot>;

/**
 * Renders an open path to a polygon by applying thickness.
 */
Paths render_path(const Path &path, double thickness, const coord::Format &fmt, bool square=false);

} // namespace plot
