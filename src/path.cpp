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
 * Defines shorthands for path/polygon operations.
 */

#include "gerbertools/clipper.hpp"
#include "gerbertools/path.hpp"

namespace gerbertools {
namespace path {

/**
 * Renders an open path to a polygon by applying thickness.
 */
coord::Paths render(const coord::Paths &paths, double thickness, bool square, ClipperLib::ClipperOffset &&co) {
    co.AddPaths(
        paths,
        square ? ClipperLib::jtMiter : ClipperLib::jtRound,
        square ? ClipperLib::etOpenButt : ClipperLib::etOpenRound
    );
    coord::Paths out;
    co.Execute(out, thickness * 0.5);
    return out;
}

/**
 * Append paths odd-even style.
 */
void append(coord::Paths &dest, const coord::Paths &src) {
    if (src.empty()) {
        return;
    }
    if (dest.empty()) {
        dest = src;
        return;
    }
    dest.insert(dest.end(), src.begin(), src.end());
    ClipperLib::SimplifyPolygons(dest);
}

/**
 * Perform an operation between two sets of paths.
 */
static coord::Paths path_op(const coord::Paths &lhs, const coord::Paths &rhs, ClipperLib::ClipType op) {
    ClipperLib::Clipper c;
    c.AddPaths(lhs, ClipperLib::ptSubject, true);
    c.AddPaths(rhs, ClipperLib::ptClip, true);
    coord::Paths result;
    c.Execute(op, result);
    return result;
}

/**
 * Compute the union of two sets of paths.
 */
coord::Paths add(const coord::Paths &lhs, const coord::Paths &rhs) {
    return path_op(lhs, rhs, ClipperLib::ctUnion);
}

/**
 * Compute the difference between two sets of paths.
 */
coord::Paths subtract(const coord::Paths &lhs, const coord::Paths &rhs) {
    return path_op(lhs, rhs, ClipperLib::ctDifference);
}

/**
 * Compute the intersection between two sets of paths.
 */
coord::Paths intersect(const coord::Paths &lhs, const coord::Paths &rhs) {
    return path_op(lhs, rhs, ClipperLib::ctIntersection);
}

/**
 * Offsets the paths by the given amount.
 */
coord::Paths offset(const coord::Paths &src, double amount, bool square, ClipperLib::ClipperOffset &&co) {
    co.AddPaths(src, square ? ClipperLib::jtMiter : ClipperLib::jtRound, ClipperLib::etClosedLine);
    coord::Paths result;
    if (amount < 0) {
        co.Execute(result, -amount);
        return subtract(src, result);
    } else {
        co.Execute(result, amount);
        return add(src, result);
    }
}

} // namespace path
} // namespace gerbertools
