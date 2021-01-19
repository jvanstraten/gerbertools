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

#pragma once

#include <memory>
#include "coord.hpp"
#include "clipper.hpp"

/**
 * Contains shorthands for path/polygon operations.
 */
namespace path {

/**
 * Renders an open path to a polygon by applying thickness.
 */
coord::Paths render(const coord::Path &path, double thickness, const coord::Format &fmt, bool square=false);

/**
 * Append paths odd-even style.
 */
void append(coord::Paths &dest, const coord::Paths &src);

/**
 * Compute the union of two sets of paths.
 */
coord::Paths add(const coord::Paths &lhs, const coord::Paths &rhs);

/**
 * Compute the difference between two sets of paths.
 */
coord::Paths subtract(const coord::Paths &lhs, const coord::Paths &rhs);

/**
 * Compute the intersection between two sets of paths.
 */
coord::Paths intersect(const coord::Paths &lhs, const coord::Paths &rhs);

/**
 * Offsets the paths by the given amount.
 */
coord::Paths offset(coord::Paths src, double amount);

} // namespace path
