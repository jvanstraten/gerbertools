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
 * Types and conversions for Gerber, NC drill, and internal coordinate formats.
 */

#include <cmath>
#include "coord.hpp"

namespace coord {

/**
 * Throws an exception if the coordinate format has been used to convert
 * coordinates already.
 */
void Format::try_to_reconfigure() const {
    if (used) {
        throw std::runtime_error(
            "cannot reconfigure coordinate format after "
            "coordinates have already been interpreted"
        );
    }
}

/**
 * Throws an exception if the coordinate format has not been fully
 * configured yet.
 */
void Format::try_to_use() const {
    if (!fmt_configured) {
        throw std::runtime_error(
            "cannot convert coordinates before coordinate "
            "format is configured"
        );
    }
    if (!unit_configured) {
        throw std::runtime_error(
            "cannot convert coordinates before unit is configured"
        );
    }
    used = true;
}

/**
 * Constructs a new coordinate format object with the given maximum arc
 * deviation and miter limit (both in millimeters).
 */
Format::Format(double max_deviation, double miter_limit) :
    fmt_configured(false),
    unit_configured(false),
    used(false),
    miter_limit(miter_limit),
    max_deviation(max_deviation)
{
}

/**
 * Configure the coordinate format to the given number of integer and
 * decimal digits.
 */
void Format::configure_format(int n_int, int n_dec) {
    try_to_reconfigure();
    fmt_configured = true;
    this->n_int = n_int;
    this->n_dec = n_dec;
}

/**
 * Use Freedum units.
 */
void Format::configure_inch() {
    try_to_reconfigure();
    unit_configured = true;
    factor = 25.4;
}

/**
 * Use sane units.
 */
void Format::configure_mm() {
    try_to_reconfigure();
    unit_configured = true;
    factor = 1.0;
}

/**
 * Parses a fixed-point coordinate and converts it to the internal 64-bit
 * CInt representation.
 */
CInt Format::parse_fixed(const std::string &s) const {
    try_to_use();
    return std::stoll(s) * PRECISION_MULT;
}

/**
 * Parses a floating-point coordinate and converts it to the internal 64-bit
 * CInt representation.
 */
CInt Format::parse_float(const std::string &s) const {
    try_to_use();
    return std::round(std::stod(s) * std::pow(10.0, n_dec) * PRECISION_MULT);
}

/**
 * Converts a previously parsed floating-point coordinate to the internal
 * 64-bit CInt representation.
 */
CInt Format::to_fixed(double d) const {
    try_to_use();
    return std::round(d * std::pow(10.0, n_dec) * PRECISION_MULT);
}

/**
 * Returns the maximum deviation in the internal 64-bit CInt representation.
 */
CInt Format::get_max_deviation() const {
    try_to_use();
    return to_fixed(max_deviation / factor);
}

/**
 * Returns the miter limit in the internal 64-bit CInt representation.
 */
CInt Format::get_miter_limit() const {
    try_to_use();
    return to_fixed(miter_limit / factor);
}

/**
 * Returns a ClipperOffset object with appropriate configuration.
 */
ClipperLib::ClipperOffset Format::build_clipper_offset() const {
    return ClipperLib::ClipperOffset(get_miter_limit(), get_max_deviation());
}

/**
 * Converts the internal 64-bit CInt representation to millimeters.
 */
double Format::to_mm(CInt i) const {
    try_to_use();
    return i / (std::pow(10.0, n_dec) * PRECISION_MULT) * factor;
}

} // namespace coord
