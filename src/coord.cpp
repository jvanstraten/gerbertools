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
    add_trailing_zeros(false),
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
 * Configures whether trailing zeros may be omitted.
 */
void Format::configure_trailing_zeros(bool add_trailing_zeros) {
    try_to_reconfigure();
    this->add_trailing_zeros = add_trailing_zeros;
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
 * CInt representation. Falls back to parse_float() when a period is found in
 * the string, however.
 */
CInt Format::parse_fixed(const std::string &s) const {
    try_to_use();
    if (s.find('.') != std::string::npos) {
        return parse_float(s);
    }
    size_t add_zeros = 10 - n_dec;
    if (factor == 25.4) {
        add_zeros = 9 - n_dec;
    } else if (factor != 1.0) {
        throw std::runtime_error("unknown conversion factor");
    }
    if (add_trailing_zeros) {
        size_t digits = s.size();
        if (s.at(0) == '-' || s.at(0) == '+') {
            digits--;
        }
        if (digits < n_int + n_dec) {
            add_zeros += n_int + n_dec - digits;
        }
    }
    CInt val = std::stoll(s + std::string(add_zeros, '0'));
    if (factor == 25.4) {
        val *= 254;
    }
    return val;
}

/**
 * Parses a floating-point coordinate and converts it to the internal 64-bit
 * CInt representation.
 */
CInt Format::parse_float(const std::string &s) const {
    return to_fixed(std::stod(s));
}

/**
 * Converts a previously parsed floating-point coordinate to the internal
 * 64-bit CInt representation.
 */
CInt Format::to_fixed(double d) const {
    try_to_use();
    return std::round(d * factor * 1e10);
}

/**
 * Returns the maximum deviation in the internal 64-bit CInt representation.
 */
CInt Format::get_max_deviation() const {
    return from_mm(max_deviation);
}

/**
 * Returns the miter limit in the internal 64-bit CInt representation.
 */
CInt Format::get_miter_limit() const {
    return from_mm(miter_limit);
}

/**
 * Returns a ClipperOffset object with appropriate configuration.
 */
ClipperLib::ClipperOffset Format::build_clipper_offset() const {
    return ClipperLib::ClipperOffset(get_miter_limit(), get_max_deviation());
}

/**
 * Converts millimeters to the internal 64-bit CInt representation.
 */
CInt Format::from_mm(double i) {
    return (CInt)std::round(i * 1e10);
}

/**
 * Converts the internal 64-bit CInt representation to millimeters.
 */
double Format::to_mm(CInt i) {
    return i / 1e10;
}

} // namespace coord
