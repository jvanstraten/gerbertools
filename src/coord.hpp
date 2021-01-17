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

#pragma once

#include "clipper.hpp"

namespace coord {

/**
 * Internal 64-bit fixed-point coordinate representation. Conversion from Gerber
 * units and to millimeters is handled by the CoordFormat class.
 */
using CInt = ClipperLib::cInt;

/**
 * 2D point, consisting of two CInts.
 */
using CPt = ClipperLib::IntPoint;

/**
 * Coordinate format handling. This class converts between Gerber fixed-point
 * and floating point format coordinates, an internal high-accuracy integer
 * representation used for polygon operations, and millimeters for the output.
 * It also stores maximum deviation and miter limit for the polygon operations.
 */
class Format {
private:

    /**
     * Whether the fixed point format has been configured yet.
     */
    bool fmt_configured;

    /**
     * Number of integer positions.
     */
    int n_int;

    /**
     * Number of decimal positions.
     */
    int n_dec;

    /**
     * Whether the unit (inch or millimeters) has been configured yet.
     */
    bool unit_configured;

    /**
     * Conversion factor from Gerber unit to millimeters; 25.4 for inch, 1.0 for
     * millimeters.
     */
    double factor;

    /**
     * Multiplier from Gerber fixed format to our high-accuracy 64-bit internal
     * fixed-point representation.
     */
    static const CInt PRECISION_MULT = 0x1000000ll;

    /**
     * Whether any coordinates have been converted yet. An exception is thrown
     * if the format is reconfigured after this, as the internal representation
     * changes with the Gerber format.
     */
    mutable bool used;

    /**
     * Maximum arc deviation in millimeters.
     */
    double max_deviation;

    /**
     * Miter limit in millimeters.
     */
    double miter_limit;

    /**
     * Throws an exception if the coordinate format has been used to convert
     * coordinates already.
     */
    void try_to_reconfigure() const;

    /**
     * Throws an exception if the coordinate format has not been fully
     * configured yet.
     */
    void try_to_use() const;

public:

    /**
     * Constructs a new coordinate format object with the given maximum arc
     * deviation and miter limit (both in millimeters).
     */
    explicit Format(double max_deviation=0.005, double miter_limit=1.0);

    /**
     * Configure the coordinate format to the given number of integer and
     * decimal digits.
     */
    void configure_format(int n_int, int n_dec);

    /**
     * Use Freedum units.
     */
    void configure_inch();

    /**
     * Use sane units.
     */
    void configure_mm();

    /**
     * Parses a fixed-point coordinate and converts it to the internal 64-bit
     * CInt representation.
     */
    CInt parse_fixed(const std::string &s) const;

    /**
     * Parses a floating-point coordinate and converts it to the internal 64-bit
     * CInt representation.
     */
    CInt parse_float(const std::string &s) const;

    /**
     * Converts a previously parsed floating-point coordinate to the internal
     * 64-bit CInt representation.
     */
    CInt to_fixed(double d) const;

    /**
     * Returns the maximum deviation in the internal 64-bit CInt representation.
     */
    CInt get_max_deviation() const;

    /**
     * Returns the miter limit in the internal 64-bit CInt representation.
     */
    CInt get_miter_limit() const;

    /**
     * Returns a ClipperOffset object with appropriate configuration.
     */
    ClipperLib::ClipperOffset build_clipper_offset() const;

    /**
     * Converts the internal 64-bit CInt representation to millimeters.
     */
    double to_mm(CInt i) const;

};

} // namespace coord
