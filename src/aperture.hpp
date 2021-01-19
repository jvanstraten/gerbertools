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
 * Defines objects representing apertures.
 */

// This was a triumph
// I'm making a note here:
//
//      HUGE SUCCESS
//
// It's hard to overstate my satisfaction

#pragma once

#include <memory>
#include "coord.hpp"
#include "plot.hpp"

/**
 * Namespace for objects representing apertures.
 */
namespace aperture {

/**
 * Base class for aperture objects.
 */
class Base {
protected:

    /**
     * Plot object representing the shape of the aperture.
     */
    plot::Ref plot;

public:
    virtual ~Base() = default;

    /**
     * Returns the plot that represents the shape of the aperture.
     */
    const plot::Plot &get_plot() const;

    /**
     * Returns whether this is a simple circle aperture, suitable for
     * interpolation. If this is indeed a simple circle and diameter is
     * non-null, the circle diameter is written to the pointer in addition.
     */
    virtual bool is_simple_circle(coord::CInt *diameter) const;

};

/**
 * Reference to an aperture object.
 */
using Ref = std::shared_ptr<Base>;

/**
 * Represents a custom aperture object, either built using an aperture macro or
 * a block aperture.
 */
class Custom : public Base {
public:

    /**
     * Constructs the custom aperture from the given plot.
     */
    explicit Custom(const plot::Ref &data);

};

/**
 * Base class for standard apertures.
 */
class Standard : public Base {
protected:

    /**
     * Diameter of the optional hole common to all standard apertures.
     */
    coord::CInt hole_diameter;

    /**
     * Returns the path for the hole. That is, a negatively wound circle, or
     * no path at all if the diameter is zero.
     */
    coord::Paths get_hole(const coord::Format &fmt) const;

};

/**
 * Represents a standard circle aperture.
 */
class Circle : public Standard {
private:

    /**
     * Diameter of the circle.
     */
    coord::CInt diameter;

public:

    /**
     * Constructs the circle aperture from the given parameters, reported as a
     * vector of strings. The first string is ignored; it is assumed to be "C"
     * to select this type of aperture.
     */
    explicit Circle(const std::vector<std::string> &csep, const coord::Format &fmt);

    /**
     * Returns whether this is a simple circle aperture, suitable for
     * interpolation. If this is indeed a simple circle and diameter is
     * non-null, the circle diameter is written to the pointer in addition.
     */
    bool is_simple_circle(coord::CInt *diameter) const override;

};

/**
 * Represents a standard rectangle aperture.
 */
class Rectangle : public Standard {
private:

    /**
     * Horizontal size of the aperture.
     */
    coord::CInt x_size;

    /**
     * Vertical size of the aperture.
     */
    coord::CInt y_size;

public:

    /**
     * Constructs the rectangle aperture from the given parameters, reported as
     * a vector of strings. The first string is ignored; it is assumed to be "R"
     * to select this type of aperture.
     */
    explicit Rectangle(const std::vector<std::string> &csep, const coord::Format &fmt);

};

/**
 * Represents a standard obround aperture.
 */
class Obround : public Standard {
private:

    /**
     * Horizontal size of the aperture.
     */
    coord::CInt x_size;

    /**
     * Vertical size of the aperture.
     */
    coord::CInt y_size;

public:

    /**
     * Constructs the obround aperture from the given parameters, reported as
     * a vector of strings. The first string is ignored; it is assumed to be "O"
     * to select this type of aperture.
     */
    explicit Obround(const std::vector<std::string> &csep, const coord::Format &fmt);

};

/**
 * Represents a standard regular polygon aperture.
 */
class Polygon : public Standard {
private:

    /**
     * Outer diameter of the regular polygon.
     */
    coord::CInt diameter;

    /**
     * Number of vertices. Should be at least 3.
     */
    size_t n_vertices;

    /**
     * Rotation of the polygon in radians counterclockwise from the positive X
     * axis.
     */
    double rotation;

public:

    /**
     * Constructs the polygon aperture from the given parameters, reported as
     * a vector of strings. The first string is ignored; it is assumed to be "P"
     * to select this type of aperture.
     */
    explicit Polygon(const std::vector<std::string> &csep, const coord::Format &fmt);

};

} // namespace aperture
