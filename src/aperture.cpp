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

#define _USE_MATH_DEFINES
#include <cmath>
#include "gerbertools/aperture.hpp"
#include "gerbertools/path.hpp"

namespace gerbertools {
namespace aperture {

/**
 * Returns the plot that represents the shape of the aperture.
 */
const plot::Plot &Base::get_plot() const {
    return *plot;
}

/**
 * Returns whether this is a simple circle aperture, suitable for
 * interpolation. If this is indeed a simple circle and diameter is
 * non-null, the circle diameter is written to the pointer in addition.
 */
bool Base::is_simple_circle(coord::CInt *diameter) const {
    (void)diameter;
    return false;
}

/**
 * Constructs the custom aperture from the given plot.
 */
Custom::Custom(const plot::Ref &data) {
    plot = data;
}

/**
 * Returns the path for the hole. That is, a negatively wound circle, or
 * no path at all if the diameter is zero.
 */
coord::Paths Standard::get_hole(const coord::Format &fmt) const {
    if (hole_diameter <= 0.0) {
        return {};
    }
    auto paths = path::render({{0, 0}}, hole_diameter, fmt);
    ClipperLib::ReversePaths(paths);
    return paths;
}

/**
 * Constructs the circle aperture from the given parameters, reported as a
 * vector of strings. The first string is ignored; it is assumed to be "C"
 * to select this type of aperture.
 */
Circle::Circle(const std::vector<std::string> &csep, const coord::Format &fmt) {

    // Parse the command.
    if (csep.size() < 2 || csep.size() > 3) {
        throw std::runtime_error("invalid circle aperture");
    }
    diameter = fmt.parse_float(csep.at(1));
    hole_diameter = (csep.size() > 2) ? fmt.parse_float(csep.at(2)) : 0;

    // Construct the plot.
    auto paths = path::render({{0, 0}}, diameter, fmt);
    auto hole = get_hole(fmt);
    paths.insert(paths.end(), hole.begin(), hole.end());
    plot = std::make_shared<plot::Plot>(paths);

}

/**
 * Returns whether this is a simple circle aperture, suitable for
 * interpolation. If this is indeed a simple circle and diameter is
 * non-null, the circle diameter is written to the pointer in addition.
 */
bool Circle::is_simple_circle(coord::CInt *diameter) const {
    if (hole_diameter > 0.0) return false;
    if (diameter) *diameter = this->diameter;
    return true;
}

/**
 * Constructs the rectangle aperture from the given parameters, reported as
 * a vector of strings. The first string is ignored; it is assumed to be "R"
 * to select this type of aperture.
 */
Rectangle::Rectangle(const std::vector<std::string> &csep, const coord::Format &fmt) {

    // Parse the command.
    if (csep.size() < 3 || csep.size() > 4) {
        throw std::runtime_error("invalid rectangle aperture");
    }
    x_size = std::abs(fmt.parse_float(csep.at(1)));
    y_size = std::abs(fmt.parse_float(csep.at(2)));
    hole_diameter = (csep.size() > 3) ? fmt.parse_float(csep.at(3)) : 0;

    // Construct the plot.
    coord::CInt x = x_size / 2;
    coord::CInt y = y_size / 2;
    coord::Paths paths{{
         {x, y},
         {x, -y},
         {-x, -y},
         {-x, y}
    }};
    auto hole = get_hole(fmt);
    paths.insert(paths.end(), hole.begin(), hole.end());
    plot = std::make_shared<plot::Plot>(paths);

}

/**
 * Constructs the obround aperture from the given parameters, reported as
 * a vector of strings. The first string is ignored; it is assumed to be "O"
 * to select this type of aperture.
 */
Obround::Obround(const std::vector<std::string> &csep, const coord::Format &fmt) {

    // Parse the command.
    if (csep.size() < 3 || csep.size() > 4) {
        throw std::runtime_error("invalid obround aperture");
    }
    x_size = std::abs(fmt.parse_float(csep.at(1)));
    y_size = std::abs(fmt.parse_float(csep.at(2)));
    hole_diameter = (csep.size() > 3) ? fmt.parse_float(csep.at(3)) : 0;

    // Construct the plot.
    coord::CInt x = x_size / 2;
    coord::CInt y = y_size / 2;
    coord::CInt r = std::min(x, y);
    x -= r;
    y -= r;
    auto paths = path::render({{-x, -y}, {x, y}}, r * 2.0, fmt);
    auto hole = get_hole(fmt);
    paths.insert(paths.end(), hole.begin(), hole.end());
    plot = std::make_shared<plot::Plot>(paths);

}

/**
 * Constructs the polygon aperture from the given parameters, reported as
 * a vector of strings. The first string is ignored; it is assumed to be "P"
 * to select this type of aperture.
 */
Polygon::Polygon(const std::vector<std::string> &csep, const coord::Format &fmt) {

    // Parse the command.
    if (csep.size() < 3 || csep.size() > 5) {
        throw std::runtime_error("invalid polygon aperture");
    }
    diameter = fmt.parse_float(csep.at(1));
    n_vertices = std::stoul(csep.at(2));
    if (n_vertices < 3) {
        throw std::runtime_error("invalid polygon aperture");
    }
    rotation = (csep.size() > 3) ? (std::stod(csep.at(3)) / 180.0 * M_PI) : 0.0;
    hole_diameter = (csep.size() > 4) ? fmt.parse_float(csep.at(4)) : 0;

    // Construct the plot.
    coord::Paths paths = {{}};
    for (size_t i = 0; i < n_vertices; i++) {
        double a = ((double)i / (double)n_vertices) * 2.0 * M_PI;
        paths.back().push_back({
            (coord::CInt)std::round(diameter * 0.5 * std::cos(a + rotation)),
            (coord::CInt)std::round(diameter * 0.5 * std::sin(a + rotation))
        });
    }
    auto hole = get_hole(fmt);
    paths.insert(paths.end(), hole.begin(), hole.end());
    plot = std::make_shared<plot::Plot>(paths);

}

} // namespace aperture
} // namespace gerbertools
