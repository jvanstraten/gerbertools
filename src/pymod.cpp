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
 * PyBind11 module for the project.
 */

#include <pybind11/pybind11.h>
#include "gerbertools/coord.hpp"
#include "gerbertools/path.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

using namespace gerbertools;

class Shape {
private:

    /**
     * Path data.
     */
    coord::Paths p;

    /**
     * Resolution, in steps per floating point unit.
     */
    double resolution;

    /**
     * Miter limit in float units.
     */
    double miter_limit;

    /**
     * Arc tolerance in float units.
     */
    double arc_tolerance;

    size_t get_index(int64_t index) const {
        if (index < p.size() || index >= p.size()) {
            throw std::domain_error("index " + std::to_string(index) + " is out of range");
        }
        if (index < 0) {
            index += p.size();
        }
        return (size_t)index;
    }

    Shape combine(const Shape &rhs, coord::Paths (*op)(const coord::Paths &, const coord::Paths &)) const {
        const Shape *rhs_conv = &rhs;
        auto rhs2 = Shape(resolution);
        if (rhs.resolution != resolution) {
            double factor = resolution / rhs.resolution;
            for (const auto &path : rhs.p) {
                rhs2.p.emplace_back();
                for (const auto &pt : path) {
                    rhs2.p.back().push_back({
                        (coord::CInt)std::round(pt.X * factor),
                        (coord::CInt)std::round(pt.Y * factor)
                    });
                }
            }
            rhs_conv = &rhs2;
        }
        auto out = Shape(resolution);
        out.p = op(p, rhs_conv->p);
        return out;
    }

    ClipperLib::ClipperOffset get_co() const {
        return ClipperLib::ClipperOffset(miter_limit * resolution, arc_tolerance * resolution);
    }

public:
    explicit Shape(
        double resolution = 1e10,
        double miter_limit = 1.0,
        double arc_tolerance = 0.005
    ) :
        resolution(resolution),
        miter_limit(miter_limit),
        arc_tolerance(arc_tolerance)
    { }

    double get_resolution() const {
        return resolution;
    }

    double get_miter_limit() const {
        return miter_limit;
    }

    double set_miter_limit(double x) {
        miter_limit = x;
    }

    double get_arc_tolerance() const {
        return arc_tolerance;
    }

    double set_arc_tolerance(double x) {
        arc_tolerance = x;
    }

    void append(const std::vector<std::tuple<double, double>> &data) {
        p.emplace_back();
        for (const auto &pt : data) {
            p.back().push_back({
                (coord::CInt)std::round(std::get<0>(pt) * resolution),
                (coord::CInt)std::round(std::get<1>(pt) * resolution)
            });
        }
    }

    void append_int(const std::vector<std::tuple<int64_t, int64_t>> &data) {
        p.emplace_back();
        for (const auto &pt : data) {
            p.back().push_back({
                (coord::CInt)std::round(std::get<0>(pt)),
                (coord::CInt)std::round(std::get<1>(pt))
            });
        }
    }

    std::vector<std::tuple<double, double>> get(int64_t index) {
        std::vector<std::tuple<double, double>> out;
        for (const auto &pt : p.at(get_index(index))) {
            out.push_back({
                pt.X / resolution,
                pt.Y / resolution
            });
        }
        return out;
    }

    std::vector<std::tuple<int64_t, int64_t>> get_int(int64_t index) {
        std::vector<std::tuple<int64_t, int64_t>> out;
        for (const auto &pt : p.at(get_index(index))) {
            out.push_back({
                pt.X,
                pt.Y
            });
        }
        return out;
    }

    void del(int64_t index) {
        p.erase(p.begin() + get_index(index));
    }

    int64_t len() const {
        return p.size();
    }

    void simplify() {
        ClipperLib::SimplifyPolygons(p);
    }

    Shape render(double thickness, bool square) const {
        auto out = Shape(resolution);
        out.p = path::render(p, thickness * resolution, square, get_co());
        return out;
    }

    Shape offset(double amount, bool square) const {
        auto out = Shape(resolution);
        out.p = path::offset(p, amount * resolution, square, get_co());
        return out;
    }

    Shape union_(const Shape &rhs) const {
        return combine(rhs, path::add);
    }

    Shape difference(const Shape &rhs) const {
        return combine(rhs, path::subtract);
    }

    Shape intersect(const Shape &rhs) const {
        return combine(rhs, path::intersect);
    }
};

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(_gerbertools, m) {
    py::class_<Shape>(m, "Shape", "A shape described via zero or more odd-even-wound polygons.")
        .def(py::init<double, double, double>(), py::arg("resolution") = 1e10, py::arg("miter_limit") = 1.0, py::arg("arc_tolerance") = 0.005,
             "Creates a new shape. Coordinates are 64-bit integers internally; resolution specifies how many integer steps there should be per floating point unit. "
             "miter_limit and arc_tolerance specify the tolerances to use for render() and offset() in floating point units.")
        .def("get_resolution", &Shape::get_resolution, "Returns the resolution.")
        .def("get_miter_limit", &Shape::get_miter_limit, "Returns the current miter limit.")
        .def("set_miter_limit", &Shape::set_miter_limit, "Reconfigures the miter limit.")
        .def("get_arc_tolerance", &Shape::get_arc_tolerance, "Returns the current arc tolerance.")
        .def("set_arc_tolerance", &Shape::set_arc_tolerance, "Reconfigures the arc tolerance.")
        .def("append", &Shape::append, "Adds a path to the shape.")
        .def("append_int", &Shape::append_int, "Adds a path to the shape using the internal integer format.")
        .def("__getitem__", &Shape::get, "Reads a path from the shape.")
        .def("get_int", &Shape::get_int, "Reads a path from the shape in the internal integer format.")
        .def("__del__", &Shape::del, "Deletes a path from the shape.")
        .def("__len__", &Shape::len, "Returns the number of paths in this shape.")
        .def("simplify", &Shape::simplify, "Simplifies the paths, interpreting the input with odd-even winding.")
        .def("render", &Shape::render, py::arg("amount"), py::arg("square") = false, "Interprets the path components as open polygons and adds thickness to them. Square selects mitered vs rounded corners and butt vs round endpoints.")
        .def("offset", &Shape::offset, py::arg("amount"), py::arg("square") = false, "Offsets polygons by some amount. Square selects mitered vs rounded corners.")
        .def("union", &Shape::union_, "Computes the union between two shapes.")
        .def("__add__", &Shape::union_, "Computes the union between two shapes.")
        .def("__or__", &Shape::union_, "Computes the union between two shapes.")
        .def("difference", &Shape::difference, "Computes the difference between two shapes.")
        .def("__sub__", &Shape::difference, "Computes the difference between two shapes.")
        .def("intersect", &Shape::intersect, "Computes the intersection of two shapes.")
        .def("__and__", &Shape::intersect, "Computes the intersection of two shapes.");

    m.def("add", &add, R"pbdoc(
        Add two numbers
        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers
        Some other explanation about the subtract function.
    )pbdoc");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
