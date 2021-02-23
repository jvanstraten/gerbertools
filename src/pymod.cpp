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
#include <pybind11/stl.h>
#include "gerbertools/coord.hpp"
#include "gerbertools/path.hpp"
#include "gerbertools/pcb.hpp"

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
        if (index < -(int64_t)p.size() || index >= (int64_t)p.size()) {
            throw std::domain_error("index " + std::to_string(index) + " is out of range, size is " + std::to_string(p.size()));
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

    void set_miter_limit(double x) {
        miter_limit = x;
    }

    double get_arc_tolerance() const {
        return arc_tolerance;
    }

    void set_arc_tolerance(double x) {
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

using Color = std::tuple<float, float, float, float>;

Color color_to_tuple(const color::Color &c) {
    return {c.r, c.g, c.b, c.a};
}

color::Color tuple_to_color(const Color &c) {
    return {std::get<0>(c), std::get<1>(c), std::get<2>(c), std::get<3>(c)};
}

class Netlist {
private:
    netlist::Netlist nl;
    coord::CInt annular_ring;
public:

    Netlist() {
        throw std::runtime_error("Netlists cannot be constructed directly. Use CircuitBoard.build_netlist() instead.");
    }

    Netlist(netlist::Netlist &&nl, coord::CInt annular_ring) : nl(std::move(nl)), annular_ring(annular_ring) {
    }

    std::list<std::string> drc() {
        return nl.perform_drc(annular_ring);
    }

};

class CircuitBoard {
private:
    pcb::CircuitBoard pcb;
public:

    CircuitBoard(
        const std::string &basename,
        const std::string &outline,
        const std::string &drill,
        const std::string &drill_nonplated,
        const std::string &mill,
        double plating_thickness
    ) : pcb(basename, outline, drill, drill_nonplated, mill, plating_thickness) {}

    void add_mask_layer(const std::string &mask, const std::string &silk) {
        pcb.add_mask_layer(mask, silk);
    }

    void add_copper_layer(const std::string &gerber, double thickness) {
        pcb.add_copper_layer(gerber, thickness);
    }

    void add_substrate_layer(double thickness) {
        pcb.add_substrate_layer(thickness);
    }

    void add_surface_finish() {
        pcb.add_surface_finish();
    }

    Netlist build_netlist(const std::list<std::tuple<std::tuple<double, double>, int, std::string>> &nets, double clearance, double annular_ring) {
        auto nb = pcb.get_netlist_builder();
        for (const auto &net : nets) {
            auto coord = std::get<0>(net);
            auto layer = std::get<1>(net);
            auto name = std::get<2>(net);
            auto x = coord::Format::from_mm(std::get<0>(coord));
            auto y = coord::Format::from_mm(std::get<1>(coord));
            nb.net({x, y}, layer, name);
        }
        return Netlist(nb.build(coord::Format::from_mm(clearance)), coord::Format::from_mm(annular_ring));
    }

    std::tuple<double, double, double, double> get_bounds() {
        auto bounds = pcb.get_bounds();
        return {
            coord::Format::to_mm(bounds.left),
            coord::Format::to_mm(bounds.top),
            coord::Format::to_mm(bounds.right),
            coord::Format::to_mm(bounds.bottom)
        };
    }

    std::string get_svg(
        bool flipped,
        const Color &soldermask,
        const Color &silkscreen,
        const Color &finish,
        const Color &substrate,
        const Color &copper,
        const std::string &id_mask
    ) {
        return pcb.get_svg(flipped, pcb::ColorScheme(
            tuple_to_color(soldermask),
            tuple_to_color(silkscreen),
            tuple_to_color(finish),
            tuple_to_color(substrate),
            tuple_to_color(copper)
        ), id_mask);
    }

    void write_svg(
        const std::string &fname,
        bool flipped,
        double scale,
        const Color &soldermask,
        const Color &silkscreen,
        const Color &finish,
        const Color &substrate,
        const Color &copper
    ) {
        pcb.write_svg(fname, flipped, scale, pcb::ColorScheme(
            tuple_to_color(soldermask),
            tuple_to_color(silkscreen),
            tuple_to_color(finish),
            tuple_to_color(substrate),
            tuple_to_color(copper)
        ));
    }

    void write_obj(const std::string &fname) {
        pcb.write_obj(fname);
    }
};

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

    py::class_<Netlist>(m, "Netlist", "A representation of the copper part of a circuit board and its relationship with the logical netlist.")
        .def("drc", &Netlist::drc, "Runs a design-rule check. Returns a list of violation messages. If the list is empty, the DRC passed.");

    py::class_<CircuitBoard>(m, "CircuitBoard", "A complete circuit board.")
        .def(py::init<std::string, std::string, std::string, std::string, std::string, double>(),
            py::arg("basename"),
            py::arg("outline") = ".GKO",
            py::arg("drill") = ".TXT",
            py::arg("drill_nonplated") = "",
            py::arg("mill") = "",
            py::arg("plating_thickness") = 0.5 * pcb::COPPER_OZ,
            R"doc(
                Constructs a circuit board. The following files are to be specified:
                 - basename: prefix for all filenames.
                 - outline: the board outline Gerber file (GKO, GM1, etc). May also include milling information.
                 - drill: the NC drill file (TXT).
                 - drill_nonplated: if specified, non-plated holes will be added from this auxiliary NC drill file.
                 - mill: if specified, adds another Gerber-based milling layer. Interpreted in the same way as outline.
                 - plating_thickness: thickness of the hole plating in millimeters.
            )doc")
        .def("add_mask_layer", &CircuitBoard::add_mask_layer, py::arg("mask"), py::arg("silk")="", "Adds a soldermask/silkscreen layer. Layers are added bottom-up.")
        .def("add_copper_layer", &CircuitBoard::add_copper_layer, py::arg("copper"), py::arg("thickness")=pcb::COPPER_OZ, "Adds a copper layer. Layers are added bottom-up.")
        .def("add_substrate_layer", &CircuitBoard::add_substrate_layer, py::arg("thickness")=1.5, "Adds a substrate layer. Layers are added bottom-up.")
        .def("add_surface_finish", &CircuitBoard::add_surface_finish, "Derives the surface finish layer for all exposed copper. Call after adding all layers.")
        .def("build_netlist", &CircuitBoard::build_netlist, py::arg("nets"), py::arg("clearance")=0.0, py::arg("annular_ring")=0.0, "Builds a netlist from the PCB and the given list of coordinate, layer index, and net name tuples.")
        .def("get_bounds", &CircuitBoard::get_bounds, "Returns the axis-aligned bounds of the PCB in millimeters, ordered left, top, right, bottom.")
        .def("get_svg", &CircuitBoard::get_svg,
             py::arg("flipped")=false,
             py::arg("soldermask")=color_to_tuple(color::MASK_GREEN),
             py::arg("silkscreen")=color_to_tuple(color::SILK_WHITE),
             py::arg("finish")=color_to_tuple(color::FINISH_TIN),
             py::arg("substrate")=color_to_tuple(color::SUBSTRATE),
             py::arg("copper")=color_to_tuple(color::COPPER),
             py::arg("id_prefix")="",
             "Renders the circuit board to SVG, returning only its body as a string."
        )
        .def("write_svg", &CircuitBoard::write_svg,
             py::arg("fname"),
             py::arg("flipped")=false,
             py::arg("scale")=1.0,
             py::arg("soldermask")=color_to_tuple(color::MASK_GREEN),
             py::arg("silkscreen")=color_to_tuple(color::SILK_WHITE),
             py::arg("finish")=color_to_tuple(color::FINISH_TIN),
             py::arg("substrate")=color_to_tuple(color::SUBSTRATE),
             py::arg("copper")=color_to_tuple(color::COPPER),
             "Renders the circuit board to an SVG file."
         )
        .def("write_obj", &CircuitBoard::write_obj,
             py::arg("fname"),
             "Renders the circuit board to a Wavefront OBJ file."
         );

    auto m2 = m.def_submodule("color", "Defines some default colors.");
    m2.def("copper", []() { return color_to_tuple(color::COPPER); }, "Returns the default color for copper.");
    m2.def("finish_tin", []() { return color_to_tuple(color::FINISH_TIN); }, "Returns the default color for tin/HASL surface finish.");
    m2.def("substrate", []() { return color_to_tuple(color::SUBSTRATE); }, "Returns the default color for an FR-4 substrate.");
    m2.def("mask_green", []() { return color_to_tuple(color::MASK_GREEN); }, "Returns the default color for green soldermask.");
    m2.def("mask_white", []() { return color_to_tuple(color::MASK_WHITE); }, "Returns the default color for white soldermask.");
    m2.def("silk_white", []() { return color_to_tuple(color::SILK_WHITE); }, "Returns the default color for white silkscreen.");
    m2.def("silk_black", []() { return color_to_tuple(color::SILK_BLACK); }, "Returns the default color for black silkscreen.");

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
