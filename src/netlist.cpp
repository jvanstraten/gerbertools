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
 * Contains the logic for constructing physical netlists from coord::Paths
 * representations for copper layers and a list of via locations.
 */

#include <algorithm>
#include <sstream>
#include <cmath>
#include <limits>
#include <iostream> // TODO removeme
#include "gerbertools/netlist.hpp"
#include "gerbertools/clipper.hpp"
#include "gerbertools/path.hpp"

namespace gerbertools {
namespace netlist {

/**
 * Resolves a layer index where negative numbers are allowed Python-style to
 * a complete layer index using the total number of layers.
 */
static size_t resolve_layer_index(int layer, size_t num_layers) {
    int actual_layer = layer;
    if (actual_layer < 0) {
        actual_layer += (int)num_layers;
    }
    if (actual_layer < 0 || actual_layer >= (int)num_layers) {
        throw std::out_of_range("layer index " + std::to_string(layer) + " is out of range");
    }
    return (size_t)actual_layer;
}

/**
 * Constructs a new via object.
 */
Via::Via(
    const coord::Path &path,
    coord::CInt finished_hole_size,
    coord::CInt plating_thickness,
    int lower_layer,
    int upper_layer
) :
    path(path),
    finished_hole_size(finished_hole_size),
    plating_thickness(plating_thickness),
    lower_layer(lower_layer),
    upper_layer(upper_layer)
{}

/**
 * Returns the coordinate of the via.
 */
coord::CPt Via::get_coordinate() const {
    return path.at(0);
}

/**
 * Returns the entire milled path for the via.
 */
coord::Path Via::get_path() const {
    return path;
}

/**
 * Returns the size of the finished hole.
 */
coord::CInt Via::get_finished_hole_size() const {
    return finished_hole_size;
}

/**
 * Returns the size of the hole made in the dielectric substrate.
 */
coord::CInt Via::get_substrate_hole_size() const {
    return finished_hole_size + 2 * plating_thickness;
}

/**
 * Returns the plating thickness, such that the substrate hole size is
 * hole_size + 2 * thickness.
 */
coord::CInt Via::get_plating_thickness() const {
    return plating_thickness;
}

/**
 * Index of the first layer that this via can reach. Layer indices are from
 * 0 to N-1 for bottom to top.
 */
size_t Via::get_lower_layer(size_t num_layers) const {
    return resolve_layer_index(lower_layer, num_layers);
}

/**
 * Index of the last layer that this via can reach. Layer indices are from
 * 0 to N-1 for bottom to top.
 */
size_t Via::get_upper_layer(size_t num_layers) const {
    return resolve_layer_index(upper_layer, num_layers);
}

/**
 * Constructs a new copper shape.
 */
Shape::Shape(
    const coord::Path &outline,
    const coord::Paths &holes,
    size_t layer
) :
    outline(outline),
    holes(holes),
    layer(layer)
{
    bounding_box.left = bounding_box.right = outline.at(0).X;
    bounding_box.bottom = bounding_box.top = outline.at(0).Y;
    for (const auto &coord : outline) {
        bounding_box.left = std::min(bounding_box.left, coord.X);
        bounding_box.bottom = std::min(bounding_box.bottom, coord.Y);
        bounding_box.right = std::max(bounding_box.right, coord.X);
        bounding_box.top = std::max(bounding_box.top, coord.Y);
    }
}

/**
 * Outline of the copper shape.
 */
const coord::Path &Shape::get_outline() const {
    return outline;
}

/**
 * Zero or more holes in the copper shape.
 */
const coord::Paths &Shape::get_holes() const {
    return holes;
}

/**
 * Returns the layer index for this shape. Layer indices are from 0 to N-1
 * for bottom to top.
 */
size_t Shape::get_layer() const {
    return layer;
}

/**
 * Determines whether the given point is inside the shape.
 */
bool Shape::contains(coord::CPt point) const {
    if (point.X < bounding_box.left) return false;
    if (point.X > bounding_box.right) return false;
    if (point.Y < bounding_box.bottom) return false;
    if (point.Y > bounding_box.top) return false;
    if (ClipperLib::PointInPolygon(point, outline) == 0) return false;
    for (const auto &hole : holes) {
        if (ClipperLib::PointInPolygon(point, hole) == 1) return false;
    }
    return true;
}

/**
 * Constructs a physical net from the given initial shape.
 */
PhysicalNet::PhysicalNet(const ShapeRef &shape) : shapes({shape}) {
}

/**
 * Adds a via to this net.
 */
void PhysicalNet::add_via(const ViaRef &via) {
    vias.push_back(via);
}

/**
 * Moves the objects from the given net into this net in order to merge
 * them.
 */
void PhysicalNet::merge_with(const PhysicalNetRef &net) {
    shapes.splice(shapes.end(), net->shapes);
    vias.splice(vias.end(), net->vias);
    for (const auto &logical_net : net->logical_nets) {
        logical_nets.insert(logical_net);
    }
}

/**
 * Determines whether the given point on the given layer is part of this
 * net.
 */
bool PhysicalNet::contains(coord::CPt point, size_t layer) const {
    for (const auto &shape : shapes) {
        if (shape->get_layer() == layer) {
            if (shape->contains(point)) {
                return true;
            }
        }
    }
    return false;
}

/**
 * Assigns a logical net to this physical net.
 */
void PhysicalNet::assign_logical(const LogicalNetRef &logical_net) {
    logical_nets.insert(logical_net);
}

/**
 * Returns the list of copper shapes connected to this net.
 */
const std::list<ShapeRef> &PhysicalNet::get_shapes() const {
    return shapes;
}

/**
 * Returns the list of vias connected to this net.
 */
const std::list<ViaRef> &PhysicalNet::get_vias() {
    return vias;
}

/**
 * Returns the set of logical net names associated with this physical net.
 */
const LogicalNets &PhysicalNet::get_logical_nets() const {
    return logical_nets;
}

/**
 * Constructs a logical net.
 */
LogicalNet::LogicalNet(const std::string &name) : name(name) {}

/**
 * Returns the name of the net.
 */
const std::string &LogicalNet::get_name() const {
    return name;
}

/**
 * Assigns physical nets from both physical netlists to this logical net.
 */
void LogicalNet::assign_physical(const PhysicalNetRef &connected, const PhysicalNetRef &clearance) {
    connected_nets.insert(connected);
    clearance_nets.insert(clearance);
}

/**
 * Returns the set of physical nets from the connection netlist associated
 * with this logical net.
 */
const PhysicalNets &LogicalNet::get_connected_nets() const {
    return connected_nets;
}

/**
 * Returns the set of physical nets from the clearance netlist associated
 * with this logical net.
 */
const PhysicalNets &LogicalNet::get_clearance_nets() const {
    return clearance_nets;
}

/**
 * Creates a new physical netlist.
 */
PhysicalNetlist::PhysicalNetlist() : vias_added(false) {
}

/**
 * Adds a copper shape to the netlist. The netlist is updated accordingly.
 */
void PhysicalNetlist::register_shape(const ShapeRef &shape) {
    if (vias_added) throw std::runtime_error("cannot add shapes after the first via has been added");
    nets.push_back(std::make_shared<PhysicalNet>(shape));
}

/**
 * Helper function for PhysicalNetlist::register_paths(). Converts a PolyNode to
 * a Shape, adds it to the netlist, and calls itself for fully contained shapes.
 */
static void nodes_to_physical_netlist(const ClipperLib::PolyNodes &nodes, PhysicalNetlist &pnl, size_t layer) {
    for (const auto &node : nodes) {
        if (node->IsHole()) {
            throw std::runtime_error("shape is a hole?");
        }
        coord::Paths holes;
        for (const auto &hole : node->Childs) {
            if (!hole->IsHole()) {
                throw std::runtime_error("hole is not a hole?");
            }
            holes.push_back(hole->Contour);
            nodes_to_physical_netlist(hole->Childs, pnl, layer);
        }
        pnl.register_shape(std::make_shared<Shape>(node->Contour, holes, layer));
    }
}

/**
 * Like register_shape(), but using a set of closed paths interpreted using
 * odd-even winding as input.
 */
void PhysicalNetlist::register_paths(const coord::Paths &paths, size_t layer) {
    ClipperLib::Clipper cl;
    cl.StrictlySimple(true);
    cl.AddPaths(paths, ClipperLib::ptSubject, true);
    ClipperLib::PolyTree tree;
    cl.Execute(ClipperLib::ctUnion, tree);
    nodes_to_physical_netlist(tree.Childs, *this, layer);
}

/**
 * Adds a via to the netlist. The netlist is updated accordingly. Returns
 * true if successful, or false if any of the specified layers do not have
 * copper at the center point of the via.
 */
bool PhysicalNetlist::register_via(const ViaRef &via, size_t num_layers) {
    size_t lower_layer = via->get_lower_layer(num_layers);
    size_t upper_layer = via->get_upper_layer(num_layers);
    if (lower_layer >= upper_layer) {
        throw std::runtime_error("via has null layer range or only includes one layer");
    }
    vias_added = true;
    bool ok = true;
    PhysicalNetRef target = {};
    for (size_t layer = lower_layer; layer <= upper_layer; layer++) {
        auto source = find_net(via->get_coordinate(), layer);
        if (!source) {
            ok = false;
            continue;
        } else if (!target) {
            target = source;
        } else {
            target->merge_with(source);
            nets.remove(source);
        }
    }
    if (target) {
        target->add_via(via);
    }
    return ok;
}

/**
 * Returns which net belongs to the given point. The returned pointer will
 * be null if there is no (virtual) copper at the given point.
 */
PhysicalNetRef PhysicalNetlist::find_net(coord::CPt point, size_t layer) const {
    for (const auto &net : nets) {
        if (net->contains(point, layer)) {
            return net;
        }
    }
    return {};
}

/**
 * Returns the list of all physical nets.
 */
const std::list<PhysicalNetRef> &PhysicalNetlist::get_nets() const {
    return nets;
}

/**
 * Constructs a connection point.
 */
ConnectionPoint::ConnectionPoint(
    coord::CPt coordinate,
    int layer,
    const LogicalNetRef &net
) :
    coordinate(coordinate),
    layer(layer),
    net(net)
{}

/**
 * Returns the coordinate of the connection point.
 */
coord::CPt ConnectionPoint::get_coordinate() const {
    return coordinate;
}

/**
 * Layer index for the connection point. Layer indices are from 0 to N-1 for
 * bottom to top.
 */
size_t ConnectionPoint::get_layer(size_t num_layers) const {
    return resolve_layer_index(layer, num_layers);
}

/**
 * Returns the logical net for this connection point.
 */
LogicalNetRef ConnectionPoint::get_net() const {
    return net;
}

/**
 * Adds a copper layer. The layers are added bottom-up. All layers must be
 * added before vias are added.
 */
NetlistBuilder &NetlistBuilder::layer(const coord::Paths &paths) {
    layers.push_back(paths);
    return *this;
}

/**
 * Adds a via defined by a milling path, finished hole size, and plating
 * thickness, reaching from and to the given layers. Layer indices are from
 * 0 to N-1 or from -N to -1 for bottom to top.
 */
NetlistBuilder &NetlistBuilder::via(
    const coord::Path &path,
    coord::CInt finished_hole_size,
    coord::CInt plating_thickness,
    int lower_layer,
    int upper_layer
) {
    vias.push_back(std::make_shared<Via>(
        path,
        finished_hole_size,
        plating_thickness,
        lower_layer,
        upper_layer
    ));
    return *this;
}

/**
 * Associates a point on the PCB with a logical net name.
 */
NetlistBuilder &NetlistBuilder::net(coord::CPt point, int layer, const std::string &net_name) {
    auto it = nets.find(net_name);
    if (it == nets.end()) {
        it = nets.insert({net_name, std::make_shared<LogicalNet>(net_name)}).first;
    }
    connections.emplace_back(point, layer, it->second);
    return *this;
}

/**
 * Builds the netlist with the given clearance. The builder should not be
 * used after this point.
 */
Netlist NetlistBuilder::build(coord::CInt clearance) {
    Netlist nl;

    // Convert and add the copper shapes for all layers.
    nl.num_layers = 0;
    for (const auto &paths : layers) {
        nl.connected_netlist.register_paths(paths, nl.num_layers);
        auto extended_paths = path::offset(paths, (double)clearance / 2, false);
        nl.clearance_netlist.register_paths(extended_paths, nl.num_layers);
        nl.num_layers++;
    }

    // Register vias.
    for (const auto &via : vias) {
        if (!nl.connected_netlist.register_via(via, nl.num_layers)) {
            std::ostringstream ss;
            ss << "via at coordinate (";
            ss << coord::Format::to_mm(via->get_coordinate().X);
            ss << ", ";
            ss << coord::Format::to_mm(via->get_coordinate().Y);
            ss << ") is not connected to copper on one or more layers";
            nl.builder_violations.push_back(ss.str());
        }
        nl.clearance_netlist.register_via(via, nl.num_layers);
    }

    // Register logical nets.
    nl.logical_nets = std::move(nets);

    // Register connection points to connect the logical and physical nets.
    for (const auto &connection : connections) {
        const auto coord = connection.get_coordinate();
        const auto layer = connection.get_layer(nl.num_layers);
        const auto logical_net = connection.get_net();
        auto connected_net = nl.connected_netlist.find_net(coord, layer);
        if (!connected_net) {
            std::ostringstream ss;
            ss << "connection at coordinate (";
            ss << coord::Format::to_mm(connection.get_coordinate().X);
            ss << ", ";
            ss << coord::Format::to_mm(connection.get_coordinate().Y);
            ss << ") on layer ";
            ss << coord::Format::to_mm(connection.get_layer(nl.num_layers));
            ss << " should be connected to logical net ";
            ss << logical_net->get_name();
            ss << ", but there is no copper here";
            nl.builder_violations.push_back(ss.str());
            continue;
        }
        auto clearance_net = nl.clearance_netlist.find_net(coord, layer);
        if (!clearance_net) {
            throw std::logic_error("point maps to connected netlist but not to clearance netlist");
        }
        connected_net->assign_logical(logical_net);
        clearance_net->assign_logical(logical_net);
        logical_net->assign_physical(connected_net, clearance_net);
    }

    return nl;
}

/**
 * Returns the square of the distance between the given point and the given line
 * segment.
 */
static double point_to_line_distance_sqr(coord::CPt point, coord::CPt a, coord::CPt b) {
    double A = point.X - a.X;
    double B = point.X - a.Y;
    double C = b.X - a.X;
    double D = b.Y - a.Y;

    double dot = A * C + B * D;
    double len_sq = C * C + D * D;
    double param = dot / len_sq;

    double xx, yy;

    if (param < 0 || (a.X == b.X && a.Y == b.Y)) {
        xx = a.X;
        yy = a.Y;
    } else if (param > 1) {
        xx = b.X;
        yy = b.Y;
    } else {
        xx = a.X + param * C;
        yy = a.Y + param * D;
    }

    double dx = point.X - xx;
    double dy = point.Y - yy;

    return dx * dx + dy * dy;
}

/**
 * Returns the square of the distance between the given point and a closed path.
 *
 * NOTE: this is not exactly the most efficient implementation.
 */
static double point_to_path_distance_sqr(coord::CPt point, coord::Path p, bool closed) {
    double r_sqr_min = std::numeric_limits<double>::infinity();
    if (closed) {
        point_to_line_distance_sqr(point, p.front(), p.back());
    }
    for (size_t i = 1; i < p.size(); i++) {
        r_sqr_min = std::min(r_sqr_min, point_to_line_distance_sqr(point, p.at(i-1), p.at(i)));
    }
    return r_sqr_min;
}

/**
 * Returns the minimum annular ring for the given via.
 */
static double compute_annular_ring(const ViaRef &via, const PhysicalNetRef &net) {
    double r_sqr_min = std::numeric_limits<double>::infinity();
    for (const auto &vc : via->get_path()) {
        for (const auto &shape : net->get_shapes()) {
            r_sqr_min = std::min(r_sqr_min, point_to_path_distance_sqr(vc, shape->get_outline(), true));
            for (const auto &hole : shape->get_holes()) {
                r_sqr_min = std::min(r_sqr_min, point_to_path_distance_sqr(vc, hole, true));
            }
        }
    }
    if (via->get_path().size() > 1) {
        for (const auto &shape : net->get_shapes()) {
            for (const auto &pt : shape->get_outline()) {
                r_sqr_min = std::min(r_sqr_min, point_to_path_distance_sqr(pt, via->get_path(), false));
            }
            for (const auto &hole : shape->get_holes()) {
                for (const auto &pt : hole) {
                    r_sqr_min = std::min(r_sqr_min, point_to_path_distance_sqr(pt, via->get_path(), false));
                }
            }
        }
    }
    return std::sqrt(r_sqr_min) - (double)via->get_finished_hole_size() / 2;
}

/**
 * Runs the design-rule check. A list of violation messages is returned. If
 * any message is returned, the DRC failed.
 */
std::list<std::string> Netlist::perform_drc(coord::CInt annular_ring) const {
    auto violations = builder_violations;

    // Check for open circuits.
    for (const auto &it : logical_nets) {
        auto net = it.second;
        auto pnets = net->get_connected_nets();
        if (pnets.empty()) {
            std::ostringstream ss;
            ss << "logical net " << net->get_name() << " is completely unrouted";
            violations.push_back(ss.str());
        } else if (pnets.size() > 1) {
            std::ostringstream ss;
            ss << "logical net " << net->get_name() << " is divided up into ";
            ss << pnets.size() << " islands";
            violations.push_back(ss.str());
        }
    }

    // Check for short circuits.
    std::set<std::pair<std::string, std::string>> shorts_reported;
    for (const auto &net : connected_netlist.get_nets()) {
        auto lnets = net->get_logical_nets();
        if (lnets.size() < 2) continue;
        for (const auto &a : lnets) {
            auto net_a = a.lock()->get_name();
            for (const auto &b : lnets) {
                auto net_b = b.lock()->get_name();
                if (net_a == net_b) continue;
                if (shorts_reported.count({net_a, net_b})) continue;
                shorts_reported.insert({net_a, net_b});
                shorts_reported.insert({net_b, net_a});
                std::ostringstream ss;
                ss << "logical nets " << net_a << " and " << net_b << " are short-circuited";
                violations.push_back(ss.str());
            }
        }
    }

    // Check for clearance violations.
    for (const auto &net : clearance_netlist.get_nets()) {
        auto lnets = net->get_logical_nets();
        if (lnets.size() < 2) continue;
        for (const auto &a : lnets) {
            auto net_a = a.lock()->get_name();
            for (const auto &b : lnets) {
                auto net_b = b.lock()->get_name();
                if (net_a == net_b) continue;
                if (shorts_reported.count({net_a, net_b})) continue;
                shorts_reported.insert({net_a, net_b});
                shorts_reported.insert({net_b, net_a});
                std::ostringstream ss;
                ss << "clearance violation between " << net_a << " and " << net_b;
                violations.push_back(ss.str());
            }
        }
    }

    // Check for annular ring violations.
    for (const auto &net : connected_netlist.get_nets()) {
        for (const auto &via : net->get_vias()) {
            auto via_annular_ring = compute_annular_ring(via, net);
            if (via_annular_ring < annular_ring) {
                std::ostringstream ss;
                ss << "via at coordinate (";
                ss << coord::Format::to_mm(via->get_coordinate().X);
                ss << ", ";
                ss << coord::Format::to_mm(via->get_coordinate().Y);
                ss << ") has annular ring ";
                ss << coord::Format::to_mm(via_annular_ring);
                ss << ", less than the minimum ";
                ss << coord::Format::to_mm(annular_ring);
                violations.push_back(ss.str());
            }
        }
    }

    return violations;
}

/**
 * Returns the physical nets in the netlist. That is, the pieces of
 * connected copper and their shape, as well as references to the logical
 * nets they are connected to (if they are electrically connected to
 * anything but themselves at all).
 */
const PhysicalNetlist &Netlist::get_physical_netlist() const {
    return connected_netlist;
}

/**
 * Returns the map from logical netname (as in from the circuit) to objects
 * that store which physical nets are mapped to them.
 */
const std::map<std::string, LogicalNetRef> &Netlist::get_logical_netlist() const {
    return logical_nets;
}

} // namespace netlist
} // namespace gerbertools
