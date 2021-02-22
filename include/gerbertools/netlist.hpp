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

#pragma once

#include <memory>
#include <string>
#include <list>
#include <set>
#include <map>
#include "gerbertools/coord.hpp"

namespace gerbertools {

/**
 * Namespace for netlist construction classes.
 */
namespace netlist {

/**
 * Represents a via.
 */
class Via {
private:

    /**
     * The path that the drill took to make the via or slot.
     */
    coord::Path path;

    /**
     * Size of the finished hole.
     */
    coord::CInt finished_hole_size;

    /**
     * Plating thickness, such that the substrate hole size is
     * hole_size + 2 * thickness.
     */
    coord::CInt plating_thickness;

    /**
     * Index of the lowest layer that this via can reach. Layer indices are from
     * 0 to N-1 or from -N to -1 for bottom to top.
     */
    int lower_layer;

    /**
     * Index of the topmost layer that this via can reach. Layer indices are
     * from 0 to N-1 or from -N to -1 for bottom to top.
     */
    int upper_layer;

public:

    /**
     * Constructs a new via object.
     */
    Via(
        const coord::Path &path,
        coord::CInt finished_hole_size,
        coord::CInt plating_thickness,
        int lower_layer=0,
        int upper_layer=-1
    );

    /**
     * Returns the start coordinate of the via.
     */
    coord::CPt get_coordinate() const;

    /**
     * Returns the entire milled path for the via.
     */
    coord::Path get_path() const;

    /**
     * Returns the size of the finished hole.
     */
    coord::CInt get_finished_hole_size() const;

    /**
     * Returns the size of the hole made in the dielectric substrate.
     */
    coord::CInt get_substrate_hole_size() const;

    /**
     * Returns the plating thickness, such that the substrate hole size is
     * hole_size + 2 * thickness.
     */
    coord::CInt get_plating_thickness() const;

    /**
     * Index of the first layer that this via can reach. Layer indices are from
     * 0 to N-1 for bottom to top.
     */
    size_t get_lower_layer(size_t num_layers) const;

    /**
     * Index of the last layer that this via can reach. Layer indices are from
     * 0 to N-1 for bottom to top.
     */
    size_t get_upper_layer(size_t num_layers) const;

};

/**
 * Represents a copper shape on one of the copper layers of the PCB, not
 * connected to any other bits of copper in that layer without involvement of
 * other layers.
 */
class Shape {
private:

    /**
     * Outline of the copper shape.
     */
    coord::Path outline;

    /**
     * Zero or more holes in the copper shape.
     */
    coord::Paths holes;

    /**
     * Axis-aligned bounding box of the shape, for accelerating intersection
     * tests.
     */
    coord::CRect bounding_box;

    /**
     * Layer index for this shape. Layer indices are from 0 to N-1 for bottom
     * to top.
     */
    size_t layer;

public:

    /**
     * Constructs a new copper shape.
     */
    Shape(
        const coord::Path &outline,
        const coord::Paths &holes,
        size_t layer
    );

    /**
     * Outline of the copper shape.
     */
    const coord::Path &get_outline() const;

    /**
     * Zero or more holes in the copper shape.
     */
    const coord::Paths &get_holes() const;

    /**
     * Returns the layer index for this shape. Layer indices are from 0 to N-1
     * for bottom to top.
     */
    size_t get_layer() const;

    /**
     * Determines whether the given point is inside the shape.
     */
    bool contains(coord::CPt point) const;

};

/**
 * Compare function for weak pointers by what they point to.
 */
template <typename T>
struct WeakPointerCompare {
    bool operator() (const T &lhs, const T &rhs) const {
        auto rptr = rhs.lock();
        if (!rptr) return false;
        auto lptr = lhs.lock();
        if (!lptr) return true;
        return (unsigned long)lptr.get() < (unsigned long)rptr.get();
    }
};

/**
 * Compare function for shared pointers by what they point to.
 */
template <typename T>
struct SharedPointerCompare {
    bool operator() (const T &lptr, const T &rptr) const {
        if (!rptr) return false;
        if (!lptr) return true;
        return (unsigned long)lptr.get() < (unsigned long)rptr.get();
    }
};

// Forward declarations for physical and logical nets such that we can make
// reference typedefs for them.
class PhysicalNet;
class LogicalNet;

/**
 * Reference to a shape object.
 */
using ShapeRef = std::shared_ptr<Shape>;

/**
 * Reference to a via object.
 */
using ViaRef = std::shared_ptr<Via>;

/**
 * Reference to a physical net.
 */
using PhysicalNetRef = std::shared_ptr<PhysicalNet>;

/**
 * Set of physical nets.
 */
using PhysicalNets = std::set<PhysicalNetRef, SharedPointerCompare<PhysicalNetRef>>;

/**
 * Reference to a logical net.
 */
using LogicalNetRef = std::shared_ptr<LogicalNet>;

/**
 * Weak reference to a logical net.
 */
using LogicalNetWeakRef = std::weak_ptr<LogicalNet>;

/**
 * Set of logical nets.
 */
using LogicalNets = std::set<LogicalNetWeakRef, WeakPointerCompare<LogicalNetWeakRef>>;

/**
 * Represents a physical net, i.e. a set of connected copper shapes and vias
 * that are physically connected on the PCB.
 */
class PhysicalNet {

    /**
     * The copper shapes connected to this net.
     */
    std::list<ShapeRef> shapes;

    /**
     * The vias connected to this net.
     */
    std::list<ViaRef> vias;

    /**
     * Weak references to the logical nets that are associated to this net.
     */
    LogicalNets logical_nets;

public:

    /**
     * Constructs a physical net from the given initial shape.
     */
    PhysicalNet(const ShapeRef &shape);

    /**
     * Adds a via to this net.
     */
    void add_via(const ViaRef &via);

    /**
     * Moves the objects from the given net into this net in order to merge
     * them.
     */
    void merge_with(const PhysicalNetRef &net);

    /**
     * Determines whether the given point on the given layer is part of this
     * net.
     */
    bool contains(coord::CPt point, size_t layer) const;

    /**
     * Assigns a logical net to this physical net.
     */
    void assign_logical(const LogicalNetRef &logical_net);

    /**
     * Returns the list of copper shapes connected to this net.
     */
    const std::list<ShapeRef> &get_shapes() const;

    /**
     * Returns the list of vias connected to this net.
     */
    const std::list<ViaRef> &get_vias();

    /**
     * Returns the set of logical net names associated with this physical net.
     */
    const LogicalNets &get_logical_nets() const;

};

/**
 * Represents a logical net from the netlist from the circuit design tool, of
 * which the correspondence with the PCB is to be verified.
 */
class LogicalNet {
private:

    /**
     * Name for this net.
     */
    std::string name;

    /**
     * Weak references to the physical nets that are associated to this net by
     * direct connection.
     */
    PhysicalNets connected_nets;

    /**
     * Weak references to the physical nets that are associated to this net due
     * to clearance violations.
     */
    PhysicalNets clearance_nets;

public:

    /**
     * Constructs a logical net.
     */
    LogicalNet(const std::string &name);

    /**
     * Returns the name of the net.
     */
    const std::string &get_name() const;

    /**
     * Assigns physical nets from both physical netlists to this logical net.
     */
    void assign_physical(const PhysicalNetRef &connected, const PhysicalNetRef &clearance);

    /**
     * Returns the set of physical nets from the connection netlist associated
     * with this logical net.
     */
    const PhysicalNets &get_connected_nets() const;

    /**
     * Returns the set of physical nets from the clearance netlist associated
     * with this logical net.
     */
    const PhysicalNets &get_clearance_nets() const;

};

/**
 * Represents a complete physical netlist derived from the PCB. Two of these are
 * used for electrical DRC: the connection netlist is built from the PCB as-is,
 * while the clearance netlist is built from the PCB after expansion of all
 * copper shapes by 1/2 the minimum clearance.
 */
class PhysicalNetlist {
private:

    /**
     * The list of nets in this netlist.
     */
    std::list<PhysicalNetRef> nets;

    /**
     * Records whether we've started adding vias.
     */
    bool vias_added;

public:

    /**
     * Creates a new physical netlist.
     */
    PhysicalNetlist();

    /**
     * Adds a copper shape to the netlist. The netlist is updated accordingly.
     * All shapes must be added before vias are added.
     */
    void register_shape(const ShapeRef &shape);

    /**
     * Like register_shape(), but using a set of closed paths for a given copper
     * layer interpreted using odd-even winding as input.
     */
    void register_paths(const coord::Paths &paths, size_t layer);

    /**
     * Adds a via to the netlist. The netlist is updated accordingly. Returns
     * true if successful, or false if any of the specified layers do not have
     * copper at the center point of the via.
     */
    bool register_via(const ViaRef &via, size_t num_layers);

    /**
     * Returns which net belongs to the given point. The returned pointer will
     * be null if there is no (virtual) copper at the given point.
     */
    PhysicalNetRef find_net(coord::CPt point, size_t layer) const;

    /**
     * Returns the list of all physical nets.
     */
    const std::list<PhysicalNetRef> &get_nets() const;

};

/**
 * Represents a connection point for a logical net.
 */
class ConnectionPoint {
private:

    /**
     * The coordinate for the connection point.
     */
    coord::CPt coordinate;

    /**
     * Layer index for the connection point. Layer indices are from 0 to N-1 or
     * from -N to -1 for bottom to top.
     */
    int layer;

    /**
     * The logical net that should be connected to this point.
     */
    LogicalNetRef net;

public:

    /**
     * Constructs a connection point.
     */
    ConnectionPoint(coord::CPt coordinate, int layer, const LogicalNetRef &net);

    /**
     * Returns the coordinate of the connection point.
     */
    coord::CPt get_coordinate() const;

    /**
     * Layer index for the connection point. Layer indices are from 0 to N-1 for
     * bottom to top.
     */
    size_t get_layer(size_t num_layers) const;

    /**
     * Returns the logical net for this connection point.
     */
    LogicalNetRef get_net() const;

};

/**
 * Forward reference for Netlist for the builder.
 */
class Netlist;

/**
 * Builder class for a complete netlist.
 */
class NetlistBuilder {
private:

    /**
     * The copper layers.
     */
    std::list<coord::Paths> layers;

    /**
     * The vias connecting the layers together.
     */
    std::list<ViaRef> vias;

    /**
     * Mapping from net names to logical nets.
     */
    std::map<std::string, LogicalNetRef> nets;

    /**
     * All the connection points.
     */
    std::list<ConnectionPoint> connections;

public:

    /**
     * Adds a copper layer. The layers are added bottom-up.
     */
    NetlistBuilder &layer(const coord::Paths &paths);

    /**
     * Adds a via defined by a milling path, finished hole size, and plating
     * thickness, reaching from and to the given layers. Layer indices are from
     * 0 to N-1 or from -N to -1 for bottom to top.
     */
    NetlistBuilder &via(
        const coord::Path &path,
        coord::CInt finished_hole_size,
        coord::CInt plating_thickness,
        int lower_layer=0,
        int upper_layer=-1
    );

    /**
     * Associates a point on the PCB with a logical net name.
     */
    NetlistBuilder &net(coord::CPt point, int layer, const std::string &net_name);

    /**
     * Builds the netlist with the given clearance. The builder should not be
     * used after this point.
     */
    Netlist build(coord::CInt clearance);

};

/**
 * A complete netlist.
 */
class Netlist {
private:
    friend class NetlistBuilder;

    /**
     * The number of copper layers on the board.
     */
    size_t num_layers;

    /**
     * The physical netlist constructed from the PCB as-is.
     */
    PhysicalNetlist connected_netlist;

    /**
     * The physical netlist constructed from the PCB with 1/2 clearance
     * expansion of all copper shapes.
     */
    PhysicalNetlist clearance_netlist;

    /**
     * List of DRC violations detected while the netlist was built.
     */
    std::list<std::string> builder_violations;

    /**
     * The list of nets in this netlist.
     */
    std::map<std::string, LogicalNetRef> logical_nets;

public:

    /**
     * Runs the design-rule check. A list of violation messages is returned. If
     * any message is returned, the DRC failed.
     */
    std::list<std::string> perform_drc(coord::CInt annular_ring) const;

    /**
     * Returns the physical nets in the netlist. That is, the pieces of
     * connected copper and their shape, as well as references to the logical
     * nets they are connected to (if they are electrically connected to
     * anything but themselves at all).
     */
    const PhysicalNetlist &get_physical_netlist() const;

    /**
     * Returns the map from logical netname (as in from the circuit) to objects
     * that store which physical nets are mapped to them.
     */
    const std::map<std::string, LogicalNetRef> &get_logical_netlist() const;

};

} // namespace netlist
} // namespace gerbertools
