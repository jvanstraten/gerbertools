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
 * Contains the CircuitBoard class, representing a complete circuit board as
 * reconstructed from Gerber files and an NC drill file.
 */

#define _USE_MATH_DEFINES
#include <cmath>
#include "gerbertools/gerber.hpp"
#include "gerbertools/pcb.hpp"
#include "gerbertools/path.hpp"

namespace gerbertools {
namespace pcb {

/**
 * Creates a new color scheme.
 */
ColorScheme::ColorScheme(
    const color::Color &soldermask,
    const color::Color &silkscreen,
    const color::Color &finish,
    const color::Color &substrate,
    const color::Color &copper
) :
    soldermask(soldermask),
    silkscreen(silkscreen),
    finish(finish),
    substrate(substrate),
    copper(copper)
{}

/**
 * Constructs a layer.
 */
Layer::Layer(const std::string &name, double thickness) : name(name), thickness(thickness) {
}

/**
 * Returns a name for the layer.
 */
std::string Layer::get_name() const {
    return name;
}

/**
 * Returns the thickness of this layer.
 */
double Layer::get_thickness() const {
    return thickness;
}

/**
 * Constructs a substrate layer.
 */
SubstrateLayer::SubstrateLayer(
    const std::string &name,
    const coord::Paths &shape,
    const coord::Paths &dielectric,
    const coord::Paths &plating,
    double thickness
) : Layer(name, thickness), shape(shape), dielectric(dielectric), plating(plating) {
}

/**
 * Returns the surface finish mask for this layer.
 */
coord::Paths SubstrateLayer::get_mask() const {
    return dielectric;
}

/**
 * Renders the layer to an SVG layer.
 */
svg::Layer SubstrateLayer::to_svg(const ColorScheme &colors, bool flipped, const std::string &id_prefix) const {
    auto layer = svg::Layer(id_prefix + get_name());
    layer.add(dielectric, colors.substrate);
    layer.add(plating, colors.finish);
    return layer;
}

/**
 * Renders the layer to an OBJ file.
 */
void SubstrateLayer::to_obj(obj::ObjFile &obj, size_t layer_index, double z, const std::string &id_prefix) const {
    obj.add_object(
        "layer" + std::to_string(layer_index) + "_" + get_name(),
        "substrate"
    ).add_sheet(
        dielectric,
        z,
        z + get_thickness()
    );
}

/**
 * Constructs a copper layer.
 */
CopperLayer::CopperLayer(
    const std::string &name,
    const coord::Paths &board_shape,
    const coord::Paths &board_shape_excl_pth,
    const coord::Paths &copper_layer,
    double thickness
) : Layer(name, thickness) {
    layer = copper_layer;
    copper = path::intersect(board_shape, copper_layer);
    copper_excl_pth = path::intersect(board_shape_excl_pth, copper_layer);
}

/**
 * Returns the surface finish mask for this layer.
 */
coord::Paths CopperLayer::get_mask() const {
    return copper;
}

/**
 * Returns the copper for this layer.
 */
const coord::Paths &CopperLayer::get_copper() const {
    return copper;
}

/**
 * Returns the copper for this layer, with only cutouts for board shape and
 * non-plated holes, not for plated holes. This is needed for annular ring DRC.
 */
const coord::Paths &CopperLayer::get_copper_excl_pth() const {
    return copper_excl_pth;
}

/**
 * Returns the original layer, without board outline intersection.
 */
const coord::Paths &CopperLayer::get_layer() const {
    return layer;
}

/**
 * Renders the layer to an SVG layer.
 */
svg::Layer CopperLayer::to_svg(const ColorScheme &colors, bool flipped, const std::string &id_prefix) const {
    auto layer = svg::Layer(id_prefix + get_name());
    layer.add(copper, colors.copper);
    return layer;
}

/**
 * Renders the layer to an OBJ file.
 */
void CopperLayer::to_obj(obj::ObjFile &obj, size_t layer_index, double z, const std::string &id_prefix) const {
    // No-operation. Copper is added via the physical netlist, such that each
    // bit of connected copper gets its own object.
}

/**
 * Constructs a solder mask.
 */
MaskLayer::MaskLayer(
    const std::string &name,
    const coord::Paths &board_outline,
    const coord::Paths &mask_layer,
    const coord::Paths &silk_layer,
    bool bottom
) : Layer(name, 0.01), bottom(bottom) {
    mask = path::subtract(board_outline, mask_layer);
    silk = path::intersect(mask, silk_layer);
}

/**
 * Returns the surface finish mask for this layer.
 */
coord::Paths MaskLayer::get_mask() const {
    return mask;
}

/**
 * Renders the layer to an SVG layer.
 */
svg::Layer MaskLayer::to_svg(const ColorScheme &colors, bool flipped, const std::string &id_prefix) const {
    auto layer = svg::Layer(id_prefix + get_name());
    if (bottom == flipped) {
        layer.add(mask, colors.soldermask);
        layer.add(silk, colors.silkscreen);
    } else {
        layer.add(silk, colors.silkscreen);
        layer.add(mask, colors.soldermask);
    }
    return layer;
}

/**
 * Renders the layer to an OBJ file.
 */
void MaskLayer::to_obj(obj::ObjFile &obj, size_t layer_index, double z, const std::string &id_prefix) const {
    double mask_z1, mask_z2, silk_z;
    std::string mask_name, silk_name;
    if (bottom) {
        mask_name = "_GBS";
        silk_name = "_GBO";
        silk_z = z;
    } else {
        mask_name = "_GTS";
        silk_name = "_GTO";
        silk_z = z + get_thickness();
    }

    // The soldermask is probably transparent, which makes Z-fighting a
    // potential issue. Just work around it by making the soldermask slightly
    // thinner than it should be. The same applies for the silkscreen layer,
    // which is modelled without thickness.
    mask_z1 = z + get_thickness() * 0.01;
    mask_z2 = z + get_thickness() * 0.99;

    obj.add_object(
        "layer" + std::to_string(layer_index) + mask_name,
        "soldermask"
    ).add_sheet(
        mask,
        mask_z1,
        mask_z2
    );
    obj.add_object(
        "layer" + std::to_string(layer_index) + silk_name,
        "silkscreen"
    ).add_surface(
        silk,
        silk_z
    );
}

/**
 * Returns an open file input stream for the given filename.
 */
std::ifstream CircuitBoard::read_file(const std::string &fname) {
    std::ifstream f(basename + fname);
    if (!f.is_open()) {
        throw std::runtime_error("file not found");
    }
    f.exceptions(std::ifstream::badbit);
    return f;
}

/**
 * Reads a Gerber file.
 */
coord::Paths CircuitBoard::read_gerber(const std::string &fname, bool outline) {
    if (fname.empty()) {
        return {};
    }
    //std::cout << "reading Gerber file " << fname << "..." << std::endl;
    auto f = read_file(fname);
    auto g = gerber::Gerber(f);
    auto paths = outline ? g.get_outline_paths() : g.get_paths();
    //std::cout << "finished reading " << fname << std::endl;
    return paths;
}

/**
 * Reads an NC drill file.
 */
void CircuitBoard::read_drill(const std::string &fname, bool plated, coord::Paths &pth, coord::Paths &npth) {
    if (fname.empty()) {
        return;
    }
    //std::cout << "reading drill file " << fname << "..." << std::endl;
    auto f = read_file(fname);
    auto d = ncdrill::NCDrill(f, plated);
    auto l = d.get_paths(true, false);
    if (pth.empty()) {
        pth = l;
    } else {
        pth.insert(pth.end(), l.begin(), l.end());
        ClipperLib::SimplifyPolygons(pth);
    }
    l = d.get_paths(false, true);
    if (npth.empty()) {
        npth = l;
    } else {
        npth.insert(npth.end(), l.begin(), l.end());
        ClipperLib::SimplifyPolygons(npth);
    }
    auto new_vias = d.get_vias();
    vias.insert(vias.end(), new_vias.begin(), new_vias.end());
}

/**
 * Constructs a circuit board. The following files are to be specified:
 *  - basename: prefix for all filenames.
 *  - outline: the board outline Gerber file (GKO, GM1, etc). May also
 *    include milling information.
 *  - drill: the NC drill file (TXT).
 *  - drill_nonplated: if specified, non-plated holes will be added from
 *    this auxiliary NC drill file.
 *  - mill: if specified, adds another Gerber-based milling layer.
 *    Interpreted in the same way as outline.
 *  - plating_thickness: thickness of the hole plating in millimeters.
 */
CircuitBoard::CircuitBoard(
    const std::string &basename,
    const std::string &outline,
    const std::string &drill,
    const std::string &drill_nonplated,
    const std::string &mill,
    double plating_thickness
) : basename(basename), num_substrate_layers(0), plating_thickness(coord::Format::from_mm(plating_thickness)) {

    // Load board outline.
    board_outline = read_gerber(outline, true);
    path::append(board_outline, read_gerber(mill, true));

    // Load drill files.
    coord::Paths pth, npth;
    read_drill(drill, true, pth, npth);
    if (!drill_nonplated.empty()) {
        read_drill(drill_nonplated, false, pth, npth);
    }

    // Make board shape.
    auto holes = path::add(pth, npth);
    board_shape = path::subtract(board_outline, holes);
    board_shape_excl_pth = path::subtract(board_outline, npth);

    // Build plating.
    coord::Paths pth_drill = path::offset(pth, this->plating_thickness, true);

    // Make substrate shape.
    substrate_dielectric = path::subtract(board_outline, path::add(pth_drill, npth));
    substrate_plating = path::subtract(pth_drill, pth);

}

/**
 * Adds a mask layer to the board. Layers are added bottom-up.
 */
void CircuitBoard::add_mask_layer(const std::string &mask, const std::string &silk) {
    layers.push_back(std::make_shared<MaskLayer>(
        "mask" + mask, board_outline, read_gerber(mask), read_gerber(silk), layers.empty()
    ));
}

/**
 * Adds a copper layer to the board. Layers are added bottom-up.
 */
void CircuitBoard::add_copper_layer(const std::string &gerber, double thickness) {
    layers.push_back(std::make_shared<CopperLayer>(
        "copper" + gerber, board_shape, board_shape_excl_pth, read_gerber(gerber), thickness
    ));
}

/**
 * Adds a substrate layer. Layers are added bottom-up.
 */
void CircuitBoard::add_substrate_layer(double thickness) {
    layers.push_back(std::make_shared<SubstrateLayer>(
        "substrate" + std::to_string(++num_substrate_layers), board_shape, substrate_dielectric, substrate_plating, thickness
    ));
}

/**
 * Derives the surface finish layer for all exposed copper.
 */
void CircuitBoard::add_surface_finish() {
    coord::Paths mask;
    for (auto it = layers.begin(); it != layers.end(); ++it) {
        auto copper = std::dynamic_pointer_cast<CopperLayer>(*it);
        if (copper) {
            bottom_finish = path::subtract(copper->get_copper(), mask);
            break;
        }
        mask = path::add(mask, (*it)->get_mask());
    }
    mask.clear();
    for (auto it = layers.rbegin(); it != layers.rend(); ++it) {
        auto copper = std::dynamic_pointer_cast<CopperLayer>(*it);
        if (copper) {
            top_finish = path::subtract(copper->get_copper(), mask);
            break;
        }
        mask = path::add(mask, (*it)->get_mask());
    }
}

/**
 * Returns a netlist builder initialized with the vias and copper regions of
 * this PCB.
 */
netlist::NetlistBuilder CircuitBoard::get_netlist_builder() const {
    netlist::NetlistBuilder nb;
    for (const auto &layer : layers) {
        auto copper = std::dynamic_pointer_cast<CopperLayer>(layer);
        if (copper) {
            nb.layer(copper->get_copper_excl_pth());
        }
    }
    for (const auto &via : vias) {
        nb.via(via.get_path(), via.get_finished_hole_size(), plating_thickness);
    }
    return nb;
}

/**
 * Returns the physical netlist for this PCB.
 */
netlist::PhysicalNetlist CircuitBoard::get_physical_netlist() const {
    netlist::PhysicalNetlist pn;
    size_t layer_index = 0;
    for (const auto &layer : layers) {
        auto copper = std::dynamic_pointer_cast<CopperLayer>(layer);
        if (copper) {
            pn.register_paths(copper->get_copper_excl_pth(), layer_index++);
        }
    }
    for (const auto &via : vias) {
        pn.register_via(
            std::make_shared<netlist::Via>(
                via.get_path(),
                via.get_finished_hole_size(),
                plating_thickness
            ),
            layer_index
        );
    }
    return pn;
}

/**
 * Returns the axis-aligned boundary coordinates of the PCB.
 */
coord::CRect CircuitBoard::get_bounds() const {
    coord::CRect bounds;
    bounds.left = bounds.bottom = INT64_MAX;
    bounds.right = bounds.top = INT64_MIN;
    for (const auto &path : board_outline) {
        for (const auto &point : path) {
            bounds.left = std::min(bounds.left, point.X);
            bounds.right = std::max(bounds.right, point.X);
            bounds.bottom = std::min(bounds.bottom, point.Y);
            bounds.top = std::max(bounds.top, point.Y);
        }
    }
    return bounds;
}

/**
 * Renders the circuit board to SVG, returning only the body of it, allowing it
 * to be composited into a larger image.
 */
std::string CircuitBoard::get_svg(bool flipped, const ColorScheme &colors, const std::string &id_prefix) const {
    std::ostringstream ss;

    if (flipped) {
        for (auto it = layers.rbegin(); it != layers.rend(); ++it) {
            ss << (*it)->to_svg(colors, flipped, id_prefix);
        }
    } else {
        for (auto it = layers.begin(); it != layers.end(); ++it) {
            ss << (*it)->to_svg(colors, flipped, id_prefix);
        }
    }

    auto finish = svg::Layer(id_prefix + "finish");
    finish.add(flipped ? bottom_finish : top_finish, colors.finish);
    ss << finish;

    return ss.str();
}

/**
 * Renders the circuit board to an SVG.
 */
void CircuitBoard::write_svg(
    const std::string &fname,
    bool flipped,
    double scale,
    const ColorScheme &colors
) const {
    auto bounds = get_bounds();

    auto width = bounds.right - bounds.left + coord::Format::from_mm(20.0);
    auto height = bounds.top - bounds.bottom + coord::Format::from_mm(20.0);
    svg::File svg{fname, {0, 0, width, height}, scale};

    std::ostringstream strm;
    auto tx = coord::Format::from_mm(10.0) - (flipped ? -bounds.right : bounds.left);
    auto ty = coord::Format::from_mm(10.0) + bounds.top;
    strm << "<g transform=\"";
    strm << "translate(" << coord::Format::to_mm(tx) << " " << coord::Format::to_mm(ty) << ") ";
    strm << "scale(" << (flipped ? "-1" : "1") << " -1) ";
    strm << "\" filter=\"drop-shadow(0 0 1 rgba(0, 0, 0, 0.2))\">\n";
    svg.add(strm.str());

    svg << get_svg(flipped, colors);

    svg << "</g>\n";
}

/**
 * Generates a path that approximates a circle of the given size with the given
 * center point.
 */
static void render_circle(coord::CPt center, coord::CInt diameter, coord::Path &output) {
    double epsilon = coord::Format().get_max_deviation();
    double r = diameter * 0.5;
    double x = (r > epsilon) ? (1.0 - epsilon / r) : 0.0;
    double th = std::acos(2.0 * x * x - 1.0) + 1e-3;
    auto n_vertices = (size_t)std::ceil(2.0 * M_PI / th);
    if (n_vertices < 3) n_vertices = 3;
    output.clear();
    output.reserve(n_vertices);
    for (size_t i = 0; i < n_vertices; i++) {
        auto a = 2.0 * M_PI * i / n_vertices;
        output.emplace_back(
            center.X + (coord::CInt)std::round(std::cos(a) * r),
            center.Y + (coord::CInt)std::round(std::sin(a) * r)
        );
    }
}

/**
 * Adds named copper shapes to the given Wavefront OBJ file manager.
 */
static void render_copper(
    obj::ObjFile &obj,
    const netlist::PhysicalNetlist &netlist,
    const std::vector<std::pair<double, double>> &copper_z
) {
    size_t name_counter = 1;
    for (const auto &net : netlist.get_nets()) {

        // Figure out a unique name for the copper object.
        std::string name;
        if (net->get_logical_nets().empty()) {
            name = "net_" + std::to_string(name_counter);
        } else {
            name = net->get_logical_nets().begin()->lock()->get_name() + "_" + std::to_string(name_counter);
        }
        name_counter++;
        auto &ob = obj.add_object(name, "copper");

        // Enumerate the vias connected to this net.
        struct Via {
            coord::CPt center;
            coord::Path inner;
            coord::Path outer;
            size_t lower_layer;
            size_t upper_layer;
        };
        std::vector<Via> vias;
        vias.reserve(net->get_vias().size());
        for (const auto &via : net->get_vias()) {
            auto center = via->get_coordinate();
            auto diameter = via->get_finished_hole_size();
            auto lower_layer = via->get_lower_layer(copper_z.size());
            auto upper_layer = via->get_upper_layer(copper_z.size());

            // Render the inner ring. This goes all the way from the bottom of
            // the lowest layer to the top of the upper layer.
            coord::Path inner;
            render_circle(center, diameter, inner);
            ob.add_ring(inner, copper_z.at(lower_layer).first, copper_z.at(upper_layer).second);

            // Render the outer rings. There is one for each layer that is
            // connected by the via, from the top of the lower of the two layers
            // to the bottom of the top of the two.
            coord::Path outer;
            render_circle(center, diameter + 2 * via->get_plating_thickness(), outer);
            for (size_t layer = lower_layer; layer < upper_layer; layer++) {
                ob.add_ring(outer, copper_z.at(layer).second, copper_z.at(layer + 1).first);
            }

            // Everything else consists of flat copper surfaces, along with
            // their own outline. We store the outer ring for this via in a
            // record to keep track of these, so we can punch holes in these
            // surfaces.
            vias.push_back({center, std::move(inner), std::move(outer), lower_layer, upper_layer});

        }

        // Enumerate and render the planar copper shapes connected to this net.
        for (const auto &shape : net->get_shapes()) {
            auto zs = copper_z.at(shape->get_layer());

            // Add the rings.
            ob.add_ring(shape->get_outline(), zs.first, zs.second);
            for (const auto &path : shape->get_holes()) {
                ob.add_ring(path, zs.first, zs.second);
            }

            // Add the surfaces.
            size_t layer = shape->get_layer();
            for (int side = 0; side < 2; side++) {

                // Side 0 is the bottom of the sheet, side 1 is the top.
                double z = side ? zs.second : zs.first;

                // Figure out the holes in this shape, including those made by
                // vias.
                coord::Paths holes = shape->get_holes();
                for (const auto &via : vias) {
                    if (layer < via.lower_layer) {
                        // Via is above this side/layer.
                        continue;
                    }
                    if (layer > via.upper_layer) {
                        // Via is below this side/layer.
                        continue;
                    }
                    if (!shape->contains(via.center)) {
                        // Via is in a different part of the PCB in-plane.
                        continue;
                    }
                    // Via punches through this side/layer, so we need to add an
                    // extra hole for it.
                    if ((layer == via.lower_layer && side == 0) || (layer == via.upper_layer && side == 1)) {
                        holes.push_back(via.inner);
                    } else {
                        holes.push_back(via.outer);
                    }
                }

                // Add the copper surface.
                ob.add_surface(shape->get_outline(), holes, z);

            }
        }
    }
}

/**
 * Renders the circuit board to a Wavefront OBJ file. Optionally, a netlist
 * can be supplied, of which the logical net names will then be used to
 * name the copper objects.
 */
void CircuitBoard::write_obj(const std::string &fname, const netlist::Netlist *netlist) const {
    obj::ObjFile obj;
    double z = 0.0;
    size_t index = 0;
    std::vector<std::pair<double, double>> copper_z;
    for (const auto &layer : layers) {
        layer->to_obj(obj, index++, z, "");
        if (std::dynamic_pointer_cast<CopperLayer>(layer)) {
            copper_z.emplace_back(z, z + layer->get_thickness());
        }
        z += layer->get_thickness();
    }
    if (netlist != nullptr) {
        render_copper(obj, netlist->get_physical_netlist(), copper_z);
    } else {
        render_copper(obj, get_physical_netlist(), copper_z);
    }
    obj.to_file(fname);
}

} // namespace pcb
} // namespace gerbertools
