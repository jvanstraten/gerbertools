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

#pragma once

#include <string>
#include <memory>
#include <fstream>
#include "gerbertools/coord.hpp"
#include "gerbertools/color.hpp"
#include "gerbertools/svg.hpp"
#include "gerbertools/obj.hpp"
#include "gerbertools/netlist.hpp"
#include "gerbertools/ncdrill.hpp"

namespace gerbertools {

/**
 * Contains the CircuitBoard class, representing a complete circuit board as
 * reconstructed from Gerber files and an NC drill file.
 */
namespace pcb {

/**
 * Copper thickness in millimeters per ounce.
 */
static const double COPPER_OZ = 0.0348;

/**
 * A PCB color scheme.
 */
struct ColorScheme {
public:

    /**
     * Soldermask color.
     */
    color::Color soldermask;

    /**
     * Silkscreen color.
     */
    color::Color silkscreen;

    /**
     * Surface-finished copper color.
     */
    color::Color finish;

    /**
     * Substrate color.
     */
    color::Color substrate;

    /**
     * Unexposed copper color.
     */
    color::Color copper;

    /**
     * Creates a new color scheme.
     */
    ColorScheme(
        const color::Color &soldermask=color::MASK_GREEN,
        const color::Color &silkscreen=color::SILK_WHITE,
        const color::Color &finish=color::FINISH_TIN,
        const color::Color &substrate=color::SUBSTRATE,
        const color::Color &copper=color::COPPER
    );

};

/**
 * Represents any PCB layer type.
 */
class Layer {
private:

    /**
     * The layer name.
     */
    std::string name;

    /**
     * Thickness of the layer.
     */
    double thickness;

protected:

    /**
     * Constructs a layer.
     */
    explicit Layer(const std::string &name, double thickness);

public:
    virtual ~Layer() = default;

    /**
     * Returns a name for the layer.
     */
    std::string get_name() const;

    /**
     * Returns the thickness of this layer.
     */
    double get_thickness() const;

    /**
     * Returns the surface finish mask for this layer.
     */
    virtual coord::Paths get_mask() const = 0;

    /**
     * Renders the layer to an SVG layer.
     */
    virtual svg::Layer to_svg(const ColorScheme &colors, bool flipped, const std::string &id_prefix) const = 0;

    /**
     * Renders the layer to an OBJ file.
     */
    virtual void to_obj(obj::ObjFile &obj, size_t layer_index, double z, const std::string &id_prefix) const = 0;

};

/**
 * Reference to a layer.
 */
using LayerRef = std::shared_ptr<Layer>;

/**
 * Represents a dielectric substrate layer.
 */
class SubstrateLayer : public Layer {
private:

    /**
     * The board shape, minus holes after plating.
     */
    coord::Paths shape;

    /**
     * The board outline, minus holes before plating.
     */
    coord::Paths dielectric;

    /**
     * Substrate plating; i.e. copper on the sides of the substrate dielectric.
     */
    coord::Paths plating;

public:

    /**
     * Constructs a substrate layer.
     */
    explicit SubstrateLayer(
        const std::string &name,
        const coord::Paths &shape,
        const coord::Paths &dielectric,
        const coord::Paths &plating,
        double thickness
    );

    /**
     * Returns the surface finish mask for this layer.
     */
    coord::Paths get_mask() const override;

    /**
     * Renders the layer to an SVG layer.
     */
    svg::Layer to_svg(const ColorScheme &colors, bool flipped, const std::string &id_prefix) const override;

    /**
     * Renders the layer to an OBJ file.
     */
    void to_obj(obj::ObjFile &obj, size_t layer_index, double z, const std::string &id_prefix) const override;

};

/**
 * Represents a copper layer.
 */
class CopperLayer : public Layer {
private:

    /**
     * Shape of the copper as specified in the Gerber file.
     */
    coord::Paths layer;

    /**
     * Actual shape of the copper. That is, layer minus board outline and all
     * finished holes.
     */
    coord::Paths copper;

    /**
     * As above, but without cutouts for plated holes.
     */
    coord::Paths copper_excl_pth;

public:

    /**
     * Constructs a copper layer.
     */
    CopperLayer(
        const std::string &name,
        const coord::Paths &board_shape,
        const coord::Paths &board_shape_excl_pth,
        const coord::Paths &copper_layer,
        double thickness
    );

    /**
     * Returns the surface finish mask for this layer.
     */
    coord::Paths get_mask() const override;

    /**
     * Returns the copper for this layer.
     */
    const coord::Paths &get_copper() const;

    /**
     * Returns the copper for this layer, with only cutouts for board shape and
     * non-plated holes, not for plated holes. This is needed for annular ring DRC.
     */
    const coord::Paths &get_copper_excl_pth() const;

    /**
     * Returns the original layer, without board outline intersection.
     */
    const coord::Paths &get_layer() const;

    /**
     * Renders the layer to an SVG layer.
     */
    svg::Layer to_svg(const ColorScheme &colors, bool flipped, const std::string &id_prefix) const override;

    /**
     * Renders the layer to an OBJ file.
     */
    void to_obj(obj::ObjFile &obj, size_t layer_index, double z, const std::string &id_prefix) const override;

};

/**
 * Represents a soldermask layer.
 */
class MaskLayer : public Layer {
private:

    /**
     * Shape of the mask. Intersection of the solder mask Gerber file and the
     * board outline.
     */
    coord::Paths mask;

    /**
     * Shape of the silkscreen. Intersection of the above solder mask shape and
     * the silkscreen layer.
     */
    coord::Paths silk;

    /**
     * Whether this mask layer is at the bottom (true) or top (false).
     * Determines draw order of the two sublayers.
     */
    bool bottom;

public:

    /**
     * Constructs a solder mask.
     */
    MaskLayer(
        const std::string &name,
        const coord::Paths &board_outline,
        const coord::Paths &mask_layer,
        const coord::Paths &silk_layer,
        bool bottom
    );

    /**
     * Returns the surface finish mask for this layer.
     */
    coord::Paths get_mask() const override;

    /**
     * Renders the layer to an SVG layer.
     */
    svg::Layer to_svg(const ColorScheme &colors, bool flipped, const std::string &id_prefix) const override;

    /**
     * Renders the layer to an OBJ file.
     */
    void to_obj(obj::ObjFile &obj, size_t layer_index, double z, const std::string &id_prefix) const override;

};

/**
 * Represents a circuit board.
 */
class CircuitBoard {
private:

    /**
     * Prefix for all filenames.
     */
    std::string basename;

    /**
     * The board outline, without removal of holes.
     */
    coord::Paths board_outline;

    /**
     * The board shape, minus holes after plating.
     */
    coord::Paths board_shape;

    /**
     * The board shape, with cutouts for non-plated holes, but not for plated
     * holes/vias.
     */
    coord::Paths board_shape_excl_pth;

    /**
     * The board outline, minus holes before plating.
     */
    coord::Paths substrate_dielectric;

    /**
     * Substrate plating; i.e. copper on the sides of the substrate dielectric.
     */
    coord::Paths substrate_plating;

    /**
     * Copper surface finish on the bottom of the PCB.
     */
    coord::Paths bottom_finish;

    /**
     * Copper surface finish on the top of the PCB.
     */
    coord::Paths top_finish;

    /**
     * Points representing vias.
     */
    std::list<ncdrill::Via> vias;

    /**
     * Plating thickness for vias.
     */
    coord::CInt plating_thickness;

    /**
     * Layers that constitute the board.
     */
    std::list<LayerRef> layers;

    /**
     * The total number of substrate layers added thus far, used for generating
     * unique layer IDs.
     */
    size_t num_substrate_layers;

    /**
     * Returns an open file input stream for the given filename.
     */
    std::ifstream read_file(const std::string &fname);

    /**
     * Reads a Gerber file.
     */
    coord::Paths read_gerber(const std::string &fname, bool outline=false);

    /**
     * Reads an NC drill file.
     */
    void read_drill(const std::string &fname, bool plated, coord::Paths &pth, coord::Paths &npth);

public:

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
    CircuitBoard(
        const std::string &basename,
        const std::string &outline,
        const std::string &drill,
        const std::string &drill_nonplated = "",
        const std::string &mill = "",
        double plating_thickness = 0.5 * COPPER_OZ
    );

    /**
     * Adds a mask layer to the board. Layers are added bottom-up.
     */
    void add_mask_layer(const std::string &mask, const std::string &silk="");

    /**
     * Adds a copper layer to the board. Layers are added bottom-up.
     */
    void add_copper_layer(const std::string &gerber, double thickness=COPPER_OZ);

    /**
     * Adds a substrate layer. Layers are added bottom-up.
     */
    void add_substrate_layer(double thickness=1.5);

    /**
     * Derives the surface finish layer for all exposed copper.
     */
    void add_surface_finish();

    /**
     * Returns a netlist builder initialized with the vias and copper regions of
     * this PCB.
     */
    netlist::NetlistBuilder get_netlist_builder() const;

    /**
     * Returns the physical netlist for this PCB.
     */
    netlist::PhysicalNetlist get_physical_netlist() const;

    /**
     * Returns the axis-aligned boundary coordinates of the PCB.
     */
    coord::CRect get_bounds() const;

    /**
     * Renders the circuit board to SVG, returning only the body of it, allowing it
     * to be composited into a larger image.
     */
    std::string get_svg(bool flipped, const ColorScheme &colors, const std::string &id_prefix="") const;

    /**
     * Renders the circuit board to an SVG.
     */
    void write_svg(const std::string &fname, bool flipped=false, double scale=1.0, const ColorScheme &colors={}) const;

    /**
     * Renders the circuit board to a Wavefront OBJ file. Optionally, a netlist
     * can be supplied, of which the logical net names will then be used to
     * name the copper objects.
     */
    void write_obj(const std::string &fname, const netlist::Netlist *netlist = nullptr) const;

};

} // namespace pcb
} // namespace gerbertools
