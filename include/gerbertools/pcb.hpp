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
 * Helper class for rendering a PCB to an SVG.
 */
class Svg {
private:

    /**
     * Output stream for the SVG.
     */
    std::ofstream f;

    /**
     * Boundbox of the PCB in millimeters.
     */
    coord::CRect bounds;

    /**
     * Number of SVG units per millimeter.
     */
    double scale;

    /**
     * Size of the SVG in SVG units.
     */
    double width, height;

public:

    /**
     * Starts rendering an SVG with the given filename, PCB bounds in
     * millimeter, and SVG units per millimeter.
     */
    Svg(const std::string &fname, const coord::CRect &bounds, double scale=1.0);

    /**
     * Destroys this SVG writer, finishing the SVG first.
     */
    ~Svg();

    /**
     * Draws a path to the SVG. If the color is fully clear, the path is
     * omitted.
     */
    void draw(coord::Paths paths, bool flipped, color::Color c = color::BLACK, bool blurred=false);

    /**
     * Finishes writing the SVG.
     */
    void close();

};

/**
 * Represents any PCB layer type.
 */
class Layer {
private:

    /**
     * Thickness of the layer.
     */
    double thickness;

protected:

    /**
     * Constructs a layer.
     */
    explicit Layer(double thickness);

public:
    virtual ~Layer() = default;

    /**
     * Returns the thickness of this layer.
     */
    double get_thickness() const;

    /**
     * Returns the surface finish mask for this layer.
     */
    virtual coord::Paths get_mask() const = 0;

    /**
     * Renders the layer to an SVG.
     */
    virtual void render(Svg &svg, bool flipped, const ColorScheme &colors) const = 0;

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
     * Renders the layer to an SVG.
     */
    void render(Svg &svg, bool flipped, const ColorScheme &colors) const override;

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
     * Actual shape of the copper. That is, layer minus board outline.
     */
    coord::Paths copper;

public:

    /**
     * Constructs a copper layer.
     */
    CopperLayer(
        const coord::Paths &board_shape,
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
     * Renders the layer to an SVG.
     */
    void render(Svg &svg, bool flipped, const ColorScheme &colors) const override;

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
     * Renders the layer to an SVG.
     */
    void render(Svg &svg, bool flipped, const ColorScheme &colors) const override;

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
    std::list<coord::CPt> via_points;

    /**
     * Layers that constitute the board.
     */
    std::list<LayerRef> layers;

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
     * Renders the circuit board to an SVG.
     */
    void write_svg(const std::string &fname, double scale=1.0, const ColorScheme &colors={});

};

} // namespace pcb
} // namespace gerbertools
