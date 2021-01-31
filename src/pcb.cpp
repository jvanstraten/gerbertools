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

#include "gerbertools/gerber.hpp"
#include "gerbertools/ncdrill.hpp"
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
 * Starts rendering an SVG with the given filename, PCB bounds in
 * millimeter, and SVG units per millimeter.
 */
Svg::Svg(
    const std::string &fname,
    const coord::CRect &bounds,
    double scale
) : f(fname), bounds(bounds), scale(scale) {
    if (!f.is_open()) {
        throw std::runtime_error("failed to open " + fname + " for writing");
    }
    width = coord::Format::to_mm(2 * (bounds.right - bounds.left)) * scale;
    height = coord::Format::to_mm(bounds.top - bounds.bottom) * scale;
    f << "<svg width=\"" << width << "\" height=\"" << height << "\" ";
    f << "xmlns=\"http://www.w3.org/2000/svg\">)\n)";
    f << R"(<defs><filter id="f1"><feGaussianBlur in="SourceGraphic" stdDeviation=")" << scale << R"(" /></filter></defs>)" << "\n";
    /*f << R"(<path fill="mediumturquoise" d=")";
    f << "M " << 0 << " " << 0 << " ";
    f << "L " << width << " " << 0 << " ";
    f << "L " << width << " " << height << " ";
    f << "L " << 0 << " " << height << " ";
    f << R"(" />)" << "\n";*/
}

/**
 * Destroys this SVG writer, finishing the SVG first.
 */
Svg::~Svg() {
    close();
}

/**
 * Draws a path to the SVG. If the color is fully clear, the path is
 * omitted.
 */
void Svg::draw(coord::Paths paths, bool flipped, color::Color color, bool blurred) {
    if (color.a == 0.0) return;

    f << "<path fill=\"rgb(";
    f << (int)(color.r*255) << ",";
    f << (int)(color.g*255) << ",";
    f << (int)(color.b*255);
    f << ")\" fill-opacity=\"" << color.a << "\" ";
    if (blurred) {
        f << "filter=\"url(#f1)\" ";
    }
    f << "d=\"";
    for (const auto &p : paths) {
        if (flipped) {
            f << "M " << width - coord::Format::to_mm(p.back().X - bounds.left) * scale;
        } else {
            f << "M " << coord::Format::to_mm(p.back().X - bounds.left) * scale;
        }
        f << " " << coord::Format::to_mm(bounds.top - p.back().Y) * scale;
        f << " ";
        for (const auto &c : p) {
            if (flipped) {
                f << "L " << width - coord::Format::to_mm(c.X - bounds.left) * scale;
            } else {
                f << "L " << coord::Format::to_mm(c.X - bounds.left) * scale;
            }
            f << " " << coord::Format::to_mm(bounds.top - c.Y) * scale;
            f << " ";
        }
    }
    f << "\"/>\n";
}

/**
 * Finishes writing the SVG.
 */
void Svg::close() {
    if (f.is_open()) {
        f << R"(</svg>)" << std::endl;
        f.close();
    }
}

/**
 * Constructs a layer.
 */
Layer::Layer(double thickness) : thickness(thickness) {
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
    const coord::Paths &shape,
    const coord::Paths &dielectric,
    const coord::Paths &plating,
    double thickness
) : Layer(thickness), shape(shape), dielectric(dielectric), plating(plating) {
}

/**
 * Returns the surface finish mask for this layer.
 */
coord::Paths SubstrateLayer::get_mask() const {
    return dielectric;
}

/**
 * Renders the layer to an SVG.
 */
void SubstrateLayer::render(Svg &svg, bool flipped, const ColorScheme &colors) const {
    svg.draw(dielectric, flipped, colors.substrate);
    svg.draw(plating, flipped, colors.finish);
}

/**
 * Constructs a copper layer.
 */
CopperLayer::CopperLayer(
    const coord::Paths &board_shape,
    const coord::Paths &copper_layer,
    double thickness
) : Layer(thickness) {
    layer = copper_layer;
    copper = path::intersect(board_shape, copper_layer);
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
 * Renders the layer to an SVG.
 */
void CopperLayer::render(Svg &svg, bool flipped, const ColorScheme &colors) const {
    svg.draw(copper, flipped, colors.copper);
}

/**
 * Constructs a solder mask.
 */
MaskLayer::MaskLayer(
    const coord::Paths &board_outline,
    const coord::Paths &mask_layer,
    const coord::Paths &silk_layer,
    bool bottom
) : Layer(0.01) {
    mask = path::subtract(board_outline, mask_layer);
    silk = path::intersect(mask, silk_layer);
    bottom = bottom;
}

/**
 * Returns the surface finish mask for this layer.
 */
coord::Paths MaskLayer::get_mask() const {
    return mask;
}

/**
 * Renders the layer to an SVG.
 */
void MaskLayer::render(Svg &svg, bool flipped, const ColorScheme &colors) const {
    if (bottom == flipped) {
        svg.draw(silk, flipped, colors.silkscreen);
        svg.draw(mask, flipped, colors.soldermask);
    } else {
        svg.draw(mask, flipped, colors.soldermask);
        svg.draw(silk, flipped, colors.silkscreen);
    }
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
        pth.insert(pth.end(), l.begin(), l.end());
        ClipperLib::SimplifyPolygons(pth);
    }
    auto vias = d.get_vias();
    via_points.insert(via_points.end(), vias.begin(), vias.end());
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
) : basename(basename) {

    // Load board outline.
    board_outline = read_gerber(outline, true);
    path::append(board_outline, read_gerber(mill, true));

    // Load drill files.
    coord::Paths pth, npth;
    read_drill(drill, true, pth, npth);
    if (!drill_nonplated.empty()) {
        read_drill(drill, false, pth, npth);
    }

    // Make board shape.
    auto holes = path::add(pth, npth);
    board_shape = path::subtract(board_outline, holes);

    // Build plating.
    coord::Paths pth_drill = path::offset(pth, plating_thickness);

    // Make substrate shape.
    substrate_dielectric = path::subtract(board_outline, path::add(pth_drill, npth));
    substrate_plating = path::subtract(pth_drill, pth);

}

/**
 * Adds a mask layer to the board. Layers are added bottom-up.
 */
void CircuitBoard::add_mask_layer(const std::string &mask, const std::string &silk) {
    layers.push_back(std::make_shared<MaskLayer>(
        board_outline, read_gerber(mask), read_gerber(silk), layers.empty()
    ));
}

/**
 * Adds a copper layer to the board. Layers are added bottom-up.
 */
void CircuitBoard::add_copper_layer(const std::string &gerber, double thickness) {
    layers.push_back(std::make_shared<CopperLayer>(
        board_shape, read_gerber(gerber), thickness
    ));
}

/**
 * Adds a substrate layer. Layers are added bottom-up.
 */
void CircuitBoard::add_substrate_layer(double thickness) {
    layers.push_back(std::make_shared<SubstrateLayer>(
        board_shape, substrate_dielectric, substrate_plating, thickness
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
 * Renders the circuit board to an SVG.
 */
void CircuitBoard::write_svg(
    const std::string &fname,
    double scale,
    const ColorScheme &colors
) {
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
    bounds.left -= coord::Format::from_mm(10.0);
    bounds.right += coord::Format::from_mm(5.0);
    bounds.bottom -= coord::Format::from_mm(10.0);
    bounds.top += coord::Format::from_mm(10.0);
    Svg svg{fname, bounds, scale};
    svg.draw(board_shape, false, {0.0f, 0.0f, 0.0f, 0.2f}, true);
    for (auto it = layers.begin(); it != layers.end(); ++it) {
        (*it)->render(svg, false, colors);
    }
    svg.draw(top_finish, false, colors.finish, false);
    svg.draw(board_shape, true, {0.0f, 0.0f, 0.0f, 0.2f}, true);
    for (auto it = layers.rbegin(); it != layers.rend(); ++it) {
        (*it)->render(svg, true, colors);
    }
    svg.draw(bottom_finish, true, colors.finish, false);
}

} // namespace pcb
} // namespace gerbertools
