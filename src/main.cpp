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

#include <iostream>
#include <fstream>
#include "gerber.hpp"
#include "ncdrill.hpp"

/**
 * Append paths odd-even style.
 */
static void path_append(plot::Paths &dest, const plot::Paths &src) {
    if (src.empty()) {
        return;
    }
    if (dest.empty()) {
        dest = src;
        return;
    }
    dest.insert(dest.end(), src.begin(), src.end());
    ClipperLib::SimplifyPolygons(dest);
}

/**
 * Perform an operation between two sets of paths.
 */
static plot::Paths path_op(const plot::Paths &lhs, const plot::Paths &rhs, ClipperLib::ClipType op) {
    ClipperLib::Clipper c;
    c.AddPaths(lhs, ClipperLib::ptSubject, true);
    c.AddPaths(rhs, ClipperLib::ptClip, true);
    plot::Paths result;
    c.Execute(op, result);
    return result;
}

/**
 * Compute the union of two sets of paths.
 */
static plot::Paths path_union(const plot::Paths &lhs, const plot::Paths &rhs) {
    return path_op(lhs, rhs, ClipperLib::ctUnion);
}

/**
 * Compute the difference between two sets of paths.
 */
static plot::Paths path_subtract(const plot::Paths &lhs, const plot::Paths &rhs) {
    return path_op(lhs, rhs, ClipperLib::ctDifference);
}

/**
 * Compute the intersection between two sets of paths.
 */
static plot::Paths path_intersect(const plot::Paths &lhs, const plot::Paths &rhs) {
    return path_op(lhs, rhs, ClipperLib::ctIntersection);
}

/**
 * Offsets the paths by the given amount.
 */
static plot::Paths path_offset(const plot::Paths &src, double amount) {
    auto co = coord::Format().build_clipper_offset();
    co.AddPaths(src, ClipperLib::jtMiter, ClipperLib::etOpenRound);
    plot::Paths result;
    co.Execute(result, coord::Format::from_mm(amount));
    ClipperLib::SimplifyPolygons(result, ClipperLib::pftPositive);
    return result;
}

struct Color {
    float r;
    float g;
    float b;
    float a;
};

static const Color COLOR_TIN = {0.7f, 0.7f, 0.7f, 1.0f};

/*
static const Color COLOR_COPPER = {0.8f, 0.59f, 0.3f, 1.0f};
static const Color COLOR_SUBSTRATE = {0.8f, 0.7f, 0.5f, 0.95f};

// Green soldermask.
static const Color COLOR_MASK = {0.1f, 0.6f, 0.3f, 0.6f};
static const Color COLOR_SILK = {0.9f, 0.9f, 0.9f, 0.9f};
*/

// White soldermask.
//static const Color COLOR_MASK = {0.85f, 0.9f, 0.95f, 0.8f};
//static const Color COLOR_SILK = {0.1f, 0.1f, 0.1f, 0.9f};

static const Color COLOR_COPPER    = {0.8f, 0.7f, 0.3f, 1.0f};
static const Color COLOR_SUBSTRATE = {0.6f, 0.5f, 0.3f, 0.95f};
static const Color COLOR_MASK      = {0.1f, 0.6f, 0.3f, 0.6f};
static const Color COLOR_SILK      = {0.9f, 0.9f, 0.9f, 0.9f};

class Svg {
private:
    std::ofstream f;
    coord::CRect bounds;
    double scale;
public:
    Svg(
        const std::string &fname,
        const coord::CRect &bounds,
        double scale = 1.0
    ) : f(fname), bounds(bounds), scale(scale) {
        if (!f.is_open()) {
            throw std::runtime_error("failed to open " + fname + " for writing");
        }
        f << "<svg width=\"" << coord::Format::to_mm(bounds.right - bounds.left) * scale;
        f << "\" height=\"" << coord::Format::to_mm(bounds.top - bounds.bottom) * scale;
        f << "\" xmlns=\"http://www.w3.org/2000/svg\">)\n)";
        f << R"(<defs><filter id="f1"><feGaussianBlur in="SourceGraphic" stdDeviation=")" << scale << R"(" /></filter></defs>)" << "\n";
    }

    void draw(plot::Paths paths, bool flipped, Color c = {0.0f, 0.0f, 0.0f, 1.0f}, bool blurred=false) {
        f << "<path fill=\"rgb(";
        f << (int)(c.r*255) << ",";
        f << (int)(c.g*255) << ",";
        f << (int)(c.b*255);
        f << ")\" fill-opacity=\"" << c.a << "\" ";
        if (blurred) {
            f << "filter=\"url(#f1)\" ";
        }
        f << "d=\"";
        for (const auto &p : paths) {
            if (flipped) {
                f << "M " << coord::Format::to_mm(bounds.right - p.back().X) * scale;
            } else {
                f << "M " << coord::Format::to_mm(p.back().X - bounds.left) * scale;
            }
            f << " " << coord::Format::to_mm(bounds.top - p.back().Y) * scale;
            f << " ";
            for (const auto &c : p) {
                if (flipped) {
                    f << "L " << coord::Format::to_mm(bounds.right - c.X) * scale;
                } else {
                    f << "L " << coord::Format::to_mm(c.X - bounds.left) * scale;
                }
                f << " " << coord::Format::to_mm(bounds.top - c.Y) * scale;
                f << " ";
            }
        }
        f << "\"/>\n";
    }

    void close() {
        if (f.is_open()) {
            f << R"(</svg>)" << std::endl;
            f.close();
        }
    }

    ~Svg() {
        close();
    }
};

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
    explicit Layer(double thickness) : thickness(thickness) {
    }

public:

    /**
     * Destroys a layer.
     * FIXME: forgot this for the other polymorphic classes!
     */
    virtual ~Layer() = default;

    /**
     * Returns the thickness of this layer.
     */
    double get_thickness() const {
        return thickness;
    }

    /**
     * Returns whether this is a substrate layer.
     */
    virtual bool is_substrate() const {
        return false;
    }

    /**
     * Renders the layer to an SVG.
     */
    virtual void render(Svg &svg, bool flipped) const {
    }

};

class SubstrateLayer : public Layer {
public:

    /**
     * Constructs a substate layer.
     */
    explicit SubstrateLayer(double thickness) : Layer(thickness) {
    }

    /**
     * Returns whether this is a substrate layer (it is).
     */
    bool is_substrate() const override {
        return true;
    }

};

class CopperLayer : public Layer {
private:

    /**
     * Shape of the copper as specified in the Gerber file.
     */
    plot::Paths layer;

    /**
     * Actual shape of the copper. That is, layer minus board outline.
     */
    plot::Paths copper;

public:

    /**
     * Constructs a copper layer.
     */
    CopperLayer(
        const plot::Paths &board_shape,
        const plot::Paths &copper_layer,
        double thickness
    ) : Layer(thickness) {
        layer = copper_layer;
        copper = path_intersect(board_shape, copper_layer);
    }

    /**
     * Renders the layer to an SVG.
     */
    void render(Svg &svg, bool flipped) const override {
        svg.draw(copper, flipped, COLOR_COPPER);
    }

};

class MaskLayer : public Layer {
private:

    /**
     * Shape of the mask. Intersection of the solder mask Gerber file and the
     * board outline.
     */
    plot::Paths mask;

    /**
     * Shape of the silkscreen. Intersection of the above solder mask shape and
     * the silkscreen layer.
     */
    plot::Paths silk;

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
        const plot::Paths &board_outline,
        const plot::Paths &mask_layer,
        const plot::Paths &silk_layer,
        bool bottom
    ) : Layer(0.01) {
        mask = path_subtract(board_outline, mask_layer);
        silk = path_intersect(mask, silk_layer);
        bottom = bottom;
    }

    /**
     * Renders the layer to an SVG.
     */
    void render(Svg &svg, bool flipped) const override {
        if (bottom == flipped) {
            svg.draw(silk, flipped, COLOR_SILK);
            svg.draw(mask, flipped, COLOR_MASK);
        } else {
            svg.draw(mask, flipped, COLOR_MASK);
            svg.draw(silk, flipped, COLOR_SILK);
        }
    }

};

/**
 * Reference to a layer.
 */
using LayerRef = std::shared_ptr<Layer>;

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
    plot::Paths board_outline;

    /**
     * The board shape, minus holes after plating.
     */
    plot::Paths board_shape;

    /**
     * The board outline, minus holes before plating.
     */
    plot::Paths substrate_dielectric;

    /**
     * Substrate plating; i.e. copper on the sides of the substrate dielectric.
     */
    plot::Paths substrate_plating;

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
    std::ifstream read_file(const std::string &fname) {
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
    plot::Paths read_gerber(const std::string &fname, bool outline=false) {
        if (fname.empty()) {
            return {};
        }
        std::cout << "reading Gerber file " << fname << "..." << std::endl;
        auto f = read_file(fname);
        auto g = gerber::Gerber(f);
        auto paths = outline ? g.get_outline_paths() : g.get_paths();
        std::cout << "finished reading " << fname << std::endl;
        return paths;
    }

    /**
     * Reads an NC drill file.
     */
    void read_drill(const std::string &fname, bool plated, plot::Paths &pth, plot::Paths &npth) {
        if (fname.empty()) {
            return;
        }
        std::cout << "reading drill file " << fname << "..." << std::endl;
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
        double plating_thickness = 0.1 // TODO look up sane default
    ) : basename(basename) {

        // Load board outline.
        board_outline = read_gerber(outline, true);
        path_append(board_outline, read_gerber(mill, true));

        // Load drill files.
        plot::Paths pth, npth;
        read_drill(drill, true, pth, npth);
        if (!drill_nonplated.empty()) {
            read_drill(drill, false, pth, npth);
        }

        // Make board shape.
        auto holes = path_union(pth, npth);
        board_shape = path_subtract(board_outline, holes);

        // Build plating.
        plot::Paths pth_drill = path_offset(pth, plating_thickness);

        // Make substrate shape.
        substrate_dielectric = path_subtract(board_outline, path_union(pth_drill, npth));
        substrate_plating = path_subtract(pth_drill, pth);

    }

    /**
     * Adds a mask layer to the board. Layers are added bottom-up.
     */
    void add_mask_layer(
        const std::string &mask,
        const std::string &silk=""
    ) {
        layers.push_back(std::make_shared<MaskLayer>(
            board_outline, read_gerber(mask), read_gerber(silk), layers.empty()
        ));
    }

    /**
     * Adds a copper layer to the board. Layers are added bottom-up.
     */
    void add_copper_layer(
        const std::string &gerber,
        double thickness = 0.04 // TODO look up sane default
    ) {
        layers.push_back(std::make_shared<CopperLayer>(
            board_shape, read_gerber(gerber), thickness
        ));
    }

    /**
     * Adds a substrate layer.
     */
    void add_substrate_layer(double thickness = 1.5) {
        layers.push_back(std::make_shared<SubstrateLayer>(thickness));
    }

    void render_layer(Svg &svg, const LayerRef &layer, bool flipped) {
        if (!layer->is_substrate()) {
            layer->render(svg, flipped);
            return;
        }
        svg.draw(board_shape, flipped, COLOR_SUBSTRATE);
    }

    void write_svg(const std::string &fname, bool flipped=false, double scale=1.0) {
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
        bounds.right += coord::Format::from_mm(10.0);
        bounds.bottom -= coord::Format::from_mm(10.0);
        bounds.top += coord::Format::from_mm(10.0);
        Svg svg{fname, bounds, scale};
        svg.draw(board_shape, flipped, {0.0f, 0.0f, 0.0f, 0.2f}, true);
        if (!flipped) {
            for (auto it = layers.begin(); it != layers.end(); ++it) {
                render_layer(svg, *it, flipped);
            }
        } else {
            for (auto it = layers.rbegin(); it != layers.rend(); ++it) {
                render_layer(svg, *it, flipped);
            }
        }
    }

};

int main(int argc, char *argv[]) {

    CircuitBoard pcb(
        "/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/fts/fts", ".GKO", ".TXT"
        //"O100030117 10by10 Green 1.6mm HASL 10pcs/mcu", ".GKO", ".TXT"
        //"pcb", ".GM3", ""
    );
    pcb.add_mask_layer(".GBS", ".GBO");
    pcb.add_copper_layer(".GBL");
    pcb.add_substrate_layer();
    pcb.add_copper_layer(".GTL");
    pcb.add_mask_layer(".GTS", ".GTO");
    pcb.write_svg("kek.svg", true, 20);

    /*std::ifstream f("/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/fts/fts.GBL");
    //std::ifstream f("O100030117 10by10 Green 1.6mm HASL 10pcs/mcu.GTL");
    if (!f.is_open()) {
        throw std::runtime_error("file not found");
    }
    f.exceptions(std::ifstream::badbit);
    auto gerber = gerber::Gerber(f);

    std::ifstream f2("/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/fts/fts.TXT");
    //std::ifstream f2("O100030117 10by10 Green 1.6mm HASL 10pcs/mcu.TXT");
    if (!f.is_open()) {
        throw std::runtime_error("file not found");
    }
    f2.exceptions(std::ifstream::badbit);
    auto drill = ncdrill::NCDrill(f2);

    ClipperLib::Clipper c;
    c.StrictlySimple(true);
    c.AddPaths(gerber.get_paths(), ClipperLib::ptSubject, true);
    c.AddPaths(drill.get_paths(), ClipperLib::ptClip, true);
    plot::Paths paths;
    c.Execute(ClipperLib::ctDifference, paths);

    std::ofstream svg("kek.svg");
    svg << R"(<svg width="10000" height="10000" xmlns="http://www.w3.org/2000/svg">))" << std::endl;
    int color = 0;
    for (const auto &p : paths) {
        if (ClipperLib::Orientation(p)) {
            color++;
            int r = 128 + 127*std::sin(color);
            int g = 128 + 127*std::sin(color + 2*M_PI / 3);
            int b = 128 + 127*std::sin(color + 4*M_PI / 3);
            svg << "<path fill=\"rgb("<<r<<","<<g<<","<<b<<")\" fill-opacity=\"0.7\" d=\"";
        } else {
            svg << R"(<path fill="white" fill-opacity="0.7" d=")";
        }
        svg << "M " << coord::Format::to_mm(p.back().X) * 25 + 5000
            << " " << coord::Format::to_mm(-p.back().Y) * 25 + 5000 << " ";
        for (const auto &c : p) {
            svg << "L " << coord::Format::to_mm(c.X) * 25 + 5000
                << " " << coord::Format::to_mm(-c.Y) * 25 + 5000 << " ";
        }
        svg << R"("/>)" << std::endl;
    }
    svg << R"(</svg>)" << std::endl;*/

}
