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
 * Contains tools for rendering SVG files.
 */

#include <gerbertools/svg.hpp>

namespace gerbertools {
namespace svg {

/**
 * Build function for adding a key/value pair.
 */
Attributes Attributes::with(const std::string &key, const std::string &value) && {
    data << " " << key << "=\"" << value << "\"";
    return std::move(*this);
}

/**
 * Writes the data to a stream.
 */
std::ostream &operator<<(std::ostream &os, const Attributes &attr) {
    os << attr.data.str();
    return os;
}

/**
 * Constructs a layer with the given SVG identifier and an optional set of
 * attributes for the group.
 */
Layer::Layer(const std::string &identifier, const Attributes &attr) {
    data << "<g id=\"" << identifier << "\"" << attr << ">\n";
}

/**
 * Adds a filled path to the layer. If the color is fully clear, the path is
 * omitted. attr may be used to add extra attributes to the path.
 */
void Layer::add(coord::Paths paths, color::Color color, const Attributes &attr) {
    if (color.a == 0.0) return;
    data << "<path fill=\"rgb(";
    data << (int)(color.r*255) << ",";
    data << (int)(color.g*255) << ",";
    data << (int)(color.b*255);
    data << ")\"";
    if (color.a < 1.0) {
        data << " fill-opacity=\"" << color.a << "\"";
    }
    data << attr;
    data << " d=\"";
    for (const auto &p : paths) {
        data << "M " << coord::Format::to_mm(p.back().X);
        data << " " << coord::Format::to_mm(p.back().Y);
        data << " ";
        for (const auto &c : p) {
            data << "L " << coord::Format::to_mm(c.X);
            data << " " << coord::Format::to_mm(c.Y);
            data << " ";
        }
    }
    data << "\"/>\n";
}

/**
 * Adds SVG data directly to the layer, in the layer coordinate system.
 */
void Layer::add(const std::string &svg_data) {
    data << svg_data << "\n";
}

/**
 * Adds SVG data directly to the layer, in the layer coordinate system.
 */
Layer &operator<<(Layer &layer, const std::string &svg_data) {
    layer.add(svg_data);
    return layer;
}

/**
 * Writes the data to a stream.
 */
std::ostream &operator<<(std::ostream &os, const Layer &layer) {
    os << layer.data.str() << "\n";
    os << "</g>\n";
    return os;
}

/**
 * Starts rendering an SVG with the given filename, PCB bounds in
 * millimeter, and scale factor.
 */
File::File(const std::string &fname, const coord::CRect &bounds, double scale) : data(fname) {
    if (!data.is_open()) {
        throw std::runtime_error("failed to open " + fname + " for writing");
    }
    auto min_x = coord::Format::to_mm(bounds.left);
    auto min_y = coord::Format::to_mm(bounds.top);
    auto width = coord::Format::to_mm(bounds.right - bounds.left);
    auto height = coord::Format::to_mm(bounds.bottom - bounds.top);
    data << "<svg viewBox=\"" << min_x << " " << min_y << " " << width << " " << height << "\"";
    data << " width=\"" << width * scale << "\" height=\"" << height * scale << "\"";
    data << " xmlns=\"http://www.w3.org/2000/svg\">\n";
}

/**
 * Destroys this SVG writer, finishing the SVG first.
 */
File::~File() {
    close();
}

/**
 * Adds a layer to the SVG.
 */
void File::add(const Layer &layer) {
    data << layer;
}

/**
 * Adds a layer to the SVG.
 */
File &operator<<(File &file, const Layer &layer) {
    file.add(layer);
    return file;
}

/**
 * Adds SVG data directly to the SVG.
 */
void File::add(const std::string &svg_data) {
    data << svg_data;
}

/**
 * Adds SVG data directly to the SVG.
 */
File &operator<<(File &file, const std::string &svg_data) {
    file.add(svg_data);
    return file;
}

/**
 * Finishes writing the SVG.
 */
void File::close() {
    if (data.is_open()) {
        data << R"(</svg>)" << std::endl;
        data.close();
    }
}

} // namespace svg
} // namespace gerbertools
