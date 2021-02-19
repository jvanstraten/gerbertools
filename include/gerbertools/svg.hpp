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

#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include "gerbertools/coord.hpp"
#include "gerbertools/color.hpp"

namespace gerbertools {

/**
 * Contains tools for rendering SVG files.
 */
namespace svg {

/**
 * Additional attributes for SVG/XML tags.
 */
class Attributes {
private:

    /**
     * The SVG/XML data.
     */
    std::ostringstream data;

public:

    /**
     * Build function for adding a key/value pair.
     */
    Attributes with(const std::string &key, const std::string &value) &&;

    /**
     * Writes the data to a stream.
     */
    friend std::ostream &operator<<(std::ostream &os, const Attributes &attr);

};

/**
 * Represents a layer of an SVG file. This becomes a <g> node in the SVG.
 */
class Layer {
private:

    /**
     * Everything written to the SVG thus far.
     */
    std::ostringstream data;

public:

    /**
     * Constructs a layer with the given SVG identifier and an optional set of
     * attributes for the group.
     */
    Layer(const std::string &identifier, const Attributes &attr = {});

    /**
     * Adds a filled path to the layer. If the color is fully clear, the path is
     * omitted. attr may be used to add extra attributes to the path.
     */
    void add(
        coord::Paths paths,
        color::Color color = color::BLACK,
        const Attributes &attr = {}
    );

    /**
     * Adds SVG data directly to the layer, in the layer coordinate system.
     */
    void add(const std::string &svg_data);

    /**
     * Adds SVG data directly to the layer, in the layer coordinate system.
     */
    friend Layer &operator<<(Layer &layer, const std::string &svg_data);

    /**
     * Writes the data to a stream.
     */
    friend std::ostream &operator<<(std::ostream &os, const Layer &layer);

};

/**
 * A complete SVG file.
 */
class File {
private:

    /**
     * Output stream for the SVG.
     */
    std::ofstream data;

public:

    /**
     * Starts rendering an SVG with the given filename, PCB bounds in
     * millimeter, and scale factor.
     */
    File(const std::string &fname, const coord::CRect &bounds, double scale=1.0);

    /**
     * Destroys this SVG writer, finishing the SVG first.
     */
    ~File();

    /**
     * Adds a layer to the SVG.
     */
    void add(const Layer &layer);

    /**
     * Adds a layer to the SVG.
     */
    friend File &operator<<(File &file, const Layer &layer);

    /**
     * Adds SVG data directly to the SVG.
     */
    void add(const std::string &svg_data);

    /**
     * Adds SVG data directly to the SVG.
     */
    friend File &operator<<(File &file, const std::string &svg_data);

    /**
     * Finishes writing the SVG.
     */
    void close();

};

} // namespace svg
} // namespace gerbertools
