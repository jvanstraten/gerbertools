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
 * Contains tools for writing Wavefront OBJ and MTL files.
 */

#pragma once

#include <array>
#include <list>
#include <map>
#include <string>
#include "gerbertools/coord.hpp"

namespace gerbertools {

/**
 * Contains tools for writing Wavefront OBJ and MTL files.
 */
namespace obj {

/**
 * A 3-dimensional vertex.
 */
using Vertex3 = std::array<double, 3>;

/**
 * A 2-dimensional vertex.
 */
using Vertex2 = std::array<double, 2>;

/**
 * Data management class for representing things that are stored in the OBJ file
 * indirectly, referred to via 1-based list indices.
 */
template <typename T>
class Indexed {
private:

    /**
     * List of all elements.
     */
    std::list<T> elements;

    /**
     * Lookup from an element to its index.
     */
    std::map<T, size_t> indices;

public:

    /**
     * Registers an element, returning its index. If the element already exists,
     * the original index is returned, and no new element is added.
     */
    size_t add(const T &element) {
        auto it = indices.find(element);
        if (it != indices.end()) {
            return it->second;
        }
        elements.push_back(element);
        auto index = elements.size();
        indices.insert({element, index});
        return index;
    }

    typename std::list<T>::const_iterator begin() const {
        return elements.cbegin();
    }

    typename std::list<T>::const_iterator end() const {
        return elements.cend();
    }

};

// Forward declaration for an object file manager.
class ObjFile;

/**
 * Represents the corner of a face.
 */
class Corner {
private:

    /**
     * Vertex index.
     */
    size_t vertex_index;

    /**
     * UV texture coordinate index.
     */
    size_t uv_coordinate_index;

public:

    /**
     * Constructs a new corner for the given vertex. UV texture coordinates are
     * copied from the X and Y coordinates of the vertex. The vertices and UV
     * coordinates are registered with the given OBJ file manager, and their
     * indices are stored.
     */
    Corner(const Vertex3 &vertex, ObjFile &obj_file);

    /**
     * Returns the vertex index within the managing OBJ file.
     */
    size_t get_vertex_index() const;

    /**
     * Returns the texture coordinate index within the managing OBJ file.
     */
    size_t get_uv_coordinate_index() const;

};

/**
 * Represents a named object, represented within the OBJ file as a group.
 */
class Object {
private:

    /**
     * The ObjFile that constructed this object.
     */
    ObjFile &owner;

    /**
     * Object name.
     */
    std::string name;

    /**
     * Material name.
     */
    std::string material;

    /**
     * The faces that the object is made of.
     */
    std::list<std::vector<Corner>> faces;

public:

    /**
     * Creates a new object.
     */
    Object(ObjFile &owner, const std::string &name, const std::string &material);

    /**
     * Adds a face to the object with the given vertices.
     */
    void add_face(const std::vector<Vertex3> &vertices);

    /**
     * Adds a surface to the object via tesselation of an odd-even wound polygon.
     * All vertices will have the same Z coordinate.
     */
    void add_surface(const coord::Paths &polygon, double z);

    /**
     * Adds a surface to the object via tesselation of a simple polygon with
     * holes. All vertices will have the same Z coordinate.
     */
    void add_surface(const coord::Path &outline, const coord::Paths &holes, double z);

    /**
     * Adds a curved vertical ring that interconnects two surfaces.
     */
    void add_ring(const coord::Path &outline, double z1, double z2);

    /**
     * Adds a flat sheet with contour specified via an odd-even wound polygon.
     */
    void add_sheet(const coord::Paths &polygon, double z1, double z2);

    /**
     * Returns the object name.
     */
    const std::string &get_name() const;

    /**
     * Returns the material name.
     */
    const std::string &get_material() const;

    /**
     * Returns the list of faces that the object is made of.
     */
    const std::list<std::vector<Corner>> &get_faces() const;

};

/**
 * Builder for object files.
 */
class ObjFile {
private:
    friend class Object;
    friend class Corner;

    /**
     * List of all vertices.
     */
    Indexed<Vertex3> vertices;

    /**
     * List of texture coordinates.
     */
    Indexed<Vertex2> uv_coordinates;

    /**
     * List of objects added to the file. The individual objects are separated
     * by group commands.
     */
    std::list<Object> objects;

public:

    /**
     * Adds a new object to the file. A reference that can be used to add faces
     * to the object is returned.
     */
    Object &add_object(const std::string &name, const std::string &material);

    /**
     * Writes the contained OBJ file data to a file.
     */
    void to_file(const std::string &fname) const;

};

} // namespace obj
} // namespace gerbertools
