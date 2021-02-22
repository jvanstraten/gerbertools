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

#include <fstream>
#include <limits>
#include <algorithm>
#include "gerbertools/obj.hpp"
#include "gerbertools/earcut.hpp"
#include "gerbertools/clipper.hpp"

namespace gerbertools {
namespace obj {

/**
 * Constructs a new corner for the given vertex. UV texture coordinates are
 * copied from the X and Y coordinates of the vertex. The vertices and UV
 * coordinates are registered with the given OBJ file manager, and their
 * indices are stored.
 */
Corner::Corner(
    const Vertex3 &vertex,
    ObjFile &obj_file
) :
    vertex_index(obj_file.vertices.add(vertex)),
    uv_coordinate_index(obj_file.uv_coordinates.add({vertex.at(0), vertex.at(1)}))
{}

/**
 * Returns the vertex index within the managing OBJ file.
 */
size_t Corner::get_vertex_index() const {
    return vertex_index;
}

/**
 * Returns the texture coordinate index within the managing OBJ file.
 */
size_t Corner::get_uv_coordinate_index() const {
    return uv_coordinate_index;
}

/**
 * Creates a new object.
 */
Object::Object(
    ObjFile &owner,
    const std::string &name,
    const std::string &material
) :
    owner(owner),
    name(name),
    material(material)
{}

/**
 * Adds a face to the object with the given vertices.
 */
void Object::add_face(const std::vector<Vertex3> &vertices) {
    if (vertices.size() < 3) {
        throw std::invalid_argument("a face needs at least 3 corners");
    }
    auto corners = std::vector<Corner>();
    corners.reserve(vertices.size());
    for (const auto &vertex : vertices) {
        corners.emplace_back(vertex, owner);
    }
    faces.emplace_back(std::move(corners));
}

/**
 * Recursively iterates over PolyNodes to add them as surfaces to the given obj.
 */
static void poly_nodes_to_surfaces(const ClipperLib::PolyNodes &nodes, Object &obj, double z) {
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
            poly_nodes_to_surfaces(hole->Childs, obj, z);
        }
        obj.add_surface(node->Contour, holes, z);
    }
}

/**
 * Adds a surface to the object via tesselation of an odd-even wound polygon.
 * All vertices will have the same Z coordinate.
 */
void Object::add_surface(const coord::Paths &polygon, double z) {
    ClipperLib::Clipper cl;
    cl.StrictlySimple(true);
    cl.AddPaths(polygon, ClipperLib::ptSubject, true);
    ClipperLib::PolyTree tree;
    cl.Execute(ClipperLib::ctUnion, tree);
    poly_nodes_to_surfaces(tree.Childs, *this, z);
}

/**
 * Adds a surface to the object via tesselation of a simple polygon with
 * holes. All vertices will have the same Z coordinate.
 */
void Object::add_surface(const coord::Path &outline, const coord::Paths &holes, double z) {
    auto shape_data = std::vector<std::vector<Vertex2>>(1 + holes.size());
    auto vertex_data = std::vector<Vertex2>();
    size_t size = outline.size();
    for (const auto &hole : holes) {
        size += hole.size();
    }
    auto &outline_data = shape_data.at(0);
    outline_data.reserve(outline.size());
    for (const auto &coord : outline) {
        Vertex2 vert = {coord::Format::to_mm(coord.X), coord::Format::to_mm(coord.Y)};
        outline_data.push_back(vert);
        vertex_data.push_back(vert);
    }
    for (size_t i = 0; i < holes.size(); i++) {
        const auto &hole = holes.at(i);
        auto &hole_data = shape_data.at(i+1);
        hole_data.reserve(hole.size());
        for (const auto &coord : hole) {
            Vertex2 vert = {coord::Format::to_mm(coord.X), coord::Format::to_mm(coord.Y)};
            hole_data.push_back(vert);
            vertex_data.push_back(vert);
        }
    }
    auto indices = mapbox::earcut<size_t>(shape_data);
    for (size_t i = 0; i < indices.size(); i += 3) {
        auto a = vertex_data.at(indices.at(i));
        auto b = vertex_data.at(indices.at(i + 1));
        auto c = vertex_data.at(indices.at(i + 2));
        add_face({
            {a.at(0), a.at(1), z},
            {b.at(0), b.at(1), z},
            {c.at(0), c.at(1), z}
        });
    }
}

/**
 * Adds a curved vertical ring that interconnects two surfaces.
 */
void Object::add_ring(const coord::Path &outline, double z1, double z2) {
    if (outline.size() < 3) {
        throw std::invalid_argument("an outline needs at least 3 coordinates");
    }
    auto x1 = coord::Format::to_mm(outline.back().X);
    auto y1 = coord::Format::to_mm(outline.back().Y);
    for (const auto &coord : outline) {
        auto x2 = coord::Format::to_mm(coord.X);
        auto y2 = coord::Format::to_mm(coord.Y);
        add_face({
            {x1, y1, z1},
            {x2, y2, z1},
            {x2, y2, z2},
        });
        add_face({
            {x1, y1, z2},
            {x1, y1, z1},
            {x2, y2, z2},
        });
        x1 = x2;
        y1 = y2;
    }
}

/**
 * Recursively iterates over PolyNodes to add the contours as rings.
 */
static void poly_nodes_to_rings(const ClipperLib::PolyNodes &nodes, Object &obj, double z1, double z2) {
    for (const auto &node : nodes) {
        obj.add_ring(node->Contour, z1, z2);
        poly_nodes_to_rings(node->Childs, obj, z1, z2);
    }
}

/**
 * Adds a flat sheet with contour specified via an odd-even wound polygon.
 */
void Object::add_sheet(const coord::Paths &polygon, double z1, double z2) {
    ClipperLib::Clipper cl;
    cl.StrictlySimple(true);
    cl.AddPaths(polygon, ClipperLib::ptSubject, true);
    ClipperLib::PolyTree tree;
    cl.Execute(ClipperLib::ctUnion, tree);
    poly_nodes_to_surfaces(tree.Childs, *this, z1);
    poly_nodes_to_surfaces(tree.Childs, *this, z2);
    poly_nodes_to_rings(tree.Childs, *this, z1, z2);
}

/**
 * Returns the object name.
 */
const std::string &Object::get_name() const {
    return name;
}

/**
 * Returns the material name.+
 */
const std::string &Object::get_material() const {
    return material;
}

/**
 * Returns the list of faces that the object is made of.
 */
const std::list<std::vector<Corner>> &Object::get_faces() const {
    return faces;
}

/**
 * Adds a new object to the file. A reference that can be used to add faces
 * to the object is returned.
 */
Object &ObjFile::add_object(const std::string &name, const std::string &material) {
    objects.emplace_back(*this, name, material);
    return objects.back();
}

/**
 * Writes the contained OBJ file data to a file.
 */
void ObjFile::to_file(const std::string &fname) const {
    std::ofstream file(fname);
    if (!file.is_open()) {
        throw std::runtime_error("failed to open " + fname + " for writing");
    }
    for (const auto &vertex : vertices) {
        file << "v " << vertex.at(0) << " " << vertex.at(1) << " " << vertex.at(2) << "\n";
    }
    double u_min = std::numeric_limits<double>::infinity();
    double u_max = -std::numeric_limits<double>::infinity();
    double v_min = std::numeric_limits<double>::infinity();
    double v_max = -std::numeric_limits<double>::infinity();
    for (const auto &uv : uv_coordinates) {
        u_min = std::min(u_min, uv.at(0));
        u_max = std::max(u_max, uv.at(0));
        v_min = std::min(v_min, uv.at(1));
        v_max = std::max(v_max, uv.at(1));
    }
    double u_scale = 1.0 / (u_max - u_min);
    double v_scale = 1.0 / (v_max - v_min);
    for (const auto &uv : uv_coordinates) {
        file << "vt " << (uv.at(0) - u_min) * u_scale << " " << (uv.at(1) - v_min) * v_scale << "\n";
    }
    for (const auto &object : objects) {
        file << "g " << object.get_name() << "\n";
        file << "usemtl " << object.get_material() << "\n";
        for (const auto &face : object.get_faces()) {
            file << "f";
            for (const auto &corner : face) {
                file << " " << corner.get_vertex_index() << "/" << corner.get_uv_coordinate_index();
            }
            file << "\n";
        }
    }
}

} // namespace obj
} // namespace gerbertools
