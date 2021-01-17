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
#include <string>
#include <sstream>
#include <cmath>
#include <cctype>
#include <map>
#include <vector>
#include <stdexcept>
#include <memory>
#include "clipper.hpp"

/**
 * Shorthand for the Clipper namespace.
 */
namespace CL { using namespace ClipperLib; }

/**
 * Internal 64-bit fixed-point coordinate representation. Conversion from Gerber
 * units and to millimeters is handled by the CoordFormat class.
 */
using CInt = CL::cInt;

/**
 * 2D point, consisting of two CInts.
 */
using CPt = CL::IntPoint;

/**
 * Coordinate format handling. This class converts between Gerber fixed-point
 * and floating point format coordinates, an internal high-accuracy integer
 * representation used for polygon operations, and millimeters for the output.
 * It also stores maximum deviation and miter limit for the polygon operations.
 */
class CoordFormat {
private:

    /**
     * Whether the fixed point format has been configured yet.
     */
    bool fmt_configured;

    /**
     * Number of integer positions.
     */
    int n_int;

    /**
     * Number of decimal positions.
     */
    int n_dec;

    /**
     * Whether the unit (inch or millimeters) has been configured yet.
     */
    bool unit_configured;

    /**
     * Conversion factor from Gerber unit to millimeters; 25.4 for inch, 1.0 for
     * millimeters.
     */
    double factor;

    /**
     * Multiplier from Gerber fixed format to our high-accuracy 64-bit internal
     * fixed-point representation.
     */
    static const CInt PRECISION_MULT = 0x1000000ll;

    /**
     * Whether any coordinates have been converted yet. An exception is thrown
     * if the format is reconfigured after this, as the internal representation
     * changes with the Gerber format.
     */
    mutable bool used;

    /**
     * Maximum arc deviation in millimeters.
     */
    double max_deviation;

    /**
     * Miter limit in millimeters.
     */
    double miter_limit;

    /**
     * Throws an exception if the coordinate format has been used to convert
     * coordinates already.
     */
    void try_to_reconfigure() const {
        if (used) {
            throw std::runtime_error(
                "cannot reconfigure coordinate format after "
                "coordinates have already been interpreted"
            );
        }
    }

    /**
     * Throws an exception if the coordinate format has not been fully
     * configured yet.
     */
    void try_to_use() const {
        if (!fmt_configured) {
            throw std::runtime_error(
                "cannot convert coordinates before coordinate "
                "format is configured"
            );
        }
        if (!unit_configured) {
            throw std::runtime_error(
                "cannot convert coordinates before unit is configured"
            );
        }
        used = true;
    }

public:

    /**
     * Constructs a new coordinate format object with the given maximum arc
     * deviation and miter limit (both in millimeters).
     */
    explicit CoordFormat(double max_deviation=0.005, double miter_limit=1.0) :
        fmt_configured(false),
        unit_configured(false),
        used(false),
        miter_limit(miter_limit),
        max_deviation(max_deviation)
    {
    }

    /**
     * Configure the coordinate format to the given number of integer and
     * decimal digits.
     */
    void configure_format(int n_int, int n_dec) {
        try_to_reconfigure();
        fmt_configured = true;
        this->n_int = n_int;
        this->n_dec = n_dec;
    }

    /**
     * Use Freedum units.
     */
    void configure_inch() {
        try_to_reconfigure();
        unit_configured = true;
        factor = 25.4;
    }

    /**
     * Use sane units.
     */
    void configure_mm() {
        try_to_reconfigure();
        unit_configured = true;
        factor = 1.0;
    }

    /**
     * Parses a fixed-point coordinate and converts it to the internal 64-bit
     * CInt representation.
     */
    CInt parse_fixed(const std::string &s) const {
        try_to_use();
        return std::stoll(s) * PRECISION_MULT;
    }

    /**
     * Parses a floating-point coordinate and converts it to the internal 64-bit
     * CInt representation.
     */
    CInt parse_float(const std::string &s) const {
        try_to_use();
        return std::round(std::stod(s) * std::pow(10.0, n_dec) * PRECISION_MULT);
    }

    /**
     * Converts a previously parsed floating-point coordinate to the internal
     * 64-bit CInt representation.
     */
    CInt to_fixed(double d) const {
        try_to_use();
        return std::round(d * std::pow(10.0, n_dec) * PRECISION_MULT);
    }

    /**
     * Returns the maximum deviation in the internal 64-bit CInt representation.
     */
    CInt get_max_deviation() const {
        try_to_use();
        return to_fixed(max_deviation / factor);
    }

    /**
     * Returns the miter limit in the internal 64-bit CInt representation.
     */
    CInt get_miter_limit() const {
        try_to_use();
        return to_fixed(miter_limit / factor);
    }

    /**
     * Returns a ClipperOffset object with appropriate configuration.
     */
    CL::ClipperOffset build_clipper_offset() const {
        return CL::ClipperOffset(get_miter_limit(), get_max_deviation());
    }

    /**
     * Converts the internal 64-bit CInt representation to millimeters.
     */
    double to_mm(CInt i) const {
        try_to_use();
        return i / (std::pow(10.0, n_dec) * PRECISION_MULT) * factor;
    }

};

/**
 * Gerber interpolation mode.
 */
enum class InterpolationMode {

    /**
     * Mode has not been set yet.
     */
    UNDEFINED,

    /**
     * Linear interpolation (G01).
     */
    LINEAR,

    /**
     * Clockwise circular interpolation (G02).
     */
    CIRCULAR_CW,

    /**
     * Counterclockwise circular interpolation (G03).
     */
    CIRCULAR_CCW
};

/**
 * Gerber quadrant mode for circular interpolation.
 */
enum class QuadrantMode {

    /**
     * Mode has not been set yet.
     */
    UNDEFINED,

    /**
     * Insane single-quadrant mode (G74).
     */
    SINGLE,

    /**
     * Sort of reasonable multi-quadrant mode (G75).
     */
    MULTI
};

/**
 * Represents a vector image, tracking which parts were explicitly cleared as
 * well as which parts were made dark.
 */
class Plot {
private:

    /**
     * Accumulator for incoming paths. It is assumed that all paths have
     * positive winding and that the intended shape is the union of the paths
     * (nonzero/positive fill rule) unless otherwise specified. When different
     * winding modes are needed, commit_paths() must be called priot to adding
     * the new paths to the accumulator, and then commit_paths() must be called
     * again with the desired fill rule.
     */
    mutable CL::Paths accum_paths;

    /**
     * Whether the paths in the accumulator are to be interpreted as making the
     * plot dark (true) or clear (false).
     */
    mutable bool accum_polarity;

    /**
     * Set of committed paths that make the plot dark. Doesn't intersect clear.
     */
    mutable CL::Paths dark;

    /**
     * Set of committed paths that make the plot clear. Doesn't intersect dark.
     */
    mutable CL::Paths clear;

    /**
     * Whether dark and clear have been simplified yet.
     */
    mutable bool simplified = false;

    /**
     * Commits paths in the accumulator to dark/clear using the given fill type.
     * No-op if the accumulator is empty.
     */
    void commit_paths(CL::PolyFillType poly_fill_type = CL::pftNonZero) const {
        if (accum_paths.empty()) return;
        CL::SimplifyPolygons(accum_paths, poly_fill_type);
        CL::Clipper cld, clc;
        cld.AddPaths(dark, CL::ptSubject, true);
        clc.AddPaths(clear, CL::ptSubject, true);
        cld.AddPaths(accum_paths, CL::ptClip, true);
        clc.AddPaths(accum_paths, CL::ptClip, true);
        cld.Execute(
            accum_polarity ? CL::ctUnion : CL::ctDifference,
            dark, CL::pftNonZero, poly_fill_type
        );
        clc.Execute(
            accum_polarity ? CL::ctDifference : CL::ctUnion,
            clear, CL::pftNonZero, poly_fill_type
        );
        simplified = false;
        accum_paths.clear();
    }

    /**
     * Simplifies the dark/clear paths. No-op if they have already been
     * simplified.
     */
    void simplify() const {
        if (simplified) return;
        CL::SimplifyPolygons(dark, CL::pftNonZero);
        CL::SimplifyPolygons(clear, CL::pftNonZero);
        simplified = true;
    }

public:

    /**
     * Constructs a plot, optionally with initial conditions for the dark and
     * clear surfaces.
     */
    explicit Plot(
        const CL::Paths &dark = {},
        const CL::Paths &clear = {}
    ) : accum_polarity(true), dark(dark), clear(clear), simplified(false) {
    }

    /**
     * Adds paths to the plot with the given polarity. It is assumed that the
     * paths are all positively wound, and that nonzero/positive winding rules
     * apply.
     */
    void draw_paths(const CL::Paths &ps, bool polarity = true) {
        if (ps.empty()) return;

        // If the polarity is not the same as the accumulator, we have to commit
        // the accumulator first.
        if (polarity != accum_polarity) commit_paths();
        accum_polarity = polarity;

        // Simply add to the accumulator.
        accum_paths.insert(accum_paths.end(), ps.begin(), ps.end());
    }

    /**
     * Advanced method for adding paths, allowing coordinate transformations and
     * the fill rule to be specified. The transformation order is translate,
     * rotate, mirror/scale.
     */
    void draw_paths(
        const CL::Paths &ps,
        bool polarity,
        CInt translate_x,
        CInt translate_y = 0,
        bool mirror_x = false,
        bool mirror_y = false,
        double rotate = 0.0,
        double scale = 1.0,
        bool special_fill_type = false,
        CL::PolyFillType fill_type = CL::pftNonZero
    ) {
        if (ps.empty()) return;

        // If we need to apply a special fill rule, commit first, so previously
        // drawn paths won't interfere.
        if (special_fill_type) commit_paths();

        // Draw the paths normally.
        draw_paths(ps, polarity);

        // Compute transformation matrix.
        double ixx = mirror_x ? -scale : scale;
        double iyy = mirror_y ? -scale : scale;
        double si = std::sin(rotate);
        double co = std::cos(rotate);
        double xx = ixx * co;
        double xy = ixx * si;
        double yx = iyy * -si;
        double yy = iyy * co;

        // Transform the paths we just added.
        for (auto it = accum_paths.end() - ps.size(); it != accum_paths.end(); ++it) {
            for (auto &c : *it) {
                double cx = c.X * xx + c.Y * yx;
                double cy = c.X * xy + c.Y * yy;
                c.X = std::round(cx) + translate_x;
                c.Y = std::round(cy) + translate_y;
            }
        }

        // Maintain winding direction in spite of mirroring.
        if (mirror_x != mirror_y) {
            for (auto it = accum_paths.end() - ps.size(); it != accum_paths.end(); ++it) {
                CL::ReversePath(*it);
            }
        }

        // If we need to apply a special fill rule, commit immediately with said
        // fill rules.
        if (special_fill_type) commit_paths(fill_type);

    }

    /**
     * Add an entire subplot to this plot with the given transformation. The
     * transformation order is translate, rotate, mirror/scale.
     */
    void draw_plot(
        const Plot &plt,
        bool polarity = true,
        CInt translate_x = 0,
        CInt translate_y = 0,
        bool mirror_x = false,
        bool mirror_y = false,
        double rotate = 0.0,
        double scale = 1.0
    ) {
        draw_paths(plt.get_dark(), polarity, translate_x, translate_y, mirror_x, mirror_y, rotate, scale, true, CL::pftNonZero);
        draw_paths(plt.get_clear(), !polarity, translate_x, translate_y, mirror_x, mirror_y, rotate, scale, true, CL::pftNonZero);
    }

    /**
     * Returns the surface that was made dark as a simplified
     * nonzero/odd-even-filled polygon.
     */
    const CL::Paths &get_dark() const {
        commit_paths();
        simplify();
        return dark;
    }

    /**
     * Returns the surface that was explicitly made clear as a simplified
     * nonzero/odd-even-filled polygon.
     */
    const CL::Paths &get_clear() const {
        commit_paths();
        simplify();
        return clear;
    }
};

/**
 * Base class for aperture objects.
 */
class Aperture {
protected:

    /**
     * Plot object representing the shape of the aperture.
     */
    std::shared_ptr<Plot> plot;

public:

    /**
     * Returns the plot that represents the shape of the aperture.
     */
    const Plot &get_plot() const {
        return *plot;
    }

    /**
     * Returns whether this is a simple circle aperture, suitable for
     * interpolation. If this is indeed a simple circle and diameter is
     * non-null, the circle diameter is written to the pointer in addition.
     */
    virtual bool is_simple_circle(CInt *diameter) const {
        return false;
    }

};

/**
 * Represents a custom aperture object, either built using an aperture macro or
 * a block aperture.
 */
class CustomAperture : public Aperture {
public:

    /**
     * Constructs the custom aperture from the given plot.
     */
    explicit CustomAperture(const std::shared_ptr<Plot> &data) {
        plot = data;
    }

};

/**
 * Base class for standard apertures.
 */
class StandardAperture : public Aperture {
protected:

    /**
     * Diameter of the optional hole common to all standard apertures.
     */
    CInt hole_diameter;

    /**
     * Returns the path for the hole. That is, a negatively wound circle, or
     * no path at all if the diameter is zero.
     */
    CL::Paths get_hole(const CoordFormat &fmt) const {
        CL::Paths ps;
        if (hole_diameter > 0.0) {
            auto co = fmt.build_clipper_offset();
            co.AddPath({{0, 0}}, CL::jtRound, CL::etOpenRound);
            co.Execute(ps, hole_diameter * 0.5);
            CL::ReversePaths(ps);
        }
        return ps;
    }
};

/**
 * Represents a standard circle aperture.
 */
class CircleAperture : public StandardAperture {
private:

    /**
     * Diameter of the circle.
     */
    CInt diameter;

public:

    /**
     * Constructs the circle aperture from the given parameters, reported as a
     * vector of strings. The first string is ignored; it is assumed to be "C"
     * to select this type of aperture.
     */
    explicit CircleAperture(const std::vector<std::string> &csep, const CoordFormat &fmt) {

        // Parse the command.
        if (csep.size() < 2 || csep.size() > 3) {
            throw std::runtime_error("invalid circle aperture");
        }
        diameter = fmt.parse_float(csep.at(1));
        hole_diameter = (csep.size() > 2) ? fmt.parse_float(csep.at(2)) : 0;

        // Construct the plot.
        CL::Paths ps;
        auto co = fmt.build_clipper_offset();
        co.AddPath({{0, 0}}, CL::jtRound, CL::etOpenRound);
        co.Execute(ps, diameter * 0.5);
        auto hole = get_hole(fmt);
        ps.insert(ps.end(), hole.begin(), hole.end());
        plot = std::make_shared<Plot>(ps);

    }

    /**
     * Returns whether this is a simple circle aperture, suitable for
     * interpolation. If this is indeed a simple circle and diameter is
     * non-null, the circle diameter is written to the pointer in addition.
     */
    bool is_simple_circle(CInt *diameter) const override {
        if (hole_diameter > 0.0) return false;
        if (diameter) *diameter = this->diameter;
        return true;
    }

};

/**
 * Represents a standard rectangle aperture.
 */
class RectangleAperture : public StandardAperture {
private:

    /**
     * Horizontal size of the aperture.
     */
    CInt x_size;

    /**
     * Vertical size of the aperture.
     */
    CInt y_size;

public:

    /**
     * Constructs the rectangle aperture from the given parameters, reported as
     * a vector of strings. The first string is ignored; it is assumed to be "R"
     * to select this type of aperture.
     */
    explicit RectangleAperture(const std::vector<std::string> &csep, const CoordFormat &fmt) {

        // Parse the command.
        if (csep.size() < 3 || csep.size() > 4) {
            throw std::runtime_error("invalid rectangle aperture");
        }
        x_size = std::abs(fmt.parse_float(csep.at(1)));
        y_size = std::abs(fmt.parse_float(csep.at(2)));
        hole_diameter = (csep.size() > 3) ? fmt.parse_float(csep.at(3)) : 0;

        // Construct the plot.
        CInt x = x_size / 2;
        CInt y = y_size / 2;
        CL::Paths ps{{
             {x, y},
             {x, -y},
             {-x, -y},
             {-x, y}
        }};
        auto hole = get_hole(fmt);
        ps.insert(ps.end(), hole.begin(), hole.end());
        plot = std::make_shared<Plot>(ps);

    }
};

/**
 * Represents a standard obround aperture.
 */
class ObroundAperture : public StandardAperture {
private:

    /**
     * Horizontal size of the aperture.
     */
    CInt x_size;

    /**
     * Vertical size of the aperture.
     */
    CInt y_size;

public:

    /**
     * Constructs the obround aperture from the given parameters, reported as
     * a vector of strings. The first string is ignored; it is assumed to be "O"
     * to select this type of aperture.
     */
    explicit ObroundAperture(const std::vector<std::string> &csep, const CoordFormat &fmt) {

        // Parse the command.
        if (csep.size() < 3 || csep.size() > 4) {
            throw std::runtime_error("invalid obround aperture");
        }
        x_size = std::abs(fmt.parse_float(csep.at(1)));
        y_size = std::abs(fmt.parse_float(csep.at(2)));
        hole_diameter = (csep.size() > 3) ? fmt.parse_float(csep.at(3)) : 0;

        // Construct the plot.
        CInt x = x_size / 2;
        CInt y = y_size / 2;
        CInt r = std::min(x, y);
        x -= r;
        y -= r;
        CL::Paths ps;
        auto co = fmt.build_clipper_offset();
        co.AddPath({{-x, -y}, {x, y}}, CL::jtRound, CL::etOpenRound);
        co.Execute(ps, r);
        auto hole = get_hole(fmt);
        ps.insert(ps.end(), hole.begin(), hole.end());
        plot = std::make_shared<Plot>(ps);

    }
};

/**
 * Represents a standard regular polygon aperture.
 */
class PolygonAperture : public StandardAperture {
private:

    /**
     * Outer diameter of the regular polygon.
     */
    CInt diameter;

    /**
     * Number of vertices. Should be at least 3.
     */
    size_t n_vertices;

    /**
     * Rotation of the polygon in radians counterclockwise from the positive X
     * axis.
     */
    double rotation;

public:

    /**
     * Constructs the polygon aperture from the given parameters, reported as
     * a vector of strings. The first string is ignored; it is assumed to be "P"
     * to select this type of aperture.
     */
    explicit PolygonAperture(const std::vector<std::string> &csep, const CoordFormat &fmt) {

        // Parse the command.
        if (csep.size() < 3 || csep.size() > 5) {
            throw std::runtime_error("invalid polygon aperture");
        }
        diameter = fmt.parse_float(csep.at(1));
        n_vertices = std::stoul(csep.at(2));
        if (n_vertices < 3) {
            throw std::runtime_error("invalid polygon aperture");
        }
        rotation = (csep.size() > 3) ? (std::stod(csep.at(3)) / 180.0 * M_PI) : 0.0;
        hole_diameter = (csep.size() > 4) ? fmt.parse_float(csep.at(4)) : 0;

        // Construct the plot.
        CL::Paths ps = {{}};
        for (size_t i = 0; i < n_vertices; i++) {
            double a = ((double)i / (double)n_vertices) * 2.0 * M_PI;
            ps.back().push_back({
                (CInt)std::round(diameter * 0.5 * std::cos(a + rotation)),
                (CInt)std::round(diameter * 0.5 * std::sin(a + rotation))
            });
        }
        auto hole = get_hole(fmt);
        ps.insert(ps.end(), hole.begin(), hole.end());
        plot = std::make_shared<Plot>(ps);

    }
};

namespace aperture_macro {

class Expression;
using ExpressionRef = std::shared_ptr<Expression>;
using ExpressionRefs = std::list<ExpressionRef>;
using Variables = std::map<size_t, double>;

/**
 * Aperture macro expression tree. Can also take the shape of a single token
 * while parsing.
 */
class Expression {
private:
    static std::string debug(
        const ExpressionRefs &expr,
        ExpressionRefs::iterator expr_begin,
        ExpressionRefs::iterator expr_end
    );
    static ExpressionRef reduce(
        ExpressionRefs &expr,
        ExpressionRefs::iterator expr_begin,
        ExpressionRefs::iterator expr_end
    );
public:
    static ExpressionRef parse(std::string expr);
    
    /**
     * Evaluate this expression with the given set of variables.
     */
    virtual double eval(const Variables &vars) const = 0;
    
    /**
     * If this is a character token, return the character it represents.
     * Otherwise retun '\0'.
     */
    virtual char get_token() const { return 0; }
    
    /**
     * Returns a debug representation of this expression node. 
     */
    virtual std::string debug() const = 0;
};

/**
 * Represents a literal.
 */
class LiteralExpression : public Expression {
private:

    /**
     * Value of the literal.
     */
    double value;

public:

    /**
     * Constructs a literal node with the given value.
     */
    explicit LiteralExpression(double value) : value(value) {
    }

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override {
        return value;
    }

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override {
        return "<" + std::to_string(value) + ">";
    }

};

/**
 * Represents a variable reference.
 */
class VariableExpression : public Expression {
private:

    /**
     * Index of the variable.
     */
    size_t index;

public:

    /**
     * Constructs a variable reference node with the given variable index.
     */
    explicit VariableExpression(size_t index) : index(index) {
    }

    /**
     * Constructs a variable node for the given variable index.
     */
    size_t get_index() const {
        return index;
    }

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override {
        auto it = vars.find(index);
        if (it != vars.end()) {
            return it->second;
        } else {
            return 0.0;
        }
    }

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override {
        return "<$" + std::to_string(index) + ">";
    }

};

/**
 * Represents a unary operation, either -x or +x.
 */
class UnaryExpression : public Expression {
private:

    /**
     * The operation, either '-' or '+'.
     */
    char oper;

    /**
     * The operand expression.
     */
    ExpressionRef expr;

public:

    /**
     * Constructs a unary expression node from the given operand character
     * (either '+' or '-') and operand expression.
     */
    UnaryExpression(
        char oper,
        const ExpressionRef &expr
    ) : oper(oper), expr(expr) {
    }

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override {
        if (oper == '+') {
            return expr->eval(vars);
        } else if (oper == '-') {
            return -expr->eval(vars);
        } else {
            throw std::runtime_error("invalid unary operator in aperture macro");
        }
    }

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override {
        return "<" + std::string(1, oper) + expr->debug() + ">";
    }

};

/**
 * Represents a unary operation, either x+y, x-y, x*y, or x/y.
 */
class BinaryExpression : public Expression {
private:

    /**
     * The operation, either '-' or '+'.
     */
    char oper;

    /**
     * Left-hand-side operand expression.
     */
    ExpressionRef lhs;

    /**
     * Right-hand-side operand expression.
     */
    ExpressionRef rhs;

public:

    /**
     * Constructs a binary expression node from the given operand character
     * (either '+', '-', 'x', or '/') and operand expressions.
     */
    BinaryExpression(
        char oper,
        const ExpressionRef &lhs,
        const ExpressionRef &rhs
    ) : oper(oper), lhs(lhs), rhs(rhs) {
    }

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override {
        if (oper == '+') {
            return lhs->eval(vars) + rhs->eval(vars);
        } else if (oper == '-') {
            return lhs->eval(vars) - rhs->eval(vars);
        } else if (oper == 'x') {
            return lhs->eval(vars) * rhs->eval(vars);
        } else if (oper == '/') {
            return lhs->eval(vars) / rhs->eval(vars);
        } else {
            throw std::runtime_error("invalid binary operator in aperture macro");
        }
    }

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override {
        return "<" + lhs->debug() + std::string(1, oper) + rhs->debug() + ">";
    }

};

/**
 * Represents an unparsed token.
 */
class Token : public Expression {
private:

    /**
     * The only tokens we need to represent are single characters, so a char is
     * sufficient for identification. These tokens are '(', ')', '+', '-', 'x',
     * and '/'.
     */
    char token;

public:

    /**
     * Constructs a token character. The character must be '(', ')', '+', '-',
     * 'x', or '/'.
     */
    Token(char token) : token(token) {
    }

    /**
     * Evaluate this expression with the given set of variables.
     */
    double eval(const Variables &vars) const override {
        throw std::runtime_error("cannot evaluate token");
    }
    char get_token() const override {
        return token;
    }

    /**
     * Returns a debug representation of this expression node.
     */
    std::string debug() const override {
        return "<" + std::string(1, token) + ">";
    }

};

/**
 * Constructs a string representation of the given expression/token list and
 * context.
 */
std::string Expression::debug(
    const ExpressionRefs &expr,
    ExpressionRefs::iterator expr_begin,
    ExpressionRefs::iterator expr_end
) {
    std::ostringstream ss;
    ss << "[";
    for (auto it = expr.begin(); it != expr.end(); ++it) {
        if (it == expr_end) {
            ss << "}";
        }
        if (it == expr_begin) {
            ss << "{";
        }
        ss << " " << (*it)->debug() << " ";
    }
    if (expr_end == expr.end()) {
        ss << "}";
    }
    if (expr_begin == expr.end()) {
        ss << "{";
    }
    ss << "]";
    return ss.str();
}

/**
 * Reduces the given expression/token list to a single expression tree by
 * applying reduction rules. Yes, it might've been better to just include a
 * proper parser lib and use that. This was easier (by some metrics).
 */
ExpressionRef Expression::reduce(
    ExpressionRefs &expr,
    ExpressionRefs::iterator expr_begin,
    ExpressionRefs::iterator expr_end
) {

    // At least one expression/token is needed.
    if (expr_begin == expr_end) {
        throw std::runtime_error("empty aperture macro (sub)expression");
    }

    // Replace parentheses.
    for (auto begin = expr_begin; begin != expr_end; ++begin) {
        if ((*begin)->get_token() == '(') {
            size_t level = 1;
            for (auto end = std::next(begin); end != expr_end; ++end) {
                char t = (*end)->get_token();
                if (t == '(') {
                    level++;
                } else if (t == ')') {
                    level--;
                    if (!level) {
                        auto e = reduce(expr, std::next(begin), end);
                        *begin = e;
                        expr.erase(
                            std::next(begin),
                            std::next(end)
                        );
                        break;
                    }
                }
            }
            if (level) {
                throw std::runtime_error("unmatched ( in aperture macro expression");
            }
        }
    }
    for (auto it = expr_begin; it != expr_end; ++it) {
        char t = (*it)->get_token();
        if (t == ')') {
            throw std::runtime_error("unmatched ) in aperture macro expression");
        } else if (t == '(') {
            throw std::runtime_error(
                "internal error in aperture macro expression parser; state = "
                + debug(expr, expr_begin, expr_end));
        }
    }

    // Handle unary operators.
    auto begin = expr_begin;
    bool ok = false;
    for (auto end = expr_begin; end != expr_end; ++end) {
        if ((*end)->get_token() == 0) {

            // Range from begin (incl) to end (excl) contains the unary
            // operators. end itself is the expression to apply the unaries to.
            // Apply in right to left order.
            while (begin != end) {
                auto un = std::prev(end);
                char t = (*un)->get_token();
                if (t != '-' && t != '+') {
                    throw std::runtime_error("invalid unary operator in aperture macro expression");
                }
                *un = std::make_shared<UnaryExpression>(t, *end);
                expr.erase(end--);
            }

            // end+1 should be a binary operator and end+2 should be the next
            // start of a unary. It's also possible that end+1 is the end of the
            // expression, if this was the last one. end+2 being the last is an
            // invalid expression.
            if (++end == expr_end) {
                ok = true;
                break;
            }
            if (std::next(end) == expr_end) {
                throw std::runtime_error("unterminated binary operator in aperture macro expression");
            }
            begin = std::next(end);
        }
    }
    if (!ok) {
        throw std::runtime_error("unterminated unary operator in aperture macro expression");
    }

    // Handle binary operators.
    for (int prec = 0; prec < 2; prec++) {
        for (auto it = expr_begin; it != expr_end; ++it) {
            char t = (*it)->get_token();
            if (
                (prec == 0 && (t == 'x' || t == '/')) ||
                (prec == 1 && (t == '+' || t == '-'))
            ) {
                if (it == expr_begin || std::next(it) == expr_end) {
                    throw std::runtime_error(
                        "internal error in aperture macro expression parser; state = "
                        + debug(expr, expr_begin, expr_end));
                }
                auto lhs = std::prev(it);
                auto rhs = std::next(it);
                *lhs = std::make_shared<BinaryExpression>(t, *lhs, *rhs);
                expr.erase(it);
                expr.erase(rhs);
                it = lhs;
            }
        }
    }

    // Expression list should now be a single entry that isn't a token.
    if (expr_begin == expr_end || next(expr_begin) != expr_end) {
        throw std::runtime_error(
            "internal error in aperture macro expression parser; state = "
            + debug(expr, expr_begin, expr_end));
    }
    return *expr_begin;

}

/**
 * Parse the given expression.
 */
ExpressionRef Expression::parse(std::string expr) {

    // Tokenize input.
    expr += " ";
    ExpressionRefs tokens;
    std::string cur;
    enum Mode {IDLE, LIT, VAR};
    Mode mode = IDLE;
    for (auto c : expr) {
        if (mode == LIT) {
            if ((c >= '0' && c <= '9') || (c == '.')) {
                cur += c;
                continue;
            } else {
                tokens.push_back(std::make_shared<LiteralExpression>(std::stod(cur)));
                cur.resize(0);
                mode = IDLE;
            }
        } else if (mode == VAR) {
            if (c >= '0' && c <= '9') {
                cur += c;
                continue;
            } else {
                tokens.push_back(std::make_shared<VariableExpression>(std::stoul(cur)));
                cur.resize(0);
                mode = IDLE;
            }
        }
        if ((c >= '0' && c <= '9') || (c == '.')) {
            cur += c;
            mode = LIT;
        } else if (c == '$') {
            mode = VAR;
        } else if (c == '-' || c == '+' || c == 'x' || c == '/' || c == '(' || c == ')') {
            tokens.push_back(std::make_shared<Token>(c));
        }
    }

    // Reduce expression.
    return reduce(tokens, tokens.begin(), tokens.end());
}

/**
 * Represents an aperture macro command... poorly.
 */
using ApertureMacroCommand = std::vector<ExpressionRef>;

/**
 * Represents an aperture macro (before being instantiated into an aperture).
 */
class ApertureMacro {
private:

    /**
     * The commands that make up the macro.
     */
    std::list<ApertureMacroCommand> cmds;

public:

    /**
     * Parses and appends an aperture macro command.
     */
    void append(const std::string &cmd) {
        if (cmd.empty()) {
            throw std::runtime_error("empty aperture macro command");
        } else if (cmd.at(0) == '0') {
            // Comment, ignore.
        } else if (cmd.at(0) == '$') {
            // Assignment command.
            auto pos = cmd.find('=');
            if (pos == std::string::npos) {
                throw std::runtime_error("invalid aperture macro assignment command");
            }
            cmds.push_back({
               Expression::parse(cmd.substr(0, pos)),
               Expression::parse(cmd.substr(pos + 1))
            });
        } else {
            size_t pos = 0;
            cmds.push_back({});
            while (true) {
                auto pos2 = cmd.find(',', pos);
                cmds.back().push_back(Expression::parse(cmd.substr(pos, pos2 - pos)));
                if (pos2 == std::string::npos) {
                    break;
                }
                pos = pos2 + 1;
            }
        }
    }

    /**
     * Executes the macro to construct an aperture using the given parameters,
     * reported as a vector of strings. The first string is ignored; it is
     * assumed to be the name of the aperture macro.
     */
    std::shared_ptr<Aperture> build(const std::vector<std::string> &csep, const CoordFormat &fmt) {
        Variables vars;
        for (size_t i = 1; i < csep.size(); i++) {
            vars[i] = std::stod(csep.at(i));
        }
        Plot plot;
        for (const auto &cmd : cmds) {

            // Handle variable assignment commands.
            if (cmd.size() == 2) {
                auto v = std::dynamic_pointer_cast<VariableExpression>(cmd.at(0));
                if (v) {
                    vars[v->get_index()] = cmd.at(1)->eval(vars);
                    continue;
                }
            }

            // Handle draw commands.
            // TODO: this is terrible design. Commands should be classes that
            //  check the command syntax upon construction, and have a
            //  polymorphic execute function. The actual rendering code is also
            //  messy, but whatever, it works.
            auto code = (size_t)std::round(cmd.at(0)->eval(vars));
            switch (code) {
                case 1: {

                    // Circle.
                    if (cmd.size() < 5 || cmd.size() > 6) {
                        throw std::runtime_error("invalid circle command in aperture macro");
                    }
                    bool exposure = cmd.at(1)->eval(vars) > 0.5;
                    double diameter = std::abs(cmd.at(2)->eval(vars));
                    double center_x = cmd.at(3)->eval(vars);
                    double center_y = cmd.at(4)->eval(vars);
                    double rotation = (cmd.size() > 5) ? cmd.at(5)->eval(vars) : 0.0;
                    CL::Paths ps;
                    auto co = fmt.build_clipper_offset();
                    co.AddPath({
                        {
                            fmt.to_fixed(center_x),
                            fmt.to_fixed(center_y)
                        }
                    }, CL::jtRound, CL::etOpenRound);
                    co.Execute(ps, fmt.to_fixed(diameter * 0.5));
                    plot.draw_paths(
                        ps, exposure,
                        0, 0,
                        false, false,
                        rotation / 180.0 * M_PI
                    );
                    break;

                }
                case 20: {

                    // Vector line.
                    if (cmd.size() < 7 || cmd.size() > 8) {
                        throw std::runtime_error("invalid vector line command in aperture macro");
                    }
                    bool exposure = cmd.at(1)->eval(vars) > 0.5;
                    double width = std::abs(cmd.at(2)->eval(vars));
                    double start_x = cmd.at(3)->eval(vars);
                    double start_y = cmd.at(4)->eval(vars);
                    double end_x = cmd.at(5)->eval(vars);
                    double end_y = cmd.at(6)->eval(vars);
                    double rotation = (cmd.size() > 7) ? cmd.at(7)->eval(vars) : 0.0;
                    CL::Paths ps;
                    auto co = fmt.build_clipper_offset();
                    co.AddPath({
                        {
                            fmt.to_fixed(start_x),
                            fmt.to_fixed(start_y)
                        }, {
                            fmt.to_fixed(end_x),
                            fmt.to_fixed(end_y)
                        }
                    }, CL::jtMiter, CL::etOpenButt);
                    co.Execute(ps, fmt.to_fixed(width * 0.5));
                    plot.draw_paths(
                        ps, exposure,
                        0, 0,
                        false, false,
                        rotation / 180.0 * M_PI
                    );
                    break;

                }
                case 21: {

                    // Center line (aka rectangle...).
                    if (cmd.size() < 6 || cmd.size() > 7) {
                        throw std::runtime_error("invalid center line command in aperture macro");
                    }
                    bool exposure = cmd.at(1)->eval(vars) > 0.5;
                    double width = std::abs(cmd.at(2)->eval(vars));
                    double height = std::abs(cmd.at(3)->eval(vars));
                    double center_x = cmd.at(4)->eval(vars);
                    double center_y = cmd.at(5)->eval(vars);
                    double rotation = (cmd.size() > 6) ? cmd.at(6)->eval(vars) : 0.0;
                    CL::Paths ps = {{
                        {
                            fmt.to_fixed(center_x + width * 0.5),
                            fmt.to_fixed(center_y + height * 0.5)
                        },
                        {
                            fmt.to_fixed(center_x - width * 0.5),
                            fmt.to_fixed(center_y + height * 0.5)
                        },
                        {
                            fmt.to_fixed(center_x - width * 0.5),
                            fmt.to_fixed(center_y - height * 0.5)
                        },
                        {
                            fmt.to_fixed(center_x + width * 0.5),
                            fmt.to_fixed(center_y - height * 0.5)
                        }
                    }};
                    plot.draw_paths(
                        ps, exposure,
                        0, 0,
                        false, false,
                        rotation / 180.0 * M_PI
                    );
                    break;

                }
                case 4: {

                    // Outline.
                    if (cmd.size() < 3) {
                        throw std::runtime_error("invalid outline command in aperture macro");
                    }
                    bool exposure = cmd.at(1)->eval(vars) > 0.5;
                    auto n_vertices = (size_t)std::round(cmd.at(2)->eval(vars));
                    size_t rotation_index = 5 + 2*n_vertices;
                    if (n_vertices < 3 || cmd.size() < rotation_index || cmd.size() > rotation_index + 1) {
                        throw std::runtime_error("invalid outline command in aperture macro");
                    }
                    double rotation = (cmd.size() > rotation_index) ? cmd.at(rotation_index)->eval(vars) : 0.0;
                    CL::Paths ps = {{}};
                    for (size_t i = 0; i < n_vertices; i++) {
                        ps.back().push_back({
                            fmt.to_fixed(cmd.at(3 + 2*i)->eval(vars)),
                            fmt.to_fixed(cmd.at(4 + 2*i)->eval(vars))
                        });
                    };
                    plot.draw_paths(
                        ps, exposure,
                        0, 0,
                        false, false,
                        rotation / 180.0 * M_PI, 1.0,
                        true, ClipperLib::pftNonZero
                    );
                    break;

                }
                case 5: {

                    // Polygon.
                    if (cmd.size() < 6 || cmd.size() > 7) {
                        throw std::runtime_error("invalid polygon command in aperture macro");
                    }
                    bool exposure = cmd.at(1)->eval(vars) > 0.5;
                    auto n_vertices = (size_t)std::round(std::abs(cmd.at(2)->eval(vars)));
                    double center_x = cmd.at(3)->eval(vars);
                    double center_y = cmd.at(4)->eval(vars);
                    double diameter = std::abs(cmd.at(5)->eval(vars));
                    double rotation = (cmd.size() > 6) ? cmd.at(6)->eval(vars) : 0.0;
                    CL::Paths ps = {{}};
                    for (size_t i = 0; i < n_vertices; i++) {
                        double a = ((double)i / (double)n_vertices) * 2.0 * M_PI;
                        ps.push_back({
                            fmt.to_fixed(center_x + diameter * 0.5 * std::cos(a)),
                            fmt.to_fixed(center_y + diameter * 0.5 * std::cos(a))
                        });
                    };
                    plot.draw_paths(
                        ps, exposure,
                        0, 0,
                        false, false,
                        rotation / 180.0 * M_PI
                    );
                    break;

                }
                case 6: {

                    // Moire.
                    if (cmd.size() < 9 || cmd.size() > 10) {
                        throw std::runtime_error("invalid moire command in aperture macro");
                    }
                    double center_x = cmd.at(1)->eval(vars);
                    double center_y = cmd.at(2)->eval(vars);
                    double diameter = std::abs(cmd.at(3)->eval(vars));
                    double thickness = std::abs(cmd.at(4)->eval(vars));
                    double gap = std::abs(cmd.at(5)->eval(vars));
                    auto max_rings = (size_t)std::round(std::abs(cmd.at(6)->eval(vars)));
                    double ch_thickness = std::abs(cmd.at(7)->eval(vars));
                    double ch_length = std::abs(cmd.at(8)->eval(vars));
                    double rotation = (cmd.size() > 9) ? cmd.at(9)->eval(vars) : 0.0;
                    CL::Paths ps = {};
                    for (size_t i = 0; i < max_rings*2 && diameter > 0.0; i++) {
                        auto co = fmt.build_clipper_offset();
                        co.AddPath({
                            {
                                fmt.to_fixed(center_x),
                                fmt.to_fixed(center_y)
                            }
                        }, CL::jtRound, CL::etOpenRound);
                        CL::Paths cps;
                        co.Execute(cps, fmt.to_fixed(diameter * 0.5));
                        if (i & 1) {
                            CL::ReversePaths(cps);
                            diameter -= gap * 2.0;
                        } else {
                            diameter -= thickness * 2.0;
                        }
                        ps.insert(ps.end(), cps.begin(), cps.end());
                    };
                    if (ch_thickness > 0.0 && ch_length > 0.0) {
                        ps.push_back({
                            {
                                fmt.to_fixed(center_x + ch_thickness * 0.5),
                                fmt.to_fixed(center_y + ch_length * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x - ch_thickness * 0.5),
                                fmt.to_fixed(center_y + ch_length * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x - ch_thickness * 0.5),
                                fmt.to_fixed(center_y - ch_length * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x + ch_thickness * 0.5),
                                fmt.to_fixed(center_y - ch_length * 0.5)
                            }
                        });
                        ps.push_back({
                            {
                                fmt.to_fixed(center_x + ch_length * 0.5),
                                fmt.to_fixed(center_y + ch_thickness * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x - ch_length * 0.5),
                                fmt.to_fixed(center_y + ch_thickness * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x - ch_length * 0.5),
                                fmt.to_fixed(center_y - ch_thickness * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x + ch_length * 0.5),
                                fmt.to_fixed(center_y - ch_thickness * 0.5)
                            }
                        });
                    }
                    plot.draw_paths(
                        ps, true,
                        0, 0,
                        false, false,
                        rotation / 180.0 * M_PI, 1.0,
                        true, ClipperLib::pftPositive
                    );
                    break;

                }
                case 7: {

                    // Thermal.
                    if (cmd.size() < 6 || cmd.size() > 7) {
                        throw std::runtime_error("invalid moire command in aperture macro");
                    }
                    double center_x = cmd.at(1)->eval(vars);
                    double center_y = cmd.at(2)->eval(vars);
                    double outer = std::abs(cmd.at(3)->eval(vars));
                    double inner = std::abs(cmd.at(4)->eval(vars));
                    double gap = std::abs(cmd.at(5)->eval(vars));
                    double rotation = (cmd.size() > 6) ? cmd.at(6)->eval(vars) : 0.0;
                    CL::Paths ps = {};

                    auto co1 = fmt.build_clipper_offset();
                    co1.AddPath({
                        {
                            fmt.to_fixed(center_x),
                            fmt.to_fixed(center_y)
                        }
                    }, CL::jtRound, CL::etOpenRound);
                    CL::Paths cps;
                    co1.Execute(cps, fmt.to_fixed(outer * 0.5));
                    ps.insert(ps.end(), cps.begin(), cps.end());

                    auto co2 = fmt.build_clipper_offset();
                    co2.AddPath({
                        {
                            fmt.to_fixed(center_x),
                            fmt.to_fixed(center_y)
                        }
                    }, CL::jtRound, CL::etOpenRound);
                    cps.clear();
                    co2.Execute(cps, fmt.to_fixed(inner * 0.5));
                    CL::ReversePaths(cps);
                    ps.insert(ps.end(), cps.begin(), cps.end());

                    if (gap > 0.0) {
                        ps.push_back({
                            {
                                fmt.to_fixed(center_x + gap * 0.5),
                                fmt.to_fixed(center_y + outer * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x + gap * 0.5),
                                fmt.to_fixed(center_y - outer * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x - gap * 0.5),
                                fmt.to_fixed(center_y - outer * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x - gap * 0.5),
                                fmt.to_fixed(center_y + outer * 0.5)
                            }
                        });
                        ps.push_back({
                            {
                                fmt.to_fixed(center_x + outer * 0.5),
                                fmt.to_fixed(center_y + gap * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x + outer * 0.5),
                                fmt.to_fixed(center_y - gap * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x - outer * 0.5),
                                fmt.to_fixed(center_y - gap * 0.5)
                            },
                            {
                                fmt.to_fixed(center_x - outer * 0.5),
                                fmt.to_fixed(center_y + gap * 0.5)
                            }
                        });
                    }
                    plot.draw_paths(
                        ps, true,
                        0, 0,
                        false, false,
                        rotation / 180.0 * M_PI, 1.0,
                        true, ClipperLib::pftPositive
                    );
                    break;

                }
                default: {
                    throw std::runtime_error("invalid aperture macro primitive code");
                }

            }

        }

        auto ap_plot = std::make_shared<Plot>();
        ap_plot->draw_paths(plot.get_dark());
        return std::make_shared<CustomAperture>(ap_plot);
    }
};

} // namespace aperture_macro

/**
 * Helper class for circular interpolations.
 */
class CircularInterpolationHelper {
private:

    /**
     * Center point for the arc.
     */
    double center_x, center_y;

    /**
     * Initial radius and angle.
     */
    double r1, r2;

    /**
     * Final radius and angle.
     */
    double a1, a2;

    /**
     * Converts the given coordinate to polar.
     */
    void to_polar(double x, double y, double &r, double &a) const {
        x -= center_x;
        y -= center_y;
        r = std::hypot(x, y);
        a = std::atan2(y, x);
    }

public:

    /**
     * Constructs a circular interpolation from one coordinate to another, given
     * a center coordinate and interpolation mode. Because 6 parameters is
     * over-constrained for an arc, the radius is interpolated as well.
     */
    CircularInterpolationHelper(CPt start, CPt end, CPt center, bool ccw, bool multi) {
        center_x = center.X;
        center_y = center.Y;
        to_polar(start.X, start.Y, r1, a1);
        to_polar(end.X, end.Y, r2, a2);
        if (multi) {
            if (ccw) {
                if (a2 <= a1) a2 += 2.0 * M_PI;
            } else {
                if (a1 <= a2) a1 += 2.0 * M_PI;
            }
        } else {
            if (ccw) {
                if (a2 < a1) a2 += 2.0 * M_PI;
            } else {
                if (a1 < a2) a1 += 2.0 * M_PI;
            }
        }
    }

    /**
     * Returns whether this arc qualifies for single-quadrant mode.
     */
    bool is_single_quadrant() const {
        return std::abs(a1 - a2) <= M_PI_2 + 1e-3;
    }

    /**
     * Returns a metric for comparing round-off error between different possible
     * interpolations. Because we're interpolating radius to solve the
     * over-constrained problem, the error is just the difference between the
     * two radii.
     */
    double error() const {
        return std::max(r1, r2);
    }

    /**
     * Interpolates the arc to a path with a given maximum deviation.
     */
    CL::Path to_path(double epsilon) const {
        double r = (r1 + r2) * 0.5;
        double x = (r > epsilon) ? (1.0 - epsilon / r) : 0.0;
        double th = std::acos(2.0 * x * x - 1.0) + 1e-3;
        auto n_vertices = (size_t)std::ceil(std::abs(a2 - a1) / th);
        CL::Path p;
        for (size_t i = 0; i <= n_vertices; i++) {
            double f2 = (double)i / (double)n_vertices;
            double f1 = 1.0 - f2;
            double vr = f1*r1 + f2*r2;
            double va = f1*a1 + f2*a2;
            double vx = center_x + vr * std::cos(va);
            double vy = center_y + vr * std::sin(va);
            p.push_back({(CInt)std::round(vx), (CInt)std::round(vy)});
        }
        return p;
    }

};

/**
 * Main class for parsing Gerber files. The input is parsed during construction
 * of the object.
 */
class GerberFile {
private:

    /**
     * Mapping from aperture index (10..infinity) to Aperture objects. Each
     * aperture is either one of the StandardAperture specializations or a
     * CustomAperture, representing either an aperture macro or a block aperture
     * depending on how it was built.
     */
    std::map<size_t, std::shared_ptr<Aperture>> apertures;

    /**
     * Mapping from name to aperture macro data. When an aperture macro is
     * instantiated, a plot is built on-the-fly to replay the macro instructions
     * with the given parameters.
     */
    std::map<std::string, std::shared_ptr<aperture_macro::ApertureMacro>> aperture_macros;

    /**
     * When non-null, an aperture macro is being constructed. Commands are
     * expected to be in the comma-separated aperture macro command syntax, to
     * be variable assignment statements of the form $<var>=<expr>, or to be
     * aperture macro comments. The attribute termination character (%) always
     * resets this to nullptr.
     */
    std::shared_ptr<aperture_macro::ApertureMacro> am_builder;

    /**
     * Stack of Plot images. The first entry is the actual image; this always
     * exists. An additional plot is pushed onto the stack when a block aperture
     * is started. Graphics objects are always pushed onto the last entry.
     */
    std::vector<std::shared_ptr<Plot>> plot_stack;

    /**
     * Coordinate format information, including expected digit counts and
     * mm/inch switch. Note that only leading zeros may be stripped; stripping
     * trailing zeros is a deprecated Gerber feature and is not supported.
     */
    CoordFormat fmt;

    /**
     * Stores linear vs. circular interpolation mode.
     */
    InterpolationMode imode;

    /**
     * Stores circular interpolation quadrant mode.
     */
    QuadrantMode qmode;

    /**
     * The currently selected aperture, or nullptr if none have been selected
     * yet.
     */
    std::shared_ptr<Aperture> aperture;

    /**
     * Current coordinate.
     */
    CPt pos;

    /**
     * Current polarity. True for dark, false for clear. Essentially, when
     * polarity is false, drawn objects remove previously drawn features.
     */
    bool polarity;

    /**
     * Whether apertures should be mirrored on the X axis.
     */
    bool ap_mirror_x;

    /**
     * Whether apertures should be mirrored on the Y axis.
     */
    bool ap_mirror_y;

    /**
     * Current rotation for apertures, in radians counter-clockwise from the
     * positive X axis.
     */
    double ap_rotate;

    /**
     * Current scale factor for apertures.
     */
    double ap_scale;

    /**
     * When set, we're inside a G36/G37 region block. Segments drawn via D01
     * commands are added to region_accum instead of being rendered as lines.
     */
    bool region_mode;

    /**
     * Path accumulator for the region command. When we're inside a G36/G37
     * block, segments drawn via D01 are added to the rather than being
     * rendered. When a G37 or D02 is encountered, the path is immediately
     * committed to the current plot and the accumulator is cleared via
     * commit_region().
     */
    CL::Path region_accum;

    /**
     * Accumulates all linear or circular interpolations drawn onto the topmost
     * plot via D01 commands, be they inside a region or not. After the Gerber
     * has been read, the contents of this can be used to look for loops, which
     * may then be used to reconstruct board outline and milling data when the
     * Gerber in question is the board outline and/or milling layer.
     */
    CL::Paths outline;

    /**
     * Render the current aperture to the current plot, taking into
     * consideration all configured aperture transformations.
     */
    void draw_aperture() {
        if (!aperture) {
            throw std::runtime_error("flash command before aperture set");
        }
        plot_stack.back()->draw_plot(
            aperture->get_plot(),
            polarity,
            pos.X, pos.Y,
            ap_mirror_x, ap_mirror_y,
            ap_rotate, ap_scale
        );
    }

    /**
     * Render an interpolation of the current aperture to the current plot, or
     * add the interpolation to the current region if we're inside a G36/G37
     * block.
     */
    void interpolate(CPt dest, CPt center) {

        // Interpolate a path from current position (pos) to dest using the
        // current interpolation mode.
        CL::Path path;
        if (imode == InterpolationMode::UNDEFINED) {
            throw std::runtime_error("interpolate command before mode set");
        } else if (imode == InterpolationMode::LINEAR) {

            // Handle linear interpolations. This is easy.
            path = {pos, dest};

        } else {

            // Handle circular interpolation. This is so esoteric we need a
            // helper class.
            std::shared_ptr<CircularInterpolationHelper> h;
            bool ccw = imode == InterpolationMode::CIRCULAR_CCW;

            if (qmode == QuadrantMode::UNDEFINED) {
                throw std::runtime_error("arc command before quadrant mode set");
            } else if (qmode == QuadrantMode::MULTI) {

                // Multi-quadrant mode is kind of sane, if over-constrained. The
                // over-constrainedness is solved by linearly interpolating
                // radius as well as angle.
                h = std::make_shared<CircularInterpolationHelper>(
                    CPt(pos.X, pos.Y),      // Start coordinate.
                    CPt(dest.X, dest.Y),    // End coordinate.
                    CPt(                    // Center point.
                        pos.X + center.X,
                        pos.Y + center.Y
                    ),
                    ccw, true               // Direction and mode.
                );

            } else {

                // Single-quadrant mode is esoteric-language-level insane. The
                // signs of the centerpoint offset are missing, so we have to
                // evaluate all four candidates, discard those that result in a
                // 90+ degree angle, and then choose the one with the least
                // radius difference from start to end. I'm sure this sounded
                // like a good idea at the time.
                for (unsigned int k = 0; k < 4; k++) {
                    auto h2 = std::make_shared<CircularInterpolationHelper>(
                        CPt(pos.X, pos.Y),   // Start coordinate.
                        CPt(dest.X, dest.Y), // End coordinate.
                        CPt(                 // Center point.
                            pos.X + ((k&1u) ? center.X : -center.X),
                            pos.Y + ((k&2u) ? center.Y : -center.Y)
                        ),
                        ccw, false          // Direction and mode.
                    );
                    if (h2->is_single_quadrant()) {
                        if (!h || h->error() > h2->error()) {
                            h = h2;
                        }
                    }
                }

            }
            if (!h) {
                throw std::runtime_error("failed to make circular interpolation");
            }
            path = h->to_path(fmt.get_max_deviation());

        }

        // Push all interpolation paths to outline, so we can reconstruct board
        // outline and/or milling data from such layers after processing.
        if (polarity && plot_stack.size() == 1) {
            outline.push_back(path);
        }

        // Handle region mode (G36/G37). Instead of drawing lines, we add to the
        // current region. We skip the start point of each path, as it matches
        // the end point of the previous path (we could equivalently have
        // decided to skip the end, and it would have worked the same).
        if (region_mode) {
            region_accum.insert(region_accum.end(), std::next(path.begin()), path.end());
            return;
        }

        // We only support interpolation for circular apertures; the Gerber
        // spec is a bit vague about whether different apertures are even legal
        // (I suspect they used to be, but were deprecated at some point).
        if (!aperture) {
            throw std::runtime_error("interpolate command before aperture set");
        }
        CInt diameter;
        if (!aperture->is_simple_circle(&diameter)) {
            throw std::runtime_error("only simple circle apertures without a hole are supported for interpolation");
        }
        double thickness = diameter * ap_scale;
        if (thickness == 0) {
            return;
        }

        // Use Clipper to add thickness to the path.
        CL::Paths ps;
        auto co = fmt.build_clipper_offset();
        co.AddPath(path, CL::jtRound, CL::etOpenRound);
        co.Execute(ps, thickness * 0.5);

        // Add the path to the plot.
        plot_stack.back()->draw_paths(ps, polarity);

    }

    /**
     * Commits the region path accumulated in region_accum via a G36/G37 block
     * to the current plot.
     */
    void commit_region() {

        // Check region size.
        if (region_accum.empty()) {
            return;
        }
        if (region_accum.size() < 3) {
            throw std::runtime_error("encountered region with less than 3 vertices");
        }

        // Gerber file paths can have any winding direction, but in order to
        // reduce clipper calls we want all paths to have positive winding if
        // we can help it. So if it's negative, reverse it.
        if (!CL::Orientation(region_accum)) {
            CL::ReversePath(region_accum);
        }

        // Render the region with the current polarity.
        plot_stack.back()->draw_paths({{region_accum}}, polarity);

        // Clear the accumulator to prepare for the next region.
        region_accum.clear();

    }

    /**
     * Handles a Gerber command. Returns true to continue, false if the command
     * marks the end, or throws a runtime error if the command is not
     * recognized.
     */
    bool command(const std::string &cmd, bool is_attrib) {
        if (am_builder) {

            // Handle aperture macro subcommands.
            am_builder->append(cmd);
            return true;

        } else if (is_attrib) {

            // Format specification.
            if (cmd.rfind("FS", 0) == 0) {
                if (cmd.size() != 10 || cmd.substr(2, 3) != "LAX" || cmd.substr(7, 1) != "Y" || cmd.substr(5, 2) != cmd.substr(8, 2)) {
                    throw std::runtime_error("invalid or deprecated and unsupported format specification: " + cmd);
                }
                fmt.configure_format(std::stoi(cmd.substr(5, 1)), std::stoi(cmd.substr(6, 1)));
                return true;
            }
            if (cmd.rfind("MO", 0) == 0) {
                if (cmd.substr(2, 2) == "IN") {
                    fmt.configure_inch();
                } else if (cmd.substr(2, 2) == "MM") {
                    fmt.configure_mm();
                } else {
                    throw std::runtime_error("invalid unit specification: " + cmd);
                }
                return true;
            }

            // Aperture definition.
            if (cmd.rfind("AD", 0) == 0) {
                if (cmd.size() < 3 || cmd.at(2) != 'D') {
                    throw std::runtime_error("invalid aperture definition: " + cmd);
                }
                size_t i = 3;
                size_t start = i;
                while (i < cmd.size()) {
                    if (!std::isdigit(cmd.at(i))) {
                        break;
                    }
                    i++;
                }
                int index = std::stoi(cmd.substr(start, i - start));
                if (index < 10) {
                    throw std::runtime_error("aperture index out of range: " + cmd);
                }
                std::vector<std::string> csep;
                start = i;
                while (i < cmd.size()) {
                    if (cmd.at(i) == ',' || (!csep.empty() && cmd.at(i) == 'X')) {
                        csep.push_back(cmd.substr(start, i - start));
                        start = i + 1;
                    }
                    i++;
                }
                csep.push_back(cmd.substr(start, i - start));
                if (csep.empty()) {
                    throw std::runtime_error("invalid aperture definition: " + cmd);
                }
                if (csep.at(0) == "C") {
                    apertures[index] = std::make_shared<CircleAperture>(csep, fmt);
                } else if (csep.at(0) == "R") {
                    apertures[index] = std::make_shared<RectangleAperture>(csep, fmt);
                } else if (csep.at(0) == "O") {
                    apertures[index] = std::make_shared<ObroundAperture>(csep, fmt);
                } else if (csep.at(0) == "P") {
                    apertures[index] = std::make_shared<PolygonAperture>(csep, fmt);
                } else {
                    auto it = aperture_macros.find(csep.at(0));
                    if (it == aperture_macros.end()) {
                        throw std::runtime_error("unsupported aperture type: " + csep.at(0));
                    }
                    apertures[index] = it->second->build(csep, fmt);
                }
                return true;
            }
            if (cmd.rfind("AM", 0) == 0) {
                auto name = cmd.substr(2);
                am_builder = std::make_shared<aperture_macro::ApertureMacro>();
                aperture_macros[name] = am_builder;
                return true;
            }
            if (cmd.rfind("AB", 0) == 0) {
                if (cmd == "AB") {
                    if (plot_stack.size() <= 1) {
                        throw std::runtime_error("unmatched aperture block close command");
                    }
                    plot_stack.pop_back();
                } else {
                    if (cmd.size() < 4 || cmd.at(2) != 'D') {
                        throw std::runtime_error("invalid aperture definition: " + cmd);
                    }
                    int index = std::stoi(cmd.substr(3));
                    if (index < 10) {
                        throw std::runtime_error("aperture index out of range: " + cmd);
                    }
                    auto plot = std::make_shared<Plot>();
                    plot_stack.push_back(plot);
                    apertures[index] = std::make_shared<CustomAperture>(plot);
                }
                return true;
            }

            // Unsupported attributes that can be ignored.
            if (cmd.rfind("TF", 0) == 0) {
                return true;
            }

            // Load polarity.
            if (cmd.rfind("LP", 0) == 0) {
                if (cmd.size() != 3 || !(cmd.at(2) == 'C' || cmd.at(2) == 'D')) {
                    throw std::runtime_error("invalid polarity command: " + cmd);
                }
                polarity = cmd.at(2) == 'D';
                return true;
            }

            // Aperture transformation commands.
            if (cmd == "LMN") {
                ap_mirror_x = false;
                ap_mirror_y = false;
                return true;
            }
            if (cmd == "LMX") {
                ap_mirror_x = true;
                ap_mirror_y = false;
                return true;
            }
            if (cmd == "LMY") {
                ap_mirror_x = false;
                ap_mirror_y = true;
                return true;
            }
            if (cmd == "LMXY") {
                ap_mirror_x = true;
                ap_mirror_y = true;
                return true;
            }
            if (cmd.rfind("LR", 0) == 0) {
                ap_rotate = std::stod(cmd.substr(2)) / 180.0 * M_PI;
                return true;
            }
            if (cmd.rfind("LS", 0) == 0) {
                ap_scale = std::stod(cmd.substr(2));
                return true;
            }

        } else {

            // Comment and deprecated commands with no effect.
            if (cmd.rfind("G04", 0) == 0) {
                return true;
            }
            if (cmd == "G54") {
                return true;
            }
            if (cmd == "G55") {
                return true;
            }

            // Interpolation mode.
            if (cmd == "G01") {
                imode = InterpolationMode::LINEAR;
                return true;
            }
            if (cmd == "G02") {
                imode = InterpolationMode::CIRCULAR_CW;
                return true;
            }
            if (cmd == "G03") {
                imode = InterpolationMode::CIRCULAR_CCW;
                return true;
            }
            if (cmd == "G74") {
                qmode = QuadrantMode::SINGLE;
                return true;
            }
            if (cmd == "G75") {
                qmode = QuadrantMode::MULTI;
                return true;
            }

            // Aperture selection.
            if (cmd.rfind('D', 0) == 0 && cmd.rfind("D0", 0) != 0) {
                auto it = apertures.find(std::stoi(cmd.substr(1)));
                if (it == apertures.end()) {
                    throw std::runtime_error("undefined aperture selected");
                } else {
                    aperture = it->second;
                    return true;
                }
            }

            // Move/draw commands.
            if (cmd.at(0) == 'X' || cmd.at(0) == 'Y' || cmd.at(0) == 'I' || cmd.at(0) == 'D') {
                std::map<char, CInt> params = {
                    {'X', pos.X},
                    {'Y', pos.Y},
                    {'I', 0},
                    {'J', 0},
                };
                int d = -1;
                char code = 0;
                size_t start = 0;
                for (size_t i = 0; i <= cmd.size(); i++) {
                    char c = (i < cmd.size()) ? cmd.at(i) : 'Z';
                    if (i == cmd.size() || isalpha(c)) {
                        if (code == 'D') {
                            d = std::stoi(cmd.substr(start, i - start));
                        } else if (code) {
                            params[code] = fmt.parse_fixed(cmd.substr(start, i - start));
                        }
                        code = c;
                        start = i + 1;
                    }
                }
                switch (d) {
                    case 1: // interpolate
                        interpolate(
                            {params.at('X'), params.at('Y')},
                            {params.at('I'), params.at('J')}
                        );
                        pos.X = params.at('X');
                        pos.Y = params.at('Y');
                        break;
                    case 2: // move
                        if (region_mode) {
                            commit_region();
                        }
                        pos.X = params.at('X');
                        pos.Y = params.at('Y');
                        break;
                    case 3: // flash
                        if (region_mode) {
                            throw std::runtime_error("cannot flash in region mode");
                        }
                        pos.X = params.at('X');
                        pos.Y = params.at('Y');
                        draw_aperture();
                        break;
                    default:
                        throw std::runtime_error("invalid draw/move command: " + std::to_string(d));
                }
                return true;
            }

            // Region mode.
            if (cmd == "G36") {
                if (region_mode) {
                    throw std::runtime_error("already in region mode");
                }
                region_mode = true;
                return true;
            }
            if (cmd == "G37") {
                if (!region_mode) {
                    throw std::runtime_error("not in region mode");
                }
                commit_region();
                region_mode = false;
                return true;
            }

            // Deprecated format specification; support anyway.
            if (cmd == "G70") {
                fmt.configure_inch();
                return true;
            }
            if (cmd == "G71") {
                fmt.configure_mm();
                return true;
            }
            if (cmd == "G90") {
                return true;
            }
            if (cmd == "G91") {
                throw std::runtime_error("incremental mode is not supported");
            }

            // Program stop.
            if (cmd == "M00") {
                return false;
            }
            if (cmd == "M01") {
                return false;
            }
            if (cmd == "M02") {
                return false;
            }

        }
        throw std::runtime_error("unknown command: " + cmd);
    }

    void end_attrib() {
        am_builder.reset();
    }

public:
    explicit GerberFile(const std::string &fname) {
        std::ifstream f(fname);
        if (!f.is_open()) {
            throw std::runtime_error("file not found");
        }
        f.exceptions(std::ifstream::badbit);
        imode = InterpolationMode::UNDEFINED;
        qmode = QuadrantMode::UNDEFINED;
        pos = {0, 0};
        polarity = true;
        ap_mirror_x = false;
        ap_mirror_y = false;
        ap_rotate = 0.0;
        ap_scale = 1.0;
        plot_stack = {std::make_shared<Plot>()};
        region_mode = false;

        bool is_attrib = false;
        std::ostringstream ss{};
        while (!f.eof()) {
            char c;
            f.get(c);
            if (isspace(c)) {
                continue;
            } else if (c == '%') {
                if (!ss.str().empty()) throw std::runtime_error("attribute mid-command");
                if (is_attrib) end_attrib();
                is_attrib = !is_attrib;
            } else if (c == '*') {
                if (ss.str().empty()) throw std::runtime_error("empty command");
                if (!command(ss.str(), is_attrib)) {
                    break;
                }
                ss = std::ostringstream();
            } else {
                ss << c;
            }
        }
        if (is_attrib) throw std::runtime_error("unterminated attribute");

        std::ofstream svg("kek.svg");
        svg << R"(<svg width="10000" height="10000" xmlns="http://www.w3.org/2000/svg">))" << std::endl;
        int color = 0;
        for (const auto &p : plot_stack.back()->get_dark()) {
            if (CL::Orientation(p)) {
                color++;
                int r = 128 + 127*sin(color);
                int g = 128 + 127*sin(color + 2*M_PI / 3);
                int b = 128 + 127*sin(color + 4*M_PI / 3);
                svg << "<path fill=\"rgb("<<r<<","<<g<<","<<b<<")\" fill-opacity=\"0.7\" d=\"";
            } else {
                svg << R"(<path fill="white" fill-opacity="0.7" d=")";
            }
            svg << "M " << fmt.to_mm(p.back().X) * 25 + 5000
                << " " << fmt.to_mm(-p.back().Y) * 25 + 5000 << " ";
            for (const auto &c : p) {
                svg << "L " << fmt.to_mm(c.X) * 25 + 5000
                    << " " << fmt.to_mm(-c.Y) * 25 + 5000 << " ";
            }
        svg << R"("/>)" << std::endl;
        }
        svg << R"(</svg>)" << std::endl;
    }
};

int main(int argc, char *argv[]) {
    GerberFile("/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/fts/fts.GBL");
}
