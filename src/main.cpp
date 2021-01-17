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

namespace CL { using namespace ClipperLib; }
using CInt = CL::cInt;
using CPt = CL::IntPoint;

const CInt PRECISION_MULT = 0x1000000ll;
const CInt MITER_LIMIT = 0x10000000ll;
const CInt MAX_DEVIATION = 0x4000000ll;

class CoordUnit {
private:
    bool configured;
    double factor;
public:
    CoordUnit() : configured(false), factor(0.0) {}
    explicit CoordUnit(double factor) : configured(true), factor(factor) {}
    double unit_to_mm(double unit) const {
        if (!configured) {
            throw std::runtime_error("coordinate parsed before unit configured");
        }
        return unit * factor;
    }
    double mm_to_unit(double mm) const {
        if (!configured) {
            throw std::runtime_error("coordinate parsed before unit configured");
        }
        return mm / factor;
    }
};

class CoordFormat {
private:
    bool configured;
    int n_int;
    int n_dec;
public:
    CoordUnit unit;
    CoordFormat() : configured(false), n_int(0), n_dec(0) {}
    CoordFormat(int n_int, int n_dec) : configured(true), n_int(n_int), n_dec(n_dec) {}
    CInt parse_fixed(const std::string &s) const {
        if (!configured) {
            throw std::runtime_error("coordinate parsed before format configured");
        }
        return std::stoll(s) * PRECISION_MULT;
    }
    CInt parse_float(const std::string &s) const {
        if (!configured) {
            throw std::runtime_error("coordinate parsed before format configured");
        }
        return std::round(std::stod(s) * std::pow(10.0, n_dec) * PRECISION_MULT);
    }
    CInt to_fixed(double d) const {
        return std::round(d * std::pow(10.0, n_dec) * PRECISION_MULT);
    }
    double to_unit(CInt i) const {
        return i / (std::pow(10.0, n_dec) * PRECISION_MULT);
    }
    double miter_limit() const {
        return to_fixed(unit.mm_to_unit(0.005));
    }
    double max_deviation() const {
        return to_fixed(unit.mm_to_unit(0.005));
    }
};

enum class InterpolationMode {
    UNDEFINED, LINEAR, CIRCULAR_CW, CIRCULAR_CCW
};

enum class QuadrantMode {
    UNDEFINED, SINGLE, MULTI
};

/**
 * Represents a vector image, remembering which parts were explicitly cleared
 * as well as which parts are dark.
 */
class Plot {
private:
    mutable CL::Paths accum_paths;
    mutable bool accum_polarity;
    mutable CL::Paths dark;
    mutable CL::Paths clear;
    mutable bool simplified = false;
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
    void simplify() const {
        if (simplified) return;
        CL::SimplifyPolygons(dark, CL::pftNonZero);
        CL::SimplifyPolygons(clear, CL::pftNonZero);
        simplified = true;
    }
public:
    explicit Plot(
        const CL::Paths &dark = {},
        const CL::Paths &clear = {}
    ) : accum_polarity(true), dark(dark), clear(clear) {
    }
    void draw_paths(const CL::Paths &ps, bool polarity = true) {
        if (ps.empty()) return;
        if (polarity != accum_polarity) commit_paths();
        accum_polarity = polarity;
        accum_paths.insert(accum_paths.end(), ps.begin(), ps.end());
    }
    void draw_paths(
        const CL::Paths &ps,
        bool polarity,
        CInt x,
        CInt y = 0,
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
                c.X = std::round(cx) + x;
                c.Y = std::round(cy) + y;
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
    void draw_plot(
        const Plot &plt,
        bool polarity = true,
        CInt x = 0,
        CInt y = 0,
        bool mirror_x = false,
        bool mirror_y = false,
        double rotate = 0.0,
        double scale = 1.0
    ) {
        draw_paths(plt.get_dark(), polarity, x, y, mirror_x, mirror_y, rotate, scale);
        draw_paths(plt.get_clear(), !polarity, x, y, mirror_x, mirror_y, rotate, scale);
    }
    const CL::Paths &get_dark() const {
        commit_paths();
        simplify();
        return dark;
    }
    const CL::Paths &get_clear() const {
        commit_paths();
        simplify();
        return clear;
    }
};

class Aperture {
protected:
    std::shared_ptr<Plot> plot;
public:
    const Plot &get_plot() const {
        return *plot;
    }
    virtual bool is_simple_circle(CInt *diameter) const {
        return false;
    }
};

class CustomAperture : public Aperture {
public:
    explicit CustomAperture(const std::shared_ptr<Plot> &data) {
        plot = data;
    }
};

class StandardAperture : public Aperture {
protected:
    CInt hole_diameter = 0;
    CL::Paths get_hole(const CoordFormat &fmt) const {
        CL::Paths ps;
        if (hole_diameter > 0.0) {
            CL::ClipperOffset co{fmt.miter_limit(), fmt.max_deviation()};
            co.AddPath({{0, 0}}, CL::jtRound, CL::etOpenRound);
            co.Execute(ps, hole_diameter * 0.5);
            CL::ReversePaths(ps);
        }
        return ps;
    }
};

class CircleAperture : public StandardAperture {
private:
    CInt diameter = 0;
public:
    explicit CircleAperture(const std::vector<std::string> &csep, const CoordFormat &fmt) {
        if (csep.size() < 2 || csep.size() > 3) {
            throw std::runtime_error("invalid circle aperture");
        }
        diameter = fmt.to_fixed(std::stod(csep.at(1)));
        hole_diameter = (csep.size() > 2) ? fmt.to_fixed(std::stod(csep.at(2))) : 0;

        CL::Paths ps;
        CL::ClipperOffset co{fmt.miter_limit(), fmt.max_deviation()};
        co.AddPath({{0, 0}}, CL::jtRound, CL::etOpenRound);
        co.Execute(ps, diameter * 0.5);
        auto hole = get_hole(fmt);
        ps.insert(ps.end(), hole.begin(), hole.end());
        plot = std::make_shared<Plot>(ps);
    }
    bool is_simple_circle(CInt *diameter) const override {
        if (hole_diameter > 0.0) return false;
        if (diameter) *diameter = this->diameter;
        return true;
    }
};

class RectangleAperture : public StandardAperture {
private:
    CInt x_size;
    CInt y_size;
public:
    explicit RectangleAperture(const std::vector<std::string> &csep, const CoordFormat &fmt) {
        if (csep.size() < 3 || csep.size() > 4) {
            throw std::runtime_error("invalid rectangle aperture");
        }
        x_size = fmt.to_fixed(std::abs(std::stod(csep.at(1))));
        y_size = fmt.to_fixed(std::abs(std::stod(csep.at(2))));
        hole_diameter = (csep.size() > 3) ? fmt.to_fixed(std::stod(csep.at(3))) : 0;

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

class ObroundAperture : public StandardAperture {
private:
    CInt x_size;
    CInt y_size;
public:
    explicit ObroundAperture(const std::vector<std::string> &csep, const CoordFormat &fmt) {
        if (csep.size() < 3 || csep.size() > 4) {
            throw std::runtime_error("invalid obround aperture");
        }
        x_size = fmt.to_fixed(std::abs(std::stod(csep.at(1))));
        y_size = fmt.to_fixed(std::abs(std::stod(csep.at(2))));
        hole_diameter = (csep.size() > 3) ? fmt.to_fixed(std::stod(csep.at(3))) : 0;

        CInt x = x_size / 2;
        CInt y = y_size / 2;
        CInt r = std::min(x, y);
        x -= r;
        y -= r;
        CL::Paths ps;
        CL::ClipperOffset co{fmt.miter_limit(), fmt.max_deviation()};
        co.AddPath({{-x, -y}, {x, y}}, CL::jtRound, CL::etOpenRound);
        co.Execute(ps, r);
        auto hole = get_hole(fmt);
        ps.insert(ps.end(), hole.begin(), hole.end());
        plot = std::make_shared<Plot>(ps);
    }
};

class PolygonAperture : public StandardAperture {
private:
    CInt diameter = 0;
    size_t n_vertices;
    double rotation;
public:
    explicit PolygonAperture(const std::vector<std::string> &csep, const CoordFormat &fmt) {
        if (csep.size() < 3 || csep.size() > 5) {
            throw std::runtime_error("invalid polygon aperture");
        }
        diameter = fmt.to_fixed(std::stod(csep.at(1)));
        n_vertices = std::stoul(csep.at(2));
        rotation = (csep.size() > 3) ? (std::stod(csep.at(3)) / 180.0 * M_PI) : 0.0;
        hole_diameter = (csep.size() > 4) ? fmt.to_fixed(std::stod(csep.at(4))) : 0;
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

class ApertureMacroExpression {
private:
    static std::string debug(
        const std::list<std::shared_ptr<ApertureMacroExpression>> &expr,
        std::list<std::shared_ptr<ApertureMacroExpression>>::iterator expr_begin,
        std::list<std::shared_ptr<ApertureMacroExpression>>::iterator expr_end
    );
    static std::shared_ptr<ApertureMacroExpression> reduce(
        std::list<std::shared_ptr<ApertureMacroExpression>> &expr,
        std::list<std::shared_ptr<ApertureMacroExpression>>::iterator expr_begin,
        std::list<std::shared_ptr<ApertureMacroExpression>>::iterator expr_end
    );
public:
    static std::shared_ptr<ApertureMacroExpression> parse(std::string expr);
    virtual double eval(const std::map<size_t, double> &vars) const = 0;
    virtual char get_token() const { return 0; }
    virtual std::string debug() const = 0;
};

class ApertureMacroLiteral : public ApertureMacroExpression {
private:
    double value;
public:
    explicit ApertureMacroLiteral(double value) : value(value) {
    }
    double eval(const std::map<size_t, double> &vars) const override {
        return value;
    }
    std::string debug() const override {
        return "<" + std::to_string(value) + ">";
    }
};

class ApertureMacroVariable : public ApertureMacroExpression {
private:
    size_t index;
public:
    explicit ApertureMacroVariable(size_t index) : index(index) {
    }
    size_t get_index() const {
        return index;
    }
    double eval(const std::map<size_t, double> &vars) const override {
        auto it = vars.find(index);
        if (it != vars.end()) {
            return it->second;
        } else {
            return 0.0;
        }
    }
    std::string debug() const override {
        return "<$" + std::to_string(index) + ">";
    }
};

class ApertureMacroUnary : public ApertureMacroExpression {
private:
    char oper;
    std::shared_ptr<ApertureMacroExpression> expr;
public:
    ApertureMacroUnary(
        char oper,
        const std::shared_ptr<ApertureMacroExpression> &expr
    ) : oper(oper), expr(expr) {
    }
    double eval(const std::map<size_t, double> &vars) const override {
        if (oper == '+') {
            return expr->eval(vars);
        } else if (oper == '-') {
            return -expr->eval(vars);
        } else {
            throw std::runtime_error("invalid unary operator in aperture macro");
        }
    }
    std::string debug() const override {
        return "<" + std::string(1, oper) + expr->debug() + ">";
    }
};

class ApertureMacroBinary : public ApertureMacroExpression {
private:
    char oper;
    std::shared_ptr<ApertureMacroExpression> lhs;
    std::shared_ptr<ApertureMacroExpression> rhs;
public:
    ApertureMacroBinary(
        char oper,
        const std::shared_ptr<ApertureMacroExpression> &lhs,
        const std::shared_ptr<ApertureMacroExpression> &rhs
    ) : oper(oper), lhs(lhs), rhs(rhs) {
    }
    double eval(const std::map<size_t, double> &vars) const override {
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
    std::string debug() const override {
        return "<" + lhs->debug() + std::string(1, oper) + rhs->debug() + ">";
    }
};

class ApertureMacroToken : public ApertureMacroExpression {
private:
    char token;
public:
    ApertureMacroToken(char token) : token(token) {
    }
    double eval(const std::map<size_t, double> &vars) const override {
        throw std::runtime_error("cannot evaluate token");
    }
    char get_token() const override {
        return token;
    }
    std::string debug() const override {
        return "<" + std::string(1, token) + ">";
    }
};

std::string ApertureMacroExpression::debug(
    const std::list<std::shared_ptr<ApertureMacroExpression>> &expr,
    std::list<std::shared_ptr<ApertureMacroExpression>>::iterator expr_begin,
    std::list<std::shared_ptr<ApertureMacroExpression>>::iterator expr_end
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

std::shared_ptr<ApertureMacroExpression> ApertureMacroExpression::reduce(
    std::list<std::shared_ptr<ApertureMacroExpression>> &expr,
    std::list<std::shared_ptr<ApertureMacroExpression>>::iterator expr_begin,
    std::list<std::shared_ptr<ApertureMacroExpression>>::iterator expr_end
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
                *un = std::make_shared<ApertureMacroUnary>(t, *end);
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
                *lhs = std::make_shared<ApertureMacroBinary>(t, *lhs, *rhs);
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

std::shared_ptr<ApertureMacroExpression> ApertureMacroExpression::parse(
    std::string expr
) {

    // Tokenize input.
    expr += " ";
    std::list<std::shared_ptr<ApertureMacroExpression>> tokens;
    std::string cur;
    enum Mode {IDLE, LIT, VAR};
    Mode mode = IDLE;
    for (auto c : expr) {
        if (mode == LIT) {
            if ((c >= '0' && c <= '9') || (c == '.')) {
                cur += c;
                continue;
            } else {
                tokens.push_back(std::make_shared<ApertureMacroLiteral>(std::stod(cur)));
                cur.resize(0);
                mode = IDLE;
            }
        } else if (mode == VAR) {
            if (c >= '0' && c <= '9') {
                cur += c;
                continue;
            } else {
                tokens.push_back(std::make_shared<ApertureMacroVariable>(std::stoul(cur)));
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
            tokens.push_back(std::make_shared<ApertureMacroToken>(c));
        }
    }

    // Reduce expression.
    return reduce(tokens, tokens.begin(), tokens.end());
}

using ApertureMacroCommand = std::vector<std::shared_ptr<ApertureMacroExpression>>;

class ApertureMacro {
private:
    std::vector<ApertureMacroCommand> cmds;
public:
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
                ApertureMacroExpression::parse(cmd.substr(0, pos)),
                ApertureMacroExpression::parse(cmd.substr(pos + 1))
            });
        } else {
            size_t pos = 0;
            cmds.push_back({});
            while (true) {
                auto pos2 = cmd.find(',', pos);
                cmds.back().push_back(ApertureMacroExpression::parse(cmd.substr(pos, pos2 - pos)));
                if (pos2 == std::string::npos) {
                    break;
                }
                pos = pos2 + 1;
            }
        }
    }
    std::shared_ptr<Aperture> build(const std::vector<std::string> &csep, const CoordFormat &fmt) {
        std::map<size_t, double> vars;
        for (size_t i = 1; i < csep.size(); i++) {
            vars[i] = std::stod(csep.at(i));
        }
        Plot plot;
        for (const auto &cmd : cmds) {

            // Handle variable assignment commands.
            if (cmd.size() == 2) {
                auto v = std::dynamic_pointer_cast<ApertureMacroVariable>(cmd.at(0));
                if (v) {
                    vars[v->get_index()] = cmd.at(1)->eval(vars);
                    continue;
                }
            }

            // Handle draw commands.
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
                    CL::ClipperOffset co{fmt.miter_limit(), fmt.max_deviation()};
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
                    CL::ClipperOffset co{fmt.miter_limit(), fmt.max_deviation()};
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
                        CL::ClipperOffset co{fmt.miter_limit(), fmt.max_deviation()};
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

                    CL::ClipperOffset co1{fmt.miter_limit(), fmt.max_deviation()};
                    co1.AddPath({
                        {
                            fmt.to_fixed(center_x),
                            fmt.to_fixed(center_y)
                        }
                    }, CL::jtRound, CL::etOpenRound);
                    CL::Paths cps;
                    co1.Execute(cps, fmt.to_fixed(outer * 0.5));
                    ps.insert(ps.end(), cps.begin(), cps.end());

                    CL::ClipperOffset co2{fmt.miter_limit(), fmt.max_deviation()};
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

/**
 * Helper class for circular interpolations.
 */
class CircularInterpolationHelper {
private:
    double center_x, center_y;
    double r1, r2;
    double a1, a2;
    void to_polar(double x, double y, double &r, double &a) const {
        x -= center_x;
        y -= center_y;
        r = std::hypot(x, y);
        a = std::atan2(y, x);
    }
public:
    CircularInterpolationHelper(CInt x1, CInt y1, CInt x2, CInt y2, CInt xc, CInt yc, bool ccw, bool multi) {
        center_x = xc;
        center_y = yc;
        to_polar(x1, y1, r1, a1);
        to_polar(x2, y2, r2, a2);
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
    bool is_single_quadrant() const {
        return std::abs(a1 - a2) <= M_PI_2 + 1e-3;
    }
    double error() const {
        return std::max(r1, r2);
    }
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
    std::map<std::string, std::shared_ptr<ApertureMacro>> aperture_macros;

    /**
     * When non-null, an aperture macro is being constructed. Commands are
     * expected to be in the comma-separated aperture macro command syntax, to
     * be variable assignment statements of the form $<var>=<expr>, or to be
     * aperture macro comments. The attribute termination character (%) always
     * resets this to nullptr.
     */
    std::shared_ptr<ApertureMacro> am_builder;

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
                    pos.X, pos.Y,           // Start coordinate.
                    dest.X, dest.Y,         // End coordinate.
                    pos.X + center.X,       // Center point, X.
                    pos.Y + center.Y,       // Center point, Y.
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
                        pos.X, pos.Y,                            // Start coordinate.
                        dest.X, dest.Y,                          // End coordinate.
                        pos.X + ((k&1u) ? center.X : -center.X), // Center point, X.
                        pos.Y + ((k&2u) ? center.Y : -center.Y), // Center point, Y.
                        ccw, false                               // Direction and mode.
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
            path = h->to_path(fmt.max_deviation());

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
        CL::ClipperOffset co{fmt.miter_limit(), fmt.max_deviation()};
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
                CoordUnit unit = fmt.unit;
                fmt = CoordFormat(std::stoi(cmd.substr(5, 1)), std::stoi(cmd.substr(6, 1)));
                fmt.unit = unit;
                return true;
            }
            if (cmd.rfind("MO", 0) == 0) {
                if (cmd.substr(2, 2) == "IN") {
                    fmt.unit = CoordUnit(25.4);
                } else if (cmd.substr(2, 2) == "MM") {
                    fmt.unit = CoordUnit(1.0);
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
                am_builder = std::make_shared<ApertureMacro>();
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
                            d = std::stoull(cmd.substr(start, i - start));
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
                fmt.unit = CoordUnit(25.4);
                return true;
            }
            if (cmd == "G71") {
                fmt.unit = CoordUnit(1.0);
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
            svg << "M " << fmt.unit.unit_to_mm(fmt.to_unit(p.back().X)) * 25 + 5000
                << " " << fmt.unit.unit_to_mm(fmt.to_unit(-p.back().Y)) * 25 + 5000 << " ";
            for (const auto &c : p) {
                svg << "L " << fmt.unit.unit_to_mm(fmt.to_unit(c.X)) * 25 + 5000
                    << " " << fmt.unit.unit_to_mm(fmt.to_unit(-c.Y)) * 25 + 5000 << " ";
            }
        svg << R"("/>)" << std::endl;
        }
        svg << R"(</svg>)" << std::endl;
    }
};

int main(int argc, char *argv[]) {
    GerberFile("/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/fts/fts.GBL");
}
