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
 * Handles all the complexity related to aperture macros.
 */

#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <memory>
#include <list>
#include <map>
#include <vector>
#include <sstream>
#include "gerbertools/aperture_macro.hpp"
#include "gerbertools/path.hpp"

namespace aperture_macro {

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
 * If this is a character token, return the character it represents.
 * Otherwise retun '\0'.
 */
char Expression::get_token() const {
    return 0;
}

/**
 * Constructs a literal node with the given value.
 */
LiteralExpression::LiteralExpression(double value) : value(value) {
}

/**
 * Evaluate this expression with the given set of variables.
 */
double LiteralExpression::eval(const Variables &vars) const {
    return value;
}

/**
 * Returns a debug representation of this expression node.
 */
std::string LiteralExpression::debug() const {
    return "<" + std::to_string(value) + ">";
}

/**
 * Constructs a variable reference node with the given variable index.
 */
VariableExpression::VariableExpression(size_t index) : index(index) {
}

/**
 * Constructs a variable node for the given variable index.
 */
size_t VariableExpression::get_index() const {
    return index;
}

/**
 * Evaluate this expression with the given set of variables.
 */
double VariableExpression::eval(const Variables &vars) const {
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
std::string VariableExpression::debug() const {
    return "<$" + std::to_string(index) + ">";
}

/**
 * Constructs a unary expression node from the given operand character
 * (either '+' or '-') and operand expression.
 */
UnaryExpression::UnaryExpression(
    char oper,
    const ExpressionRef &expr
) : oper(oper), expr(expr) {
}

/**
 * Evaluate this expression with the given set of variables.
 */
double UnaryExpression::eval(const Variables &vars) const {
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
std::string UnaryExpression::debug() const {
    return "<" + std::string(1, oper) + expr->debug() + ">";
}

/**
 * Constructs a binary expression node from the given operand character
 * (either '+', '-', 'x', or '/') and operand expressions.
 */
BinaryExpression::BinaryExpression(
    char oper,
    const ExpressionRef &lhs,
    const ExpressionRef &rhs
) : oper(oper), lhs(lhs), rhs(rhs) {
}

/**
 * Evaluate this expression with the given set of variables.
 */
double BinaryExpression::eval(const Variables &vars) const {
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
std::string BinaryExpression::debug() const {
    return "<" + lhs->debug() + std::string(1, oper) + rhs->debug() + ">";
}

/**
 * Constructs a token character. The character must be '(', ')', '+', '-',
 * 'x', or '/'.
 */
Token::Token(char token) : token(token) {
}

/**
 * Evaluate this expression with the given set of variables.
 */
double Token::eval(const Variables &vars) const {
    throw std::runtime_error("cannot evaluate token");
}

/**
 * Returns the token's character representation.
 */
char Token::get_token() const {
    return token;
}

/**
 * Returns a debug representation of this expression node.
 */
std::string Token::debug() const {
    return "<" + std::string(1, token) + ">";
}

/**
 * Parses and appends an aperture macro command.
 */
void ApertureMacro::append(const std::string &cmd) {
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
aperture::Ref ApertureMacro::build(const std::vector<std::string> &csep, const coord::Format &fmt) {
    Variables vars;
    for (size_t i = 1; i < csep.size(); i++) {
        vars[i] = std::stod(csep.at(i));
    }
    plot::Plot plot;
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
                auto paths = path::render({
                    {
                        fmt.to_fixed(center_x),
                        fmt.to_fixed(center_y)
                    }
                }, fmt.to_fixed(diameter), fmt);
                plot.draw_paths(
                    paths, exposure,
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
                auto paths = path::render({
                    {
                        fmt.to_fixed(start_x),
                        fmt.to_fixed(start_y)
                    }, {
                        fmt.to_fixed(end_x),
                        fmt.to_fixed(end_y)
                    }
                }, fmt.to_fixed(width), fmt, true);
                plot.draw_paths(
                    paths, exposure,
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
                coord::Paths paths = {{
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
                    paths, exposure,
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
                coord::Paths paths = {{}};
                for (size_t i = 0; i < n_vertices; i++) {
                    paths.back().push_back({
                        fmt.to_fixed(cmd.at(3 + 2*i)->eval(vars)),
                        fmt.to_fixed(cmd.at(4 + 2*i)->eval(vars))
                    });
                };
                plot.draw_paths(
                    paths, exposure,
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
                coord::Paths ps = {{}};
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
                coord::Paths paths = {};
                for (size_t i = 0; i < max_rings*2 && diameter > 0.0; i++) {
                    auto circle_paths = path::render({
                        {
                            fmt.to_fixed(center_x),
                            fmt.to_fixed(center_y)
                        }
                    }, fmt.to_fixed(diameter), fmt);
                    if (i & 1) {
                        ClipperLib::ReversePaths(circle_paths);
                        diameter -= gap * 2.0;
                    } else {
                        diameter -= thickness * 2.0;
                    }
                    paths.insert(paths.end(), circle_paths.begin(), circle_paths.end());
                };
                if (ch_thickness > 0.0 && ch_length > 0.0) {
                    paths.push_back({
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
                    paths.push_back({
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
                    paths, true,
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

                auto paths = path::render({
                    {
                        fmt.to_fixed(center_x),
                        fmt.to_fixed(center_y)
                    }
                }, fmt.to_fixed(outer), fmt);

                auto inner_paths = path::render({
                    {
                        fmt.to_fixed(center_x),
                        fmt.to_fixed(center_y)
                    }
                }, fmt.to_fixed(inner), fmt);
                ClipperLib::ReversePaths(inner_paths);
                paths.insert(paths.end(), inner_paths.begin(), inner_paths.end());

                if (gap > 0.0) {
                    paths.push_back({
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
                    paths.push_back({
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
                    paths, true,
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

    auto ap_plot = std::make_shared<plot::Plot>();
    ap_plot->draw_paths(plot.get_dark());
    return std::make_shared<aperture::Custom>(ap_plot);
}

} // namespace aperture_macro
