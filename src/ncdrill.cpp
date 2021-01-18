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
 * Handles parsing NC drill files.
 */

#include <sstream>
#include <cmath>
#include "ncdrill.hpp"

namespace ncdrill {

/**
 * Constructs a new tool.
 */
Tool::Tool(coord::CInt diameter, bool plated) : diameter(diameter), plated(plated) {
}

/**
 * Returns the diameter of this tool.
 */
coord::CInt Tool::get_diameter() const {
    return diameter;
}

/**
 * Returns whether this tool is plated.
 */
bool Tool::is_plated() const {
    return plated;
}

/**
 * Commits the path in the path field to the plots and to vias based on the
 * current tool.
 */
void NCDrill::commit_path() {
    if (!tool) {
        throw std::runtime_error("tool use before any tool is selected");
    }
    if (tool->is_plated()) {
        plot_pth.draw_paths(plot::render_path(path, tool->get_diameter(), fmt));
        vias.insert(vias.end(), path.begin(), path.end());
    } else {
        plot_npth.draw_paths(plot::render_path(path, tool->get_diameter(), fmt));
    }
    path.clear();
}

/**
 * Adds an arc to the current path.
 */
void NCDrill::add_arc(
    coord::CPt start,
    coord::CPt end,
    coord::CInt radius,
    bool ccw
) {

    // Convert to double for the math.
    double x0 = start.X;
    double y0 = start.Y;
    double x1 = end.X;
    double y1 = end.Y;
    double r = radius;

    // Find the center.
    // https://math.stackexchange.com/a/88067
    double d = std::hypot(x0 - x1, y0 - y1);
    double e = 2.0 * r / d;
    if (e < 1.0) {
        e = 0.0;
    } else {
        e *= e;
        e = std::sqrt(e - 1.0);
    }
    if (!ccw) {
        e = -e;
    }
    double ax = (x0 - x1) * 0.5;
    double ay = (y0 - y1) * 0.5;
    double xc = ax + ay * e;
    double yc = ay - ax * e;

    // Find the start and end angle.
    double a0 = std::atan2(y0 - yc, x0 - xc);
    double a1 = std::atan2(y1 - yc, x1 - xc);
    if (ccw) {
        if (a1 < a0) {
            a1 += 2.0 * M_PI;
        }
    } else {
        if (a0 < a1) {
            a0 += 2.0 * M_PI;
        }
    }

    // Find the number of vertices we need to make.
    double epsilon = fmt.get_max_deviation();
    double f = (r > epsilon) ? (1.0 - epsilon / r) : 0.0;
    double th = std::acos(2.0 * f * f - 1.0) + 1e-3;
    auto n_vertices = (size_t)std::ceil(std::abs(a1 - a0) / th);

    // Actually add the path.
    for (size_t i = 1; i <= n_vertices; i++) {
        double f1 = (double)i / (double)n_vertices;
        double f0 = 1.0 - f1;
        double va = f0*a0 + f1*a1;
        double vx = xc + r * std::cos(va);
        double vy = yc + r * std::sin(va);
        path.push_back({(coord::CInt)std::round(vx), (coord::CInt)std::round(vy)});
    }

}

/**
 * Parses a regular command, consisting of concatenated letter-number pairs.
 * Throws a std::runtime_error if the command does not conform to that
 * syntax.
 */
std::map<char, std::string> NCDrill::parse_regular_command(const std::string &cmd) {
    std::map<char, std::string> params;
    char code = 0;
    size_t start = 0;
    for (size_t i = 0; i <= cmd.size(); i++) {
        char c = (i < cmd.size()) ? cmd.at(i) : '\0';
        if (i == cmd.size() || isalpha(c)) {
            if (i > 0) {
                if (i <= start) {
                    throw std::runtime_error("unknown/unexpected NC drill command: " + cmd);
                }
                params[code] = cmd.substr(start, i - start);
            }
            code = c;
            start = i + 1;
        }
    }
    return params;
}

/**
 * Processes a command. Returns whether processing is complete.
 */
bool NCDrill::command(const std::string &cmd) {

    // Ignore empty lines.
    if (cmd.empty()) {
        return true;
    }

    // Parse commands state-dependently.
    if (parse_state == ParseState::PRE_HEADER) {

        // Comment.
        if (cmd.at(0) == ';') {
            return true;
        }

        // Start of header.
        if (cmd == "M48") {
            parse_state = ParseState::HEADER;
            return true;
        }

    } else if (parse_state == ParseState::HEADER) {

        // Comment.
        if (cmd.at(0) == ';') {

            // Coordinate format comment as generated by Altium.
            if (cmd.size() == 16 &&cmd.substr(0, 13) == ";FILE_FORMAT=" && cmd.at(14) == ':') {
                fmt.configure_format(std::stoi(cmd.substr(13, 1)), std::stoi(cmd.substr(15, 1)));
            }

            // Plating types as generated by Altium.
            if (cmd == ";TYPE=PLATED") {
                plated = true;
            }
            if (cmd == ";TYPE=NON_PLATED") {
                plated = false;
            }

            return true;
        }

        // Unit & coordinate format.
        if (cmd.substr(0, 6) == "METRIC") {
            fmt.configure_mm();

            // Interpret leading/trailing zero data as generated by
            // Altium.
            if (cmd.substr(6) == ",LZ") {
                fmt.configure_trailing_zeros(true);
            } else {
                fmt.configure_trailing_zeros(false);
            }

            return true;
        }
        if (cmd.substr(0, 4) == "INCH") {
            fmt.configure_inch();

            // Interpret leading/trailing zero data as generated by
            // Altium.
            if (cmd.substr(4) == ",LZ") {
                fmt.configure_trailing_zeros(true);
            } else {
                fmt.configure_trailing_zeros(false);
            }

            return true;
        }

        // End of header.
        if (cmd == "%") {
            parse_state = ParseState::BODY;
            return true;
        }

        // Anything not parsed yet through the above exceptions should
        // be a regular command.
        auto params = parse_regular_command(cmd);

        // Handle tool definition commands.
        auto it = params.find('T');
        if (it != params.end()) {
            size_t tool_no = std::stoull(it->second);
            it = params.find('C');
            if (it == params.end()) {
                throw std::runtime_error("missing tool diameter in " + cmd);
            }
            tools[tool_no] = std::make_shared<Tool>(fmt.parse_float(it->second), plated);
            return true;
        }

    } else if (parse_state == ParseState::BODY) {

        // Ignore comments.
        if (cmd.at(0) == ';') {
            return true;
        }

        // Parse the command as a regular command.
        auto params = parse_regular_command(cmd);

        // Handle coordinates.
        coord::CPt start_point = pos;
        bool coord_set = false;
        auto it = params.find('X');
        if (it != params.end()) {
            pos.X = fmt.parse_fixed(it->second);
            coord_set = true;
        }
        it = params.find('Y');
        if (it != params.end()) {
            pos.Y = fmt.parse_fixed(it->second);
            coord_set = true;
        }
        coord::CPt end_point = pos;

        // Handle T (tool change) commands.
        it = params.find('T');
        if (it != params.end()) {
            size_t t = std::stoull(it->second);
            if (rout_mode == RoutMode::ROUT_TOOL_DOWN) {
                throw std::runtime_error("unexpected tool change; tool is down");
            }
            auto it2 = tools.find(t);
            if (it2 == tools.end()) {
                throw std::runtime_error("attempting to change to undefined tool");
            }
            tool = it2->second;
            return true;
        }

        // Handle G commands.
        it = params.find('G');
        if (it != params.end()) {
            int g = std::stoi(it->second);

            // Set rout mode.
            if (g == 0) {
                if (rout_mode != RoutMode::DRILL) {
                    throw std::runtime_error("unexpected G00; already routing");
                }
                rout_mode = RoutMode::ROUT_TOOL_UP;
                return true;
            }

            // Linear rout.
            if (g == 1) {
                if (rout_mode == RoutMode::DRILL) {
                    throw std::runtime_error("unexpected G01; not routing");
                } else if (rout_mode == RoutMode::ROUT_TOOL_DOWN) {
                    path.push_back(end_point);
                }
                return true;
            }

            // Circular clockwise rout.
            if (g == 2) {
                if (rout_mode == RoutMode::DRILL) {
                    throw std::runtime_error("unexpected G02; not routing");
                } else if (rout_mode == RoutMode::ROUT_TOOL_DOWN) {
                    it = params.find('A');
                    if (it == params.end()) {
                        throw std::runtime_error("arc radius is missing for G02");
                    }
                    add_arc(start_point, end_point, fmt.parse_fixed(it->second), false);
                }
                return true;
            }

            // Circular counterclockwise rout.
            if (g == 3) {
                if (rout_mode == RoutMode::DRILL) {
                    throw std::runtime_error("unexpected G03; not routing");
                } else if (rout_mode == RoutMode::ROUT_TOOL_DOWN) {
                    it = params.find('A');
                    if (it == params.end()) {
                        throw std::runtime_error("arc radius is missing for G03");
                    }
                    add_arc(start_point, end_point, fmt.parse_fixed(it->second), true);
                }
                return true;
            }

            // Enter drill mode.
            if (g == 5) {
                if (rout_mode == RoutMode::ROUT_TOOL_DOWN) {
                    throw std::runtime_error("unexpected G05; cannot exit route mode with tool down");
                } else if (rout_mode == RoutMode::DRILL) {
                    throw std::runtime_error("unexpected G05; already in drill mode");
                }
                rout_mode = RoutMode::DRILL;
                return true;
            }

            throw std::runtime_error("unsupported G command: " + cmd);
        }

        // Handle M commands.
        it = params.find('M');
        if (it != params.end()) {
            int m = std::stoi(it->second);

            // Routing tool down.
            if (m == 15) {
                if (rout_mode == RoutMode::ROUT_TOOL_DOWN) {
                    throw std::runtime_error("unexpected M15; tool already down");
                } else if (rout_mode == RoutMode::DRILL) {
                    throw std::runtime_error("unexpected M15; not in rout mode");
                }
                rout_mode = RoutMode::ROUT_TOOL_DOWN;
                path.push_back(end_point);
                return true;
            }

            // Routing tool up.
            if (m == 16) {
                if (rout_mode == RoutMode::ROUT_TOOL_UP) {
                    throw std::runtime_error("unexpected M16; tool already up");
                } else if (rout_mode == RoutMode::DRILL) {
                    throw std::runtime_error("unexpected M16; not in rout mode");
                }
                rout_mode = RoutMode::ROUT_TOOL_UP;
                commit_path();
                return true;
            }

            // End of file.
            if (m == 30) {
                if (rout_mode == RoutMode::ROUT_TOOL_DOWN) {
                    throw std::runtime_error("end of file with routing tool down");
                }
                return false;
            }

        }

        // Remaining commands should be tool hits.
        if (coord_set) {
            path.push_back(end_point);
            commit_path();
            return true;
        }

    }

    throw std::runtime_error("unknown/unexpected command: " + cmd);
}

/**
 * Parses an NC drill file.
 */
NCDrill::NCDrill(std::istream &s, bool default_plated) {
    parse_state = ParseState::PRE_HEADER;
    plated = default_plated;
    fmt.configure_format(4, 3);
    fmt.configure_mm();
    pos = {0, 0};
    rout_mode = RoutMode::DRILL;

    bool terminated = false;
    std::ostringstream ss{};
    while (!s.eof()) {
        char c;
        s.get(c);
        if (c == '\n') {
            terminated = !command(ss.str());
            if (terminated) break;
            ss = std::ostringstream();
        } else if (isspace(c)) {
            continue;
        } else {
            ss << c;
        }
    }
    if (!terminated) {
        throw std::runtime_error("unterminated NC drill file");
    }
}

/**
 * Returns the cutout paths for this NC drill file as negatively-wound
 * polygons.
 */
plot::Paths NCDrill::get_paths(bool plated, bool unplated) const {
    plot::Paths paths;
    if (plated) {
        if (unplated) {
            ClipperLib::Clipper c;
            c.AddPaths(plot_pth.get_dark(), ClipperLib::ptSubject, true);
            c.AddPaths(plot_npth.get_dark(), ClipperLib::ptSubject, true);
            c.Execute(ClipperLib::ctUnion, paths, ClipperLib::pftPositive);
        } else {
            paths = plot_pth.get_dark();
        }
    } else {
        if (unplated) {
            paths = plot_npth.get_dark();
        } else {
            return {};
        }
    }
    ClipperLib::ReversePaths(paths);
    return paths;
}

/**
 * Returns the center coordinates of all plated holes.
 */
const std::list<coord::CPt> &NCDrill::get_vias() const {
    return vias;
}

} // namespace ncdrill
