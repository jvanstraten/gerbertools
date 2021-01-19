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
 * Defines the top-level class for parsing Gerber files.
 */

#include <cmath>
#include <vector>
#include "gerber.hpp"
#include "path.hpp"

namespace gerber {

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
    CircularInterpolationHelper(
        coord::CPt start,
        coord::CPt end,
        coord::CPt center,
        bool ccw, bool multi
    ) {
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
    coord::Path to_path(double epsilon) const {
        double r = (r1 + r2) * 0.5;
        double x = (r > epsilon) ? (1.0 - epsilon / r) : 0.0;
        double th = std::acos(2.0 * x * x - 1.0) + 1e-3;
        auto n_vertices = (size_t)std::ceil(std::abs(a2 - a1) / th);
        coord::Path p;
        for (size_t i = 0; i <= n_vertices; i++) {
            double f2 = (double)i / (double)n_vertices;
            double f1 = 1.0 - f2;
            double vr = f1*r1 + f2*r2;
            double va = f1*a1 + f2*a2;
            double vx = center_x + vr * std::cos(va);
            double vy = center_y + vr * std::sin(va);
            p.push_back({(coord::CInt)std::round(vx), (coord::CInt)std::round(vy)});
        }
        return p;
    }

};

/**
 * Render the current aperture to the current plot, taking into
 * consideration all configured aperture transformations.
 */
void Gerber::draw_aperture() {
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
void Gerber::interpolate(coord::CPt dest, coord::CPt center) {

    // Interpolate a path from current position (pos) to dest using the
    // current interpolation mode.
    coord::Path path;
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
                coord::CPt(pos.X, pos.Y),   // Start coordinate.
                coord::CPt(dest.X, dest.Y), // End coordinate.
                coord::CPt(                 // Center point.
                    pos.X + center.X,
                    pos.Y + center.Y
                ),
                ccw, true                   // Direction and mode.
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
                    coord::CPt(pos.X, pos.Y),   // Start coordinate.
                    coord::CPt(dest.X, dest.Y), // End coordinate.
                    coord::CPt(                 // Center point.
                        pos.X + ((k&1u) ? center.X : -center.X),
                        pos.Y + ((k&2u) ? center.Y : -center.Y)
                    ),
                    ccw, false                  // Direction and mode.
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
    coord::CInt diameter;
    if (!aperture->is_simple_circle(&diameter)) {
        throw std::runtime_error("only simple circle apertures without a hole are supported for interpolation");
    }
    double thickness = diameter * ap_scale;
    if (thickness == 0) {
        return;
    }

    // Use Clipper to add thickness to the path.
    coord::Paths paths = path::render(path, thickness, fmt);

    // Add the path to the plot.
    plot_stack.back()->draw_paths(paths, polarity);

}

/**
 * Commits the region path accumulated in region_accum via a G36/G37 block
 * to the current plot.
 */
void Gerber::commit_region() {

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
    if (!ClipperLib::Orientation(region_accum)) {
        ClipperLib::ReversePath(region_accum);
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
bool Gerber::command(const std::string &cmd, bool is_attrib) {
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
                apertures[index] = std::make_shared<aperture::Circle>(csep, fmt);
            } else if (csep.at(0) == "R") {
                apertures[index] = std::make_shared<aperture::Rectangle>(csep, fmt);
            } else if (csep.at(0) == "O") {
                apertures[index] = std::make_shared<aperture::Obround>(csep, fmt);
            } else if (csep.at(0) == "P") {
                apertures[index] = std::make_shared<aperture::Polygon>(csep, fmt);
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
                auto plot = std::make_shared<plot::Plot>();
                plot_stack.push_back(plot);
                apertures[index] = std::make_shared<aperture::Custom>(plot);
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
            std::map<char, coord::CInt> params = {
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

/**
 * Handles the end of an attribute.
 */
void Gerber::end_attrib() {
    am_builder.reset();
}

/**
 * Loads a gerber file from the given stream. This reads until the end
 * command; the stream may not be EOF if there is more data at the end of
 * the stream.
 */
Gerber::Gerber(std::istream &s) {
    imode = InterpolationMode::UNDEFINED;
    qmode = QuadrantMode::UNDEFINED;
    pos = {0, 0};
    polarity = true;
    ap_mirror_x = false;
    ap_mirror_y = false;
    ap_rotate = 0.0;
    ap_scale = 1.0;
    plot_stack = {std::make_shared<plot::Plot>()};
    region_mode = false;
    outline_constructed = false;

    bool terminated = false;
    bool is_attrib = false;
    std::ostringstream ss{};
    while (!s.eof()) {
        char c;
        s.get(c);
        if (isspace(c)) {
            continue;
        } else if (c == '%') {
            if (!ss.str().empty()) throw std::runtime_error("attribute mid-command");
            if (is_attrib) end_attrib();
            is_attrib = !is_attrib;
        } else if (c == '*') {
            if (ss.str().empty()) throw std::runtime_error("empty command");
            if (!command(ss.str(), is_attrib)) {
                terminated = true;
                break;
            }
            ss = std::ostringstream();
        } else {
            ss << c;
        }
    }
    if (is_attrib) {
        throw std::runtime_error("unterminated attribute");
    }
    if (!terminated) {
        throw std::runtime_error("unterminated gerber file");
    }
    if (plot_stack.size() != 1) {
        throw std::runtime_error("unterminated block aperture");
    }
    if (region_mode) {
        throw std::runtime_error("unterminated region block");
    }
}

/**
 * Returns the paths representing the Gerber file.
 */
const coord::Paths &Gerber::get_paths() const {
    return plot_stack.back()->get_dark();
}

/**
 * Attempts to interpret the Gerber file data as the board outline and/or
 * milling data, returning polygons that follow the center of closed loops
 * of traces/regions in the file rather than the traces themselves. This is
 * a bit sensitive to round-off error and probably not work right if the
 * file isn't a proper outline; your mileage may vary.
 */
const coord::Paths &Gerber::get_outline_paths() {

    // Return immediately when the outline has already been constructed.
    if (outline_constructed) {
        return outline;
    }

    // Make a list of path indices that end in a particular coordinate for each
    // coordinate (within tolerance).
    using EndPoints = std::shared_ptr<std::list<size_t>>;
    std::map<std::pair<double, double>, EndPoints> point_map;
    std::set<EndPoints> points;
    std::vector<std::pair<EndPoints, EndPoints>> end_points;
    double eps = fmt.get_max_deviation();
    double eps_sqr = eps * eps;
    for (const auto &path : outline) {
        std::pair<EndPoints, EndPoints> endpts;
        for (int endpt = 0; endpt < 2; endpt++) {
            auto c_fix = endpt ? path.back() : path.front();
            auto c = std::make_pair<double, double>(c_fix.X, c_fix.Y);

            // Look for exact match first.
            auto it = point_map.find(c);
            if (it == point_map.end()) {

                // Nope. Try inexact.
                auto start = point_map.lower_bound({c.first - eps, c.second - eps});
                auto end = point_map.upper_bound({c.first + eps, c.second + eps});
                double least_error_sqr = INFINITY;
                for (auto it2 = start; it2 != end; ++it2) {
                    auto c2 = it2->first;
                    double dx = c.first - c2.first;
                    double dy = c.second - c2.second;
                    double error_sqr = dx*dx + dy*dy;
                    if (error_sqr < eps_sqr && error_sqr < least_error_sqr) {
                        least_error_sqr = error_sqr;
                        it = it2;
                    }
                }
                if (it == point_map.end()) {

                    // Nope. Add a new record for this point.
                    EndPoints ep = std::make_shared<std::list<size_t>>();
                    it = point_map.insert({c, ep}).first;
                    points.insert(ep);

                }
            }

            // Record the path as starting or ending in the point we found or
            // constructed.
            it->second->push_back(end_points.size());

            // Also record the other direction; from path to either point.
            if (endpt == 0) {
                endpts.first = it->second;
            } else {
                endpts.second = it->second;
            }
        }
        end_points.push_back(endpts);
    }

    coord::Paths paths;
    while (!points.empty()) {
        auto cur = *(points.begin());

        // Any point with something other than two paths ending in it is
        // non-manifold and will be ignored.
        if (cur->size() != 2) {
            points.erase(cur);
            continue;
        }

        // We have a potential start point. Now we just need to see if it's a
        // valid cycle or not. Even if it isn't, we can remove all the points
        // we find from the set.
        bool is_loop = true;
        coord::Path path;
        size_t start_index = cur->front();
        size_t cur_index = cur->back();
        while (true) {
            points.erase(cur);
            auto endpts = end_points.at(cur_index);
            auto section = outline.at(cur_index);
            if (endpts.first == cur) {
                path.insert(path.end(), section.begin(), std::prev(section.end()));
                cur = endpts.second;
            } else if (endpts.second == cur) {
                path.insert(path.end(), section.rbegin(), std::prev(section.rend()));
                cur = endpts.first;
            } else {
                throw std::runtime_error("this should never happen");
            }
            if (cur_index == start_index) {
                break;
            }
            if (cur->size() != 2) {
                is_loop = false;
                break;
            }
            if (cur->front() == cur_index) {
                cur_index = cur->back();
            } else if (cur->back() == cur_index) {
                cur_index = cur->front();
            } else {
                throw std::runtime_error("this should never happen");
            }
        }

        // If this is a valid loop, add the polygon.
        if (is_loop) {
            if (!ClipperLib::Orientation(path)) {
                ClipperLib::ReversePath(path);
            }
            paths.push_back(std::move(path));
        }

    }

    // Simplify the polygons with the correct fill rule.
    ClipperLib::SimplifyPolygons(paths, ClipperLib::pftEvenOdd);

    // Cache and return the construct outline.
    outline_constructed = true;
    outline = paths;
    return outline;
}

} // namespace gerber
