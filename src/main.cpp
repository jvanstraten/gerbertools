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
#include <cmath>
#include "gerber.hpp"

int main(int argc, char *argv[]) {
    //std::ifstream f("/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/fts/fts.GKO");
    std::ifstream f("pcb.GTO");
    if (!f.is_open()) {
        throw std::runtime_error("file not found");
    }
    f.exceptions(std::ifstream::badbit);
    auto gerber = gerber::Gerber(f);

    std::ofstream svg("kek.svg");
    svg << R"(<svg width="10000" height="10000" xmlns="http://www.w3.org/2000/svg">))" << std::endl;
    int color = 0;
    for (const auto &p : gerber.get_outline_paths()) {
        if (ClipperLib::Orientation(p)) {
            color++;
            int r = 128 + 127*std::sin(color);
            int g = 128 + 127*std::sin(color + 2*M_PI / 3);
            int b = 128 + 127*std::sin(color + 4*M_PI / 3);
            svg << "<path fill=\"rgb("<<r<<","<<g<<","<<b<<")\" fill-opacity=\"0.7\" d=\"";
        } else {
            svg << R"(<path fill="white" fill-opacity="0.7" d=")";
        }
        svg << "M " << gerber.to_mm(p.back().X) * 25 + 5000
            << " " << gerber.to_mm(-p.back().Y) * 25 + 5000 << " ";
        for (const auto &c : p) {
            svg << "L " << gerber.to_mm(c.X) * 25 + 5000
                << " " << gerber.to_mm(-c.Y) * 25 + 5000 << " ";
        }
        svg << R"("/>)" << std::endl;
    }
    svg << R"(</svg>)" << std::endl;

}
