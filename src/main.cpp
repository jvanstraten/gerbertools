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

#include "gerbertools/pcb.hpp"

int main(int argc, char *argv[]) {
    (void)argc;
    (void)argv;

    /*pcb::CircuitBoard pcb(
        "/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/ecu-bottom/ecu-bottom", ".GKO", ".TXT"
        //"/mnt/d/Jeroen/hobby/house/monitoring/cv-controller/hw/vert/JLCPCB/vert", ".GM2", ".TXT"
        //"O100030117 10by10 Green 1.6mm HASL 10pcs/mcu", ".GKO", ".TXT"
        //"pcb", ".GM3", ""
    );
    pcb.add_mask_layer(".GBS", ".GBO");
    pcb.add_copper_layer(".GBL");
    pcb.add_substrate_layer();
    pcb.add_copper_layer(".GTL");
    pcb.add_mask_layer(".GTS", ".GTO");
    pcb.add_surface_finish();
    pcb.write_svg("kek.svg", 20, {color::MASK_WHITE, color::SILK_BLACK});*/


    /*pcb::CircuitBoard pcb("/mnt/d/Jeroen/hobby/clock/components/cross/output/", "gbr.GM4", "ncd.TXT", "", "gbr.GM3");
    pcb.add_mask_layer("gbr.GBS", "gbr.GBO");
    pcb.add_copper_layer("gbr.GBL");
    pcb.add_substrate_layer(0.2);
    pcb.add_copper_layer("gbr.G2");
    pcb.add_substrate_layer(1.0);
    pcb.add_copper_layer("gbr.G1");
    pcb.add_substrate_layer(0.2);
    pcb.add_copper_layer("gbr.GTL");
    pcb.add_mask_layer("gbr.GTS", "gbr.GTO");
    pcb.add_surface_finish();
    pcb.write_svg("kek.svg", 50);*/

    pcb::CircuitBoard pcb("/mnt/d/Jeroen/hobby/clock/gated-clock/", "kek.GM1", "kek.TXT");
    pcb.add_mask_layer("kek.GBS", "kek.GBO");
    pcb.add_copper_layer("kek.GBL");
    pcb.add_substrate_layer(0.2);
    pcb.add_copper_layer("kek.G2");
    pcb.add_substrate_layer(1.0);
    pcb.add_copper_layer("kek.G1");
    pcb.add_substrate_layer(0.2);
    pcb.add_copper_layer("kek.GTL");
    pcb.add_mask_layer("kek.GTS", "kek.GTO");
    pcb.add_surface_finish();
    pcb.write_svg("/mnt/d/Jeroen/hobby/clock/gated-clock/kek.svg", 20, {color::MASK_WHITE, color::SILK_BLACK});


    /*std::ifstream f("/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/fts/fts.GBL");
    //std::ifstream f("O100030117 10by10 Green 1.6mm HASL 10pcs/mcu.GTL");
    if (!f.is_open()) {
        throw std::runtime_error("file not found");
    }
    f.exceptions(std::ifstream::badbit);
    auto gerber = gerber::Gerber(f);

    std::ifstream f2("/mnt/e/git/DARE subrepos/projects/stratos2plus/orders/2015-06-16/fts/fts.TXT");
    //std::ifstream f2("O100030117 10by10 Green 1.6mm HASL 10pcs/mcu.TXT");
    if (!f.is_open()) {
        throw std::runtime_error("file not found");
    }
    f2.exceptions(std::ifstream::badbit);
    auto drill = ncdrill::NCDrill(f2);

    ClipperLib::Clipper c;
    c.StrictlySimple(true);
    c.AddPaths(gerber.get_paths(), ClipperLib::ptSubject, true);
    c.AddPaths(drill.get_paths(), ClipperLib::ptClip, true);
    coord::Paths paths;
    c.Execute(ClipperLib::ctDifference, paths);

    std::ofstream svg("kek.svg");
    svg << R"(<svg width="10000" height="10000" xmlns="http://www.w3.org/2000/svg">))" << std::endl;
    int color = 0;
    for (const auto &p : paths) {
        if (ClipperLib::Orientation(p)) {
            color++;
            int r = 128 + 127*std::sin(color);
            int g = 128 + 127*std::sin(color + 2*M_PI / 3);
            int b = 128 + 127*std::sin(color + 4*M_PI / 3);
            svg << "<path fill=\"rgb("<<r<<","<<g<<","<<b<<")\" fill-opacity=\"0.7\" d=\"";
        } else {
            svg << R"(<path fill="white" fill-opacity="0.7" d=")";
        }
        svg << "M " << coord::Format::to_mm(p.back().X) * 25 + 5000
            << " " << coord::Format::to_mm(-p.back().Y) * 25 + 5000 << " ";
        for (const auto &c : p) {
            svg << "L " << coord::Format::to_mm(c.X) * 25 + 5000
                << " " << coord::Format::to_mm(-c.Y) * 25 + 5000 << " ";
        }
        svg << R"("/>)" << std::endl;
    }
    svg << R"(</svg>)" << std::endl;*/

}
