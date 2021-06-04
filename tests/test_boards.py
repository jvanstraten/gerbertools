import unittest
import os

class TestClass(unittest.TestCase):

    def test_einsy_1_1(self):
        import gerbertools
        pcb = gerbertools.read('Einsy-Rambo-1.2b/Einsy Rambo_1.1a', 'GM88')
        os.makedirs('output', exist_ok=True)
        pcb.write_svg('output/rambo_1.1.svg', scale=20)

    def test_einsy_1_2(self):
        import gerbertools
        pcb = gerbertools.read('Einsy-Rambo-1.2b/Einsy Rambo_1.2a', 'GM88')
        os.makedirs('output', exist_ok=True)
        pcb.write_svg('output/rambo_1.2.svg', scale=20)

    def test_system76(self):
        import gerbertools
        pcb = gerbertools.CircuitBoard('system76-launch-1.3/launch-', 'Edge_Cuts.gbr', 'PTH.drl', 'NPTH.drl');
        pcb.add_mask_layer('B_Mask.gbr', 'B_SilkS.gbr');
        pcb.add_copper_layer('B_Cu.gbr')
        pcb.add_substrate_layer(0.2)
        pcb.add_copper_layer('In2_Cu.gbr')
        pcb.add_substrate_layer(1.1)
        pcb.add_copper_layer('In1_Cu.gbr')
        pcb.add_substrate_layer(0.2)
        pcb.add_copper_layer('F_Cu.gbr')
        pcb.add_mask_layer('F_Mask.gbr', 'F_SilkS.gbr');
        pcb.add_surface_finish()
        os.makedirs('output', exist_ok=True)
        pcb.write_svg('output/system76.svg', scale=20)

if __name__ == '__main__':
    unittest.main()
