import unittest
import os

class TestClass(unittest.TestCase):

    def test_einsy_1_1(self):
        import gerbertools
        pcb = gerbertools.read('Einsy-Rambo-1.2b/Einsy Rambo_1.1a', 'GM88')
        os.makedirs('output', exist_ok=True)
        pcb.write_svg('output/rambo_1.1.svg', scale=20)
        self.assertTrue(os.path.isfile('output/rambo_1.1.svg'))

    def test_einsy_1_2(self):
        import gerbertools
        os.makedirs('output', exist_ok=True)
        result, code = gerbertools.main([
            'Einsy-Rambo-1.2b/Einsy Rambo_1.1a',
            '--outline', 'GM88',
            '--bounds',
            '--size',
            '--front', 'output/rambo_1.2.front.svg',
            '--back', 'output/rambo_1.2.back.svg',
            '--svg-scale', '20',
            '--obj', 'output/rambo_1.2.obj'
        ])
        self.assertEqual(code, 0)
        self.assertEqual(len(result), 2)
        self.assertEqual(result[0], '-22.485,-57.106,82.518,13.887')
        self.assertEqual(result[1], '105.004,70.993')
        self.assertTrue(os.path.isfile('output/rambo_1.2.front.svg'))
        self.assertTrue(os.path.isfile('output/rambo_1.2.back.svg'))
        self.assertTrue(os.path.isfile('output/rambo_1.2.obj'))

    def test_system76(self):
        import gerbertools
        result, code = gerbertools.main([
            'system76-launch-1.3/launch-',
            '--stackup', ';'.join([
                'Edge_Cuts.gbr',
                'PTH.drl:NPTH.drl',
                'm:B_Mask.gbr:B_SilkS.gbr',
                'c:B_Cu.gbr',
                's:0.2',
                'c:In2_Cu.gbr',
                's:1.1',
                'c:In1_Cu.gbr',
                's:0.2',
                'c:F_Cu.gbr',
                'm:F_Mask.gbr:F_SilkS.gbr'
            ]),
            '--bounds',
            '--size',
            '--front', 'output/system76.front.svg',
            '--back', 'output/system76.back.svg',
            '--svg-scale', '20',
            '--obj', 'output/system76.obj',
            '--drc', 'test.csv'
        ])
        self.assertEqual(code, 1)
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], '1.000,-114.000,300.000,13.000')
        self.assertEqual(result[1], '299.000,127.000')
        self.assertRegex(result[2], 'logical nets short[12] and short[12] are short-circuited')
        self.assertTrue(os.path.isfile('output/system76.front.svg'))
        self.assertTrue(os.path.isfile('output/system76.back.svg'))
        self.assertTrue(os.path.isfile('output/system76.obj'))

if __name__ == '__main__':
    unittest.main()
