import unittest
from os import chdir

from mcsteel import mc


class CFASTScenarioTestCase(unittest.TestCase):
    def setUp(self):
        chdir('files/CaseCFAST')

    def test_unify_characters(self):
        self.assertEqual(
            mc.CFASTScenario._unify_characters('&MISC LOWER_OXYGEN_LIMIT = 0.15, MAX_TIME_STEP =  0.1/').split(),
            ['&MISC', 'LOWER_OXYGEN_LIMIT', '0.15', 'MAX_TIME_STEP', '0.1']
        )
        
    def test_get_fire_location(self):
        passages = [
            ["&FIRE ID   'f19' COMP_ID   'r19' FIRE_ID   'f19' LOCATION   5.63  6.7 ", [5.63, 6.7]],
            ["&FIRE ID   'f19' COMP_ID   'r19' LOCATION   5.63  6.7 FIRE_ID   'f19' ", [5.63, 6.7]],
        ]
        for p in passages:
            self.assertEqual(mc.CFASTScenario._get_fire_location(p[0]), p[1])

    def test_get_fire_id(self):
        passages = [
            ["&FIRE ID   'f19' COMP_ID   'r19' FIRE_ID   'f19' LOCATION   5.63  6.7 ", 'f19'],
            ["&FIRE ID   'f19' COMP_ID   'r19' LOCATION   5.63  6.7 FIRE_ID   'f19' ", 'f19'],
        ]
        for p in passages:
            self.assertEqual(mc.CFASTScenario._get_fire_id(p[0]), p[1])

    def test_get_hr_area_record(self):
        passages = [
            ["&TABL ID   'f19' DATA   30  3216  1.63  1.99  0.903  0.113  0.0 ", (3216, 1.99)],
            ["&TABL  DATA   30  3216  1.63  1.99  0.903  0.113  0.0 ID   'f19' ", (3216, 1.99)]
            ]
        for p in passages:
            self.assertEqual(mc.CFASTScenario._get_hrr_area_record(p[0]), p[1])

    def test_get_room_height(self):
        passages = [
            ["&COMP ID   'r4'  WIDTH   18.76  DEPTH   23.84  HEIGHT   2.5  CEILING_MATL_ID   'concrete'  WALL_MATL_ID   'concrete'  FLOOR_MATL_ID   'concrete'  CEILING_THICKNESS   0.3  WALL_THICKNESS   0.3  FLOOR_THICKNESS   0.3  ORIGIN   10.0  10.96  0.0  LEAK_AREA   0.00035  5.2e-05  GRID   50  50  50  ", ('r4', 2.5)],
            ["&COMP  WIDTH   18.76  DEPTH   23.84  HEIGHT   2.5  CEILING_MATL_ID  'concrete' ID   'r4'  WALL_MATL_ID   'concrete'  FLOOR_MATL_ID   'concrete'  CEILING_THICKNESS   0.3  WALL_THICKNESS   0.3  FLOOR_THICKNESS   0.3  ORIGIN   10.0  10.96  0.0  LEAK_AREA   0.00035  5.2e-05  GRID   50  50  50  ", ('r4', 2.5)],
        ]
        for p in passages:
            self.assertEqual(mc.CFASTScenario._get_hrr_area_record(p[0]), p[1])

    def test_calculate_hrrpua(self):
        passages = [
            [[(1000, 2), (2000, 2), (3000, 2)], 1500],
            [[(1000, 2), (2000, 2), (3000, 2), (2500, 2), (1500, 2)], 1500]
        ]
        for p in passages:
            self.assertEqual(mc.CFASTScenario._calculate_hrrpua(p[0]), p[1])

    def test_cfast_extract_in(self):
        cfast_scenario = mc.CFASTScenario(mc.Config('setup/test.user'), 'cfast')
        self.assertEqual(cfast_scenario._cfast_extract_in(), ((None, 1620.1, [5.63, 6.7] ), 2.5))

    def test_init(self):
        cfast_scenario = mc.CFASTScenario(mc.Config('setup/test.user'), 'cfast')
        self.assertEqual(cfast_scenario.hrrpua, 1620.1)
        self.assertEqual(cfast_scenario.alpha, None)
        self.assertEqual(cfast_scenario.fire_curve, [[], []])
        self.assertEqual(cfast_scenario.ceiling, 2.5)
        self.assertEqual(cfast_scenario.fire_type, 'cfast')
        self.assertEqual(cfast_scenario.sprinklers, None)

    def test_create_fire_curve(self):
        cfast_scenario = mc.CFASTScenario(mc.Config('setup/test.user'), 'cfast')
        cfast_scenario.create_fire_curve()
        # complex to test for assertEquals :(







class MultisimulationTestCase(unittest.TestCase):
    pass
