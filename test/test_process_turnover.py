
import unittest, os, time
from brian2.units import second
import numpy as np

# get tested
from methods.process_turnover_pd import extract_lifetimes


class Test_extract_lifetimes(unittest.TestCase):

    def test_undying_synapses_added_to_wthsrv_but_not_dthonly(self):
        t_cut, Tmax = 10*second, 100*second
        turnover_data = np.array([[1, 20.85, 0, 0]])

        lts_wthsrv, lts_dthonly, _ = extract_lifetimes(turnover_data, 10,
                                                    t_cut, Tmax)

        self.assertEqual(lts_wthsrv, [Tmax/second-20.85])
        self.assertEqual(lts_dthonly, [])


    def test_synases_below_grown_pre_tcut_are_not_inclded(self):
        t_cut, Tmax = 10*second, 100*second
        turnover_data = np.array([[1, 2.85, 0, 0],
                                  [0, 7.85, 0, 0],
                                  [0, 0.85, 1, 0],
                                  [1, -50.85, 0, 12]])

        lts_wthsrv, lts_dthonly, _ = extract_lifetimes(turnover_data, 10,
                                                    t_cut, Tmax)

        self.assertEqual(lts_wthsrv, [])
        self.assertEqual(lts_dthonly, [])


    def test_synpase_first_died_started_but_did_not_die_again(self):
        t_cut, Tmax = 10*second, 100*second
        turnover_data = np.array([[0, 25.85, 0, 0],
                                  [1, 32.85, 0, 0]])

        lts_wthsrv, lts_dthonly, _ = extract_lifetimes(turnover_data, 10,
                                                    t_cut, Tmax)

        self.assertEqual(lts_wthsrv, [Tmax/second-32.85])
        self.assertEqual(lts_dthonly, [])


    def test_even_number_growth_death_events(self):
        t_cut, Tmax = 10*second, 100*second
        turnover_data = np.array([[0, 25.85, 0, 0],
                                  [1, 32.85, 0, 0],
                                  [0, 38.85, 0, 0],
                                  [1, 82.85, 0, 0],
                                  [0, 88.85, 0, 0]])

        lts_wthsrv, lts_dthonly, _ = extract_lifetimes(turnover_data, 10,
                                                    t_cut, Tmax)

        self.assertEqual(lts_wthsrv, [6.,6.])
        self.assertEqual(lts_dthonly, [6.,6.])


    def test_odd_number_growth_death_events(self):
        t_cut, Tmax = 10*second, 100*second
        turnover_data = np.array([[0, 25.85, 0, 0],
                                  [1, 32.85, 0, 0],
                                  [0, 38.85, 0, 0],
                                  [1, 82.85, 0, 0],
                                  [0, 88.85, 0, 0],
                                  [1, 91.02, 0, 0]])

        lts_wthsrv, lts_dthonly, _ = extract_lifetimes(turnover_data, 10,
                                                    t_cut, Tmax)

        np.testing.assert_array_almost_equal(lts_wthsrv, [6.,6.,8.98])
        np.testing.assert_array_almost_equal (lts_dthonly, [6.,6.])


    # def test_turnover_data_set1(self):
    #     t_split, t_cut, bin_w = 5*second, 2*second, 1*second
    #     turnover_data = [[1, 2.5, 0, 0]]

    #     full_t, ex_ids = extract_survival(np.array(turnover_data),
    #                                       bin_w, 10, t_split, t_cut)

    #     self.assertEqual(len(full_t),1)
    #     self.assertEqual(full_t[0],t_split/second)
        

    # def test_turnover_data_set2(self):
    #     t_split, t_cut, bin_w = 5*second, 2*second, 1*second
    #     turnover_data = [[1, 2.5, 0, 0],
    #                      [0, 2.7, 0, 0],
    #                      [1, 3.8, 0, 0],
    #                      [0, 4.592, 0, 0]]

    #     full_t, ex_ids = extract_survival(np.array(turnover_data),
    #                                       bin_w, 10, t_split, t_cut)

    #     np.testing.assert_array_almost_equal(full_t, [0.2, 0.792])

        
    # def test_turnover_data_set3(self):
    #     t_split, t_cut, bin_w = 5*second, 2*second, 1*second
    #     turnover_data = [[0, 2.3, 0, 0],
    #                      [1, 2.5, 0, 0],
    #                      [0, 2.7, 0, 0],
    #                      [1, 3.8, 0, 0],
    #                      [0, 4.592, 0, 0]]

    #     full_t, ex_ids = extract_survival(np.array(turnover_data),
    #                                       bin_w, 10, t_split, t_cut)

    #     np.testing.assert_array_almost_equal(full_t, [0.2, 0.792])


    # def test_turnover_data_set4(self):
    #     t_split, t_cut, bin_w = 5*second, 2*second, 1*second
    #     turnover_data = [[1, 2.5, 0, 0],
    #                      [0, 2.7, 0, 0],
    #                      [1, 3.8, 0, 0],
    #                      [0, 8.6, 0, 0]]

    #     full_t, ex_ids = extract_survival(np.array(turnover_data),
    #                                       bin_w, 10, t_split, t_cut)

    #     # test correct lifetimes
    #     np.testing.assert_array_almost_equal(full_t, [0.2, 4.8])


    # def test_turnover_data_set5(self):
    #     t_split, t_cut, bin_w = 5*second, 2*second, 1*second
    #     turnover_data = [[1, 2.5, 0, 0],
    #                      [0, 2.7, 0, 0],
    #                      [1, 3.8, 0, 0],
    #                      [1, 8.6, 0, 0]]

    #     full_t, ex_ids = extract_survival(np.array(turnover_data),
    #                                       bin_w, 10, t_split, t_cut)

    #     # test excluded ids
    #     self.assertEqual(len(ex_ids), 1)
    #     self.assertEqual(ex_ids[0], 0)

        
    # def test_turnover_data_set6(self):
    #     t_split, t_cut, bin_w = 5*second, 2*second, 1*second
    #     turnover_data = [[1, 2.5, 0, 0],
    #                      [0, 2.7, 0, 0],
    #                      [1, 3.8, 0, 0]]

    #     full_t, ex_ids = extract_survival(np.array(turnover_data),
    #                                       bin_w, 10, t_split, t_cut)

    #     np.testing.assert_array_almost_equal(full_t, [0.2, t_split/second])

        
    # def test_turnover_data_set7(self):
    #     t_split, t_cut, bin_w = 5*second, 2*second, 1*second
    #     turnover_data = [[1, 5.9, 0, 0]]

    #     full_t, ex_ids = extract_survival(np.array(turnover_data),
    #                                       bin_w, 10, t_split, t_cut)

    #     np.testing.assert_array_almost_equal(full_t, [5.])


    # def test_turnover_data_set8(self):
    #     t_split, t_cut, bin_w = 5*second, 2*second, 1*second
    #     turnover_data = [[1, 7.9, 0, 0]]

    #     full_t, ex_ids = extract_survival(np.array(turnover_data),
    #                                       bin_w, 10, t_split, t_cut)

    #     np.testing.assert_array_almost_equal(full_t, [])
        
    # def test_turnover_data_speed(self):
    #     t_split, t_cut, bin_w = 4*second, 2*second, 1*second
    #     turnover_data = np.genfromtxt('test/test_sets/turnover_test_set1',
    #                                   delimiter=',')

    #     a = time.time()
    #     full_t, ex_ids = extract_survival(turnover_data,
    #                                       bin_w, 1000,
    #                                       t_split, t_cut)
    #     b = time.time()

    #     print('Test Set 1 took :', b-a, ' s')
        


if __name__ == '__main__':
    unittest.main()
