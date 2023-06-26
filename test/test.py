#!/usr/bin/env python

import unittest
import os
import sys
import subprocess
import IMP
import pickle


TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..',
                                      'rnapolii'))

class Tests(unittest.TestCase):

    def test_complete(self):
        """Test modeling and analysis"""
        # Run modeling
        os.chdir(os.path.join(TOPDIR, 'modeling'))
        p = subprocess.check_call(["python", 'modeling.py', "--test"])
        self.assertTrue(os.path.exists('output/rmfs/0.rmf3'))

        # Run clustering
        os.chdir(os.path.join(TOPDIR, 'analysis'))
        p = subprocess.check_call(["python", 'clustering.py', "--test"])
        self.assertTrue(os.path.exists('kmeans_5_1/dist_matrix.pdf'))
        self.assertTrue(os.path.exists('kmeans_5_1/cluster.0/0.rmf3'))

        # Test analysis
        p = subprocess.check_call(["python", 'precision_rmsf.py'])
        self.assertTrue(os.path.exists('kmeans_5_1/precision.0.0.out'))

        p = subprocess.check_call(["python", 'accuracy.py'])

    def test_pickle(self):
        """Test that pickled ReplicaExchange object works"""
        # Run modeling
        os.chdir(os.path.join(TOPDIR, 'modeling'))
        with open('modeling.py') as fh:
            contents = fh.read().replace('mc1.execute_macro()', '')
        g = {}
        exec(contents, g)
        mc1 = g['mc1']
        mc1.vars['number_of_frames'] = 5

        dump = pickle.dumps((mc1.model, mc1))

        # Run the original ReplicaExchange and get the final score
        IMP.random_number_generator.seed(99)
        mc1.execute_macro()
        rs = IMP.pmi.tools.get_restraint_set(mc1.model)
        original_score = rs.evaluate(False)
        del mc1, rs

        # With the same random seed, we should get the exact same trajectory
        # with the pickled object
        newm, newmc1 = pickle.loads(dump)
        IMP.random_number_generator.seed(99)
        newmc1.execute_macro()
        rs = IMP.pmi.tools.get_restraint_set(newmc1.model)
        new_score = rs.evaluate(False)
        self.assertAlmostEqual(original_score, new_score, delta=1e-4)


if __name__ == '__main__':
    unittest.main()
