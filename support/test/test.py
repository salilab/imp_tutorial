#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..',
                                      '..', 'rnapolii'))

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

if __name__ == '__main__':
    unittest.main()
