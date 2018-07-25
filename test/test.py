#!/usr/bin/env python

import unittest
import os
import sys
import subprocess

TOPDIR = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), '..',
                                      '.'))

class Tests(unittest.TestCase):

    def test_complete(self):
        """Test modeling and analysis"""
        # Run modeling
        os.chdir(os.path.join(TOPDIR, 'modeling'))
        p = subprocess.check_call(["python", 'IMP_example.py'])

        p = subprocess.check_call(["python", 'modeling.py'])
        self.assertTrue(os.path.exists('output/rmfs/0.rmf3'))

        # Run short analysis
        p = subprocess.check_call(["python", 'short_analysis.py'])

        # Test long analysis
        p = subprocess.check_call(["python", 'long_analysis.py'])

if __name__ == '__main__':
    unittest.main()
