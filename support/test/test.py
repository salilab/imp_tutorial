#!/usr/bin/env python

import unittest
import os
import glob
import shutil
import subprocess
import tarfile
import pickle
import urllib.request
import IMP


TOPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..',
                                      '..', 'rnapolii'))
RESULTS = "https://salilab.org/ftp/tutorials/imp/rnapolii/results.tar.gz"
ANALYSIS = "https://salilab.org/ftp/tutorials/imp/rnapolii/analysis.tar.gz"


class Tests(unittest.TestCase):

    def clean_output(self):
        for subdir in ('pdbs', 'rmfs'):
            fp = os.path.join(TOPDIR, 'modeling', 'output', subdir)
            if os.path.exists(fp):
                shutil.rmtree(fp, ignore_errors=False)
        for name in ('best.scores.rex.py', 'initial.0.rmf3',
                     'stat.0.out', 'stat_replica.0.out'):
            fp = os.path.join(TOPDIR, 'modeling', 'output', name)
            if os.path.exists(fp):
                os.unlink(fp)
        for subdir in glob.glob(os.path.join(TOPDIR, 'analysis',
                                             'kmeans_*_*')):
            shutil.rmtree(subdir, ignore_errors=False)

    def test_analysis_on_precomputed(self):
        """Make sure that analysis works on precomputed results"""
        self.clean_output()
        os.chdir(os.path.join(TOPDIR, 'analysis'))
        # Get and extract precomputed results
        with urllib.request.urlopen(RESULTS) as fh:
            with tarfile.open(fileobj=fh, mode='r:gz') as tf:
                if hasattr(tarfile, 'data_filter'):
                    tf.extractall(filter='data')
                else:
                    tf.extractall()
        for subdir in os.listdir('results/output'):
            shutil.move('results/output/%s' % subdir, '../modeling/output/')

        # Run clustering
        subprocess.check_call(["python", 'clustering.py'])
        self.assertTrue(os.path.exists('kmeans_5_1/dist_matrix.pdf'))
        self.assertTrue(os.path.exists('kmeans_5_1/cluster.0/0.rmf3'))

        # Test analysis
        subprocess.check_call(["python", 'precision_rmsf.py'])
        self.assertTrue(os.path.exists('kmeans_5_1/precision.0.0.out'))

        subprocess.check_call(["python", 'accuracy.py'])

    def test_accuracy_of_precomputed_analysis(self):
        """Make sure that precomputed analysis can be checked for accuracy"""
        self.clean_output()
        os.chdir(os.path.join(TOPDIR, 'analysis'))
        # Get and extract precomputed analysis
        with urllib.request.urlopen(ANALYSIS) as fh:
            with tarfile.open(fileobj=fh, mode='r:gz') as tf:
                if hasattr(tarfile, 'data_filter'):
                    tf.extractall(filter='data')
                else:
                    tf.extractall()

        subprocess.check_call(["python", 'accuracy.py'])

    def test_complete(self):
        """Test modeling and analysis"""
        self.clean_output()
        # Run modeling
        os.chdir(os.path.join(TOPDIR, 'modeling'))
        subprocess.check_call(["python", 'modeling.py', "--test"])
        self.assertTrue(os.path.exists('output/rmfs/0.rmf3'))

        # Run clustering
        os.chdir(os.path.join(TOPDIR, 'analysis'))
        subprocess.check_call(["python", 'clustering.py', "--test"])
        self.assertTrue(os.path.exists('kmeans_5_1/dist_matrix.pdf'))
        self.assertTrue(os.path.exists('kmeans_5_1/cluster.0/0.rmf3'))

        # Test analysis
        subprocess.check_call(["python", 'precision_rmsf.py'])
        self.assertTrue(os.path.exists('kmeans_5_1/precision.0.0.out'))

        subprocess.check_call(["python", 'accuracy.py'])

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
