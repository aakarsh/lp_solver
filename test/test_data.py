#!/usr/bin/env python3.4

import glob
import os
import sys
import unittest
import itertools
import functools

# what is the best way to get
sys.path.append("../src/")

import matrix
import tableau
import simplex
import scipy_lprog
import main

class TestDataRunner(unittest.TestCase):

    test_directory =  "/home/aakarsh/src/MOOC/coursera/lp_solver/test_data"

    def assertValidBounded(self,fname,tolerance=1e-4):
        s = simplex.Simplex.parse(file=fname)
        anst,ansx=s.solve()
        assert anst == 0
        self.assertTrue(s.verify_bounds(ansx,tolerance))
        self.assertTrue(s.verify_scipy(tolerance))

    def assertAnswerType(self,expected,fname,predicates=[]):
        try:
            self.assertEqual(anst,expected)
            for p in predicates:
                self.assertTrue(p(s,anst,ansx))
        except AssertionError as err:
            print(err)
            if anst != expected:
                raise AssertionError("Checking %s expected %s got %s " \
                                     % (fname,
                                        simplex.Simplex.answer_type_str(expected),
                                        simplex.Simplex.answer_type_str(anst)))
            raise err

    def flatten(self,ls):
        return list(itertools.chain(*ls))

    def test_compare_algorithms(self):
        patterns = ["%s/bounded/*","%s/inf/*","%s/no/*"]
        ignore = lambda fname: not fname.endswith(".a")
        test_dirs  = [ p % TestDataRunner.test_directory for p in patterns ]
        test_files = self.flatten([filter(ignore,glob.iglob(d)) for d in test_dirs])
        debug=False

        failures_by_file = {}
        
        for tf in test_files:
            
            (A,b,c,n,m) = main.Reader.parse(file=tf)
            
            algorithms  = {"Tableau" : tableau.Tableau(A,b,c,n,m,debug),
                           "Simplex" : simplex.Simplex(A,b,c,n,m,debug),
                           "Scipy"   : scipy_lprog.SciPy(A,b,c,n,m,debug)}
            results  = {}
            anst,ansx = None,None
            failed = False

            for (name,algo) in algorithms.items():
                results[name] = algo.solve()

                if anst is None and ansx is None:
                    anst,ansx = results[name]
                    
                cur_t,cur_x   = results[name]
                
                if cur_t != anst:
                    failed = True

            if not tf in failures_by_file and failed:
                failures_by_file[tf] = results


        for (tf,results) in failures_by_file.items():
            tf_dir,tf_file = tf.split(os.path.sep)[-2:]
            print("%s/%s:" %(tf_dir,tf_file))
                    
            for key in results.keys():
                anst,ansx = results[key]
                print("| %10s | %10s | %20s |" % (key,anst,matrix.map_optional(float,ansx)))


if __name__ == '__main__':
    unittest.main()


# def test_unbounded_files(self):
#     pat = ("%s/inf/[0-9]*"% TestDataRunner.test_directory)
#     not_answer_p = lambda fname: not fname.endswith(".a")
#     for fname in filter(not_answer_p,glob.iglob(pat)):
#         self.assertAnswerType(1,fname)
#
# def test_nosolution_files(self):
#     pat = ("%s/no/[0-9]*" % TestDataRunner.test_directory)
#     not_answer_p = lambda fname: not fname.endswith(".a")
#     for fname in filter(not_answer_p,glob.iglob(pat)):
#         self.assertAnswerType(-1,fname)
#
# def test_bounded_solution(self):
#     pat = ("%s/bounded/[0-9]*" % TestDataRunner.test_directory)
#     not_answer_p = lambda fname: not fname.endswith(".a")
#     for fname in filter(not_answer_p,glob.iglob(pat)):
#         self.assertAnswerType(-1,fname)
