#!/usr/bin/env python3.4

import glob
import os
import sys
import unittest
import itertools
import functools
import simplex
from simplex import *

class SimplexTest(unittest.TestCase):

    def setUp(self): pass

    def test_01(self):
        debug=True
        A = [[-1 ,-1 ],
             [ 1 , 0 ],
             [ 0 , 1 ]]
        b = [-1 ,2 ,2]
        c = [-1 ,2]
        n = 3
        m = 2
        s = simplex.Simplex(A,b,c,n,m)
        opt,ass = s.solve()
        self.assertEqual(s.optimum,4.0)

    def test_feasible_start(self):
        A = [[1 ,1 ,3 ],
             [2 ,2 ,5 ],
             [4 ,1 ,2 ]]
        b = [30 ,24 ,36]
        c = [3 ,1 ,2]
        n = 3
        m = 3
        s = simplex.Simplex(A,b,c,n,m)
        opt,ass = s.solve()
        self.assertEqual(s.optimum,28)

    def test_infeasible_start(self):
        A = [[ 2 ,-1],
             [ 1 ,-5]]
        b = [2, -4]
        c = [2, -1]
        n = 2
        m = 2
        s = simplex.Simplex(A,b,c,n,m)
        opt,ass = s.solve()
        self.assertEqual(s.optimum,2.0)

    def test_read_01(self):
        s = simplex.Simplex.parse(file="./tests/bounded/01")
        self.assertIsNotNone(s)
        self.assertEqual(3,s.n)
        self.assertEqual(2,s.m)
        self.assertEqual(s.n,len(s.A))

    def test_read_04(self):
        s = simplex.Simplex.parse(file="./tests/bounded/04")
        anst,ansx = s.solve()
        if debug: print("%s"%ansx)

    def test_solve_01(self):
        s = simplex.Simplex.parse(file="./tests/bounded/01")
        anst,ansx = s.solve()
        if debug: print("anst:%s ansx:%s"%(anst,ansx))
        self.assertIsNotNone(anst)
        self.assertIsNotNone(ansx)
        
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

    def test_unbounded_files(self):
        pat = "./tests/inf/[0-9]*"
        not_answer_p = lambda fname: not fname.endswith(".a")
        for fname in filter(not_answer_p,glob.iglob(pat)):
            self.assertAnswerType(1,fname)
    
    def test_nosolution_files(self):
        pat = "./tests/no/[0-9]*"
        not_answer_p = lambda fname: not fname.endswith(".a")
        for fname in filter(not_answer_p,glob.iglob(pat)):
            self.assertAnswerType(-1,fname)

    # def test_tolerance_failure(self):
    #     fname = "./fail/11"
    #     simplex.debug=True            
    #     self.assertValidBounded(fname)


            
def run_tests():
    unittest.main(module='test_simplex',exit=False)

if __name__ == "__main__":
    if len(sys.argv) > 1 and (sys.argv[1] == "-d" or sys.argv[1] == "--debug") :
        simplex.debug = True    
    run_tests()
