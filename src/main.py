#!/usr/bin/env python3.4

import glob
import os
import sys
import unittest
import itertools
import functools
import argparse
import decimal
from decimal import Decimal
from itertools import chain

from matrix import *
from simplex import *
from scipy_lprog import *
from tableau import *

class Reader:
    
    @staticmethod
    def parse(stream=sys.stdin,file=None):
        if file:
            stream = open(file,"r")            
        with stream as stream:
            n, m = list(map(int, stream.readline().split()))
            A = []
            for i in range(n):
                A += [list(map(int, stream.readline().split()))]
            b = list(map(int, stream.readline().split()))
            c = list(map(int, stream.readline().split()))
            return (A,b,c,n,m)


if __name__ =="__main__":
    parser = argparse.ArgumentParser(prog="simplex.py",
                                     description='Run simplex on matrix from standard input.')
    parser.add_argument("-d","--debug",action="count",help="enable debug level log")
    parser.add_argument("-t","--tolerance",help="floating point tolerance to tolerate in an intolerable world")
    parser.add_argument("-s","--scipy",action='count',help="Use sci-py instead to compare answers")
    parser.add_argument("-v","--verify",action='count',help="Verify sci-py instead to compare answers")

    args = parser.parse_args()

    if args.debug:
        debug = True
    if args.tolerance:
        global_tolerance = float(args.tolerance)

    (A,b,c,n,m) = Reader.parse()

    if debug:
        print("-------------------- Tableau --------------------")
    tableau = Tableau(A,b,c,n,m,debug)
    t_anst, t_ansx =  tableau.solve()

    compare = False
    if compare:
        if debug: print("-------------------------------------------------")
        if debug: print("-------------------- Scipy --------------------")        
        scipy = SciPy(A,b,c,n,m,debug)    
        l_anst, l_ansx = scipy.solve()
        if debug: print("---------------------End:Scipy--------------------")

        if debug: print("-------------------- Simplex --------------------")
        simplex = Simplex(A,b,c,n,m,debug)
        anst, ansx = simplex.solve()
        if debug: print("-----------------------------------------------")

        if debug:
            print("---------------------------------------------------------------")
            solvers = [("Tableau", t_anst, list(map(float,t_ansx)) if t_ansx else None),
                       ("Simplex",anst,list(map(float,ansx)) if ansx else None),
                       ("Scipy",l_anst,l_ansx)]
            for s in solvers:
                print("%-15s | %15s | %4s " % s)
            print("---------------------------------------------------------------")
    

    def answer_type_str(anst):
        try:
            return ["No solution","Bounded solution","Infinity"][anst+1]
        except :
            return "Unrecognized Answer : %s" % anst

    print(answer_type_str(t_anst))

    if t_anst == 0:
        print(' '.join(list( map( lambda x : '%.18f' % x, t_ansx))))
        
        # simplex.verify_bounds(tolerance = global_tolerance)
        # if args.verify:
        #     simplex.verify_scipy(tolerance = global_tolerance)

    # if args.scipy:
    #     print(simplex.solve_scipy())

