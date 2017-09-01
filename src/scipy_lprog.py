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

import sys
import numpy as np
import lprog

class SciPy:

    def __init__(self,A,b,c,n,m,debug=False):
        (self.A,self.b,self.c) =  (A,b,c)
        self.n = n
        self.m = m
        self.debug = debug

    def solve(self,tolerance = global_tolerance):

        linprog_res = lprog.linprog([ -x for x in self.c ],
                                    A_ub = self.A,
                                    b_ub = self.b,
                                    options  = { 'tol': tolerance },
                                    callback = lprog.linprog_verbose_callback if self.debug else None) 

        if self.debug:
            print("%s" % linprog_res)
            print("x--:%s" % linprog_res.x)
            print("status--:%s" % linprog_res.status)
        
        if linprog_res.status == 0:
            return (0,linprog_res.x)
        if linprog_res.status == 3: # unbounded
            return (1,None)
        if linprog_res.status == 2: # no solution
            return (-1,None)                
        else:
            return (linprog_res.status,linprog_res.x)
        

