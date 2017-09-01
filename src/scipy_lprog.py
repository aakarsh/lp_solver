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


class SciPy:

    import sys
    import numpy as np
    import lprog

    def __init__(self,A,b,c,n,m):
        (self.A,self.b,self.c) =  (A,b,c)
        self.n = n
        self.m = m

    def solve(self,tolerance = global_tolerance):
        import sys
        import numpy as np
        import lprog
        linprog_res = lprog.linprog([ -x for x in self.c ],
                                    A_ub = self.A,
                                    b_ub = self.b,
                                    options  = { 'tol': tolerance },
                                    callback = lprog.linprog_verbose_callback)
        print("%s" % linprog_res)
        if linprog_res.status == 0:
            return (1,linprog_res.x)

