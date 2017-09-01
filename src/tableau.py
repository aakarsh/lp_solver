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

# TODO: Remove unecessary features in order to keep it simple 
# TODO : Write tests for Infeasible and Unbounded solutions, other stuff here.

class Tableau:
    
    """ Tableau based solver which is based on lp solver in simplex. """
    
    def __init__(self,A,b,c,n,m):
        (self.A,self.b,self.c,
         self.n,self.m,self.v) = (A,b,c,n,m,0)
            
        self.v = 0
        v = self.v
            
        # what is difference between decimal and numpy 64
        c = list(map(Decimal,c))

        # A represents the table of size total x total with None possible entries
        nvars = n + m

        # One slack variable for each upper bound constraint
        self.n_slack = m

        # Artificial variables for every negative value of b.
        self.n_artificial = (list_wrapper(b) < 0).count_nonzero()

        self.T = Matrix.zeros([ m + 2 , self.n + self.n_slack + self.n_artificial + 1])

        # copy objective
        self.T[-2,:n]  = -list_wrapper(c)
        self.T[-2,-1]  = -v

        # copy constants
        self.T[0:-2,-1] = b

        # copy A
        self.T[0:m,:n]   = A

        self.T.fill_diagonal([slice(0,m,1),
                              slice(n,n + self.n_slack, 1)],1)

        if debug:
            print("After filling in the diagonals")
            print("T:")
            print(self.T)

        slcount = 0 # count of slack variables
        avcount = 0 # count of artificial

        self.basis = [0] * m
        artificial = [0] * self.n_artificial

        for i in range(m) :
            # Negative constant need to Introduce a artificial variables
            if self.T[i,-1] < 0 :
                # Namespace of artificial variables starts just beyond
                self.basis[i] = n + self.n_slack + avcount
                artificial[avcount] = i

                self.T[i,:-1]  *= -1
                self.T[i ,-1]  *= -1

                self.T[ i, self.basis[i]]  = 1
                self.T[-1, self.basis[i]]  = 1

                avcount += 1
            else:
                self.basis[i] = n + slcount
                slcount += 1

        if debug:
            print("T:")            
            print(self.T)
            print("basis variables: %s " % self.basis)

        for r in artificial:
            self.T[-1,:] = self.T[-1,:] - self.T[r,:]

        if debug:
            print("T:")
            print(self.T)



    def pivot_column(self,T,tol=global_tolerance):
        """Go through the objective row and find the minimum entry above
           tolerance"""
        # Ignore all positive values where: positive is defined as
        # anything greater than -tol
        ignored = lambda e: (e is None) or (e >= -tol)
        objective = T[-1,:-1].masked_where(ignored)
        return objective.min_index()

    def pivot_row(self,T,pivcol,tol,phase=1):
        """ Find the appropriate pivot row. """
        # Skip objective rows, in first phase we have two objective rows
        # including a psuedo-objective row.

        skip_rows = 2 if (phase == 1) else 1

        # Mask values less than tolerance
        ignored = lambda e: (e is None) or ( e <= tol )

        # print(">> T\n %s"%T)
        # print(">> T[:-skip_rows,pivcol] : \n %s",T[:-skip_rows,pivcol])
        # print(">> T[:-skip_rows] :\n %s",T[:-skip_rows])
        # only seem to be getting back a single integer
        # need to get all rows except the skipped rows

        ma = T[:-skip_rows,pivcol].masked_where(ignored)

        # All pivot column entries
        if ma.count_notnone() == 0 :
            return (False,None)

        mb = self.T[:-skip_rows,-1].masked_where(ignored)

        q = mb / ma
        print("q: %s" % q)
        return q.min_index()


    def do_pivot(self, T, pivrow, pivcol, basis, phase):
        "Perform the pivot operation using pivot row and pivot column and basis"

        if debug:
            print("Before Pivot[%d,%d] " % (pivrow,pivcol))
            print("T : \n%s" % T)
            print("Basic Variables : %s\n" % basis)

        basis[pivrow] = pivcol
        T[pivrow,:]   = T[pivrow,:] / T[pivrow][pivcol]

        for irow in range(T.shape()[0]):
            if irow != pivrow:
                T[irow,:] = T[irow,:] - T[pivrow,:] * T[irow,pivcol]

        if debug:
            print("After Pivot[%d,%d] " % (pivrow,pivcol))
            print("T : \n%s" % T)
            print("Basic Variables : %s\n" % basis)
            print("Current Objective Value:\n%s\n" % -T[-1,-1])


    def simplex_solve(self, T, n, basis, phase =2, tol = global_tolerance):
        # Ignore original and new objectives.
        if phase == 1:
            m = T.shape()[0] - 2
        elif phase == 2:
            m = T.shape()[0] - 1

        self.nvars = (T.shape()[1]-1)
        complete = False
        solution = [0] * self.nvars

        max_iterations = 5
        nit = 0

        # Identify and substitute
        if phase == 2:
            # Identify aritificial variables still in the objective
            ncols = T.shape()[1]
            is_artificial = lambda idx : basis[idx] > ncols - 2

            # Check basis for artificial variables
            variables = list(filter(is_artificial,range(len(basis))))

            if debug:
                print("Basic Variables : %s\n" % basis)
                print("Artificial Pivot Variables : %s " % variables)

            for pivrow in variables:
                non_zero_col = lambda col: self.T[pivrow,col] != 0
                pivcols = filter(non_zero_col,range(ncols -1))

                if len(pivcols) == 0: continue
                self.do_pivot(T,pivrow,pivcols[0],basis,phase)

        while not complete and (nit < max_iterations):

            nit += 1
            pivcol_found, pivcol = self.pivot_column(T,tol)

            if debug:
                print("T:")
                print(T)

            if not pivcol_found:  # Finished with all the columns, in basic form
                status, complete = 0, True
                break

            pivrow_found, pivrow = self.pivot_row(T,pivcol,tol)

            print("pivrow_found : %s pivrow: %s" % (pivrow_found, pivrow))

            if not pivrow_found: # Not finding the pivot row is very serious.
                status, complete = 3, True
                break

            if not complete: # perform the pivot on pivot entry
                self.do_pivot(T,pivrow,pivcol, basis,phase)

        if complete:
            print("Pivot Result == Phase: %d == " % phase)
            print("Basic Variables : %s\n" % basis)
            print("T : \n%s" % T)
        return (status,complete)


    def solve(self):
        # Pivot to basic flexible.
        status, complete = self.simplex_solve(self.T,self.n,self.basis,phase=1)

        print("Pseudo-Objecive : %d " % self.T[-1,-1])

        if abs(self.T[-1,-1]) < global_tolerance:
            # Remove pseudo-objective row
            self.T.del_rows([len(self.T)-1])
            # Remove artificial variable columns
            self.T.del_columns(range(self.n + self.n_slack,self.n + self.n_slack + self.n_artificial))
        else:
            # Infeasible without starting point
            status = 2
            print("Infeasible soltion")

        if status == 2: # infeasible solution
            raise InfeasibleError()

        sys.exit(-1)
        status, complete = self.simplex_solve(self.T,self.n,self.basis,phase=2)

        if status == 0:
            obj,ansx = -self.T[-1,-1], [0] * self.nvars
        else:
            raise Exception("Expected status == 0")

        print("found objective: %f" % obj)

        return (status,ansx)
