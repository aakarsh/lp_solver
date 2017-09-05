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


    def __init__(self,A,b,c,m,n,debug=False):

        """Called with A set of coefficieant for linear equations
        along with their upper bound constraints b and the
        coefficients for the objective function to maximize."""

        self.debug = debug

        (self.A,self.b,self.c,
         self.n,self.m,self.v) = (A,b,c,n,m,0)

        if self.debug:
            print("n: %d  m:%d" % (n,m))
            print("c: %s" % c)
            print("b: %s" % b)
            print("A:\n%s" % A)

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

        self.T = Matrix.zeros([ m + 2 ,
                                self.n + self.n_slack + self.n_artificial + 1])


        if self.debug:
            print("T :")
            print(self.T)

        # copy objective
        self.T[-2,:n]  = -list_wrapper(c)
        self.T[-2,-1]  = -v

        # copy constants
        self.T[0:-2,-1] = b

        # copy A
        self.T[0:m,:n]   = A

        if self.debug:
            print("copying coefficents:")
            print("b:%s" %b)
            print("T:")
            print(self.T)


        self.T.fill_diagonal([slice(0,m,1),
                              slice(n,n + self.n_slack, 1)],1)

        if self.debug:
            print("adding slack variables:")
            print("T:")
            print(self.T)

        slcount = 0 # count of slack variables
        avcount = 0 # count of artificial

        self.basis = [0] * m
        artificial = [0] * self.n_artificial

        for i in range(m) :
            # negative constant need to introduce a artificial variables
            # maybe this comparison is not safe.
            if self.T[i,-1] < 0 :
                # namespace of artificial variables starts just beyond
                
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

        if self.debug:
            print("adding artificial variables")
            print("artificial: %s "% artificial)
            print("T:")
            print(self.T)
            print("basis variables: %s " % self.basis)

        for r in artificial:
            self.T[-1,:] = self.T[-1,:] - self.T[r,:]


        if self.debug:
            print("subtracting artifical variables from objective row  ")
            print("T:")
            print(self.T)


    def pivot_column(self,T,tol=global_tolerance):        
        """Go through the objective row and find the minimum entry
           above tolerance"""
        
        # Ignore all positive values where: positive is defined as
        # anything greater than -tol
        ignored = lambda e: (e is None) or (e >= -tol)

        objective = T[-1,:-1].masked_where(ignored)
        idx =  objective.min_index()
        
        if self.debug:
            print("pivot-column:%s index : %s\n" % (T[-1,:-1], idx[1]))
            
        return idx


    def pivot_row(self,T,pivcol,tol,phase=1):
        """ Find the appropriate pivot row. """
        # Skip objective rows, in first phase we have two objective rows
        # including a psuedo-objective row.

        skip_rows = 2 if (phase == 1) else 1

        # Mask values less than tolerance
        ignored = lambda e: (e is None) or ( e <= -tol ) # not negative tolerance

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
        if self.debug:
            print("q: %s" % q)
        # this call is not trivial or obvious
        return q.min_index()


    def do_pivot(self, T, pivrow, pivcol, basis, phase):
        "Perform the pivot operation using pivot row and pivot column and basis"

        if self.debug:
            print("Before Pivot[%d,%d] " % (pivrow,pivcol))
            print("T : \n%s" % T)
            print("Basic Variables : %s\n" % basis)

        basis[pivrow] = pivcol
        T[pivrow,:]   = T[pivrow,:] / T[pivrow][pivcol]

        for irow in range(T.shape()[0]):
            if irow != pivrow:
                T[irow,:] = T[irow,:] - T[pivrow,:] * T[irow,pivcol]

        if self.debug:
            print("After Pivot[%d,%d] " % (pivrow,pivcol))
            print("T : \n%s" % T)
            print("Basic Variables : %s\n" % basis)
            print("Current Objective Value:\n%s\n" % -T[-1,-1])


    def simplex_solve(self, T, n, basis, phase =2, tol = global_tolerance):
        # Ignore original and new objectives.

        # When called in phase-2 we are presented with both the original and the pseudo
        # objective function. Thus we must ignore the last two rows.
        
        if self.debug:
            print("simplex_solve phase: %d" % phase)
            print("T:")
            print(T)
            print("n: %d" % n)
            print("basis: %s " % basis)

        if phase == 1:
            m = T.shape()[0] - 2

        elif phase == 2:
            m = T.shape()[0] - 1

        self.nvars = (T.shape()[1]-1)

        solution = [0] * self.nvars

        max_iterations = 5000
        nit = 0

        # Identify and substitute
        if phase == 2:

            # Identify aritificial variables still in the objective
            ncols = T.shape()[1]
            is_artificial = lambda idx : basis[idx] > (ncols - 2)

            # Check basis for artificial variables
            artificial_variables = list(filter(is_artificial, range(len(basis))))

            if self.debug:
                print("basic variables : %s \n" % self.basis)
                print("artificial pivot variables : %s " % artificial_variables)

            # This should pivot out all the artificial variables
            for pivrow in artificial_variables:
                def non_zero_col(col):
                    return self.T[pivrow,col] != 0
                
                pivcols = filter(non_zero_col,range(ncols -1))

                if len(pivcols) == 0:
                    continue
                
                self.do_pivot(T,pivrow,pivcols[0],basis,phase)

        complete = False

        while not complete and (nit < max_iterations):

            nit += 1
            pivcol_found, pivcol = self.pivot_column(T,tol)

            if self.debug:
                print("T:")
                print(T)

            if not pivcol_found:  # Finished with all the columns, in basic form
                status, complete = 0, True
                break

            pivrow_found, pivrow = self.pivot_row(T,pivcol,tol,phase=phase)

            if not pivrow_found: # not finding the pivot row is very serious.
                if self.debug:
                    print("cound not find pivot row")
                # Unbounded
                status, complete = 3, True
                break
            if self.debug:
                print("pivot: (pivrow,pivcol): (%d,%d)" %(pivrow,pivcol))

            if not complete: # perform the pivot on pivot entry
                self.do_pivot(T,pivrow,pivcol, basis,phase)

        if complete:
            if self.debug:
                print("Pivot Result == Phase: %d == " % phase)
                print("Basic Variables : %s\n" % basis)
                print("Status : %s" % status)
                print("Complete : %s" % complete)
                print("T : \n%s" % T)
        else :
            # unbounded
            return (3,True)

        return (status, complete)


    def solve(self):
        # Pivot to basic flexible.
        
        status, complete = self.simplex_solve(self.T,self.n,self.basis,phase=1)
        pseudo_objective  = self.T[-1,-1]
        if self.debug:
            print("pseudo-objecive : %d " % pseudo_objective)

        if abs(pseudo_objective) < global_tolerance:

            # remove pseudo-objective row
            self.T.del_rows([len(self.T)-1])

            # Remove artificial variable columns
            self.T.del_columns(range(self.n + self.n_slack,
                                     self.n + self.n_slack + self.n_artificial))
        else:
            # Infeasible without starting point
            status = 2
            if self.debug:
                print("Infeasible soltion")

        if status == 2: # infeasible solution
            return (-1,None)


        # Move to phase 2 compuation
        status, complete = self.simplex_solve(self.T,self.n,self.basis,phase = 2)

        if status == 0:                
            solution = list_wrapper([0] * (self.n + self.n_slack + self.n_artificial))
            
            if self.debug:
                print(" m: %d, n: %d"  %  (self.m, self.n))
                print(" basis: %s"     %   self.basis[:self.m])

                # assign the value of the constant
                print("constants : %s" % list(self.T[:self.m,-1]))
            solution[self.basis[:self.m]] = list(self.T[:self.m,-1])
            solution = solution[:self.n]
            if self.debug:
                print("solution: %s"  %  solution)

            # for i in range(len(basis)):
            # solution[basis[:m]] = T[:]

            obj,ansx = -self.T[-1,-1], solution
            return (0,ansx)
        
        elif status == 3:
            return (1,None)
        else:
            return (-1,None)

        print("status : %d found objective: %f solution: %s" % (status,obj,ansx))

        #return (status,ansx)
