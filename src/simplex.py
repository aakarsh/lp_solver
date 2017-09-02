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

class UnboundedError(Exception):
    "Raised when the solution for given equations is unbounded "
    pass

class InfeasibleError(Exception):
    "Raised when the solution for given equations is infeasible(has no-solution)."
    pass


class fp_ops:

    def __init__(self,rel_tol=global_tolerance, abs_tol=global_tolerance):
        self.abs_tol = abs_tol
        self.rel_tol = rel_tol
        self.tolerance = global_tolerance

    def ispositive(self,x):
        return x and self.isgt(x,Decimal(0.0))
    #(not self.iszero(x)) and self.sub(x , 0.0) > Decimal(global_tolerance)

    def iszero(self,a): return self.iseq(a,0.0)

    def iseq(self,a, b):
        return abs(self.sub(a,b)) <= max(self.mul(self.rel_tol, max(abs(a), abs(b))), self.abs_tol)
        #return abs(self.sub(a,b)) <= global_tolerance

    def isgt(self,a,b):
        return Decimal(a)-Decimal(b) > -Decimal(self.abs_tol)  #and a > b

    def islt(self,a,b):
        return not self.iseq(a,b) and a < b

    def round(self,a):
        return Decimal(0.0) if abs(Decimal(a) - Decimal(0.0)) <= Decimal(global_tolerance) else Decimal(a)

    def add(self, a, b):
        return self.round(a)+self.round(b)

    def sub(self, a, b):
        return self.round(a)-self.round(b)

    def mul(self, a, b):
        return self.round(a)*self.round(b)

    def div(self, a, b):
        return self.round(a)/self.round(b)

class SlackForm:
    "Represents simplex in slack form ..."

    class Constraint:
        "Helper class to deal with constraints..."

        def __init__(self,constant=None,coefficients=None,fp_op=None,**kwargs):
            self.fp_op= fp_op
            if "row" in kwargs:
                self.row = kwargs["row"]
            else:
                self.row = [constant]
                self.row.extend(coefficients)

        def __repr__(self):
            return "%s" % Matrix.PrettyPrinter.format_list(self.row)
        #(Matrix.PrettyPrinter.format_number(self.row[0]),Matrix.PrettyPrinter.format_list(self.row[1:]))

        def get_row(self):             return self.row
        def get_coefficients(self):    return self.row[1:] # O(n)
        def get_constant(self):        return self.row[0]
        def get_coefficient(self,idx): return self.row[idx+1]

        def sexp(self,other,factor):
            "Perform linear transform on two vectors self and other of a+(factor*b) "
            opt    = lambda v :  v if not  v is None else 0.0
            inc    = lambda a,c : self.fp_op.add(opt(a) ,self.fp_op.mul( factor,opt(c)))
            #
            retval =  [inc(a,c)  for a,c in zip(self.get_row(),other.get_row())]
            for idx in range(len(retval)):
                if(self.fp_op.iszero(retval[idx])):
                    retval[idx] = 0.0
            return retval

        def substitute(self,entering_idx,other):
            # if debug:
            #     print("subs: %-20s @ [%3d] => %-20s " % (other, entering_idx, self))

            variable_coefficient = self.get_coefficient(entering_idx)

            # TODO maybe not correct return the entering constraints
            if not variable_coefficient or self.fp_op.iszero(variable_coefficient):
                #if debug: print("subs => : %s " % (self))
                return self

            # treat constant and coefficients as single row
            return_row = self.sexp(other,variable_coefficient)
            return_row[entering_idx+1] = None
            retval = SlackForm.Constraint(row = return_row,fp_op=self.fp_op)

            #if debug: print("subs => : %s " % (retval))

            return retval

    def replace(self,**kwargs):
        if "A" in kwargs:
            self.A = kwargs["A"]
        if "b" in kwargs:
            self.b = kwargs["b"]
        if "c" in kwargs:
            self.c = kwargs["c"]
        if "v" in kwargs:
            self.v = kwargs["v"]
        if "n" in kwargs:
            self.n = kwargs["n"]
        if "m" in kwargs:
             self.m = kwargs["m"]

    def _pivot_col(self,T,tol=global_tolerance):
        """Go through the objective row and find the minimum entry above
           tolerance"""
        # ignore all positive values where: positive is defined as
        # anything greater than -tol
        ignored = lambda e: (e is None) or (e >= -tol)
        objective = T[-1,:-1].masked_where(ignored)
        return objective.min_index()


    def _pivot_row(self,T,pivcol,tol,phase=1):
        """ Find the appropriate pivot row. """

        first_phase = phase == 1
        skip_rows = 2 if first_phase else 1

        # Mask values less than tolerance
        ignored = lambda e: (e is None) or (e <= tol)

        ma = T[:-skip_rows,pivcol].masked_where(ignored)

        # All pivot column entries
        if ma.count_notnone() == 0 : return (False,None)

        mb = self.T[:-skip_rows,-1].masked_where(ignored)

        mr = mb / ma
        print("mr: %s" % mr)

        return mr.min_index()


    def _simplex_solve(self, T, n, basis, phase =2, tol = global_tolerance):
        # Ignore original and new objectives.

        if phase == 1:
            m = T.shape()[0] - 2
        elif phase == 2:
            m = T.shape()[0] - 1

        nvars = (T.shape()[1]-1)
        complete = False
        solution = [0] * nvars

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
            pivcol_found, pivcol = self._pivot_col(T,tol)

            if debug:
                print("T:")
                print(T)

            if not pivcol_found:  # Finished with all the columns, in basic form
                status, complete = 0, True
                break

            pivrow_found, pivrow = self._pivot_row(T,pivcol,tol)

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



    def do_pivot(self, T, pivrow, pivcol, basis, phase):
        "Perform the pivot operation using pivot row and pivot column and basis"

        if debug:
            print("Before Pivot[%d,%d] " % (pivrow,pivcol))
            print("T : \n%s" % T)
            print("Basic Variables : %s\n" % basis)


        basis[pivrow] = pivcol
        pivval = T[pivrow][pivcol]
        T[pivrow,:]  = T[pivrow,:] / pivval

        for irow in range(T.shape()[0]):
            if irow != pivrow:
                T[irow,:] = T[irow,:] - T[pivrow,:] * T[irow,pivcol]

        if debug:
            print("After Pivot[%d,%d] " % (pivrow,pivcol))
            print("T : \n%s" % T)
            print("Basic Variables : %s\n" % basis)


    def __init__(self,**kwargs):
        # expand
        if ("simplex" in kwargs) and kwargs["simplex"]:
            simplex = kwargs["simplex"]
            n,m = (simplex.n,simplex.m)
            A = simplex.A
            b = list(map(Decimal,simplex.b))
            c = list(map(Decimal,simplex.c))
            v = simplex.v
            self.fp_op = simplex.fp_op
            # TODO: Add tableau here which parrallels the simplex tableau
        else:
            # TODO: Make the tableau here.
            n,m = (kwargs["n"],kwargs["m"])
            A = kwargs["A"]
            b = kwargs["b"]
            c = list(map(Decimal,(kwargs["c"])))
            v = kwargs["v"]
            self.fp_op = fp_ops(rel_tol=global_tolerance,abs_tol=global_tolerance)

        # A represents the table of size total x total
        # with None possible entreis
        nvars = n+m
        self.A = Matrix(nvars,nvars)

        # We will use a numbering where
        # dependent variables will be from [x_n to x_nvars)
        # independent variables will be from [x_0 to x_n)
        # x_nvars will be a psuedo variable used in case the initial
        # simplex is not in basic form.

        for i in range(n):
            for j in range(m):
                self.A[i][j+n] = -A[i][j]

        if debug:
            print("\nA:\n%s" % self.A)
            print("Input Equations:")
            print(Matrix(matrix=A))

        if debug:
            print("orig-b:\n%s" % Matrix.PrettyPrinter.format_list(b))

        self.b = list(map(float,b))

        # Fill rest of it with None
        self.b.extend([None]*(nvars-len(b)))

        if debug:
            print("b:\n%s"% self.b)

        self.c = [None] * n
        self.c.extend(map(float,c))

        if debug:
            print("orig-c:%s" % Matrix.PrettyPrinter.format_list(c))
            print("c:%s"      % Matrix.PrettyPrinter.format_list(self.c))

        self.v = v
        # Original variables
        self.independent   = set(range(n,nvars))

        # Slack varaibles
        self.dependent     = set(range(0,n))

        if debug:
            print("dependent : %s, independent : %s" % (self.independent,self.dependent))

        # TODO - do explicit comparisons
        leaving_coeff,leaving_idx = list_wrapper(self.b).min_index()

        if debug:
            print("leaving_coeff:%d leaving_idx: %d" % (leaving_coeff,leaving_idx))

    def to_constraint(self,equation_idx = None,type='constraint'):
        "Particular equation-id from simplex wrapped  equation"
        if  equation_idx is None: # return objective constant
            return SlackForm.Constraint(self.v,self.c,fp_op=self.fp_op)
        return SlackForm.Constraint(self.b[equation_idx],self.A[equation_idx],fp_op=self.fp_op)

    def store(self,idx,new_eq):
        "Store constraint back into slack form"
        self.b[idx],self.A[idx] = new_eq.get_constant(), new_eq.get_coefficients()

    def substitute_objective(self,entering_constraint,entering_idx):
        new_obj = self.to_constraint(type='objective').substitute(entering_idx,entering_constraint)
        #if debug: print("subs-obj: (%s,%s) => %s " % (self.c,self.v,new_obj))
        self.v,self.c = new_obj.get_constant(),new_obj.get_coefficients()

    def substitute_equation(self,equation_idx,entering_constraint,entering_idx):
        "Substititue enterint constraint into equation "
        new_eq  = self.to_constraint(equation_idx).substitute(entering_idx,entering_constraint)
        self.store(equation_idx,new_eq)

    def make_pivot_constraint(self,leaving_idx, entering_idx):
        "Construct constriant based on leaving and entering ids "
        pivot = self.A[leaving_idx][entering_idx]
        pivot_inverse = self.fp_op.mul(-1 ,self.fp_op.div (1.0 , pivot))
        if self.fp_op.iszero(pivot_inverse):
            pivot_inverse = 0.0
        # Treat slack row as [b,A] @ leaving_idx
        row = chain([self.b[leaving_idx]],self.A[leaving_idx])

        def pivot_scale(value,default=0.0):
            value =  self.fp_op.mul(pivot_inverse , value) if (not value is None)  else default
            if self.fp_op.iszero(value): value = 0.0
            return value

        scaled_row = list(map(pivot_scale ,row))

        # leaving idx will get 1/pivot
        scaled_row[entering_idx+1],scaled_row[leaving_idx +1]  = None, self.fp_op.mul((-1) , (pivot_inverse))
        if self.fp_op.iszero(scaled_row[leaving_idx+1] ): scaled_row[leaving_idx+1]  = 0.0
        return SlackForm.Constraint(row=scaled_row,fp_op=self.fp_op)


    def pivot(self,entering_idx,leaving_idx,phase=0):
        "Pivot on the entering and leaving ids"
        #
        entering_constraint = self.make_pivot_constraint(leaving_idx,entering_idx)
        self.store(entering_idx,entering_constraint)
        #
        if debug:
            print("(Eq:%d) => :%s" % (entering_idx,entering_constraint))
            print("\nA:\n%s" % Matrix(matrix=self.A))
        #
        for idx in self.dependent:
            if idx == leaving_idx: continue
            self.substitute_equation(idx,entering_constraint,entering_idx)

        self.substitute_objective(entering_constraint,entering_idx)
        #
        # move leaving_idx  from dependent   to independent set
        # move entering_idx from independent to dependent set
        if debug:
            print("dependent: %s   -  %d + %d " % (self.dependent,leaving_idx,entering_idx))

        def set_replace(s,rem=None,add=None):
            s.remove(rem)
            s.add(add)

        set_replace(self.dependent,add=entering_idx,rem=leaving_idx)

        if debug:
            print("dependent=> %s" % self.dependent)

        if debug:
            print("independent: %s -  %d + %d" % (self.independent,entering_idx,leaving_idx))

        set_replace(self.independent,add=leaving_idx,rem=entering_idx)


        if debug:
            print("independent: %s" % (self.independent))

    def objective_has_positive_coefficients(self):
        return list_wrapper(self.c).contains(predicate = lambda x: x and x >= -global_tolerance)
        #return list_wrapper(self.c).contains(predicate = self.fp_op.ispositive)

    # TODO: This is a deviation since since we use min for simplex
    # but first negative index here.
    def pick_entering_idx(self):
        idx = 0
        min_idx,min_value = -1,None
        for c_v in self.c:
            if c_v is not None and c_v >= global_tolerance and ((min_value is None) or c_v < min_value ):
                min_value = c_v
                min_idx   = idx
            idx+=1

        return min_idx if ((not min_value is None ) and min_value >= global_tolerance) else None

    #list_wrapper(self.c).find_first_idx(in_set = self.independent,
    #by     = lambda x: x and x >= -global_tolerance)
    #self.fp_op.ispositive)

    def pick_leaving_idx(self,entering_idx):
        if debug:
            print("pick_leaving_idx: dependent:%s entering:%d" % (self.dependent,entering_idx))

        slack = [None] * len(self.b)

        # TODO: Seeing non-None  values in b where it is not a depenedent

        if debug:
            print("b: %s" %  Matrix.PrettyPrinter.format_list(self.b))
            print("c: %s" %  Matrix.PrettyPrinter.format_list(self.c))
            print("\nA:\n %s\n" % Matrix(matrix=self.A))

        for idx in self.dependent:
            constant = self.b[idx]
            coeff    = self.fp_op.mul(-1 , self.A[idx][entering_idx])

            # if debug:
            #     print("slack[%d]=:: [%s]/[%s] " % (idx,constant,coeff))

            if self.fp_op.ispositive(coeff):
                slack[idx]  = self.fp_op.div(constant,coeff)
                if self.fp_op.iszero(slack[idx]): slack[idx] = 0.0

            if constant <= global_tolerance or coeff <= global_tolerance:
                slack[idx] = 0.0

                # if debug:
                #     print("slack[%d]=> [%f]/[%f] " % (idx,constant,coeff))
                #     if slack[idx]:
                #         print("%2.2f",slack[idx])
                #     else: print()

        # if debug: print("A col[%d] %s" % (entering_idx,[self.A[idx][entering_idx] for idx in self.dependent]))
        # if debug: print("b %s" % [self.b[idx] for idx in self.dependent])

        #TODO : Perform explicit Comparisons Instead
        min_slack,slack_idx = list_wrapper(slack).min_index()

        if debug:
            if min_slack:
                print("slack: %s \nmin_slack[%2.5f] :%2.5f " % (Matrix.PrettyPrinter.format_list(slack),slack_idx,min_slack))
            else:
                print("min_slack is None slacks: %s",   Matrix.PrettyPrinter.format_list(slack ))
            if slack_idx:
                print("leave: %d" % slack_idx)

        return slack_idx if min_slack else None


    def __repr__(self):
        return "A:\n%s\nb:\n%s" % (Matrix(matrix=self.A), Matrix.PrettyPrinter.format_list(self.b))

    def solve(self,phase):
        cnt = 0
        opt_d = lambda v: Decimal(v) if not v is None else None
        self.b = list(map(opt_d,self.b))

        while self.objective_has_positive_coefficients():
            #
            # TODO: Print out the value of objective function at each
            # we are here and need to solve this.
            #
            if debug: print(self)
            cnt += 1

            if debug:
                print("=========== Pivot Iteration : %d Phase %d ==========" % (cnt,phase) )
                print("%s" % Matrix(matrix=self.A))
                # why is this like this
                print("objective:\n%s\n" % Matrix.PrettyPrinter.format_list(self.c))
                print("constant:\n%s " % self.v)


            entering_idx = self.pick_entering_idx()
            leaving_idx  = self.pick_leaving_idx(entering_idx)

            if (leaving_idx is None)  or  (entering_idx is None):  # unbounded
                raise UnboundedError()

            if debug:
                print("pivot:  %s >> (dep) >> %s" % (entering_idx,leaving_idx))

            self.pivot(entering_idx , leaving_idx,phase)

        if debug:
            print("\nA:\n%s\n"   % Matrix(matrix=self.A))
            print("\nb:\n%s\n"   % Matrix.PrettyPrinter.format_list(self.b))
            #print("objective after pivots: %s" % Matrix.PrettyPrinter.format_list(self.c))

        # Variable assignment after pivot stops.
        assignment = [0] * len(self.b)
        #
        # All dependent variables get value of constant
        for var_idx in self.dependent:
            assignment[var_idx] = self.b[var_idx]
        #
        opt_value = self.v
        for i in range(0,len(self.c)):
            coefficient = self.c[i]
            value = assignment[i]
            if not coefficient is  None and not value is None  and not self.fp_op.iszero(coefficient) and not self.fp_op.iszero(value):
                opt_value = self.fp_op.add(opt_value,self.fp_op.mul(coefficient, value))

        if debug:
            print("opt-value:%d ,assignment:%s" % (opt_value,Matrix.PrettyPrinter.format_list(assignment)))
        # does that mean all assigned variables were explicity assigned.
        # how does it correspond to initial non-slack assignmentx
        return (opt_value, assignment)


class Simplex:

    def __init__(self,A,b,c,n,m,debug=False): # called here so modify in solve
        self.A = A
        self.b = b
        self.c = c
        # start with (no-constant)
        # in objective function
        self.v = 0
        self.n = n
        self.m = m
        self.debug  = debug
        self.fp_op = fp_ops(rel_tol=global_tolerance,abs_tol=global_tolerance)





    def find_basic_feasible(self, min_idx):
        "Find basic feasible solution."

        slackform = SlackForm(A=self.A,b=self.b,c=self.c,v=0,n=self.n,m=self.m)
        aux_simplex = self.make_auxiliary_form()
        aux_sf = SlackForm(simplex = aux_simplex)
        nvars = self.n + self.m

        if self.debug:
            print("%d => (dependent) => %d " % (nvars,min_idx))
        aux_sf.phase = 0
        aux_sf.pivot( self.n + self.m , min_idx,phase=0)
        if self.debug: print("========== Search for Auxiliar Form ========== ")
        (opt,ansx) = aux_sf.solve(phase=0)
        if self.debug: print("========== End Search for Auxiliar Form ========== ")

        if not self.fp_op.iszero(opt):
            if self.debug: print("Raising Infeasible %f %s" % (opt,self.fp_op.iszero(opt)))
            raise InfeasibleError()

        if self.debug:
            print("dependent : %s" % aux_sf.dependent)

        new_obj = SlackForm.Constraint(slackform.v, slackform.c, fp_op = self.fp_op)

        for idx in aux_sf.dependent:
            constant, equation = aux_sf.b[idx],aux_sf.A[idx]
            entering_constraint = SlackForm.Constraint(constant,equation,fp_op = self.fp_op)
            new_obj = SlackForm.Constraint.substitute(new_obj,idx,entering_constraint)

        # if self.debug: print("sub-obj: %s" % new_obj)
        if self.debug: print("END:subsitite-auxiliary-form-objective ")

        list_wrapper(new_obj.get_coefficients()).set_values( None, aux_sf.dependent)
        list_wrapper(aux_sf.b).set_values( None, aux_sf.independent )

        # set all unused rows  and columns to None
        aux_sf.A.set_value(None,rows= aux_sf.independent,columns= aux_sf.dependent)

        # using new objective and return new slack form need to remove
        # x_nvars
        if self.debug:
            print("Objective: %s" % new_obj)
            print("A: \n%s" % aux_sf.A)
            print("b: \n%s" % Matrix.PrettyPrinter.format_list(aux_sf.b))
            print("depenedents:   \n%s"  % aux_sf.dependent)
            print("indepenedents: \n%s"  % aux_sf.independent)

        # remove last entry from rows
        aux_sf.A.pop_column()
        aux_sf.A.pop_row()

        # remove last entry from constant
        aux_sf.b.pop()

        # replace the new objective
        aux_sf.replace(c = new_obj.get_coefficients(),
                       v = new_obj.get_constant())

        # remove one independent variable
        aux_sf.independent.remove(nvars)
        aux_sf.replace(m=len(aux_sf.independent)-1)

        if self.debug:
            print("objective: %s " % new_obj)
            print("A: \n%s" % aux_sf.A)
            print("b: \n%s" % Matrix.PrettyPrinter.format_list(aux_sf.b))
            print("depenedents:   \n%s"  % aux_sf.dependent)
            print("indepenedents: \n%s"  % aux_sf.independent)

        return aux_sf

    def solve(self): # phase used to track
        if self.debug: print("\nSimplex solve");
        try:
            if self.debug:
                print("simplex:b %s"%self.b)

            min_constant, idx = list_wrapper(self.b).min_index();

            #TODO : Error-Here
            if  self.fp_op.islt(min_constant , 0): # 0-comparison
                if self.debug:
                    print(" Currently not in basic form :%d - index %d " % (min_constant,idx))

                simplex_basic_form = self.find_basic_feasible(idx)

                try:
                    if self.debug: print("-" * 100)
                    if self.debug: print("========== Solve Transformed  basic form ========== ")
                    opt,ansx = simplex_basic_form.solve( phase = 1)
                    if self.debug: print("========== End Transformed  basic form ========== ")
                    self.optimum = opt
                    self.assignment = ansx
                    ## slack variables are the first n varialbes
                    del self.assignment[:self.n]
                    self.anst = 0
                    return (0,self.assignment)
                except UnboundedError as err:
                    self.anst = 1
                    return(1,None)
                except InfeasibleError as err:
                    self.anst = -1
                    return (-1, None)

            # basic solution is already feasible
            opt,ansx = None,None
            if self.fp_op.ispositive(min_constant): # basic form already in feasible form
                opt,ansx = SlackForm(simplex=self).solve(phase=2)
                self.optimum = opt
                self.assignment = ansx
                del self.assignment[:self.n]

            if self.debug:
                print("Optimimum Value: %f" % opt)
                print("Assignment: %s "     % Matrix.PrettyPrinter.format_list(self.assignment))

            return (0,self.assignment)
        except UnboundedError as err:
            return(1,None)
        except InfeasibleError as err:
            return (-1, None)

    def make_auxiliary_form(self):
        "Make auxiliiary form"
        nvars = (self.n + self.m)
        # new size
        aux_size = (nvars + 1)
        # create copy of b and extend it for x_{nvars}
        #
        aux_b = self.b[:]
        if self.debug:
            print("make_auxiliary_form: b:\n%s\naux-b:\n%s" % (self.b,aux_b))
        #
        # set objective function: -x_{nvars}
        aux_c = ([0] * self.m)
        aux_c.append(-1)
        aux_v = 0
        #
        aux_dep = set(range(self.n))
        aux_indep = set(range(self.n,nvars))
        #
        # simplex will have one extra column and row
        # to represent x_{nvars}
        aux_A = Matrix(self.n,self.m+1)

        # copy An
        aux_A.copy(self.A)

        # set nvars column to -1
        aux_A.set_column_value(self.m,-1)

        return Simplex(aux_A,aux_b,aux_c,self.n,self.m+1)


    def verify_bounds(self,tolerance=1e-3):
        r = 0
        for row in self.A:
            dot = [ row[i] * ansx[i] for i in range(len(row))]
            sum = 0
            for x in dot:
                sum += x
            assert(sum <= (self.b[r]) + (1e-3))
            r+=1

    def solve_scipy(self,tolerance=global_tolerance):
        import sys
        import numpy as np
        import lprog

        linprog_res = lprog.linprog([-x for x in simplex.c],
                                    A_ub = simplex.A,
                                    b_ub = simplex.b,
                                    options  = { 'tol': tolerance, },
                                    callback = lprog.linprog_verbose_callback)
        return linprog_res


    def verify_scipy(self,tolerance=1e-3):
        import sys
        import numpy as np
        from scipy.optimize import linprog

        linprog_res = self.solve_scipy()

        if self.debug:
            print("max_myprog = %f " % self.optimum)
            print ("linprog_res.status : %f " % linprog_res.status)
        if self.anst == 0:
            assert (linprog_res.status == 0)

        max_ref = np.dot(simplex.c, linprog_res.x)

        if self.debug:
            print('max_ref = %f' % max_ref)
            print('tolerance = %f' % tolerance)
            print("delta= %f " % abs(self.fp_op.sub(simplex.optimum , max_ref)))
            print(" %s " % linprog_res)

        assert(abs(self.fp_op.sub(simplex.optimum , max_ref)) <= tolerance)

    def __repr__(self):
        return "A:\n%s\nb:\n%s\nc:\n%s\n" % (self.A,self.b,self.c)

    def  __str__(self):
        return self.__repr__()

