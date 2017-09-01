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

debug = False
global_tolerance = 1e-12
decimal.getcontext().prec = 64

class Matrix:

    "Simplistic matrix class."

    def __init__(self,nrows = None, ncols = None,**kwargs):
        if "matrix" in kwargs and kwargs["matrix"]:
            self.matrix = kwargs["matrix"]
        else:
            initial = None
            if "initial" in kwargs :
                initial = kwargs["initial"]
            self.matrix =  [ [initial for x in range(ncols)] for y in range(nrows) ]


    def _slice(self,i,type="row"):
        if type == "row": max = len(self.matrix)
        else:             max = len(self.matrix[0])

        if isinstance(i,slice):
            start ,stop,step = i.start,i.stop,i.step
            if start is None: start = 0
            if start < 0: start = max + start
            if stop  is None: stop  = max
            if stop < 0 : stop = max + stop
            if step  is None: step  = 1
            return slice(start,stop,step)
        else:
            if i < 0 : i = max + i
            return slice(i,i+1,1)

    def _map(self,row_slice,col_slice,action = lambda rs,cs,matrix: None ):
        r = 0
        count = 0
        loop_ctx={}
        for rs in range(row_slice.start,row_slice.stop,row_slice.step):
            c = 0
            for cs in range(col_slice.start,col_slice.stop,col_slice.step):
                loop_ctx["r"],loop_ctx["c"] = r,c
                loop_ctx["count"] = count
                loop_ctx["element"] = self.matrix[rs][cs]
                e = self.matrix[rs][cs]
                action(rs,cs,loop_ctx)
                c += 1
                count +=1
            r += 1


    def __getitem__(self,idx):
        if hasattr(idx,'__getitem__'):
            i,j = idx
            is_cell_ref = not isinstance(i,slice) and not isinstance(j,slice)

            # print("Matrix::__getitem__(i,j) : (%s,%s) "  % idx)
            row_slice = self._slice(i, type="row")
            # print("Matrix::__getitem__(row_slice: %s)  " % row_slice)
            col_slice = self._slice(j, type="column")
            # print("Matrix::__getitem__(col_slice: %s)  " % col_slice)

            retval = []

            def save_cell(rs,cs,ctx):
                v  = self.matrix[rs][cs]
                retval.append(v)

            self._map(row_slice,col_slice, save_cell)

            # TODO: Why do I have this here ?
            if is_cell_ref and len(retval) == 1 :
                return retval[0]

            return list_wrapper(retval)
        else:
            return self.matrix[idx]


    def __setitem__(self,idx,value):
        # tuplie
        if isinstance(idx,tuple):
            sl = idx
            row_slice = self._slice(sl[0],type="row")
            col_slice = self._slice(sl[1],type="col")

            value_iterator = iter(value) if hasattr(value,"__iter__") else None

            def do_assign(rs,cs,ctx):
                r,c,count = ctx["r"],ctx["c"],ctx["count"]

                if hasattr(value,"__getitem__"):
                    if r < len(value) and  hasattr(value[r],"__getitem__"):
                        if c < len(value[r]):
                            self.matrix[rs][cs] = value[r][c]
                    elif count < len(value) :
                        self.matrix[rs][cs] = value[count]
                elif not value_iterator is None :
                    self.matrix[rs][cs] = next(value_iterator)
                else:
                    self.matrix[rs][cs] = value

            self._map(row_slice,col_slice,do_assign)

        # assume it is just a list of numbers?
        # elif (isinstance(idx,list) or isinstance(idx,list_wrapper)) and hasattr(value,"__iter__"):            
        #     value_iterator = iter(value)            
        #     for i in idx:
        #         self.matrix[i] = next(value_iterator)
                
        else:
            self.matrix[idx] = value

    def fill_diagonal(self,sl,value):
        row_slice = self._slice(sl[0],type="row")
        col_slice = self._slice(sl[1],type="column")
        step_size = 0
        for rs in range(row_slice.start,row_slice.stop,row_slice.step):
            cs = step_size + col_slice.start
            if cs > col_slice.stop:
                break;
            self.matrix[rs][cs] = value
            step_size += 1

    def __delitem__(self,idx):  del self.matrix[idx]
    def __iter__(self):         return iter(self.matrix)
    def __len__(self):          return len(self.matrix)

    def pop(self):              self.matrix.pop()
    def pop_row(self):          self.pop()

    def pop_column(matrix):
        for row in matrix: row.pop()

    def del_columns(self,col_idxs):
        "Delete_columns: (matrix, col, idxs)"
        for row in self:
            deleted = 0
            for idx in col_idxs:
                del row[idx - deleted]
                deleted += 1

    def del_rows(self, row_idxs):
        for idx in row_idxs:
            del self[idx]

    def __repr__(self):
        return Matrix.PrettyPrinter.format_table(self.matrix)

    def shape(self):
        return (len(self.matrix),len(self.matrix[0]))

    @staticmethod
    def zeros(shape):
        return Matrix(shape[0],shape[1],initial=0)

    class PrettyPrinter:

        @staticmethod
        def format_number(elem,width=15,precision=6):
            nfmt = "%%+%d.%df" %(width,precision)
            sfmt = "%%%ds" %(width)
            output=""
            if not elem is None:
                if  isinstance(elem,float) or isinstance(elem,Decimal):

                    output += nfmt % float(elem)
                else:
                    output += sfmt % elem
            else:
                output +=sfmt % "None"
            return output

        @staticmethod
        def format_list(ls,sep="|",end="\n"):
            output=""
            first = True
            for elem in ls:
                if first:
                    output += sep + " "
                    first=False
                output += Matrix.PrettyPrinter.format_number(elem,width = 14)
                output +=" %s " % sep
            output += end
            return output

        @staticmethod
        def format_table(matrix):
            "Format matrix as a table"
            output =""
            for row in matrix:
                output += Matrix.PrettyPrinter.format_list(row,sep="|",end="\n")
            return output

    def set_value_all_rows(matrix,row_idxs, value):
        for row in row_idxs:
            matrix[row] = [value] * len(matrix[row])

    def set_column_value(matrix,col,value):
        for row in range(len(matrix)):
            matrix[row][col]= value

    def set_value_all_columns(matrix,col_idxs, value):
        for row in range(len(matrix)):
            for col in col_idxs:
                matrix[row][col] = value

    def set_value(self,value,**kwargs):
        if "columns" in kwargs:
            self.set_value_all_columns(kwargs["columns"],value)
        if "rows" in kwargs:
            self.set_value_all_rows(kwargs["rows"],value)

    def copy(m1,m2,column_offset = 0,row_offset = 0):
        for r in range(len(m2)):
            for c in range(len(m2[r])):
                m1[r+row_offset][c+column_offset] = m2[r][c]

class list_wrapper():
    "Wrapper class around built in list providing some method"

    def __init__(self,*args):
        if len(args) == 1 and hasattr(args[0],'__iter__'):
            self.ls = list(args[0])
        else:
            self.ls = args

    def __repr__(self):         return Matrix.PrettyPrinter.format_list(self.ls,end="\n")
    def __getitem__(self,idx):  return self.ls[idx]

    def __setitem__(self,idx,value):
        if isinstance(idx,list):
            i = 0
            for id in idx:
                self.ls[id] = value[i]
                i += 1
        else:
            self.ls[idx] = value

    def __delitem__(self,idx):  del    self.ls[idx]

    def masked_where(self,condition):
        return list_wrapper([elem if not condition(elem) else None for elem in self.ls ])

    def _isnum(self,v): return isinstance(v,Decimal) or isinstance(v,float) or isinstance(v,int)
    def _iseq(self,v) : return isinstance(v,list)    or isinstance(v,list_wrapper)

    def min_index(self):

        idx = 0
        min = None
        min_idx = None
        found = False

        for elem in self.ls:

            if elem is None:
                idx += 1
                continue

            if min is None:
                min = elem
                min_idx = idx

            if min > elem:
                found=True
                min = pivot_elem
                min_idx  = idx
            idx += 1
        return (found,min_idx)

    def apply_op(self,value,operator):
        op = self.wrap_optinal(operator)
        if self._isnum(value):
            return list_wrapper([op(element,value) for element in self.ls])
        elif self._iseq(value):
            ret = [None] * len(self.ls)
            idx = 0
            for v in value:
                if not self.ls[idx] is  None and not v is None:
                    ret[idx] = op(self.ls[idx], v)
                idx += 1
            return list_wrapper(ret)

    def wrap_optinal(self,operator):
        return lambda a,b: None if (a is None) or (b is None) else operator(Decimal(a),Decimal(b))

    def __mul__(self,value):      return self.apply_op(value, lambda a,b: a * b )
    def __truediv__(self, value): return self.apply_op(value, lambda a,b: a / b if not (b is None) and not (abs(b) < global_tolerance) else None )
    def __add__(self,value):      return self.apply_op(value, lambda a,b: a+b)
    def __sub__(self,value):      return self.apply_op(value, lambda a,b: a-b)
    def __neg__(self):            return self.apply_op(-1,    lambda a,b: a*b)

    def __iter__(self):         return iter(self.ls)
    def __len__(self):          return len(self.ls)

    def __lt__(self,other):
        if isinstance(other,int):
            return self.bools(lambda e: e < other)

    def __gt__(self,other):
        if isinstance(other,int):
            return self.bools(lambda e: e > other)

    def bools(self,predicate = lambda v : 1 if not v is None else 0):
        "Retrun a 1/0 list with predicate applied to each element in the list"
        retval = [0] * len(self.ls)
        idx = 0
        for e in self.ls:
            if predicate(e):
                retval[idx] = 1
                idx += 1
        return list_wrapper(retval)

    def count(self,predicate):
        count = 0
        for e in self.ls:
            if predicate(e):
                count += 1
        return count

    def count_nonzero(self):
        def non_zero(n): return n!= 0
        return self.count(non_zero)

    def count_notnone(self):
        def not_none(n): return not n is None
        return self.count(not_none)

    def zip_set(l,member_set,by = lambda x: True):
        return [(i,l[i]) for i in member_set if by(l[i])]

    # def find_first_idx(l, in_set = None, by = lambda x: True):
    #     "Find first element of list set in by predicate "
    #     if not in_set:
    #         in_set = set(range(len(l)))
    #     for (idx, value) in l.zip_set(in_set,by):
    #         return idx
    #     return None

    def contains(ls , predicate = lambda x: x):
        return any(predicate(v) for v in ls)

    def set_values(l, value,idxs ):
        for idx in idxs: l[idx] = value

    def min_index(l,max=float('inf'),custom_min=min):
        "Return a pair of (value,index) of elment with minimum index"
        idx = 0
        min_idx,min_value = None,None
        for c_v in l:
            if c_v is not None and abs(c_v) >= global_tolerance and  ((min_value is None) or (min_value-c_v  >= global_tolerance)):
                min_value = c_v
                min_idx   = idx
            idx+=1
        return (min_value,min_idx) if ((not min_value is None ) ) else (None,None)

        # zl = zip(l,range(len(l)))
        # # hiding a comparison
        # return custom_min(zl, key = lambda z: z[0] if z[0] else max)


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

# from matrix import *

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
        """Go through the objective row and find the minimum entry above
           tolerance"""
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
        if self.debug:
            print("q: %s" % q)
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

        max_iterations = 500
        nit = 0

        # Identify and substitute
        if phase == 2:

            # Identify aritificial variables still in the objective
            ncols = T.shape()[1]
            is_artificial = lambda idx : basis[idx] > ncols - 2

            # Check basis for artificial variables
            variables = list(filter(is_artificial,range(len(basis))))

            if self.debug:
                print("basic variables : %s\n" % basis)
                print("artificial pivot variables : %s " % variables)

            for pivrow in variables:
                non_zero_col = lambda col: self.T[pivrow,col] != 0
                pivcols = filter(non_zero_col,range(ncols -1))

                if len(pivcols) == 0: continue
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
            # raise InfeasibleError()

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

# from matrix import *
# from simplex import *
# from scipy_lprog import *
# from tableau import *

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
    if debug:
        print("-------------------------------------------------")

    if debug:
        print("-------------------- Scipy --------------------")
    # scipy   = SciPy(A,b,c,n,m,debug)
    # l_anst, l_ansx = scipy.solve()
    if debug:
        print("---------------------End:Scipy--------------------")

    if debug: print("-------------------- Simplex --------------------")
    # simplex = Simplex(A,b,c,n,m,debug)
    # anst, ansx = simplex.solve()
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

