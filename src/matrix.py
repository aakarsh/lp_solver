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

def map_optional(func,ls):
    "Maps of the list skipping over but preserving None objects "
    if ls is None: return None    
    result =[]
    for l in ls:
        if l:
            result.append(func(l))
        else:
            result.append(None)
    return result

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
