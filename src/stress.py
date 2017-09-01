#!/usr/bin/env python3.4
#
# Author :  https://gist.github.com/dtrebbien/396973cccf906d3b9afe949bad8c0613
#
#
from __future__ import print_function
import sys
import numpy as np
from scipy.optimize import linprog
from scipy import optimize
import subprocess

def print_lp(n,m,A,b,c,file=sys.stdout):
  print(n, m,file=file)
  for i in range(n):
    print(*(A[i][j] for j in range(m)), file=file)
  print(*b,file=file)
  print(*c,file=file)
  print(file=file)
  
log = open("stress.log","a")

while True:

  n, m = np.random.randint(1, 8, size = (2,))
  A    = np.random.randint(-100, 100, size = (n, m,))
  b    = np.random.randint(-1000000, 1000000, size = (n,))
  c    = np.random.randint(-100, 100, size = (m,))

  print_lp(n,m,A,b,c,file = sys.stdout)

  proc = subprocess.Popen(['./simplex.py'],
                          stdin = subprocess.PIPE,
                          stdout = subprocess.PIPE,
                          universal_newlines = True)
  
  print_lp(n,m,A,b,c,file = proc.stdin)

  stdoutdata, _ = proc.communicate()
  assert proc.returncode == 0

  stdout_lines = stdoutdata.splitlines()

  print(stdoutdata, end = '')

  # Try to reduce "false positive" results, where linprog() returns an answer that is
  # not correct. Use the idea for "a slightly modified procedure" described at:
  # https://www.coursera.org/learn/advanced-algorithms-and-complexity/discussions/all/threads/XBz2qmB5EeaqYRKO7-Ax0Q/replies/GdzNnHYvEealrxI7520HyQ/comments/AAAGGXbSEeaF8w5uWHT1BQ
  
  linprog_res = linprog(-c, A_ub = A, b_ub = b, options = { 'tol': 1e-4 })
  
  if linprog_res.status != 2:
    prev_linprog_res = linprog_res
    tolerances = [1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13]
    for tol in tolerances:
      linprog_res = linprog(-c, A_ub = A, b_ub = b, options = { 'tol': tol })
      # (no-soultion)
      if linprog_res.status == 2: # no solution
        linprog_res = prev_linprog_res
        break
      # ( unbounded |infinite ) previous liner program
      prev_linprog_res = linprog_res

  try:
    if linprog_res.status == 3:
      assert stdout_lines[0] == 'Infinity'
    elif linprog_res.status == 2:
      assert stdout_lines[0] == 'No solution'
    elif linprog_res.status == 0:
      x_ref = linprog_res.x
      print('x_ref =', ' '.join(list(map(lambda x: '%.18f' % float(x), x_ref))))
      assert stdout_lines[0] == 'Bounded solution'
      x_myprog = np.array([float(num_str) for num_str in stdout_lines[1].split(' ')])
      assert len(x_myprog) == len(x_ref)

      # Verify that all inequalities are satisfied.
      assert (np.dot(A, x_myprog) <= b + 1e-3).all()
    
      max_myprog = np.dot(c, x_myprog)
      print('max_myprog =', max_myprog)
      max_ref = np.dot(c, x_ref)
      print('max_ref =', max_ref)

      # Total pleasure differs from the optimum by at most 1e-3.
      assert (abs(max_myprog - max_ref) <= 1e-3)
      print()
  except AssertionError as err:
    print("%s" % err,file=log)
    print_lp(n,m,A,b,c,file = log)
    

  
