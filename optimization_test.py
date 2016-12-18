# -*- coding: utf-8 -*-
"""
Created on Fri Oct  9 15:57:18 2015

@author: hazem.soliman
"""

from scipy.optimize import minimize

def sum_points(x):
    N=3
    sum_p = []
    for i in range(N):
        for j in range(i+1,N,1):
            print('i=',i)
            print('j=',j)
            sum_p.append( -(x[i] - x[j])**2)
    print('end')
    return min(sum_p)

N=3
si = [2,2,3]

fun = lambda x: sum([sum([-(x[i] - x[j])**2 for j in range(i+1,N,1)]) for i in range(N)])

fun2 = lambda x: -x[N]


cons3 = ({'type': 'ineq', 'fun': lambda x:  x[0] - 2},
        {'type': 'ineq', 'fun': lambda x:  -x[0] + 8},
        {'type': 'ineq', 'fun': lambda x:  x[1] - 3},
        {'type': 'ineq', 'fun': lambda x:  -x[1] + 7})
        
#cons1 = tuple({'type': 'ineq', 'fun': lambda x:  x[i] - si[i]/2} for i in range(N))
#cons2 = tuple({'type': 'ineq', 'fun': lambda x:  -x[i] + (10-si[i]/2)} for i in range(N))
#cons = cons1 + cons2

cons = []
for i in range(N):
    print(i)
    cons.append({'type': 'ineq', 'fun': lambda x, ind = i:  -x[ind] + (10-si[ind]/2)})
    cons.append({'type': 'ineq', 'fun': lambda x, ind = i:  x[ind] - si[ind]/2})
    
for i in range(N):
    for j in range(i+1,N,1):
        cons.append({'type': 'ineq', 'fun': lambda x, ind1 = i, ind2 = j:  (x[ind1] - x[ind2])**2 -x[N]})

cons = tuple(cons)



bnds = [(0, 10)]*(N) + [(-100,100)]

res = minimize(fun2, (2, 0, 3, 5), method='SLSQP', bounds=bnds, constraints=cons)
print(res)