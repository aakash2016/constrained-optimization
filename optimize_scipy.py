import numpy as np
from scipy.optimize import minimize

def objective(x):
	return(pow((pow(x[0],2) + x[1] - 11),2) + pow((pow(x[1],2) + x[0] - 7),2))

def constraint(x):
    return(-26.0 + pow((x[0]-5.0), 2) + pow(x[1],2))

x0 = np.array([0.11, 0.1])

# show initial objective
print('Initial SSE Objective: ' + str(objective(x0)))

# optimize
bnds = ((0,10), (0,10))
con1 = {'type': 'ineq', 'fun': constraint}  
cons = ([con1])
solution = minimize(objective,x0,method='SLSQP', bounds=bnds, constraints=cons)
x = solution.x

# show final objective
print('Final SSE Objective: ' + str(objective(x)))

# print solution
print('x1 = ' + str(x[0]))
print('x2 = ' + str(x[1]))
