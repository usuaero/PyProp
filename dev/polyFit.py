#Code written by Zachary Montgomery.
#Altered by Cory Goates to handle "NaN" values.
#Altered by Cory Goates to allow the user to specify if some power coefficients
#should be set to zero.

import numpy as np

def poly_fit(n,xx,yy,forcezero=None):
	good = np.where(~(np.isnan(xx)) & ~(np.isnan(yy)))
	x = np.copy(xx[good])
	y = np.copy(yy[good])
 
	if forcezero is None:
		# initialize variables
		A = np.zeros((n,n))
		b = np.zeros(n)
		powers = range(n)
	else:
		size = n - len(forcezero)
		A = np.zeros((size,size))
		b = np.zeros(size)
		powers = np.delete(range(n), forcezero)
  
	for i,e in enumerate(powers):
		for j,f in enumerate(powers):
			A[i,j] = sum( x ** float(e + f) )
		b[i] = sum( y * x ** float(e) )
  
	a = np.linalg.solve(A,b)
 
	if forcezero is not None:
		for zero in forcezero:
			a = np.insert(a,zero,0.)
	r = r2(a,x,y)
	return a, r

def poly_func(a,x):
	f = 0.
	n = len(a)
	for i in range(n):
		f += a[i] * x ** float(i)
	return f

def r2(aa,xx,yy):
	a = np.copy(aa)
	x = np.copy(xx)
	y = np.copy(yy)
	n = len(y)
	y_ = sum(y) / float(n)
	SSt = sum((y-y_)**2.)
	f = np.zeros(n)
	for i in range(n):
		f[i] = poly_func(a,x[i])
	SSr = sum((y-f)**2.)
	return 1. - SSr / SSt

def force_good_poly_fit(n,xx,yy,tol,sym=False):
	x = np.copy(xx)
	y = np.copy(yy)
	
	a, r = poly_fit(n,x,y,sym=sym)
	x_bad = []
	y_bad = []
	print('R2 is {:7.5f}'.format(r))
	while r < tol:
		
		e_loc = 0
		l = y.size
		e_max = 0.
		for i in range(l):
			e = abs(poly_func(n,a,x[i]) - y[i])
			if e > e_max:
				e_max = e
				e_loc = i
		x_bad.append(x[e_loc])
		y_bad.append(y[e_loc])
		x = np.delete(x,e_loc)
		y = np.delete(y,e_loc)
		
		a, r = poly_fit(n,x,y)
		print('R2 is {:7.5f}'.format(r))
	
	mx = max(y)
	
	return a, x, y, x_bad, y_bad, r, mx
