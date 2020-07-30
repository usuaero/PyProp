import numpy as np

def poly_fit(n,xx,yy,sym=False):
	x = np.copy(xx)
	y = np.copy(yy)
	start = 0
	if not sym:
		# initialize variables
		A = np.zeros((n,n))
		b = np.zeros(n)
		step = 1
	else:
		step = 2
		if n % 2 == 0:
			if n > 2:
				size = int(n / 2)
				start = 1
			else:
				size = n
				step = 1
		else:
			size = int(n / 2) + 1
		A = np.zeros( (size,size) )
		b = np.zeros( size )
	for i,e in enumerate(range(start,n,step)):
		for j,f in enumerate(range(start,n,step)):
			A[i,j] = sum( x ** float(e + f) )
		b[i] = sum( y * x ** float(e) )
	a = np.linalg.solve(A,b)
	if sym:
		start = 1
		if n % 2 == 0:
			start = 0
		a = np.insert(a,range(start,a.size),0.)
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

def rms(raw_x, raw_y, a):
	x = np.copy( raw_x )
	y = np.copy( raw_y )
	avg = np.mean(abs(y))
	l = len(x)
	func = np.zeros(l)
	e = np.zeros(l)
	e_per = np.zeros(l)
	for i in range(l):
		func[i] = poly_func(a, x[i])
		e[i] = (y[i] - func[i]) ** 2.
		e_per[i] = ((y[i] - func[i])/avg) ** 2.
	return np.sqrt(np.mean(e)), np.sqrt(np.mean(e_per))

