## Penalty Function method.
## Loading the libraries.
import numpy as np

eps = 0.000001  #Global Var.
delta = 0.001  #Global Var.
nc = 1 ## number of constraints. 
g = np.zeros(nc) 

## This function will give the functional value at a point,
def multi_f(x):
	## Himmel Blau function!
	sum_ = pow((pow(x[0],2) + x[1] - 11),2) + pow((pow(x[1],2) + x[0] - 7),2)
	g[0] = -26.0 + pow((x[0]-5.0), 2) + pow(x[1],2);#constraints.

	for i in range(nc):
		if(g[i] < 0.0): ## meaning that the constraint is violatd.
			sum_ = sum_ + r*g[i]*g[i];

	return(sum_)

## This function will give us the gradient vector..
def grad_multi_f(grad, x_ip):
	d1 = np.zeros(M);
	d2 = np.zeros(M);
	delta_=0.001;
	
	for i in range(M):
		for j in range(M):
			d1[j]=x_ip[j]; d2[j]=x_ip[j];	

		d1[i] = d1[i] + delta_;
		d2[i] = d2[i] - delta_;

		fdiff = multi_f(d1); bdiff = multi_f(d2);		
		grad[i] = (fdiff - bdiff)/(2*delta_);
	
	return(grad)

## This function finds the alphastar <<-- UNIDIRECTIONAL SEARCH,,
def uni_search(guess, x_0, s_0, a1, b1):
	xd = np.zeros(M);	
	xw = a1 + (b1-a1)*guess;
	for i in range(M):
		xd[i] = x_0[i] + xw*s_0[i]
	return multi_f(xd)

## BOUNDING PHASE METHOD,, TO BRACKET THE ALPHASTAR.
def bracketing_(x_0, s_0):
	k=0;

	## Step__1
	## "choose an initial guess and increment(delta) for the Bounding phase method".
	del_=0.0000001; w_0=0.5;#initial guess for bounding phase.
	
	while(1):
		f0 = uni_search(w_0, x_0, s_0, 0.0, 1.0);
		fp = uni_search(w_0+del_, x_0, s_0, 0.0, 1.0);
		fn = uni_search(w_0-del_, x_0, s_0, 0.0, 1.0);
		
		if(fn >= f0):
			if(f0 >= fp):
				del_ = abs(del_); ## step__2
				break;
	
			else:
				a=w_0-del_; 
				b=w_0+del_;			
		
		elif((fn <= f0) & (f0 <= fp)):
			del_ = -1*abs(del_); 
			break;

		else:
			del_ /= 2.0;

	wn = w_0 - del_;## wn = lower limit and w_1 = upper limit.

	## Step__3
	w_1=w_0 + del_*pow(2,k); ## << exponential purterbation
	f1=uni_search(w_1, x_0, s_0, 0.0, 1.0);	

	## Step__4
	while(f1 < f0):
		k+=1;
		wn=w_0;
		fn=f0;
		w_0=w_1;
		f0=f1;
		w_1=w_0 + del_*pow(2,k);
		f1=uni_search(w_1, x_0, s_0, 0.0, 1.0);

	a=wn;
	b=w_1;
	## if del_ is + we will bracket the minima in (a, b);
	## If del_ was -,
	if(b < a):
		temp=a;a=b;b=temp;#swapp
	return(a, b)

# In order to apply Secant Method(In order to compute z). We need to also calculate the derivative of the F(alpha) function,
def f_dash(x_0, s_0, xm):
	xd = np.zeros(M);

	# using central difference.	
	for i in range(M):
		xd[i] = x_0[i] + (xm+delta)*s_0[i]; 
		fdif = multi_f(xd);

	for i in range(M): 
		xd[i] = x_0[i] + (xm-delta)*s_0[i]; 
		bdif = multi_f(xd);	
	
	f_ = (fdif - bdif)/(2*delta);
	return(f_)

## This is a function to compute the "z" used in the formula of SECANT Method..
def compute_z(x_0, s_0, x1, x2):
	z_ = x2 - ((x2-x1)*f_dash(x_0, s_0, x2))/(f_dash(x_0, s_0, x2)-f_dash(x_0, s_0,x1));
	return(z_)

def secant_minima(a, b, x_0, s_0):
	## Step__1	
	x1=a;
	x2=b;

	##Step__2 --> Compute New point.
	z = compute_z(x_0, s_0, x1, x2);
	
	##Step__3 
	while(abs(f_dash(x_0, s_0, z)) > eps):
		if(f_dash(x_0, s_0, z) >= 0):
			x2=z; # i.e. eliminate(z, x2) 		
		else:
			x1=z; # i.e. elimintae(x1, z)
		z = compute_z(x_0, s_0, x1, x2); #Step__2

	return(z)

def DFP(x_ip):
	eps1 = 0.0001;
	eps2 = 0.0001;
		
	grad0=np.zeros(M);
	grad1=np.zeros(M);
	xdiff=np.zeros(M);
	Ae=np.zeros(M);
	graddiff=np.zeros(M);	
	x_1=np.zeros(M);
	x_2=np.zeros(M);
	
	## step_1
	k=0;
	x_0=x_ip #This is the initial Guess_.  	

	## step_2	
	s_0=np.zeros(M); ##Search Direction.
	grad0 = grad_multi_f(grad0, x_ip);
	s_0 = -grad0
	
	## step_3 --> unidirectional search along s_0 from x_0. Aim -- to find alphastar and eventually x(1).
	## find alphastar such that f(x + alphastar*s) is minimum.
	a,b = bracketing_(x_0, s_0); 
	alphastar = secant_minima(a, b, x_0, s_0);
	
	for i in range(M):
		x_1[i] = x_0[i] + alphastar*s_0[i]; 
	
	grad1 = grad_multi_f(grad1, x_1); # Now computing grad(x(1))..	

	## step_4
	## Lets initialize the matrix A,
	A = np.eye(M) #Identity Matrix.
	dA = np.zeros((M,M))
	A1 = np.zeros((M,M))

	for j in range(50): # 50 is the total iterations we will perform. 		
		k+=1;
		#print("Iteration Number -- ", k)
		for i in range(M):
			xdiff[i]=x_1[i]-x_0[i]; 
			graddiff[i]= grad1[i]-grad0[i];
		
		# Now, Applying that GIANT FORMULA..
		dr=np.inner(xdiff, graddiff);
		for i in range(M):
			Ae[i]=0;
			for j in range(M):
				dA[i][j] = (xdiff[i]*xdiff[j])/dr;
				Ae[i] += (A[i][j]*graddiff[j]);
		
		dr=np.inner(graddiff, Ae)
		for i in range(M):
			for j in range(M):
				A1[i][j] = A[i][j] + dA[i][j] - (Ae[i]*Ae[j])/dr;
				A[i][j] =A1[i][j];

		## Finding the New Search Direction s_1.. -A_1*grad(x_1)
		s_1 =np.zeros(M);
		for i in range(M):
			s_1[i]=0;
			for j in range(M):
				s_1[i] += (-A1[i][j]*grad1[j]);

		## Unit vector along s_1;
		unitv = np.linalg.norm(s_1);
		for j in range(M):
			s_1[j] /= unitv;
		
		## Step__5
		# Performing Unidirectional Search along s_1 from x_1; to find x_2;)
		a,b=bracketing_(x_1, s_1); 
		alphastar = secant_minima(a, b, x_1, s_1);
		alphastar = (int)(alphastar * 1000.0)/1000.0; #Rounding alphastar to three decimal places.
		for i in range(M):
			x_0[i] = x_1[i];
			x_1[i] = x_1[i] + alphastar*s_1[i]; 

		## Step__6 --->> TERMINATION CONDITION,,
		grad1 = grad_multi_f(grad1, x_1);
		grad0 = grad_multi_f(grad0, x_0);
		
		for j in range(M):
			x_2[j] = x_1[j] - x_0[j];
		
		if((np.linalg.norm(grad1) <= eps1) or (np.linalg.norm(x_2)/np.linalg.norm(x_0)) <= eps2):
			break;

	return(x_1)

if __name__ == "__main__":
	M = 2 ## Specifies the total dimension we are working with.. //Global Var.
	c=1.55; ## factor for updating r.
	numseq = 20; ## number of sequences for penalty function method.
	r = 0.1

	x_ip = np.array([0.11,0.1])## Our initial guess,
	grad = np.zeros(M)
	sol = np.zeros(M); 
	x_2 = np.zeros(M);	
		
	## BRACKET OPERATOR PENALTY METHOD..1.1,1.1
	for k in range(numseq):
		print("\n")
		print("sequence number ", k)
		print("current function value is ", multi_f(x_ip))
		x_1 = np.copy(x_ip);		
		sol = DFP(x_ip);
		x_ip = np.copy(x_1);

		for j in range(M):
			x_2[j] = sol[j] - x_ip[j];
		
		track = np.linalg.norm(x_2)/np.linalg.norm(x_ip);
		for i in range(M): 
			print(sol[i]); 
			x_ip[i] = sol[i]

		r *= c;
		if(track <= 0.001):# TERMINATION CONDITION.
			break;
		
