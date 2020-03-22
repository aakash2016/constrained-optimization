/*
AAKASH AGRAWAL 160122001.
Phase3 - ME609
### CONSTRAINED OPTIMIZATION ####
Unidirectional Search = Bounding Phase + Secant Method. 
Multivariable Optimization = Variable Metric Method/DFP METHOD.
Constrained Optimization = Penalty function method. <Bracket Operator Penalty Method>.
*/

#include<iostream>
#include<math.h>
#include<tuple>
#include <numeric>
#include<vector>
#include<fstream>

using namespace std;

//arrays cant be returned from a function in cpp so we will use vector throughout.
double eps = 0.000001;//Global Var.
double delta = 0.001;//Global Var.
double del, r=1.15;//GLobal Var.
int M_, M, nc=2;//M Specifies the total dimension we are working with.. //Global Var.,,//nc number of constrains.
int restart;
vector<double> sigma(nc,0);//Global vector.

vector<double> constraint(const vector<double> &x){
	vector<double> g(nc);
	//basic himmelblau function.
	//g[0] = -26.0 + pow((x[0]-5.0), 2) + pow(x[1],2);//constrints.
	
	//problem -- II
	g[0] = -(pow(x[0],2) - x[1] + 1.0);
	g[1] = -(1 - x[0] + pow((x[1]-4),2));
	/*
	//problem -- I
	sum = pow((x[0]-10),3) + pow((x[1]-20),3);
	g[0]= pow((x[0]-5),2) + pow((x[1]-5),2) - 100.0;
	g[1]= -(pow((x[0]-5),2) + pow((x[1]-5),2) - 82.81);



	
	//problem -- III // solve for 8 variables. M=8;
	g[0] = 1 - (x[3]+ x[5])*0.0025;
	g[1]= 1 - (x[4]+ x[6] - x[3])*0.0025;
	g[2]=1 - (x[7]- x[5])*0.01;
	g[3]=-(100*x[0] - x[0]*x[5] + 833.33252*x[3] - 83333.333);
	g[4]=-(x[1]*x[3] - x[1]*x[6] - 1250*x[3] + 1250*x[4]);
	g[5]=-(x[2]*x[4] - x[2]*x[7] - 2500*x[4] + 1250000);
	*/
	return g;
}

//This Function will give us the Function value at a point//
double multi_f(const vector<double> &x){
	double sum;
	vector<double> g(nc);
	// Check for Himmel Blau function 
	//sum = pow((pow(x[0],2) + x[1] - 11),2) + pow((pow(x[1],2) + x[0] - 7),2);
	
	//problem -- II
	double pi = atan(1)*4;
	sum = -(pow(sin(2*pi*x[0]), 3)*sin(2*pi*x[1]))/pow(x[0],3)*(x[0] + x[1]);//duality consition.	

	/*
	//problem -- I
	//sum = pow((x[0]-10),3) + pow((x[1]-20),3);
	
		
	//problem -- III // solve for 8 variables. M=8;
	for(int i=0; i<M ;i++){
		sum += x[i];	
	}
	*/	
	
	g = constraint(x);
	//for method of multipliers/
	for(int i=0; i<nc; i++){	
		if(g[i] < 0.0){
			sum = sum + r*pow((g[i]+sigma[i]),2) - r*pow(sigma[i],2);
		}
	}

	return sum;
}

//This function will give us the gradient vector..
void grad_multi_f(vector<double> &grad, const vector<double> &x_ip){
	double fdiff, bdiff;
	vector<double> dumb1(M);
	vector<double> dumb2(M);
	double delta_=0.000001;
	for (int i = 1; i <= M; ++i){dumb1[i-1]=0; dumb2[i-1]=0;}
	
	for(int i=1; i<=M; i++){
		for (int j=0; j<M; j++){dumb1[j]=x_ip[j]; dumb2[j]=x_ip[j];}	

		dumb1.at(i-1) = dumb1[i-1] + delta_;
		dumb2.at(i-1) = dumb2[i-1] - delta_;

		fdiff = multi_f(dumb1); bdiff = multi_f(dumb2);		
		grad[i-1] = (fdiff - bdiff)/(2*delta_);
	}
}

// In order to apply Secant Method(In order to compute z). We need to also calculate the derivative of the F(alpha) function,
double f_dash(const vector<double> &x_0, const vector<double> &s_0, double xm){
	double f_;// using central difference.
	double fdif, bdif;
	vector<double> xd(M);	
	for (int i = 0; i < M; ++i){xd[i] = x_0[i] + (xm+delta)*s_0[i]; fdif = multi_f(xd);}
	for (int i = 0; i < M; ++i){xd[i] = x_0[i] + (xm-delta)*s_0[i]; bdif = multi_f(xd);}	
	
	f_ = (fdif - bdif)/(2*delta);
	return f_;
}

// This is a function to compute the "z" used in the formula of SECANT Method..
double compute_z(const vector<double> &x_0, const vector<double> &s_0, double x1, double x2){
	double z_;
	z_ = x2 - ((x2-x1)*f_dash(x_0, s_0, x2))/(f_dash(x_0, s_0, x2)-f_dash(x_0, s_0,x1));
	return z_;
}

//This function finds the alphastar <<-- UNIDIRECTIONAL SEARCH,,
double uni_search(double guess, const vector<double> &x_0, const vector<double> &s_0, double a1, double b1){
	vector<double> xd(M);	
	double xw = a1 + (b1-a1)*guess;
	for (int i = 0; i < M; ++i){xd[i] = x_0[i] + xw*s_0[i];}
	return multi_f(xd);
}

//BOUNDING PHASE METHOD,, TO BRACKET THE ALPHASTAR.
std::tuple<double, double> bracketing_(vector<double> x_0, vector<double> s_0){
	double init_guess, del; // wn = lower limit and w_1 = upper limit.
	double w_0, w_1, temp, wn;
	double f0, f1, fp, fn;
	double a,b;
	int k=0;

	//Step__1
	//"choose an initial guess and increment(delta) for the Bounding phase method".
	//cout << "Enter the initial guess for bracketing method and an increment delta" << endl; 	
	//cin >> init_guess; 
	//cin >> del;
	del=0.0000001; w_0=0.5;	
	

	while(1){
		f0 = uni_search(w_0, x_0, s_0, 0.0, 1.0);
		fp = uni_search(w_0+del, x_0, s_0, 0.0, 1.0);
		fn = uni_search(w_0-del, x_0, s_0, 0.0, 1.0);
		
		if(fn >= f0){
			if(f0 >= fp){del = abs(del); break;}	
			else{a = w_0-del; b=w_0+del;}			
		}
		
		else if((fn <= f0) && (f0 <= fp)){del = -1*abs(del); break;}
		else{del /= 2.0;}
	}

	/*
	//calculate the function values,
	f0 = uni_search(init_guess, x_0, s_0);
	fp = uni_search(init_guess+del, x_0, s_0);
	fn = uni_search(init_guess-del, x_0, s_0);

	//Step__2
	if(fn > fp){del = abs(del);}// the increment is positive.
	else if(fn < fp){del = -1*abs(del);} // the increment is negative.
	*/

	wn = w_0 - del;

	//Step__3
	w_1=w_0 + del*pow(2,k); // << exponential purterbation
	f1=uni_search(w_1, x_0, s_0, 0.0, 1.0);	

	//Step__4
	while(f1 < f0){
		k+=1;
		wn=w_0;
		fn=f0;
		w_0=w_1;
		f0=f1;
		w_1=w_0 + del*pow(2,k);
		f1=uni_search(w_1, x_0, s_0, 0.0, 1.0);
	}

	//cout << k << endl;
	a=wn;
	b=w_1;
	// if del is + we will bracket the minima in (a, b);
	// If del was -,
	if(b < a){temp=a;a=b;b=temp;}//swapp
	return make_tuple(a, b);
}

double secant_minima(double a, double b, const vector<double> &x_0, const vector<double> &s_0){
	double z,x1,x2;
	double t_factor=eps;

	//Step__1	
	x1=a;
	x2=b;

	//Step__2 --> Compute New point.
	z = compute_z(x_0, s_0, x1, x2);
	
	//Step__3 
	while(abs(f_dash(x_0, s_0, z)) > t_factor){
		if(f_dash(x_0, s_0, z) >= 0){x2=z;} // i.e. eliminate(z, x2) 		
		else{x1=z;} // i.e. elimintae(x1, z)
		z = compute_z(x_0, s_0, x1, x2); //Step__2
	}
	
	return z;
}

double accurate(double a, double b, const vector<double> &x_0, const vector<double> &s_0){
	double xstar = a;
	int maxfun = 10000;
	double gold = (sqrt(5.0) - 1.0)/2.0;
	double aw = 0.0; 
	double bw =1.0;
	double lw=1.0;
	int k=1, ic =0;;
	double w1prev = gold;
	double w2prev = 1.0- gold;
	double fw1 = uni_search(w1prev, x_0, s_0, a, b);
	double fw2 = uni_search(w2prev, x_0, s_0, a, b);
	double w1;
	double w2;
	while(1){
		w1 = w1prev;
		w2 = w2prev;
		if(ic == 1){fw2 = fw1; fw1 = uni_search(w1, x_0, s_0, a, b);}
		else if(ic == 2){fw1 = fw2; fw2 = uni_search(w2, x_0, s_0, a, b);}
		else if(ic == 3){fw1 = uni_search(w1, x_0, s_0, a, b); fw2 = uni_search(w2, x_0, s_0, a, b);}

		//Region Elimination Rule.
		if(fw1 < fw2){ic =1; aw=w2; lw=bw-aw; w1prev = aw + gold*lw; w2prev = w1;}
		else if(fw2 < fw1){ic =2; bw = w1; lw=bw-aw; w1prev =w2; w2prev = bw - gold*lw;}
		else{ic=3; aw=w2;bw=w1;lw=bw-aw; w1prev = aw + gold*lw; w2prev = bw - gold*lw;}
		if(abs(lw) <= 0.0001){
			xstar= a + (b-a)*(aw+bw)/2.0; 
			break;
		} 
	}
	return xstar;
}

double norm(const vector<double> &x){
	return pow(inner_product(begin(x), end(x), begin(x), 0.0), 0.5);
}

vector<double> DFP(const vector<double> &x_ip){
	double eps1 = 0.000001;
	double eps2 = 0.000001;
	double x, a, b, alphastar, dr;	
	vector<double> grad0(M);
	vector<double> grad1(M);
	vector<double> xdiff(M);
	vector<double> Ae(M);
	vector<double> graddiff(M);	
	vector<double> x_1(M);
	vector<double> x_2(M);
	int k;
	
	//step_1
	k=0;
	vector<double> x_0(x_ip); //This is the initial Guess_. //copy of x_ip 	

	//step_2	
	vector<double> s_0(M); //Search Direction.
	grad_multi_f(grad0, x_ip);
	for (int i = 0; i < M; ++i){s_0[i] = -grad0[i];}
	
	//step_3 --> unidirectional search along s_0 from x_0. Aim -- to find alphastar and eventually x(1).
	// find alphastar such that f(x + alphastar*s) is minimum.
	//cout << "Performing Undirectional Search along s_0 from x_0" << endl;	
	tie(a,b)=bracketing_(x_0, s_0); 

	//cout << endl;
	//cout << "Iteration Number -- " << k << endl;
	alphastar = secant_minima(a, b, x_0, s_0);
	//alphastar = (int)(alphastar * 1000.0)/1000.0; //Rounding alphastar to three decimal places.
	for (int i=0; i < M; ++i){x_1[i] = x_0[i] + alphastar*s_0[i];} //cout << x_1[i] << endl;
	grad_multi_f(grad1, x_1);//Now computing grad(x(1))..	

	//step_4
	//Lets initialize the matrix A,
	vector<vector<double>> A(M, vector<double>(M,0));
	for(unsigned int t = 0; t < M; t++)
		A[t][t] = 1;//Identity Matrix.

	vector<vector<double>> dA(M, vector<double>(M,0));
	vector<vector<double>> A1(M, vector<double>(M,0));
	
	for(int j=0; j<50; j++){// 50 is the total iterations we will perform. 		
		k+=1;
		//cout << endl;
		//cout << "Iteration Number -- " << k << endl;
		for (int i=0; i < M; ++i){
			xdiff[i]=x_1[i]-x_0[i]; 
			graddiff[i]= grad1[i]-grad0[i];
		}
		
		//Now, Applying that GIANT FORMULA..
		dr=inner_product(begin(xdiff), end(xdiff), begin(graddiff), 0.0);
		for(int i = 0; i< M; i++){
			Ae[i]=0;
			for(int j = 0; j < M; j++){
				dA[i][j] = (xdiff[i]*xdiff[j])/dr;
				Ae[i] += (A[i][j]*graddiff[j]);
			}
		}
		
		dr=inner_product(begin(graddiff), end(graddiff), begin(Ae), 0.0);
		for(int i = 0; i< M; i++){
			for(int j = 0; j < M; j++){
				A1[i][j] = A[i][j] + dA[i][j] - (Ae[i]*Ae[j])/dr;
				A[i][j] =A1[i][j];
				//cout << A[i][j] << endl;
			}
		}

		//Finding the New Search Direction s_1.. -A_1*grad(x_1)
		vector<double> s_1(M);
		for(int i = 0; i< M; i++){
			s_1[i]=0;
			for(int j = 0; j < M; j++){
				s_1[i] += (-A1[i][j]*grad1[j]);
			}
			//cout << s_1[i] << endl;
		}

		//Unit vector along s_1;
		double unitv = norm(s_1);
		for(int j = 0; j < M; j++){
			s_1[j] /= unitv;
		}	
		
		//Step__5
		//Performing Unidirectional Search along s_1 from x_1; to find x_2;)
		tie(a,b)=bracketing_(x_1, s_1); 
		//cout << a << "  " << b << endl;
		alphastar = secant_minima(a, b, x_1, s_1);
		alphastar = (int)(alphastar * 1000.0)/1000.0; //Rounding alphastar to three decimal places.
		for (int i=0; i < M; ++i){
			x_0[i] = x_1[i];
			x_1[i] = x_1[i] + alphastar*s_1[i]; 
			//cout << x_1[i] << endl;
		}

		//Step__6 --->> TERMINATION CONDITION,,
		grad_multi_f(grad1, x_1);
		grad_multi_f(grad0, x_0);
		
		for(int j = 0; j < M; j++){x_2[j] = x_1[j] - x_0[j];}
		if((norm(grad1) <= eps1) || (norm(x_2)/norm(x_0)) <= eps2){break;}
		
		//RESTART
		if((restart == 1) && (j ==3)){//RESTART OPTION. --> restarting after 4th iteration.
			cout << endl << "We are Restarting the DFP" << endl;
			for(int i = 0; i< M; i++){
				for(int j = 0; j < M; j++){
					A[i][j] = 0;
				}
			}
			
			for(int t = 0; t < M; t++)
				A[t][t] = 1;//Identity Matrix.	
			
			k=0;		
		}
	}

	/*
	// Printing the A matrix -- HESSIAN,
	for(int i = 0; i< M; i++){
		for(int j = 0; j < M; j++){
			cout << (int)(A[i][j]*1000.0)/1000.0 << "\t";
		}
		cout << endl;		
	}

	cout << endl;
	cout << multi_f(x_1) << endl;
	cout << k << endl;
	*/

	return x_1;
}

int main(){
	ifstream file_("read_p3.txt"); //ifstream :: stream class to read from files,
	cout << "#################### DFP - Varible Metric Method ##################" << endl;
	double track=0, pi = atan(1)*4;
	double epss = 0.00001;
	int n, M_, numseq; //numseq = number of sequences for penalty function method.
	int restart; // 0 or 1(TRUE).
	cout << "Enter the number of variables you want to work with --> ";
	cin >> n; 
	::M = n;

	//Restart Option.
	::restart =0;
	::M_=50;//total number of iterations.
	numseq = 10;

	vector<double> x_ip(M);
	if(file_.is_open()){
		for(int i = 1; i <= M; ++i){file_ >> x_ip[i-1];} // xi can take values -10 to 10.
		file_.close();
	}
	
	vector<double> sol(M);//0.836, 2.940 //0.829, 2.933
	vector<double> x_2(M);	
	vector<double> g(nc);

	//METHOD OF MULTIPLIERS.
	// r remains same throughout simulation.
	ofstream f("phase_me609_10.txt"); //ofstream :: stream class to write on files,
	if(f.is_open()){
		for(int k = 1; k <= numseq; k++){
			cout << endl;
			cout << "Sequence number " << k << endl;
			cout << "initial function value is " << multi_f(x_ip) << endl;//P(x(0), sigma(0))
			f << multi_f(x_ip) << endl;
			g = constraint(x_ip);
			sol = DFP(x_ip);

			for(int j = 0; j < M; j++){
				cout << sol[j] << endl;
				x_ip[j] = sol[j]; //1st sequence.
			}

			if(abs(track-multi_f(sol)) <= epss){break;}// TERMINATION CONDITION.
			track = multi_f(sol);//P(x(1), sigma(0))//P(x(2), sigma(1)).
			
			//step--5 //updating the multipliers.
			for(int i=0; i<nc; i++){	
				if(g[i] < 0.0){
					::sigma[i] += g[i];
				}
			}
		}
		f.close();
	}
	return 0;
}


