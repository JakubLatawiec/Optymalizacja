//Do not edit the code below (unless you know what you are doing)

#include"ode_solver.h"
#include <fstream>

matrix *solve_ode(double t0, double dt, double tend, const matrix &Y0, matrix P)
{
	int N = static_cast<int>(floor((tend - t0) / dt) + 1);
	if (N < 2)
		throw "The time interval is defined incorrectly";
	int *s = get_size(Y0);
	if (s[1] != 1)
		throw "Initial condition must be a vector";
	int n = s[0];
	delete[]s;
	matrix *S = new matrix[2]{ matrix(N), matrix(n,N) };
	S[0](0) = t0;
	for (int i = 0; i < n; ++i)
		S[1](i, 0) = Y0(i);
	matrix k1(n), k2(n), k3(n), k4(n);
	for (int i = 1; i < N; ++i)
	{
		S[0](i) = S[0](i - 1) + dt;
		k1 = dt*diff(S[0](i - 1), S[1][i - 1], P);
		k2 = dt*diff(S[0](i - 1) + 0.5*dt, S[1][i - 1] + 0.5*k1, P);
		k3 = dt*diff(S[0](i - 1) + 0.5*dt, S[1][i - 1] + 0.5*k2, P);
		k4 = dt*diff(S[0](i - 1) + dt, S[1][i - 1] + k3, P);
		for (int j = 0; j < n; ++j)
			S[1](j, i) = S[1](j, i - 1) + (k1(j) + 2 * k2(j) + 2 * k3(j) + k4(j)) / 6;
	}


	///
	cout << P(0) << ";" << P(1) << endl;

	S[1] = trans(S[1]);
	return S;
}

//You can edit the following code

matrix diff(double t, const matrix &Y, matrix P)
{
#if LAB_NO==1 
	double m = 5, d = 1.5, k = 1, F = 0;
	matrix dY(Y);
	dY(0) = Y(1);
	dY(1) = (F - d * Y(1) - k * Y(0)) / m;
	return dY;

#elif LAB_NO==2

	matrix dY(Y);
	//Y(O)=VA, Y(1)=VB, Y(2)=TB

	double a = 0.96;			//wsp. - lepkoœæ
	double b = 0.63;			//wsp. - zwê¿enie strumienia
	double g = 9.81;			//przysp. ziemskie
	double PA = 1;				//pole pow. zbiornika A
	double PB = 1;				//pole pow. zbiornika B
	double DB = 0.00365665;		//wielkoœæ otworu w zbiorniku B
	double Fin = 0.01;			//iloœæ wody która wp³ywa do B (z rury)
	double Tin = 10;			//temp  wody ^
	double TA = 90.0;			//temp wody w zbiorniku A
	double DA = P(0);			//wielkosc otworu w zbiorniku A 

	//ile wody wylewa sie z A:
	double FAout = Y(0) > 0 ? a * b * DA * sqrt(2 * g * Y(0) / PA) : 0;
	//ile wody wylewa sie z B:
	double FBout = Y(1) > 0 ? a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;
	
	//VA':
	dY(0) = -FAout;
	//VB':
	dY(1) = FAout + Fin - FBout;
	//TB':
	dY(2) = FAout / Y(1) * (TA - Y(2)) + Fin / Y(1) * (Tin - Y(2));
	return dY;

#elif LAB_NO == 3 && LAB_PART == 2

	double mr = 1.0;				//masa ramienia
	double mc = 10.0;				//masa ciezarka
	double l = 0.5;					//dl. ramienia
	double alfa_ref = 3.14159265358979323846;		//pi rad
	double omega_ref = 0.0;			//0 rad/s
	double b = 0.5;					//wsp. tarcia

	//moment bezwladnosci:
	double I = (mr * l * l) / 3 + mc * l * l;

	double k1 = P(0);			//wsp. wzmocnienia
	double k2 = P(1);			//przesylane w P

	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (k1 * (alfa_ref - Y(0)) + k2 * (omega_ref - Y(1)) - b * Y(1)) / I;
	//przmieszczenie i predkosc:
	//cout << Y(0) << " " << Y(1) << endl;
	return dY;

#else

	matrix dY;
	return dY;

#endif
}