#include"opt_alg.h"
//#include"matrix.cpp"
#if LAB_NO>1
double* expansion(double x0, double d, double alfa, int Nmax, matrix O)//, matrix O
{
	double* p = new double[2];
	solution X0(x0);
	solution X1(x0 + d);
	X0.fit_fun();
	X1.fit_fun();
	if (X0.y(0) == X1.y(0))
	{
		p[0] = X0.x(0);
		p[1] = X1.x(0);
		return p;
	}
	if (X1.y(0) > X0.y(0))
	{
		d *= -1;
		X1.x(0) = X0.x(0) + d;
		X1.fit_fun(X1.x);
		if (X1.y(0) >= X0.y(0))
		{
			p[0] = X1.x(0);
			p[1] = X0.x(0);
			return p;
		}
	}
	solution X2;
	int i = 1;
	while (true)
	{
		X2.x(0) = x0 + pow(alfa, i) * d;
		X2.fit_fun();
		if (i > Nmax || X0.y(0) <= X2.y(0))
			break;
		X0.x = X1.x;
		X0.fit_fun();
		X1.x = X2.x;
		++i;
	}
	if (X0.x(0) < X2.x(0))
	{
		p[0] = X0.x(0);
		p[1] = X2.x(0);
	}
	else {
		p[0] = X2.x(0);
		p[1] = X0.x(0);
	}
	return p;
}

solution fib(double a, double b, double epsilon, matrix O)
{

	int n = 2; //(b-a)/epsilon +1
	int Fib = 0;
	while (true) {
		Fib = (1 / sqrt(5)) * (pow(((1 + sqrt(5)) / (2)), n) - pow(((1 - sqrt(5)) / (2)), n));
		if (Fib > (b - a) / epsilon)
			break;
		else
			n++;
	}
	int* F = new int[n] {1, 1};
	for (int i = 2; i < n; ++i)
		F[i] = F[i - 2] + F[i - 1];
	solution A(a), B(b), C, D;
	double B0 = (double)((F[(n - 2)]) / (double)(F[(n - 1)]));
	C.x = A.x(0) - B0 * (B.x(0) - A.x(0));
	D.x = A.x(0) + B.x(0) - C.x(0);
	C.fit_fun();
	D.fit_fun();
	for (int i = 0; i <= n - 3; ++i)
	{
		cout << "Fib iteracja i: " << i << "  b-a= " << B.x(0) - A.x(0) << endl;
		if (C.y(0) < D.y(0))
		{
			B.x = D.x(0);
		}
		else {
			A.x = C.x(0);
		}
		B0 = (double)((F[n - i - 2]) / (double)(F[n - i - 1]));
		C.x = B.x(0) - B0 * (B.x(0) - A.x(0));
		D.x = A.x(0) + B.x(0) - C.x(0);
		C.fit_fun();
		D.fit_fun();
	}
	//cout << "C.x znalezione: " << C.x(0) << endl;
	//cout << "C.y znalezione: " << C.y << endl;
	cout << "Fib iteracja koncowa: " << "  b-a= " << B.x(0) - A.x(0) << endl;
	return C;
}

solution lag(double a, double b, double epsilon, double gamma, int Nmax, matrix O)
{
	solution A(a), B(b), C, D, D0;
	C.x = (a + b) / 2;
	A.fit_fun();
	B.fit_fun();
	C.fit_fun();
	double l, m;
	int i = 0;
	while (true)
	{
		cout << "Lag iteracja i: " << i << "  b-a= " << B.x(0) - A.x(0) << endl;
		/*cout << "A.x znalezione: " << A.x(0) << "  y: " << A.y(0) << endl;
		cout << "B.x znalezione: " << B.x(0) << "  y: " << B.y(0) << endl;
		cout << "C.x znalezione: " << C.x(0) << "  y: " << C.y(0) << endl;*/
		l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
		m = A.y(0) * (B.x(0) - C.x(0)) + B.y(0) * (C.x(0) - A.x(0)) + C.y(0) * (A.x(0) - B.x(0));


		if (m <= 0)
		{
			C.x = NAN;
			C.y = NAN;
			return C;
		}
		D0.x = D.x;
		D.x = l / (2 * m);
		D.fit_fun();
		//cout << "D.x znalezione: " << D.x(0) << "  y: " << D.y(0) << endl;
		if (A.x(0) < D.x(0) && D.x(0) < C.x(0))
		{
			if (D.y(0) < C.y(0))
			{
				B.x = C.x;
				C.x = D.x;
				B.fit_fun();
				C.fit_fun();
			}
			else
				A.x = D.x;
			A.fit_fun();
		}
		else if (C.x(0) < D.x(0) && D.x(0) < B.x(0))
		{
			if (D.y(0) < C.y(0))
			{
				A.x = C.x;
				C.x = D.x;
				A.fit_fun();
				C.fit_fun();
			}
			else
				B.x = D.x;
			B.fit_fun();
		}
		else
		{
			C.x = NAN;
			C.y = NAN;
			return C;
		}

		if (i > Nmax || B.x(0) - A.x(0) < epsilon || abs(D.x(0) - D0.x(0)) <= gamma)
		{

			break;
		}


		i++;
	}
	cout << "Lag iteracja koncowa: " << "  b-a= " << B.x(0) - A.x(0) << endl;
	//cout << "Lag C.x znalezione: " << C.x(0) << "  y: " << C.y(0) << endl;
	return C;
}

#endif
#if LAB_NO>2
solution HJ(matrix x0, double s, double alfa, double epsilon, int Nmax, matrix O)
{

	solution XB, XB_old, X;
	XB.x = x0;
	XB.fit_fun();
	while (true)
	{
		X = HJ_trial(XB, s);
		if (X.y(0) < XB.y(0))
		{
			while (true)
			{
				cout << "Krok HJ" << XB.x(0) << " " << XB.x(1) << endl;
				XB_old = XB;
				XB = X;
				//cout << " XB.x1: " << XB.x(0) << " XB.x1: " << XB.x(1) << endl;
				X.x = 2.0 * XB.x - XB_old.x;
				//cout << "NEW XB.x1: " << XB.x(0) << " XB.x1: " << XB.x(1) << endl;
				//X.x =  XB.x - XB_old.x;
				X.fit_fun();
				X = HJ_trial(X, s);
				if (X.y(0) >= XB.y(0))
					break;
				if (solution::f_calls > Nmax)
					return XB;
			}
		}
		else
			s *= alfa;
		if (s< epsilon || solution::f_calls>Nmax)
			return XB;
	}
}

solution HJ_trial(solution XB, double s, matrix O)
{
	int* n = get_size(XB.x);
	matrix D(n[0], n[0]);
	for (int i = 0; i < n[0]; i++)
		D(i, i) = 1;
	solution X;
	for (int i = 0; i < n[0]; ++i)
	{
		X.x = XB.x + s * D[i];
		X.fit_fun();
		if (X.y < XB.y)
			XB = X;
		else
		{
			X.x = XB.x - s * D[i];
			X.fit_fun();
			if (X.y < XB.y)
				XB = X;
		}
	}
	return XB;
}

solution Rosen(matrix x0, matrix s0, double alfa, double beta, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	matrix l(n[0], 1), p(n[0], 1), s(s0);
	matrix D(n[0], n[0]);
	for (int i = 0; i < n[0]; i++)
		D(i, i) = 1;

	solution X, Xt;
	X.x = x0;
	X.fit_fun();
	while (true)
	{
		for (int i = 0; i < n[0]; ++i)
		{
			cout << "Krok Ros" << X.x(0) << " " << X.x(1) << endl;
			Xt.x = X.x + s(i) * D[i];
			Xt.fit_fun();
			if (Xt.y(0) < X.y(0))
			{
				X = Xt;
				l(i) += s(i);
				s(i) *= alfa;
			}
			else
			{
				p(i) = p(i) + 1;
				s(i) *= -beta;
			}
		}
		bool change = true;
		for (int i = 0; i < n[0]; ++i)
			if (p(i) == 0 || l(i) == 0)
			{
				change = false;
				break;
			}
		if (change)
		{
			matrix Q(n[0], n[0]), v(n[0], 1);
			for (int i = 0; i < n[0]; ++i)
				for (int j = 0; j <= i; ++j)
					Q(i, j) = l(i);
			Q = D * Q;
			v = Q[0] / norm(Q[0]);
			D = set_col(D, v, 0);
			for (int i = 1; i < n[0]; ++i)
			{
				matrix temp(n[0], 1);
				for (int j = 0; j < i; ++j)
					temp = temp + (trans(Q[i]) * D[j]) * D[j];
				v = Q[i] - temp;
				D = set_col(D, v, i);
			}
			s = s0;
			l = matrix(n[0], 1);
			p = matrix(n[0], 1);
		}
		double max_s = abs(s(0));
		for (int i = 1; i < n[0]; ++i)
			if (max_s < abs(s(i)))
				max_s = abs(s(i));
		if (max_s < epsilon || solution::f_calls > Nmax)
			return X;
	}
}
#endif
#if LAB_NO>3
solution pen(matrix x0, double c0, double dc, double epsilon, int Nmax, matrix O)
{
	double alfa = 1, beta = 0.5, gama = 2, delta = 0.5, s = 0.5;
	matrix A(new double[2]{ c0,O(0) }, 2);
	solution X, X1;
	X.x = x0;
	while (true)
	{
		X1 = sym_NM(X.x, s, alfa, beta, gama, delta, epsilon, Nmax, A);
		if (A(0) < epsilon || solution::f_calls > Nmax || A(0) > 1 / epsilon)
			return X1;
		A(0) *= dc;
		X = X1;
	}
}

solution sym_NM(matrix x0, double s, double alfa, double beta, double gama, double delta, double epsilon, int Nmax, matrix O)
{
	/**/
	int* n = get_size(x0);
	matrix D(n[0], n[0]);
	for (int i = 0; i < n[0]; i++)
		D(i, i) = 1;


	int N = n[0] + 1;
	solution* S = new solution[N];
	S[0].x = x0;
	S[0].fit_fun(O);
	for (int i = 1; i < N; ++i)
	{
		S[i].x = S[0].x + s;
		S[i].fit_fun(O);
	}
	solution p_o, p_e, p_z;
	matrix p_sr;
	int i_min, i_max;
	while (true)
	{
		i_min = 0;
		i_max = 0;
		for (int i = 1; i < N; ++i)
		{
			if (S[i].y <= S[i_min].y)
				i_min = i;
			if (S[i].y >= S[i_max].y)
				i_max = i;
		}
		p_sr = matrix(n[0], 1);
		for (int i = 0; i < N; ++i)
			if (i != i_max)
				p_sr = p_sr + S[i].x;
		p_sr = p_sr / N;
		p_o.x = p_sr + alfa * (p_sr - S[i_max].x);
		p_o.fit_fun(O);
		if (S[i_min].y <= p_o.y && S[i_max].y > p_o.y)
			S[i_max] = p_o;
		else if (S[i_min].y > p_o.y)
		{
			p_e.x = p_sr + gama * (p_o.x - p_sr);
			p_e.fit_fun(O);
			if (p_e.y < p_o.y)
				S[i_max] = p_e;
			else
				S[i_max] = p_o;
		}
		else
		{
			p_z.x = p_sr + beta * (S[i_max].x - p_sr);
			p_z.fit_fun(O);
			if (p_z.y < S[i_max].y)
				S[i_max] = p_z;
			else
			{
				for (int i = 0; i < N; ++i)
					if (i != i_min)
					{
						S[i].x = delta * (S[i].x + S[i_min].x);
						S[i].fit_fun(O);
					}
			}
		}
		double max_s = norm(S[0].x - S[i_min].x);
		for (int i = 1; i < N; ++i)
			if (max_s < norm(S[i].x - S[i_min].x))
				max_s = norm(S[i].x - S[i_min].x);
		if (solution::f_calls > Nmax || max_s < epsilon)
			return S[i_min];
	}
}
#endif
#if LAB_NO>4
solution SD(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b;
	while (true)
	{
		X.grad();
		d = -X.g;
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		if (h0 < 0)
		{
			b = compute_b(X.x, d, limits);
			h = golden(0, b, epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
			X1.x = X.x + h0 * d;
		if (norm(X1.x - X.x) < epsilon ||
			solution::g_calls > Nmax ||
			solution::f_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution CG(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b, beta;
	X.grad();
	d = -X.g;
	while (true)
	{
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		if (h0 < 0)
		{
			b = compute_b(X.x, d, limits);
			h = golden(0, b, epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
			X1.x = X.x + h0 * d;
		if (norm(X1.x - X.x) < epsilon ||
			solution::g_calls > Nmax ||
			solution::f_calls > Nmax)
		{
			X1.fit_fun();
			return X1;
		}
		X1.grad();
		beta = pow(norm(X1.g), 2) / pow(norm(X.g), 2);
		d = -X1.g + beta * d;
		X = X1;
	}
}

solution Newton(matrix x0, double h0, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	solution X, X1;
	X.x = x0;
	matrix d(n[0], 1), P(n[0], 2), limits = O;
	solution h;
	double b;
	while (true)
	{
		X.grad();
		X.hess();
		d = -inv(X.H) * X.g;
		P = set_col(P, X.x, 0);
		P = set_col(P, d, 1);
		if (h0 < 0)
		{
			b = compute_b(X.x, d, limits);
			h = golden(0, b, epsilon, Nmax, P);
			X1.x = X.x + h.x * d;
		}
		else
			X1.x = X.x + h0 * d;
		if (norm(X1.x - X.x) < epsilon ||
			solution::g_calls > Nmax ||
			solution::f_calls > Nmax ||
			det(X.H) == 0)
		{
			X1.fit_fun();
			return X1;
		}
		X = X1;
	}
}

solution golden(double a, double b, double epsilon, int Nmax, matrix O)
{
	double alfa = (sqrt(5) - 1) / 2;
	solution A, B, C, D;
	A.x = a;
	B.x = b;
	C.x = B.x - alfa * (B.x - A.x);
	C.fit_fun(O);
	D.x = A.x + alfa * (B.x - A.x);
	D.fit_fun(O);
	while (true)
	{
		if (C.y < D.y)
		{
			B = D;
			D = C;
			C.x = B.x - alfa * (B.x - A.x);
			C.fit_fun(O);
		}
		else
		{
			A = C;
			C = D;
			D.x = A.x + alfa * (B.x - A.x);
			D.fit_fun(O);
		}
		if (solution::f_calls > Nmax || B.x - A.x < epsilon)
		{
			A.x = (A.x + B.x) / 2.0;
			A.fit_fun(O);
			return A;
		}
	}
}

double compute_b(matrix x, matrix d, matrix limits)
{
	int* n = get_size(x);
	double b = 1e9, bi;
	for (int i = 0; i < n[0]; ++i)
	{
		if (d(i) == 0)
			bi = 1e9;
		else if (d(i) > 0)
			bi = (limits(i, 1) - x(i)) / d(i);
		else
			bi = (limits(i, 0) - x(i)) / d(i);
		if (b != bi)
			b = bi;
	}
	return b;
}
#endif
#if LAB_NO>5
solution Powell(matrix x0, double epsilon, int Nmax, matrix O)
{
	int* n = get_size(x0);
	matrix D = ident_mat(n[0]), A(n[0], 3), limits(n[0], 2);
	limits = set_col(limits, O[0], 0);
	limits = set_col(limits, O[1], 1);
	A(0, 2) = O(0, 2);
	solution X, P, h;
	X.x = x0;
	double* ab;
	while (true)
	{
		P = ? ;
		for (int i = 0; i < ? ; ++i)
		{
			A = set_col(A, P.x, 0);
			A = set_col(A, D[i], 1);
			ab = compute_ab(? , ? , limits);
			h = golden(? , ? , epsilon, Nmax, A);
			P.x = ? ;
		}
		if (? )
		{
			P.fit_fun();
			return P;
		}
		for (int i = 0; i < n[0] - 1; ++i)
			D = ? ;
		D = ? ;
		A = set_col(A, P.x, 0);
		A = set_col(A, D[n[0] - 1], 1);
		ab = compute_ab(? , ? , limits);
		h = golden(? , ? , epsilon, Nmax, A);
		X.x = ? ;
	}
}

double* compute_ab(matrix x, matrix d, matrix limits)
{
	int* n = get_size(x);
	double* ab = new double[2]{ -1e9,1e9 };
	double ai, bi;
	for (int i = 0; i < n[0]; ++i)
	{
		if (d(i) == 0)
		{
			ai = ? ;
			bi = ? ;
		}
		else if (d(i) > 0)
		{
			ai = ? ;
			bi = ? ;
		}
		else
		{
			ai = ? ;
			bi = ? ;
		}
		if (? )
			ab[0] = ai;
		if (? )
			ab[1] = bi;
	}
	return ab;
}
#endif
#if LAB_NO>6
solution EA(int N, matrix limits, double epsilon, int Nmax, matrix O)
{
	int mi = 20, lambda = 40;
	solution* P = new solution[mi + lambda];
	solution* Pm = new solution[mi];
	random_device rd;
	default_random_engine gen;
	gen.seed(static_cast<unsigned int>(chrono::system_clock::now().time_since_epoch().count()));
	normal_distribution<double> distr(0.0, 1.0);
	matrix IFF(mi, 1), temp(N, 2);
	double r, s, s_IFF;
	double tau = ? , tau1 = ? ;
	int j_min;
	for (int i = 0; i < ? ; ++i)
	{
		P[i].x = matrix(N, 2);
		for (int j = 0; j < N; ++j)
		{
			P[i].x(j, 0) = ? ;
			P[i].x(j, 1) = ? ;
		}
		P[i].fit_fun();
		if (P[i].y < epsilon)
			return P[i];
	}
	while (true)
	{
		s_IFF = 0;
		for (int i = 0; i < ? ; ++i)
		{
			IFF(i) = 1 / P[i].y(0);
			s_IFF += IFF(i);
		}
		for (int i = 0; i < ? ; ++i)
		{
			r = ? ;
			s = 0;
			for (int j = 0; j < ? ; ++j)
			{
				s += ? ;
				if (? )
				{
					P[mi + i] = ? ;
					break;
				}
			}
		}
		for (int i = 0; i < ? ; ++i)
		{
			r = distr(gen);
			for (int j = 0; j < N; ++j)
			{
				P[mi + i].x(j, 1) *= ? ;
				P[mi + i].x(j, 0) += ? ;
			}
		}
		for (int i = 0; i < ? ; i += 2)
		{
			r = ? ;
			temp = P[mi + i].x;
			P[mi + i].x = ? ;
			P[mi + i + 1].x = ? ;
		}
		for (int i = 0; i < ? ; ++i)
		{
			P[mi + i].fit_fun();
			if (P[mi + i].y < epsilon)
				return P[mi + i];
		}
		for (int i = 0; i < ? ; ++i)
		{
			j_min = 0;
			for (int j = 1; j < ? ; ++j)
				if (P[j_min].y > P[j].y)
					j_min = j;
			Pm[i] = P[j_min];
			P[j_min].y = 1e10;
		}
		for (int i = 0; i < ? ; ++i)
			P[i] = Pm[i];
		if (? )
			return P[0];
	}
}
#endif
