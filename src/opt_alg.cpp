#include "configuration.h"
#include"opt_alg.h"
#include "environment.h"


solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int i = 0;
		solution xi_sol, xi_next_sol;
		double xi, xi_next;

		xi = x0;
		xi_next = xi + d;

		xi_sol.x = xi;
		xi_sol.fit_fun(ff, ud1);

		xi_next_sol.x = xi_next;
		xi_next_sol.fit_fun(ff, ud1);

		if (xi_next_sol.y == xi_sol.y)
			return new double[3] {xi, xi_next, (double)solution::f_calls};

		if (xi_next_sol.y > xi_sol.y)
		{
			d = -d;
			xi_next = xi + d;
			xi_next_sol.x = xi_next;
			xi_next_sol.fit_fun(ff, ud1);

			if (xi_next_sol.y >= xi_sol.y)
				return new double[3] {xi_next, xi - d, (double)solution::f_calls};
		}

		solution::clear_calls();
		double xi_prev{};
		double f_xi = m2d(xi_sol.y);
		do
		{
			if (solution::f_calls > Nmax)
			{
				xi_next_sol.flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

			++i;
			xi_next = xi + pow(alpha, i) * d;

			xi_next_sol.x = xi_next;
			xi_next_sol.fit_fun(ff, ud1);

			if (!(f_xi > xi_next_sol.y))
				break;

			xi_prev = xi;
			xi = xi_next;
			f_xi = m2d(xi_next_sol.y);

		} while (true);

		if (d > 0)
			return new double[3] {xi_prev, xi_next, (double)solution::f_calls};

		return new double[3] {xi_next, xi_prev, (double)solution::f_calls};
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		std::stringstream ss;
    
		std::vector<double> sigma = { 1, 1 };
		double ratio = (b - a) / epsilon;
		while (true)
		{
			if (sigma.back() > ratio)
				break;

			sigma.push_back(sigma[sigma.size() - 1] + sigma[sigma.size() - 2]);
		}	
		int k = sigma.size() - 1;
		
		double a0 = a;
		double b0 = b;
		double c0 = b0 - sigma[k - 1] / sigma[k] * (b0 - a0);
		double d0 = a0 + b0 - c0;

		solution c_sol, d_sol;
		for (int i = 0; i <= k - 3; ++i)
		{
			c_sol.x = c0;
			c_sol.fit_fun(ff, ud1);

			d_sol.x = d0;
			d_sol.fit_fun(ff, ud1);

			if (c_sol.y < d_sol.y)
				b0 = d0;
			else
				a0 = c0;

			c0 = b0 - sigma[k - i - 2] / sigma[k - i - 1] * (b0 - a0);
			d0 = a0 + b0 - c0;
      
			if (SAVE_CHART_DATA)
				ss << i << ";" << b0 - a0 << ";\n";
		}

		solution Xopt;
		Xopt.x = c0;
		Xopt.fit_fun(ff, ud1);

		if (SAVE_CHART_DATA)
			save_to_file("fibonacci_chart.csv", ss.str());

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		std::stringstream ss;
		solution Xopt;

		double ai = a;
		double bi = b;
		double ci = (a + b) / 2;
		double di{};

		int i = 0;
		double l{}, m{};
		solution ai_sol, bi_sol, ci_sol, di_sol;
		double l_prev{}, m_prev{}, di_prev{};
		do
		{
			ai_sol.x = ai;
			ai_sol.fit_fun(ff, ud1);

			bi_sol.x = bi;
			bi_sol.fit_fun(ff, ud1);

			ci_sol.x = ci;
			ci_sol.fit_fun(ff, ud1);

			l = m2d(ai_sol.y) * (pow(bi, 2) - pow(ci, 2)) + m2d(bi_sol.y) * (pow(ci, 2) - pow(ai, 2)) + m2d(ci_sol.y) * (pow(ai, 2) - pow(bi, 2));
			m = m2d(ai_sol.y) * (bi - ci) + m2d(bi_sol.y) * (ci - ai) + m2d(ci_sol.y) * (ai - bi);

			if (m <= 0)
			{
				Xopt.flag = 0;
				break;
			}

			di = 0.5 * l / m;
			di_sol.x = di;
			di_sol.fit_fun(ff, ud1);

			if (ai < di && di < ci)
			{
				if (di_sol.y < ci_sol.y)
				{
					bi = ci;
					ci = di;
				}
				else
					ai = di;
			}
			else
			{
				if (ci < di && di < bi)
				{
					if (di_sol.y < ci_sol.y)
					{
						ai = ci;
						ci = di;
					}
					else
						bi = di;
				}
				else
				{
					Xopt.flag = 0;
					break;
				}
			}

			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

			if (i > 0)
			{
				di_prev = 0.5 * l_prev / m_prev;
			}

			l_prev = l;
			m_prev = m;

			if (SAVE_CHART_DATA)
				ss << i << ";" << bi - ai << ";\n";
      
			++i;
		} while (!(bi - ai < epsilon || abs(di - di_prev) < gamma));

		Xopt.x = di;
		Xopt.fit_fun(ff, ud1);

		if (SAVE_CHART_DATA)
			save_to_file("lagrange_chart.csv", ss.str());

		return Xopt;
	}
	catch (string& ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//String builder do zapisu danych
		std::stringstream ss{};

		solution XT;

		//Wyliczanie punktu bazowego
		solution XB;
		XB.x = x0;
		XB.fit_fun(ff,ud1, ud2);

		while (s > epsilon)
		{
			//Poszukiwanie lepszego punktu
			XT = HJ_trial(ff, XB, s, ud1, ud2);

			//Znalezienie lepszego punktu
			if (XT.y < XB.y)
			{
				//Szukanie minimum przy sta³ym kroku do momentu roœniêcia funkcji
				do
				{
					//Zapis danych do wykresu
					if (SAVE_CHART_DATA)
						ss << m2d(XB.x(0)) << ";" << m2d(XB.x(1)) << "\n";

					solution XB_old = XB;
					XB = XT;
					XT.x = 2.0 * XB.x - XB_old.x;
					XT.fit_fun(ff, ud1, ud2);
					XT = HJ_trial(ff, XT, s, ud1, ud2);
					if (solution::f_calls > Nmax)
					{
						XT.flag = 0;
						throw std::string("Maximum amount of f_calls reached!");
					}
				} while (XT.y < XB.y);

				//Zast¹pienie wyiczonego punktu poprzednim (XT to by³ punkt w którym funkcja ju¿ ros³a)
				XT = XB;
			}
			//Brak lepszego punktu
			else
				//Zmniejszenie kroku poszukiwañ
				s = s * alpha;

			if (solution::f_calls > Nmax)
			{
				XT.flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}
		}

		//Zapis danych do pliku
		if (SAVE_CHART_DATA)
			save_to_file("HJ_chart.csv", ss.str());

		return XT;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_len(XB.x);

		//Baza szukania minimum
		matrix d = matrix(n, n);
		for (int i = 0; i < n; ++i)
			d(i, i) = 1.0;

		//Poszukiwanie lepszego minimum po ka¿dym kierunku funkcji
		for (int j = 0; j < n; ++j)
		{
			//Wyliczanie nowego punktu
			solution X_trial;
			X_trial.x = XB.x + s * d[j];
			X_trial.fit_fun(ff, ud1, ud2);

			//Znalezienie lepszego punktu w kierunku dodatnim
			if (X_trial.y < XB.y)
				XB = X_trial;
			//Sprawdzenie kierunku ujemnego
			else
			{
				X_trial.x = XB.x - s * d[j];
				X_trial.fit_fun(ff, ud1, ud2);
				//Znalezienie lepszego punktu w kierunku ujemnym
				if (X_trial.y < XB.y)
					XB = X_trial;
			}
		}
		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		std::stringstream ss;

		solution Xopt;

		//Funkcja pomocnicza szukaj¹ca maksymalnej wartoœci bezwzglêdnej w wektorze pionowym
		auto max = [&](matrix m) -> double
		{
			int len = get_len(m);
			double result = 0.0;
			for (int i = 0; i < len; ++i)
				if (abs(result) < abs(m(i)))
					result = abs(m(i));

			return result;
		};

		int i = 0;
		//Iloœæ zmiennych funkcji
		int n = get_len(x0);

		//Baza szukania kierunku
		matrix d = matrix(n, n);
		for (int j = 0; j < n; ++j)
			d(j, j) = 1.0;

		//Wektory pomocnicze
		matrix lambda = matrix(n, new double[n] {0.0});
		matrix p = matrix(n, new double[n] {0.0});

		//Wektor pionowy kroku
		matrix s = s0;

		//Wyliczanie punktu bazowego
		solution XB;
		XB.x = x0;
		XB.fit_fun(ff, ud1, ud2);

		if (SAVE_CHART_DATA)
			ss << m2d(XB.x(0)) << ";" << m2d(XB.x(1)) << "\n";

		while (max(s) >= epsilon)
		{
			//Poszukiwanie lepszego minimum po ka¿dym kierunku funkcji
			for (int j = 0; j < n; ++j)
			{
				//Wyliczanie nowego punktu
				solution XT;
				XT.x = XB.x + s(j) * d[j];
				XT.fit_fun(ff, ud1, ud2);

				//Znalezienie lepszego punktu
				if (XT.y < XB.y)
				{
					XB = XT;
					lambda(j) += s(j);
					s(j) *= alpha; //Ekspansja
				}
				//Brak lepszego punktu
				else
				{
					s(j) *= -beta; //Kontrakcja
					p(j) += 1;
				}
			}
			++i;

			//Zapis danych do wykresu
			if (SAVE_CHART_DATA)
				ss << m2d(XB.x(0)) << ";" << m2d(XB.x(1)) << "\n";

			//Znaleziony punkt przypisywany do wynikowego
			Xopt = XB;

			//Sprawdzanie czy wszystkie mo¿liwe kierunki zosta³y sprawdzone
			bool changeDirectionBase = true;
			for (int j = 0; j < n; ++j)
			{
				if (lambda(j) == 0 || p(j) == 0)
				{
					changeDirectionBase = false;
					break;
				}
			}

			//Zmiana bazy kierunku
			if (changeDirectionBase)
			{
				//Macierz trójk¹tna dolna lambdy
				matrix lambdaMatrix(n, n);
				int l = 0;
				for (int k = 0; k < n; ++k)
				{
					for (int j = 0; j <= k; ++j)
					{
						lambdaMatrix(k, j) = lambda(l);
					}
					++l;
				}

				//Wyliczanie macierzy Q
				matrix Q = d * lambdaMatrix;

				//Wyznaczanie sk³adowych wektora V
				matrix V = matrix(n, 1);
				V = Q[0];
				d[0] = V / norm(V);
				for (int j = 1; j < n; ++j)
				{
					matrix sum = matrix(n, new double[n] {0.0});
					for (int k = 0; k < j; ++k)
					{
						sum = sum + (trans(Q[j]) * d[k]) * d[k];
					}
					V = Q[j] - sum;
					d[j] = V / norm(V);
				}

				//Zerowanie wektorów pomocniczych
				lambda = matrix(n, new double[n] {0.0});
				p = matrix(n, new double[n] {0.0});

				//Zmiana d³ugoœci kroku na pocz¹tkowy
				s = s0;
			}

			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}
		}

		//Zapis danych do pliku
		if (SAVE_CHART_DATA)
			save_to_file("Rosen_chart.csv", ss.str());

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution XB;
		XB.x = x0;
		XB.fit_fun(ff, ud1, c);

		solution XT;
		XT = XB;

		double s = 0.5; //D³ugoœæ boku trójk¹ta
		double alpha = 1.0; //Wspó³czynnik odbicia
		double beta = 0.5; //Wspó³czynnik zwê¿enia
		double gamma = 2.0; //Wspó³czynnik ekspansji
		double delta = 0.5; //Wspó³czynnik redukcji

		do
		{
			XT = sym_NM(ff, XB.x, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1, c);
			c *= dc;

			if (solution::f_calls > Nmax)
			{
				XT.flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

			if (norm(XT.x - XB.x) < epsilon)
				break;

			XB = XT;
		} while (true);

		return XT;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		//Funkcja pomocnicza do znajdywania maksymum normy
		auto max = [&](std::vector<solution> sim, int i_min) -> double
		{
			double result = 0.0;
			for (int i = 0; i < sim.size(); ++i)
			{
				double normal = norm(sim[i_min].x - sim[i].x);
				if (result < normal)
					result = normal;
			}
			return result;
		};
		
		int n = get_len(x0);

		//Tworzenie bazy ortogonalnej
		matrix d = matrix(n, n);
		for (int i = 0; i < n; ++i)
			d(i, i) = 1.0;

		//Tworzenie simplexu i uzupe³nianie go danymi
		std::vector<solution> simplex;
		simplex.resize(n + 1);
		simplex[0].x = x0;
		simplex[0].fit_fun(ff, ud1, ud2);
		for (int i = 1; i < simplex.size(); ++i)
		{
			simplex[i].x = simplex[0].x + s * d[i - 1];
			simplex[i].fit_fun(ff, ud1, ud2);
		}

		//Indeks najmniejszej wartoœci wierzcho³ka simplexu
		int i_min{};
		//Indeks najwiêkszej wartoœci wierzcho³ka simplexu
		int i_max{};

		while (max(simplex, i_min) >= epsilon)
		{
			//Wyznaczanie maksymalnego i minimalnego indeksu
			i_min = 0;
			i_max = 0;
			for (int i = 1; i < simplex.size(); ++i)
			{
				if (simplex[i].y < simplex[i_min].y)
					i_min = i;
				if (simplex[i].y > simplex[i_max].y)
					i_max = i;
			}

			//Wyznaczenie œrodka ciê¿koœci
			matrix simplex_CoG{};
			for (int i = 0; i < simplex.size(); ++i)
			{
				if (i == i_max)
					continue;
				simplex_CoG = simplex_CoG + simplex[i].x;
			}
			simplex_CoG = simplex_CoG / simplex.size();

			//Obliczanie wartoœci funkcji odbitego simplexu
			solution simplex_reflected{};
			simplex_reflected.x = simplex_CoG + alpha * (simplex_CoG - simplex[i_max].x);
			simplex_reflected.fit_fun(ff, ud1, ud2);

			if (simplex_reflected.y < simplex[i_min].y)
			{
				//Obliczanie wartoœci funkcji powiêkszonego simplexu
				solution simplex_expansion{};
				simplex_expansion.x = simplex_CoG + gamma * (simplex_reflected.x - simplex_CoG);
				simplex_expansion.fit_fun(ff, ud1, ud2);
				if (simplex_expansion.y < simplex_reflected.y)
					simplex[i_max] = simplex_expansion;
				else
					simplex[i_max] = simplex_reflected;
			}
			else
			{
				if (simplex[i_min].y <= simplex_reflected.y && simplex_reflected.y < simplex[i_max].y)
					simplex[i_max] = simplex_reflected;
				else
				{
					//Obliczanie wartoœci funkcji pomniejszonego simplexu
					solution simplex_narrowed{};
					simplex_narrowed.x = simplex_CoG + beta * (simplex[i_max].x - simplex_CoG);
					simplex_narrowed.fit_fun(ff, ud1, ud2);
					if (simplex_narrowed.y >= simplex[i_max].y)
					{
						for (int i = 0; i < simplex.size(); ++i)
						{
							if (i == i_min)
								continue;
							simplex[i].x = delta * (simplex[i].x + simplex[i_min].x);
							simplex[i].fit_fun(ff, ud1, ud2);
						}
					}
					else
						simplex[i_max] = simplex_narrowed;
				}
			}

			if (solution::f_calls > Nmax)
			{
				simplex[i_min].flag = 0;
				throw std::string("Maximum amount of f_calls reached!");
			}

		}

		return simplex[i_min];
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		std::stringstream ss{};

		solution XB;
		XB.x = x0;

		solution XT;
		XT = XB;

		matrix d;

		while (true)
		{
			if (SAVE_CHART_DATA)
				ss << XB.x(0) << ";" << XB.x(1) << "\n";

			//Wyliczanie kierunku spadku funkcji
			XB.grad(gf, ud1, ud2);
			d = -XB.g;

			//Metoda zmiennokrokowa
			if (h0 <= 0)
			{
				matrix h_fun_data(2, 2);
				h_fun_data.set_col(XB.x, 0);
				h_fun_data.set_col(d, 1);
				solution h_sol = golden(ff, 0, 1, epsilon, Nmax, ud1, h_fun_data);
				matrix h = h_sol.x;
				XT.x = XB.x + h * d;
			}
			//Metoda sta³okrokowa
			else
			{
				XT.x = XB.x + h0 * d;
			}
			

			if (solution::g_calls > Nmax)
			{
				if (SAVE_CHART_DATA)
					save_to_file("SD_chart_h_" + std::to_string(h0) + ".csv", ss.str());
				XT.fit_fun(ff, ud1, ud2);
				return XT;
			}
				

			if (norm(XT.x - XB.x) <= epsilon)
				break;

			XB = XT;
		}

		if (SAVE_CHART_DATA)
			save_to_file("SD_chart_h_" + std::to_string(h0) + ".csv", ss.str());

		XT.fit_fun(ff, ud1, ud2);
		return XT;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		std::stringstream ss{};

		solution XB;
		XB.x = x0;

		solution XT;
		matrix g_prev;
		matrix g_curr;
		matrix d;

		XB.grad(gf, ud1, ud2);
		g_prev = XB.g;
		d = -g_prev;

		while (true)
		{
			ss << XB.x(0) << ";" << XB.x(1) << "\n";
			//Metoda zmiennokrokowa
			if (h0 <= 0)
			{
				matrix h_fun_data(2, 2);
				h_fun_data.set_col(XB.x, 0);
				h_fun_data.set_col(d, 1);
				solution h_sol = golden(ff, 0, 1, epsilon, Nmax, ud1, h_fun_data);
				matrix h = h_sol.x;
				XT.x = XB.x + h * d;
			}
			//Metoda sta³okrokowa
			else
			{
				XT.x = XB.x + h0 * d;
			}

			if (solution::g_calls > Nmax)
			{
				if (SAVE_CHART_DATA)
					save_to_file("CG_chart_h_" + std::to_string(h0) + ".csv", ss.str());

				XT.fit_fun(ff, ud1, ud2);
				return XT;
			}

			if (norm(XT.x - XB.x) <= epsilon)
				break;

			XT.grad(gf, ud1, ud2);
			g_curr = XT.g;

			double beta = pow(norm(g_curr), 2) / pow(norm(g_prev), 2);
			d = -g_curr + beta * d;

			g_prev = g_curr;
			XB = XT;
		}

		if (SAVE_CHART_DATA)
			save_to_file("CG_chart_h_" + std::to_string(h0) + ".csv", ss.str());

		XT.fit_fun(ff, ud1, ud2);
		return XT;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		std::stringstream ss{};

		solution XB;
		XB.x = x0;

		solution XT;

		matrix d;

		while (true)
		{
			if (SAVE_CHART_DATA)
				ss << XB.x(0) << ";" << XB.x(1) << "\n";

			XB.grad(gf, ud1, ud2);
			XB.hess(Hf, ud1, ud2);

			d = -inv(XB.H) * XB.g;

			//Metoda zmiennokrokowa
			if (h0 <= 0)
			{
				matrix h_fun_data(2, 2);
				h_fun_data.set_col(XB.x, 0);
				h_fun_data.set_col(d, 1);
				solution h_sol = golden(ff, 0, 1, epsilon, Nmax, ud1, h_fun_data);
				matrix h = h_sol.x;
				XT.x = XB.x + h * d;
			}
			//Metoda sta³okrokowa
			else
			{
				XT.x = XB.x + h0 * d;
			}

			if (solution::g_calls > Nmax)
			{
				if (SAVE_CHART_DATA)
					save_to_file("Newton_chart_h_" + std::to_string(h0) + ".csv", ss.str());

				XT.fit_fun(ff, ud1, ud2);
				return XT;
			}

			if (norm(XT.x - XB.x) <= epsilon)
				break;

			XB = XT;
		}

		if (SAVE_CHART_DATA)
			save_to_file("Newton_chart_h_" + std::to_string(h0) + ".csv", ss.str());

		XT.fit_fun(ff, ud1, ud2);
		return XT;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double alpha = (pow(5, 0.5) - 1) / 2;
		double a0 = a;
		double b0 = b;
		double c0 = b0 - alpha * (b0 - a0);
		double d0 = a0 + alpha * (b0 - a0);

		do
		{
			solution c0_sol;
			c0_sol.x = c0;
			c0_sol.fit_fun(ff, ud1, ud2);

			solution d0_sol;
			d0_sol.x = d0;
			d0_sol.fit_fun(ff, ud1, ud2);

			if (c0_sol.y < d0_sol.y)
			{
				b0 = d0;
				d0 = c0;
				c0 = b0 - alpha * (b0 - a0);
			}
			else
			{
				a0 = c0;
				c0 = d0;
				d0 = a0 + alpha * (b0 - a0);
			}

			if (solution::f_calls > Nmax)
				throw std::string("Maximum amount of f_calls reached!: ");

		} while (b0 - a0 > epsilon);

		Xopt.x = (a0 + b0) / 2;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

static double* expansion_Powell(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		double* p = new double[2] { 0, 0 };
		double x1 = x0 + d;
		int fcalls = 1;
		matrix x0Mat(x0);
		matrix x1Mat(x1);

		if (ff(x1Mat, ud1, ud2) == ff(x0Mat, ud1, ud2))
		{
			p[0] = m2d(x0Mat);
			p[1] = m2d(x1Mat);
			return p;
		}

		x1Mat = x1;
		if (ff(x1Mat, ud1, ud2) > ff(x0Mat, ud1, ud2))
		{
			d = -d;
			x1 = x0 + d;
			fcalls++;
			x1Mat = x1;
			if (ff(x1Mat, ud1, ud2) >= ff(x0Mat, ud1, ud2))
			{
				p[0] = m2d(x1Mat);
				p[1] = m2d(x0Mat) - d;
				return p;
			}
		}

		int i = 1;
		double xi_next = 0;
		matrix xtempMat(0);

		while (fcalls <= Nmax)
		{
			xi_next = x0 + alpha * i * d;
			fcalls++;
			x1Mat = xi_next;
			xtempMat = (x0 + alpha * (i - 1) * d);
			if (ff(x1Mat, ud1, ud2) >= ff(xtempMat, ud1, ud2))
				break;

			i++;
		}

		if (fcalls > Nmax)
		{
			throw string("Przekroczono liczbe wywolan\n");
		}
		else
		{
			if (d > 0)
			{
				p[0] = m2d(x0Mat) + alpha * (i - 2) * d;
				p[1] = m2d(x1Mat);
				return p;
			}
			else
			{
				p[0] = m2d(x1Mat);
				p[1] = m2d(x0Mat) + alpha * (i - 2) * d;
				return p;
			}
		}
	}
	catch (string ex_info)
	{
		throw("double* expansion(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;

		int n = get_len(x0);
		matrix d = ident_mat(n);
		matrix p;

		matrix h_sol_data(n, 2);
		solution h;

		matrix XB;
		XB = x0;

		double* range;

		while (true)
		{
			p = XB;

			for (int j = 0; j < n; ++j)
			{
				h_sol_data.set_col(p, 0);
				h_sol_data.set_col(d[j], 1);
				range = expansion(ff, 0.0, 1.0, 1.2, Nmax, ud1, h_sol_data);
				h = golden(ff, range[0], range[1], epsilon, Nmax, ud1, h_sol_data);
				p = p + h.x * d[j];
			}

			if (norm(p - XB) < epsilon)
			{
				Xopt.x = p;
				Xopt.fit_fun(ff, ud1, ud2);

				return Xopt;
			}

			if (solution::f_calls > Nmax)
				throw std::string("Maximum amount of f_calls reached!");

			for (int j = 0; j < n - 1; ++j)
				d.set_col(d[j + 1], j);
			d.set_col(p - XB, n - 1);

			h_sol_data.set_col(p, 0);
			h_sol_data.set_col(d[n - 1], 1);
			range = expansion(ff, 0.0, 1.0, 1.2, Nmax, ud1, h_sol_data);
			h = golden(ff, range[0], range[1], epsilon, Nmax, ud1, h_sol_data);

			XB = p + h.x * d[n - 1];
		}

		free(range);
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution* population = new solution[mi + lambda]; // Populacja (rodzice + potomkowie)
		solution* best_mu = new solution[mi]; // Najlepsi osobnicy

		matrix inv_fitness(mi, 1), temp_individual(N, 2); 

		double r, cumulative_sum, total_inv_fitness;
		double tau = 1.0 / sqrt(2 * N), tau_prime = 1.0 / sqrt(2 * sqrt(N)); // Wspó³czynniki mutacji

		int worst_idx; // Indeks najgorszego osobnika

		// Inicjalizacja populacji pocz¹tkowej
		for (int i = 0; i < mi; ++i)
		{
			population[i].x = matrix(N, 2);

			for (int j = 0; j < N; ++j)
			{
				population[i].x(j, 0) = (ub(j) - lb(j)) * m2d(rand_mat()) + lb(j);
				population[i].x(j, 1) = sigma0(0);
			}

			population[i].fit_fun(ff, ud1, ud2);

			if (population[i].y < epsilon)
			{
				population[i].flag = 1;
				return population[i];
			}
		}

		// G³ówna pêtla ewolucji
		while (true)
		{
			total_inv_fitness = 0;

			for (int i = 0; i < mi; ++i)
			{
				inv_fitness(i) = 1 / population[i].y(0);
				total_inv_fitness += inv_fitness(i);
			}

			// Selekcja rodziców
			for (int i = 0; i < lambda; ++i)
			{
				r = total_inv_fitness * m2d(rand_mat());
				cumulative_sum = 0;
				for (int j = 0; j < mi; ++j)
				{
					cumulative_sum += inv_fitness(j);
					if (r <= cumulative_sum)
					{
						population[mi + i] = population[j];
						break;
					}
				}
			}

			// Mutacja potomków
			for (int i = 0; i < lambda; ++i)
			{
				r = m2d(randn_mat());
				for (int j = 0; j < N; ++j)
				{
					population[mi + i].x(j, 1) *= exp(tau_prime * r + tau * m2d(randn_mat()));
					population[mi + i].x(j, 0) += population[mi + i].x(j, 1) * m2d(randn_mat());
				}
			}

			// Krzy¿owanie
			for (int i = 0; i < lambda; i += 2)
			{
				r = m2d(rand_mat());
				temp_individual = population[mi + i].x;
				population[mi + i].x = r * population[mi + i].x + (1 - r) * population[mi + i + 1].x;
				population[mi + i + 1].x = r * population[mi + i + 1].x + (1 - r) * temp_individual;
			}

			// Ocena funkcji celu dla potomków
			for (int i = 0; i < lambda; ++i)
			{
				population[mi + i].fit_fun(ff, ud1, ud2);
				if (population[mi + i].y < epsilon)
				{
					population[mi + i].flag = 1;
					return population[mi + i];
				}
			}

			// Selekcja mi najlepszych osobników
			for (int i = 0; i < mi; ++i)
			{
				worst_idx = 0;
				for (int j = 1; j < mi + lambda; ++j)
				{
					if (population[worst_idx].y > population[j].y)
						worst_idx = j;
				}
				best_mu[i] = population[worst_idx];
				population[worst_idx].y = 1e10;
			}

			// Aktualizacja populacji bazowej
			for (int i = 0; i < mi; ++i)
				population[i] = best_mu[i];

			if (solution::f_calls > Nmax)
				break;
		}

		population[0].flag = 0;
		return population[0];
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
