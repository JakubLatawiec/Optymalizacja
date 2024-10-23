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

			if (ai_sol.f_calls > Nmax)
			{
				Xopt.flag = 0;
				throw std::string("Error message!");
				break;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
		//Tu wpisz kod funkcji

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
		solution Xopt;
		//Tu wpisz kod funkcji

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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
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
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

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
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution EA(...):\n" + ex_info);
	}
}
