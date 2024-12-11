#define _USE_MATH_DEFINES
#include <cmath>

#include"user_funs.h"

//LAB 0
matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m * pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
	return dY;
}

//LAB 1
matrix ff1T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * m2d(x)) * exp(-1.0 * pow((0.1 * m2d(x) - 2 * M_PI), 2)) + 0.002 * pow(0.1 * m2d(x), 2);
	return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)
{
	//Warunki pocz�tkowe
	matrix Y0 = matrix(3, new double[3] {
		5.0, //Va_0
		1.0,  //Vb_0
		20.0 //Tb_0
	});

	//Symulacja dla podanych warunk�w pocz�tkowych, danych zadania i szukanego Da
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, x);

	//Szykanie maksymalnej temperatury
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 0; i < n; ++i)
	{
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	}

	//Odchylenie maksymalnej temperatury od 50 stopni Celsujsza
	matrix y;
	y = abs(max - 50.0);

	Y[0].~matrix();
	Y[1].~matrix();

	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	//Wektor zmian po czasie
	matrix dY(3, 1);

	//Zmienne zadania
	double Va = Y(0);
	double Vb = Y(1);
	double Tb = Y(2);

	//Dane z zadania
	double Pa = ud1(0);
	double Ta = ud1(1);
	double Pb = ud1(2);
	double Db = ud1(3);
	double F_in = ud1(4);
	double T_in = ud1(5);
	double a = ud1(6);
	double b = ud1(7);
	double g = ud1(8);

	//Pole otworu w zbiorniku A
	double Da = ud2(0);

	//Wyliczanie wylanego strumienia ze zbiornika A
	double Fa_out{};
	if (Va > 0.0)
		Fa_out = a * b * Da * sqrt(2 * g * Va / Pa);

	//Wyliczanie wylanego strumienia ze zbiornika B
	double Fb_out{};
	if (Vb > 0.0)
		Fb_out = a * b * Db * sqrt(2 * g * Vb / Pb);


	//Ustalanie zmiany obj�to�ci w zbiornku A
	if (Y(0) + dY(0) < 0)
		dY(0) = -Y(0); //Wylanie reszty je�li strumie� wi�kszy od obj�to�ci wody
	else
		dY(0) = -Fa_out; //Wylanie strumienia

	//Ustalanie zmien obj�to�ci w zbiorniku B
	if (Y(1) + dY(1) < 0)
		dY(1) = -Y(1); //Wylanie reszty je�li strumie� wi�kszy od obj�to�ci wody
	else
		dY(1) = Fa_out + F_in - Fb_out; //Wylanie strumienia oraz wlanie wody z kranu i ze zbiornika A

	//Ustalenie zmian temperatury w zbiorniku B
	if (Vb > 0)
		dY(2) = (F_in / Vb) * (T_in - Tb) + (Fa_out / Vb) * (Ta - Tb); //Formu�a je�li zbiornik B nie jest pusty
	else
		dY(2) = 0; //Pusty zbiornik B

	//Zwracanie zmian po czasie
	return dY;
}

//LAB 2

matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
}

matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = 0;

	//Warunki początkowe
	matrix Y0 = matrix(2, 1);

	//Symulacja
	matrix* Y = solve_ode(df2, 0.0, 0.1, 100, Y0, ud1, x);

	//Dane referencyjne
	double alpha_ref = ud1(4);
	double omega_ref = ud1(5);

	//Obliczanie funkjconału jakości metodą prostokątów
	int n = get_len(Y[0]);
	for (int i = 0; i < n; ++i)
	{
		y = y + 10 * pow(alpha_ref - Y[1](i, 0), 2) + pow(omega_ref - Y[1](i, 1), 2) + pow(x(0) * (alpha_ref - Y[1](i, 0)) + x(1) * (omega_ref - Y[1](i, 1)), 2);
	}
	y = 0.1 * y;

	return y;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{
	//Wektor zmian po czasie
	matrix dY(2, 1);

	//Zmienne zadania
	double alpha = Y(0);
	double omega = Y(1);

	//Dane z zadania
	double l = ud1(0);
	double m_r = ud1(1);
	double m_c = ud1(2);
	double b = ud1(3);
	double alpha_ref = ud1(4);
	double omega_ref = ud1(5);

	//Współczynniki wzmocnienia
	double k1 = ud2(0);
	double k2 = ud2(1);

	//Moment bezwładności
	double I = (1.0 / 3.0) * m_r * pow(l, 2) + m_c * pow(l, 2);

	//Moment siły
	double Mt = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);

	dY(0) = Y(1);
	dY(1) = (Mt - b * omega) / I;

	return dY;
}

//LAB 3

matrix ff3T_outside(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = (sin(M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)));

	double a = m2d(ud1);
	double c = m2d(ud2);

	//g1
	if (-x(0) + 1 > 0)
		y = y + c * pow(-x(0) + 1, 2);
	//g2
	if (-x(1) + 1 > 0)
		y = y + c * pow(-x(1) + 1, 2);
	//g3
	if (norm(x) - a > 0)
		y = y + c * pow(norm(x) - a, 2);

	return y;
}

matrix ff3T_inside(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = (sin(M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)))) / (M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2)));

	double a = m2d(ud1);
	double c = m2d(ud2);

	//g1
	if (-x(0) + 1 >= 0)
		y = 1E10;
	else
		y = y - c / (-x(0) + 1);
	//g2
	if (-x(1) + 1 >= 0)
		y = 1E10;
	else
		y = y - c / (-x(1) + 1);
	//g3
	if (norm(x) - a > 0)
		y = 1E10;
	else
		y = y - c / (norm(x) - a);

	return y;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0(4, new double[4] {0.0, x(0), 100, 0});
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, x(1));
	int n = get_len(Y[0]);
	int i50 = 0;
	int i_0 = 0;
	for (int i = 0; i < n; ++i)
	{
		if (abs(Y[1](i, 2) - 50.0) < abs(Y[1](i50, 2) - 50.0))
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i_0, 2)))
			i_0 = i;
	}

	y = -Y[1](i_0, 0);

	if (abs(x(0)) - 10 > 0)
		y = y + ud2 * pow(abs(x(0)) - 10, 2);
	if (abs(x(1)) - 15 > 0)
		y = y + ud2 * pow(abs(x(1)) - 15, 2);
	if (abs(Y[1](i50, 0) - 5.0) - 0.5 > 0)
		y = y + ud2 * pow(abs(Y[1](i50, 0) - 5.0) - 0.5, 2);

	return y;
}

matrix df3(double t, matrix Y, matrix ud1, matrix ud2)
{
	//Wektor zmian po czasie
	matrix dY(4, 1);

	//Zmienne dane
	double x = Y(0);
	double v_x = Y(1);
	double y = Y(2);
	double v_y = Y(3);

	//Dane zadania
	double C = ud1(0);
	double rho = ud1(1);
	double r = ud1(2);
	double m = ud1(3);
	double g = ud1(4);

	double s = M_PI * pow(r, 2);

	double D_x = (1.0 / 2.0) * C * rho * s * v_x * abs(v_x);
	double D_y = (1.0 / 2.0) * C * rho * s * v_y * abs(v_y);
	double FM_x = rho * v_y * m2d(ud2) * M_PI * pow(r, 3);
	double FM_y = rho * v_x * m2d(ud2) * M_PI * pow(r, 3);

	dY(0) = v_x;
	dY(1) = (-D_x - FM_x) / m;
	dY(2) = v_y;
	dY(3) = ((- m * g) - D_y - FM_y) / m;

	return dY;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;

	if (isnan(ud2(0, 0)))
		y = pow((x(0) + 2 * x(1) - 7), 2) + pow((2 * x(0) + x(1) - 5), 2);
	else
		y = ff4T(ud2[0] + x * ud2[1]);

	return y;
}

matrix gf4T(matrix x, matrix ud1, matrix ud2)
{
	matrix y(2, 1);

	y(0) = -34.0 + 10.0 * x(0) + 8.0 * x(1);
	y(1) = -38.0 + 8.0 * x(0) + 10.0 * x(1);

	return y;
}

double sigmoid(matrix theta, matrix x)
{
	return 1.0 / (1.0 + exp(-1.0 * m2d(trans(theta) * x)));
}

matrix hf4T(matrix x, matrix ud1, matrix ud2)
{
	matrix y(2, 2);

	y(0, 0) = 10;
	y(0, 1) = 8;
	y(1, 0) = 8;
	y(1, 1) = 10;

	return y;
}

matrix ff4R(matrix theta, matrix X, matrix Y)
{
	matrix y;

	double sum = 0.0;
	for (int i = 0; i < 100; ++i)
	{
		double y_i = Y[i](0);
		matrix x_i = X[i];
		sum += y_i * log(sigmoid(theta, x_i)) + (1.0 - y_i) * log(1.0 - sigmoid(theta, x_i));
	}

	y = (-1.0 / 100.0) * sum;

	return y;
}

matrix gf4R(matrix theta, matrix X, matrix Y)
{
	matrix y(3, 1);

	for (int j = 0; j < 3; ++j)
	{
		double sum = 0.0;
		for (int i = 0; i < 100; ++i)
		{
			double y_i = Y[i](0);
			matrix x_i = X[i];

			sum += (sigmoid(theta, x_i) - y_i) * x_i(j);
		}
		y(j) = (1.0 / 100.0) * sum;
	}

	return y;
}


