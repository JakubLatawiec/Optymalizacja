#define _USE_MATH_DEFINES
#include <cmath>
#include"user_funs.h"

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
	matrix result(1, 1);  
	result(0, 0) = x(0, 0) * x(0, 0) + x(1, 0) * x(1, 0) - cos(2.5 * M_PI * x(0, 0)) - cos(2.5 * M_PI * x(1, 0)) + 2;
	return result;
}

matrix ff1(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = -cos(0.1 * m2d(x)) * exp(-1.0 * pow((0.1 * m2d(x) - 2 * M_PI), 2)) + 0.002 * pow(0.1 * m2d(x), 2);
	return y;
}

matrix ff0T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
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
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}


//lab_2

const double l = 1.0;       // D³ugoœæ ramienia (m)
const double m = 1.0;       // Masa ramienia (kg)
const double mc = 5.0;      // Masa ciê¿arka (kg)
const double b = 0.5;       // Wspó³czynnik tarcia (Nms)
const double alpha_target = M_PI;   // Docelowy k¹t (rad)
const double omega_target = 0.0;    // Docelowa prêdkoœæ k¹towa (rad/s)
const double I = m * l * l + mc * l * l; // Moment bezw³adnoœci


double calculate_torque(double k1, double k2, double alpha, double omega)
{

	return (k1 * (alpha_target - alpha) + k2 * (omega_target - omega));
}


double simulate(double k1, double k2) {
	double dt = 0.1;      // Krok czasowy (s)
	double t_max = 100.0;  // Czas symulacji (s)

	double alpha = 0.0;    // Pocz¹tkowy k¹t (rad)
	double omega = 0.0;    // Pocz¹tkowa prêdkoœæ k¹towa (rad/s)
	double Q = 0.0;  // Wartoœæ funkcjona³u jakoœci

	// Symulacja dynamiki
	for (double t = 0; t <= t_max; t += dt) {
		double M = calculate_torque(k1, k2, alpha, omega); // Moment si³y
		double alpha_ddot = (M - b * omega) / I;           // Przyspieszenie k¹towe

		// Aktualizacja prêdkoœci i pozycji (metoda Eulera)
		omega += alpha_ddot * dt;
		alpha += omega * dt;

		// Aktualizacja funkcjona³u jakoœci
		Q += dt * (10 * std::pow(alpha_target - alpha, 2) +
			std::pow(omega_target - omega, 2) +
			std::pow(M, 2));
	}

	return Q;
}

// Funkcja celu wywo³ywana przez algorytmy optymalizacji
matrix Q_function(matrix k, matrix, matrix) {
	double k1 = k(0, 0);
	double k2 = k(1, 0);
	double Q_value = simulate(k1, k2);
	matrix result(1, 1);
	result(0, 0) = Q_value;
	return result;
}


