/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include "environment.h"
#include"opt_alg.h"
#include "configuration.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		//lab0();
		//lab1();
		//lab2();
		//lab3();
		lab4();
	}
	catch (const string& EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

//void lab0()
//{
//
//#ifdef SAVE_TO_FILE
//	create_environment("lab00");
//#endif
//
//#ifdef CALC_TEST
//	//Funkcja testowa
//	double epsilon = 1e-2;
//	int Nmax = 10000;
//	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
//	solution opt;
//	a(0) = -1;
//	a(1) = 2;
//	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
//	cout << opt << endl << endl;
//	solution::clear_calls();
//#endif
//
//	
//#ifdef CALC_SIMULATION
//	//Wahadlo
//	Nmax = 1000;
//	epsilon = 1e-2;
//	lb = 0;
//	ub = 5;
//	double teta_opt = 1;
//	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
//	cout << opt << endl << endl;
//	solution::clear_calls();
//
//	//Zapis symulacji do pliku csv
//	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });
//	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
//
//#ifdef SAVE_TO_FILE
//	save_to_file("simulation.csv", hcat(Y[0], Y[1]));
//#endif
//
//	Y[0].~matrix();
//	Y[1].~matrix();
//#endif
//}

void lab1()
{
#ifdef SAVE_TO_FILE
	create_environment("lab01");
#endif

	//Dane dok³adnoci wyników
	double epsilon = 1e-18;
	double gamma = 1e-30;
	int Nmax = 200;
	double d = 0.01;
	double alpha = 1.1;

	//Funkcja testowa
#ifdef CALC_TEST

	//Generator losowania liczb
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(0.0, 100.0);

	//Stringstream do zapisu danych
	std::stringstream test_ss;

	//Solution dla testów
	solution test_opt;

	double test_alpha = 1.5;
	for (int j = 0; j < 3; ++j)
	{
		//Liczenie ekspansji, fibonacciego i lagrange'a dla danego wspó³czynnika alpha
		for (int i = 0; i < 100; ++i)
		{
			double x0 = x0_dist(gen);
			double* bounds = expansion(ff1T, x0, d, test_alpha, Nmax);
			test_ss << x0 << ";" << bounds[0] << ";" << bounds[1] << ";" << bounds[2] << ";";
			solution::clear_calls();

			test_opt = fib(ff1T, bounds[0], bounds[1], epsilon);
			test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";";
			solution::clear_calls();

			test_opt = lag(ff1T, bounds[0], bounds[1], epsilon, gamma, Nmax);
			test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";\n";
			solution::clear_calls();
		}

		//Zapis do pliku
#ifdef SAVE_TO_FILE
		save_to_file("test_alpha_" + std::to_string(test_alpha) + ".csv", test_ss.str());
#endif

		//Czyszczenie zawartoci ss
		test_ss.str(std::string());

		//Zmiana alfy
		test_alpha += 1.3;
	}

	//Czyszczenie zawartoci ss
	test_ss.str(std::string());

	//Obliczanie minimum metod¹ Fibonacci'ego
	SAVE_CHART_DATA = true;

	test_opt = fib(ff1T, -100.0, 100.0, epsilon);
	test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";";
	solution::clear_calls();

	//Obliczanie minimum metod¹ Lagrange'a
	test_opt = lag(ff1T, -100, 100, epsilon, gamma, Nmax);
	test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";\n";
	solution::clear_calls();

	SAVE_CHART_DATA = false;

	//Zapis do pliku
#ifdef SAVE_TO_FILE
	save_to_file("test_no_expansion.csv", test_ss.str());
#endif

#endif

#ifdef CALC_SIMULATION

	//Dane do problemu rzeczywistego
	matrix ud1 = matrix(9, 1);
	ud1(0) = 0.5; //Pole zbiornika A (Pa) [m^2]
	ud1(1) = 90.0; //Temperatura wody w zbiorniku A (Ta) [Celsjusz]
	ud1(2) = 1.0; //Pole zbiornika B (Pb) [m^2]
	ud1(3) = 36.5665 * 0.0001; //Otwór w zbiorniku B (Db) [m^2]
	ud1(4) = 10 * 0.001; //Strumieñ wody z kranu (F_in) [m^3/s]
	ud1(5) = 20.0; //Temperatura wody z kranu (T_in) [Celsjusz]
	ud1(6) = 0.98; //Wspó³czynnik lepkoci cieczy (a) [-]
	ud1(7) = 0.63; //Wspó³czynnik zwê¿ania strumienia cieczy (b) [-]
	ud1(8) = 9.81; //Przypieszenie ziemskie (g) [m/s^2]

	//Zakres szukania Da
	double Da_0_s = 1.0 * 0.0001;
	double Da_0_f = 100 * 0.0001;

	//Szukanie minimum
	solution opt = fib(ff1R, Da_0_s, Da_0_f, epsilon, ud1);
	std::cout << opt;
	solution::clear_calls();

	//Warunki pocz¹tkowe
	matrix Y0 = matrix(3, 1);
	Y0(0) = 5.0; //Poczatkowa objetosc w a
	Y0(1) = 1.0; //Poczatkowa objetosc w b
	Y0(2) = 20.0;//Poczatkowa temperatura w b

	//Symulacja
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, opt.x);

#ifdef SAVE_TO_FILE
	save_to_file("simulation_fibonacci.csv", hcat(Y[0], Y[1]));
#endif

	//Szukanie maksymalnej temperatury
	int n = get_len(Y[0]);
	double Tb_max = Y[1](0, 2);
	for (int i = 0; i < n; ++i)
	{
		if (Tb_max < Y[1](i, 2))
			Tb_max = Y[1](i, 2);
	}
	std::cout << "Tb_max (fibonacci): " << Tb_max << "\n\n";

	//Szukanie minimum
	opt = lag(ff1R, Da_0_s, Da_0_f, epsilon, epsilon, Nmax, ud1);
	std::cout << opt;
	solution::clear_calls();

	//Symulacja
	Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, opt.x);

#ifdef SAVE_TO_FILE
	save_to_file("simulation_lagrange.csv", hcat(Y[0], Y[1]));
#endif

	//Szukanie maksymalnej temperatury
	n = get_len(Y[0]);
	Tb_max = Y[1](0, 2);
	for (int i = 0; i < n; ++i)
	{
		if (Tb_max < Y[1](i, 2))
			Tb_max = Y[1](i, 2);
	}
	std::cout << "Tb_max (lagrange): " << Tb_max << "\n";

	Y[0].~matrix();
	Y[1].~matrix();
#endif
}

void lab2()
{

#ifdef SAVE_TO_FILE
	create_environment("lab02");
#endif

	//Dane dokładności wyników
	double s = 0.1;
	double alpha = 0.1;
	double beta = 0.1;
	double epsilon = 1E-6;
	double Nmax = 2000;
	double alphaRosen = 1.2;

#ifdef CALC_TEST

	//Generator liczb losowych
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(-1.0, 1.0);

	//Stringstream do zapisu danych
	std::stringstream test_ss;

	//Rozwiązanie dla wyników testowych
	solution test_sol;

	//Dane dokładności dla testów
	double test_s = 0.1;

	//Punty startowe dla testów
	matrix test_x0{};

	//Liczenie testów
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 100; ++j)
		{
			test_x0 = matrix(2, new double[2] {x0_dist(gen), x0_dist(gen)});
			test_sol = HJ(ff2T, test_x0, s, alpha, epsilon, Nmax);
			test_ss << test_x0(0) << ";" << test_x0(1) << ";" << m2d(test_sol.x(0)) << ";" << m2d(test_sol.x(1)) << ";" << m2d(test_sol.y) << ";" << test_sol.f_calls << ";" << (abs(m2d(test_sol.x(0))) < 0.01 && abs(m2d(test_sol.x(1))) < 0.01 ? "TAK" : "NIE") << ";";
			solution::clear_calls();
			test_sol = Rosen(ff2T, test_x0, matrix(2, new double[2] {s, s}), alphaRosen, beta, epsilon, Nmax);
			test_ss << m2d(test_sol.x(0)) << ";" << m2d(test_sol.x(1)) << ";" << m2d(test_sol.y) << ";" << test_sol.f_calls << ";" << (abs(m2d(test_sol.x(0))) < 0.01 && abs(m2d(test_sol.x(1))) < 0.01 ? "TAK" : "NIE") << "\n";
			solution::clear_calls();
		}

		//Zapis do pliku
#ifdef SAVE_TO_FILE
		save_to_file("test_s_" + std::to_string(test_s) + ".csv", test_ss.str());
#endif

		//Czyszczenie zawartoci ss
		test_ss.str(std::string());

		//Zmiana długości kroku
		test_s += 0.6;
	}


	//Liczenie danych do wykresów
	SAVE_CHART_DATA = true;

	//Wspólny punkt dla obu metod
	test_x0 = matrix(2, new double[2] {0.45, 0.45});

	HJ(ff2T, test_x0, s, alpha, epsilon, Nmax);
	solution::clear_calls();

	Rosen(ff2T, test_x0, matrix(2, new double[2] {s, s}), alphaRosen, beta, epsilon, Nmax);
	solution::clear_calls();

	SAVE_CHART_DATA = false;
#endif

#ifdef CALC_SIMULATION

	//Dane symulacji
	matrix ud1 = matrix(6, new double[6] {
		1.0, //Długość ramienia (l) [m]
		1.0, //Masa ramienia (m_r) [kg]
		5.0, //Masa ciężarka (m_c) [kg]
		0.5, //Współczynnik tarcia (b) [Nms]
		3.14, //Referencyjny kąt (alpha_ref) [rad]
		0.0 //Referencyjna prędkość kątowa (omega_ref) [rad/s] 
	});

	//Początkowe wartości współczynników
	matrix k_0 = matrix(2, new double[2] {1.0, 1.0});

	//Warunki początkowe rónwe 0.0
	matrix Y0 = matrix(2, 1);

	//Wyliczanie optymalnych wartości współczynników metodą HJ
	solution opt = HJ(ff2R, k_0, s, alpha, epsilon, Nmax, ud1);
	std::cout << opt << "\n";
	solution::clear_calls();
	matrix* Y = solve_ode(df2, 0.0, 0.1, 100.0, Y0, ud1, opt.x);
#ifdef SAVE_TO_FILE
	save_to_file("simulation_HJ.csv", hcat(Y[0], Y[1]));
#endif

	//Wyliczanie optymalnych wartości współczynników metodą Rosenbrocka
	opt = Rosen(ff2R, k_0, matrix(2, new double[2]{s, s}), alphaRosen, beta, epsilon, Nmax, ud1);
	std::cout << opt << "\n";
	solution::clear_calls();
	Y = solve_ode(df2, 0.0, 0.1, 100.0, Y0, ud1, opt.x);
#ifdef SAVE_TO_FILE
	save_to_file("simulation_Rosen.csv", hcat(Y[0], Y[1]));
#endif

#endif
}




std::string solution_to_string(const solution& sol) {
	std::ostringstream oss;

	oss << "x = ";
	for (int i = 0; i < get_len(sol.x); ++i) {
		oss << sol.x(i, 0) << (i < get_len(sol.x) - 1 ? ", " : "");
	}
	oss << "\n";

	oss << "y = " << sol.y(0, 0) << "\n";
	oss << "Exit flag = " << sol.flag << "\n";

	return oss.str();
}




void lab3()
{
	double epsilon = 1e-3;
	double s = 0.5;
	double alpha = 1.0;
	double beta = 0.5;
	double gamma = 2.0;
	double delta = 0.5;
	int Nmax = 1000;

	std::vector<double> a_values = { 4.0, 4.4934, 5.0 };

	for (double a : a_values) {
		std::cout << "Optimizations for a = " << a << ":\n";

		matrix ud1 = matrix(1, 1, a); // Set parameter a for constraints
		for (int i = 1; i <= 100; ++i) {
			matrix x0 = matrix(2, new double[2] {(rand() % 100) / 100.0, (rand() % 100) / 100.0});

			solution result = sym_NM(ff3T, x0, s, alpha, beta, gamma, delta, epsilon, Nmax, ud1);

			double r = norm(result.x);
			std::cout << "Iteration " << i << ": x* = [" << result.x(0, 0) << ", " << result.x(1, 0) << "], f(x*) = " << result.y(0, 0) << ", r = " << r << "\n";
		}
	}
}

void lab4() {
#ifdef SAVE_TO_FILE
	create_environment("lab04");
#endif

	double epsilon = 1e-6;
	double h0 = 0.05;
	double a = 0.0;
	double b = 1.0;
	int Nmax = 1000;

	matrix x0 = matrix(2, new double[2] {0.0, 0.0});
	matrix ud1 = NAN;
	matrix ud2 = NAN;

	solution grad_result = CG(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2);
	std::cout << "Gradient Calculation Result:\n" << grad_result << "\n";

	solution golden_result = golden(ff4T, a, b, epsilon, Nmax, ud1, ud2);
	std::cout << "Golden Section Search Result:\n" << golden_result << "\n";

	solution sd_result = SD(ff4T, gf4T, x0, h0, epsilon, Nmax, ud1, ud2);
	std::cout << "Steepest Descent Result:\n" << sd_result << "\n";

#ifdef SAVE_TO_FILE
	save_to_file("lab4_gradient.csv", solution_to_string(grad_result));
	save_to_file("lab4_golden.csv", solution_to_string(golden_result));
	save_to_file("lab4_SD.csv", solution_to_string(sd_result));
#endif
}





void lab5()
{

}

void lab6()
{

}




