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
		lab3();
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

void lab3()
{

#ifdef SAVE_TO_FILE
	create_environment("lab03");
#endif

	//Dane dokładnościowe
	double epsilon = 1E-3;
	int Nmax = 10000;
	double c_inside = 100;
	double dc_inside = 0.2;
	double c_outside = 1.0;
	double dc_outside = 1.5;

#ifdef CALC_TEST	
	//Generator liczb losowych
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(1.5, 5.5);

	//Stringstream do zapisu danych
	std::stringstream test_ss;

	//Rozwiązanie dla wyników testowych
	solution test_sol;

	//Dane a dla testów
	matrix a = matrix(4.0);

	//Punty startowe dla testów
	matrix test_x0{};
	
	for (int i = 0; i < 3; ++i)
	{
		if (i == 0)
			a = matrix(4.0);
		else if (i == 1)
			a = matrix(4.4934);
		else
			a = matrix(5.0);

		for (int j = 0; j < 100; ++j)
		{
			test_x0 = matrix(2, new double[2] {x0_dist(gen), x0_dist(gen)});
			test_ss << test_x0(0) << ";" << test_x0(1) << ";";
			test_sol = pen(ff3T_outside, test_x0, c_outside, dc_outside, epsilon, Nmax, a);
			test_ss << test_sol.x(0) << ";" << test_sol.x(1) << ";" << sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << ";" << test_sol.y << test_sol.f_calls << ";";
			solution::clear_calls();
			test_sol = pen(ff3T_inside, test_x0, c_inside, dc_inside, epsilon, Nmax, a);
			test_ss << test_sol.x(0) << ";" << test_sol.x(1) << ";" << sqrt(pow(test_sol.x(0), 2) + pow(test_sol.x(1), 2)) << ";" << test_sol.y << test_sol.f_calls << "\n";
			solution::clear_calls();
		}

#ifdef SAVE_TO_FILE
		save_to_file("test_a_" + std::to_string(m2d(a)) + ".csv", test_ss.str());
#endif

		//Czyszczenie zawartoci ss
		test_ss.str(std::string());
	}

#endif

#ifdef CALC_SIMULATION
	//Dane zadania
	matrix ud1 = matrix(5, new double[5] {
		0.47, //Współczynnik oporu (C) [-]
		1.2, //Gęstość powietrza (rho) [kg/m^3]
		0.12, //Promień piłki (r) [m]
		0.6, //Masa piłki (m) [kg]
		9.81 //Przyśpieszenie ziemskie (g) [m/s^2]
	});

	//Początkowe wartości szukania minimum
	matrix x0 = matrix(2, new double[2] {-5.0, 5.0});

	//Szukanie optymalnej prędkości początkowej po osi x i początkowej prędkości obrotowej
	solution opt = pen(ff3R, x0, c_outside, dc_outside, epsilon, Nmax, ud1);
	std::cout << opt << "\n";

	//Symulacja lotu piłki dla wyznaczonych ograniczeń
	matrix Y0(4, new double[4] {0.0, opt.x(0), 100, 0});
	matrix* Y = solve_ode(df3, 0.0, 0.01, 7.0, Y0, ud1, opt.x(1));

#ifdef SAVE_TO_FILE
	save_to_file("simulation.csv", hcat(Y[0], Y[1]));
#endif

#endif
}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
