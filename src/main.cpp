/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
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
		lab1();
	}
	catch (const string& EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{

#ifdef SAVE_TO_FILE
	create_environment("lab00");
#endif

#ifdef CALC_TEST
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();
#endif

	
#ifdef CALC_SIMULATION
	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);

#ifdef SAVE_TO_FILE
	save_to_file("simulation.csv", hcat(Y[0], Y[1]));
#endif

	Y[0].~matrix();
	Y[1].~matrix();
#endif
}

void lab1()
{

#ifdef SAVE_TO_FILE
	create_environment("lab01");
#endif

	//Dane dok�adno�ci wynik�w
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

	//Solution dla test�w
	solution test_opt;

	double test_alpha = 1.5;
	for (int j = 0; j < 3; ++j)
	{
		//Liczenie ekspansji, fibonacciego i lagrange'a dla danego wsp�czynnika alpha
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

		//Czyszczenie zawarto�ci ss
		test_ss.str(std::string());

		//Zmiana alfy
		test_alpha += 1.3;
	}

	//Czyszczenie zawarto�ci ss
	test_ss.str(std::string());

	//Obliczanie minimum metod� Fibonacci'ego
	SAVE_CHART_DATA = true;

	test_opt = fib(ff1T, -100.0, 100.0, epsilon);
	test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";";
	solution::clear_calls();

	//Obliczanie minimum metod� Lagrange'a
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
	ud1(3) = 36.5665 * 0.0001; //Otw�r w zbiorniku B (Db) [m^2]
	ud1(4) = 10 * 0.001; //Strumie� wody z kranu (F_in) [m^3/s]
	ud1(5) = 20.0; //Temperatura wody z kranu (T_in) [Celsjusz]
	ud1(6) = 0.98; //Wsp�czynnik lepko�ci cieczy (a) [-]
	ud1(7) = 0.63; //Wsp�czynnik zw�ania strumienia cieczy (b) [-]
	ud1(8) = 9.81; //Przy�pieszenie ziemskie (g) [m/s^2]

	//Zakres szukania Da
	double Da_0_s = 1.0 * 0.0001;
	double Da_0_f = 100 * 0.0001;

	//Szukanie minimum
	solution opt = fib(ff1R, Da_0_s, Da_0_f, epsilon, ud1);
	std::cout << opt;
	solution::clear_calls();

	//Warunki pocz�tkowe
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

}

void lab3()
{

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
