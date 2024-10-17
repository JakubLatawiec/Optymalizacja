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
		lab0();
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

	//Dane dok³adnoœci wyników
	double epsilon = 1e-20;
	double gamma = 1e-20;
	int Nmax = 200;
	double d = 0.01;
	double alpha = 1.1;

	//Funkcja testowa
#ifdef CALC_TEST

	//Generator losowania liczb
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0_dist(0.0, 100.0);
	std::uniform_real_distribution<> alpha_dist(1.1, 2.0);

	//Stringstream do zapisu danych
	std::stringstream test_ss;

	//Solution dla testów
	solution test_opt;

	//
	for (int j = 0; j < 3; ++j)
	{
		double test_alpha = alpha_dist(gen);
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

		//Czyszczenie zawartoœci ss
		test_ss.str(std::string());
	}



	//Czyszczenie zawartoœci ss
	test_ss.str(std::string());

	//Obliczanie minimum metod¹ Fibonacci'ego
	test_opt = fib(ff1T, -100, 100, epsilon);
	test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";";
	solution::clear_calls();

	//Obliczanie minimum metod¹ Lagrange'a
	test_opt = lag(ff1T, -100, 100, epsilon, gamma, Nmax);
	test_ss << m2d(test_opt.x) << ";" << m2d(test_opt.y) << ";" << test_opt.f_calls << ";" << (test_opt.x > -1 && test_opt.x < 1 ? "lokalne" : "globalne") << ";\n";
	solution::clear_calls();

	//Zapis do pliku
#ifdef SAVE_TO_FILE
	save_to_file("test_no_expansion.csv", test_ss.str());
#endif

#endif

#ifdef CALC_SIMULATION

	matrix ud1 = matrix(6, 1);
	ud1(0) = 0.5; //Pa
	ud1(1) = 90.0; //Ta
	ud1(2) = 1.0; //Pb
	ud1(3) = 36.5665 * 0.0001; //Db
	ud1(4) = 10 * 0.001; //F_in
	ud1(5) = 20.0; //T_in

	
	double Da_0 = 1.0 * 0.0001;

	
	double* bounds = expansion(ff1R, Da_0, d, alpha, Nmax, ud1);
	solution::clear_calls();

	std::cout << bounds[0] << ", " << bounds[1] << "\n\n";

	solution opt = fib(ff1R, bounds[0], bounds[1], epsilon, ud1);
	std::cout << opt;
	solution::clear_calls();

	matrix Y0 = matrix(3, 1);
	Y0(0) = 5.0; //Poczatkowa objetosc w a
	Y0(1) = 1.0; //Poczatkowa objetosc w b
	Y0(2) = 20.0;//Poczatkowa temperatura w b


	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, opt.x);

#ifdef SAVE_TO_FILE
	save_to_file("simulation_fibonacci.csv", hcat(Y[0], Y[1]));
#endif

	int n = get_len(Y[0]);
	double Tb_max = Y[1](0, 2);
	for (int i = 0; i < n; ++i)
	{
		if (Tb_max < Y[1](i, 2))
			Tb_max = Y[1](i, 2);
	}

	std::cout << "Tb_max (fibonacci): " << Tb_max << "\n\n";

	
	opt = lag(ff1R, bounds[0], bounds[1], epsilon, epsilon, Nmax, ud1);
	std::cout << opt;
	solution::clear_calls();

	Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, opt.x);

#ifdef SAVE_TO_FILE
	save_to_file("simulation_lagrange.csv", hcat(Y[0], Y[1]));
#endif

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
