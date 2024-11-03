/*********************************************
Kod stanowi uzupełnienie materiałów do ćwiczeń
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"
#include <cmath>
#include <vector>
void lab0();
void lab1();

void lab2();
void tabela1();
void tabela2();
void tabela3();
void tabela4();
void tabela5();


void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		tabela2();
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
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	double x0 = 1;
	double d = 0.1;
	double alpha = 2.0;
	int nMax = 100;

	//for na 100 (randowanie x0, i wypis do pliku)
	double* reach = expansion(ff1, x0, d, alpha, nMax);
	std::cout << "Res: " << reach[0] << ", " << reach[1] << "\n\n";
	std::cout << solution::f_calls << "\n";
	solution::clear_calls();

	double a = reach[0];
	double b = reach[1];
	double epsilon = 0.0001;
	double gamma = 0.0001;

	solution sol = fib(ff1, a, b, epsilon);
	std::cout << sol << std::endl;

	solution::clear_calls();

	solution sol2 = lag(ff1, a, b, epsilon, gamma, 100);
	std::cout << sol2;

	solution::clear_calls();


	//Problem rzeczywisty
	double Pa = 0.5;
	double Va_0 = 5.0;
}

void lab2()
{
	double step_size = 0.5;  
	matrix s0(2, 1, 0.5);    
	double alpha = 1.5;     
	double beta = 0.5;       
	double epsilon = 1e-6;   
	int Nmax = 1000;         
	matrix ud1, ud2;         

	std::cout << "Długość kroku  Lp.  x1^(0)  x2^(0)  x1* (HJ)  x2* (HJ)  y* (HJ)  Liczba wyw. HJ  Min globalne HJ  x1* (R)  x2* (R)  y* (R)  Liczba wyw. R  Min globalne R" << std::endl;

	for (int iter = 1; iter <= 200; ++iter) {
		double start_values[2] = { 1.0 + iter * 0.1, 1.0 + iter * 0.1 };  
		matrix x0(2, start_values);  

		solution result_hj = HJ(ff2T, x0, step_size, alpha, epsilon, Nmax, ud1, ud2);

		solution result_rosen = Rosen(ff2T, x0, s0, alpha, beta, epsilon, Nmax, ud1, ud2);

		std::string global_min_hj = (result_hj.y(0, 0) <= 1e-3) ? "TAK" : "NIE";
		std::string global_min_rosen = (result_rosen.y(0, 0) <= 1e-3) ? "TAK" : "NIE";

		std::cout << step_size << "  "
			<< iter << "  "
			<< x0(0, 0) << "  "
			<< x0(1, 0) << "  "
			<< result_hj.x(0, 0) << "  "
			<< result_hj.x(1, 0) << "  "
			<< result_hj.y(0, 0) << "  "
			<< result_hj.f_calls << "  "
			<< global_min_hj << "  "
			<< result_rosen.x(0, 0) << "  "
			<< result_rosen.x(1, 0) << "  "
			<< result_rosen.y(0, 0) << "  "
			<< result_rosen.f_calls << "  "
			<< global_min_rosen
			<< std::endl;
	
	}

}


void imulation()
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


void tabela1() {
	std::cout << "Długość kroku  Lp.  x1^(0)  x2^(0)  x1* (HJ)  x2* (HJ)  y* (HJ)  Wywołania HJ  Globalne HJ  "
		"x1* (R)  x2* (R)  y* (R)  Wywołania R  Globalne R\n";
	for (int iter = 1; iter <= 100; ++iter) {
		double start_values[2] = { 1.0 + iter * 0.1, 1.0 + iter * 0.1 };
		matrix x0(2, start_values);

		solution result_hj = HJ(Q_function, x0, 0.5, 1.5, 1e-6, 1000, matrix(), matrix());
		solution result_rosen = Rosen(Q_function, x0, matrix(2, 1, 0.5), 1.5, 0.5, 1e-6, 1000, matrix(), matrix());

		std::string global_hj = (result_hj.y(0, 0) <= 1e-3) ? "TAK" : "NIE";
		std::string global_rosen = (result_rosen.y(0, 0) <= 1e-3) ? "TAK" : "NIE";

		std::cout << "0.5  " << iter << "  " << x0(0, 0) << "  " << x0(1, 0) << "  "
			<< result_hj.x(0, 0) << "  " << result_hj.x(1, 0) << "  " << result_hj.y(0, 0) << "  "
			<< result_hj.f_calls << "  " << global_hj << "  "
			<< result_rosen.x(0, 0) << "  " << result_rosen.x(1, 0) << "  " << result_rosen.y(0, 0) << "  "
			<< result_rosen.f_calls << "  " << global_rosen << "\n";
	}
}


void tabela2() {
	int global_min_hj = 0, global_min_rosen = 0;
	for (int iter = 1; iter <= 100; ++iter) {
		std::cout << iter;
		double start_values[2] = { 1.0 + iter * 0.1, 1.0 + iter * 0.1 };
		matrix x0(2, start_values);

		solution result_hj = HJ(Q_function, x0, 0.5, 1.5, 1e-6, 1000, matrix(), matrix());
		solution result_rosen = Rosen(Q_function, x0, matrix(2, 1, 0.5), 1.5, 0.5, 1e-6, 1000, matrix(), matrix());

		if (result_hj.y(0, 0) <= 1e-3) global_min_hj++;
		if (result_rosen.y(0, 0) <= 1e-3) global_min_rosen++;
	}

	std::cout << "Liczba globalnych minimów\n";
	std::cout << "Metoda Hooke'a-Jeevesa: " << global_min_hj << "\n";
	std::cout << "Metoda Rosenbrocka: " << global_min_rosen << "\n";
}

void tabela3() {
	std::cout << "Nr iteracji  x1* (HJ)  x2* (HJ)  x1* (R)  x2* (R)\n";
	for (int iter = 0; iter < 100; ++iter) {
		double start_values[2] = { 1.0 + iter * 0.1, 1.0 + iter * 0.1 };
		matrix x0(2, start_values);

		solution result_hj = HJ(Q_function, x0, 0.5, 1.5, 1e-6, 1000, matrix(), matrix());
		solution result_rosen = Rosen(Q_function, x0, matrix(2, 1, 0.5), 1.5, 0.5, 1e-6, 1000, matrix(), matrix());

		std::cout << iter << "  " << result_hj.x(0, 0) << "  " << result_hj.x(1, 0) << "  "
			<< result_rosen.x(0, 0) << "  " << result_rosen.x(1, 0) << "\n";
	}
}

void tabela4() {
	std::cout << "Długość kroku  k1* (HJ)  k2* (HJ)  Q* (HJ)  Wywołania HJ  k1* (R)  k2* (R)  Q* (R)  Wywołania R\n";
	double start_values[2] = { 5.0, 5.0 };
	matrix x0(2, start_values);

	solution result_hj = HJ(Q_function, x0, 0.5, 1.5, 1e-6, 1000, matrix(), matrix());
	solution result_rosen = Rosen(Q_function, x0, matrix(2, 1, 0.5), 1.5, 0.5, 1e-6, 1000, matrix(), matrix());

	std::cout << "0.5  " << result_hj.x(0, 0) << "  " << result_hj.x(1, 0) << "  " << result_hj.y(0, 0) << "  "
		<< result_hj.f_calls << "  " << result_rosen.x(0, 0) << "  " << result_rosen.x(1, 0) << "  "
		<< result_rosen.y(0, 0) << "  " << result_rosen.f_calls << "\n";
}

void tabela5(double k1_hj, double k2_hj, double k1_rosen, double k2_rosen) {
	std::cout << "t  α (HJ)  ω (HJ)  α (R)  ω (R)\n";
	double dt = 0.1;
	double alpha_hj = 0.0, omega_hj = 0.0;
	double alpha_rosen = 0.0, omega_rosen = 0.0;
}