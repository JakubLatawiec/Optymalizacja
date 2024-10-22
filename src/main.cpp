/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

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
		lab2();
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
	// Parametry pocz¹tkowe
	double step_size = 0.5;  // Pocz¹tkowy krok (dla Hooke'a-Jeevesa)
	matrix s0(2, 1, 0.5);    // Wektor pocz¹tkowych d³ugoœci kroków (dla Rosenbrocka)
	double alpha = 1.5;      // Wspó³czynnik ekspansji (dla Rosenbrocka)
	double beta = 0.5;       // Wspó³czynnik kontrakcji (dla Rosenbrocka)
	double epsilon = 1e-6;   // Dok³adnoœæ
	int Nmax = 1000;         // Maksymalna liczba wywo³añ funkcji celu
	matrix ud1, ud2;         // Dodatkowe dane u¿ytkownika (nieu¿ywane tutaj)

	// Wypisanie nag³ówka tabeli
	std::cout << "D³ugoœæ kroku  Lp.  x1^(0)  x2^(0)  x1* (HJ)  x2* (HJ)  y* (HJ)  Liczba wyw. HJ  Min globalne HJ  x1* (R)  x2* (R)  y* (R)  Liczba wyw. R  Min globalne R" << std::endl;

	// Przeprowadzamy 100 iteracji
	for (int iter = 1; iter <= 100; ++iter) {
		// Zmieniamy punkt startowy dla ka¿dej iteracji
		double start_values[2] = { 1.0 + iter * 0.1, 1.0 + iter * 0.1 };  // Przyk³adowa zmiana
		matrix x0(2, start_values);  // Punkt pocz¹tkowy jako macierz 2x1

		// Wywo³anie algorytmu Hooke'a-Jeevesa (HJ)
		solution result_hj = HJ(ff2T, x0, step_size, alpha, epsilon, Nmax, ud1, ud2);

		// Wywo³anie algorytmu Rosenbrocka
		solution result_rosen = Rosen(ff2T, x0, s0, alpha, beta, epsilon, Nmax, ud1, ud2);

		// Sprawdzamy, czy znaleziono minimum globalne
		std::string global_min_hj = (result_hj.y(0, 0) <= 1e-3) ? "TAK" : "NIE";
		std::string global_min_rosen = (result_rosen.y(0, 0) <= 1e-3) ? "TAK" : "NIE";

		// Wypisujemy dane bez dodatkowego formatowania
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
