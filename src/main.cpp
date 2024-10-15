#include<iostream>
#include<random>
#include<chrono>
#include<fstream>
#include"opt_alg.h"
#include"ode_solver.h"
#include <iomanip>

using namespace std;

int main()
{
	try
	{
		cout << "LAB NUMBER " << LAB_NO << endl;
		cout << "LAB PART " << LAB_PART << endl << endl;

#if LAB_NO==2
		solution F, L;
		random_device R;
		double x0 = 200.0 * R() / R.max() - 100;
		double d = 6.2;
		double tabx0[100];
		double taba[100];
		double tabb[100];
		int tabExpfcalls[100];
		double tabFx[100];
		double tabFy[100];
		int tabFibfcalls[100];
		double tabLx[100];
		double tabLy[100];
		int tabLagfcalls[100];
		//cout << "x0 = " << x0 << endl;
		double alfa = 2.2;
		int Nmax = 45;
		double* p = expansion(x0, d, alfa, Nmax);
		solution::clear_calls();
		double a, b, epsilon;
		epsilon = 0.01;
		solution C();


		F = fib(-100, 100, epsilon);
		cout << "Fibonaciego x=" << F.x << "  y=" << F.y << endl;
		cout << "Fibonaciego fCalls=" << solution::f_calls << endl;
		solution::clear_calls();
		//double a, b, epsilon,
		double	gamma;
		epsilon = 0.01;
		gamma = 0, 0001;

		//lag(a, b, epsilon, gamma, Nmax);
		L = lag(-100, 100, epsilon, gamma, Nmax);
		cout << "Lagrange x=" << L.x << "  y=" << L.y << endl;
		cout << "Lagrange fCalls=" << solution::f_calls << endl;
		cout << endl << endl << endl;



		for (int i = 0; i < 100; i++) {
			cout << "Iteracja :" << i << endl;
			solution::clear_calls();
			x0 = 200.0 * R() / R.max() - 100;
			cout << "x0 = " << x0 << endl;
			tabx0[i] = x0;
			p = expansion(x0, d, alfa, Nmax);
			cout << "przedzial a=" << p[0] << "  b=" << p[1] << endl;
			taba[i] = p[0];
			tabb[i] = p[1];
			cout << "liczba wywolan: " << solution::f_calls << endl;
			tabExpfcalls[i] = solution::f_calls;
			solution::clear_calls();
			F = fib(p[0], p[1], epsilon);
			cout << "Fibonaciego x=" << F.x << "  y=" << F.y << endl;
			tabFx[i] = F.x(0);
			tabFy[i] = F.y(0);
			cout << "liczba wywolan: " << solution::f_calls << endl;
			tabFibfcalls[i] = solution::f_calls;
			solution::clear_calls();
			L = lag(p[0], p[1], epsilon, gamma, Nmax);
			cout << "Lagrange x=" << L.x << "  y=" << L.y << endl;
			tabLx[i] = L.x(0);
			tabLy[i] = L.y(0);
			cout << "liczba wywolan: " << solution::f_calls << endl;
			tabLagfcalls[i] = solution::f_calls;
		}
		cout << "tak";
		ofstream myfile;
		myfile.open("X0.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx0[i] << "\n";
		myfile.close();


		myfile.open("Expa.txt");
		for (int i = 0; i < 100; i++)
			myfile << taba[i] << "\n";
		myfile.close();

		myfile.open("Expb.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabb[i] << "\n";
		myfile.close();

		myfile.open("tabExpfcalls.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabExpfcalls[i] << "\n";
		myfile.close();

		myfile.open("tabFx.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabFx[i] << "\n";
		myfile.close();

		myfile.open("tabFy.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabFy[i] << "\n";
		myfile.close();

		myfile.open("tabFibfcalls.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabFibfcalls[i] << "\n";
		myfile.close();

		myfile.open("tabLx.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabLx[i] << "\n";
		myfile.close();

		myfile.open("tabLy.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabLy[i] << "\n";
		myfile.close();

		myfile.open("tabLagfcalls.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabLagfcalls[i] << "\n";
		myfile.close();

#endif
#if LAB_NO==3
		random_device R;
		double step = 0.001;
		matrix Step = matrix(new double[2]{ 0.001,0.001 }, 2);
		matrix X = matrix(new double[2]{ 2.0,2.0 }, 2);
		X(0) = 2.0 * R() / R.max() - 1.0;
		X(1) = 2.0 * R() / R.max() - 1.0;
		cout << "start" << X(0) << " " << X(1) << endl;
		solution HJs = HJ(X, step, 0.5, 0.01, 1000);
		cout << "End HJ" << HJs.x(0) << " " << HJs.x(1) << endl;
		cout << HJs.y << endl;
		cout << "F calls: " << solution::f_calls << endl;
		solution::clear_calls;



		solution Ros = Rosen(X, Step, 2, 0.6, 0.01, 1000);
		cout << "Krok Ros" << Ros.x(0) << " " << Ros.x(1) << endl;
		cout << Ros.y << endl;
		cout << "F calls: " << solution::f_calls << endl;

		double tabx10[100];
		double tabx20[100];

		double tabHJx1[100];
		double tabHJx2[100];
		double tabHJy[100];
		double tabHJfcalls[100];

		double tabRosx1[100];
		double tabRosx2[100];
		double tabRosy[100];
		double tabRosfcalls[100];



		for (int i = 0; i < 100; i++) {
			cout << "\n Iteracja :" << i << endl;
			solution::clear_calls();
			X(0) = 2.0 * R() / R.max() - 1.0;
			X(1) = 2.0 * R() / R.max() - 1.0;
			tabHJ[0] = 2.0 * R() / R.max() - 1.0;
			tabHJ[1] = 2.0 * R() / R.max() - 1.0;
			matrix x0 = (tabHJ, 2);
			cout << X(0) << " " << X(1) << endl;
			tabx10[i] = X(0);
			tabx20[i] = X(1);

			HJs = HJ(X, step, 0.5, 0.01, 1000);
			cout << "HJ: " << HJs.x(0) << " " << HJs.x(1) << endl;
			cout << HJs.y << endl;
			cout << "F calls: " << solution::f_calls << endl;
			tabHJx1[i] = HJs.x(0);
			tabHJx2[i] = HJs.x(1);
			tabHJy[i] = HJs.y(0);
			tabHJfcalls[i] = solution::f_calls;

			solution::clear_calls;


			Ros = Rosen(X, Step, 2, 0.6, 0.01, 1000);
			cout << "Rosen: " << Ros.x(0) << " " << Ros.x(1) << endl;
			cout << Ros.y << endl;
			cout << "F calls: " << solution::f_calls << endl;
			tabRosx1[i] = Ros.x(0);
			tabRosx2[i] = Ros.x(1);
			tabRosy[i] = Ros.y(0);
			tabRosfcalls[i] = solution::f_calls;
		}


		ofstream myfile;
		myfile.open("X10.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx10[i] << "\n";
		myfile.close();


		myfile.open("X20.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx20[i] << "\n";
		myfile.close();


		myfile.open("tabHJx1.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabHJx1[i] << "\n";
		myfile.close();


		myfile.open("tabHJx2.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabHJx2[i] << "\n";
		myfile.close();


		myfile.open("tabHJy.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabHJy[i] << "\n";
		myfile.close();


		myfile.open("tabHJfcalls.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabHJfcalls[i] << "\n";
		myfile.close();


		myfile.open("tabRosx1.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabRosx1[i] << "\n";
		myfile.close();


		myfile.open("tabRosx2.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabRosx2[i] << "\n";
		myfile.close();


		myfile.open("tabRosy.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabRosy[i] << "\n";
		myfile.close();


		myfile.open("tabRosfcalls.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabRosfcalls[i] << "\n";
		myfile.close();
#endif
#if LAB_NO==4
		solution::g_calls = 0;
		double c = 10.0;
		double epsilon = 0.00001;
		int Nmax = 1000000;
		double a = 5;
		random_device r;
		double dc = 1.7;
		double tabx0[100];
		double tabx1[100];
		double tabx1zew[100];
		double tabx2zew[100];
		double tabyzew[100];
		double tabFcallszew[100];
		matrix X = matrix(new double[2]{ 2.0,2.0 }, 2);
		for (int i = 0; i < 100; i++) {
			cout << "\n ZEW Iteracja :" << i << endl;
			X(0) = (5.0 * r()) / r.max() + 1;
			X(1) = (5.0 * r()) / r.max() + 1;
			tabx0[i] = X(0);
			tabx1[i] = X(1);
			solution::clear_calls();
			matrix O(a);
			solution Pen = pen(X, c, dc, epsilon, Nmax, O);
			cout << "x1= " << X(0) << "   x2= " << X(1) << endl;
			cout << Pen << endl;
			tabx1zew[i] = Pen.x(0);
			tabx2zew[i] = Pen.x(1);
			tabyzew[i] = Pen.y(0);
			tabFcallszew[i] = solution::f_calls;
			solution::g_calls++;
		}

		ofstream myfile;
		myfile.open("tabx0.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx0[i] << "\n";
		myfile.close();


		myfile.open("tabx1.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx1[i] << "\n";
		myfile.close();


		myfile.open("tabx1zew.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx1zew[i] << "\n";
		myfile.close();


		myfile.open("tabx2zew.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx2zew[i] << "\n";
		myfile.close();


		myfile.open("tabyzew.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabyzew[i] << "\n";
		myfile.close();


		myfile.open("tabFcallszew.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabFcallszew[i] << "\n";
		myfile.close();




		cout << "LAB PART " << LAB_PART << endl << endl;
#define LAB_PART 2
		dc = 0.5;
		cout << "LAB PART " << LAB_PART << endl << endl;
		double tabx1wew[100];
		double tabx2wew[100];
		double tabywew[100];
		double tabFcallswew[100];

		for (int i = 0; i < 100; i++) {
			cout << "\n WEW Iteracja :" << i << endl;
			solution::clear_calls();
			matrix O(a);
			X(0) = tabx0[i];
			X(1) = tabx1[i];
			solution Pen = pen(X, c, dc, epsilon, Nmax, O);
			cout << "x1= " << X(0) << "   x2= " << X(1) << endl;
			cout << Pen << endl;
			tabx1wew[i] = Pen.x(0);
			tabx2wew[i] = Pen.x(1);
			tabywew[i] = Pen.y(0);
			tabFcallswew[i] = solution::f_calls;
			solution::g_calls++;
		}



		myfile.open("tabx1wew.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx1wew[i] << "\n";
		myfile.close();


		myfile.open("tabx2wew.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx2wew[i] << "\n";
		myfile.close();


		myfile.open("tabywew.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabywew[i] << "\n";
		myfile.close();


		myfile.open("tabFcallswew.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabFcallswew[i] << "\n";
		myfile.close();



#endif
#if LAB_NO==5 && LAB_PART == 1

		double tabx0[10];
		double tabx1[10];

		double tabx1SD[300];
		double tabx2SD[300];
		double tabySD[300];
		double tabFcallSD[300];
		double tabGcallSD[300];

		double tabx1CG[300];
		double tabx2CG[300];
		double tabyCG[300];
		double tabFcallCG[300];
		double tabGcallCG[300];


		double tabx1Newton[300];
		double tabx2Newton[300];
		double tabyNewton[300];
		double tabFcallNewton[300];
		double tabGcallNewton[300];
		double tabHcallNewton[300];

		random_device R;
		double epsilon = 1e-5;
		int Nmax = 10000;
		matrix O(2, 2);
		O(0, 0) = -10;
		O(0, 1) = 10;
		O(1, 0) = -10;
		O(1, 1) = 10;
		matrix x0(new double[2]{ 2, 2 }, 2);
		double h0;
		solution Newtone;
		solution CGe;
		solution SDe;
		int j = 0;
		for (int i = 0; i < 300; i++)
		{
			cout << "\n WEW Iteracja :" << i << endl;

			h0 = 0.05;
			cout << "h0" << h0 << endl;
			x0(0) = 20.0 * R() / R.max() - 10;
			x0(1) = 20.0 * R() / R.max() - 10;
			tabx0[j] = x0(0);
			tabx1[j] = x0(1);

			cout << "pkt startowy: " << x0(0) << ", " << x0(1) << endl;
			cout << "SDe:\n";

			solution::clear_calls();
			SDe = SD(x0, h0, epsilon, Nmax, O);
			cout << SDe << "\n\n";
			/// ////////////////////////////////////////////////////////////////////////////////////////
			tabx1SD[i] = SDe.x(0);
			tabx2SD[i] = SDe.x(1);
			tabySD[i] = SDe.y(0);
			tabFcallSD[i] = solution::f_calls;
			tabGcallSD[i] = solution::g_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////

			cout << "Newton:\n";
			solution::clear_calls();
			Newtone = Newton(x0, h0, epsilon, Nmax, O);
			cout << Newtone << "\n\n";


			/// ////////////////////////////////////////////////////////////////////////////////////////
			tabx1Newton[i] = Newtone.x(0);
			tabx2Newton[i] = Newtone.x(1);
			tabyNewton[i] = Newtone.y(0);
			tabFcallNewton[i] = solution::f_calls;
			tabGcallNewton[i] = solution::g_calls;
			tabHcallNewton[i] = solution::H_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////
			cout << "CGe:\n";
			solution::clear_calls();
			CGe = CG(x0, h0, epsilon, Nmax, O);
			cout << CGe << "\n\n";
			/// ////////////////////////////////////////////////////////////////////////////////////////
			tabx1CG[i] = CGe.x(0);
			tabx2CG[i] = CGe.x(1);
			tabyCG[i] = CGe.y(0);
			tabFcallCG[i] = solution::f_calls;
			tabGcallCG[i] = solution::g_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////
			i++;
			h0 = 0.12;
			cout << "h0" << h0 << endl;
			x0(0) = 20.0 * R() / R.max() - 10;
			x0(1) = 20.0 * R() / R.max() - 10;
			cout << "pkt startowy: " << x0(0) << ", " << x0(1) << endl;
			cout << "SDe:\n";


			solution::clear_calls();
			SDe = SD(x0, h0, epsilon, Nmax, O);
			cout << SDe << "\n\n";

			/// ////////////////////////////////////////////////////////////////////////////////////////
			tabx1SD[i] = SDe.x(0);
			tabx2SD[i] = SDe.x(1);
			tabySD[i] = SDe.y(0);
			tabFcallSD[i] = solution::f_calls;
			tabGcallSD[i] = solution::g_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////

			cout << "Newton:\n";
			solution::clear_calls();
			Newtone = Newton(x0, h0, epsilon, Nmax, O);
			cout << Newtone << "\n\n";

			/// ////////////////////////////////////////////////////////////////////////////////////////
			tabx1Newton[i] = Newtone.x(0);
			tabx2Newton[i] = Newtone.x(1);
			tabyNewton[i] = Newtone.y(0);
			tabFcallNewton[i] = solution::f_calls;
			tabGcallNewton[i] = solution::g_calls;
			tabHcallNewton[i] = solution::H_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////

			cout << "CGe:\n";
			solution::clear_calls();
			CGe = CG(x0, h0, epsilon, Nmax, O);
			cout << CGe << "\n\n";

			/// ////////////////////////////////////////////////////////////////////////////////////////
			tabx1CG[i] = CGe.x(0);
			tabx2CG[i] = CGe.x(1);
			tabyCG[i] = CGe.y(0);
			tabFcallCG[i] = solution::f_calls;
			tabGcallCG[i] = solution::g_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////
			i++;
			h0 = -0.12;
			cout << "h0" << h0 << endl;
			x0(0) = 20.0 * R() / R.max() - 10;
			x0(1) = 20.0 * R() / R.max() - 10;
			cout << "pkt startowy: " << x0(0) << ", " << x0(1) << endl;
			cout << "SDe:\n";


			solution::clear_calls();
			SDe = SD(x0, h0, epsilon, Nmax, O);
			cout << SDe << "\n\n";


			/// ////////////////////////////////////////////////////////////////////////////////////////
			tabx1SD[i] = SDe.x(0);
			tabx2SD[i] = SDe.x(1);
			tabySD[i] = SDe.y(0);
			tabFcallSD[i] = solution::f_calls;
			tabGcallSD[i] = solution::g_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////

			cout << "Newton:\n";
			solution::clear_calls();
			Newtone = Newton(x0, h0, epsilon, Nmax, O);
			cout << Newtone << "\n\n";

			/// ////////////////////////////////////////////////////////////////////////////////////////
			tabx1Newton[i] = Newtone.x(0);
			tabx2Newton[i] = Newtone.x(1);
			tabyNewton[i] = Newtone.y(0);
			tabFcallNewton[i] = solution::f_calls;
			tabGcallNewton[i] = solution::g_calls;
			tabHcallNewton[i] = solution::H_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////

			cout << "CGe:\n";
			solution::clear_calls();
			CGe = CG(x0, h0, epsilon, Nmax, O);
			cout << CGe << "\n\n";


			/// ////////////////////////////////////////////////////////////////////////////////////////

			tabx1CG[i] = CGe.x(0);
			tabx2CG[i] = CGe.x(1);
			tabyCG[i] = CGe.y(0);
			tabFcallCG[i] = solution::f_calls;
			tabGcallCG[i] = solution::g_calls;
			/// ////////////////////////////////////////////////////////////////////////////////////////

			j++;
		}
		ofstream myfile;

		myfile.open("tabx0.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx0[i] << "\n";
		myfile.close();

		myfile.open("tabx1.txt");
		for (int i = 0; i < 100; i++)
			myfile << tabx1[i] << "\n";
		myfile.close();





		myfile.open("tabx1SD.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabx1SD[i] << "\n";
		myfile.close();

		myfile.open("tabx2SD.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabx2SD[i] << "\n";
		myfile.close();

		myfile.open("tabySD.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabySD[i] << "\n";
		myfile.close();

		myfile.open("tabFcallSD.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabFcallSD[i] << "\n";
		myfile.close();

		myfile.open("tabGcallSD.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabGcallSD[i] << "\n";
		myfile.close();





		myfile.open("tabx1CG.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabx1CG[i] << "\n";
		myfile.close();

		myfile.open("tabx2CG.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabx2CG[i] << "\n";
		myfile.close();

		myfile.open("tabyCG.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabyCG[i] << "\n";
		myfile.close();

		myfile.open("tabFcallCG.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabFcallCG[i] << "\n";
		myfile.close();

		myfile.open("tabGcallCG.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabGcallCG[i] << "\n";
		myfile.close();



		myfile.open("tabx1Newton.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabx1Newton[i] << "\n";
		myfile.close();

		myfile.open("tabx2Newton.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabx2Newton[i] << "\n";
		myfile.close();

		myfile.open("tabyNewton.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabyNewton[i] << "\n";
		myfile.close();

		myfile.open("tabFcallNewton.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabFcallNewton[i] << "\n";
		myfile.close();

		myfile.open("tabGcallNewton.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabGcallNewton[i] << "\n";
		myfile.close();


		myfile.open("tabHcallNewton.txt");
		for (int i = 0; i < 300; i++)
			myfile << tabHcallNewton[i] << "\n";
		myfile.close();

		/*double epsilon = 0.00001;
		int Nmax = 1000000;
		matrix X = matrix(new double[2]{ 2.0,2.0 }, 2);
		double h = 0.05;
		X(0) = -8;
		X(1) = 8;
		solution SDe = SD(X, h, epsilon, Nmax);
		cout << SDe << endl;


		X(0) = -8;
		X(1) = 8;
		solution CGe = CG(X, h, epsilon, Nmax);
		cout << CGe << endl;



		X(0) = -8;
		X(1) = 8;
		solution Newtone = Newton(X, h, epsilon, Nmax);
		cout << Newtone << endl;*/

#endif

#if LAB_NO == 5 && LAB_PART == 2
		double epsilon = 0.0000001;
		int Nmax = 10000;
		matrix x0(new double[3]{ 0,0,0, }, 3);
		cout << "= = = = = = = = = = = = = = =\n";
		solution::clear_calls();
		solution CG1 = CG(x0, 0.01, epsilon, Nmax);
		cout << CG1;
		cout << "= = = = = = = = = = = = = = =\n";
		solution::clear_calls();
		solution CG2 = CG(x0, 0.001, epsilon, Nmax);
		cout << CG2;
		cout << "= = = = = = = = = = = = = = =\n";
		solution::clear_calls();
		solution CG3 = CG(x0, 0.0001, epsilon, Nmax);
		cout << CG3;

#endif

	}
	catch (char* EX_INFO)
	{
		cout << EX_INFO << endl;
	}
	system("pause");
	return 0;
}
