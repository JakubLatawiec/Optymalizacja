#pragma once

#include"ode_solver.h"

matrix ff1(matrix, matrix = NAN, matrix = NAN);
matrix ff0T(matrix, matrix = NAN, matrix = NAN);
matrix ff0R(matrix, matrix = NAN, matrix = NAN);
matrix ff2T(matrix, matrix = NAN, matrix = NAN);
matrix df0(double, matrix, matrix = NAN, matrix = NAN);

//lab_2
double simulate(double k1, double k2);
double calculate_torque(double k1, double k2, double alpha, double omega);
matrix Q_function(matrix k, matrix, matrix);
