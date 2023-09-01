#pragma once
#ifndef SIMPLEX_H_
#define SIMPLEX_H_
#endif

#include <vector>
#include <iostream>
#include<map>
#include <algorithm>

class Simplex
{
private:
	// Original Variables
	int ori_var;
	// Target Variables
	int tar_var;
	// Loose Variables
	int loo_var;
	// Manual Variables
	int man_var;
	// Optimal Matrix
	std::vector<std::vector<double>> M;
	// Base and Other Variables
	std::vector<bool> variables;
	// Solution Vector
	std::vector<double> sol_vec;
	// Coefficients of Target Function
	std::vector<double> tar_coe;

public:
	Simplex(int,
			int,
			int,
			std::vector<std::vector<double>>);
	~Simplex();
	// Optimality Test
	bool Opt_Test();
	// Get Value of Base Variables
	void Compute_Base();
	// Find Wanted Value
	double Find_Val();

	// Gaussian Elemination
	std::pair<std::vector<std::vector<double>>, bool> Gaussian(std::vector<std::vector<double>>);
	// Compute Det
	double Det(std::vector<std::vector<double>>);
	// Solve Equation Groups
	std::vector<double> Solve(std::vector<std::vector<double>>);
	// Gaussian-Jordan Elemination
	std::vector<double> Gauss_Jordan(std::vector<std::vector<double>>);
	// Output the Matrix
	void output_matrix(std::vector<std::vector<double>>);
};
