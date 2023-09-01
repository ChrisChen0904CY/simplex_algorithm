#include "Simplex.h"
#include <vector>
#include <iostream>
#include <map>
#include <algorithm>

// Initialize
Simplex::Simplex(int ori_var,
				 int loo_var = 0,
				 int man_var = 0,
				 std::vector<std::vector<double>> M = {}) {
	this->ori_var = ori_var;
	this->tar_var = loo_var;
	this->loo_var = loo_var;
	this->man_var = man_var;
	this->variables = std::vector<bool>(ori_var+loo_var+man_var, false);
	this->M = M;
	int r = M.size();
	int c = M[0].size();
	// Initialize the Coefficients of Target Function
	for (int i = 0; i < c - 1; i++) {
		this->tar_coe.push_back(this->M[r - 1][i]);
	}
	// Generate the First Feasible Solution
	for (int i = ori_var; i < ori_var+loo_var; i++) {
		// The First x variable will be a element of base solution
		// x here refers to the value of target variables
		this->variables[i]=true;
	}
	this->Compute_Base();
}
// Destruct
Simplex::~Simplex() {

}

void Simplex::Compute_Base() {
	int c = this->M[0].size();
	std::vector<std::vector<double>> S;
	// Only pull target variables in to S to get their Value
	int i, j;
	for (i = 0; i < this->M.size()-1.; i++) {
		std::vector<double> row = {};
		for (j = 0; j < c-1; j++) {
			if (this->variables[j]==true) {
				row.push_back(this->M[i][j]);
			}
		}
		// Get the values column
		row.push_back(this->M[i][c - 1]);
		S.push_back(row);
	}
	// this->output_matrix(S);
	// Get the value of base variables
	std::vector<double> solutions = this->Gauss_Jordan(S);
	j = 0;
	this->sol_vec = {};
	for (i = 0; i < c-1; i++) {
		if (this->variables[i]==true) {
			this->sol_vec.push_back(solutions[j++]);
		}
		else {
			this->sol_vec.push_back(0);
		}
		// std::cout << sol_vec[i] << " ";
	}
	// std::cout << std::endl;
}

// Optimality Test
bool Simplex::Opt_Test() {
	int r = this->M.size();
	int c = this->M[0].size();
	for (int i = 0; i < c - 1; i++) {
		if (this->variables[i]) {
			if (this->M.back()[i] != 0) {
				return false;
			}
		}
		else {
			if (this->M.back()[i] > 0) {
				return false;
			}
		}
	}
	return true;
}

// Find Wanted Value
double Simplex::Find_Val() {
	int r = this->M.size();
	int c = this->M[0].size();
	// While current Matrix can't pass the optimality
	while (!this->Opt_Test()) {
		this->output_matrix(this->M);
		// Index of the Enter Variable and the Leave Variable
		int inInd, outInd;
		inInd = outInd = 0;
		int maxx = INT_MIN;
		int minn = INT_MAX;
		// Find the inInd
		for (int i = 0; i < c - 1; i++) {
			if (!this->variables[i]) {
				if (this->M[r - 1][i] > maxx) {
					maxx = this->M[r - 1][i];
					inInd = i;
				}
			}
		}
		// Find the outInd
		for (int i = 0; i < r-1; i++) {
			if (this->M[i][inInd] > 0) {
				if (this->M[i][c-1] / this->M[i][inInd] < minn) {
					minn = this->M[i][c-1] / this->M[i][inInd];
					outInd = i;
				}
			}
		}
		// Disable the Prior Base Variable who's about to leave
		int i = 0;
		for (auto v : this->variables) {
			if (v) {
				i++;
			}
			if (i == outInd + 1) {
				v = false;
				break;
			}
		}
		// Enter & Leave on Variables
		this->variables[inInd] = true;
		// Update M
		for (int i = 0; i < r; i++) {
			double coe = this->M[i][inInd] / this->M[outInd][inInd];
			if (i != outInd) {
				for (int j = 0; j < c; j++) {
					this->M[i][j] -= coe* this->M[outInd][j];
				}
				this->M[i][inInd] = 0;
			}
		}
		for (int j = 0; j < c; j++) {
			if (j != inInd) {
				this->M[outInd][j] /= this->M[outInd][inInd];
			}
		}
		this->M[outInd][inInd] = 1;
		// Update the values of Base Viriables
		this->Compute_Base();
		// Update the target-value
		double Z = 0;
		for (int i = 0; i < c - 1; i++) {
			if (this->variables[i]==true) {
				Z += this->tar_coe[i]*this->sol_vec[i];
			}
		}
		this->M.back().back() = Z;
	}
	// Print the Base Variables
	for (int i = 0; i < c - 1; i++) {
		if (this->variables[i]) {
			std::cout << "X" << i + 1 << "> " << this->sol_vec[i] << std::endl;
		}
	}
	return this->M.back().back();
}

// Gaussian Elemination
std::pair<std::vector<std::vector<double>>, bool> Simplex::Gaussian(std::vector<std::vector<double>> mat) {
	std::map<int, std::vector<double>> pivots_row;
	std::vector<int> pivots;
	std::vector<std::vector<double>> ans;
	int rev = 0;
	int n = mat.size();
	int pivot = 0;
	for (int i = 0; i < n-1; i++) {
		pivot = 0;
		for (int j = 0; j < n; j++) {
			if (mat[i][j] != 0) {
				// Find the first non-zero element
				pivot = j;
				pivots.push_back(pivot);
				break;
			}
		}
		// Eliminates by rows
		for (int k = i + 1; k < n; k++) {
			// Compute the coeffient
			double coe = mat[k][pivot]/ mat[i][pivot];
			for (int j = 0; j < n; j++) {
				mat[k][j] -= coe * mat[i][j];
			}
			mat[k][pivot] = 0;
		}
		// Deal with Current row
		for (int j = 0; j < n; j++) {
			if (j!=pivot) {
				mat[i][j] = 0;
			}
		}
		pivots_row[pivot] = mat[i];
	}
	// Deal with the last row
	for (int i = 0; i < n; i++) {
		if (mat[n - 1][i] != 0) {
			pivot = i;
			pivots.push_back(pivot);
			break;
		}
	}
	pivots_row[pivot] = mat[n - 1];
	// Combine the rows sorted by pivots
	for (int i = 0; i < n; i++) {
		ans.push_back(pivots_row[i]);
	}
	// Compute the value of reverse
	for (int i = 0; i < n-1; i++) {
		for (int j = i + 1; j < n; j++) {
			if (pivots[j] < pivots[i]) {
				rev++;
			}
		}
	}
	return { ans, (rev &1)==0 };
}

// Compute Det
double Simplex::Det(std::vector<std::vector<double>> mat) {
	std::pair<std::vector<std::vector<double>>, bool> Gass = this->Gaussian(mat);
	double ans = 1.0;
	for (int i = 0; i < mat.size(); i++) {
		ans *= Gass.first[i][i];
		ans *= Gass.second ? 1 : (-1);
	}
	return ans;
}

// Solve Equation Groups
std::vector<double> Simplex::Solve(std::vector<std::vector<double>> mat) {
	std::vector<double> solutions;
	int m = mat.size();
	int n = mat[0].size();
	// Use to det computation
	std::vector<std::vector<double>> comp = {};
	std::vector<double> tmp = {};
	// Slice the wanted part of matrix
	std::vector<double>::const_iterator left;
	std::vector<double>::const_iterator right;
	// Sliced by colums
	for (int i = 0; i < m; i++) {
		tmp = {};
		left = mat[i].begin();
		right = mat[i].end() - 1;
		tmp.assign(left, right);
		comp.push_back(tmp);
	}
	double D = this->Det(comp);
	// Compute for Each variable
	for (int i = 0; i < m; i++) {
		// Slice First
		comp = {};
		for (int r = 0; r < m; r++) {
			tmp = {};
			// Left Part
			if (i != 0) {
				left = mat[r].begin();
				right = mat[r].begin() + i;
				tmp.insert(tmp.end(), left, right);
			}
			// Current Column
			tmp.push_back(mat[r][n-1]);
			// Right Part
			left = mat[r].begin() + i + 1;
			right = mat[r].end()-1;
			tmp.insert(tmp.end(), left, right);
			comp.push_back(tmp);
		}
		// Get the value
		solutions.push_back(this->Det(comp)/D);
	}
	return solutions;
}

// Gauss-Jordan Elemination
std::vector<double> Simplex::Gauss_Jordan(std::vector<std::vector<double>> mat) {
	int r = mat.size();
	int c = mat[0].size();
	std::vector<std::vector<double>> tri_mat;
	if (r < c - 1) {
		std::cout << "Equations are less to solve." << std::endl;
	}
	std::map<int, std::vector<double>> pivots_row;
	std::vector<int> pivots;
	std::vector<std::vector<double>> ans;
	std::vector<double> solutions;
	int pivot = 0;
	// Pivot Matching
	for (int i = 0; i < r - 1; i++) {
		pivot = 0;
		for (int j = 0; j < c; j++) {
			if (mat[i][j] != 0) {
				// Find the first non-zero element
				pivot = j;
				pivots.push_back(pivot);
				break;
			}
		}
		// Eliminates by rows
		for (int k = i + 1; k < r; k++) {
			// Compute the coeffient
			double coe = mat[k][pivot] / mat[i][pivot];
			for (int j = 0; j < c; j++) {
				mat[k][j] -= coe * mat[i][j];
			}
			mat[k][pivot] = 0;
		}
		pivots_row[pivot] = mat[i];
	}
	// Deal with the last row
	for (int i = 0; i < c-1; i++) {
		if (mat[r - 1][i] != 0) {
			pivot = i;
			pivots.push_back(pivot);
			break;
		}
	}
	pivots_row[pivot] = mat[r-1];
	// Combine the rows sorted by pivots
	for (int i = 0; i < r; i++) {
		tri_mat.push_back(pivots_row[i]);
	}
	// this->output_matrix(tri_mat);
	// Triangle --> Unit
	for (int i = r-1; i >= 0; i--) {
		for (int j = c-2; j > i; j--) {
			tri_mat[i][c-1] -= tri_mat[i][j]* tri_mat[j][c-1];
			tri_mat[i][j] = 0;
		}
		tri_mat[i][c - 1] /= tri_mat[i][i];
		solutions.push_back(tri_mat[i][c - 1]);
		// std::cout << tri_mat[i][c - 1] << " ";
		tri_mat[i][i] = 1.0;
	}
	// this->output_matrix(tri_mat);
	reverse(solutions.begin(), solutions.end());
	return solutions;
}

// Output the Matrix
void Simplex::output_matrix(std::vector<std::vector<double>> mat) {
	for (auto r : mat) {
		for (auto x : r) {
			std::cout << x << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
