#include "Simplex.h"
#include <iostream>
#include <vector>


int main() {
	std::vector<std::vector<double>> equations = { {0, 5, 1, 0, 0, 15},
												   {6, 2, 0, 1, 0, 24},
												   {1, 1, 0, 0, 1, 5},
												   {2, 1, 0, 0, 0, 0} };
	// Get a Simplex Object
	Simplex s(2, 3, 0, equations);
	// Output
	std::cout<<s.Find_Val()<<std::endl;
}