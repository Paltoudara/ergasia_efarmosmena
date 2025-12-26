#include"header.h"

int main() {
	//elegxos
	solution obj{};
	double a{ 0.0 }, b{ 1.0 }, h{ 0.1 }, n{2};
	std::vector<std::vector<double>>matrix{};
	std::vector<std::function<double(double)>>functions{};
	functions.push_back([](double x)->double{return 2.0 * x*x - 2.0 * x - 2.0; });
	functions.push_back([](double x)->double {return 2.0 * x*x*x - 8.0 * x; });
	/*
	functions[0]->lf1
	functions[1]->lf2
	*/
	std::vector<double>B{};
	for (std::size_t i = 0; i < n; i++) {
		std::vector<double>A{};
		for (std::size_t j = 0; j < n; j++) {
			A.push_back(obj.trapezoid_integral([&](double x) {return functions[i](x) * functions[j](x); },a,b,h));
		}
		matrix.push_back(A);
	}
	
	for (std::size_t i = 0; i < n; i++) {
		B.push_back(obj.trapezoid_integral([&](double x) {return obj.func_f(x)*functions[i](x); },a,b,h));
	}
	std::vector<double>sol = obj.gauss_elim(matrix, B);
	if (sol.empty()) {
		std::cout << "this system has no solution\n";
		return 0;
	}
	for (std::size_t i = 0; i < sol.size(); i++) {
		std::cout << sol[i] << '\n';
	}
	std::cout << '\n';
	obj.error(sol,a,b,h);
}
