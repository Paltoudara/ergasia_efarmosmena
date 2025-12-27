//more generic code next time
#pragma once
#include<iostream>
#include<cmath>
#include<numbers>
#include<vector>
#include<functional>
#include<optional>
//
class solution {
public:
	//func_f done
	double func_f(const double x) {
		return (4 * std::pow(std::numbers::pi, 2) + 2) * std::sin(2 * std::numbers::pi * x);
	}
	//u_exact done
	double u_exact(const double x) {
		return std::sin(2 * std::numbers::pi * x);
	}
	//trapezoid_integral done
	double trapezoid_integral(const std::function<double(double)>& func,const double a,const double b,const double h) {
		std::size_t n = static_cast<std::size_t>((b - a) / h);
		double sum{ (func(a) + func(b)) / 2 };
		for (std::size_t i = 1; i < n; i++) {
			sum += func(a + i * h);
		}
		return sum * h;
	}
    double fi1(const double x) {
        return x * x - x;
    }
    double fi2(const double x) {
        return x * x * x - x;
    }
    //done 
    std::vector<double> gauss_elim(std::vector<std::vector<double>>& A, const std::vector<double>& b) {//tes reference edo giati kaneis copy 
        const size_t n = A.size(); // Number of equations (rows)
        if (n == 0)return {};
        // Augment the coefficient matrix with the right-hand side vector
        for (size_t i = 0; i < n; i++) {
            A[i].emplace_back(b[i]); // Appending the corresponding element from b to each row of A
        }//epaujhmenos
        std::cout << '\n';
        // Perform Gaussian elimination
        for (size_t i = 0; i < n; i++) { // Loop over each row (equation)
            // Find the row with the maximum absolute value in the ith column and swap rows
            size_t max_row = i;
            for (size_t j = i + 1; j < n; j++) {
                if (std::abs(A[j][i]) > std::abs(A[max_row][i])) {
                    max_row = j; // Update the index of the row with the maximum absolute value
                }
            }
            std::swap(A[i], A[max_row]); // Swap the current row with the row with the maximum absolute value
            // Perform row operations to eliminate coefficients below the pivot element
            for (size_t j = i + 1; j < n; j++) { // Loop over rows below the pivot row
                if (std::abs(A[i][i]-0.0) < 1E-8)
                    return {};
                double factor = A[j][i] / A[i][i]; // Compute the factor by which the pivot row will be multiplied
                for (size_t k = i; k < n + 1; k++) { // Loop over columns including the augmented column
                    A[j][k] -= factor * A[i][k]; // Perform row operation to eliminate coefficients below the pivot element
                }
            }
        }
        if (std::abs(A[n - 1][n - 1]- 0.0)<1E-8) {
            return {};
        }
         // Back-substitution to solve for x
        std::vector<double> x(n, 0.0); // Initialize the solution vector x with zeros
        for (size_t i = n - 1; i > 0; i--) { // Start from the last equation and move upwards
            double sum = 0.0;
            for (size_t j = i + 1; j < n; j++) { // Loop over elements to the right of the diagonal in the current row
                sum += A[i][j] * x[j]; // Compute the sum of products of coefficients and corresponding elements of x
            }
            if (std::abs(A[i][i] - 0.0) < 1E-8) {
                return {};
            }
            x[i] = (A[i][n] - sum) / A[i][i]; // Compute the value of x[i] using back-substitution
        }
        size_t i = 0;
        while (i == 0) {
            double sum = 0.0;
            for (size_t j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            if (std::abs(A[i][i] - 0.0) < 1E-8) {
                return {};
            }
            x[i] = (A[i][n] - sum) / A[i][i];
            break;
        }
        return  x; // Return the solution vector x
    }
    void error(const  std::vector<double>&sol,const double a,const double b,const double h) {
        //return un
        std::vector<std::function<double(double)>>fi{ };
        fi.push_back([&](double x)->double {return fi1(x); });
        fi.push_back([&](double x)->double {return fi2(x); });
        std::size_t n = static_cast<std::size_t>((b - a) / h);
        for (std::size_t i = 0; i <= n; i++) {
            double sum{ 0.0 };
            for (std::size_t j = 0; j < 2; j++) {
                sum += sol[j] * fi[j](a + i * h);
            }
           /* sum += sol[0] * fi1(a + i * h);
            sum += sol[1] * fi2(a + i * h);*/
            std::cout << "un: " << sum << " u_exact:" << u_exact(a + i * h) << " error: " << std::abs(sum - u_exact(a + i * h)) << '\n';
            //std::cout << std::abs(sum - u_exact(a + i*h))<<'\n';
        }
       
    }

};
