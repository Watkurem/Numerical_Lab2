/* Copyright 2016 Alexander Melnyk / Олександр Мельник
 *
 * This file is part of Numerical_Lab2 package.
 *
 * Numerical_Lab2 is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * Numerical_Lab2 is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * Numerical_Lab2. If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include "table.h"

/* This comment was for my personal reference, is not very useful and is saved
   here only for historical purposes
 */

/* x*x+2*sin(x)-1 = 0
 * x*x = 1-2*sin(x)
 *
 * f(x) = x*x+2*sin(x)-1;
 * Графічно визначаю, що корені знаходяться на (-2, -1.5) і (0.5, 1)
 *
 * Для методу ітерацій знаходжу похідну функції: 2(x+cos(x)).
 * На відрізку [-2, -1.5] похідна зростає; в такому випадку
 *  M1 = f'(-1.5) = -2.859;
 *  m1 = f'(-2) = -4.832;
 * На відрізку [0.5, 1] похідна зростає; в такому випадку
 *  M2 = f'(1) = 3.081;
 *  m2 = f'(0.5) = 2.755;
 * λ1 = 1/M1 = -0.3498;
 * λ2 = 1/M2 = 0.3246;
 * q1 = 1-m1/M1 = -0.6905l
 * q2 = 1-m2/M2 = 0.1056;
 *
 * Отже, φ1(x) = x - λ1*f(x);
 *       φ2(x) = x - λ2*f(x).
 *
 * Для методу дотичних m1 визначаю на [-2, -1.5] як |f'(-2)|
 *                     m1 на [0.5, 1] як |f'(0.5)|
 */

#define COLUMNS 3
#define ROWS 5

#define X1_A (-2)
#define X1_B (-1.5)
#define X2_A 0.5
#define X2_B 1

#define LAMBDA1 -0.34983069634458664308
#define LAMBDA2  0.32461160260238120922
#define Q1 -0.69048466060011397130
#define Q2  0.10564143373534422962
#define CONTROL_X1 -1.72517120542893012713
#define CONTROL_X2  0.42302818188516042885

#define X 1 /* Which root to find, 1st or 2nd */

/* This crazy stuff is needed to allow the program to be easily recompiled for
   finding first or second roots. This way is, most definitely, not the best
   way (and probably not even a good way), but it allowed me to experience
   several preprocessor quirks.
 */
#define XCONTROL_X(N) CONTROL_X ## N
#define CONTROL_X(N) XCONTROL_X(N)
#define XLAMBDA(N) LAMBDA ## N
#define LAMBDA(N) XLAMBDA(N)
#define XQ(N) Q ## N
#define Q(N) XQ(N)
#define XX_A(N) X ## N ## _A
#define X_A(N) XX_A(N)
#define XX_B(N) X ## N ## _B
#define X_B(N) XX_B(N)

/* Left part of eq */
double f(const double x){
	return x*x+2*sin(x)-1;
}

/* First derivative */
double df(const double x){
	return 2*(x+cos(x));
}

/* Second derivative */
double ddf(const double x){
	return 2 - 2*sin(x);
}

/* Formula implementation */
double phi(const double x, const double lambda){
	return x - lambda*f(x);
}

/* Another formula implementation (iterative method), recursive */
double iterate(const double x, const double eps, const double lambda,
               const double q, unsigned long *n){
	double n_x = phi(x, lambda);
	(*n)++;
	if(fabs(n_x - x) <= fabs(eps*(1-q)/q))
		return n_x;
	return iterate(n_x, eps, lambda, q, n);
}

/* Another formula implementation (Newton method), recursive part */
double newton_recur(const double x, const double eps, const double m1, unsigned long *n){
	double n_x = x - f(x)/df(x);
	(*n)++;
	if(fabs(f(n_x)/m1) <= eps)
		return n_x;
	return newton_recur(n_x, eps, m1, n);
}

/* Another formula implementation (Newton method), non-recursive part */
double newton_method(const double eps, const double a, const double b, const double m1, unsigned long *n){
	if(f(a)*ddf(a) > 0)
		return newton_recur(a, eps, m1, n);
	else if(f(b)*ddf(b) > 0)
		return newton_recur(b, eps, m1, n);
	else return NAN;
}

int main(){
	int i, j;
	char *header1[COLUMNS] = {" eps ",
	                          "  Значення кореня  ",
	                          "Оцінка точності кореня за методом І"};
	char *header2[COLUMNS] = {" eps ",
	                          "  Значення кореня  ",
	                          "Оцінка точності кореня за методом Д"};
	char *header3[COLUMNS] = {" eps ",
	                          "Кількість ітерацій за методом І",
	                          "Кількість ітерацій за методом Д"};
	char ***rows;
	double epses[ROWS] = {0.01};
	double x0 = (X_A(X)+X_B(X))/2;
	double answer = 0;
	unsigned long n_iteration[ROWS];
	unsigned long n_newton[ROWS];

	setlocale(LC_ALL, ""); /* For mbstowcs */

	/* Init output tables */
	rows = calloc(ROWS, sizeof(char **));
	for(j = 0; j < ROWS; j++){
		rows[j] = calloc(COLUMNS, sizeof(char *));
		for(i = 0; i < COLUMNS; i++){
			/* This here makes sure that there will be no overflows when initing
			 * headers yet a sane amount of memory will be used.
			 * mbstowcs is used so that extra bytes will not be allocated for the
			 * non-header rows (which, as far as my concept goes, should never
			 * contain non-ascii strings, while headers can be non-ascii)
			 */
			rows[j][i] = calloc(mbstowcs(NULL, header1[i], 0), 1);
		}
	}

	/* Fill output tables with computed data */
	for(i = 0; i < ROWS; i++){
		if(i)
			epses[i] = epses[i-1]*0.001;
		n_iteration[i] = 0;
		answer = iterate(x0, epses[i], LAMBDA(X), Q(X), &(n_iteration[i]));
		sprintf(rows[i][0], "%.0e", epses[i]);
		sprintf(rows[i][1], "%19.16f", answer);
		sprintf(rows[i][2], "%35.29e", fabs(answer-CONTROL_X(X)));
	}

	print_table(header1, COLUMNS, rows, ROWS);

	/* Free used memory */
	for(j = 0; j < ROWS; j++){
		for(i = 0; i < COLUMNS; i++){
			free(rows[j][i]);
		}
		free(rows[j]);
	}
	free(rows);

	/* Init output tables for second table */
	rows = calloc(ROWS, sizeof(char **));
	for(j = 0; j < ROWS; j++){
		rows[j] = calloc(COLUMNS, sizeof(char *));
		for(i = 0; i < COLUMNS; i++){
			rows[j][i] = calloc(mbstowcs(NULL, header2[i], 0), 1);
		}
	}

	/* Fill output tables with computed data */
	for(i = 0; i < ROWS; i++){
		n_newton[i] = 0;
		answer = newton_method(epses[i], X_A(X), X_B(X), df(X_A(X)), &(n_newton[i]));
		sprintf(rows[i][0], "%.0e", epses[i]);
		sprintf(rows[i][1], "%19.16f", answer);
		sprintf(rows[i][2], "%35.29e", fabs(answer-CONTROL_X(X)));
	}

	print_table(header2, COLUMNS, rows, ROWS);

	/* Free used memory */
	for(j = 0; j < ROWS; j++){
		for(i = 0; i < COLUMNS; i++){
			free(rows[j][i]);
		}
		free(rows[j]);
	}
	free(rows);

	/* Init output tables for third table */
	rows = calloc(ROWS, sizeof(char **));
	for(j = 0; j < ROWS; j++){
		rows[j] = calloc(COLUMNS, sizeof(char *));
		for(i = 0; i < COLUMNS; i++){
			rows[j][i] = calloc(mbstowcs(NULL, header3[i], 0), 1);
		}
	}

	/* Fill output tables with computed data */
	for(i = 0; i < ROWS; i++){
		sprintf(rows[i][0], "%.0e", epses[i]);
		sprintf(rows[i][1], "%31lu", n_iteration[i]);
		sprintf(rows[i][2], "%31lu", n_newton[i]);
	}

	print_table(header3, COLUMNS, rows, ROWS);

	/* Free used memory */
	for(j = 0; j < ROWS; j++){
		for(i = 0; i < COLUMNS; i++){
			free(rows[j][i]);
		}
		free(rows[j]);
	}
	free(rows);

	/* I should've probably done error detection too, but this program can
	 * crash only if there was not enough memory, and if you don't have enough
	 * memory for something like this then you're clearly doing something
	 * really wrong.*/
	return 0;
}
