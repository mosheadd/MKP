#pragma once
#include "iostream"
#include "cmath"
#include "vector"
#include "gnuplot-iostream.h"


namespace mkp
{

	
	const float PI = 3.14159;
	const float kappa = 1.618;


	class SpaceObject
	{

		double e;
		double n;

		float a;

		int nyu;
		int T;

	public:


		SpaceObject(double e, float a, int nyu);

		double get_e();
		double get_n();

		float get_a();

		int get_nyu();
		int get_T();

	};


	double KeplerEquation(double M, double e, double E);

	double trueAnomaly(double e, double E);


	double bisection(double M, double e, double epsilon, int max_it);

	double goldensection(double M, double e, double epsilon, int max_it);

	double fixedpoint(double M, double e, double epsilon, int max_it);

	double newton(double M, double e, double epsilon, int max_it);

	void rootsOfKepEq(SpaceObject& object, double(*func)(double, double, double, int), double epsilon, int max_it);


}
