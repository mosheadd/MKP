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

		float e;
		float a;
		float n;

		int nyu;
		int T;

	public:


		SpaceObject(float e, float a, int nyu);

		float get_e();
		float get_a();
		float get_n();

		int get_nyu();
		int get_T();

	};


	float KeplerEquation(float M, float e, float E);

	float trueAnomaly(float e, float E);


	float bisection(float M, float e, float epsilon, int max_it);

	float goldensection(float M, float e, float epsilon, int max_it);

	float fixedpoint(float M, float e, float epsilon, int max_it);

	float newton(float M, float e, float epsilon, int max_it);

	void rootsOfKepEq(SpaceObject& object, float(*func)(float, float, float, int), float epsilon, int max_it);


}
