#pragma once
#include "iostream"
#include "cmath"
#include "vector"


namespace mkp
{


	const float PI = 3.14159;
	const float kappa = 1.618;
	const float G = 6.67e-11;


	class SpaceObject
	{

		float e;
		float a;
		float n;

		unsigned _int64 mass;
		int T;


	public:
		SpaceObject(unsigned _int64 mass, float e, float a);

		float get_e();
		float get_a();
		float get_n();

		int get_T();

	};


	float KeplerEquation(float M, float e, float E);


	float bisection(float M, float e, float epsilon, int max_it);

}
