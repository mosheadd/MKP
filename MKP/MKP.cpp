#include "MKP.h"


namespace mkp
{

	//Реализация класса космического объекта

	SpaceObject::SpaceObject(unsigned _int64 mass, float e, float a)
	{

		this->mass = mass;
		this->e = e;
		this->a = a;

		float nyu = G * mass;

		T = 2 * PI * sqrt(pow(a, 3) / nyu);

		std::cout << T;

		n = 2 * PI / T;

	}

	float SpaceObject::get_e() { return e; }
	float SpaceObject::get_a() { return a; }
	float SpaceObject::get_n() { return n; }

	int SpaceObject::get_T() { return T; }


	//Остальный функции

	float KeplerEquation(float M, float e, float E) { return E - e * sin(E) - M; }

	float trueAnomaly(float e, float E)
	{

		float V = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));

		return V;

	}


}
