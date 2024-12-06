#include "MKP.h"


namespace mkp
{

	//Реализация класса космического объекта

	SpaceObject::SpaceObject(double e, float a, int nyu)
	{

		this->e = e;
		this->a = a;
		this->nyu = nyu;


		T = 2 * PI * sqrt(pow(a, 3) / nyu);

		n = 2 * PI / T;

		std::cout << "T: " << T << "\n";
		std::cout << "n: " << n << "\n\n";

	}

	double SpaceObject::get_e() { return e; }
	double SpaceObject::get_n() { return n; }

	float SpaceObject::get_a() { return a; }

	int SpaceObject::get_nyu() { return nyu; }
	int SpaceObject::get_T() { return T; }


	//Остальный функции

	double KeplerEquation(double M, double e, double E) { return E - e * sin(E) - M; }

	double trueAnomaly(double e, double E)
	{

		double V = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));

		if (V < 0)
			V += 2 * PI;

		return V;

	}


}
