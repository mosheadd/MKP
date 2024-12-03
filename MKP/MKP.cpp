#include "MKP.h"


namespace mkp
{


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

	
	//Алгоритмы вычисления корня уравнения Кеплеа

	float bisection(float M, float e, float epsilon, int max_it)
	{

		struct abPares
		{

			float a;
			float b;

			abPares(float a, float b) : a(a), b(b) {}

		};


		float a = 0;
		float b = 2 * PI;

		std::vector<abPares> seqOfSegments_ab;


		for (int i = 0; i < max_it; i++)
		{

			float E = (a + b) / 2;

			float f = KeplerEquation(M, e, E);

			if (f == 0)
				return E;

			float fa = KeplerEquation(M, e, a);
			float fb = KeplerEquation(M, e, b);

			if (fa * fb < 0)
				seqOfSegments_ab.push_back(abPares(a, b));


			if (fa * f < 0)
				b = E;
			else
				a = E;


			if (abs(b - a) < 2 * epsilon)
			{

				std::cout << "Не достигнуто нужного количества итераций.\n";
				break;

			}

		}

		abPares retValue = seqOfSegments_ab[seqOfSegments_ab.size() - 1];

		return (retValue.a + retValue.b) / 2;

	}


	//Остальный функции

	float KeplerEquation(float M, float e, float E) { return E - e * sin(E) - M; }

	float trueAnomaly(float e, float E)
	{

		float V = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));

		return V;

	}


}
