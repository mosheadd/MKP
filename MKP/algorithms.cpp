#include "MKP.h"


using namespace mkp;


//Реализация метода половинного деления
float mkp::bisection(float M, float e, float epsilon, int max_it)
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




//Функция нахожждения решения 
void mkp::findRootsOfKepEq(SpaceObject& object)
{



}
