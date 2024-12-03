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

			std::cout << "Not reached\n";
			break;

		}

	}

	abPares retValue = seqOfSegments_ab[seqOfSegments_ab.size() - 1];

	return (retValue.a + retValue.b) / 2;

}


//Реализации метода золотого сечения
float mkp::goldenratio(float M, float e, float epsilon, int max_it)
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

		float E = a + (b - a) / kappa;

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

			std::cout << "Not reached\n";
			break;

		}

	}

	abPares retValue = seqOfSegments_ab[seqOfSegments_ab.size() - 1];

	return retValue.a + (retValue.b - retValue.a) / kappa;

}

float mkp::fixedpoint(float M, float e, float epsilon, int max_it)
{
	
	float E0 = M;
	float E = e * sin(E0) + M;

	for (int i = 0; i < max_it; i++)
	{

		if (abs(E - E0) < epsilon)
			return E;

	}

	std::cout << "Max number of iterations reached\n";

	return E;

}



//Функция нахождения решения 
void mkp::findRootsOfKepEq(SpaceObject& object, float(*func)(float, float, float, int), float epsilon)
{

	//Данные графика
	std::vector<float> M_axis;
	std::vector<float> E_axis;
	std::vector<float> V_axis;

	std::vector<int> t_axis;


	//Данные объекта
	float e = object.get_e();
	float a = object.get_a();
	float n = object.get_n();

	int T = object.get_T();


	//Сам алгоритм
	int iteration_delta = 500;

	float M, E, V;

	for (int t = 500; t < T; t += iteration_delta)
	{

		M = n * t;

		M_axis.push_back(M);


		E = func(M, e, epsilon, 10000);

		E_axis.push_back(E);


		V = trueAnomaly(e, E);

		V_axis.push_back(E);


		std::cout << V << " " << E << " " << M << std::endl;


	}

	M = n * T;
	
	M_axis.push_back(M);


	E = func(M, e, epsilon, 10000);

	E_axis.push_back(E);


	V = trueAnomaly(e, E);

	V_axis.push_back(V);

	std::cout << V << " " << E << " " << M << std::endl;


}
