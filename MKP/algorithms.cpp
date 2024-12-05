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

			//std::cout << "Not reached\n";
			break;

		}

	}

	abPares retValue = seqOfSegments_ab[seqOfSegments_ab.size() - 1];

	return (retValue.a + retValue.b) / 2;

}


//Реализации метода золотого сечения
float mkp::goldensection(float M, float e, float epsilon, int max_it)
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

			//std::cout << "Not reached\n";
			break;

		}

	}

	abPares retValue = seqOfSegments_ab[seqOfSegments_ab.size() - 1];

	return retValue.a + (retValue.b - retValue.a) / kappa;

}


//Реализация метода итераций
float mkp::fixedpoint(float M, float e, float epsilon, int max_it)
{
	
	float E0 = M;
	float E = e * sin(E0) + M;

	for (int i = 0; i < max_it; i++)
	{

		if (abs(E - E0) < epsilon)
			return E;

		E0 = E;
		E = e * sin(E0) + M;

	}

	std::cout << "Max number of iterations reached\n";

	return E;

}


//Реализация метода Ньютона
float mkp::newton(float M, float e, float epsilon, int max_it)
{
	
	float E = M;

	for (int i = 0; i < max_it; i++)
	{

		float f = KeplerEquation(M, e, E);
		float df = 1 - e * cos(E);


		float dE = -(f / df);

		E += dE;


		if (abs(dE) < epsilon)
			return E;

	}


	std::cout << "Max number of iterations reached\n";

	return E;


}


//Функция нахождения всех параметров
void mkp::rootsOfKepEq(SpaceObject& object, float(*func)(float, float, float, int), float epsilon, int max_it)
{

	//Для построения графиков
	std::vector<std::pair<float, float>> M_graph;
	std::vector<std::pair<float, float>> E_graph;
	std::vector<std::pair<float, float>> v_graph;

	M_graph.emplace_back(0.f, 0.f);
	E_graph.emplace_back(0.f, 0.f);
	v_graph.emplace_back(0.f, 0.f);


	std::vector <std::pair<float, float>> r_graph;


	std::vector<std::pair<float, float>> Vr_graph;
	std::vector<std::pair<float, float>> Vn_graph;
	std::vector<std::pair<float, float>> V_graph;


	//Данные объекта
	float e = object.get_e();
	float a = object.get_a();
	float n = object.get_n();

	int nyu = object.get_nyu();
	int T = object.get_T();


	//Фокальный праметр
	float p = a * (1 - e * e);


	float r0 = p / (1 + e);
	r_graph.emplace_back(0.f, r0);


	float Vn0 = sqrt(nyu / p) * (1 + e);

	Vr_graph.emplace_back(0.f, 0.f);
	Vn_graph.emplace_back(0.f, Vn0);
	V_graph.emplace_back(0.f, Vn0);


	//Расчет
	float M, E, v, r, Vr, Vn, V;

	int iteration_delta = 500;

	for (int t = iteration_delta; t < T; t += iteration_delta)
	{

		M = n * t;

		M_graph.emplace_back((float)t, M);


		E = func(M, e, epsilon, max_it);

		E_graph.emplace_back((float)t, E);


		v = trueAnomaly(e, E);

		v_graph.emplace_back((float)t, v);


		r = p / (1 + e * cos(v));

		r_graph.emplace_back(t, r);


		Vr = sqrt(nyu / p) * e * sin(v);

		Vr_graph.emplace_back((float)t, Vr);


		Vn = sqrt(nyu / p) * (1 + e * cos(v));

		Vn_graph.emplace_back((float)t, Vn);


		V = sqrt(Vr * Vr + Vn * Vn);

		V_graph.emplace_back((float)t, V);


	}

	M = n * T;

	M_graph.emplace_back((float)T, M);


	E = func(M, e, epsilon, 10000);

	E_graph.emplace_back((float)T, E);


	v = trueAnomaly(e, E);

	v_graph.emplace_back((float)T, v);


	r = p / (1 + e * cos(v));

	r_graph.emplace_back(T, r);


	Vr = sqrt(nyu / p) * e * sin(v);

	Vr_graph.emplace_back((float)T, Vr);


	Vn = sqrt(nyu / p) * (1 + e * cos(v));

	Vn_graph.emplace_back((float)T, Vn);


	V = sqrt(Vr * Vr + Vn * Vn);

	V_graph.emplace_back((float)T, V);


	Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");


	while (1)
	{

		std::cout << "1.Plot anomalies.\n";
		std::cout << "2.Plot position.\n";
		std::cout << "3.Plot velocities.\n";
		std::cout << "4.Output all in files.\n";
		std::cout << "(Other num).Exit.\n";

		int choice;

		std::cin >> choice;


		switch(choice)
		{
		case 1:

			gp << "set title 'График аномалий объекта'\n";

			gp << "plot ";
			gp << "'-' with lines smooth mcsplines title 'v', ";
			gp << "'-' with lines smooth mcsplines title 'E', ";
			gp << "'-' with lines smooth mcsplines title 'M'\n";

			gp.send1d(v_graph);
			gp.send1d(E_graph);
			gp.send1d(M_graph);

			break;
		case 2:

			gp << "set title 'График радиус-вектора'\n";

			gp << "plot '-' with lines smooth mcsplines title 'r'\n";

			gp.send1d(r_graph);
			break;
		case 3:

			gp << "set title 'Графики скоростей'\n";

			gp << "plot ";
			gp << "'-' with lines smooth mcsplines title 'Vr', ";
			gp << "'-' with lines smooth mcsplines title 'Vn', ";
			gp << "'-' with lines smooth mcsplines title 'V'\n";

			gp.send1d(Vr_graph);
			gp.send1d(Vn_graph);
			gp.send1d(V_graph);
			
			break;
		case 4:

			break;
		default:
			return;
		}


	}

	std::cout << "Success...\n\n";


}

