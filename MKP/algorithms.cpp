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


//Функция нахождения самого решения 
void mkp::findRootsOfKepEq(SpaceObject& object, float(*func)(float, float, float, int), float epsilon, int max_it)
{


	//Данные графика
	std::vector<std::pair<float, float>> M_graph;
	std::vector<std::pair<float, float>> E_graph;
	std::vector<std::pair<float, float>> v_graph;
	M_graph.emplace_back(0.f, 0.f);
	E_graph.emplace_back(0.f, 0.f);
	v_graph.emplace_back(0.f, 0.f);


	//Данные объекта
	float e = object.get_e();
	float a = object.get_a();
	float n = object.get_n();

	int T = object.get_T();


	//Сам алгоритм
	int iteration_delta = 500;

	float M, E, v;

	for (int t = 500; t < T; t += iteration_delta)
	{

		M = n * t;

		M_graph.emplace_back((float)t, M);


		E = func(M, e, epsilon, max_it);

		E_graph.emplace_back((float)t, E);


		v = trueAnomaly(e, E);

		v_graph.emplace_back((float)t, v);


		std::cout << v << " " << E << " " << M << std::endl;


	}

	M = n * T;
	
	M_graph.emplace_back((float)T, M);


	E = func(M, e, epsilon, 10000);

	E_graph.emplace_back((float)T, E);


	v = trueAnomaly(e, E);

	v_graph.emplace_back((float)T, v);

	std::cout << v << " " << E << " " << M << "\n\n";


	std::cout << "Ploting...\n";


	try
	{

		Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");

		gp << "set title 'График аномалий объекта'\n";

		gp << "plot ";
		gp << "'-' with lines smooth mcsplines title 'V', ";
		gp << "'-' with lines smooth mcsplines title 'E', ";
		gp << "'-' with lines smooth mcsplines title 'M'\n";

		gp.send1d(v_graph);
		gp.send1d(E_graph);
		gp.send1d(M_graph);

		std::cout << "Success...\n";

		system("pause");

	}
	catch (...)
	{

		std::cout << "Error ploting...\n\n";

	}


}


//Реализация графика радиус вектора
void mkp::findRadVec(SpaceObject& object, float(*func)(float, float, float, int), float epsilon, int max_it)
{

	//Данные графика для радиус-вектора
	std::vector <std::pair<float, float>> r_graph;


	//Данные объекта
	float e = object.get_e();
	float a = object.get_a();
	float n = object.get_n();

	int T = object.get_T();


	//Фокальный праметр
	float p = a * (1 - e * e);
	std::cout <<"\n" << p << "\n";


	float r0 = p / (1 + e);
	r_graph.emplace_back(0.f, r0);


	float M, E, v, r;

	//Расчет

	int iteration_delta = 500;

	for (int t = 500; t < T; t+=iteration_delta)
	{

		M = n * t;

		E = func(M, e, epsilon, max_it);

		v = trueAnomaly(e, E);

		r = p / (1 + e * cos(v));

		//std::cout << V << "\n";

		r_graph.emplace_back(t, r);

	}

	M = n * T;

	E = func(M, e, epsilon, max_it);

	v = trueAnomaly(e, E);

	r = p / (1 + e * cos(v));

	r_graph.emplace_back(T, r);


	std::cout << "Ploting...\n";

	try
	{

		Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");

		
		gp << "set title 'График радиус-вектора'\n";

		gp << "plot '-' with lines smooth mcsplines title 'r'\n";

		
		gp.send1d(r_graph);


		std::cout << "Success...\n";

		system("pause");


	}
	catch (...)
	{

		std::cout << "Error...\n";

	}
	

}

void mkp::findVelocities(SpaceObject& object, float(*func)(float, float, float, int), float epsilon, int max_it, int nyu)
{

	//Данные графиков скоростей
	std::vector<std::pair<float, float>> Vr_graph;
	std::vector<std::pair<float, float>> Vn_graph;
	std::vector<std::pair<float, float>> V_graph;


	//Данные объекта
	float e = object.get_e();
	float a = object.get_a();
	float n = object.get_n();

	int T = object.get_T();


	//Фокальный праметр
	float p = a * (1 - e * e);

	Vr_graph.emplace_back(0.f, 0.f);
	Vn_graph.emplace_back(0.f, 0.f);
	V_graph.emplace_back(0.f, 0.f);
	

	float M, E, v, V, Vr, Vn;

	int iteration_delta = 500;


	for (int t = 500; t < T; t+= iteration_delta)
	{

		M = n * t;

		E = func(M, e, epsilon, max_it);

		v = trueAnomaly(e, E);

		Vr = sqrt(nyu / p) * e * sin(v);

		Vr_graph.emplace_back((float)t, Vr);

		Vn = sqrt(nyu / p) * (1 + e * cos(v));

		Vn_graph.emplace_back((float)t, Vr);

		V = sqrt(Vr * Vr + Vn * Vn);

		V_graph.emplace_back((float)t, V);

		std::cout << Vr << " " << Vn << " " << V << "\n";

	}

	M = n * T;

	E = func(M, e, epsilon, max_it);

	v = trueAnomaly(e, E);

	Vr = sqrt(nyu / p) * e * sin(v);

	Vr_graph.emplace_back((float)T, Vr);

	Vn = sqrt(nyu / p) * (1 + e * cos(v));

	Vn_graph.emplace_back((float)T, Vr);

	V = sqrt(Vr * Vr + Vn * Vn);

	V_graph.emplace_back((float)T, V);


	std::cout << "Ploting...\n";

	try
	{

		Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");

		gp << "set title 'Графики скоростей'\n";

		gp << "plot ";
		gp << "'-' with lines smooth mcsplines title 'Vr', ";
		gp << "'-' with lines smooth mcsplines title 'Vn', ";
		gp << "'-' with lines smooth mcsplines title 'V'\n";

		gp.send1d(Vr_graph);
		gp.send1d(Vn_graph);
		gp.send1d(V_graph);

		std::cout << "Success...\n";

		system("pause");

	}
	catch (...)
	{

		std::cout << "Error...\n";

	}


}
