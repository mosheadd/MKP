#include "MKP.h"
#include "fstream"


using namespace mkp;


//Реализация метода половинного деления
double mkp::bisection(double M, double e, double epsilon, int max_it)
{

	int iterationsCount = 0;

	struct abPares
	{

		double a;
		double b;

		abPares(double a, double b) : a(a), b(b) {}

	};


	double a = 0;
	double b = 2 * PI;

	std::vector<abPares> seqOfSegments_ab;


	for (int i = 0; i < max_it; i++)
	{

		double E = (a + b) / 2;

		double f = KeplerEquation(M, e, E);

		if (f == 0)
		{
			std::cout << iterationsCount << "\n";
			return E;
		}

		double fa = KeplerEquation(M, e, a);
		double fb = KeplerEquation(M, e, b);

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

		iterationsCount++;

	}


	std::cout << iterationsCount << "\n";

	abPares retValue = seqOfSegments_ab[seqOfSegments_ab.size() - 1];

	return (retValue.a + retValue.b) / 2;

}


//Реализации метода золотого сечения
double mkp::goldensection(double M, double e, double epsilon, int max_it)
{

	int iterationscount = 0;

	struct abPares
	{

		double a;
		double b;

		abPares(double a, double b) : a(a), b(b) {}

	};


	double a = 0;
	double b = 2 * PI;

	std::vector<abPares> seqOfSegments_ab;


	for (int i = 0; i < max_it; i++)
	{

		double E = a + (b - a) / kappa;

		double f = KeplerEquation(M, e, E);

		if (f == 0)
		{

			std::cout << iterationscount << "\n";
			return E;
		}

		double fa = KeplerEquation(M, e, a);
		double fb = KeplerEquation(M, e, b);

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

		iterationscount++;

	}


	std::cout << iterationscount << "\n";


	abPares retValue = seqOfSegments_ab[seqOfSegments_ab.size() - 1];

	return retValue.a + (retValue.b - retValue.a) / kappa;

}


//Реализация метода итераций
double mkp::fixedpoint(double M, double e, double epsilon, int max_it)
{
	
	double E0 = M;
	double E = e * sin(E0) + M;


	int iterationsCount = 0;

	for (int i = 0; i < max_it; i++)
	{

		if (abs(E - E0) < epsilon)
		{
			std::cout << iterationsCount << "\n";
			return E;
		}

		E0 = E;
		E = e * sin(E0) + M;


		iterationsCount++;

	}

	std::cout<<iterationsCount << "\n";

	return E;

}


//Реализация метода Ньютона
double mkp::newton(double M, double e, double epsilon, int max_it)
{
	
	double E = M;

	int iterationsCount = 0;

	for (int i = 0; i < max_it; i++)
	{

		double f = KeplerEquation(M, e, E);
		double df = 1 - e * cos(E);


		double dE = -(f / df);

		E += dE;


		if (abs(dE) < epsilon)
		{
			std::cout << iterationsCount << "\n";
			return E;
		}

	}


	std::cout<<iterationsCount << "\n";

	return E;


}


//Функция нахождения всех параметров
void mkp::rootsOfKepEq(SpaceObject& object, double(*func)(double, double, double, int), double epsilon, int max_it)
{

	//Для построения графиков
	std::vector<std::pair<double, double>> M_graph;
	std::vector<std::pair<double, double>> E_graph;
	std::vector<std::pair<double, double>> v_graph;

	M_graph.emplace_back(0.f, 0.f);
	E_graph.emplace_back(0.f, 0.f);
	v_graph.emplace_back(0.f, 0.f);


	std::vector <std::pair<double, double>> r_graph;


	std::vector<std::pair<double, double>> Vr_graph;
	std::vector<std::pair<double, double>> Vn_graph;
	std::vector<std::pair<double, double>> V_graph;


	//Данные объекта
	double e = object.get_e();
	double a = object.get_a();
	double n = object.get_n();

	int nyu = object.get_nyu();
	int T = object.get_T();


	//Фокальный праметр
	double p = a * (1 - e * e);


	double r0 = p / (1 + e);
	r_graph.emplace_back(0.f, r0);


	double Vn0 = sqrt(nyu / p) * (1 + e);

	Vr_graph.emplace_back(0.f, 0.f);
	Vn_graph.emplace_back(0.f, Vn0);
	V_graph.emplace_back(0.f, Vn0);


	//Расчет
	double M, E, v, r, Vr, Vn, V;

	int iteration_delta = 50;

	for (int t = iteration_delta; t < T; t += iteration_delta)
	{

		M = n * t;

		M_graph.emplace_back((double)t, M);


		E = func(M, e, epsilon, max_it);

		E_graph.emplace_back((double)t, E);


		v = trueAnomaly(e, E);

		v_graph.emplace_back((double)t, v);


		r = p / (1 + e * cos(v));

		r_graph.emplace_back(t, r);


		Vr = sqrt(nyu / p) * e * sin(v);

		Vr_graph.emplace_back((double)t, Vr);


		Vn = sqrt(nyu / p) * (1 + e * cos(v));

		Vn_graph.emplace_back((double)t, Vn);


		V = sqrt(Vr * Vr + Vn * Vn);

		V_graph.emplace_back((double)t, V);


	}

	M = n * T;

	M_graph.emplace_back((double)T, M);


	E = func(M, e, epsilon, 10000);

	E_graph.emplace_back((double)T, E);


	v = trueAnomaly(e, E);

	v_graph.emplace_back((double)T, v);


	r = p / (1 + e * cos(v));

	r_graph.emplace_back(T, r);


	Vr = sqrt(nyu / p) * e * sin(v);

	Vr_graph.emplace_back((double)T, Vr);


	Vn = sqrt(nyu / p) * (1 + e * cos(v));

	Vn_graph.emplace_back((double)T, Vn);


	V = sqrt(Vr * Vr + Vn * Vn);

	V_graph.emplace_back((double)T, V);


	Gnuplot gp("\"C:\\Program Files\\gnuplot\\bin\\gnuplot.exe\"");


	std::ofstream outAnomalies("anomalies.txt", std::ios::out);
	std::ofstream outPosition("position.txt", std::ios::out);
	std::ofstream outVelocities("velocities.txt", std::ios::out);


	while (1)
	{

		std::cout << "1.Plot anomalies.\n";
		std::cout << "2.Plot position.\n";
		std::cout << "3.Plot velocities.\n";
		std::cout << "4.Output all in files.\n";
		std::cout << "(Other num).Exit.\n";

		int choice;

		std::cin >> choice;


		int iterationCount = T / iteration_delta + 2;


		std::string vstr, Mstr, Estr, rstr, Vrstr, Vnstr, Vstr, tstr;

		switch(choice)
		{
		case 1:

			gp << "set title 'График аномалий объекта'\n";

			gp << "set ytics pi\n";
			gp << "set format y '%.0Ppi'\n";

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

			
			if (outAnomalies.is_open() && outPosition.is_open() && outVelocities.is_open())
			{

				for (int i = 0; i < iterationCount; i++)
				{

					
					vstr = std::to_string(v_graph[i].second);
					Estr = std::to_string(E_graph[i].second);
					Mstr = std::to_string(M_graph[i].second);

					rstr = std::to_string(r_graph[i].second);

					Vrstr = std::to_string(Vr_graph[i].second);
					Vnstr = std::to_string(Vn_graph[i].second);
					Vstr = std::to_string(V_graph[i].second);

					tstr = std::to_string(i * iteration_delta);

					if(i * iteration_delta > T)
						tstr = std::to_string(T);

					outAnomalies << vstr << " " << Estr << " " << Mstr << " " << tstr << "\n";
					outPosition << rstr << " " << tstr << "\n";
					outVelocities << Vrstr << " " << Vnstr << " " << Vstr << " " << tstr << "\n";


				}

			}
			

			break;
		default:
			outAnomalies.close();
			outPosition.close();
			outVelocities.close();

			return;
		}


	}

	std::cout << "Success...\n\n";


}

