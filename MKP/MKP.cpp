#include "MKP.h"


namespace mkp
{


	SpaceObject::SpaceObject(int mass, float e, float a)
	{

		this->mass = mass;
		this->e = e;
		this->a = a;

		float nyu = G * mass;

		T = 2 * PI * sqrt(pow(a, 3) / nyu);

		n = 2 * PI / T;

	}


}
