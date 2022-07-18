#include "typedef.h"

double lerp(double t, double x1, double x2)
{
	return (1 - t) * x1 + t * x2;
}
