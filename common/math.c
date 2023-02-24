#include "typedef.h"

Float lerp(Float t, Float x1, Float x2)
{
	return (1 - t) * x1 + t * x2;
}
