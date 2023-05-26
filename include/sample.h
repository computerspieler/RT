#ifndef _SAMPLE_H_
#define _SAMPLE_H_

#include "config.h"
#include "typedef.h"

typedef struct Sample Sample;
struct Sample {
	Float mean;
	long int count;
};

typedef struct TaskOutput TaskOutput;
struct TaskOutput {
	Float values[MAX_DEPTH];
	int triangles_associated[MAX_DEPTH];
	int length;
	bool went_out_of_bounds;
};

#endif
