#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "sphere.h"

Primitive *createSphere(SphereMetadata metadata)
{
	Primitive *output;
	SphereMetadata *output_metadata;

	output = (Primitive*) malloc(sizeof(Primitive));
	assert(output);

	output_metadata = (SphereMetadata*) malloc(sizeof(SphereMetadata));
	assert(output_metadata);

	memcpy(output_metadata, &metadata, sizeof(SphereMetadata));
	output->metadata = output_metadata;
	return output;
}