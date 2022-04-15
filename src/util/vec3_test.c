#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "vec3.h"
#include "vec3_test.h"

int vec3_check_equality(Vec3 v1, Vec3 v2)
{
	return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
}

void vec3_test_create()
{
	Vec3 valid_response = (Vec3) {.x = 1, .y = 2, .z = 3};
	assert(vec3_check_equality(vec3_create(1, 2, 3), valid_response));
}

void vec3_test_add()
{
	Vec3 v1 = vec3_create(1, 2, 3);
	Vec3 v2 = vec3_create(2, 3, 4);
	Vec3 valid_response = (Vec3) {.x = 3, .y = 5, .z = 7};

	assert(vec3_check_equality(vec3_add(v1, v2), valid_response));
}

void vec3_test_diff()
{
	Vec3 v1 = vec3_create(1, 2, 3);
	Vec3 v2 = vec3_create(2, 4, 4);
	Vec3 valid_response = (Vec3) {.x = -1, .y = -2, .z = -1};

	assert(vec3_check_equality(vec3_diff(v1, v2), valid_response));
}

void vec3_test_cross()
{
	Vec3 v1 = vec3_create(1, 2, 3);
	Vec3 v2 = vec3_create(4, 5, 6);
	Vec3 valid_response = (Vec3) {.x = -3, .y = 6, .z = -3};

	assert(vec3_check_equality(vec3_cross(v1, v2), valid_response));
}

void vec3_test_dot()
{
	Vec3 v1 = vec3_create(1, 2, 3);
	Vec3 v2 = vec3_create(2, 4, 4);

	assert(vec3_dot(v1, v2) == 22);
}

void vec3_test_smul()
{
	Vec3 v = vec3_create(1, 2, 3);
	Vec3 valid_response = (Vec3) {.x = 3, .y = 6, .z = 9};

	assert(vec3_check_equality(vec3_smul(3, v), valid_response));
}

void vec3_test_norm()
{
	Vec3 v = vec3_create(1, 1, 1);
	assert(vec3_norm(v) == sqrt(3));
}

void vec3_test_normalize()
{
	Vec3 v1 = vec3_create(-4, 0, 0);
	Vec3 null_v = vec3_create(0, 0, 0);

	assert(vec3_check_equality(vec3_normalize(v1), vec3_create(-1, 0, 0)));
	assert(vec3_check_equality(vec3_normalize(null_v), vec3_create(1, 0, 0)));
}

void test()
{
	vec3_test_create();
	vec3_test_add();
	vec3_test_diff();
	vec3_test_cross();
	vec3_test_dot();
	vec3_test_smul();
	vec3_test_norm();
	vec3_test_normalize();
}
