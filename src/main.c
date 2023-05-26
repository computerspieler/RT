#include <CL/cl.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <sys/time.h>

#include <GL/glut.h>

#include "bvhtree.h"
#include "cl_handler.h"
#include "material.h"
#include "scene_loader.h"
#include "scene.h"
#include "transform.h"
#include "host.h"
#include "sample.h"
#include "kernel/main.cl.h"

long long get_current_time()
{
	struct timeval te;
	gettimeofday(&te, NULL);
	return te.tv_sec * 1000LL + te.tv_usec / 1000;
}

vec3 rotation, position;
float *viewer_vertices;
float *viewer_colors;
HostContext ctx;

void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glClearColor(.4, .58, .94, 1);

	glLoadIdentity();
	glRotatef(rotation.x, 1, 0, 0);
	glRotatef(rotation.y, 0, 1, 0);
	glRotatef(rotation.z, 0, 0, 1);
	glTranslatef(position.x, position.y, position.z);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, viewer_vertices);
	
	glEnableClientState(GL_COLOR_ARRAY);
	glColorPointer(3, GL_FLOAT, 0, viewer_colors);
	
	glDrawArrays(GL_TRIANGLES, 0, 3 * ctx.scene_meta.triangles_count);
	
	glDisableClientState(GL_VERTEX_ARRAY);
	glFlush();
	glutSwapBuffers();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat) w / (GLfloat) h, 0.1, 100.0);
	glMatrixMode(GL_MODELVIEW);
}

const double speed = .25;
void keyboard(unsigned char key, int x, int y)
{
	switch (key) {
	case 'q': case 'Q': position.x += speed; break;
	case 'd': case 'D': position.x -= speed; break;
	case  8 :			position.y += speed; break;
	case ' ':			position.y -= speed; break;
	case 'z': case 'Z': position.z += speed; break;
	case 's': case 'S': position.z -= speed; break;

	case 27: exit(0); return;
	default: return;
	}
	glutPostRedisplay();
}

void keyboard_special(int key, int x, int y)
{
	switch (key) {
	case 101: rotation.x +=  speed * 2. * M_PI; break;
	case 103: rotation.x += -speed * 2. * M_PI; break;
	case 100: rotation.y += -speed * 2. * M_PI; break;
	case 102: rotation.y +=  speed * 2. * M_PI; break;
	case 104: rotation.z += -speed * 2. * M_PI; break;
	case 105: rotation.z +=  speed * 2. * M_PI; break;

	default: return;
	}
	glutPostRedisplay();
}

Scene obj;

Float MRSE;
Float mean_time;
Float total_area;
long long first_start;
long int lost_photons;

size_t origin[3] = {0};
size_t region[3] = {SAMPLE_PER_ITERATION, 1, 1};

Sample *density_per_tri;
TaskOutput *task_output;
OpenCL_GeneralContext cl_gen;
OpenCL_ProgramContext cl_prg;

long int count = 0;

void update(int iteration_count)
{
	long long start = get_current_time();
	ctx.seed = (uint) rand();
	opencl_run(&cl_gen, &cl_prg, origin, region);
	opencl_postrun(&cl_gen, &cl_prg, origin, region);
	
	for(int i = 0; i < SAMPLE_PER_ITERATION; i ++) {
		for(int j = 0; j < task_output[i].length; j ++) {
			int triangle_id = task_output[i].triangles_associated[j];

			density_per_tri[triangle_id].count ++;
			Float count_f = (Float) density_per_tri[triangle_id].count;
			
			density_per_tri[triangle_id].mean *= (count_f - 1) / count_f;
			density_per_tri[triangle_id].mean += task_output[i].values[j] / count_f;
		}
		count += task_output[i].length-1;
		if(!task_output[i].went_out_of_bounds)
			count ++;
	}

	mean_time *= iteration_count - 1;
	mean_time += get_current_time() - start;
	mean_time /= iteration_count;

	printf(
		"Milliseconds spent: %lld, Mean time per sampling: %f, iteration count: %d\r",
		get_current_time() - start,
		mean_time,
		iteration_count
	);

	if(iteration_count % 10 == 0) {
		MRSE = 0;
		printf("\n");

		Float max = 0;
		Float tmp = 0;

		for(size_t i = 0; i < ctx.scene_meta.triangles_count; i ++) {
			tmp  = density_per_tri[i].mean;
			tmp *= (Float) density_per_tri[i].count / (Float) count;
			MRSE += tmp;
			
			for(int j = 0; j < 9; j ++)
				viewer_colors[9 * i + j] = density_per_tri[i].mean / max;
		}
		
		MRSE /= total_area;
		printf("MRSE: " FLOAT_FMT "\n", MRSE);

		glutPostRedisplay();
	}

	glutTimerFunc(2, update, iteration_count + 1);
}

int createBuffers()
{
	cl_int err;

	viewer_vertices = malloc(ctx.scene_meta.triangles_count * 9 * sizeof(float));
	for(size_t i = 0; i < ctx.scene_meta.triangles_count; i ++) {
		for(int j = 0; j < 3; j ++) {
			vec3 p = obj.vertices[obj.triangles[i].vertices[j]];

			viewer_vertices[9 * i + 3 * j    ] = p.x;
			viewer_vertices[9 * i + 3 * j + 1] = p.y;
			viewer_vertices[9 * i + 3 * j + 2] = p.z;
		}
	}

	viewer_colors   = malloc(ctx.scene_meta.triangles_count * 9 * sizeof(float));
	bzero(viewer_colors, ctx.scene_meta.triangles_count * 9 * sizeof(float));

	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.materials, ctx.scene_meta.materials_count * sizeof(Material));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.vertices, ctx.scene_meta.vertices_count * sizeof(vec3));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.vertices_tex, ctx.scene_meta.vertices_tex_count * sizeof(vec2));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.triangles, ctx.scene_meta.triangles_count * sizeof(Triangle));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, ctx.scene_meta.tree.nodes, ctx.scene_meta.tree.nodes_count * sizeof(BVHNode));

	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, &ctx, sizeof(HostContext));
	
	TRY(opencl_add_input_buffer, exit,
		&cl_gen, &cl_prg, obj.lights, ctx.scene_meta.lights_count * sizeof(Light));

	density_per_tri = (Sample*) malloc(ctx.scene_meta.triangles_count * sizeof(Sample));
	bzero(density_per_tri, ctx.scene_meta.triangles_count * sizeof(Sample));

	task_output = (TaskOutput*) malloc(SAMPLE_PER_ITERATION * sizeof(TaskOutput));
	TRY(opencl_add_input_output_buffer, exit,
		&cl_gen, &cl_prg, task_output, SAMPLE_PER_ITERATION * sizeof(TaskOutput));
		
exit:
	return err;
}

void freeBuffers()
{
	free(task_output);
	free(density_per_tri);
	free(viewer_colors);
	free(viewer_vertices);

	scene_delete(&obj, &ctx.scene_meta);
}

int main(int argc, char* argv[])
{
	FILE *f;
	cl_int err;

	if(argc != 2) {
		printf("Usage: program [SCENE]\n");
		return -1;
	}

	srand(time(NULL));
	
	// Open the context
	err = opencl_init_general_context(&cl_gen);
	if(err != CL_SUCCESS)
		goto exit;

	// Load the kernel
	err = opencl_init_program_context(&cl_prg);
	if(err != CL_SUCCESS)
		goto exit;

	err = opencl_add_program_source(&cl_prg, kernel_main_cl_src, sizeof(kernel_main_cl_src));
	if(err != CL_SUCCESS)
		goto exit;

	// Load the scene
	f = fopen(argv[1], "r");
	if(loadFromObj(f, &obj, &ctx.scene_meta, false))
		goto exit;
	fclose(f);

	printObjectInfo(ctx.scene_meta);


	char arguments[1024];
	printf("Compile the program\n");
	sprintf(
		arguments, "-DMAX_TREE_DEPTH=%u",
		BVHTree_depth(ctx.scene_meta.tree.nodes, ctx.scene_meta.tree.root)
	);

	TRY(opencl_build_program, exit,
		&cl_gen, &cl_prg, "compute_sample", arguments);
	
	TRY(createBuffers, exit);

	printf("Run the program\n");
	opencl_prerun(&cl_prg);

	first_start = get_current_time();
	mean_time = 0.f;
	total_area = 0;
	for(size_t i = 0; i < ctx.scene_meta.triangles_count; i ++)
		total_area += obj.triangles[i].area;
	printf("Area: " FLOAT_FMT "\n", total_area);
	rotation = position = VEC3(0, 0, 0);
	
	glutInitWindowSize(200, 200);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutCreateWindow(argv[0]);

	glClearColor(.4, .58, .94, 1);
	glShadeModel(GL_FLAT);
	glEnable(GL_DEPTH_TEST);

	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutDisplayFunc(display);
	glutSpecialFunc(keyboard_special);
	glutTimerFunc(0, update, 1);
	glutMainLoop();

	freeBuffers();

exit:
	printf("Error code : %i\n", err);
	printf("Exit !\n");

	return 0;
}
