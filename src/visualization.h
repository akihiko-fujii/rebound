#include <stdio.h>
#include <GLFW/glfw3.h>
#include <OpenGL/glext.h>
#include <OpenGL/gl3ext.h>
#include <stdlib.h>
#include <math.h>

void createshader();

extern GLuint shader_program;
extern float *points;
extern GLFWwindow *window1, *window2;

void draw_particles();

void initialize_glfw();

