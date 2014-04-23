#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "particle.h"

#ifdef OPENGL
#ifndef USEGLUT
#include <OpenGL/glext.h>
#include <OpenGL/gl3ext.h>
#include <GLFW/glfw3.h>
#include "visualization.h"

GLuint shader_program;

static const GLchar *vertex_shader_source= 
  {
    "attribute vec4 position;"
    "void main()\n"
    "{\n"
    "gl_Position = position;"
    "gl_PointSize = 2.;"
    "}\n"
  };

static const GLchar *fragment_shader_source=
  {
    "void main()\n"
    "{\n"
    "gl_FragColor = vec4(0.,1.,0.,1.);"
    "}\n"
  };


void initialize_glfw(){

  points = (float *)calloc(N*6+1,sizeof(float));

  glfwInit();

  /* display_init(argc, argv); */

  int width, height; width=800;height=650;
  window1 = glfwCreateWindow(width, height, "Rebound (Fork)", NULL, NULL);
  /* window2 = glfwCreateWindow(width, height/2, "Rebound (2)", NULL, NULL); */

  glfwGetFramebufferSize(window1, &width, &height);
  /* glViewport(0, height / 2, width / 2, height / 2); */

  glfwMakeContextCurrent(window1);
  /* glScissor(0, height / 2, width / 2, height / 2); */

  glfwSetWindowPos(window1,30,50);

  /* glfwSetWindowPos(window2,30,height+70); */
  /* glViewport(0, 0, width, height); */

  createshader();

  glfwSwapInterval(1);

  glClearColor(0.,0.,0.5,1.);	/* blue background */

  draw_particles();

}


static GLboolean printShaderInfoLog(GLuint shader, const char *str)
{
  GLint status;
  glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
  if (status == GL_FALSE) printf("shader compile ERROR\n");
  
  GLsizei bufSize;
  glGetShaderiv(shader, GL_INFO_LOG_LENGTH , &bufSize);
  
  if (bufSize > 1)
    {
      GLchar infoLog[1024];
      GLsizei length;
      glGetShaderInfoLog(shader, bufSize, &length, infoLog);
      printf("ERROR:%s\n",infoLog);
    }
  
  return (GLboolean)status;
}

static GLboolean printProgramInfoLog(GLuint program)
{
  GLint status;
  glGetProgramiv(program, GL_LINK_STATUS, &status);
  if (status == GL_FALSE) printf("Link error\n");
  
  GLsizei bufSize;
  glGetProgramiv(program, GL_INFO_LOG_LENGTH , &bufSize);
  
  if (bufSize > 1)
  {
    GLchar infoLog[1024];
    GLsizei length;
    glGetProgramInfoLog(program, bufSize, &length, infoLog);
    fprintf(stderr,"%s",infoLog);
  }
  return (GLboolean)status;
}

void createshader(){

  GLuint vshader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vshader,1,&vertex_shader_source, 0); /* third argument is pointer to shader source */
  glCompileShader(vshader);
  printShaderInfoLog(vshader, "vertex shader");

  GLuint fshader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fshader,1,&fragment_shader_source, 0); /* third argument is pointer to shader source */
  glCompileShader(fshader);
  printShaderInfoLog(fshader, "fragment shader");


  shader_program = glCreateProgram();
  glAttachShader(shader_program,vshader);
  glDeleteShader(vshader);
  glAttachShader(shader_program,fshader);
  glDeleteShader(fshader);

  glBindAttribLocation(shader_program, 0, "position");

  glLinkProgram(shader_program);
  printProgramInfoLog(shader_program);
}

GLuint myVBO_id;

void draw_particles_xy(){

  int i;
  for(i=0;i<N;i++){
    points[6*i] = (float)particles[i].x/boxsize_x*1.8;
    points[6*i+1] = (float)particles[i].y/boxsize_y*1.8;
    points[6*i+2] = 0.;
#ifndef COLLISIONS_NONE
    points[6*i+3] = (float)particles[i].r+1.;
#else
    points[6*i+3] = 0.;
#endif
    points[6*i+4] = 0.;
    points[6*i+5] = 0.;
  }


  glEnable(GL_SCISSOR_TEST);
  glViewport(30,200,350,350);
  glScissor(30,200,350,350);

  glBindBuffer(GL_ARRAY_BUFFER, myVBO_id);  
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  glEnableClientState(GL_VERTEX_ARRAY); 
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 6*sizeof(float), points); /* 3 float values for each vertex, offset 0*/
  glColorPointer(3, GL_FLOAT, 6*sizeof(float), &points[3]); /* 3 float values for each vertex, offset 0*/

  glDrawArrays(GL_POINTS,0,N);
  glDisableClientState(GL_COLOR_ARRAY); 
  glDisableClientState(GL_VERTEX_ARRAY); 

  glDisable(GL_SCISSOR_TEST);
}

void draw_particles_xz(){

  int i;
  for(i=0;i<N;i++){
    points[6*i] = (float)particles[i].x/boxsize_x*1.8;
    points[6*i+1] = (float)particles[i].z/boxsize_x*1.8;
    points[6*i+2] = 0.;
#ifndef COLLISIONS_NONE
    points[6*i+3] = (float)particles[i].r+1.;
#else
    points[6*i+3] = 0.;
#endif
    points[6*i+4] = 0.;
    points[6*i+5] = 0.;
  }

  glEnable(GL_SCISSOR_TEST);
  glViewport(30,0,350,200);
  glScissor(30,0,350,200);

  glBindBuffer(GL_ARRAY_BUFFER, myVBO_id);  
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  glEnableClientState(GL_VERTEX_ARRAY); 
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 6*sizeof(float), points); /* 3 float values for each vertex, offset 0*/
  glColorPointer(3, GL_FLOAT, 6*sizeof(float), &points[3]); /* 3 float values for each vertex, offset 0*/

  glDrawArrays(GL_POINTS,0,N);
  glDisableClientState(GL_COLOR_ARRAY); 
  glDisableClientState(GL_VERTEX_ARRAY); 

  glDisable(GL_SCISSOR_TEST);
}

void draw_particles_yz(){

  int i;
  for(i=0;i<N;i++){
    points[6*i] = (float)particles[i].z/boxsize_y*1.8;
    points[6*i+1] = (float)particles[i].y/boxsize_y*1.8;
    points[6*i+2] = 0.;
#ifndef COLLISIONS_NONE
    points[6*i+3] = (float)particles[i].r+1.;
#else
    points[6*i+3] = 0.;
#endif
    points[6*i+4] = 0.;
    points[6*i+5] = 0.;
  }

  glEnable(GL_SCISSOR_TEST);
  glViewport(350,200,200,350);
  glScissor(350,200,200,350);

  glBindBuffer(GL_ARRAY_BUFFER, myVBO_id);  
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  glEnableClientState(GL_VERTEX_ARRAY); 
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 6*sizeof(float), points); /* 3 float values for each vertex, offset 0*/
  glColorPointer(3, GL_FLOAT, 6*sizeof(float), &points[3]); /* 3 float values for each vertex, offset 0*/

  glDrawArrays(GL_POINTS,0,N);
  glDisableClientState(GL_COLOR_ARRAY); 
  glDisableClientState(GL_VERTEX_ARRAY); 

  glDisable(GL_SCISSOR_TEST);
}

void draw_particles_vxvy(){

  int i;
  for(i=0;i<N;i++){
    points[6*i] = (float)particles[i].vx/(boxsize_x*1e-4)*1.8;
    points[6*i+1] = ((float)particles[i].vy+1.5*0.00013143*(float)particles[i].x)/(boxsize_x*1e-4)*1.8;
    points[6*i+2] = 0.;
#ifndef COLLISIONS_NONE
    points[6*i+3] = (float)particles[i].r+1.;
#else
    points[6*i+3] = 0.;
#endif
    points[6*i+4] = 0.;
    points[6*i+5] = 0.;
  }

  glEnable(GL_SCISSOR_TEST);
  glViewport(450,200,350,350);
  glScissor(450,200,350,350);

  glBindBuffer(GL_ARRAY_BUFFER, myVBO_id);  
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

  glEnableClientState(GL_VERTEX_ARRAY); 
  glEnableClientState(GL_COLOR_ARRAY);

  glVertexPointer(3, GL_FLOAT, 6*sizeof(float), points); /* 3 float values for each vertex, offset 0*/
  glColorPointer(3, GL_FLOAT, 6*sizeof(float), &points[3]); /* 3 float values for each vertex, offset 0*/

  glDrawArrays(GL_POINTS,0,N);
  glDisableClientState(GL_COLOR_ARRAY); 
  glDisableClientState(GL_VERTEX_ARRAY); 

  glDisable(GL_SCISSOR_TEST);
}

void draw_particles(){

  draw_particles_xy();

  draw_particles_xz();

  draw_particles_yz();

  draw_particles_vxvy();
}

#endif // USEGLUT
#endif // OPENGL
