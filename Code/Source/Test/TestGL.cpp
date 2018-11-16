#include "../Common.hpp"
#include "../Application.hpp"
#include "../FilePath.hpp"
#include "../Plugin.hpp"
#include "../System.hpp"
#include "../Graphics/RenderSystem.hpp"
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <ctime>

#ifdef THEA_MAC
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

using namespace std;
using namespace Thea;
using namespace Graphics;

void testGL(int argc, char * argv[]);
int cleanup(int status);
void draw();
void reshape(int width, int height);
void normalKey(unsigned char key, int x, int y);
void specialKey(int key, int x, int y);
void init();

int
main(int argc, char * argv[])
{
  // Do the OpenGL tests
  try
  {
    testGL(argc, argv);
  }
  THEA_STANDARD_CATCH_BLOCKS(return cleanup(-1);, ERROR, "%s", "An error occurred")

  // Hooray, all tests passed
  cout << "Test completed" << endl;

  return cleanup(0);
}

void
testGL(int argc, char * argv[])
{
  // Get the path containing the executable
  string bin_path = FilePath::parent(argv[0]);

  // Try to load the OpenGL plugin from the same parent directory as the executable
#ifdef THEA_DEBUG_BUILD
  string plugin_path = FilePath::concat(bin_path, "../lib/libTheaPluginGLd");
#else
  string plugin_path = FilePath::concat(bin_path, "../lib/libTheaPluginGL");
#endif

  cout << "Loading plugin: " << plugin_path << endl;
  Plugin * gl_plugin = Application::getPluginManager().load(plugin_path);

  // Create a GL context via a GLUT window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(640, 480);
  glutCreateWindow("Thea OpenGL Test");

  // Start up the plugin
  gl_plugin->startup();

  // We should now have a GL rendersystem factory
  RenderSystemFactory * factory = Application::getRenderSystemManager().getFactory("OpenGL");

  // Create a rendersystem
  RenderSystem * rs = factory->createRenderSystem("My OpenGL rendersystem");

  cout << "GL_RENDERER  =  " << (char *)glGetString(GL_RENDERER) << endl;
  cout << "GL_VERSION   =  " << (char *)glGetString(GL_VERSION)  << endl;
  cout << "GL_VENDOR    =  " << (char *)glGetString(GL_VENDOR)   << endl;

  init();

  glutDisplayFunc(draw);
  glutIdleFunc(draw);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(normalKey);
  glutSpecialFunc(specialKey);
  glutMainLoop();

  // Destroy the rendersystem
  factory->destroyRenderSystem(rs);

  // Cleanup and quit
  gl_plugin->shutdown();
}

int
cleanup(int status)
{
  Application::getPluginManager().unloadAllPlugins();
  return status;
}

// The following functions are lifted from the Mesa source code (http://www.mesa3d.org), courtesy Brian Paul.

#ifndef M_PI
#  define M_PI 3.14159265
#endif

GLfloat view_rotx = 20.0, view_roty = 30.0, view_rotz = 0.0;
GLint gear1, gear2, gear3;
GLfloat angle = 0.0;

void
gear(GLfloat inner_radius, GLfloat outer_radius, GLfloat width, GLint teeth, GLfloat tooth_depth)
{
   GLint i;
   GLfloat r0, r1, r2;
   GLfloat angle, da;
   GLfloat u, v, len;

   r0 = inner_radius;
   r1 = outer_radius - tooth_depth / 2.0;
   r2 = outer_radius + tooth_depth / 2.0;

   da = 2.0 * M_PI / teeth / 4.0;

   glShadeModel(GL_FLAT);

   glNormal3f(0.0, 0.0, 1.0);

   /* draw front face */
   glBegin(GL_QUAD_STRIP);
   for (i = 0; i <= teeth; i++) {
      angle = i * 2.0 * M_PI / teeth;
      glVertex3f(r0 * cos(angle), r0 * sin(angle), width * 0.5);
      glVertex3f(r1 * cos(angle), r1 * sin(angle), width * 0.5);
      if (i < teeth) {
         glVertex3f(r0 * cos(angle), r0 * sin(angle), width * 0.5);
         glVertex3f(r1 * cos(angle + 3 * da), r1 * sin(angle + 3 * da),
                    width * 0.5);
      }
   }
   glEnd();

   /* draw front sides of teeth */
   glBegin(GL_QUADS);
   da = 2.0 * M_PI / teeth / 4.0;
   for (i = 0; i < teeth; i++) {
      angle = i * 2.0 * M_PI / teeth;

      glVertex3f(r1 * cos(angle), r1 * sin(angle), width * 0.5);
      glVertex3f(r2 * cos(angle + da), r2 * sin(angle + da), width * 0.5);
      glVertex3f(r2 * cos(angle + 2 * da), r2 * sin(angle + 2 * da),
                 width * 0.5);
      glVertex3f(r1 * cos(angle + 3 * da), r1 * sin(angle + 3 * da),
                 width * 0.5);
   }
   glEnd();

   glNormal3f(0.0, 0.0, -1.0);

   /* draw back face */
   glBegin(GL_QUAD_STRIP);
   for (i = 0; i <= teeth; i++) {
      angle = i * 2.0 * M_PI / teeth;
      glVertex3f(r1 * cos(angle), r1 * sin(angle), -width * 0.5);
      glVertex3f(r0 * cos(angle), r0 * sin(angle), -width * 0.5);
      if (i < teeth) {
         glVertex3f(r1 * cos(angle + 3 * da), r1 * sin(angle + 3 * da),
                    -width * 0.5);
         glVertex3f(r0 * cos(angle), r0 * sin(angle), -width * 0.5);
      }
   }
   glEnd();

   /* draw back sides of teeth */
   glBegin(GL_QUADS);
   da = 2.0 * M_PI / teeth / 4.0;
   for (i = 0; i < teeth; i++) {
      angle = i * 2.0 * M_PI / teeth;

      glVertex3f(r1 * cos(angle + 3 * da), r1 * sin(angle + 3 * da),
                 -width * 0.5);
      glVertex3f(r2 * cos(angle + 2 * da), r2 * sin(angle + 2 * da),
                 -width * 0.5);
      glVertex3f(r2 * cos(angle + da), r2 * sin(angle + da), -width * 0.5);
      glVertex3f(r1 * cos(angle), r1 * sin(angle), -width * 0.5);
   }
   glEnd();

   /* draw outward faces of teeth */
   glBegin(GL_QUAD_STRIP);
   for (i = 0; i < teeth; i++) {
      angle = i * 2.0 * M_PI / teeth;

      glVertex3f(r1 * cos(angle), r1 * sin(angle), width * 0.5);
      glVertex3f(r1 * cos(angle), r1 * sin(angle), -width * 0.5);
      u = r2 * cos(angle + da) - r1 * cos(angle);
      v = r2 * sin(angle + da) - r1 * sin(angle);
      len = sqrt(u * u + v * v);
      u /= len;
      v /= len;
      glNormal3f(v, -u, 0.0);
      glVertex3f(r2 * cos(angle + da), r2 * sin(angle + da), width * 0.5);
      glVertex3f(r2 * cos(angle + da), r2 * sin(angle + da), -width * 0.5);
      glNormal3f(cos(angle), sin(angle), 0.0);
      glVertex3f(r2 * cos(angle + 2 * da), r2 * sin(angle + 2 * da),
                 width * 0.5);
      glVertex3f(r2 * cos(angle + 2 * da), r2 * sin(angle + 2 * da),
                 -width * 0.5);
      u = r1 * cos(angle + 3 * da) - r2 * cos(angle + 2 * da);
      v = r1 * sin(angle + 3 * da) - r2 * sin(angle + 2 * da);
      glNormal3f(v, -u, 0.0);
      glVertex3f(r1 * cos(angle + 3 * da), r1 * sin(angle + 3 * da),
                 width * 0.5);
      glVertex3f(r1 * cos(angle + 3 * da), r1 * sin(angle + 3 * da),
                 -width * 0.5);
      glNormal3f(cos(angle), sin(angle), 0.0);
   }

   glVertex3f(r1 * cos(0), r1 * sin(0), width * 0.5);
   glVertex3f(r1 * cos(0), r1 * sin(0), -width * 0.5);

   glEnd();

   glShadeModel(GL_SMOOTH);

   /* draw inside radius cylinder */
   glBegin(GL_QUAD_STRIP);
   for (i = 0; i <= teeth; i++) {
      angle = i * 2.0 * M_PI / teeth;
      glNormal3f(-cos(angle), -sin(angle), 0.0);
      glVertex3f(r0 * cos(angle), r0 * sin(angle), -width * 0.5);
      glVertex3f(r0 * cos(angle), r0 * sin(angle), width * 0.5);
   }
   glEnd();
}

void
draw()
{
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

   glPushMatrix();
   glRotatef(view_rotx, 1.0, 0.0, 0.0);
   glRotatef(view_roty, 0.0, 1.0, 0.0);
   glRotatef(view_rotz, 0.0, 0.0, 1.0);

   glPushMatrix();
   glTranslatef(-3.0, -2.0, 0.0);
   glRotatef(angle, 0.0, 0.0, 1.0);
   glCallList(gear1);
   glPopMatrix();

   glPushMatrix();
   glTranslatef(3.1, -2.0, 0.0);
   glRotatef(-2.0 * angle - 9.0, 0.0, 0.0, 1.0);
   glCallList(gear2);
   glPopMatrix();

   glPushMatrix();
   glTranslatef(-3.1, 4.2, 0.0);
   glRotatef(-2.0 * angle - 25.0, 0.0, 0.0, 1.0);
   glCallList(gear3);
   glPopMatrix();

   glPopMatrix();

   glutSwapBuffers();

   // next frame
   angle += 2.0;

   // calc framerate
   {
      static double t0 = -1;
      static int frames = 0;
      double t = System::time();

      if (t0 < 0)
         t0 = t;

      frames++;

      if (t - t0 >= 5.0)
      {
         double seconds = t - t0;
         double fps = frames / seconds;
         printf("%d frames in %3.1lf seconds = %6.3lf FPS\n", frames, seconds, fps);
         t0 = t;
         frames = 0;
      }
   }
}

void
reshape(int width, int height)
{
   GLfloat h = (GLfloat) height / (GLfloat) width;

   glViewport(0, 0, (GLint) width, (GLint) height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glFrustum(-1.0, 1.0, -h, h, 5.0, 60.0);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glTranslatef(0.0, 0.0, -40.0);
   glutPostRedisplay();
}

void
normalKey(unsigned char key, int x, int y)
{
  if (key == 27) exit(0);
}

void
specialKey(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_LEFT:  view_roty += 5.0; glutPostRedisplay(); break;
    case GLUT_KEY_RIGHT: view_roty -= 5.0; glutPostRedisplay(); break;
    case GLUT_KEY_UP:    view_rotx += 5.0; glutPostRedisplay(); break;
    case GLUT_KEY_DOWN:  view_rotx -= 5.0; glutPostRedisplay(); break;
  }
}

void
init()
{
   static GLfloat pos[4] = { 5.0, 5.0, 10.0, 0.0 };
   static GLfloat red[4] = { 0.8, 0.1, 0.0, 1.0 };
   static GLfloat green[4] = { 0.0, 0.8, 0.2, 1.0 };
   static GLfloat blue[4] = { 0.2, 0.2, 1.0, 1.0 };

   glLightfv(GL_LIGHT0, GL_POSITION, pos);
   glEnable(GL_CULL_FACE);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);
   glEnable(GL_DEPTH_TEST);

   /* make the gears */
   gear1 = glGenLists(1);
   glNewList(gear1, GL_COMPILE);
   glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, red);
   gear(1.0, 4.0, 1.0, 20, 0.7);
   glEndList();

   gear2 = glGenLists(1);
   glNewList(gear2, GL_COMPILE);
   glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, green);
   gear(0.5, 2.0, 2.0, 10, 0.7);
   glEndList();

   gear3 = glGenLists(1);
   glNewList(gear3, GL_COMPILE);
   glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, blue);
   gear(1.3, 2.0, 0.5, 10, 0.7);
   glEndList();

   glEnable(GL_NORMALIZE);
}
