// #define USE_GENERAL_MESH

#include "../Common.hpp"
#include "../AffineTransform3.hpp"
#include "../Application.hpp"
#include "../Colors.hpp"
#include "../FilePath.hpp"
#include "../Math.hpp"
#include "../MatVec.hpp"
#include "../IPlugin.hpp"
#include "../Stopwatch.hpp"
#include "../Graphics/IRenderSystem.hpp"

#ifdef USE_GENERAL_MESH
#  include "../Graphics/GeneralMesh.hpp"
#  include "../Graphics/GraphicsAttributes.hpp"
#else
#  include "../Graphics/DisplayMesh.hpp"
#endif

#include "../Graphics/MeshGroup.hpp"
#include <cstdlib>
#include <cstdio>

#ifdef THEA_MAC
#  include <GLUT/glut.h>
#else
#  include <GL/glut.h>
#endif

using namespace std;
using namespace Thea;
using namespace Graphics;

int testDisplayMesh(int argc, char * argv[]);
void load(string const & filename);
int cleanup(int status);
void draw();
void reshape(int width, int height);
void normalKey(unsigned char key, int x, int y);
void specialKey(int key, int x, int y);
void init();

// Global variables
#ifdef USE_GENERAL_MESH
  typedef GeneralMesh< ColorAttribute<ColorRgba> > Mesh;
#else
  typedef DisplayMesh Mesh;
#endif

typedef MeshGroup<Mesh> MG;
MG mg;
GLfloat color[4] = { 0.7, 0.7, 0.7, 1.0 };
Vector3 translate(0, 0, 0);
GLfloat scale = 1;
IRenderSystem * render_system = nullptr;
GLfloat view_rotx = 20.0, view_roty = 30.0, view_rotz = 0.0;

int
main(int argc, char * argv[])
{
  int status = 0;
  try
  {
    status = testDisplayMesh(argc, argv);
  }
  THEA_STANDARD_CATCH_BLOCKS(return cleanup(-1);, ERROR, "%s", "An error occurred")

  // Check if tests passed
  if (status == 0) cout << "Test completed" << endl;

  return cleanup(0);
}

bool
initMesh(Mesh & mesh)
{
  cout << "Loaded mesh with " << mesh.numVertices() << " vertices, " << mesh.numTriangles() << " triangles and "
       << mesh.numQuads() << " quads" << endl;

// #define WIREFRAME

#ifdef USE_GENERAL_MESH

  if (mesh.hasVertexColors())
    cout << "Mesh has vertex colors" << endl;
  else
    cout << "Mesh does not have vertex colors" << endl;

  if (mesh.hasVertexTexCoords())
    cout << "Mesh has vertex texture coordinates" << endl;
  else
    cout << "Mesh does not have vertex texture coordinates" << endl;

  for (Mesh::VertexIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi)
    vi->attr().setColor(ColorRgb::random());

  mesh.setGpuBufferedRendering(true);

#ifdef WIREFRAME
  mesh.setGpuBufferedWireframe(true);
#endif

#else

  if (!mesh.hasNormals())
  {
    cout << "Computing normals for mesh " << mesh.getName() << endl;
    mesh.computeAveragedVertexNormals();
  }

#ifdef WIREFRAME
  mesh.setWireframeEnabled(true);
#endif

#endif

  return false;
}

void
load(string const & filename)
{
  mg.load(filename);

  cout << "Mesh group has " << mg.numMeshes() << " mesh(es)" << endl;
  mg.forEachMeshUntil(initMesh);

  mg.updateBounds();
  translate = -mg.getBounds().getCenter();
  Vector3 ext = mg.getBounds().getExtent();
  scale = 1.0f / max(ext.x(), max(ext.y(), ext.z()));
}

int
testDisplayMesh(int argc, char * argv[])
{
  // Load the mesh
  if (argc > 1)
    load(argv[1]);
  else
  {
    cerr << "Usage: " << argv[0] << " <mesh-filename> [--color=R,G,B]" << endl;
    return -1;
  }

  // Set the color, if specified
  if (argc > 2)
  {
    bool parsed = false;
    if (beginsWith(argv[2], "--color"))
    {
      float r = 0, g = 0, b = 0;
      if (sscanf(argv[2], "--color = %f , %f , %f", &r, &g, &b) == 3)
      {
        color[0] = r;
        color[1] = g;
        color[2] = b;
        parsed = true;
      }
    }

    if (!parsed)
    {
      cerr << "Could not parse argument " << argv[2] << endl;
      return -1;
    }
  }

  // Get the path containing the executable
  string bin_path = FilePath::parent(argv[0]);

  // Try to load the OpenGL plugin from the same parent directory as the executable
#ifdef THEA_DEBUG_BUILD
  string gl_plugin_path = FilePath::concat(bin_path, "../lib/libTheaPluginGLd");
#else
  string gl_plugin_path = FilePath::concat(bin_path, "../lib/libTheaPluginGL");
#endif

  cout << "Loading plugin: " << gl_plugin_path << endl;
  IPlugin * gl_plugin = Application::getPluginManager().load(gl_plugin_path);

  // Create a GL context via a GLUT window
  glutInit(&argc, argv);

#ifdef THEA_MAC
  glutInitDisplayString("double rgba depth samples");
#else
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
#endif

  glutInitWindowSize(640, 480);
  glutCreateWindow("Thea Mesh Test");

  // Start up the plugin
  gl_plugin->startup();

  // We should now have a GL rendersystem factory
  IRenderSystemFactory * factory = Application::getRenderSystemManager().getFactory("OpenGL");

  // Create a rendersystem
  render_system = factory->createRenderSystem("My OpenGL rendersystem");

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
  factory->destroyRenderSystem(render_system);

  // Cleanup and quit
  gl_plugin->shutdown();
  return 0;
}

int
cleanup(int status)
{
  Application::getPluginManager().unloadAllPlugins();
  return status;
}

#ifdef USE_GENERAL_MESH

// Just for some stress tests to see how fast we can push stuff to the GPU
bool
invalidateGpuBuffers(Mesh & mesh)
{
  mesh.invalidateGpuBuffers();
  return false;
}

#endif

void
draw()
{
// #define FRAME_TIMER
#ifdef FRAME_TIMER
  Stopwatch timer;
  timer.tick();
#endif

  render_system->setClearColor(ColorRgba::white().data());
  render_system->clear();

  Matrix3 rotx = Math::rotationAxisAngle(Vector3(1, 0, 0), Math::degreesToRadians(view_rotx));
  Matrix3 roty = Math::rotationAxisAngle(Vector3(0, 1, 0), Math::degreesToRadians(view_roty));
  Matrix3 rotz = Math::rotationAxisAngle(Vector3(0, 0, 1), Math::degreesToRadians(view_rotz));

  // std::cout << "rotx = " << toString(rotx) << std::endl;
  // std::cout << "roty = " << toString(roty) << std::endl;
  // std::cout << "rotz = " << toString(rotz) << std::endl;

  Matrix4 mat = (AffineTransform3(rotx * roty * rotz, Vector3::Zero())
               * AffineTransform3::scaling(scale)
               * AffineTransform3::translation(translate)).homogeneous();

  glShadeModel(GL_SMOOTH);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  static GLfloat const pos[4] = { 5.0, 5.0, 10.0, 0.0 };
  glLightfv(GL_LIGHT0, GL_POSITION, pos);

  RenderOptions render_opts = *RenderOptions::defaults();
#ifdef WIREFRAME
  render_opts.drawEdges() = true;
#endif

#ifdef USE_GENERAL_MESH
  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
  render_opts.sendColors() = true;
#else
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, color);
#endif

  render_system->setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW);
  render_system->pushMatrix();
    render_system->setMatrix(&asLvalue(Math::wrapMatrix(mat)));

#ifdef USE_GENERAL_MESH
    // Just a stress test to see how fast we can push stuff to the GPU
    // mg.forEachMeshUntil(invalidateGpuBuffers);
#endif

    mg.draw(render_system, &render_opts);

  render_system->popMatrix();

  glutSwapBuffers();

#ifdef FRAME_TIMER
  timer.tock();
  cout << "Frame rendered in " << timer.elapsedTime() << " seconds" << endl;
#endif
}

void
reshape(int width, int height)
{
  glViewport(0, 0, (GLint)width, (GLint)height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  GLfloat h = (GLfloat)height / (GLfloat)width;
  glOrtho(-1.0, 1.0, -h, h, 1.0, 100.0);
  glTranslatef(0.0, 0.0, -5.0);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glutPostRedisplay();
}

void
normalKey(unsigned char key, int x, int y)
{
  if (key == 27 || key == 'q' || key == 'Q')
  {
    Application::getPluginManager().unloadAllPlugins();
    exit(0);
  }
}

void
specialKey(int key, int x, int y)
{
  int shift = glutGetModifiers() & GLUT_ACTIVE_SHIFT;

  switch (key)
  {
    case GLUT_KEY_LEFT:
      if (shift) view_rotz -= 2.0;
      else       view_roty -= 2.0;
      glutPostRedisplay(); break;

    case GLUT_KEY_RIGHT:
      if (shift) view_rotz += 2.0;
      else       view_roty += 2.0;
      glutPostRedisplay(); break;

    case GLUT_KEY_UP:    view_rotx -= 2.0; glutPostRedisplay(); break;
    case GLUT_KEY_DOWN:  view_rotx += 2.0; glutPostRedisplay(); break;
  }
}

void
init()
{
  glDisable(GL_CULL_FACE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
}
