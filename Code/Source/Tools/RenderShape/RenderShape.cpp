#include "../../Common.hpp"
#include "../../Algorithms/MeshTriangles.hpp"
#include "../../Graphics/Camera.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Application.hpp"
#include "../../Ball3.hpp"
#include "../../CoordinateFrame3.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../Math.hpp"
#include "../../Plugin.hpp"
#include "../../Vector3.hpp"
#include <cstdio>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

typedef DisplayMesh Mesh;
typedef MeshGroup<Mesh> MG;

RenderSystem * render_system = NULL;

bool loadPlugins(int argc, char * argv[]);
bool averageNormals(Mesh & mesh);
bool initShader(Shader & shader);
Ball3 modelBSphere(MG const & mg);
Camera fitCameraToModel(MG const & mg, int width, int height);
bool render(MG const & mg, bool is_overlay);

int
main(int argc, char * argv[])
{
  if (argc < 5)
  {
    THEA_CONSOLE << "Usage: " << argv[0] << " <mesh> <output-image> <width> <height> [<overlay0> [<overlay1> ...]]";
    return -1;
  }

  if (!loadPlugins(argc, argv))
    return -1;

  string mesh_path = argv[1];
  string out_path = argv[2];

  int width, height;
  if (sscanf(argv[3], "%d", &width) != 1 || sscanf(argv[4], "%d", &height) != 1 || width <= 0 || height <= 0)
  {
    THEA_ERROR << "Could not parse output image dimensions: " << argv[3] << " x " << argv[4];
    return -1;
  }

  // Load the mesh
  MG mg;
  try
  {
    mg.load(mesh_path);
    mg.forEachMeshUntil(averageNormals);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not load mesh from '%s'", mesh_path.c_str())

  // Set up framebuffer for offscreen drawing
  Texture * color_tex = NULL;
  try
  {
    Texture::Options tex_opts = Texture::Options::defaults();
    tex_opts.interpolateMode = Texture::InterpolateMode::NEAREST_NO_MIPMAP;
    color_tex = render_system->createTexture("Color", width, height, 1, Texture::Format::RGBA8(), Texture::Dimension::DIM_2D,
                                             tex_opts);
    if (!color_tex)
    {
      THEA_ERROR << "Could not create color buffer";
      return -1;
    }

    Texture * depth_tex = render_system->createTexture("Depth", width, height, 1, Texture::Format::DEPTH24(),
                                                       Texture::Dimension::DIM_2D, tex_opts);
    if (!depth_tex)
    {
      THEA_ERROR << "Could not create depth buffer";
      return -1;
    }

    Framebuffer * fb = render_system->createFramebuffer("Framebuffer");
    if (!fb)
    {
      THEA_ERROR << "Could not create offscreen framebuffer";
      return -1;
    }

    fb->attach(Framebuffer::AttachmentPoint::COLOR_0, color_tex);
    fb->attach(Framebuffer::AttachmentPoint::DEPTH,   depth_tex);

    // Initialize the camera
    Camera camera = fitCameraToModel(mg, width, height);

    // Initialize the shader
    Shader * shader = render_system->createShader("Shader");
    if (!shader)
    {
      THEA_ERROR << "Could not create shader";
      return -1;
    }

    if (!initShader(*shader))
    {
      THEA_ERROR << "Could not initialize shader";
      return -1;
    }

    // Render the mesh to the offscreen framebuffer
    render_system->pushFramebuffer();
    render_system->pushShader();
    render_system->pushDepthFlags();
    render_system->pushColorFlags();
    render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
    render_system->setMatrixMode(RenderSystem::MatrixMode::PROJECTION); render_system->pushMatrix();

      render_system->setFramebuffer(fb);
      render_system->setShader(shader);
      render_system->setCamera(camera);
      render_system->setDepthTest(RenderSystem::DepthTest::LESS);
      render_system->setDepthWrite(true);
      render_system->setColorWrite(true, true, true, true);

      render_system->setColorClearValue(ColorRGB(1, 1, 1));
      render_system->clear();

      render(mg, false);

    render_system->setMatrixMode(RenderSystem::MatrixMode::PROJECTION); render_system->popMatrix();
    render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
    render_system->popColorFlags();
    render_system->popDepthFlags();
    render_system->popShader();
    render_system->popFramebuffer();
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "Could not render mesh")

  // Save the rendered image;
  try
  {
    Image image(Image::Type::RGB_8U, width, height);
    color_tex->getImage(image);
    image.save(out_path);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not save rendered image to '%s'", out_path.c_str())

  return 0;
}

bool
averageNormals(Mesh & mesh)
{
  if (!mesh.hasNormals())
    mesh.computeAveragedVertexNormals();

  return false;
}

class FarthestPoint
{
  public:
    FarthestPoint(Vector3 const & center_) : center(center_), max_sqdist(0) {}

    bool operator()(Mesh const & mesh)
    {
      Mesh::VertexArray const & vertices = mesh.getVertices();
      for (array_size_t i = 0; i < vertices.size(); ++i)
      {
        Real sqdist = (vertices[i] - center).squaredLength();
        if (sqdist > max_sqdist)
          max_sqdist = sqdist;
      }

      return false;
    }

    Real getFarthestDistance() const { return sqrt(max_sqdist); }

  private:
    Vector3 center;
    Real max_sqdist;
};

Ball3
modelBSphere(MG const & mg)
{
  MeshTriangles<Mesh> tris;
  tris.add(const_cast<MG &>(mg));

  double sum_x = 0, sum_y = 0, sum_z = 0;
  double sum_w = 0;

  MeshTriangles<Mesh>::TriangleArray const & tri_array = tris.getTriangles();
  for (array_size_t i = 0; i < tri_array.size(); ++i)
  {
    Vector3 c = tri_array[i].getCentroid();
    Real area = tri_array[i].getArea();

    sum_x += (area * c[0]);
    sum_y += (area * c[1]);
    sum_z += (area * c[2]);

    sum_w += area;
  }

  Vector3 center(0, 0, 0);
  if (sum_w > 0)
  {
    center[0] = (Real)(sum_x / sum_w);
    center[1] = (Real)(sum_y / sum_w);
    center[2] = (Real)(sum_z / sum_w);
  }

  FarthestPoint fp(center);
  mg.forEachMeshUntil(&fp);
  Real radius = fp.getFarthestDistance();

  return Ball3(center, radius);
}

Camera
fitCameraToModel(MG const & mg, int width, int height)
{
  // Orientation
  Ball3 bsphere = modelBSphere(mg);
  Vector3 center = bsphere.getCenter();
  Real diameter = bsphere.getDiameter();

  Real camera_separation = diameter > 1.0e-3f ? 2 * diameter : 1.0e-3f;
  Vector3 dir = -Vector3::unitZ();
  Vector3 up  = Vector3::unitY();

  CoordinateFrame3 cframe = CoordinateFrame3::fromViewFrame(center - camera_separation * dir,  // eye
                                                            center,                            // look-at
                                                            up);                               // up

  // Projection
  static Real const HALF_WIDTH = 0.5;
  Real hw = 0, hh = 0;
  if (height > width)
  {
    Real aspect_ratio = height / (Real)width;
    hw = HALF_WIDTH;
    hh = aspect_ratio * HALF_WIDTH;
  }
  else
  {
    Real aspect_ratio = width / (Real)height;
    hw = aspect_ratio * HALF_WIDTH;
    hh = HALF_WIDTH;
  }

  Real near_dist = 0.01 * camera_separation;
  Real far_dist  = 5 * camera_separation;

  hw *= (0.5f * near_dist);
  hh *= (0.5f * near_dist);

  return Camera(cframe,
                Camera::ProjectionType::PERSPECTIVE, -hw, hw, -hh, hh, near_dist, far_dist, Camera::ProjectedYDirection::UP);
}

bool
render(MG const & mg, bool is_overlay)
{
  if (is_overlay)
  {
    return false;
  }
  else
  {
    mg.draw(*render_system);
  }

  return true;
}

bool
loadPlugins(int argc, char * argv[])
{
  string app_path = FileSystem::resolve(argv[0]);
  string plugin_dir = FilePath::concat(FilePath::parent(FilePath::parent(app_path)), "lib");

  // Try to load the OpenGL plugin
#ifdef THEA_DEBUG_BUILD

#ifdef THEA_WINDOWS
  string debug_plugin_path    =  FilePath::concat(plugin_dir, "TheaPluginGLd");
  string release_plugin_path  =  FilePath::concat(plugin_dir, "TheaPluginGL");
#else
  string debug_plugin_path    =  FilePath::concat(plugin_dir, "libTheaPluginGLd");
  string release_plugin_path  =  FilePath::concat(plugin_dir, "libTheaPluginGL");
#endif

#ifdef THEA_WINDOWS
  string debug_plugin_path_ext = debug_plugin_path + ".dll";
#elif THEA_OSX
  string debug_plugin_path_ext = debug_plugin_path + ".dylib";
#else
  string debug_plugin_path_ext = debug_plugin_path + ".so";
#endif

  string plugin_path = FileSystem::exists(debug_plugin_path_ext) ? debug_plugin_path : release_plugin_path;
#else

#ifdef THEA_WINDOWS
  string plugin_path = FilePath::concat(plugin_dir, "TheaPluginGL");
#else
  string plugin_path = FilePath::concat(plugin_dir, "libTheaPluginGL");
#endif

#endif

  Plugin * gl_plugin = Application::getPluginManager().load(plugin_path);
  if (!gl_plugin)
  {
    THEA_ERROR << "Could not load OpenGL plugin: " << plugin_path;
    return false;
  }

  gl_plugin->startup();

  RenderSystemFactory * render_system_factory = Application::getRenderSystemManager().getFactory("OpenGL");
  if (!render_system_factory)
  {
    THEA_ERROR << "Could not get OpenGL rendersystem factory";
    return false;
  }

  render_system = render_system_factory->createRenderSystem("OpenGL RenderSystem");
  if (!render_system)
  {
    THEA_ERROR << "Could not create OpenGL rendersystem";
    return false;
  }

  return true;
}

bool
initShader(Shader & shader)
{
  static string const VERTEX_SHADER =
"varying vec3 position;  // position in camera space\n"
"varying vec3 normal;  // normal in camera space\n"
"\n"
"void main()\n"
"{\n"
"  gl_Position = ftransform();\n"
"\n"
"  position = vec3(gl_ModelViewMatrix * gl_Vertex);  // assume rigid transform, so we can drop w\n"
"  normal = gl_NormalMatrix * gl_Normal;\n"
"\n"
"  gl_FrontColor = gl_Color;\n"
"  gl_BackColor = gl_Color;\n"
"}\n";

  static string const FRAGMENT_SHADER =
"uniform vec3 ambient_color;\n"
"uniform vec3 light_dir;  // must be specified in camera space, pointing towards object\n"
"uniform vec3 light_color;\n"
"uniform vec4 material;  // [ka, kl, <ignored>, <ignored>]\n"
"uniform float two_sided;\n"
"\n"
"varying vec3 position;  // position in camera space\n"
"varying vec3 normal;  // normal in camera space\n"
"\n"
"void main()\n"
"{\n"
"  vec3 N = normalize(normal);\n"
"  vec3 L = normalize(light_dir);\n"
"\n"
"  vec3 ambt_color = material[0] * gl_Color.rgb * ambient_color;\n"
"\n"
"  float NdL = -dot(N, L);\n"
"  vec3 lamb_color = (NdL >= -two_sided) ? material[1] * abs(NdL) * gl_Color.rgb * light_color : vec3(0.0, 0.0, 0.0);\n"
"\n"
"  gl_FragColor = vec4(ambt_color + lamb_color, 1.0);\n"
"}\n";

  try
  {
    shader.attachModuleFromString(Shader::ModuleType::VERTEX, VERTEX_SHADER.c_str());
    shader.attachModuleFromString(Shader::ModuleType::FRAGMENT, FRAGMENT_SHADER.c_str());
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "Could not attach shader module")

  shader.setUniform("light_dir", Vector3(-1, -1, -2));
  shader.setUniform("light_color", ColorRGB(1, 1, 1));
  shader.setUniform("ambient_color", ColorRGB(1, 0.8f, 0.7f));
  shader.setUniform("two_sided", 1.0f);
  shader.setUniform("material", Vector4(0.2f, 0.6f, 0.2f, 25));

  return true;
}
