#include "../../Common.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Algorithms/MeshTriangles.hpp"
#include "../../Graphics/Camera.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Plugins/GL/GLHeaders.hpp"
#include "../../Application.hpp"
#include "../../Ball3.hpp"
#include "../../ColorRGBA.hpp"
#include "../../CoordinateFrame3.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../Math.hpp"
#include "../../Matrix4.hpp"
#include "../../Plugin.hpp"
#include "../../Vector3.hpp"
#include <cstdio>
#include <fstream>
#include <sstream>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

typedef DisplayMesh Mesh;
typedef MeshGroup<Mesh> MG;

struct Model
{
  Model(bool convert_to_points_ = false) : convert_to_points(convert_to_points_), is_point_cloud(false) {}
  bool load(string const & path);
  Camera fitCamera(Matrix4 const & transform, int width, int height);
  bool render(ColorRGBA const & color);

  bool convert_to_points;
  MG mesh_group;
  bool is_point_cloud;
  TheaArray<Vector3> points;
};

enum PointUsage
{
  POINTS_NONE     = 0x0000,
  POINTS_PRIMARY  = 0x0001,
  POINTS_OVERLAY  = 0x0002,
  POINTS_ALL      = 0xFFFF,
};

RenderSystem * render_system = NULL;
TheaArray<string> model_paths;
TheaArray<Matrix4> transforms;
string out_path;
int out_width, out_height;
Vector3 view_dir(-1, -1, -1);
Vector3 view_up(0, 1, 0);
ColorRGBA primary_color(1.0f, 0.9f, 0.8f, 1.0f);
ColorRGBA background_color(1, 1, 1, 1);
int antialiasing_level = 1;
PointUsage show_points = POINTS_NONE;

bool parseArgs(int argc, char * argv[]);
bool loadPlugins(int argc, char * argv[]);
ColorRGBA getPaletteColor(long n);

int
main(int argc, char * argv[])
{
  if (!parseArgs(argc, argv))
    return -1;

  if (!loadPlugins(argc, argv))
    return -1;

  // Load the mesh
  Model model(show_points & POINTS_PRIMARY);
  if (!model.load(model_paths[0]))
    return -1;

  // Set up framebuffer for offscreen drawing
  int buffer_width   =  antialiasing_level * out_width;
  int buffer_height  =  antialiasing_level * out_height;
  Texture * color_tex = NULL;
  try
  {
    Texture::Options tex_opts = Texture::Options::defaults();
    tex_opts.interpolateMode = Texture::InterpolateMode::NEAREST_NO_MIPMAP;
    color_tex = render_system->createTexture("Color", buffer_width, buffer_height, 1, Texture::Format::RGBA8(),
                                             Texture::Dimension::DIM_2D, tex_opts);
    if (!color_tex)
    {
      THEA_ERROR << "Could not create color buffer";
      return -1;
    }

    Texture * depth_tex = render_system->createTexture("Depth", buffer_width, buffer_height, 1, Texture::Format::DEPTH24(),
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
    Camera camera = model.fitCamera(transforms[0], buffer_width, buffer_height);

    // Render the mesh to the offscreen framebuffer
    render_system->pushFramebuffer();
      render_system->setFramebuffer(fb);

      render_system->pushDepthFlags();
      render_system->pushColorFlags();
      render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
      render_system->setMatrixMode(RenderSystem::MatrixMode::PROJECTION); render_system->pushMatrix();

        render_system->setCamera(camera);
        render_system->setDepthTest(RenderSystem::DepthTest::LESS);
        render_system->setDepthWrite(true);
        render_system->setColorWrite(true, true, true, true);

        render_system->setColorClearValue(background_color);
        render_system->clear();

        // Draw primary model
        render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
          render_system->multMatrix(transforms[0]);
          model.render(primary_color);
        render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();

        // Draw overlay models
        for (array_size_t i = 1; i < model_paths.size(); ++i)
        {
          model.convert_to_points = (show_points & POINTS_OVERLAY);
          if (!model.load(model_paths[i]))
            return -1;

          ColorRGBA overlay_color = getPaletteColor((long)i - 1);
          overlay_color.a() = (model.is_point_cloud ? 1.0f : 0.5f);

          render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
            render_system->multMatrix(transforms[i]);
            model.render(overlay_color);
          render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
        }

      render_system->setMatrixMode(RenderSystem::MatrixMode::PROJECTION); render_system->popMatrix();
      render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
      render_system->popColorFlags();
      render_system->popDepthFlags();

    render_system->popFramebuffer();
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "Could not render mesh")

  // Save the rendered image
  try
  {
    Image image(Image::Type::RGB_8U, buffer_width, buffer_height);
    color_tex->getImage(image);

    if (antialiasing_level > 1 && !image.rescale(out_width, out_height, Image::Filter::BICUBIC))
    {
      THEA_ERROR << "Could not rescale image to output dimensions";
      return -1;
    }

    image.save(out_path);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not save rendered image to '%s'", out_path.c_str())

  return 0;
}

bool
usage()
{
  string app_path = FilePath::objectName(Application::programPath());

  THEA_CONSOLE << "";
  THEA_CONSOLE << "Usage: " << app_path << " [OPTIONS] <mesh> <output-image> <width> <height>";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << "  -o <overlay-shape>    (may be mesh or point set)";
  THEA_CONSOLE << "  -p <scope>            (show 'none' | 'primary' | 'overlay' | 'all' shapes";
  THEA_CONSOLE << "                         as points)";
  THEA_CONSOLE << "  -t <transform>        (row-major comma-separated 3x4 or 4x4 matrix,";
  THEA_CONSOLE << "                         applied to all subsequent shapes)";
  THEA_CONSOLE << "  -v <viewing-dir>      (3 chars, one for each coordinate, each one of";
  THEA_CONSOLE << "                         +, - or 0)";
  THEA_CONSOLE << "  -u <up-dir>           (x, y or z, optionally preceded by + or -)";
  THEA_CONSOLE << "  -c <argb>             (mesh color)";
  THEA_CONSOLE << "  -b <argb>             (background color)";
  THEA_CONSOLE << "  -a N                  (enable NxN antialiasing: 2 is normal, 4 is very";
  THEA_CONSOLE << "                         high quality)";
  THEA_CONSOLE << "";

  return false;
}

bool
parseTransform(string const & s, Matrix4 & m)
{
  TheaArray<string> fields;
  long nfields = stringSplit(trimWhitespace(s), ",;:[({<>})] \t\n\r\f", fields, true);
  if (nfields != 12 && nfields != 16)
  {
    THEA_ERROR << "Could not read row-major comma-separated matrix from '" << s << '\'';
    return false;
  }

  m = Matrix4::identity();
  for (int i = 0; i < nfields; ++i)
  {
    istringstream field_in(fields[i]);
    if (!(field_in >> m(i / 4, i % 4)))
    {
      THEA_ERROR << "Could nor parse matrix entry '" << fields[i] << '\'';
      return false;
    }
  }

  return true;
}

bool
parseViewDirection(string const & s, Vector3 & dir, Vector3 & up)
{
  if (s.length() != 3)
  {
    THEA_ERROR << "Viewing direction string must have exactly 3 characters, one for each coordinate";
    return -1;
  }

  if (s == "000")
  {
    THEA_ERROR << "View direction is zero vector";
    return -1;
  }

  for (int i = 0; i < 3; ++i)
  {
    switch (s[i])
    {
      case '+': dir[i] =  1; break;
      case '-': dir[i] = -1; break;
      case '0': dir[i] =  0; break;
      default:
        THEA_ERROR << "Invalid view direction string '" << s << '\'';
        return -1;
    }
  }

  if (s == "0-0")
    up = -Vector3::unitZ();
  else if (s == "0+0")
    up = Vector3::unitZ();
  else
    up = Vector3::unitY();

  return true;
}

bool
parseViewUp(string const & s, Vector3 & up)
{
  if (s.length() != 1 && s.length() != 2)
  {
    THEA_ERROR << "Up direction must be 'x', 'y' or 'z', optionally preceded by + or -";
    return false;
  }

  char c = s[0];
  bool neg = false;
  if (c == '+' || c == '-')
  {
    if (s.length() < 2)
    {
      THEA_ERROR << "Up direction must be 'x', 'y' or 'z', optionally preceded by + or -";
      return false;
    }

    if (c == '-')
      neg = true;

    c = s[1];
  }

  if (c == 'x' || c == 'X')
    up = Vector3::unitX();
  else if (c == 'y' || c == 'Y')
    up = Vector3::unitY();
  else if (c == 'z' || c == 'Z')
    up = Vector3::unitZ();
  else
  {
    THEA_ERROR << "Up direction must be 'x', 'y' or 'z', optionally preceded by + or -";
    return false;
  }

  if (neg)
    up = -up;

  return true;
}

bool
parseColor(string const & s, ColorRGBA & c)
{
  std::stringstream ss;
  ss << std::hex << s;

  uint32 argb;
  if (!(ss >> argb))
  {
    THEA_ERROR << "Could not parse color '" << s << '\'';
    return false;
  }

  c = ColorRGBA::fromARGB(argb);

  if (trimWhitespace(s).length() <= 6)  // alpha channel not specified
    c.a() = 1.0;

  return true;
}

bool
parseArgs(int argc, char * argv[])
{
  if (argc < 5)
    return usage();

  Matrix4 current_transform = Matrix4::identity();
  bool has_up = false;

  argv++;
  argc--;
  int pos = 0;

  while (argc > 0)
  {
    string arg = *argv;
    argv++; argc--;

    if (arg.length() <= 0)
      continue;

    if (arg[0] == '-')
    {
      if (arg.length() != 2)
        return usage();

      switch (arg[1])
      {
        case 'o':
        {
          if (argc < 1) { THEA_ERROR << "-o: Overlay path not specified"; return false; }
          model_paths.push_back(*argv);
          transforms.push_back(current_transform);
          argv++; argc--; break;
        }

        case 'p':
        {
          if (argc < 1) { THEA_ERROR << "-p: Shapes to show as points not specified"; return false; }
          string scope = toLower(*argv);
          if (scope == "all")
            show_points = POINTS_ALL;
          else if (scope == "primary")
            show_points = POINTS_PRIMARY;
          else if (scope == "overlay")
            show_points = POINTS_OVERLAY;
          else if (scope == "none")
            show_points = POINTS_NONE;
          else
          {
            THEA_ERROR << "Unknown set of shapes to show as points: '" << scope << '\'';
            return false;
          }

          argv++; argc--; break;
        }

        case 't':
        {
          if (argc < 1) { THEA_ERROR << "-t: Transform not specified"; return false; }
          if (!parseTransform(*argv, current_transform)) return false;
          argv++; argc--; break;
        }

        case 'v':
        {
          if (argc < 1) { THEA_ERROR << "-v: View direction not specified"; return false; }
          Vector3 up;
          if (!parseViewDirection(*argv, view_dir, up)) return false;
          if (!has_up) view_up = up;
          argv++; argc--; break;
        }

        case 'u':
        {
          if (argc < 1) { THEA_ERROR << "-u: Up direction not specified"; return false; }
          if (!parseViewUp(*argv, view_up)) return false;
          has_up = true;
          argv++; argc--; break;
        }

        case 'c':
        {
          if (argc < 1) { THEA_ERROR << "-c: Mesh color not specified"; return false; }
          if (!parseColor(*argv, primary_color)) return false;
          argv++; argc--; break;
        }

        case 'b':
        {
          if (argc < 1) { THEA_ERROR << "-b: Background color not specified"; return false; }
          if (!parseColor(*argv, background_color)) return false;
          argv++; argc--; break;
        }

        case 'a':
        {
          if (argc < 1) { THEA_ERROR << "-a: Antialiasing level not specified"; return false; }
          if (sscanf(*argv, " %d", &antialiasing_level) != 1)
          {
            THEA_ERROR << "Could not parse antialiasing level '" << *argv << '\'';
            return false;
          }
          if (antialiasing_level < 1)
          {
            THEA_ERROR << "Invalid antialiasing level " << antialiasing_level;
            return false;
          }
          argv++; argc--; break;
        }
      }
    }
    else
    {
      pos++;

      switch (pos)
      {
        case 1:
        {
          model_paths.insert(model_paths.begin(), arg);
          transforms.insert(transforms.begin(), current_transform);
          break;
        }

        case 2:
        {
          out_path = arg;
          break;
        }

        case 3:
        {
          if (argc < 1) { THEA_ERROR << "Width not followed by height"; return false; }

          if (sscanf(arg.c_str(), "%d", &out_width) != 1 || sscanf(*argv, "%d", &out_height) != 1
           || out_width <= 0 || out_height <= 0)
          {
            THEA_ERROR << "Could not parse output image dimensions: " << arg << " x " << *argv;
            return false;
          }

          argv++; argc--; pos++;
          break;
        }

        default:
        {
          THEA_ERROR << "Too many positional arguments";
          return false;
        }
      }
    }
  }

  if (pos < 4)
  {
    THEA_ERROR << "Too few positional arguments";
    return usage();
  }

  return true;
}

ColorRGBA
getPaletteColor(long n)
{
  static ColorRGBA PALETTE[] = {
    ColorRGBA::fromARGB(0xFFFF0000),
    ColorRGBA::fromARGB(0xFF00FF00),
    ColorRGBA::fromARGB(0xFF0000FF),
    ColorRGBA::fromARGB(0xFF00FFFF),
    ColorRGBA::fromARGB(0xFFFF00FF),
    ColorRGBA::fromARGB(0xFFFFFF00),
    ColorRGBA::fromARGB(0xFF800000),
    ColorRGBA::fromARGB(0xFF008000),
    ColorRGBA::fromARGB(0xFF000080),
    ColorRGBA::fromARGB(0xFF008080),
    ColorRGBA::fromARGB(0xFF800080),
    ColorRGBA::fromARGB(0xFF808000),
  };

  return PALETTE[n % (sizeof(PALETTE) / sizeof(ColorRGBA))];
}

bool
averageNormals(Mesh & mesh)
{
  if (!mesh.hasNormals())
    mesh.computeAveragedVertexNormals();

  return false;
}

bool
Model::load(string const & path)
{
  mesh_group.clear();
  points.clear();
  is_point_cloud = false;

  if (endsWith(toLower(path), ".pts"))
  {
    ifstream in(path.c_str());
    if (!in)
    {
      THEA_ERROR << "Could not open points file '" << path << '\'';
      return false;
    }

    string line;
    double x, y, z;
    while (getline(in, line))
    {
      if (trimWhitespace(line).empty())
        continue;

      istringstream line_in(line);
      if (!(line_in >> x >> y >> z))
      {
        THEA_ERROR << "Could not read point " << points.size() << " from '" << path << '\'';
        return false;
      }

      points.push_back(Vector3((Real)x, (Real)y, (Real)z));
    }

    is_point_cloud = true;

    THEA_CONSOLE << "Read " << points.size() << " points from '" << path << '\'';
  }
  else
  {
    try
    {
      mesh_group.load(path);
    }
    THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "Could not load model from '%s'", path.c_str())

    if (convert_to_points)
    {
      MeshSampler<Mesh> sampler(mesh_group);

      double out_scale = min(out_width, out_height);
      long max_samples = (long)ceil(10 * out_scale);
      long npoints = sampler.sampleEvenlyByArea(max_samples, points);

      THEA_CONSOLE << "Sampled " << npoints << " points from '" << path << '\'';

      mesh_group.clear();
      is_point_cloud = true;
    }
    else
      mesh_group.forEachMeshUntil(averageNormals);
  }

  return true;
}

class FarthestPoint
{
  public:
    FarthestPoint(Vector3 const & center_, Matrix4 const & transform_)
    : center(center_), transform(transform_), max_sqdist(0) {}

    bool operator()(Mesh const & mesh)
    {
      Mesh::VertexArray const & vertices = mesh.getVertices();
      for (array_size_t i = 0; i < vertices.size(); ++i)
      {
        Real sqdist = (transform * vertices[i] - center).squaredLength();
        if (sqdist > max_sqdist)
          max_sqdist = sqdist;
      }

      return false;
    }

    Real getFarthestDistance() const { return sqrt(max_sqdist); }

  private:
    Vector3 center;
    Matrix4 const & transform;
    Real max_sqdist;
};

Ball3
modelBSphere(Model const & model, Matrix4 const & transform)
{
  double sum_x = 0, sum_y = 0, sum_z = 0;
  double sum_w = 0;

  if (model.is_point_cloud)
  {
    for (array_size_t i = 0; i < model.points.size(); ++i)
    {
      sum_x += model.points[i][0];
      sum_y += model.points[i][1];
      sum_z += model.points[i][2];
      sum_w += 1;
    }
  }
  else
  {
    MeshTriangles<Mesh> tris;
    tris.add(const_cast<MG &>(model.mesh_group));

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
  }

  Vector3 center(0, 0, 0);
  if (sum_w > 0)
  {
    center[0] = (Real)(sum_x / sum_w);
    center[1] = (Real)(sum_y / sum_w);
    center[2] = (Real)(sum_z / sum_w);
  }

  center = transform * center;

  Real radius = 0;
  if (model.is_point_cloud)
  {
    for (array_size_t i = 0; i < model.points.size(); ++i)
    {
      Real sqdist = (transform * model.points[i] - center).squaredLength();
      if (i == 0 || sqdist > radius)
        radius = sqdist;
    }

    radius = sqrt(radius);
  }
  else
  {
    FarthestPoint fp(center, transform);
    model.mesh_group.forEachMeshUntil(&fp);
    radius = fp.getFarthestDistance();
  }

  return Ball3(center, radius);
}

Camera
Model::fitCamera(Matrix4 const & transform, int width, int height)
{
  // Orientation
  Ball3 bsphere = modelBSphere(*this, transform);
  Vector3 center = bsphere.getCenter();
  Real diameter = bsphere.getDiameter();

  Real camera_separation = diameter > 1.0e-10f ? 2.1 * diameter : 1.0e-10f;
  Vector3 dir = view_dir.unit();
  Vector3 up  = view_up.unit();

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
initPointShader(Shader & shader)
{
  static string const VERTEX_SHADER =
"void main()\n"
"{\n"
"  gl_Position = ftransform();\n"
"  gl_FrontColor = gl_Color;\n"
"  gl_BackColor = gl_Color;\n"
"}\n";

  static string const FRAGMENT_SHADER =
"void main()\n"
"{\n"
"  gl_FragColor = gl_Color;\n"
"}\n";

  try
  {
    shader.attachModuleFromString(Shader::ModuleType::VERTEX, VERTEX_SHADER.c_str());
    shader.attachModuleFromString(Shader::ModuleType::FRAGMENT, FRAGMENT_SHADER.c_str());
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "Could not attach point shader module")

  return true;
}

bool
initMeshShader(Shader & shader)
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
"  gl_FragColor = vec4(ambt_color + lamb_color, gl_Color.a);\n"
"}\n";

  try
  {
    shader.attachModuleFromString(Shader::ModuleType::VERTEX, VERTEX_SHADER.c_str());
    shader.attachModuleFromString(Shader::ModuleType::FRAGMENT, FRAGMENT_SHADER.c_str());
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "Could not attach mesh shader module")

  shader.setUniform("light_dir", Vector3(-1, -1, -2));
  shader.setUniform("light_color", ColorRGB(1, 1, 1));
  shader.setUniform("ambient_color", ColorRGB(1, 1, 1));
  shader.setUniform("two_sided", 1.0f);
  shader.setUniform("material", Vector4(0.2f, 0.6f, 0.2f, 25));

  return true;
}

bool
Model::render(ColorRGBA const & color)
{
  render_system->pushShader();
  render_system->pushShapeFlags();
  render_system->pushColorFlags();

  render_system->setColor(color);

  bool has_transparency = (color.a() < 1);
  if (has_transparency)
  {
    // Enable alpha-blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }

  if (is_point_cloud)
  {
    // Initialize the shader
    static Shader * point_shader = NULL;
    if (!point_shader)
    {
      point_shader = render_system->createShader("Point shader");
      if (!point_shader)
      {
        THEA_ERROR << "Could not create point shader";
        return -1;
      }

      if (!initPointShader(*point_shader))
      {
        THEA_ERROR << "Could not initialize point shader";
        return -1;
      }
    }

    render_system->setShader(point_shader);

    render_system->setPointSize(antialiasing_level);
    render_system->beginPrimitive(RenderSystem::Primitive::POINTS);

      for (array_size_t i = 0; i < points.size(); ++i)
        render_system->sendVertex(points[i]);

    render_system->endPrimitive();
  }
  else
  {
    // Initialize the shader
    static Shader * mesh_shader = NULL;
    if (!mesh_shader)
    {
      mesh_shader = render_system->createShader("Mesh shader");
      if (!mesh_shader)
      {
        THEA_ERROR << "Could not create mesh shader";
        return -1;
      }

      if (!initMeshShader(*mesh_shader))
      {
        THEA_ERROR << "Could not initialize mesh shader";
        return -1;
      }
    }

    render_system->setShader(mesh_shader);

    if (has_transparency)
    {
      // First back faces
      render_system->setCullFace(RenderSystem::CullFace::FRONT);
      mesh_group.draw(*render_system);

      // Then front faces
      render_system->setCullFace(RenderSystem::CullFace::BACK);
      mesh_group.draw(*render_system);
    }
    else
      mesh_group.draw(*render_system);
  }

  render_system->popColorFlags();
  render_system->popShapeFlags();
  render_system->popShader();

  return true;
}

bool
loadPlugins(int argc, char * argv[])
{
  string app_path = FileSystem::resolve(Application::programPath());
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
