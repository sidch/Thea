#include "ShapeRenderer.hpp"
#include "../../Algorithms/BestFitSphere3.hpp"
#include "../../Algorithms/KDTreeN.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Algorithms/MeshTriangles.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Graphics/Camera.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshCodec.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Plugins/GL/GLHeaders.hpp"
#include "../../Application.hpp"
#include "../../Ball3.hpp"
#include "../../BoundedSortedArrayN.hpp"
#include "../../ColorRGBA.hpp"
#include "../../CoordinateFrame3.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../Image.hpp"
#include "../../Math.hpp"
#include "../../Matrix4.hpp"
#include "../../Memory.hpp"
#include "../../Plugin.hpp"
#include "../../Random.hpp"
#include "../../UnorderedMap.hpp"
#include "../../Vector3.hpp"
#include "../../Vector4.hpp"
#include <boost/functional/hash.hpp>
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <utility>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

// #define DRAW_EDGES

typedef DisplayMesh Mesh;
typedef MeshGroup<Mesh> MG;
typedef std::pair<Mesh const *, long> FaceRef;
typedef TheaUnorderedMap<FaceRef, uint32> FaceIndexMap;

struct View
{
  Vector3 dir;
  bool has_eye;
  Vector3 eye;
  Vector3 up;

  View() : dir(-1, -1, -1), has_eye(false), up(0, 1, 0) {}
};

struct Model
{
  Model(bool convert_to_points_ = false) : convert_to_points(convert_to_points_), is_point_cloud(false) {}
  Camera fitCamera(Matrix4 const & transform, View const & view, Real zoom, int width, int height);

  bool convert_to_points;
  MG mesh_group;
  bool is_point_cloud;
  TheaArray<Vector3> points;
  TheaArray<ColorRGBA> point_colors;
};

enum PointUsage
{
  POINTS_NONE     = 0x0000,
  POINTS_PRIMARY  = 0x0001,
  POINTS_OVERLAY  = 0x0002,
  POINTS_ALL      = 0xFFFF,
};

class ShapeRendererImpl
{
  private:
    static AtomicInt32 has_render_system;
    static RenderSystemFactory * render_system_factory;
    static RenderSystem * render_system;
    static Shader * point_shader;
    static Shader * mesh_shader;
    static Shader * face_id_shader;
    static Texture * matcap_tex;
    static Texture * tex3d;

    TheaArray<string> model_paths;
    TheaArray<Matrix4> transforms;
    float zoom;
    string out_path;
    int out_width, out_height;
    TheaArray<View> views;
    bool has_up;
    Vector3 view_up;
    float point_size;
    bool color_by_id;
    bool color_by_label;
    bool color_by_features;
    string labels_path;
    TheaArray<int> labels;
    string features_path;
    bool accentuate_features;
    bool color_by_tex3d;
    string tex3d_path;
    string selected_mesh;
    Vector4 material;
    string matcap_path;
    ColorRGBA primary_color;
    ColorRGBA background_color;
    int antialiasing_level;
    PointUsage show_points;
    bool flat;
    bool save_depth;
    bool print_camera;
    int palette_shift;

    bool loadPlugins(int argc, char ** argv);
    bool parseArgs(int argc, char ** argv);
    bool usage();
    bool parseTransform(string const & s, Matrix4 & m);
    bool parseViewDiscrete(string const & s, View & view, bool silent = false);
    bool parseViewContinuous(string const & s, View & view, bool silent = false);
    bool parseViewFile(string const & path);
    bool parseViewUp(string const & s, Vector3 & up);
    bool parseColor(string const & s, ColorRGBA & c);
    int parseVector(string const & str, Vector4 & v);
    void resetArgs();
    bool loadModel(Model & model, string const & path);
    bool loadLabels(Model & model, FaceIndexMap const * tri_ids, FaceIndexMap const * quad_ids);
    bool loadFeatures(Model & model);
    bool renderModel(Model const & model, ColorRGBA const & color);
    void colorizeMeshSelection(MG & mg, uint32 parent_id);
    ColorRGBA getPaletteColor(long n) const;
    ColorRGBA getLabelColor(long label) const;

    friend struct FaceColorizer;

  public:
    ShapeRendererImpl(int argc, char * argv[]);  // just loads plugins and initializes variables
    ~ShapeRendererImpl();  // shuts down plugins

    int exec(string const & cmdline);
    int exec(int argc, char ** argv);

}; // class ShapeRendererImpl

ShapeRenderer::ShapeRenderer(int argc, char * argv[])
{
  impl = new ShapeRendererImpl(argc, argv);
}

ShapeRenderer::~ShapeRenderer()
{
  delete impl;
}

int
ShapeRenderer::exec(string const & cmdline)
{
  return impl->exec(cmdline);
}

int
ShapeRenderer::exec(int argc, char ** argv)
{
  return impl->exec(argc, argv);
}

AtomicInt32 ShapeRendererImpl::has_render_system(0);
RenderSystemFactory * ShapeRendererImpl::render_system_factory = NULL;
RenderSystem * ShapeRendererImpl::render_system = NULL;
Shader * ShapeRendererImpl::point_shader = NULL;
Shader * ShapeRendererImpl::mesh_shader = NULL;
Shader * ShapeRendererImpl::face_id_shader = NULL;
Texture * ShapeRendererImpl::matcap_tex = NULL;
Texture * ShapeRendererImpl::tex3d = NULL;

ShapeRendererImpl::ShapeRendererImpl(int argc, char * argv[])
{
  resetArgs();

  if (has_render_system.compareAndSet(0, 1) == 0)
  {
    if (!loadPlugins(argc, argv))
      throw Error("Could not load plugins");
  }
}

ShapeRendererImpl::~ShapeRendererImpl()
{
  if (has_render_system.compareAndSet(1, 0) == 1)
  {
    if (render_system_factory)
      render_system_factory->destroyRenderSystem(render_system);
  }

  Application::getPluginManager().unloadAllPlugins();
}

void
ShapeRendererImpl::resetArgs()
{
  model_paths.clear();
  transforms.clear();
  zoom = 1.0f;
  out_path = "";
  out_width = out_height = -1;
  views.clear();
  has_up = false;
  view_up = Vector3(0, 1, 0);
  point_size = 1.0f;
  color_by_id = false;
  color_by_label = false;
  color_by_features = false;
  labels_path = "";
  labels.clear();
  features_path = "";
  accentuate_features = false;
  color_by_tex3d = false;
  tex3d_path = "";
  selected_mesh = "";
  material = Vector4(0.3f, 0.7f, 0.2f, 25);
  matcap_path = "";
  primary_color = ColorRGBA(1.0f, 0.9f, 0.8f, 1.0f);
  background_color = ColorRGBA(1, 1, 1, 1);
  antialiasing_level = 1;
  show_points = POINTS_NONE;
  flat = false;
  save_depth = false;
  print_camera = false;
  palette_shift = 0;
}

int
ShapeRendererImpl::exec(string const & cmdline)  // cmdline should not include program name
{
  TheaArray<string> args;
  stringSplit(cmdline, " \t\n\f\r", args, true);

  TheaArray<char *> argv(args.size() + 1);

  string prog_name = FilePath::objectName(Application::programPath());
  argv[0] = new char[prog_name.length() + 1];
  strcpy(argv[0], prog_name.c_str());

  for (size_t i = 0; i < args.size(); ++i)
  {
    argv[i + 1] = new char[args[i].length() + 1];
    strcpy(argv[i + 1], args[i].c_str());
  }

  int status = exec((int)argv.size(), &argv[0]);

  for (size_t i = 0; i < argv.size(); ++i)
    delete [] argv[i];

  return status;
}

// Guaranteed to return a value between 0 and 2^32 - 1
uint32
labelHash(string const & label)
{
  boost::hash<string> hasher;
  string rev = label; reverse(rev.begin(), rev.end());
  return (uint32)((hasher(label) + hasher(rev)) & 0x7FFFFFFF);
}

ColorRGBA
ShapeRendererImpl::getPaletteColor(long n) const
{
  static ColorRGBA PALETTE[] = {
    ColorRGBA::fromARGB(0xFFFF66FF),
    ColorRGBA::fromARGB(0xFFFFFF66),
    ColorRGBA::fromARGB(0xFF66FFFF),
    ColorRGBA::fromARGB(0xFF3399FF),
    ColorRGBA::fromARGB(0xFFFF9933),
    ColorRGBA::fromARGB(0xFFFF3399),
    ColorRGBA::fromARGB(0xFF99FF33),
    ColorRGBA::fromARGB(0xFF9933FF),
    ColorRGBA::fromARGB(0xFFFFAAAA),
    ColorRGBA::fromARGB(0xFFAAFFAA),
    ColorRGBA::fromARGB(0xFFAAAAFF),
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
    ColorRGBA::fromARGB(0xFFA0A0A0),
  };

  static long const PALETTE_SIZE = (long)(sizeof(PALETTE) / sizeof(ColorRGBA));
  return PALETTE[abs((n % (PALETTE_SIZE + palette_shift)) % PALETTE_SIZE)];
}

ColorRGBA
ShapeRendererImpl::getLabelColor(long label) const
{
  return getPaletteColor(label);
}

int
ShapeRendererImpl::exec(int argc, char ** argv)
{
  if (!parseArgs(argc, argv))
    return -1;

  // Load the mesh
  Model model(show_points & POINTS_PRIMARY);
  if (!loadModel(model, model_paths[0]))
    return -1;

  // Set up framebuffer for offscreen drawing
  int buffer_width   =  antialiasing_level * out_width;
  int buffer_height  =  antialiasing_level * out_height;
  Texture * color_tex = NULL;
  Texture * depth_tex = NULL;
  Framebuffer * fb = NULL;
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

    depth_tex = render_system->createTexture("Depth", buffer_width, buffer_height, 1, Texture::Format::DEPTH16(),
                                             Texture::Dimension::DIM_2D, tex_opts);
    if (!depth_tex)
    {
      THEA_ERROR << "Could not create depth buffer";
      return -1;
    }

    fb = render_system->createFramebuffer("Framebuffer");
    if (!fb)
    {
      THEA_ERROR << "Could not create offscreen framebuffer";
      return -1;
    }

    fb->attach(Framebuffer::AttachmentPoint::COLOR_0, color_tex);
    fb->attach(Framebuffer::AttachmentPoint::DEPTH,   depth_tex);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "Could not render shape")

  matcap_tex = NULL;
  if (!matcap_path.empty())
  {
    try
    {
      Image matcap_img(matcap_path);
      matcap_tex = render_system->createTexture("Matcap", matcap_img);
    }
    THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "Could not create matcap texture")
  }

  tex3d = NULL;
  if (!tex3d_path.empty())
  {
    try
    {
      Image tex3d_img(tex3d_path);
      tex3d = render_system->createTexture("Texture3D", tex3d_img, Texture::Format::AUTO(), Texture::Dimension::DIM_3D);
    }
    THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "%s", "Could not create 3D texture")
  }

  typedef Thea::shared_ptr<Model> ModelPtr;
  TheaArray<ModelPtr> overlay_models;
  for (size_t i = 1; i < model_paths.size(); ++i)
  {
    ModelPtr overlay_model(new Model);
    overlay_model->convert_to_points = (show_points & POINTS_OVERLAY);
    if (!loadModel(*overlay_model, model_paths[i]))
      return -1;

    overlay_models.push_back(overlay_model);
  }

  // Do the rendering
  for (size_t v = 0; v < views.size(); ++v)
  {
    try
    {
      // Initialize the camera
      Camera camera = model.fitCamera(transforms[0], views[v], zoom, buffer_width, buffer_height);
      if (print_camera)
        THEA_CONSOLE << "Camera for view " << v << " is: " << camera.toString();

      // Render the mesh to the offscreen framebuffer
      render_system->pushFramebuffer();
        render_system->setFramebuffer(fb);

        render_system->pushDepthFlags();
        render_system->pushColorFlags();
        render_system->pushShapeFlags();
        render_system->pushTextures();
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
            renderModel(model, primary_color);
          render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();

          // Draw overlay models
          for (size_t i = 0; i < overlay_models.size(); ++i)
          {
            Model const & overlay = *overlay_models[i];
            ColorRGBA overlay_color = getPaletteColor((long)i - 1);
            overlay_color.a() = (!color_by_id && !overlay.is_point_cloud ? 0.5f : 1.0f);

            render_system->setPolygonOffset(true, -1.0f);  // make sure overlays appear on top of primary shape

            render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
              render_system->multMatrix(transforms[i]);
              renderModel(overlay, overlay_color);
            render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
          }

        render_system->setMatrixMode(RenderSystem::MatrixMode::PROJECTION); render_system->popMatrix();
        render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
        render_system->popTextures();
        render_system->popShapeFlags();
        render_system->popColorFlags();
        render_system->popDepthFlags();

        // Grab and save the rendered image
        Image image(Image::Type::RGB_8U, buffer_width, buffer_height);
        color_tex->getImage(image);

        if (antialiasing_level > 1 && !image.rescale(out_width, out_height, 1, Image::Filter::BICUBIC))
        {
          THEA_ERROR << "Could not rescale image to output dimensions";
          return -1;
        }

        string path = out_path;
        if (views.size() > 1)
        {
          path = FilePath::concat(FilePath::parent(path),
                                  FilePath::baseName(path) + format("_%06ld.", (long)v) + FilePath::completeExtension(path));
        }

        image.save(path);

        // Grab and save the depth image
        if (save_depth)
        {
          Image depth_image(Image::Type::LUMINANCE_16U, buffer_width, buffer_height);
          depth_tex->getImage(depth_image);

          if (antialiasing_level > 1 && !depth_image.rescale(out_width, out_height, 1, Image::Filter::BICUBIC))
          {
            THEA_ERROR << "Could not rescale depth image to output dimensions";
            return -1;
          }

          string suffix = (views.size() > 1 ? format("_%06ld.tif", (long)v) : ".tif");
          string depth_path = FilePath::concat(FilePath::parent(out_path), FilePath::baseName(out_path) + "_depth" + suffix);

          depth_image.save(depth_path);
        }

      render_system->popFramebuffer();
    }
    THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not render view %ld of shape", (long)v)
  }

  render_system->destroyTexture(tex3d);
  render_system->destroyTexture(matcap_tex);

  THEA_CONSOLE << "Rendered " << views.size() << " view(s) of the shape";

  return 0;
}

bool
ShapeRendererImpl::usage()
{
  string app_path = FilePath::objectName(Application::programPath());

  THEA_CONSOLE << "";
  THEA_CONSOLE << "Usage: " << app_path << " [OPTIONS] <mesh> <output-image> <width> <height>";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << "  -o <overlay-shape>    (may be mesh or point set)";
  THEA_CONSOLE << "  -s <scope>            (show 'none' | 'primary' | 'overlay' | 'all' shapes";
  THEA_CONSOLE << "                         as points)";
  THEA_CONSOLE << "  -t <transform>        (row-major comma-separated 3x4 or 4x4 matrix,";
  THEA_CONSOLE << "                         applied to all subsequent shapes)";
  THEA_CONSOLE << "  -z <factor>           (zoom factor, default 1)";
  THEA_CONSOLE << "  -v <arg>              (comma-separated 3-vector (viewing direction);";
  THEA_CONSOLE << "                         or 6-vector (direction + eye position);";
  THEA_CONSOLE << "                         or 9-vector (direction + eye + up);";
  THEA_CONSOLE << "                         or string of 3 chars, one for each coordinate,";
  THEA_CONSOLE << "                           each one of +, - or 0;";
  THEA_CONSOLE << "                         or file containing one of the above per line)";
  THEA_CONSOLE << "  -u <up-dir>           (x, y or z, optionally preceded by + or -)";
  THEA_CONSOLE << "  -p <pixels>           (size of points in pixels -- can be fractional)";
  THEA_CONSOLE << "  -c <argb>             (shape color, or 'id' to color faces by face ID and";
  THEA_CONSOLE << "                         points by point ID)";
  THEA_CONSOLE << "  -l <path>             (color faces/points by labels from <path>)";
  THEA_CONSOLE << "  -f <path>             (color vertices/points by features from <path>)";
  THEA_CONSOLE << "  -3 <path>             (color mesh by a 3D texture)";
  THEA_CONSOLE << "  -e                    (accentuate features)";
  THEA_CONSOLE << "  -m <name>             (render submesh with the given name as white, the";
  THEA_CONSOLE << "                         rest of the shape as black)";
  THEA_CONSOLE << "  -k 'ka kd ks ksp'     (Phong material coefficients)";
  THEA_CONSOLE << "  -k <path>             (texture to be used as matcap material)";
  THEA_CONSOLE << "  -b <argb>             (background color)";
  THEA_CONSOLE << "  -a N                  (enable NxN antialiasing: 2 is normal, 4 is very";
  THEA_CONSOLE << "                         high quality)";
  THEA_CONSOLE << "  -0                    (flat shading)";
  THEA_CONSOLE << "  -d                    (also save the depth image)";
  THEA_CONSOLE << "  -w cam                (print the camera parameters)";
  THEA_CONSOLE << "  -i <shift>            (add a shift to how indices are mapped to colors)";
  THEA_CONSOLE << "";

  return false;
}

bool
ShapeRendererImpl::parseTransform(string const & s, Matrix4 & m)
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
ShapeRendererImpl::parseViewDiscrete(string const & s, View & view, bool silent)
{
  if (s.length() != 3)
  {
    if (!silent) THEA_ERROR << "Viewing direction string must have exactly 3 characters, one for each coordinate";
    return false;
  }

  if (s == "000")
  {
    if (!silent) THEA_ERROR << "View direction is zero vector";
    return false;
  }

  view = View();

  for (int i = 0; i < 3; ++i)
  {
    switch (s[i])
    {
      case '+': view.dir[i] =  1; break;
      case '-': view.dir[i] = -1; break;
      case '0': view.dir[i] =  0; break;
      default:
        if (!silent) THEA_ERROR << "Invalid view direction string '" << s << '\'';
        return false;
    }
  }

  if (view.dir.squaredLength() <= 1e-10)
  {
    if (!silent) THEA_ERROR << "View direction is zero vector";
    return false;
  }

  view.dir.unitize();

  if (has_up)
    view.up = view_up;
  else if (s == "0-0")
    view.up = -Vector3::unitZ();
  else if (s == "0+0")
    view.up = Vector3::unitZ();
  else
    view.up = Vector3::unitY();

  return true;
}

bool
ShapeRendererImpl::parseViewContinuous(string const & s, View & view, bool silent)
{
  double dx, dy, dz;
  double ex, ey, ez;
  double ux, uy, uz;
  char trailing[2];  // to make sure there's nothing after the 9th number
  int num_params = sscanf(s.c_str(), " %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf , %lf %1s",
                          &dx, &dy, &dz, &ex, &ey, &ez, &ux, &uy, &uz, trailing);

  if (!(num_params == 9 || num_params == 6 || num_params == 3))
  {
    if (!silent) THEA_ERROR << "Invalid view string '" << s << '\'';
    return false;
  }

  view = View();

  view.dir = Vector3(dx, dy, dz);
  if (view.dir.squaredLength() <= 1e-10)
  {
    if (!silent) THEA_ERROR << "View direction is zero vector";
    return false;
  }

  view.dir.unitize();

  if (num_params == 6)
  {
    view.has_eye = true;
    view.eye = Vector3(ex, ey, ez);
  }

  if (num_params == 9)
  {
    view.up = Vector3(ux, uy, uz);

    if (view.up.squaredLength() <= 1e-10)
    {
      if (!silent) THEA_ERROR << "View up is zero vector";
      return false;
    }

    view.up.unitize();
  }
  else if (has_up)
    view.up = view_up;
  else
  {
    Real d = view.dir.dot(Vector3::unitY());
    if (Math::fuzzyEq(d, (Real)-1))
      view.up = -Vector3::unitZ();
    else if (Math::fuzzyEq(d, (Real)1))
      view.up = Vector3::unitZ();
    else
      view.up = Vector3::unitY();
  }

  return true;
}

bool
ShapeRendererImpl::parseViewFile(string const & path)
{
  ifstream in(path.c_str());
  if (!in)
  {
    THEA_ERROR << "Could not open view directions file '" << path << '\'';
    return false;
  }

  string line;
  View view;
  while (getline(in, line))
  {
    line = trimWhitespace(line);
    if (line.empty())
      continue;
    else if (line.length() == 3)
    {
      if (!parseViewDiscrete(line, view))
        return false;
    }
    else
    {
      if (!parseViewContinuous(line, view))
        return false;
    }

    views.push_back(view);
  }

  return true;
}

bool
ShapeRendererImpl::parseViewUp(string const & s, Vector3 & up)
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
ShapeRendererImpl::parseColor(string const & s, ColorRGBA & c)
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

// Returns number of components found
int
ShapeRendererImpl::parseVector(string const & str, Vector4 & v)
{
  string s = trimWhitespace(str);
  if (s.length() >= 2 && s[0] == '(' && s[s.length() - 1] == ')')
    s = s.substr(1, s.length() - 2);

  double x[4];
  size_t num_parsed = sscanf(s.c_str(), " %lf , %lf , %lf , %lf", &x[0], &x[1], &x[2], &x[3]);
  if (num_parsed < 2)
    num_parsed = sscanf(s.c_str(), " %lf %lf %lf %lf", &x[0], &x[1], &x[2], &x[3]);

  if (num_parsed > 0)
  {
    v.fill(-1);
    for (size_t i = 0; i < num_parsed; ++i)
      v[i] = (Real)x[i];
  }
  else
    THEA_ERROR << "Could not parse vector '" << s << '\'';

  return (int)num_parsed;
}

bool
ShapeRendererImpl::parseArgs(int argc, char ** argv)
{
  if (argc < 5)
    return usage();

  resetArgs();

  Matrix4 current_transform = Matrix4::identity();

  argv++;
  argc--;
  int pos = 0;

  bool explicit_primary_color = false;

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

        case 's':
        {
          if (argc < 1) { THEA_ERROR << "-s: Shapes to show as points not specified"; return false; }
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

        case 'z':
        {
          if (argc < 1) { THEA_ERROR << "-z: Zoom not specified"; return false; }
          if (sscanf(*argv, " %f", &zoom) != 1)
          {
            THEA_ERROR << "Could not parse zoom '" << *argv << '\'';
            return false;
          }
          if (zoom <= 0)
          {
            THEA_ERROR << "Invalid zoom " << zoom;
            return false;
          }
          argv++; argc--; break;
        }

        case 'v':
        {
          if (argc < 1) { THEA_ERROR << "-v: View parameters not specified"; return false; }
          View view;
          bool status = false;
          if (strlen(*argv) == 3)
            status = parseViewDiscrete(*argv, view, true);
          else
            status = parseViewContinuous(*argv, view, true);

          if (status)
            views.push_back(view);
          else
          {
            if (FileSystem::fileExists(*argv))
            {
              if (!parseViewFile(*argv))
                return false;
            }
            else
            {
              THEA_ERROR << "Could not parse view direction '" << *argv << '\'';
              return false;
            }
          }

          argv++; argc--; break;
        }

        case 'u':
        {
          if (argc < 1) { THEA_ERROR << "-u: Up direction not specified"; return false; }
          if (!parseViewUp(*argv, view_up)) return false;
          has_up = true;
          argv++; argc--; break;
        }

        case 'p':
        {
          if (argc < 1) { THEA_ERROR << "-p: Point size not specified"; return false; }
          if (sscanf(*argv, " %f", &point_size) != 1)
          {
            THEA_ERROR << "Could not parse point size '" << *argv << '\'';
            return false;
          }
          if (zoom <= 0)
          {
            THEA_ERROR << "Invalid point size " << point_size;
            return false;
          }
          argv++; argc--; break;
        }

        case 'c':
        {
          if (argc < 1) { THEA_ERROR << "-c: Color not specified"; return false; }

          if (toLower(trimWhitespace(*argv)) == "id")
            color_by_id = true;
          else
          {
            color_by_id = false;
            if (!parseColor(*argv, primary_color))
              return false;
          }

          explicit_primary_color = true;

          argv++; argc--; break;
        }

        case 'l':
        {
          if (argc < 1) { THEA_ERROR << "-l: Labels not specified"; return false; }

          color_by_label = true;
          labels_path = *argv;

          argv++; argc--; break;
        }

        case 'f':
        {
          if (argc < 1) { THEA_ERROR << "-f: Features not specified"; return false; }

          color_by_features = true;
          features_path = *argv;

          argv++; argc--; break;
        }

        case 'e':
        {
          accentuate_features = true;
          break;
        }

        case '3':
        {
          if (argc < 1) { THEA_ERROR << "-3: Texture image not specified"; return false; }

          color_by_tex3d = true;
          tex3d_path = *argv;

          argv++; argc--; break;
        }

        case 'm':
        {
          if (argc < 1) { THEA_ERROR << "-m: Selected mesh not specified"; return false; }

          selected_mesh = toLower(*argv);

          argv++; argc--; break;
        }

        case 'k':
        {
          if (argc < 1) { THEA_ERROR << "-k: Material not specified"; return false; }

          if (FileSystem::fileExists(*argv))
            matcap_path = *argv;
          else
          {
            Vector4 m;
            int num_fields = parseVector(*argv, m);
            if (num_fields <= 0)
              return false;

            for (int i = 0; i < num_fields; ++i)
              if (m[i] >= -0.001)
                material[i] = m[i];
          }

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

        case '0':
        {
          flat = true;
          break;
        }

        case 'd':
        {
          save_depth = true;
          break;
        }

        case 'w':
        {
          if (argc < 1) { THEA_ERROR << "-w: Field to write not specified"; return false; }
          if (string(*argv) == "cam")
            print_camera = true;
          else
          {
            THEA_ERROR << "Unknown field to write: '" << *argv << '\'';
            return false;
          }
          argv++; argc--; break;
        }

        case 'i':
        {
          if (argc < 1) { THEA_ERROR << "-i: Palette shift not specified"; return false; }
          if (sscanf(*argv, " %d", &palette_shift) != 1)
          {
            THEA_ERROR << "Could not parse palette shift: '" << *argv << '\'';
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
          THEA_ERROR << "Extra positional argument: '" << arg << '\'';
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

  // If no primary color was explicitly specified, set it to white by default if we're also going to multiply it by the matcap
  // or 3D texture color. The assumption is that the latter colors should be presented accurately unless the user explicitly
  // modulates them by an additional color.
  if (!explicit_primary_color && (!matcap_path.empty() || color_by_tex3d))
    primary_color = ColorRGBA(1, 1, 1, 1);

  if (views.empty())
  {
    View view;
    if (has_up)
    {
      view.up = view_up;
      int axis = view_up.maxAbsAxis();
      switch (axis)
      {
        case 0:  view.dir = Vector3(-Math::sign(view_up[0]),   Math::sign(view_up[0]),  -1); break;
        case 2:  view.dir = Vector3(-Math::sign(view_up[2]),   1,                       -Math::sign(view_up[2])); break;
        default: view.dir = Vector3(-Math::sign(view_up[1]),  -Math::sign(view_up[1]),  -1);
      }
    }
    else
    {
      view.up = Vector3(0, 1, 0);
      view.dir = Vector3(-1, -1, -1);
    }

    views.push_back(view);
  }

  return true;
}

struct MeshReadCallback : public MeshCodec<Mesh>::ReadCallback
{
  FaceIndexMap & tri_ids;
  FaceIndexMap & quad_ids;

  MeshReadCallback(FaceIndexMap & tri_ids_, FaceIndexMap & quad_ids_) : tri_ids(tri_ids_), quad_ids(quad_ids_) {}

  void faceRead(Mesh * mesh, long index, Mesh::FaceHandle face)
  {
    if (face.hasTriangles())
    {
      long base_tri = face.getFirstTriangle();
      for (long i = 0; i < face.numTriangles(); ++i)
      {
        FaceRef ref(mesh, base_tri + i);
        tri_ids[ref] = (uint32)index;

        // THEA_CONSOLE << ref.first->getName() << ": Triangle " << ref.second << " has id " << index;
      }
    }

    if (face.hasQuads())
    {
      long base_quad = face.getFirstQuad();
      for (long i = 0; i < face.numQuads(); ++i)
      {
        FaceRef ref(mesh, base_quad + i);
        quad_ids[ref] = (uint32)index;

        // THEA_CONSOLE << ref.first->getName() << ": Quad " << ref.second << " has id " << index;
      }
    }
  }
};

bool
enableWireframe(Mesh & mesh)
{
  mesh.setWireframeEnabled(true);
  return false;
}

bool
flattenFaces(Mesh & mesh)
{
  mesh.isolateFaces();
  mesh.computeAveragedVertexNormals();
  return false;
}

bool
averageNormals(Mesh & mesh)
{
  if (!mesh.hasNormals())
    mesh.computeAveragedVertexNormals();

  return false;
}

ColorRGBA8
indexToColor(uint32 index, bool is_point)
{
  ColorRGBA8 color((uint8)((index      ) & 0xFF),
                   (uint8)((index >>  8) & 0xFF),
                   (uint8)((index >> 16) & 0xFF),
                   255);

  if (is_point)
  {
    if (color.b() & 0x80)
      THEA_WARNING << "Too many points -- point IDs will overflow and not be one-to-one!";

    color.b() = (color.b() | 0x80);
  }

  return color;
}

struct FaceColorizer
{
  ShapeRendererImpl const * parent;
  FaceIndexMap const & tri_ids;
  FaceIndexMap const & quad_ids;
  TheaArray<int> const * labels;

  FaceColorizer(ShapeRendererImpl const * parent_, FaceIndexMap const & tri_ids_, FaceIndexMap const & quad_ids_,
                TheaArray<int> const * labels_ = NULL)
  : parent(parent_), tri_ids(tri_ids_), quad_ids(quad_ids_), labels(labels_) {}

  bool operator()(Mesh & mesh)
  {
    mesh.isolateFaces();
    if (labels) mesh.computeAveragedVertexNormals();
    mesh.addColors();

    Mesh::IndexArray const & tris = mesh.getTriangleIndices();
    ColorRGBA8 color;
    for (size_t i = 0; i < tris.size(); i += 3)
    {
      FaceRef face(&mesh, (long)i / 3);
      FaceIndexMap::const_iterator existing = tri_ids.find(face);
      if (existing == tri_ids.end())
        throw Error(format("Could not find index of triangle %ld in mesh '%s'", (long)i / 3, mesh.getName()));

      uint32 id = existing->second;
      if (labels)
      {
        if ((size_t)id >= labels->size())
          color = ColorRGBA8(255, 0, 0, 255);
        else
          color = parent->getLabelColor((*labels)[(size_t)id]);
      }
      else
        color = indexToColor(id, false);

      mesh.setColor((long)tris[i    ], color);
      mesh.setColor((long)tris[i + 1], color);
      mesh.setColor((long)tris[i + 2], color);
    }

    Mesh::IndexArray const & quads = mesh.getQuadIndices();
    for (size_t i = 0; i < quads.size(); i += 4)
    {
      FaceRef face(&mesh, (long)i / 4);
      FaceIndexMap::const_iterator existing = quad_ids.find(face);
      if (existing == quad_ids.end())
        throw Error(format("Could not find index of quad %ld in mesh '%s'", (long)i / 4, mesh.getName()));

      uint32 id = existing->second;
      if (labels)
      {
        if ((size_t)id >= labels->size())
          color = ColorRGBA8(255, 0, 0, 255);
        else
          color = parent->getLabelColor((*labels)[(size_t)id]);
      }
      else
        color = indexToColor(id, false);

      mesh.setColor((long)quads[i    ], color);
      mesh.setColor((long)quads[i + 1], color);
      mesh.setColor((long)quads[i + 2], color);
      mesh.setColor((long)quads[i + 3], color);
    }

    return false;
  }
};

bool
ShapeRendererImpl::loadLabels(Model & model, FaceIndexMap const * tri_ids, FaceIndexMap const * quad_ids)
{
  if (labels_path.empty())
  {
    THEA_ERROR << "No labels file specified";
    return false;
  }

  TheaArray<int> labels;
  string ext = toLower(FilePath::extension(labels_path));
  if (ext == "lab")
  {
    // <label1 name>
    // <face/point ids starting from 1 labeled according to label1>
    // <label2 name>
    // <face/point ids starting from 1 labeled according to label2>
    // ...

    ifstream in(labels_path.c_str());
    if (!in)
    {
      THEA_ERROR << "Could not open labels file '" << labels_path << '\'';
      return false;
    }

    typedef TheaUnorderedMap<string, int> LabelIndexMap;
    LabelIndexMap label_indices;

    string line;
    while (getline(in, line))
    {
      string label_name = trimWhitespace(line);
      if (label_name.empty())
        continue;

      if (!getline(in, line))
      {
        THEA_ERROR << "Could not read list of elements for label '" << label_name << '\'';
        return false;
      }

      int label_index = -1;
      LabelIndexMap::const_iterator existing = label_indices.find(label_name);
      if (existing == label_indices.end())
      {
        label_index = (int)labelHash(label_name);
        label_indices[label_name] = label_index;
      }
      else
        label_index = existing->second;

      istringstream line_in(line);
      long index;
      while (line_in >> index)
      {
        index--;  // first face has index 1 in this format
        if (index < 0)
        {
          THEA_ERROR << "Index " << index << " of element to be labeled is out of bounds";
          return false;
        }

        if (index >= (long)labels.size())
          labels.resize((size_t)(2 * (index + 1)), -1);

        labels[(size_t)index] = label_index;
      }
    }
  }
  else if (ext == "labels")
  {
    // <label1 name>
    // <representative face ids starting from 0, one per submesh labeled according to label1>
    // <label2 name>
    // <representative face ids starting from 0, one per submesh labeled according to label2>
    // ...

    if (model.is_point_cloud)
    {
      THEA_ERROR << "The .labels format is supported only for meshes";
      return false;
    }

    ifstream in(labels_path.c_str());
    if (!in)
    {
      THEA_ERROR << "Could not open labels file '" << labels_path << '\'';
      return false;
    }

    // Build mapping from face indices to their parent meshes
    TheaArray<Mesh const *> face_meshes((size_t)(tri_ids->size() + quad_ids->size()), NULL);
    for (FaceIndexMap::const_iterator fi = tri_ids->begin(); fi != tri_ids->end(); ++fi)
      face_meshes[(size_t)fi->second] = fi->first.first;

    for (FaceIndexMap::const_iterator fi = quad_ids->begin(); fi != quad_ids->end(); ++fi)
      face_meshes[(size_t)fi->second] = fi->first.first;

    // Build map from meshes to their labels
    typedef TheaUnorderedMap<Mesh const *, int> MeshLabelMap;
    MeshLabelMap mesh_labels;

    typedef TheaUnorderedMap<string, int> LabelIndexMap;
    LabelIndexMap label_indices;

    string line;
    while (getline(in, line))
    {
      string label_name = trimWhitespace(line);
      if (label_name.empty())
        continue;

      if (!getline(in, line))
      {
        THEA_ERROR << "Could not read list of representative faces for label '" << label_name << '\'';
        return false;
      }

      int label_index = -1;
      LabelIndexMap::const_iterator existing = label_indices.find(label_name);
      if (existing == label_indices.end())
      {
        label_index = (int)labelHash(label_name);
        label_indices[label_name] = label_index;
      }
      else
        label_index = existing->second;

      istringstream line_in(line);
      long rep_index;
      while (line_in >> rep_index)
      {
        if (rep_index < 0)
        {
          THEA_ERROR << "Index " << rep_index << " of element to be labeled is out of bounds";
          return false;
        }

        Mesh const * mesh = face_meshes[rep_index];
        if (!mesh)
        {
          THEA_ERROR << "Face " << rep_index << " not associated with a parent mesh";
          return false;
        }

        mesh_labels[mesh] = label_index;
      }
    }

    // Map from mesh labels to face labels
    labels.resize(face_meshes.size(), -1);
    for (size_t i = 0; i < labels.size(); ++i)
    {
      if (face_meshes[i])
      {
        MeshLabelMap::const_iterator existing = mesh_labels.find(face_meshes[i]);
        labels[i] = (existing != mesh_labels.end() ? existing->second : -1);
      }
    }
  }
  else
  {
    // <name of label of face/point 0>
    // <name of label of face/point 1>
    // <name of label of face/point 2>
    // ...

    ifstream in(labels_path.c_str());
    if (!in)
    {
      THEA_ERROR << "Could not open labels file '" << labels_path << '\'';
      return false;
    }

    typedef TheaUnorderedMap<string, int> LabelIndexMap;
    LabelIndexMap label_indices;

    string line;
    while (getline(in, line))
    {
      string label_name = trimWhitespace(line);

      int label_index = -1;
      LabelIndexMap::const_iterator existing = label_indices.find(label_name);
      if (existing == label_indices.end())
      {
        label_index = (int)labelHash(label_name);
        label_indices[label_name] = label_index;
      }
      else
        label_index = existing->second;

      labels.push_back(label_index);
    }
  }

  if (model.is_point_cloud)
  {
    if (labels.size() != model.points.size())
    {
      THEA_ERROR << "Label points are not in one-one correspondence with model points";
      return false;
    }

    model.point_colors.resize(model.points.size());
    for (size_t i = 0; i < model.points.size(); ++i)
      model.point_colors[i] = (labels[i] < 0 ? ColorRGBA(1, 0, 0, 1) : getLabelColor(labels[i]));
  }
  else
  {
    alwaysAssertM(tri_ids && quad_ids, "Face IDs necessary for coloring by label");
    FaceColorizer label_colorizer(this, *tri_ids, *quad_ids, &labels);
    model.mesh_group.forEachMeshUntil(&label_colorizer);
  }

  THEA_CONSOLE << "Read labels from '" << labels_path << '\'';

  return true;
}

istream &
getNextNonBlankLine(istream & in, string & line)
{
  while (getline(in, line))
  {
    if (!trimWhitespace(line).empty())
      break;
  }

  return in;
}

ColorRGB
featToColor(Real f0, Real const * f1, Real const * f2)
{
  if (!f2)
  {
    if (!f1)
      return ColorRGB::jetColorMap(0.2 + 0.6 * f0);
    else
      return ColorRGB(f0, *f1, 1.0f);
  }
  else
    return ColorRGB(f0, *f1, *f2);
}

typedef KDTreeN<Vector3, 3> PointKDTree;

struct VertexColorizer
{
  VertexColorizer(PointKDTree * fkdtree_, Real const * feat_vals0_, Real const * feat_vals1_, Real const * feat_vals2_)
  : fkdtree(fkdtree_), feat_vals0(feat_vals0_), feat_vals1(feat_vals1_), feat_vals2(feat_vals2_) {}

  bool operator()(Mesh & mesh)
  {
    static int const MAX_NBRS = 8;  // ok to have a few neighbors for output quality -- this is an offline renderer
    Real scale = std::max(0.2f * fkdtree->getBounds().getExtent().length(), (Real)1.0e-8);
    Real scale2 = scale * scale;

    mesh.addColors();

    Mesh::VertexArray const & vertices = mesh.getVertices();
    BoundedSortedArrayN<MAX_NBRS, PointKDTree::NeighborPair> nbrs;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      nbrs.clear();
      long num_nbrs = fkdtree->kClosestPairs<Algorithms::MetricL2>(vertices[i], nbrs, 2 * scale);
      if (num_nbrs <= 0)
        num_nbrs = fkdtree->kClosestPairs<Algorithms::MetricL2>(vertices[i], nbrs);

      if (num_nbrs > 0)
      {
        ColorRGB c(0, 0, 0);
        double sum_weights = 0;
        for (int j = 0; j < num_nbrs; ++j)
        {
          double dist = nbrs[j].getDistance<MetricL2>();
          double weight = Math::fastMinusExp(dist * dist / scale2);
          long nn_index = nbrs[j].getTargetIndex();
          sum_weights += weight;
          c += weight * featToColor(feat_vals0[nn_index],
                                    (feat_vals1 ? &feat_vals1[nn_index] : NULL),
                                    (feat_vals2 ? &feat_vals2[nn_index] : NULL));
        }

        mesh.setColor((long)i, sum_weights > 0 ? c / sum_weights : c);
      }
      else
      {
        THEA_WARNING << "No nearest neighbor found!";
        mesh.setColor((long)i, ColorRGB(1, 1, 1));
      }
    }

    return false;
  }

  PointKDTree * fkdtree;
  Real const * feat_vals0;
  Real const * feat_vals1;
  Real const * feat_vals2;
};

bool
ShapeRendererImpl::loadFeatures(Model & model)
{
  TheaArray<Vector3> feat_pts;
  TheaArray< TheaArray<Real> > feat_vals(1);
  try
  {
    ifstream in(features_path.c_str());
    if (!in)
      throw Error("Couldn't open file");

    string line;
    Vector3 p;
    Real f;
    while (getNextNonBlankLine(in, line))
    {
      istringstream line_in(line);
      if (!(line_in >> p[0] >> p[1] >> p[2] >> f))
        throw Error("Couldn't read feature");

      feat_pts.push_back(p);
      feat_vals[0].push_back(f);

      if (feat_pts.size() == 1)
      {
        while (line_in >> f)
        {
          feat_vals.push_back(TheaArray<Real>());
          feat_vals.back().push_back(f);
        }
      }
      else
      {
        for (size_t i = 1; i < feat_vals.size(); ++i)
        {
          if (!(line_in >> f))
            throw Error("Couldn't read feature");

          feat_vals[i].push_back(f);
        }
      }
    }

    if (feat_pts.empty())
      return true;

    if (accentuate_features)
    {
      if (feat_vals.size() == 3)
      {
        Real abs_max = -1;
        for (size_t i = 0; i < feat_vals.size(); ++i)
        {
          for (size_t j = 0; j < feat_vals[i].size(); ++j)
          {
            Real abs_feat_val = fabs(feat_vals[i][j]);
            if (abs_feat_val > abs_max)
              abs_max = abs_feat_val;
          }
        }

        if (abs_max > 0)
        {
          for (size_t i = 0; i < feat_vals.size(); ++i)
            for (size_t j = 0; j < feat_vals[i].size(); ++j)
              feat_vals[i][j] = Math::clamp(0.5 * (feat_vals[i][j]  / abs_max + 1), (Real)0, (Real)1);
        }
      }
      else
      {
        for (size_t i = 0; i < feat_vals.size(); ++i)
        {
          TheaArray<Real> sorted = feat_vals[i];
          sort(sorted.begin(), sorted.end());

          size_t tenth = (int)(0.1 * sorted.size());
          Real lo = *(sorted.begin() + tenth);

          size_t ninetieth = (int)(0.9 * sorted.size());
          Real hi = *(sorted.begin() + ninetieth);

          Real range = hi - lo;

          if (range < 1e-20)
          {
            lo = sorted.front();
            hi = sorted.back();
            range = hi - lo;

            if (range < 1e-20)
              continue;
          }

          if (sorted[0] >= 0)  // make a guess if this is a [0, 1] feature (e.g. SDF) or a [-1, 1] feature (e.g. curvature)
          {
            for (size_t j = 0; j < feat_vals[i].size(); ++j)
              feat_vals[i][j] = Math::clamp((feat_vals[i][j] - lo) / range, (Real)0, (Real)1);
          }
          else
          {
            Real abs_max = max(fabs(lo), fabs(hi));
            for (size_t j = 0; j < feat_vals[i].size(); ++j)
              feat_vals[i][j] = Math::clamp((feat_vals[i][j] + abs_max) / (2 * abs_max), (Real)0, (Real)1);
          }
        }
      }
    }
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, WARNING, "Couldn't load model features from '%s'", features_path.c_str())

  if (model.is_point_cloud)
  {
    if (feat_vals[0].size() != model.points.size())
    {
      THEA_ERROR << "Feature points are not in one-one correspondence with model points";
      return false;
    }

    model.point_colors.resize(model.points.size());
    for (size_t i = 0; i < model.points.size(); ++i)
    {
      model.point_colors[i] = featToColor(feat_vals[0][i],
                                          (feat_vals.size() > 1 ? &feat_vals[1][i] : NULL),
                                          (feat_vals.size() > 2 ? &feat_vals[2][i] : NULL));
    }
  }
  else
  {
    PointKDTree fkdtree(feat_pts.begin(), feat_pts.end());
    VertexColorizer visitor(&fkdtree,
                            &feat_vals[0][0],
                            feat_vals.size() > 1 ? &feat_vals[1][0] : NULL,
                            feat_vals.size() > 2 ? &feat_vals[2][0] : NULL);
    model.mesh_group.forEachMeshUntil(&visitor);
  }

  THEA_CONSOLE << "Read features from '" << features_path << '\'';

  return true;
}

void
ShapeRendererImpl::colorizeMeshSelection(MG & mg, uint32 parent_id)
{
  for (MG::MeshConstIterator mi = mg.meshesBegin(); mi != mg.meshesEnd(); ++mi)
  {
    Mesh & mesh = **mi;

    uint32 mesh_id = 0;
    if (parent_id != 0)
      mesh_id = parent_id;
    else if (!selected_mesh.empty() && toLower(mesh.getName()) == selected_mesh)
      mesh_id = 0xFFFFFF;

    ColorRGBA8 color = indexToColor(mesh_id, false);

    mesh.isolateFaces();
    mesh.addColors();

    Mesh::IndexArray const & tris = mesh.getTriangleIndices();
    for (size_t i = 0; i < tris.size(); i += 3)
    {
      mesh.setColor((long)tris[i    ], color);
      mesh.setColor((long)tris[i + 1], color);
      mesh.setColor((long)tris[i + 2], color);
    }

    Mesh::IndexArray const & quads = mesh.getQuadIndices();
    for (size_t i = 0; i < quads.size(); i += 4)
    {
      mesh.setColor((long)quads[i    ], color);
      mesh.setColor((long)quads[i + 1], color);
      mesh.setColor((long)quads[i + 2], color);
      mesh.setColor((long)quads[i + 3], color);
    }
  }

  for (MG::GroupConstIterator ci = mg.childrenBegin(); ci != mg.childrenEnd(); ++ci)
  {
    MG & child = **ci;
    uint32 child_id = 0;

    if (parent_id != 0)
      child_id = parent_id;
    else if (!selected_mesh.empty() && toLower(child.getName()) == selected_mesh)
      child_id = 0xFFFFFF;

    colorizeMeshSelection(child, child_id);
  }
}

bool
ShapeRendererImpl::loadModel(Model & model, string const & path)
{
  model.mesh_group.clear();
  model.points.clear();
  model.is_point_cloud = false;

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
        THEA_ERROR << "Could not read point " << model.points.size() << " from '" << path << '\'';
        return false;
      }

      model.points.push_back(Vector3((Real)x, (Real)y, (Real)z));
    }

    model.is_point_cloud = true;

    if (color_by_label)
    {
      if (!loadLabels(model, NULL, NULL))
        return false;
    }
    else if (color_by_features)
    {
      if (!loadFeatures(model))
        return false;
    }

    THEA_CONSOLE << "Read " << model.points.size() << " points from '" << path << '\'';
  }
  else
  {
    FaceIndexMap tri_ids, quad_ids;
    try
    {
      MeshReadCallback callback(tri_ids, quad_ids);
      model.mesh_group.load(path, Codec_AUTO(), &callback);
    }
    THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "Could not load model from '%s'", path.c_str())

    if (model.convert_to_points)
    {
      MeshSampler<Mesh> sampler(model.mesh_group);

      double out_scale = min(out_width, out_height);
      long max_samples = (long)ceil(10 * out_scale);
      long npoints = sampler.sampleEvenlyByArea(max_samples, model.points);

      THEA_CONSOLE << "Sampled " << npoints << " points from '" << path << '\'';

      model.mesh_group.clear();
      model.is_point_cloud = true;
    }
    else
    {
      bool needs_normals = true;

      if (color_by_id)
      {
        FaceColorizer id_colorizer(this, tri_ids, quad_ids);
        model.mesh_group.forEachMeshUntil(&id_colorizer);
        needs_normals = false;
      }
      else if (color_by_label)
      {
        if (!loadLabels(model, &tri_ids, &quad_ids))
          return false;

        needs_normals = false;
      }
      else if (color_by_features)
      {
        if (!loadFeatures(model))
          return false;
      }
      else if (!selected_mesh.empty())
      {
        colorizeMeshSelection(model.mesh_group, 0);
        needs_normals = false;
      }

      if (needs_normals)
      {
        if (flat)
          model.mesh_group.forEachMeshUntil(flattenFaces);
        else
          model.mesh_group.forEachMeshUntil(averageNormals);
      }

#ifdef DRAW_EDGES
      model.mesh_group.forEachMeshUntil(enableWireframe);
#endif
    }
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
      for (size_t i = 0; i < vertices.size(); ++i)
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
    for (size_t i = 0; i < model.points.size(); ++i)
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
    for (size_t i = 0; i < tri_array.size(); ++i)
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
    for (size_t i = 0; i < model.points.size(); ++i)
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
Model::fitCamera(Matrix4 const & transform, View const & view, Real zoom, int width, int height)
{
  // Orientation
  Ball3 bsphere = modelBSphere(*this, transform);
  Vector3 center = bsphere.getCenter();
  Real diameter = bsphere.getDiameter();

  // Make absolutely sure these are unit vectors
  Vector3 dir  =  view.dir.unit();
  Vector3 up   =  view.up.unit();
  Vector3 eye;

  if (view.has_eye)
  {
    eye = view.eye;
  }
  else
  {
    Real camera_separation = diameter > 1.0e-10f ? 2.1 * diameter : 1.0e-10f;
    eye = center - camera_separation * dir;
  }

  CoordinateFrame3 cframe = CoordinateFrame3::fromViewFrame(eye, eye + dir, up);

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

  Real center_distance = std::max((bsphere.getCenter() - eye).dot(dir), (Real)0);
  Real near_dist = std::max(center_distance - 1.1f * diameter, 0.01f * diameter);
  Real far_dist  = center_distance + 2 * diameter;

  hw = (hw / zoom) * (0.5f * near_dist);
  hh = (hh / zoom) * (0.5f * near_dist);

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
initMeshShader(Shader & shader, Vector4 const & material, Texture * matcap_tex = NULL, Texture * tex3d = NULL,
               AxisAlignedBox3 const & bbox = AxisAlignedBox3())
{
  static string const VERTEX_SHADER =
"varying vec3 src_pos;  // position in mesh coordinates\n"
"varying vec3 normal;  // normal in camera space\n"
"\n"
"void main()\n"
"{\n"
"  gl_Position = ftransform();\n"
"\n"
"  src_pos = vec3(gl_Vertex);\n"
"  normal = gl_NormalMatrix * gl_Normal;\n"
"\n"
"  gl_FrontColor = gl_Color;\n"
"  gl_BackColor = gl_Color;\n"
"}\n";

  static string const FRAGMENT_SHADER_HEADER_1 =
"uniform float two_sided;\n"
"varying vec3 src_pos;  // position in mesh coordinates\n"
"varying vec3 normal;  // normal in camera space\n";

  static string const FRAGMENT_SHADER_HEADER_PHONG =
"uniform vec3 ambient_color;\n"
"uniform vec3 light_dir;  // must be specified in camera space, pointing towards object\n"
"uniform vec3 light_color;\n"
"uniform vec4 material;  // [ka, kl, <ignored>, <ignored>]\n";

  static string const FRAGMENT_SHADER_HEADER_MATCAP =
"uniform sampler2D matcap_tex;\n";

  static string const FRAGMENT_SHADER_HEADER_TEX3D =
"uniform sampler3D tex3d;\n"
"uniform vec3 bbox_lo;\n"
"uniform vec3 bbox_hi;\n";

  static string const FRAGMENT_SHADER_BODY_1 =
"void main()\n"
"{\n"
"  vec3 N = normalize(normal);\n"
"  vec4 color = gl_Color;\n";

  static string const FRAGMENT_SHADER_BODY_TEX3D =
"  vec3 tex3d_p = (src_pos - bbox_lo) / (bbox_hi - bbox_lo);\n"
"  vec4 tex3d_color = texture3D(tex3d, tex3d_p);\n"
"  color.rgb = mix(color.rgb, tex3d_color.rgb, tex3d_color.a);\n";

  static string const FRAGMENT_SHADER_BODY_PHONG =
"  vec3 ambt_color = material[0] * color.rgb * ambient_color;\n"
"  vec3 L = normalize(light_dir);\n"
"  float NdL = -dot(N, L);\n"
"  vec3 lamb_color = (NdL >= -two_sided) ? material[1] * abs(NdL) * color.rgb * light_color : vec3(0.0, 0.0, 0.0);\n"
"  gl_FragColor = vec4(ambt_color + lamb_color, color.a);\n"
"}\n";

  static string const FRAGMENT_SHADER_BODY_MATCAP =
"  vec2 matcap_uv = N.xy * (two_sided > 0.5 ? sign(N.z) : 1.0);\n"
"  vec4 matcap_color = texture2D(matcap_tex, 0.5 * (matcap_uv + vec2(1.0, 1.0)));\n"
"  gl_FragColor = vec4(color.rgb * matcap_color.rgb, color.a);\n"
"}\n";

  string fragment_shader = FRAGMENT_SHADER_HEADER_1;
  if (matcap_tex)
    fragment_shader += FRAGMENT_SHADER_HEADER_MATCAP;
  else
    fragment_shader += FRAGMENT_SHADER_HEADER_PHONG;

  if (tex3d) fragment_shader += FRAGMENT_SHADER_HEADER_TEX3D;

  fragment_shader += FRAGMENT_SHADER_BODY_1;
  if (tex3d) fragment_shader += FRAGMENT_SHADER_BODY_TEX3D;
  if (matcap_tex)
    fragment_shader += FRAGMENT_SHADER_BODY_MATCAP;
  else
    fragment_shader += FRAGMENT_SHADER_BODY_PHONG;

  try
  {
    shader.attachModuleFromString(Shader::ModuleType::VERTEX, VERTEX_SHADER.c_str());
    shader.attachModuleFromString(Shader::ModuleType::FRAGMENT, fragment_shader.c_str());
  }
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "Could not attach mesh shader module")

  shader.setUniform("two_sided", 1.0f);

  if (matcap_tex)
    shader.setUniform("matcap_tex", matcap_tex);
  else
  {
    shader.setUniform("light_dir", Vector3(-1, -1, -2));
    shader.setUniform("light_color", ColorRGB(1, 1, 1));
    shader.setUniform("ambient_color", ColorRGB(1, 1, 1));
    shader.setUniform("material", material);
  }

  if (tex3d)
  {
    shader.setUniform("tex3d", tex3d);

    // Tightly fit the shape within the texture volume, without distorting the shape or the cubical voxels
    int width   =  tex3d->getWidth();
    int height  =  tex3d->getHeight();
    int depth   =  tex3d->getDepth();
    Vector3 box_ext = bbox.getExtent();
    Vector3 vol_ext(width, depth, height);
    Vector3 dim_ratio = vol_ext / box_ext;
    int fit_axis = dim_ratio.minAxis();
    box_ext = vol_ext / dim_ratio[fit_axis];
    Vector3 center = bbox.getCenter();

    shader.setUniform("bbox_lo", center - 0.5 * box_ext);
    shader.setUniform("bbox_hi", center + 0.5 * box_ext);
  }

  return true;
}

bool
initFaceIDShader(Shader & shader)
{
  static string const VERTEX_SHADER =
"\n"
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
  THEA_STANDARD_CATCH_BLOCKS(return false;, ERROR, "%s", "Could not attach face ID shader module")

  return true;
}

bool
ShapeRendererImpl::renderModel(Model const & model, ColorRGBA const & color)
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

  if (model.is_point_cloud)
  {
    // Initialize the shader
    if (!point_shader)
    {
      point_shader = render_system->createShader("Point shader");
      if (!point_shader)
      {
        THEA_ERROR << "Could not create point shader";
        return false;
      }

      if (!initPointShader(*point_shader))
      {
        THEA_ERROR << "Could not initialize point shader";
        return false;
      }
    }

    render_system->setShader(point_shader);

    render_system->setPointSize(point_size * antialiasing_level);
    render_system->beginPrimitive(RenderSystem::Primitive::POINTS);

      for (size_t i = 0; i < model.points.size(); ++i)
      {
        if (color_by_id)
          render_system->setColor(indexToColor((uint32)i, true));
        else if (color_by_label || color_by_features)
          render_system->setColor(model.point_colors[i]);

        render_system->sendVertex(model.points[i]);
      }

    render_system->endPrimitive();
  }
  else
  {
    // Initialize the shader
    if (color_by_id || !selected_mesh.empty())
    {
      if (!face_id_shader)
      {
        face_id_shader = render_system->createShader("Face ID shader");
        if (!face_id_shader)
        {
          THEA_ERROR << "Could not create face ID shader";
          return false;
        }

        if (!initFaceIDShader(*face_id_shader))
        {
          THEA_ERROR << "Could not initialize face ID shader";
          return false;
        }
      }

      render_system->setShader(face_id_shader);
    }
    else
    {
      if (!mesh_shader)
      {
        mesh_shader = render_system->createShader("Mesh shader");
        if (!mesh_shader)
        {
          THEA_ERROR << "Could not create mesh shader";
          return false;
        }

        if (!initMeshShader(*mesh_shader, material, matcap_tex, tex3d, model.mesh_group.getBounds()))
        {
          THEA_ERROR << "Could not initialize mesh shader";
          return false;
        }
      }

      render_system->setShader(mesh_shader);
    }

    RenderOptions opts = RenderOptions::defaults();
    if (color_by_id || color_by_label || color_by_features)
      opts.useVertexData() = true;

#ifdef DRAW_EDGES
    opts.drawEdges() = true;
    opts.edgeColor() = ColorRGB(0.5, 0.5, 1);
#endif

    if (has_transparency && !color_by_id)
    {
      // First back faces...
      render_system->setCullFace(RenderSystem::CullFace::FRONT);
      model.mesh_group.draw(*render_system, opts);

      // ... then front faces
      render_system->setCullFace(RenderSystem::CullFace::BACK);
      model.mesh_group.draw(*render_system, opts);
    }
    else
    {
      model.mesh_group.draw(*render_system, opts);
    }
  }

  render_system->popColorFlags();
  render_system->popShapeFlags();
  render_system->popShader();

  return true;
}

bool
ShapeRendererImpl::loadPlugins(int argc, char ** argv)
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

  render_system_factory = Application::getRenderSystemManager().getFactory("OpenGL");
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
