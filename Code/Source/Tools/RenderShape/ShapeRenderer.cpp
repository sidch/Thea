#include "ShapeRenderer.hpp"
#include "../../Algorithms/BestFitSphere3.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Algorithms/MeshTriangles.hpp"
#include "../../Graphics/Camera.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshCodec.hpp"
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
#include "../../Random.hpp"
#include "../../UnorderedMap.hpp"
#include "../../Vector3.hpp"
#include <boost/functional/hash.hpp>
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
    string labels_path;
    TheaArray<int> labels;
    string selected_mesh;
    ColorRGBA primary_color;
    ColorRGBA background_color;
    int antialiasing_level;
    PointUsage show_points;
    bool flat;
    bool save_depth;

    bool loadPlugins(int argc, char ** argv);
    bool parseArgs(int argc, char ** argv);
    bool usage();
    bool parseTransform(string const & s, Matrix4 & m);
    bool parseViewDiscrete(string const & s, View & view, bool silent = false);
    bool parseViewContinuous(string const & s, View & view, bool silent = false);
    bool parseViewFile(string const & path);
    bool parseViewUp(string const & s, Vector3 & up);
    bool parseColor(string const & s, ColorRGBA & c);
    void resetArgs();
    bool loadModel(Model & model, string const & path);
    bool loadLabels();
    bool renderModel(Model const & model, ColorRGBA const & color);
    void colorizeMeshSelection(MG & mg, uint32 parent_id);

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
  labels_path = "";
  labels.clear();
  selected_mesh = "";
  primary_color = ColorRGBA(1.0f, 0.9f, 0.8f, 1.0f);
  background_color = ColorRGBA(1, 1, 1, 1);
  antialiasing_level = 1;
  show_points = POINTS_NONE;
  flat = false;
  save_depth = false;
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

  for (array_size_t i = 0; i < args.size(); ++i)
  {
    argv[i + 1] = new char[args[i].length() + 1];
    strcpy(argv[i + 1], args[i].c_str());
  }

  int status = exec((int)argv.size(), &argv[0]);

  for (array_size_t i = 0; i < argv.size(); ++i)
    delete [] argv[i];

  return status;
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

  // Do the rendering
  for (array_size_t v = 0; v < views.size(); ++v)
  {
    try
    {
      // Initialize the camera
      Camera camera = model.fitCamera(transforms[0], views[v], zoom, buffer_width, buffer_height);

      // Render the mesh to the offscreen framebuffer
      render_system->pushFramebuffer();
        render_system->setFramebuffer(fb);

        render_system->pushDepthFlags();
        render_system->pushColorFlags();
        render_system->pushShapeFlags();
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
          for (array_size_t i = 1; i < model_paths.size(); ++i)
          {
            model.convert_to_points = (show_points & POINTS_OVERLAY);
            if (!loadModel(model, model_paths[i]))
              return -1;

            ColorRGBA overlay_color = getPaletteColor((long)i - 1);
            overlay_color.a() = (!color_by_id && !model.is_point_cloud ? 0.5f : 1.0f);

            render_system->setPolygonOffset(-1.0f);  // make sure overlays appear on top of primary shape

            render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
              render_system->multMatrix(transforms[i]);
              renderModel(model, overlay_color);
            render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
          }

        render_system->setMatrixMode(RenderSystem::MatrixMode::PROJECTION); render_system->popMatrix();
        render_system->setMatrixMode(RenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
        render_system->popShapeFlags();
        render_system->popColorFlags();
        render_system->popDepthFlags();

        // Grab and save the rendered image
        Image image(Image::Type::RGB_8U, buffer_width, buffer_height);
        color_tex->getImage(image);

        if (antialiasing_level > 1 && !image.rescale(out_width, out_height, Image::Filter::BICUBIC))
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

          if (antialiasing_level > 1 && !depth_image.rescale(out_width, out_height, Image::Filter::BICUBIC))
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
  THEA_CONSOLE << "  -p <scope>            (show 'none' | 'primary' | 'overlay' | 'all' shapes";
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
  THEA_CONSOLE << "  -s <pixels>           (size of points in pixels -- can be fractional)";
  THEA_CONSOLE << "  -c <argb>             (shape color, or 'id' to color faces by face ID and";
  THEA_CONSOLE << "                         points by point ID)";
  THEA_CONSOLE << "  -l <path>             (color faces/points by labels from <path>)";
  THEA_CONSOLE << "  -m <name>             (render submesh with the given name as white, the";
  THEA_CONSOLE << "                         rest of the shape as black)";
  THEA_CONSOLE << "  -b <argb>             (background color)";
  THEA_CONSOLE << "  -a N                  (enable NxN antialiasing: 2 is normal, 4 is very";
  THEA_CONSOLE << "                         high quality)";
  THEA_CONSOLE << "  -f                    (flat shading)";
  THEA_CONSOLE << "  -d                    (also save the depth image)";
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

        case 's':
        {
          if (argc < 1) { THEA_ERROR << "-s: Point size not specified"; return false; }
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

          argv++; argc--; break;
        }

        case 'l':
        {
          if (argc < 1) { THEA_ERROR << "-l: Labels not specified"; return false; }

          color_by_label = true;
          labels_path = trimWhitespace(*argv);

          argv++; argc--; break;
        }

        case 'm':
        {
          if (argc < 1) { THEA_ERROR << "-m: Selected mesh not specified"; return false; }

          selected_mesh = toLower(*argv);

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

        case 'f':
        {
          flat = true;
          break;
        }

        case 'd':
        {
          save_depth = true;
          break;
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

  if (views.empty())
  {
    View view;
    view.dir = Vector3(-1, -1, -1);
    view.up = Vector3(0, 1, 0);

    views.push_back(view);
  }

  return true;
}

typedef std::pair<Mesh const *, long> FaceRef;
typedef TheaUnorderedMap<FaceRef, uint32> FaceIndexMap;

struct MeshReadCallback : public MeshCodec<Mesh>::ReadCallback
{
  FaceIndexMap & tri_ids;
  FaceIndexMap & quad_ids;

  MeshReadCallback(FaceIndexMap & tri_ids_, FaceIndexMap & quad_ids_) : tri_ids(tri_ids_), quad_ids(quad_ids_) {}

  void faceAdded(Mesh * mesh, long index, IncrementalMeshBuilder<Mesh>::FaceHandle face)
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
  FaceIndexMap const & tri_ids;
  FaceIndexMap const & quad_ids;
  TheaArray<int> const * labels;

  FaceColorizer(FaceIndexMap const & tri_ids_, FaceIndexMap const & quad_ids_, TheaArray<int> const * labels_ = NULL)
  : tri_ids(tri_ids_), quad_ids(quad_ids_), labels(labels_) {}

  bool operator()(Mesh & mesh)
  {
    mesh.isolateFaces();
    if (labels) mesh.computeAveragedVertexNormals();
    mesh.addColors();

    Mesh::IndexArray const & tris = mesh.getTriangleIndices();
    ColorRGBA8 color;
    for (array_size_t i = 0; i < tris.size(); i += 3)
    {
      FaceRef face(&mesh, (long)i / 3);
      FaceIndexMap::const_iterator existing = tri_ids.find(face);
      if (existing == tri_ids.end())
        throw Error(format("Could not find index of triangle %ld in mesh '%s'", (long)i / 3, mesh.getName()));

      uint32 id = existing->second;
      if (labels)
      {
        if ((array_size_t)id >= labels->size())
          color = ColorRGBA8(255, 0, 0, 255);
        else
          color = getPaletteColor((*labels)[(array_size_t)id]);
      }
      else
        color = indexToColor(id, false);

      mesh.setColor((long)tris[i    ], color);
      mesh.setColor((long)tris[i + 1], color);
      mesh.setColor((long)tris[i + 2], color);
    }

    Mesh::IndexArray const & quads = mesh.getQuadIndices();
    for (array_size_t i = 0; i < quads.size(); i += 4)
    {
      FaceRef face(&mesh, (long)i / 4);
      FaceIndexMap::const_iterator existing = quad_ids.find(face);
      if (existing == quad_ids.end())
        throw Error(format("Could not find index of quad %ld in mesh '%s'", (long)i / 4, mesh.getName()));

      uint32 id = existing->second;
      if (labels)
      {
        if ((array_size_t)id >= labels->size())
          color = ColorRGBA8(255, 0, 0, 255);
        else
          color = getPaletteColor((*labels)[(array_size_t)id]);
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
    for (array_size_t i = 0; i < tris.size(); i += 3)
    {
      mesh.setColor((long)tris[i    ], color);
      mesh.setColor((long)tris[i + 1], color);
      mesh.setColor((long)tris[i + 2], color);
    }

    Mesh::IndexArray const & quads = mesh.getQuadIndices();
    for (array_size_t i = 0; i < quads.size(); i += 4)
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
      if (!loadLabels())
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
      if (color_by_id)
      {
        FaceColorizer id_colorizer(tri_ids, quad_ids);
        model.mesh_group.forEachMeshUntil(&id_colorizer);
      }
      else if (color_by_label)
      {
        if (!loadLabels())
          return false;

        FaceColorizer label_colorizer(tri_ids, quad_ids, &labels);
        model.mesh_group.forEachMeshUntil(&label_colorizer);
      }
      else if (!selected_mesh.empty())
      {
        colorizeMeshSelection(model.mesh_group, 0);
      }
      else
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

uint32
labelHash(string const & label)
{
  boost::hash<string> hasher;
  return (uint32)hasher(label);
}

bool
ShapeRendererImpl::loadLabels()
{
  if (labels_path.empty())
  {
    THEA_ERROR << "No labels file specified";
    return false;
  }

  labels.clear();

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
          labels.resize((array_size_t)(2 * (index + 1)), -1);

        labels[(array_size_t)index] = label_index;
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

  THEA_CONSOLE << "Read labels from '" << labels_path << '\'';

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
initMeshShader(Shader & shader)
{
  static string const VERTEX_SHADER =
// "varying vec3 position;  // position in camera space\n"
"varying vec3 normal;  // normal in camera space\n"
"\n"
"void main()\n"
"{\n"
"  gl_Position = ftransform();\n"
"\n"
// "  position = vec3(gl_ModelViewMatrix * gl_Vertex);  // assume rigid transform, so we can drop w\n"
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
// "varying vec3 position;  // position in camera space\n"
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

      for (array_size_t i = 0; i < model.points.size(); ++i)
      {
        if (color_by_id)
          render_system->setColor(indexToColor((uint32)i, true));
        else if (color_by_label)
        {
          int label = labels[(array_size_t)i];
          ColorRGBA color = (label < 0 ? ColorRGBA(1, 0, 0, 1) : getPaletteColor(label));
          render_system->setColor(color);
        }

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

        if (!initMeshShader(*mesh_shader))
        {
          THEA_ERROR << "Could not initialize mesh shader";
          return false;
        }
      }

      render_system->setShader(mesh_shader);
    }

    RenderOptions opts = RenderOptions::defaults();
    if (color_by_id || color_by_label)
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
