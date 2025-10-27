#include "ShapeRenderer.hpp"
#include "../../Algorithms/BvhN.hpp"
#include "../../Algorithms/MeshSampler.hpp"
#include "../../Algorithms/MeshTriangles.hpp"
#include "../../Algorithms/MetricL2.hpp"
#include "../../Graphics/Camera.hpp"
#include "../../Graphics/DisplayMesh.hpp"
#include "../../Graphics/MeshCodec.hpp"
#include "../../Graphics/MeshGroup.hpp"
#include "../../Plugins/GL/GlHeaders.hpp"
#include "../../Application.hpp"
#include "../../Ball3.hpp"
#include "../../BoundedSortedArrayN.hpp"
#include "../../ColorRgba.hpp"
#include "../../CoordinateFrame3.hpp"
#include "../../FilePath.hpp"
#include "../../FileSystem.hpp"
#include "../../Image.hpp"
#include "../../Math.hpp"
#include "../../MatrixWrapper.hpp"
#include "../../MatVec.hpp"
#include "../../Memory.hpp"
#include "../../IPlugin.hpp"
#include "../../Random.hpp"
#include "../../UnorderedMap.hpp"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <utility>

using namespace std;
using namespace Thea;
using namespace Algorithms;
using namespace Graphics;

typedef DisplayMesh Mesh;
typedef MeshGroup<Mesh> MG;

struct View
{
  Vector3 dir;
  bool has_eye;
  Vector3 eye;
  Vector3 up;

  bool has_view_matrix;
  Matrix4 view_matrix;

  bool has_proj_matrix;
  Matrix4 proj_matrix;

  View() : dir(-1, -1, -1), has_eye(false), up(0, 1, 0), has_view_matrix(false), has_proj_matrix(false) {}
};

struct Model
{
  Model(bool convert_to_points_ = false) : convert_to_points(convert_to_points_), is_point_cloud(false), max_id(-1) {}
  Camera fitCamera(Matrix4 const & transform, View const & view, Real zoom, int width, int height);

  bool convert_to_points;
  MG mesh_group;
  bool is_point_cloud;
  intx max_id;
  Array<Vector3> points;
  Array<ColorRgba> point_colors;
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
    static std::atomic_bool has_render_system;
    static IRenderSystemFactory * render_system_factory;
    static IRenderSystem * render_system;
    static IShader * point_shader;
    static IShader * mesh_shader;
    static IShader * face_id_shader;
    static ITexture * matcap_tex;
    static ITexture * tex2d;
    static ITexture * tex3d;

    Array<string> model_paths;
    Array<Matrix4> transforms;
    float zoom;
    string out_path;
    int out_width, out_height;
    Array<View> views;
    bool has_up;
    Vector3 view_up;
    float point_size;
    bool color_by_id;
    bool color_by_leaf;
    bool color_by_leafname;
    bool color_by_label;
    bool color_by_features;
    bool show_edges;
    ColorRgba edge_color;
    string labels_path;
    Array<int> labels;
    string features_path;
    bool accentuate_features;
    bool color_by_tex2d;
    string tex2d_path;
    bool color_by_tex3d;
    string tex3d_path;
    string selected_mesh;
    ColorRgba selected_color;
    bool selected_binary_mask;
    Vector4 material;
    string matcap_path;
    ColorRgba primary_color;
    ColorRgba background_color;
    int antialiasing_level;
    PointUsage show_points;
    bool flat;
    bool two_sided;
    Vector3 light_dir;
    bool toon_shading;
    bool out_normals;
    bool save_color;
    bool save_depth;
    bool save_hitcounts;
    bool print_camera;
    int palette_shift;

    bool loadPlugins(int argc, char ** argv);
    bool parseArgs(int argc, char ** argv);
    bool usage(int argc, char ** argv);
    bool parseTransform(string const & s, Matrix4 & m);
    bool parseViewDiscrete(string const & s, View & view, bool silent = false);
    bool parseViewContinuous(string const & s, View & view, bool silent = false);
    bool parseViewFile(string const & path);
    bool parseViewUp(string const & s, Vector3 & up);
    bool parseARGB(string const & s, uint32 & argb);
    bool parseColor(string const & s, ColorRgba & c);
    int parseVector(string const & str, Vector4 & v);
    void resetArgs();
    bool loadModel(Model & model, string const & path);
    bool loadLabels(Model & model);
    bool loadFeatures(Model & model);
    bool renderModel(Model const & model, ColorRgba const & color);
    void colorizeMeshSelection(MG & mg, uint32 parent_id);
    ColorRgba getPaletteColor(intx n) const;
    ColorRgba getLabelColor(intx label) const;

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

std::atomic_bool ShapeRendererImpl::has_render_system(false);
IRenderSystemFactory * ShapeRendererImpl::render_system_factory = nullptr;
IRenderSystem * ShapeRendererImpl::render_system = nullptr;
IShader * ShapeRendererImpl::point_shader = nullptr;
IShader * ShapeRendererImpl::mesh_shader = nullptr;
IShader * ShapeRendererImpl::face_id_shader = nullptr;
ITexture * ShapeRendererImpl::matcap_tex = nullptr;
ITexture * ShapeRendererImpl::tex2d = nullptr;
ITexture * ShapeRendererImpl::tex3d = nullptr;

ShapeRendererImpl::ShapeRendererImpl(int argc, char * argv[])
{
  resetArgs();
}

ShapeRendererImpl::~ShapeRendererImpl()
{
  if (has_render_system.exchange(false))
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
  color_by_leaf = false;
  color_by_leafname = false;
  color_by_label = false;
  color_by_features = false;
  show_edges = false;
  edge_color = ColorRgba(0.15f, 0.25f, 0.5f, 1.0f);
  labels_path = "";
  labels.clear();
  features_path = "";
  accentuate_features = false;
  color_by_tex2d = false;
  tex2d_path = "";
  color_by_tex3d = false;
  tex3d_path = "";
  selected_mesh = "";
  selected_color = ColorRgba(1, 1, 1, 1);
  selected_binary_mask = true;
  material = Vector4(0.3f, 0.7f, 0.2f, 25);
  matcap_path = "";
  primary_color = ColorRgba(1.0f, 0.9f, 0.8f, 1.0f);
  background_color = ColorRgba(1, 1, 1, 1);
  antialiasing_level = 1;
  show_points = POINTS_NONE;
  flat = false;
  two_sided = true;
  light_dir = Vector3(-1, -1, -2).stableNormalized();
  toon_shading = false;
  out_normals = false;
  save_color = true;
  save_depth = false;
  save_hitcounts = false;
  print_camera = false;
  palette_shift = 0;
}

int
ShapeRendererImpl::exec(string const & cmdline)  // cmdline should not include program name
{
  Array<string> args;
  stringSplit(cmdline, " \t\n\f\r", args, true);

  Array<char *> argv(args.size() + 1);

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
  std::hash<string> hasher;
  string rev = label; reverse(rev.begin(), rev.end());
  return (uint32)((hasher(label) + hasher(rev)) & 0x7FFFFFFF);
}

ColorRgba
ShapeRendererImpl::getPaletteColor(intx n) const
{
  static ColorRgba PALETTE[] = {
    ColorRgba::fromARGB(0xFFFF66FF),
    ColorRgba::fromARGB(0xFFFFFF66),
    ColorRgba::fromARGB(0xFF66FFFF),
    ColorRgba::fromARGB(0xFF3399FF),
    ColorRgba::fromARGB(0xFFFF9933),
    ColorRgba::fromARGB(0xFFFF3399),
    ColorRgba::fromARGB(0xFF99FF33),
    ColorRgba::fromARGB(0xFF9933FF),
    ColorRgba::fromARGB(0xFFFFAAAA),
    ColorRgba::fromARGB(0xFFAAFFAA),
    ColorRgba::fromARGB(0xFFAAAAFF),
    ColorRgba::fromARGB(0xFFFF0000),
    ColorRgba::fromARGB(0xFF00FF00),
    ColorRgba::fromARGB(0xFF0000FF),
    ColorRgba::fromARGB(0xFF00FFFF),
    ColorRgba::fromARGB(0xFFFF00FF),
    ColorRgba::fromARGB(0xFFFFFF00),
    ColorRgba::fromARGB(0xFF800000),
    ColorRgba::fromARGB(0xFF008000),
    ColorRgba::fromARGB(0xFF000080),
    ColorRgba::fromARGB(0xFF008080),
    ColorRgba::fromARGB(0xFF800080),
    ColorRgba::fromARGB(0xFF808000),
    ColorRgba::fromARGB(0xFFA0A0A0),
  };

  static intx const PALETTE_SIZE = (intx)(sizeof(PALETTE) / sizeof(ColorRgba));
  return PALETTE[abs((n % (PALETTE_SIZE + palette_shift)) % PALETTE_SIZE)];
}

ColorRgba
ShapeRendererImpl::getLabelColor(intx label) const
{
  return getPaletteColor(label);
}

int
ShapeRendererImpl::exec(int argc, char ** argv)
{
  if (!parseArgs(argc, argv))
    return -1;

  if (!has_render_system.exchange(true))
  {
    if (!loadPlugins(argc, argv))
      throw Error("Could not load plugins");
  }

  // Load the mesh
  Model model(show_points & POINTS_PRIMARY);
  if (!loadModel(model, model_paths[0]))
    return -1;

  // Set up framebuffer for offscreen drawing
  int buffer_width   =  antialiasing_level * out_width;
  int buffer_height  =  antialiasing_level * out_height;
  ITexture * color_tex = nullptr;
  ITexture * depth_tex = nullptr;
  IFramebuffer * fb = nullptr;

  TextureOptions tex_opts;
  tex_opts.setInterpolateMode(ITextureOptions::InterpolateMode::NEAREST_NO_MIPMAP);
  color_tex = render_system->createTexture("Color", buffer_width, buffer_height, 1, TextureFormat::RGBA8(),
                                           ITexture::Dimension::DIM_2D, &tex_opts);
  if (!color_tex)
  {
    THEA_ERROR << "Could not create color buffer";
    return -1;
  }

  depth_tex = render_system->createTexture("Depth", buffer_width, buffer_height, 1, TextureFormat::DEPTH16(),
                                           ITexture::Dimension::DIM_2D, &tex_opts);
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

  if (!fb->attach(IFramebuffer::AttachmentPoint::COLOR_0, color_tex)) return -1;
  if (!fb->attach(IFramebuffer::AttachmentPoint::DEPTH,   depth_tex)) return -1;

  // Load matcap material, if any
  matcap_tex = nullptr;
  if (!matcap_path.empty())
  {
    Image matcap_img(matcap_path);
    matcap_tex = render_system->createTexture("Matcap", &matcap_img, TextureFormat::AUTO(), ITexture::Dimension::DIM_2D);
    if (!matcap_tex) return -1;
  }

  // Load 2D texture, if any
  tex2d = nullptr;
  if (!tex2d_path.empty())
  {
    Image tex2d_img(tex2d_path);
    tex2d = render_system->createTexture("Texture2D", &tex2d_img, TextureFormat::AUTO(), ITexture::Dimension::DIM_2D);
    if (!tex2d) return -1;
  }

  // Load 3D texture, if any
  tex3d = nullptr;
  if (!tex3d_path.empty())
  {
    Image tex3d_img(tex3d_path);
    tex3d = render_system->createTexture("Texture3D", &tex3d_img, TextureFormat::AUTO(), ITexture::Dimension::DIM_3D);
    if (tex3d) return -1;
  }

  // Load overlay models
  typedef std::shared_ptr<Model> ModelPtr;
  Array<ModelPtr> overlay_models;
  for (size_t i = 1; i < model_paths.size(); ++i)
  {
    ModelPtr overlay_model(new Model);
    overlay_model->convert_to_points = (show_points & POINTS_OVERLAY);
    if (!loadModel(*overlay_model, model_paths[i]))
      return -1;

    overlay_models.push_back(overlay_model);
  }

  // Set up accumulator to count number of pixels rendered per face/point
  Array<intx> hitcounts;
  if (save_hitcounts)
    hitcounts.resize((size_t)(model.max_id + 1), 0);

  // Do the rendering
  for (size_t v = 0; v < views.size(); ++v)
  {
    try
    {
      // Initialize the camera
      Camera camera = model.fitCamera(transforms[0], views[v], zoom, buffer_width, buffer_height);
      if (print_camera)
      {
        THEA_CONSOLE << "Camera for view " << v << " is: " << camera.toString();
        THEA_CONSOLE << "Viewing matrix (world to camera) for view " << v << " is: "
                     << toString(camera.getWorldToCameraTransform().homogeneous());
        THEA_CONSOLE << "Projection matrix (camera to projection) for view " << v << " is: "
                     << toString(camera.getProjectionTransform());
      }

      // Render the mesh to the offscreen framebuffer
      render_system->pushFramebuffer();
        render_system->setFramebuffer(fb);

        render_system->pushDepthFlags();
        render_system->pushColorFlags();
        render_system->pushShapeFlags();
        render_system->pushTextures();
        render_system->setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
        render_system->setMatrixMode(IRenderSystem::MatrixMode::PROJECTION); render_system->pushMatrix();

          camera.makeCurrent(render_system);
          render_system->setDepthTest(IRenderSystem::DepthTest::LESS);
          render_system->setDepthWrite(true);
          render_system->setColorWrite(true, true, true, true);

          render_system->setClearColor(background_color.data());
          render_system->clear();

          // Draw primary model
          render_system->setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
          {
            render_system->multMatrix(&asLvalue(Math::wrapMatrix(transforms[0])));
            if (!renderModel(model, primary_color)) return -1;
          }
          render_system->setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();

          // Draw overlay models
          for (size_t i = 0; i < overlay_models.size(); ++i)
          {
            Model const & overlay = *overlay_models[i];
            ColorRgba overlay_color = getPaletteColor((intx)i - 1);
            overlay_color.a() = (!color_by_id && !overlay.is_point_cloud ? 0.5f : 1.0f);

            render_system->setPolygonOffset(true, -1.0f);  // make sure overlays appear on top of primary shape

            render_system->setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW); render_system->pushMatrix();
            {
              render_system->multMatrix(&asLvalue(Math::wrapMatrix(transforms[i])));
              if (!renderModel(overlay, overlay_color)) return -1;
            }
            render_system->setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
          }

        render_system->setMatrixMode(IRenderSystem::MatrixMode::PROJECTION); render_system->popMatrix();
        render_system->setMatrixMode(IRenderSystem::MatrixMode::MODELVIEW); render_system->popMatrix();
        render_system->popTextures();
        render_system->popShapeFlags();
        render_system->popColorFlags();
        render_system->popDepthFlags();

        // Grab and save the rendered image
        //
        // We have to explicitly pick a 24-bit (instead of 32-bit) buffer when no background transparency is used, else the
        // FreeImage JPEG codec (and who knows what other codecs) can't save it.
        Image image((background_color.a() <= 0.9999f ? Image::Type::RGBA_8U : Image::Type::RGB_8U),
                    buffer_width, buffer_height);
        color_tex->getImage(&image);

        if (save_hitcounts)
        {
          // We know image type is RGBA_8U because background alpha is set to zero when saving hitcounts

          for (intx r = 0; r < image.getHeight(); ++r)
          {
            uint8 const * pixel = (uint8 const *)image.getScanLine(r);
            for (intx c = 0; c < image.getWidth(); ++c, pixel += 4)
            {
              if (pixel[3] != 0xFF) continue;  // background

              size_t index = (size_t)((((uint32)pixel[2] << 16) | ((uint32)pixel[1] << 8) |  (uint32)pixel[0]) & 0x7FFFFF);
              alwaysAssertM(index < hitcounts.size(), "Rendered face index out of bounds");

              hitcounts[index]++;
            }
          }
        }

        // Grab and save the color image
        if (save_color)
        {
          if (antialiasing_level > 1 && !image.rescale(out_width, out_height, 1, Image::Filter::BICUBIC))
          {
            THEA_ERROR << "Could not rescale image to output dimensions";
            return -1;
          }

          string path = out_path;
          if (views.size() > 1)
          {
            path = FilePath::concat(FilePath::parent(path),
                                    FilePath::baseName(path) + format("_%06ld.", (intx)v) + FilePath::completeExtension(path));
          }

          image.save(path);
        }

        // Grab and save the depth image
        if (save_depth)
        {
          Image depth_image(Image::Type::LUMINANCE_16U, buffer_width, buffer_height);
          depth_tex->getImage(&depth_image);

          if (antialiasing_level > 1 && !depth_image.rescale(out_width, out_height, 1, Image::Filter::BICUBIC))
          {
            THEA_ERROR << "Could not rescale depth image to output dimensions";
            return -1;
          }

          string suffix = (views.size() > 1 ? format("_%06ld.png", (intx)v) : ".png");
          string depth_path = FilePath::concat(FilePath::parent(out_path), FilePath::baseName(out_path) + "_depth" + suffix);

          depth_image.save(depth_path);
        }

      render_system->popFramebuffer();
    }
    THEA_CATCH(return -1;, ERROR, "Could not render view %ld of shape", (intx)v)
  }

  render_system->destroyTexture(tex3d);
  render_system->destroyTexture(tex2d);
  render_system->destroyTexture(matcap_tex);

  if (save_hitcounts)
  {
    string hitcount_path = FilePath::changeCompleteExtension(out_path, "hitcount");
    ofstream out(hitcount_path, ios::binary);
    for (size_t i = 0; i < hitcounts.size(); ++i)
      out << hitcounts[i] << '\n';
  }

  THEA_CONSOLE << "Rendered " << views.size() << " view(s) of the shape";

  return 0;
}

bool
ShapeRendererImpl::usage(int argc, char ** argv)
{
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Usage: " << argv[0] << " [OPTIONS] <mesh> <output-image> <width> <height>";
  THEA_CONSOLE << "";
  THEA_CONSOLE << "Options:";
  THEA_CONSOLE << "  -o <overlay-shape>    (may be mesh or point set)";
  THEA_CONSOLE << "  -s <scope>            (show 'none' | 'primary' | 'overlay' | 'all' shapes";
  THEA_CONSOLE << "                         as points, sampled if necessary)";
  THEA_CONSOLE << "  -t <transform>        (row-major 3x4 or 4x4 matrix, applied to all";
  THEA_CONSOLE << "                         subsequent shapes)";
  THEA_CONSOLE << "  -z <factor>           (zoom factor, default 1)";
  THEA_CONSOLE << "  -v <arg>              (comma-separated 3-vector (viewing direction);";
  THEA_CONSOLE << "                         or 6-vector (direction + eye position);";
  THEA_CONSOLE << "                         or 9-vector (direction + eye + up);";
  THEA_CONSOLE << "                         or string of 3 chars, one for each coordinate,";
  THEA_CONSOLE << "                           each one of +, - or 0;";
  THEA_CONSOLE << "                         or file containing one of the above per line)";
  THEA_CONSOLE << "  -q <view> <proj>      (row-major 3x4 or 4x4 view (world to camera) and";
  THEA_CONSOLE << "                         projection (camera to projection) matrices,";
  THEA_CONSOLE << "                         alternative to -v that allows copy-pasting the";
  THEA_CONSOLE << "                         matrices from -w cam; use '-' to choose a default";
  THEA_CONSOLE << "                         setting for an argument)";
  THEA_CONSOLE << "  -u <dir>              (x, y or z, optionally preceded by + or -)";
  THEA_CONSOLE << "  -p <pixels>           (size of points in pixels -- can be fractional)";
  THEA_CONSOLE << "  -c <argb>             (shape color; or 'id' to color faces by face ID and";
  THEA_CONSOLE << "                         points by point ID; or 'leaf' to assign a random";
  THEA_CONSOLE << "                         color to each leaf submesh; or 'leafname' to derive";
  THEA_CONSOLE << "                         the leaf submesh color from its name, searching for";
  THEA_CONSOLE << "                         'RGB(r,g,b)' if it exists, else hashing the entire";
  THEA_CONSOLE << "                         name to a color)";
  THEA_CONSOLE << "  -i <shift>            (add a shift to how indices are mapped to colors)";
  THEA_CONSOLE << "  -j <argb>             (draw mesh edges in the given color)";
  THEA_CONSOLE << "  -l <path>             (color faces/points by labels from <path>)";
  THEA_CONSOLE << "  -f <path>             (color vertices/points by features from <path>)";
  THEA_CONSOLE << "  -2 <path>             (color mesh by a 2D texture)";
  THEA_CONSOLE << "  -3 <path>             (color mesh by a 3D texture)";
  THEA_CONSOLE << "  -e                    (accentuate features)";
  THEA_CONSOLE << "  -m <name>             (render submesh with the given name as white, the";
  THEA_CONSOLE << "                         rest of the shape as black, without any shading)";
  THEA_CONSOLE << "  -h <name> <argb>      (highlight submesh with the given name in the given";
  THEA_CONSOLE << "                         color, with shading)";
  THEA_CONSOLE << "  -k 'ka kd ks ksp'     (Phong material coefficients)";
  THEA_CONSOLE << "  -k <path>             (texture to be used as matcap material)";
  THEA_CONSOLE << "  -b <argb>             (background color)";
  THEA_CONSOLE << "  -a N                  (enable NxN antialiasing: 2 is normal, 4 is very";
  THEA_CONSOLE << "                         high quality)";
  THEA_CONSOLE << "  -0                    (flat shading)";
  THEA_CONSOLE << "  -1                    (one-sided lighting)";
  THEA_CONSOLE << "  -g <dir>              (light direction, as 'x y z')";
  THEA_CONSOLE << "  --toon                (render with toon/cel shading)";
  THEA_CONSOLE << "  -n                    (output a normal map instead of colors)";
  THEA_CONSOLE << "  -d                    (also save the depth image)";
  THEA_CONSOLE << "  -#                    (also save a text file with extension '.hitcount'";
  THEA_CONSOLE << "                         containing the number of pixels rendered per";
  THEA_CONSOLE << "                         face/point)";
  THEA_CONSOLE << "  -x                    (suppress color image output)";
  THEA_CONSOLE << "  -w cam                (print the camera parameters)";
  THEA_CONSOLE << "";

  return false;
}

bool
ShapeRendererImpl::parseTransform(string const & s, Matrix4 & m)
{
  auto t = trimWhitespace(s);
  if (t.size() < 2) { THEA_ERROR << "Could not parse transform from string: " << s; return false; }
  if (t.front() == '[')
  {
    if (t.back() == ']') { t = t.substr(1, t.length() - 2); }
    else                 { THEA_ERROR << "Could not parse transform from string: " << s; return false; }
  }

  m.setIdentity();
  istringstream t_in(t);
  for (int i = 0, c = 0; i < m.rows() && c != char_traits<char>::eof(); ++i)
    for (int j = 0; j < m.cols(); ++j)
    {
      if (!(t_in >> m(i, j))) { THEA_ERROR << "Could not parse transform from string: " << s; return false; }

      do
      {
        c = t_in.get();
        if (!t_in && (i < m.rows() - 2 || j != m.cols() - 1))
        { THEA_ERROR << "Could not parse transform from string: " << s; return false; }
      } while (t_in && isspace(c));

      if (t_in && c != ',' && (j != m.cols() - 1 || c != ';')) t_in.unget();
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

  if (view.dir.squaredNorm() <= 1e-10)
  {
    if (!silent) THEA_ERROR << "View direction is zero vector";
    return false;
  }

  view.dir.normalize();

  if (has_up)
    view.up = view_up;
  else if (s == "0-0")
    view.up = -Vector3::UnitZ();
  else if (s == "0+0")
    view.up = Vector3::UnitZ();
  else
    view.up = Vector3::UnitY();

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
  if (view.dir.squaredNorm() <= 1e-10)
  {
    if (!silent) THEA_ERROR << "View direction is zero vector";
    return false;
  }

  view.dir.normalize();

  if (num_params >= 6)
  {
    view.has_eye = true;
    view.eye = Vector3(ex, ey, ez);
  }

  if (num_params >= 9)
  {
    view.up = Vector3(ux, uy, uz);

    if (view.up.squaredNorm() <= 1e-10)
    {
      if (!silent) THEA_ERROR << "View up is zero vector";
      return false;
    }

    view.up.normalize();
  }
  else if (has_up)
    view.up = view_up;
  else
  {
    Real d = view.dir.dot(Vector3::UnitY());
    if (Math::fuzzyEq(d, (Real)-1))
      view.up = -Vector3::UnitZ();
    else if (Math::fuzzyEq(d, (Real)1))
      view.up = Vector3::UnitZ();
    else
      view.up = Vector3::UnitY();
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
    up = Vector3::UnitX();
  else if (c == 'y' || c == 'Y')
    up = Vector3::UnitY();
  else if (c == 'z' || c == 'Z')
    up = Vector3::UnitZ();
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
ShapeRendererImpl::parseARGB(string const & s, uint32 & argb)
{
  std::stringstream ss;
  ss << std::hex << s;

  if (!(ss >> argb))
  {
    THEA_ERROR << "Could not parse ARGB color '" << s << '\'';
    return false;
  }

  return true;
}

bool
ShapeRendererImpl::parseColor(string const & s, ColorRgba & c)
{
  uint32 argb;
  if (!parseARGB(s, argb))
    return false;

  c = ColorRgba::fromARGB(argb);

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
    return usage(argc, argv);

  resetArgs();

  Matrix4 current_transform = Matrix4::Identity();

  argv++;
  argc--;
  int pos = 0;

  bool explicit_primary_color = false;

  while (argc > 0)
  {
    string arg = *argv;
    argv++; argc--;

    if (arg.empty())
      continue;

    if (arg[0] == '-')
    {
      if (arg == "-o")
      {
        if (argc < 1) { THEA_ERROR << "-o: Overlay path not specified"; return false; }
        model_paths.push_back(*argv);
        transforms.push_back(current_transform);
        argv++; argc--;
      }
      else if (arg == "-s")
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

        argv++; argc--;
      }
      else if (arg == "-t")
      {
        if (argc < 1) { THEA_ERROR << "-t: Transform not specified"; return false; }
        if (!parseTransform(*argv, current_transform)) return false;
        argv++; argc--;
      }
      else if (arg == "-z")
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
        argv++; argc--;
      }
      else if (arg == "-v")
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
            THEA_ERROR << "Invalid or unparseable view parameters '" << *argv << '\'';
            return false;
          }
        }

        argv++; argc--;
      }
      else if (arg == "-q")
      {
        if (argc < 2) { THEA_ERROR << "-q: Camera matrices not specified"; return false; }
        View view;

        string view_str = *(argv++); argc--;
        if (view_str.empty() || view_str == "-")
          view.has_view_matrix = false;
        else
        {
          if (parseTransform(view_str, view.view_matrix))
            view.has_view_matrix = true;
          else
            return false;
        }

        string proj_str = *(argv++); argc--;
        if (proj_str.empty() || proj_str == "-")
          view.has_proj_matrix = false;
        else
        {
          if (parseTransform(proj_str, view.proj_matrix))
            view.has_proj_matrix = true;
          else
            return false;
        }

        views.push_back(view);
      }
      else if (arg == "-u")
      {
        if (argc < 1) { THEA_ERROR << "-u: Up direction not specified"; return false; }
        if (!parseViewUp(*argv, view_up)) return false;
        has_up = true;
        argv++; argc--;
      }
      else if (arg == "-p")
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
        argv++; argc--;
      }
      else if (arg == "-c")
      {
        if (argc < 1) { THEA_ERROR << "-c: Color not specified"; return false; }

        color_by_id = false;
        color_by_leaf = false;
        color_by_leafname = false;

        string arg_lc = toLower(trimWhitespace(*argv));
        if (arg_lc == "id")
          color_by_id = true;
        else if (arg_lc == "leaf")
          color_by_leaf = true;
        else if (arg_lc == "leafname")
          color_by_leafname = true;
        else
        {
          if (!parseColor(*argv, primary_color))
            return false;
        }

        explicit_primary_color = true;

        argv++; argc--;
      }
      else if (arg == "-j")
      {
        if (argc < 1) { THEA_ERROR << "-j: Edge color not specified"; return false; }
        if (!parseColor(*argv, edge_color)) return false;
        show_edges = true;
        argv++; argc--;
      }
      else if (arg == "-l")
      {
        if (argc < 1) { THEA_ERROR << "-l: Labels not specified"; return false; }

        color_by_label = true;
        labels_path = *argv;

        argv++; argc--;
      }
      else if (arg == "-f")
      {
        if (argc < 1) { THEA_ERROR << "-f: Features not specified"; return false; }

        color_by_features = true;
        features_path = *argv;

        argv++; argc--;
      }
      else if (arg == "-e")
      {
        accentuate_features = true;
      }
      else if (arg == "-2")
      {
        if (argc < 1) { THEA_ERROR << "-2: 2D texture image not specified"; return false; }

        color_by_tex2d = true;
        tex2d_path = *argv;

        argv++; argc--;
      }
      else if (arg == "-3")
      {
        if (argc < 1) { THEA_ERROR << "-3: 3D texture image not specified"; return false; }

        color_by_tex3d = true;
        tex3d_path = *argv;

        argv++; argc--;
      }
      else if (arg == "-m")
      {
        if (argc < 1) { THEA_ERROR << "-m: Selected mesh not specified"; return false; }

        selected_binary_mask = true;
        selected_mesh = toLower(*argv);
        selected_color = ColorRgba(1, 1, 1, 1);

        argv++; argc--;
      }
      else if (arg == "-h")
      {
        if (argc < 2) { THEA_ERROR << "-h: Selected mesh or color not specified"; return false; }

        selected_binary_mask = false;
        selected_mesh = toLower(*argv);
        argv++; argc--;
        if (!parseColor(*argv, selected_color))
          return false;

        argv++; argc--;
      }
      else if (arg == "-k")
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

        argv++; argc--;
      }
      else if (arg == "-b")
      {
        if (argc < 1) { THEA_ERROR << "-b: Background color not specified"; return false; }
        if (!parseColor(*argv, background_color)) return false;
        argv++; argc--;
      }
      else if (arg == "-a")
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
        argv++; argc--;
      }
      else if (arg == "-0")
      {
        flat = true;
      }
      else if (arg == "-1")
      {
        two_sided = false;
      }
      else if (arg == "-g")
      {
        Vector4 d;
        if (parseVector(*argv, d) != 3) { THEA_ERROR << "-g: Could not parse light direction"; return false; }
        light_dir = d.head<3>().stableNormalized();

        argv++; argc--;
      }
      else if (arg == "-n")
      {
        out_normals = true;
      }
      else if (arg == "--toon")
      {
        toon_shading = true;
      }
      else if (arg == "-d")
      {
        save_depth = true;
      }
      else if (arg == "-#")
      {
        save_hitcounts = true;
      }
      else if (arg == "-x")
      {
        save_color = false;
      }
      else if (arg == "-w")
      {
        if (argc < 1) { THEA_ERROR << "-w: Field to write not specified"; return false; }
        if (string(*argv) == "cam")
          print_camera = true;
        else
        {
          THEA_ERROR << "Unknown field to write: '" << *argv << '\'';
          return false;
        }
        argv++; argc--;
      }
      else if (arg == "-i")
      {
        if (argc < 1) { THEA_ERROR << "-i: Palette shift not specified"; return false; }
        if (sscanf(*argv, " %d", &palette_shift) != 1)
        {
          THEA_ERROR << "Could not parse palette shift: '" << *argv << '\'';
          return false;
        }
        argv++; argc--;
      }
      else
      {
        THEA_ERROR << "Unknown option: " << arg;
        return usage(argc, argv);
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
    return usage(argc, argv);
  }

  // If no primary color was explicitly specified, set it to white by default if we're also going to multiply it by the matcap
  // or 2D/3D texture color. The assumption is that the latter colors should be presented accurately unless the user explicitly
  // modulates them by an additional color.
  if (!explicit_primary_color && (!matcap_path.empty() || color_by_tex2d || color_by_tex3d))
    primary_color = ColorRgba(1, 1, 1, 1);

  if (views.empty())
  {
    View view;
    if (has_up)
    {
      view.up = view_up;
      int axis = Math::maxAbsAxis(view_up);
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

  if (save_hitcounts)
  {
    if (!color_by_id) { THEA_ERROR << "Saving hitcounts requires coloring elements by id"; return false; }
    if (antialiasing_level > 1) { THEA_ERROR << "Saving hitcounts requires no antialiasing"; return false; }
    if (model_paths.size() > 1) { THEA_ERROR << "Saving hitcounts requires no overlay models"; return false; }

    THEA_WARNING << "Background alpha overridden to zero to accurately measure hitcounts";
    background_color.a() = 0;
  }

  return true;
}

bool
flattenFaces(Mesh & mesh)
{
  mesh.isolateTriangles();
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

ColorRgba8
indexToColor(uint32 index, bool is_point)
{
  ColorRgba8 color((uint8)((index      ) & 0xFF),
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
  Array<int> const * labels;

  FaceColorizer(ShapeRendererImpl const * parent_, Array<int> const * labels_ = nullptr)
  : parent(parent_), labels(labels_)
  {}

  bool operator()(Mesh & mesh)
  {
    if (parent->flat || !(parent->color_by_leaf || parent->color_by_leafname))
      mesh.isolateTriangles();

    if (parent->flat || labels || parent->color_by_leaf || parent->color_by_leafname)
      mesh.computeAveragedVertexNormals();

    mesh.addColors();

    ColorRgba8 color;
    if (parent->color_by_leaf)
      color = parent->getPaletteColor(Random::common().integer());
    else if (parent->color_by_leafname)
    {
      string const & name = mesh.getName();
      auto loc = name.find("RGB(");
      float r, g, b;
      if (loc != string::npos && sscanf(name.substr(loc).c_str(), "RGB( %f , %f , %f )", &r, &g, &b) == 3)
        color = ColorRgba8(ColorRgba(r, g, b, 1.0));
      else
        color = parent->getLabelColor((intx)labelHash(mesh.getName()));
    }

    Mesh::IndexArray const & tris = mesh.getTriangleIndices();
    for (size_t i = 0; i < tris.size(); i += 3)
    {
      if (!parent->color_by_leaf && !parent->color_by_leafname)
      {
        auto id = mesh.getTriangleSourceFaceIndex((intx)i / 3);
        if (labels)
        {
          if ((size_t)id >= labels->size())
            color = ColorRgba8(255, 0, 0, 255);
          else
            color = parent->getLabelColor((*labels)[(size_t)id]);
        }
        else
          color = indexToColor((uint32)id, false);
      }

      mesh.setColor((intx)tris[i    ], color);
      mesh.setColor((intx)tris[i + 1], color);
      mesh.setColor((intx)tris[i + 2], color);
    }

    return false;
  }
};

struct FaceCounter
{
  intx num_faces;

  FaceCounter() : num_faces(0) {}
  bool operator()(Mesh & mesh) { num_faces += mesh.numFaces(); return false; }
};

struct MeshToFaceLabels
{
  UnorderedMap<intx, int> const & input_face_labels;
  Array<int> & all_face_labels;

  MeshToFaceLabels(UnorderedMap<intx, int> const & input_face_labels_, Array<int> & all_face_labels_)
  : input_face_labels(input_face_labels_), all_face_labels(all_face_labels_)
  {
    std::fill(all_face_labels_.begin(), all_face_labels_.end(), -1);
  }

  bool operator()(Mesh & mesh)
  {
    intx nt = mesh.numTriangles();
    for (intx i = 0; i < nt; ++i)
    {
      auto loc = input_face_labels.find(mesh.getTriangleSourceFaceIndex(i));
      if (loc != input_face_labels.end())
      {
        for (intx j = 0; j < nt; ++j)
          all_face_labels[(size_t)mesh.getTriangleSourceFaceIndex(j)] = loc->second;

        break;  // found the label for this mesh
      }
    }

    return false;
  }
};

bool
ShapeRendererImpl::loadLabels(Model & model)
{
  if (labels_path.empty())
  {
    THEA_ERROR << "No labels file specified";
    return false;
  }

  Array<int> labels;
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

    typedef UnorderedMap<string, int> LabelIndexMap;
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
      intx index;
      while (line_in >> index)
      {
        index--;  // first face has index 1 in this format
        if (index < 0)
        {
          THEA_ERROR << "Index " << index << " of element to be labeled is out of bounds";
          return false;
        }

        if (index >= (intx)labels.size())
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

    FaceCounter face_counter;
    model.mesh_group.forEachMeshUntil(std::ref(face_counter));

    typedef UnorderedMap<intx, int> FaceLabelMap;
    FaceLabelMap input_face_labels;

    typedef UnorderedMap<string, int> LabelIndexMap;
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
      intx rep_index;
      while (line_in >> rep_index)
      {
        if (rep_index < 0 || rep_index >= face_counter.num_faces)
        {
          THEA_ERROR << "Index " << rep_index << " of element to be labeled is out of bounds";
          return false;
        }

        input_face_labels[rep_index] = label_index;
      }
    }

    // Map from mesh labels to face labels
    labels.resize((size_t)face_counter.num_faces);
    MeshToFaceLabels m2f(input_face_labels, labels);
    model.mesh_group.forEachMeshUntil(m2f);
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

    typedef UnorderedMap<string, int> LabelIndexMap;
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
      model.point_colors[i] = (labels[i] < 0 ? ColorRgba(1, 0, 0, 1) : getLabelColor(labels[i]));
  }
  else
  {
    FaceColorizer label_colorizer(this, &labels);
    model.mesh_group.forEachMeshUntil(label_colorizer);
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

bool
parseInt(std::string const & s, intx & n)
{
  long i;
  char c;
  if (std::sscanf(s.c_str(), " %ld %c", &i, &c) != 1)
  { THEA_ERROR << "String '" << s << "' is not an integer"; return false; }

  n = (intx)i;
  return true;
}

bool
parseReal(std::string const & s, Real & x)
{
  long double d;
  char c;
  if (std::sscanf(s.c_str(), " %Lf %c", &d, &c) != 1)
  { THEA_ERROR << "String '" << s << "' is not a real number"; return false; }

  x = (Real)d;
  return true;
}

ColorRgb
featToColor(Real f0, Real const * f1, Real const * f2)
{
  if (!f2)
  {
    if (!f1)
      return ColorRgb::jetColorMap(0.2 + 0.6 * f0);
    else
      return ColorRgb(f0, *f1, 1.0f);
  }
  else
    return ColorRgb(f0, *f1, *f2);
}

typedef BvhN<Vector3, 3> PointBvh;
typedef UnorderedMap<intx, size_t> IndexMap;

struct VertexColorizer
{
  VertexColorizer(PointBvh const * fbvh_, IndexMap const * vertex_features_,
                  Real const * feat_vals0_, Real const * feat_vals1_, Real const * feat_vals2_)
  : fbvh(fbvh_), vertex_features(vertex_features_), feat_vals0(feat_vals0_), feat_vals1(feat_vals1_), feat_vals2(feat_vals2_) {}

  bool operator()(Mesh & mesh)
  {
    static int const MAX_NBRS = 8;  // ok to have a few neighbors for output quality -- this is an offline renderer
    Real scale = std::max(0.2f * fbvh->getBounds().getExtent().norm(), (Real)1.0e-8);
    Real scale2 = scale * scale;

    Mesh::VertexArray const & vertices = mesh.getVertices();
    BoundedSortedArrayN<MAX_NBRS, PointBvh::NeighborPair> nbrs;
    for (size_t i = 0; i < vertices.size(); ++i)
    {
      if (vertex_features->find(mesh.getVertexSourceIndex((intx)i)) != vertex_features->end())
        continue;

      nbrs.clear();
      intx num_nbrs = fbvh->kClosestPairs<Algorithms::MetricL2>(vertices[i], nbrs, 2 * scale);
      if (num_nbrs <= 0)
        num_nbrs = fbvh->kClosestPairs<Algorithms::MetricL2>(vertices[i], nbrs);

      if (num_nbrs > 0)
      {
        ColorRgb c(0, 0, 0);
        double sum_weights = 0;
        for (size_t j = 0; j < nbrs.size(); ++j)
        {
          double dist = nbrs[j].getDistance<MetricL2>();
          double weight = Math::fastMinusExp(dist * dist / scale2);
          intx nn_index = nbrs[j].getTargetIndex();
          sum_weights += weight;
          c += weight * featToColor(feat_vals0[nn_index],
                                    (feat_vals1 ? &feat_vals1[nn_index] : nullptr),
                                    (feat_vals2 ? &feat_vals2[nn_index] : nullptr));
        }

        mesh.setColor((intx)i, sum_weights > 0 ? c / sum_weights : c);
      }
      else
      {
        THEA_WARNING << "No nearest neighbor found!";
        mesh.setColor((intx)i, ColorRgb(1, 1, 1));
      }
    }

    return false;
  }

  PointBvh const * fbvh;
  IndexMap const * vertex_features;
  Real const * feat_vals0;
  Real const * feat_vals1;
  Real const * feat_vals2;
};

bool
ShapeRendererImpl::loadFeatures(Model & model)
{
  IndexMap vertex_features;
  Array<Vector3> feat_pts;
  Array< Array<Real> > feat_vals(1);
  try
  {
    std::ifstream in(features_path.c_str());
    if (!in)
      throw Error("Couldn't open file");

    std::string line, fields[3];
    Vector3 p;
    long double f;  // to handle underflow problems with smaller representations
    while (getNextNonBlankLine(in, line))
    {
      std::istringstream line_in(line);
      if (!(line_in >> fields[0] >> fields[1] >> fields[2] >> f))
        throw Error("Couldn't read feature from line '" + line + '\'');

      if (fields[0] == "v")  // direct vertex reference
      {
        intx vindex;
        if (!parseInt(fields[1], vindex))  // fields[2] is ignored, currently
          throw Error("Couldn't read feature from line '" + line + '\'');

        vertex_features[vindex] = feat_pts.size();
        p = Vector3::Zero();  // will be set below
      }
      else
        if (!parseReal(fields[0], p[0]) || !parseReal(fields[1], p[1]) || !parseReal(fields[2], p[2]))
          throw Error("Couldn't read feature from line '" + line + '\'');

      feat_pts.push_back(p);
      feat_vals[0].push_back(f);

      if (feat_pts.size() == 1)
      {
        while (line_in >> f)
        {
          feat_vals.push_back(Array<Real>());
          feat_vals.back().push_back((Real)Math::clamp(f, -std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()));
        }
      }
      else
      {
        for (size_t i = 1; i < feat_vals.size(); ++i)
        {
          if (!(line_in >> f))
            throw Error("Couldn't read feature from line '" + line + '\'');

          feat_vals[i].push_back((Real)Math::clamp(f, -std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max()));
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
          Array<Real> sorted = feat_vals[i];
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
  THEA_CATCH(return false;, WARNING, "Couldn't load model features from '%s'", features_path.c_str())

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
                                          (feat_vals.size() > 1 ? &feat_vals[1][i] : nullptr),
                                          (feat_vals.size() > 2 ? &feat_vals[2][i] : nullptr));
    }
  }
  else
  {
    model.mesh_group.forEachMeshUntil([&](Mesh & mesh) {
      mesh.addColors();
      auto const & vertices = mesh.getVertices();
      for (size_t i = 0; i < vertices.size(); ++i)
      {
        auto existing = vertex_features.find(mesh.getVertexSourceIndex((intx)i));
        if (existing != vertex_features.end())
          mesh.setColor((intx)i, featToColor(feat_vals[0][existing->second],
                                             (feat_vals.size() > 1 ? &feat_vals[1][existing->second] : nullptr),
                                             (feat_vals.size() > 2 ? &feat_vals[2][existing->second] : nullptr)));
      }

      return false;
    });

    PointBvh fbvh(feat_pts.begin(), feat_pts.end());
    VertexColorizer visitor(&fbvh, &vertex_features,
                            &feat_vals[0][0],
                            feat_vals.size() > 1 ? &feat_vals[1][0] : nullptr,
                            feat_vals.size() > 2 ? &feat_vals[2][0] : nullptr);
    model.mesh_group.forEachMeshUntil(visitor);
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

    ColorRgba8 color;
    if (selected_binary_mask)
      color = (mesh_id == 0 ? ColorRgba8(0, 0, 0, 255) : ColorRgba8(255, 255, 255, 255));
    else
      color = (mesh_id == 0 ? ColorRgba8(primary_color) : ColorRgba8(selected_color));

    if (selected_binary_mask) mesh.isolateTriangles();
    mesh.addColors();

    Mesh::IndexArray const & tris = mesh.getTriangleIndices();
    for (size_t i = 0; i < tris.size(); i += 3)
    {
      mesh.setColor((intx)tris[i    ], color);
      mesh.setColor((intx)tris[i + 1], color);
      mesh.setColor((intx)tris[i + 2], color);
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

struct MeshReadCallback : public MeshCodec<Mesh>::ReadCallback
{
  intx max_id;

  MeshReadCallback() : max_id(-1) {}
  void faceRead(Mesh * mesh, intx index, Mesh::FaceHandle face) { max_id = std::max(max_id, index); }
};

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
    model.max_id = (intx)model.points.size() - 1;

    if (color_by_label)
    {
      if (!loadLabels(model))
        return false;
    }
    else if (color_by_features)
    {
      if (!loadFeatures(model))
        return false;
    }
    else if (color_by_id)
    {
      model.point_colors.resize(model.points.size());
      for (size_t i = 0; i < model.points.size(); ++i)
        model.point_colors[i] = indexToColor((uint32)i, true);
    }

    THEA_CONSOLE << "Read " << model.points.size() << " points from '" << path << '\'';
  }
  else
  {
    try
    {
      MeshReadCallback callback;
      model.mesh_group.load(path, CodecAuto(), &callback);
      model.max_id = callback.max_id;
    }
    THEA_CATCH(return false;, ERROR, "Could not load model from '%s'", path.c_str())

    if (model.convert_to_points)
    {
      typedef MeshSampler<Mesh>::Triangle MeshTriangle;
      MeshSampler<Mesh> sampler(model.mesh_group);
      Array<MeshTriangle const *> sampled_tris;

      double out_scale = min(out_width, out_height);
      intx max_samples = (intx)ceil(10 * out_scale);
      intx npoints = sampler.sampleEvenlyByArea(max_samples, model.points, nullptr, (color_by_id ? &sampled_tris : nullptr));

      THEA_CONSOLE << "Sampled " << npoints << " points from '" << path << '\'';

      model.mesh_group.clear();
      model.is_point_cloud = true;

      if (color_by_id)
      {
        alwaysAssertM(sampled_tris.size() == model.points.size(), "Number of source triangles != number of sampled points");

        model.point_colors.resize(model.points.size());
        for (size_t i = 0; i < model.points.size(); ++i)
        {
          MeshTriangle::VertexTriple const & tverts = sampled_tris[i]->getVertices();
          intx face_id = tverts.getMesh()->getTriangleSourceFaceIndex(tverts.getMeshTriangleIndex());
          model.point_colors[i] = indexToColor((uint32)face_id, true);
        }
      }
    }
    else
    {
      bool needs_normals = true;

      if (color_by_id || color_by_leaf || color_by_leafname)
      {
        FaceColorizer colorizer(this);
        model.mesh_group.forEachMeshUntil(colorizer);
        needs_normals = false;  // FaceColorizer has already computed as needed
      }
      else if (color_by_label)
      {
        if (!loadLabels(model))
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
        needs_normals = !selected_binary_mask;
      }

      if (needs_normals)
      {
        if (flat)
          model.mesh_group.forEachMeshUntil(flattenFaces);
        else
          model.mesh_group.forEachMeshUntil(averageNormals);
      }
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
        Real sqdist = (Math::hmul(transform, vertices[i]) - center).squaredNorm();
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

  center = Math::hmul(transform, center);

  Real radius = 0;
  if (model.is_point_cloud)
  {
    for (size_t i = 0; i < model.points.size(); ++i)
    {
      Real sqdist = (Math::hmul(transform, model.points[i]) - center).squaredNorm();
      if (i == 0 || sqdist > radius)
        radius = sqdist;
    }

    radius = sqrt(radius);
  }
  else
  {
    FarthestPoint fp(center, transform);
    model.mesh_group.forEachMeshUntil(std::ref(fp));
    radius = fp.getFarthestDistance();
  }

  return Ball3(center, radius);
}

Camera
Model::fitCamera(Matrix4 const & transform, View const & view, Real zoom, int width, int height)
{
  Ball3 bsphere;
  Vector3 center;
  Real diameter;
  if (!view.has_view_matrix || !view.has_proj_matrix)
  {
    bsphere = modelBSphere(*this, transform);
    center = bsphere.getCenter();
    diameter = bsphere.getDiameter();
  }

  CoordinateFrame3 cframe;
  if (view.has_view_matrix)
  {
    Matrix3 rot = Math::orthonormalBasis(Matrix3(view.view_matrix.topLeftCorner<3, 3>()));  // world to camera coordinates
    Vector3 trn = view.view_matrix.topRightCorner<3, 1>();
    cframe = CoordinateFrame3::_fromAffine(AffineTransform3(rot, trn)).inverse();  // camera to world
  }
  else
  {
    // Make absolutely sure these are unit vectors
    Vector3 dir  =  view.dir.normalized();
    Vector3 up   =  view.up.normalized();
    Vector3 eye;

    if (view.has_eye)
      eye = view.eye;
    else
    {
      Real camera_separation = diameter > 1.0e-10f ? 2.1 * diameter : 1.0e-10f;
      eye = center - camera_separation * dir;
    }

    cframe = CoordinateFrame3::fromViewFrame(eye, eye + dir, up);
  }

  Camera camera;
  if (view.has_proj_matrix)
  {
    camera.setFrame(cframe);
    if (!camera.setProjection(view.proj_matrix))
      throw Error("Could not infer projection parameters from projection matrix");

    alwaysAssertM(camera.getProjectionTransform().isApprox(view.proj_matrix),
                  "Inferred projection parameters do not reconstruct projection matrix");
  }
  else
  {
    // Projection
    Vector3 eye = cframe.getTranslation();
    Vector3 dir = cframe.lookVector();
    Real center_distance = std::max((center - eye).dot(dir), (Real)0);
    Real near_dist = std::max(center_distance - 1.1f * diameter, 0.01f * diameter);
    Real far_dist  = center_distance + 2 * diameter;

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

    hw = (hw / zoom) * (0.5f * near_dist);
    hh = (hh / zoom) * (0.5f * near_dist);

    camera.set(cframe,
               Camera::ProjectionType::PERSPECTIVE, -hw, hw, -hh, hh, near_dist, far_dist, Camera::ProjectedYDirection::UP);
  }

  return camera;
}

bool
initPointShader(IShader & shader)
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
    shader.attachModuleFromString(IShader::ModuleType::VERTEX, VERTEX_SHADER.c_str());
    shader.attachModuleFromString(IShader::ModuleType::FRAGMENT, FRAGMENT_SHADER.c_str());
  }
  THEA_CATCH(return false;, ERROR, "%s", "Could not attach point shader module")

  return true;
}

bool
initMeshShader(IShader & shader, Vector4 const & material, Vector3 light_dir, bool two_sided = true, bool out_normals = false,
               bool toon_shading = false, ITexture * matcap_tex = nullptr, ITexture * tex2d = nullptr,
               ITexture * tex3d = nullptr, AxisAlignedBox3 const & bbox = AxisAlignedBox3())
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
"\n"
"  gl_TexCoord[0] = gl_MultiTexCoord0;\n"
"}\n";

  static string const FRAGMENT_SHADER_HEADER_1 =
"varying vec3 src_pos;  // position in mesh coordinates\n"
"varying vec3 normal;  // normal in camera space\n";

  static string const FRAGMENT_SHADER_HEADER_PHONG =
"uniform float two_sided;\n"
"uniform vec3 ambient_color;\n"
"uniform vec3 light_dir;  // must be specified in camera space, pointing towards object\n"
"uniform vec3 light_color;\n"
"uniform vec4 material;  // [ka, kl, <ignored>, <ignored>]\n";

  static string const FRAGMENT_SHADER_HEADER_TOON =
"uniform float two_sided;\n"
"uniform vec3 light_dir;  // must be specified in camera space, pointing towards object\n";

  static string const FRAGMENT_SHADER_HEADER_MATCAP =
"uniform float two_sided;\n"
"uniform sampler2D matcap_tex;\n"
"\n"
"float lightness(vec3 color)\n"
"{\n"
"  return 0.5 * (max(max(color.r, color.g), color.b) + min(min(color.r, color.g), color.b));\n"
"}\n";

  static string const FRAGMENT_SHADER_HEADER_TEX2D =
"uniform sampler2D tex2d;\n";

  static string const FRAGMENT_SHADER_HEADER_TEX3D =
"uniform sampler3D tex3d;\n"
"uniform vec3 bbox_lo;\n"
"uniform vec3 bbox_hi;\n";

  static string const FRAGMENT_SHADER_BODY_1 =
"void main()\n"
"{\n"
"  vec3 N = normalize(normal);\n"
"  vec4 color = gl_Color;\n";

  static string const FRAGMENT_SHADER_BODY_TEX2D =
"  vec4 tex2d_color = texture2D(tex2d, gl_TexCoord[0].st);\n"
"  color = vec4(color.rgb * tex2d_color.rgb, tex2d_color.a);\n";

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

  static string const FRAGMENT_SHADER_BODY_TOON =
"  vec3 L = normalize(light_dir);\n"
"  float NdL = -dot(N, L);\n"
"  if (NdL >= -two_sided)\n"
"  {\n"
"    float intensity = abs(NdL);\n"
"    if (intensity > 0.95)\n"
"      gl_FragColor = vec4(color.rgb, 1.0);\n"
"    else if (intensity > 0.5)\n"
"      gl_FragColor = vec4(0.6 * color.rgb, 1.0);\n"
"    else if (intensity > 0.25)\n"
"      gl_FragColor = vec4(0.4 * color.rgb, 1.0);\n"
"    else\n"
"      gl_FragColor = vec4(0.2 * color.rgb, 1.0);\n"
"  }\n"
"  else\n"
"    gl_FragColor = vec4(0.0, 0.0, 0.0, 1.0);\n"
"}\n";

  static string const FRAGMENT_SHADER_BODY_MATCAP =
"  vec2 matcap_uv = (two_sided < 0.5 && N.z < 0.0 ? normalize(N.xy) : N.xy);\n"
"  vec4 matcap_color = texture2D(matcap_tex, 0.495 * matcap_uv + 0.5);\n"
"  float lite = lightness(matcap_color.rgb);\n"
"  float highlight_blend = lite * lite * lite * lite;\n"
"  gl_FragColor = vec4(mix(color.rgb, vec3(1.0, 1.0, 1.0), highlight_blend) * matcap_color.rgb, color.a);\n"
"}\n";

  static string const FRAGMENT_SHADER_BODY_NORMAL =
"  gl_FragColor = vec4(0.5 * N + 0.5, 1.0);\n"
"}\n";

  string fragment_shader = FRAGMENT_SHADER_HEADER_1;
  if (matcap_tex)
    fragment_shader += FRAGMENT_SHADER_HEADER_MATCAP;
  else if (toon_shading)
    fragment_shader += FRAGMENT_SHADER_HEADER_TOON;
  else if (!out_normals)
    fragment_shader += FRAGMENT_SHADER_HEADER_PHONG;

  if (tex2d) fragment_shader += FRAGMENT_SHADER_HEADER_TEX2D;
  if (tex3d) fragment_shader += FRAGMENT_SHADER_HEADER_TEX3D;

  fragment_shader += FRAGMENT_SHADER_BODY_1;
  if (tex2d) fragment_shader += FRAGMENT_SHADER_BODY_TEX2D;
  if (tex3d) fragment_shader += FRAGMENT_SHADER_BODY_TEX3D;
  if (matcap_tex)
    fragment_shader += FRAGMENT_SHADER_BODY_MATCAP;
  else if (out_normals)
    fragment_shader += FRAGMENT_SHADER_BODY_NORMAL;
  else if (toon_shading)
    fragment_shader += FRAGMENT_SHADER_BODY_TOON;
  else
    fragment_shader += FRAGMENT_SHADER_BODY_PHONG;

  if (!shader.attachModuleFromString(IShader::ModuleType::VERTEX, VERTEX_SHADER.c_str())
   || !shader.attachModuleFromString(IShader::ModuleType::FRAGMENT, fragment_shader.c_str()))
    return false;

  if (!out_normals)
    shader.setUniform("two_sided", (two_sided ? 1.0f : 0.0f));

  if (matcap_tex)
    shader.setUniform("matcap_tex", matcap_tex);
  else if (toon_shading)
    shader.setUniform("light_dir", &asLvalue(Math::wrapMatrix(light_dir)));
  else if (!out_normals)
  {
    Vector3 light_color(1, 1, 1);
    Vector3 ambient_color(1, 1, 1);

    shader.setUniform("light_dir", &asLvalue(Math::wrapMatrix(light_dir)));
    shader.setUniform("light_color", &asLvalue(Math::wrapMatrix(light_color)));
    shader.setUniform("ambient_color", &asLvalue(Math::wrapMatrix(ambient_color)));
    shader.setUniform("material", &asLvalue(Math::wrapMatrix(const_cast<Vector4 &>(material))));
  }

  if (tex2d)
    shader.setUniform("tex2d", tex2d);

  if (tex3d)
  {
    shader.setUniform("tex3d", tex3d);

    // Tightly fit the shape within the texture volume, without distorting the shape or the cubical voxels
    int width   =  tex3d->getWidth();
    int height  =  tex3d->getHeight();
    int depth   =  tex3d->getDepth();
    Vector3 box_ext = bbox.getExtent();
    Vector3 vol_ext(width, depth, height);
    Vector3 dim_ratio = vol_ext.cwiseQuotient(box_ext);
    int fit_axis = (int)Math::minAxis(dim_ratio);
    box_ext = vol_ext / dim_ratio[fit_axis];
    Vector3 center = bbox.getCenter();

    Vector3 bbox_lo = center - 0.5 * box_ext;
    Vector3 bbox_hi = center + 0.5 * box_ext;

    shader.setUniform("bbox_lo", &asLvalue(Math::wrapMatrix(bbox_lo)));
    shader.setUniform("bbox_hi", &asLvalue(Math::wrapMatrix(bbox_hi)));
  }

  return true;
}

bool
initFaceIdShader(IShader & shader)
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

  if (!shader.attachModuleFromString(IShader::ModuleType::VERTEX, VERTEX_SHADER.c_str())
   || !shader.attachModuleFromString(IShader::ModuleType::FRAGMENT, FRAGMENT_SHADER.c_str()))
    return false;

  return true;
}

bool
ShapeRendererImpl::renderModel(Model const & model, ColorRgba const & color)
{
  render_system->pushShader();
  render_system->pushShapeFlags();
  render_system->pushColorFlags();

  render_system->setColor(color.data());

  render_system->setPolygonSmooth(false);  // smoothing causes blending halos around primitives with OSMesa
  render_system->setLineSmooth(false);
  render_system->setPointSmooth(!color_by_id);  // this also can cause blending halos, but no other way to get circular points

  bool has_transparency = (color.a() <= 0.9999f || background_color.a() <= 0.9999f || color_by_tex2d);
  if (has_transparency && !color_by_id)
  {
    // Enable alpha-blending
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  }
  else
    glDisable(GL_BLEND);

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
    render_system->beginPrimitive(IRenderSystem::Primitive::POINTS);

      for (size_t i = 0; i < model.points.size(); ++i)
      {
        if (color_by_id)
        {
          ColorRgba c = indexToColor((uint32)i, true);
          render_system->setColor(c.data());
        }
        else if (color_by_label || color_by_features)
          render_system->setColor(model.point_colors[i].data());

        render_system->sendVertex(3, model.points[i].data());
      }

    render_system->endPrimitive();
  }
  else
  {
    // Initialize the shader
    if (color_by_id || (selected_binary_mask && !selected_mesh.empty()))
    {
      if (!face_id_shader)
      {
        face_id_shader = render_system->createShader("Face ID shader");
        if (!face_id_shader)
        {
          THEA_ERROR << "Could not create face ID shader";
          return false;
        }

        if (!initFaceIdShader(*face_id_shader))
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

        if (!initMeshShader(*mesh_shader, material, light_dir, two_sided, out_normals, toon_shading, matcap_tex, tex2d, tex3d,
                            model.mesh_group.getBounds()))
        {
          THEA_ERROR << "Could not initialize mesh shader";
          return false;
        }
      }

      render_system->setShader(mesh_shader);
    }

    RenderOptions opts;
    opts.setSendNormals(true);
    opts.setSendColors(true);
    opts.setSendTexCoords((bool)tex2d);

    if (show_edges)
      opts.setDrawEdges(true).setEdgeColor(edge_color.data());
    else
      opts.setDrawEdges(false);

    if (has_transparency && !color_by_id)
    {
      // First back faces...
      render_system->setCullFace(IRenderSystem::CullFace::FRONT);
      if (!model.mesh_group.draw(render_system, &opts)) return false;

      // ... then front faces
      render_system->setCullFace(IRenderSystem::CullFace::BACK);
      if (!model.mesh_group.draw(render_system, &opts)) return false;
    }
    else
    {
      if (!model.mesh_group.draw(render_system, &opts)) return false;
    }
  }

  render_system->popColorFlags();
  render_system->popShapeFlags();
  render_system->popShader();

  if (char const * err = render_system->getLastError())
  { THEA_ERROR << "Rendering error (" << err << ')'; return false; }

  return true;
}

bool
ShapeRendererImpl::loadPlugins(int argc, char ** argv)
{
  std::string plugin_path = Application::getPluginPath("TheaPluginGL");
  if (plugin_path.empty())
    throw Error("Could not locate OpenGL plugin 'TheaPluginGL'");

  IPlugin * gl_plugin = Application::getPluginManager().load(plugin_path);
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
