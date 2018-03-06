#include "../../Common.hpp"
#include "../../Array.hpp"
#include "../../FileSystem.hpp"
#include "../../Image.hpp"
#include "../../List.hpp"
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <utility>

#define __NO_STD_VECTOR // Use cl::vector instead of STL version
#ifdef __APPLE__
#  include <OpenCL/opencl.h>
#else
#  include <CL/cl.h>
#endif

// #define TESTING

using namespace Thea;

//=============================================================================================================================
// Basic classes
//=============================================================================================================================

// Operations
struct Op
{
  enum Value
  {
    GREYSCALE,
    INVERT,
    THIN,
    THRESHOLD,
  };

  THEA_ENUM_CLASS_BODY(Op)
};

// Map from operations to arguments
typedef TheaArray<std::string> Args;
typedef std::pair<Op, Args> OpArgs;
typedef TheaList<OpArgs> OpArgsList;

// Declare functions and map from operations to functions
#define GIP_DECL_OP(name) bool name(cl_mem & inbuf, cl_mem & outbuf, Args const & args);

GIP_DECL_OP(greyscale)
GIP_DECL_OP(invert)
GIP_DECL_OP(thin)
GIP_DECL_OP(threshold)

#undef GIP_DECL_OP

typedef bool (*Func)(cl_mem &, cl_mem &, Args const &);
typedef std::pair<Op, Func> OpFunc;

OpFunc OP_TABLE[] = {
  OpFunc(Op::GREYSCALE,       greyscale),
  OpFunc(Op::INVERT,          invert),
  OpFunc(Op::THIN,            thin),
  OpFunc(Op::THRESHOLD,       threshold),
};

// Program options
struct Config
{
  Config()
  : device_type(CL_DEVICE_TYPE_DEFAULT),
    verbose(false),
    alignment(4)
  {}

  std::string in_path;
  std::string out_path;
  cl_device_type device_type;
  bool verbose;
  int alignment;

  OpArgsList ops;

}; // struct Config

Config config;

// Global variables
struct Globals
{
  Globals()
  : image_width(0),
    image_height(0),
    image_aligned_width(0),
    image_aligned_height(0)
  {}

  int image_width;
  int image_height;
  int image_aligned_width;
  int image_aligned_height;

}; // struct Globals

Globals globals;

//=============================================================================================================================
// Parse options
//=============================================================================================================================

bool
usage(int argc, char * argv[])
{
  std::cout << "\nGIP: (G)PU (I)mage (P)rocessor\n"
            << "Siddhartha Chaudhuri, 2013\n\n"
            << "Usage: " << argv[0] << " [options] <infile> <outfile>\n\n"
            << "Options:\n"
            << "  --cpu                     Force execution on CPU\n"
            << "  --gpu                     Force execution on GPU\n"
            << "  --acc                     Force execution on accelerator\n"
            << "  --verbose                 Print debugging messages\n"
            << '\n'
            << "  --greyscale               Convert image to greyscale\n"
            << "  --invert                  Invert image colors\n"
            << "  --thin <iters>            Perform repeated thinning (e.g. to skeletonize)\n"
            << "  --threshold <value>       Convert to binary image by thresholding luminance\n"
            << std::endl;
  return false;
}

//=============================================================================================================================
// Parse options
//=============================================================================================================================

int
parseOpts(int argc, char * argv[])
{
  int num_positionals = 0;
  for (int i = 1; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (beginsWith(arg, "--"))
    {
      if (arg == "--greyscale")
        config.ops.push_back(OpArgs(Op::GREYSCALE, Args()));
      else if (arg == "--invert")
        config.ops.push_back(OpArgs(Op::INVERT, Args()));
      else if (arg == "--threshold")
      {
        if (i >= argc - 1)
          return usage(argc, argv);

        Args op_args; op_args.push_back(argv[++i]);
        config.ops.push_back(OpArgs(Op::THRESHOLD, op_args));
      }
      else if (arg == "--thin")
      {
        if (i >= argc - 1)
          return usage(argc, argv);

        Args op_args; op_args.push_back(argv[++i]);
        config.ops.push_back(OpArgs(Op::THIN, op_args));
      }
      else if (arg == "--cpu")
        config.device_type = CL_DEVICE_TYPE_CPU;
      else if (arg == "--gpu")
        config.device_type = CL_DEVICE_TYPE_GPU;
      else if (arg == "--acc")
        config.device_type = CL_DEVICE_TYPE_ACCELERATOR;
      else if (arg == "--verbose")
        config.verbose = true;
      else if (arg == "--help")
      {
        usage(argc, argv);
        std::exit(0);
      }
      else
        return usage(argc, argv);
    }
    else
    {
      if (num_positionals == 0)
        config.in_path = arg;
      else if (num_positionals == 1)
        config.out_path = arg;
      else
        return usage(argc, argv);

      num_positionals++;
    }
  }

  if (num_positionals < 2)
    return usage(argc, argv);

  return true;
}

//=============================================================================================================================
// Basic image functions
//=============================================================================================================================

int
alignSize(int s, int alignment = 4)
{
  int padding = (alignment - s % alignment) % alignment;
  return s + padding;
}

template <typename T, int NUM_CHANNELS, uint32 MAX_VALUE> struct PixelToFloats {};
template <typename T, uint32 MAX_VALUE> struct PixelToFloats<T, 1, MAX_VALUE>
{
  inline static void conv(T const * pixel, float32 * out)
  {
    out[0] = out[1] = out[2] = *pixel / (float32)MAX_VALUE;
    out[3] = 1.0f;
  }
};

template <typename T, uint32 MAX_VALUE> struct PixelToFloats<T, 3, MAX_VALUE>
{
  inline static void conv(T const * pixel, float32 * out)
  {
    out[0] = pixel[Image::Channel::RED  ] / (float32)MAX_VALUE;
    out[1] = pixel[Image::Channel::GREEN] / (float32)MAX_VALUE;
    out[2] = pixel[Image::Channel::BLUE ] / (float32)MAX_VALUE;
    out[3] = 1.0f;
  }
};

template <typename T, uint32 MAX_VALUE> struct PixelToFloats<T, 4, MAX_VALUE>
{
  inline static void conv(T const * pixel, float32 * out)
  {
    out[0] = pixel[Image::Channel::RED  ] / (float32)MAX_VALUE;
    out[1] = pixel[Image::Channel::GREEN] / (float32)MAX_VALUE;
    out[2] = pixel[Image::Channel::BLUE ] / (float32)MAX_VALUE;
    out[3] = pixel[Image::Channel::ALPHA] / (float32)MAX_VALUE;
  }
};

template <typename T, int NUM_CHANNELS, uint32 MAX_VALUE>
struct FlattenToFloatArray
{
  static bool conv(Image const & image, float32 * buffer)
  {
    int width = image.getWidth();
    int height = image.getHeight();
    int aligned_width = alignSize(width, config.alignment);

    for (int i = 0; i < height; ++i)
    {
      T const * in_pixel = (T const *)image.getScanLine(i);
      float32 * out_pixel = &buffer[i * aligned_width * 4];
      for (int j = 0; j < width; ++j, in_pixel += NUM_CHANNELS, out_pixel += 4)
        PixelToFloats<T, NUM_CHANNELS, MAX_VALUE>::conv(in_pixel, out_pixel);
    }

    return true;
  }
};

template <>
struct FlattenToFloatArray<float32, 4, 1>
{
  static bool conv(Image const & image, float32 * buffer)
  {
    int width = image.getWidth();
    int height = image.getHeight();
    int aligned_width = alignSize(width, config.alignment);

    for (int i = 0; i < height; ++i)
    {
      float32 const * in_scanline = (float32 const *)image.getScanLine(i);
      float32 * out_scanline = &buffer[i * aligned_width * 4];
      std::memcpy(out_scanline, in_scanline, width * 4 * sizeof(float32));
    }

    return true;
  }
};

bool
flattenToFloatArray(Image const & image, float32 * buffer)
{
  switch (image.getType())
  {
    case Image::Type::LUMINANCE_8U:    return FlattenToFloatArray<uint8,    1,  0xFF      >::conv(image, buffer);
    case Image::Type::LUMINANCE_16:    return FlattenToFloatArray<int16,    1,  0x7FFF    >::conv(image, buffer);
    case Image::Type::LUMINANCE_16U:   return FlattenToFloatArray<uint16,   1,  0xFFFF    >::conv(image, buffer);
    case Image::Type::LUMINANCE_32:    return FlattenToFloatArray<int32,    1,  0x7FFFFFFF>::conv(image, buffer);
    case Image::Type::LUMINANCE_32U:   return FlattenToFloatArray<uint32,   1,  0xFFFFFFFF>::conv(image, buffer);
    case Image::Type::RGB_8U:          return FlattenToFloatArray<uint8,    3,  0xFF      >::conv(image, buffer);
    case Image::Type::RGB_16U:         return FlattenToFloatArray<uint16,   3,  0xFFFF    >::conv(image, buffer);
    case Image::Type::RGB_32F:         return FlattenToFloatArray<float32,  3,  1         >::conv(image, buffer);
    case Image::Type::RGBA_8U:         return FlattenToFloatArray<uint8,    4,  0xFF      >::conv(image, buffer);
    case Image::Type::RGBA_16U:        return FlattenToFloatArray<uint16,   4,  0xFFFF    >::conv(image, buffer);
    case Image::Type::RGBA_32F:        return FlattenToFloatArray<float32,  4,  1         >::conv(image, buffer);
    default: return false;
  }

  return true;
}

template <typename T, uint32 MAX_VALUE> T convertValue(float32 val) { return static_cast<T>(MAX_VALUE * val + 0.5); }
template <> float32 convertValue<float32, 1>(float32 val) { return val; }
template <> float64 convertValue<float64, 1>(float32 val) { return val; }

template <typename T, int NUM_CHANNELS, uint32 MAX_VALUE> struct PixelFromFloats {};
template <typename T, uint32 MAX_VALUE> struct PixelFromFloats<T, 1, MAX_VALUE>
{
  inline static void conv(float32 const * in, T * pixel)
  {
    *pixel = convertValue<T, MAX_VALUE>(0.2126 * in[0] + 0.7152 * in[1] + 0.0722  * in[2]);
  }
};

template <typename T, uint32 MAX_VALUE> struct PixelFromFloats<T, 3, MAX_VALUE>
{
  inline static void conv(float32 const * in, T * pixel)
  {
    pixel[Image::Channel::RED  ] = convertValue<T, MAX_VALUE>(in[0]);
    pixel[Image::Channel::GREEN] = convertValue<T, MAX_VALUE>(in[1]);
    pixel[Image::Channel::BLUE ] = convertValue<T, MAX_VALUE>(in[2]);
  }
};

template <typename T, uint32 MAX_VALUE> struct PixelFromFloats<T, 4, MAX_VALUE>
{
  inline static void conv(float32 const * in, T * pixel)
  {
    pixel[Image::Channel::RED  ] = convertValue<T, MAX_VALUE>(in[0]);
    pixel[Image::Channel::GREEN] = convertValue<T, MAX_VALUE>(in[1]);
    pixel[Image::Channel::BLUE ] = convertValue<T, MAX_VALUE>(in[2]);
    pixel[Image::Channel::ALPHA] = convertValue<T, MAX_VALUE>(in[3]);
  }
};

template <typename T, int NUM_CHANNELS, uint32 MAX_VALUE>
struct UnflattenFromFloatArray
{
  static bool conv(float32 const * buffer, Image & image)
  {
    int width = image.getWidth();
    int height = image.getHeight();
    int aligned_width = alignSize(width, config.alignment);

    for (int i = 0; i < height; ++i)
    {
      float32 const * in_pixel = &buffer[i * aligned_width * 4];
      T * out_pixel = (T *)image.getScanLine(i);
      for (int j = 0; j < width; ++j, in_pixel += 4, out_pixel += NUM_CHANNELS)
        PixelFromFloats<T, NUM_CHANNELS, MAX_VALUE>::conv(in_pixel, out_pixel);
    }

    return true;
  }
};

template <>
struct UnflattenFromFloatArray<float32, 4, 1>
{
  static bool conv(float32 const * buffer, Image & image)
  {
    int width = image.getWidth();
    int height = image.getHeight();
    int aligned_width = alignSize(width, config.alignment);

    for (int i = 0; i < height; ++i)
    {
      float32 const * in_scanline = &buffer[i * aligned_width * 4];
      float32 * out_scanline = (float32 *)image.getScanLine(i);
      std::memcpy(out_scanline, in_scanline, width * 4 * sizeof(float32));
    }

    return true;
  }
};

bool
unflattenFromFloatArray(float32 const * buffer, Image & image)
{
  switch (image.getType())
  {
    case Image::Type::LUMINANCE_8U:    return UnflattenFromFloatArray<uint8,    1,  0xFF      >::conv(buffer, image);
    case Image::Type::LUMINANCE_16:    return UnflattenFromFloatArray<int16,    1,  0x7FFF    >::conv(buffer, image);
    case Image::Type::LUMINANCE_16U:   return UnflattenFromFloatArray<uint16,   1,  0xFFFF    >::conv(buffer, image);
    case Image::Type::LUMINANCE_32:    return UnflattenFromFloatArray<int32,    1,  0x7FFFFFFF>::conv(buffer, image);
    case Image::Type::LUMINANCE_32U:   return UnflattenFromFloatArray<uint32,   1,  0xFFFFFFFF>::conv(buffer, image);
    case Image::Type::RGB_8U:          return UnflattenFromFloatArray<uint8,    3,  0xFF      >::conv(buffer, image);
    case Image::Type::RGB_16U:         return UnflattenFromFloatArray<uint16,   3,  0xFFFF    >::conv(buffer, image);
    case Image::Type::RGB_32F:         return UnflattenFromFloatArray<float32,  3,  1         >::conv(buffer, image);
    case Image::Type::RGBA_8U:         return UnflattenFromFloatArray<uint8,    4,  0xFF      >::conv(buffer, image);
    case Image::Type::RGBA_16U:        return UnflattenFromFloatArray<uint16,   4,  0xFFFF    >::conv(buffer, image);
    case Image::Type::RGBA_32F:        return UnflattenFromFloatArray<float32,  4,  1         >::conv(buffer, image);
    default: return false;
  }

  return true;
}

//=============================================================================================================================
// OpenCL stuff
//=============================================================================================================================

namespace CL {

cl_platform_id platform_id = NULL;
cl_device_id device_id = NULL;
cl_context context;
cl_command_queue command_queue;

}; // namespace CL

std::string
clErrorString(cl_int err)
{
  switch (err)
  {
    case CL_SUCCESS:                            return "Success!";
    case CL_DEVICE_NOT_FOUND:                   return "Device not found.";
    case CL_DEVICE_NOT_AVAILABLE:               return "Device not available";
    case CL_COMPILER_NOT_AVAILABLE:             return "Compiler not available";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE:      return "Memory object allocation failure";
    case CL_OUT_OF_RESOURCES:                   return "Out of resources";
    case CL_OUT_OF_HOST_MEMORY:                 return "Out of host memory";
    case CL_PROFILING_INFO_NOT_AVAILABLE:       return "Profiling information not available";
    case CL_MEM_COPY_OVERLAP:                   return "Memory copy overlap";
    case CL_IMAGE_FORMAT_MISMATCH:              return "Image format mismatch";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED:         return "Image format not supported";
    case CL_BUILD_PROGRAM_FAILURE:              return "Program build failure";
    case CL_MAP_FAILURE:                        return "Map failure";
    case CL_INVALID_VALUE:                      return "Invalid value";
    case CL_INVALID_DEVICE_TYPE:                return "Invalid device type";
    case CL_INVALID_PLATFORM:                   return "Invalid platform";
    case CL_INVALID_DEVICE:                     return "Invalid device";
    case CL_INVALID_CONTEXT:                    return "Invalid context";
    case CL_INVALID_QUEUE_PROPERTIES:           return "Invalid queue properties";
    case CL_INVALID_COMMAND_QUEUE:              return "Invalid command queue";
    case CL_INVALID_HOST_PTR:                   return "Invalid host pointer";
    case CL_INVALID_MEM_OBJECT:                 return "Invalid memory object";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:    return "Invalid image format descriptor";
    case CL_INVALID_IMAGE_SIZE:                 return "Invalid image size";
    case CL_INVALID_SAMPLER:                    return "Invalid sampler";
    case CL_INVALID_BINARY:                     return "Invalid binary";
    case CL_INVALID_BUILD_OPTIONS:              return "Invalid build options";
    case CL_INVALID_PROGRAM:                    return "Invalid program";
    case CL_INVALID_PROGRAM_EXECUTABLE:         return "Invalid program executable";
    case CL_INVALID_KERNEL_NAME:                return "Invalid kernel name";
    case CL_INVALID_KERNEL_DEFINITION:          return "Invalid kernel definition";
    case CL_INVALID_KERNEL:                     return "Invalid kernel";
    case CL_INVALID_ARG_INDEX:                  return "Invalid argument index";
    case CL_INVALID_ARG_VALUE:                  return "Invalid argument value";
    case CL_INVALID_ARG_SIZE:                   return "Invalid argument size";
    case CL_INVALID_KERNEL_ARGS:                return "Invalid kernel arguments";
    case CL_INVALID_WORK_DIMENSION:             return "Invalid work dimension";
    case CL_INVALID_WORK_GROUP_SIZE:            return "Invalid work group size";
    case CL_INVALID_WORK_ITEM_SIZE:             return "Invalid work item size";
    case CL_INVALID_GLOBAL_OFFSET:              return "Invalid global offset";
    case CL_INVALID_EVENT_WAIT_LIST:            return "Invalid event wait list";
    case CL_INVALID_EVENT:                      return "Invalid event";
    case CL_INVALID_OPERATION:                  return "Invalid operation";
    case CL_INVALID_GL_OBJECT:                  return "Invalid OpenGL object";
    case CL_INVALID_BUFFER_SIZE:                return "Invalid buffer size";
    case CL_INVALID_MIP_LEVEL:                  return "Invalid mip-map level";
    default:                                    return format("errorcode %d", err);
  }
}

void
checkCL(cl_int err, char const * loc = NULL)
{
  if (err != CL_SUCCESS)
  {
    if (loc)
      THEA_ERROR << "OpenCL error in " << loc << ": " << clErrorString(err);
    else
      THEA_ERROR << "OpenCL error: " << clErrorString(err);

    std::exit(-1);
  }
}

bool
initCL()
{
  // From: http://www.thebigblob.com/getting-started-with-opencl-and-gpu-computing/

  // Get platform and device information
  cl_int ret = clGetPlatformIDs(1, &CL::platform_id, NULL);
  checkCL(ret, "initCL::getPlatformIDs");

  ret = clGetDeviceIDs(CL::platform_id, config.device_type, 1, &CL::device_id, NULL);
  checkCL(ret, "initCL::getDeviceIDs");

  // Create an OpenCL context
  CL::context = clCreateContext(NULL, 1, &CL::device_id, NULL, NULL, &ret);
  checkCL(ret, "initCL::createContext");

  // Create a command queue
  CL::command_queue = clCreateCommandQueue(CL::context, CL::device_id, 0, &ret);
  checkCL(ret, "initCL::createCommandQueue");

  return true;
}

void
cleanupCL()
{
  clReleaseCommandQueue(CL::command_queue);
  clReleaseContext(CL::context);
}

cl_mem
copyImageToCL(Image & image)
{
  int aligned_width = alignSize(image.getWidth(), config.alignment);
  int aligned_height = alignSize(image.getHeight(), config.alignment);
  size_t buf_size = (size_t)(aligned_width * aligned_height * 4 * sizeof(float32));

  TheaArray<float32> host_buffer((size_t)buf_size, 0);
  if (!flattenToFloatArray(image, &host_buffer[0]))
  {
    THEA_ERROR << "Could not flatten image to a 1D floating-point buffer";
    return 0;
  }

  cl_int ret;
  cl_mem buffer = clCreateBuffer(CL::context, CL_MEM_COPY_HOST_PTR, buf_size, &host_buffer[0], &ret);
  checkCL(ret, "copyImageToCL::createBuffer");

  return buffer;
}

bool
copyImageFromCL(cl_mem buf, Image & image)
{
  int aligned_width = alignSize(image.getWidth(), config.alignment);
  int aligned_height = alignSize(image.getHeight(), config.alignment);
  size_t buf_size = (size_t)(aligned_width * aligned_height * 4 * sizeof(float32));

  TheaArray<float32> host_buffer((size_t)buf_size);
  cl_int ret = clEnqueueReadBuffer(CL::command_queue, buf, CL_TRUE, 0, buf_size, &host_buffer[0], 0, NULL, NULL);
  checkCL(ret, "copyImageFromCL::enqueueReadBuffer");

  // No need to flush, we issued a blocking call above

  if (!unflattenFromFloatArray(&host_buffer[0], image))
  {
    THEA_ERROR << "Could not unflatten image from 1D floating-point buffer";
    return false;
  }

  return true;
}

//=============================================================================================================================
// Main function
//=============================================================================================================================

int
main(int argc, char * argv[])
{
  if (!parseOpts(argc, argv))
    return -1;

  // Is there anything to do?
  if (config.ops.empty())
  {
    if (!FileSystem::fileExists(config.in_path))
    {
      THEA_ERROR << "Could not load input image " << config.in_path;
      return -1;
    }

    std::string resolved_in_path = FileSystem::resolve(config.in_path);
    std::string resolved_out_path = FileSystem::resolve(config.out_path);
    if (resolved_in_path == resolved_out_path)
      return 0;
  }

  // Load image
  Image image;
  try
  {
    image.load(config.in_path);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not load input image %s", config.in_path.c_str())

  if (!config.ops.empty())
  {
#ifdef TESTING
    // Setting this environment variable forces the OpenCL source code to be recompiled every time on NVIDIA cards
    setenv("CUDA_CACHE_DISABLE", "1", 1);
#endif

    // Initialize
    if (!initCL())
      return -1;

    // Transfer image to an OpenCL buffer
    cl_mem inbuf = copyImageToCL(image);
    if (!inbuf)
      return -1;

    // Initialize global variables
    globals.image_width = image.getWidth();
    globals.image_height = image.getHeight();
    globals.image_aligned_width = alignSize(globals.image_width, config.alignment);
    globals.image_aligned_height = alignSize(globals.image_height, config.alignment);

    // Create a second OpenCL buffer for the output
    size_t buf_size = (size_t)(globals.image_aligned_width * globals.image_aligned_height * 4 * sizeof(float32));
    cl_int ret;
    cl_mem outbuf = clCreateBuffer(CL::context, CL_MEM_READ_WRITE, buf_size, NULL, &ret);
    checkCL(ret, "main::createBuffer");

    // Apply operations
    for (OpArgsList::const_iterator oi = config.ops.begin(); oi != config.ops.end(); ++oi)
    {
      size_t NUM_OPS = sizeof(OP_TABLE) / sizeof(OpFunc);
      bool found = false;
      for (size_t j = 0; j < NUM_OPS; ++j)
      {
        if (oi->first == OP_TABLE[j].first)
        {
          if (!OP_TABLE[j].second(inbuf, outbuf, oi->second))
            return -1;

          found = true;
          break;
        }
      }

      if (!found)
      {
        THEA_ERROR << "Unsupported operation";
        return -1;
      }

      std::swap(inbuf, outbuf);
    }

    // Copy final image from OpenCL buffer
    // FIXME: Create a new image, so we can load in greyscale and save in color if a colorize filter was applied, for instance
    if (!copyImageFromCL(inbuf, image))
      return -1;

    // cleanup
    ret = clReleaseMemObject(inbuf);
    ret = ret | clReleaseMemObject(outbuf);
    checkCL(ret, "main::releaseMemObject");

    cleanupCL();
  }

  // Save output image
  try
  {
    image.save(config.out_path);
  }
  THEA_STANDARD_CATCH_BLOCKS(return -1;, ERROR, "Could not save output image %s", config.out_path.c_str())
}

//=============================================================================================================================
// Create and launch kernels
//=============================================================================================================================

void
createKernel(char const * prog, char const * caller_name, cl_program & program, cl_kernel & kernel)
{
  // Get the size of the program source
  size_t prog_size = std::strlen(prog);

  // Create a program from the source
  cl_int ret;
  program = clCreateProgramWithSource(CL::context, 1, (char const **)&prog, (size_t const *)&prog_size, &ret);
  checkCL(ret, (std::string(caller_name) + "::createProgramWithSource").c_str());

  // Build the program
  ret = clBuildProgram(program, 1, &CL::device_id, NULL, NULL, NULL);
  if (ret != CL_SUCCESS)
  {
    if (ret == CL_BUILD_PROGRAM_FAILURE)
    {
      size_t log_size;
      clGetProgramBuildInfo(program, CL::device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size);
      char * log = new char[log_size];
      clGetProgramBuildInfo(program, CL::device_id, CL_PROGRAM_BUILD_LOG, log_size, log, NULL);
      THEA_ERROR << "OpenCL error in " << caller_name << "::buildProgram: Program build failure";
      std::cerr << "Build log:\n";
      std::cerr << log << std::endl;
      delete [] log;
      std::exit(-1);
    }
    else
      checkCL(ret, (std::string(caller_name) + "::buildProgram").c_str());
  }

  // Create the OpenCL kernel
  kernel = clCreateKernel(program, caller_name, &ret);
  checkCL(ret, (std::string(caller_name) + "::createKernel").c_str());
}

void
setDefaultKernelArgs(cl_kernel kernel, cl_mem inbuf, cl_mem outbuf, char const * caller_name)
{
  // Set the arguments of the kernel
  cl_int ret;
  ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), &inbuf);
  ret = ret | clSetKernelArg(kernel, 1, sizeof(cl_mem), &outbuf);
  ret = ret | clSetKernelArg(kernel, 2, sizeof(int),    &globals.image_aligned_width);
  checkCL(ret, (std::string(caller_name) + "::setKernelArg").c_str());
}

void
runKernel(cl_kernel kernel, char const * caller_name)
{
  // Execute the OpenCL kernel on the list
  size_t global_work_size[2] = { (size_t)globals.image_width, (size_t)globals.image_height };
  cl_int ret = clEnqueueNDRangeKernel(CL::command_queue, kernel, 2, NULL, global_work_size, NULL, 0, NULL, NULL);
  checkCL(ret, (std::string(caller_name) + "::enqueueNDRangeKernel").c_str());
}

//=============================================================================================================================
// Greyscale
//=============================================================================================================================

bool
greyscale(cl_mem & inbuf, cl_mem & outbuf, Args const & args)
{
  static char const * PROG =
    "__kernel void greyscale(__global float4 * in, __global float4 * out, const int width)"
    "{"
    "  int x = get_global_id(0);"
    "  int y = get_global_id(1);"
    ""
    "  int base_index = (y * width + x);"
    "  float lum = 0.2126f * in[base_index].x + 0.7152f * in[base_index].y + 0.0722f * in[base_index].z;"
    "  lum = clamp(lum, 0.0f, 1.0f);"
    "  out[base_index] = (float4)(lum, lum, lum, in[base_index].w);"
    "}";

  cl_program program;
  cl_kernel kernel;

  createKernel(PROG, "greyscale", program, kernel);
  setDefaultKernelArgs(kernel, inbuf, outbuf, "greyscale");
  runKernel(kernel, "greyscale");

  return true;
}

//=============================================================================================================================
// Invert
//=============================================================================================================================

bool
invert(cl_mem & inbuf, cl_mem & outbuf, Args const & args)
{
  static char const * PROG =
    "__kernel void invert(__global float4 * in, __global float4 * out, const int width)"
    "{"
    "  int x = get_global_id(0);"
    "  int y = get_global_id(1);"
    ""
    "  int base_index = (y * width + x);"
    "  out[base_index] = (float4)(1, 1, 1, 1) - in[base_index];"
    "}";

  cl_program program;
  cl_kernel kernel;

  createKernel(PROG, "invert", program, kernel);
  setDefaultKernelArgs(kernel, inbuf, outbuf, "invert");
  runKernel(kernel, "invert");

  return true;
}

//=============================================================================================================================
// Thinning <iters>
//=============================================================================================================================

bool
thin(cl_mem & inbuf, cl_mem & outbuf, Args const & args)
{
  static char const * PROG =
    "typedef union {"
    "  int8    v;"
    "  int     i[8];"
    "} Int8;"
    ""
    "typedef union {"
    "  float8    v;"
    "  float     f[8];"
    "} Float8;"
    ""
    "__kernel void thin(__global float4 * in, __global float4 * out, const int width, const int height,"
    "                   const float8 structure, const float8 mask)"
    "{"
    "  int x = get_global_id(0);"
    "  int y = get_global_id(1);"
    ""
    "  Int8 rows; rows.v = (int8)(1, 1, 1, 0, 0, -1, -1, -1);"
    "  Int8 cols; cols.v = (int8)(-1, 0, 1, -1, 1, -1, 0, 1);"
    ""
    "  Float8 nbd;"
    "  for (int i = 0; i < 8; ++i)"
    "  {"
    "    int xi = x + cols.i[i];"
    "    int yi = y + rows.i[i];"
    ""
    "    nbd.f[i] = (xi < 0 || xi >= width || yi < 0 || yi >= height) ? 0.0f : in[yi * width + xi].x;"
    "  }"
    ""
    "  int base_index = (y * width + x);"
    ""
    "  float8 diff_abs = fabs(nbd.v - structure);"
    "  float err = dot(mask.lo, diff_abs.lo) + dot(mask.hi, diff_abs.hi);"
    "  float lum = (err < 0.000001f) ? 0.0f : in[base_index].x;"
    "  out[base_index] = (float4)(lum, lum, lum, in[base_index].w);"
    "}";

  int niters;
  if (args.empty() || !(std::istringstream(args[0]) >> niters))
  {
    THEA_ERROR << "Could not parse number of iterations";
    return false;
  }

  if (niters <= 0)
  {
    THEA_ERROR << "Invalid number of iterations (" << niters << ')';
    return false;
  }

  // First threshold the input image to ensure we have a binary image
  Args thresh_args; thresh_args.push_back("0.5");
  if (!threshold(inbuf, outbuf, thresh_args))
    return false;

  std::swap(inbuf, outbuf);

  // Now set up the structures and masks. Each structure and mask is defined as an 8-vector of floats.
  // The index pattern is:
  //  0  1  2
  //  3     4
  //  5  6  7
  //
  // This is the Skeleton:2 kernel in ImageMagick: http://www.imagemagick.org/Usage/morphology/#skeleton2

  static float32 const STRUCTURES[][8] = {
    { 0, 0, 1,
      0,    1,
      0, 0, 1 },

    { 0, 0, 0,
      0,    1,
      0, 1, 0 },

    { 0, 0, 0,
      0,    0,
      1, 1, 1 },

    { 0, 0, 0,
      1,    0,
      0, 1, 0 },

    { 1, 0, 0,
      1,    0,
      1, 0, 0 },

    { 0, 1, 0,
      1,    0,
      0, 0, 0 },

    { 1, 1, 1,
      0,    0,
      0, 0, 0 },

    { 0, 1, 0,
      0,    1,
      0, 0, 0 },
  };

  static float32 const MASKS[][8] = {
    { 1, 0, 1,
      1,    1,
      1, 0, 1 },

    { 1, 1, 0,
      1,    1,
      0, 1, 0 },

    { 1, 1, 1,
      0,    0,
      1, 1, 1 },

    { 0, 1, 1,
      1,    1,
      0, 1, 0 },

    { 1, 0, 1,
      1,    1,
      1, 0, 1 },

    { 0, 1, 0,
      1,    1,
      0, 1, 1 },

    { 1, 1, 1,
      0,    0,
      1, 1, 1 },

    { 0, 1, 0,
      1,    1,
      1, 1, 0 },
  };

  cl_program program;
  cl_kernel kernel;

  createKernel(PROG, "thin", program, kernel);

  for (int iter = 0; iter < niters; ++iter)
  {
    for (size_t j = 0; j < sizeof(STRUCTURES) / (8 * sizeof(float32)); ++j)
    {
      cl_int ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), &inbuf);
      ret =  ret | clSetKernelArg(kernel, 1, sizeof(cl_mem), &outbuf);
      ret =  ret | clSetKernelArg(kernel, 2, sizeof(int), &globals.image_aligned_width);
      ret =  ret | clSetKernelArg(kernel, 3, sizeof(int), &globals.image_aligned_height);
      ret =  ret | clSetKernelArg(kernel, 4, 8 * sizeof(float32), STRUCTURES[j]);
      ret =  ret | clSetKernelArg(kernel, 5, 8 * sizeof(float32), MASKS[j]);
      checkCL(ret, "thin::setKernelArg");

      runKernel(kernel, "thin");

      std::swap(inbuf, outbuf);
    }
  }

  // The final result should be in outbuf
  std::swap(inbuf, outbuf);

  return true;
}

//=============================================================================================================================
// Threshold <value>
//=============================================================================================================================

bool
threshold(cl_mem & inbuf, cl_mem & outbuf, Args const & args)
{
  static char const * PROG =
    "__kernel void threshold(__global float4 * in, __global float4 * out, const int width, const float thresh)"
    "{"
    "  int x = get_global_id(0);"
    "  int y = get_global_id(1);"
    ""
    "  int base_index = (y * width + x);"
    "  float lum = 0.2126f * in[base_index].x + 0.7152f * in[base_index].y + 0.0722f * in[base_index].z;"
    "  float bit = (lum < thresh) ? 0.0f : 1.0f;"
    "  out[base_index] = (float4)(bit, bit, bit, in[base_index].w);"
    "}";

  float32 thresh;
  if (args.empty() || !(std::istringstream(args[0]) >> thresh))
  {
    THEA_ERROR << "Could not parse thresholding value";
    return false;
  }

  cl_program program;
  cl_kernel kernel;

  createKernel(PROG, "threshold", program, kernel);

  setDefaultKernelArgs(kernel, inbuf, outbuf, "threshold");
  cl_int ret = clSetKernelArg(kernel, 3, sizeof(float32), &thresh);
  checkCL(ret, "threshold::setKernelArg");

  runKernel(kernel, "threshold");

  return true;
}
