# Thea
A toolkit for visual computing with a focus on geometry processing.

Source code is available at [https://github.com/sidch/Thea](https://github.com/sidch/Thea).<br/>
Online documentation is available at [https://sidch.github.io/Thea](https://sidch.github.io/Thea).

Author: [Siddhartha Chaudhuri](https://www.cse.iitb.ac.in/~sidch). Released under the BSD license (see `LICENSE.txt`).

If you find a bug, please let me know promptly. Thank you!

## What is Thea?

*Thea* is a library of C++ classes for computer graphics, primarily for [3D geometry processing](https://www.cse.iitb.ac.in/~cs749/spr2017). It is the core library I use for nearly all my research projects, and it is also the core library for [Adobe Fuse](https://www.adobe.com/products/fuse.html), which I originally authored. As such, it is developed for personal use and its features reflect this: please do not write to me to asking for specific features to be included. However, over time, it has become quite general-purpose. Among its features are:

* Polygon mesh classes with arbitrary per-element attributes, including heavyweight ones that store full mesh topology and a lightweight one designed only for rendering.
* General linear algebra via [Eigen](http://eigen.tuxfamily.org), geometric transformations (e.g. rigid transforms), shared-library-safe pure virtual wrappers for matrices, and easy-to-use interfaces to various solver packages (NNLS, CSPARSE, ARPACK).
* 2, 3 and N-dimensional geometric primitives, including lines, line segments, rays, (hyper)planes, triangles (+ ray-triangle and triangle-triangle intersections), balls, axis-aligned boxes, oriented boxes, polygons and spline curves (+ fast spline-fitting to points).
* An eclectic collection of algorithms, including a fast N-dimensional KD-tree (on points or mesh triangles), shortest paths in graphs, best-fit boxes and ellipsoids, singular value decomposition and PCA, iterative closest point (ICP), symmetry detection, convex hulls, connected components, surface parametrization (global and local), discrete Laplace-Beltrami operators, sampling points from meshes, mesh features (curvature, distance histogram, shape diameter, spin image), and some machine learning models.
* Basic image processing (wrapper for [FreeImage](http://freeimage.sourceforge.net/)).
* A plugin architecture and included plugins providing easy interfaces to OpenGL, [ARPACK](http://www.caam.rice.edu/software/ARPACK/) and [CSPARSE](http://people.sc.fsu.edu/~jburkardt/c_src/csparse/csparse.html). The OpenGL plugin optionally (and easily) compiles with an [OSMesa](https://www.mesa3d.org/osmesa.html) driver to automatically create a headless CPU-only context.
* A variety of utility classes for filesystem navigation, serialization, timing, synchronization, hashing, logging, string manipulation/searching, memory allocation, bounded/sorted arrays, pseudo-random numbers, mathematics (including algebraic roots of polynomials upto degree 4) etc.
* Several bundled tools for 3D file viewing and annotation (*Browse3D*); offline rendering (*RenderShape*); mesh sampling (*MeshSample*), repair (*MeshFix*), features (*MeshLocalFeatures*, *MeshGlobalFeatures*) and format conversion (*MeshConv*); rigid (*ShapeAlign*) and non-rigid (*Register*) shape registration; k-NN graphs of surface samples (*SampleGraph*) etc.

**Thea is constantly under development and many parts are incomplete. Use at your own risk!** I do not provide any support (unless you have bugs to report), and I make no correctness or robustness guarantees for any part of the code. Parts of the library are reasonably battle-tested (e.g. in Fuse), and parts are one-off inclusions rarely used in anger or tested thoroughly.

*Thea* is heavily influenced by, and borrows code from, Morgan McGuire's [G3D](https://casual-effects.com/g3d) library. It started out as an extension of G3D.

The *Thea* library is not related to the independently and contemporaneously developed [Thea Render](http://www.thearender.com) photorealistic rendering engine.

## Installation

*Thea* is written in standards-compliant C++11, and should compile with any recent compiler on Mac, Linux and Windows. It uses [CMake](https:///cmake.org) as a cross-platform buildsystem. However, I do not normally work on Windows, and do not currently provide build instructions for this platform. I have successfully done Windows builds in the past and there is no reason why it should not work with a bit of effort getting the dependencies installed. I will try to add full Windows instructions in the future, time permitting.

### Installing the dependencies

*Thea* relies on [Boost](https://www.boost.org), [Eigen](https://eigen.tuxfamily.org) (3.3 or later), [lib3ds](https://code.google.com/archive/p/lib3ds), [FreeImage](http://freeimage.sourceforge.net) and [ARPACK](http://www.caam.rice.edu/software/ARPACK). (A couple of classes also depend on [CGAL](https://www.cgal.org), but these are optional -- if CGAL is not found on the system, these classes will be omitted from the build.) A convenient script installs all of these on Unix-like systems (Mac and Linux), as follows. Both local (no root) and system-wide (needs root) installs are supported.

Assume `$basedir` is some directory where you're going to check out the source code, and `$prefix` is some directory where you'll install stuff (e.g. `$basedir/Installations` or `/usr/local`).
```shell
cd "$basedir"
git clone --recursive https://bitbucket.org/sidch/theadepsunix TheaDepsUnix
cd TheaDepsUnix/Source
```
**For a local install (no root perms needed to write to `$prefix`):**
```shell
./install-defaults.sh --with-wxwidgets --prefix "$prefix" -j4
```
**For a system-wide install:**
```shell
sudo ./install-defaults.sh --with-wxwidgets --use-root --prefix "$prefix" -j4
```
`--use-root` will try to use `apt-get` on Ubuntu/Debian, omit it if you want to build everything from scratch regardless. `--with-wxwidgets` is needed to build *Browse3D*, a bundled GUI application for viewing 3D files: it can be omitted if so desired. Add `--with-osmesa` to install OSMesa for headless CPU-only rendering (good for remote servers). Replace 4 with the actual number of hardware threads on your system, typically 2, 4, or 8.

The above step will install the necessary libraries by compiling them from source (if not `apt-get`able) and placing the result in `$prefix`. To install an individual library, call `install.sh --with-<package> ...` with the same `--use-root`, `--prefix` etc arguments as above. Carefully check for errors (warnings are generally ok). If there are errors, you probably need to explicitly install some third-party libraries/tools -- see the error messages -- and rerun the command. E.g. I have found some barebones server installs without the [zlib](https://www.zlib.net) or [Expat](https://libexpat.github.io) development packages, required by Mesa: in this particular case, the source packages are included and you can use `./install.sh --with-expat --with-zlib ...`. *Make sure there are no errors in the output before proceeding further.*

### Installing the Thea library, plugins and bundled tools

Assuming there were no errors while installing the dependencies, execute the following commands:
```shell
cd "$basedir"
git clone --recursive https://github.com/sidch/Thea
cd Thea/Code/Build
cmake -DCMAKE_INSTALL_PREFIX="$prefix" -DCMAKE_BUILD_TYPE=RelWithDebInfo .
make -j4
make install    # add sudo if necessary
```
By default, CMake looks for the dependencies in the `CMAKE_INSTALL_PREFIX` directory. If for some reason the dependencies are located somewhere else, e.g. in `$deps`, you can point CMake to it by adding `-DTHEA_DEPS_ROOT="$deps"`. The bundled tools are installed to `$prefix/bin/Thea`: add this to your executable search path (e.g. the system `PATH` variable) as needed. A quick way to check if everything has installed correctly is to run
```shell
$prefix/bin/Thea/Browse3D ../../Data/Models/teapot.obj
```
and see if a window pops up displaying a 3D teapot, or
```shell
$prefix/bin/Thea/RenderShape ../../Data/Models/teapot.obj teapot.png 800 600
```
to render the teapot to an image file.

To build with OSMesa instead of the system OpenGL driver, add `-DWITH_OSMESA=true` to the CMake line above. The *RenderShape* tool will then use OSMesa. To run some test scripts (several probably out of date), run `make test` after building. To omit building the tests altogether, pass `-DWITH_TESTS=false` to CMake. To change the build type (by default `Release`), set `-DCMAKE_BUILD_TYPE=Debug|Release|RelWithDebInfo`.

***
**C++ versions**: Thea itself follows the C++11 standard. However, other packages it depends on may enforce more recent standards: e.g. CGAL has recently moved to C++14 by default. Hence, the build script detects the latest standard (upto C++17) supported by the compiler and uses that to build. To force the compiler to operate in strict C++11 mode, pass `-DFORCE_CXX11=true` to CMake.
***

## Documentation

To generate HTML documentation for the API, run [Doxygen](http://www.doxygen.org) in the `Thea/Code/Documentation` folder. Then, open `html/index.html` in a browser.

This is probably the best place to start looking at the toolkit.

Note that many convenience types, such as `Vector3` and `Matrix4`, are typedefs (for `Eigen::Matrix<Real, 3, 1, ...>` and `Eigen::Matrix<Real, 4, 4, ...>` respectively) and don't show up in the Class Index. Some will show up in Namespaces --> Namespace Members --> Typedefs. For others, you will have to look at the source code. Documenting all of these properly is work-in-progress.)

## Using the library

***
**GCC/Clang-specific**: You **MUST** compile with strict aliasing turned OFF. This is achieved with `-fno-strict-aliasing`. I also recommend `-Wall -g2 -O2` (all "W"arnings, debu"g"gable binaries, "O"ptimize for speed). `-O2` messes up the debugging a bit so turn it off temporarily if you can't track down your bug.
***

The usual command line to link your program with the library is:
```
c++ -std=c++11 -Wall -g2 -O2 -fno-strict-aliasing \
    -I"$prefix/include" -I"$prefix/include/eigen3" \
    <source-files> \
    -L"$prefix/lib" -lThea \
    -lboost_filesystem-mt -lboost_system-mt -lboost_thread-mt \
    -lm \
    [-ldl] [-framework Carbon]
```

If you're using CMake for your own code, a convenient `FindThea.cmake` module in `Thea/Code/Build/Common/CMake/Modules` (or directly from <https://github.com/sidch/CMake>) allows you to do `FIND_PACKAGE(Thea)`, including locating all the necessary dependencies.

## Sample code

Here is a simple "Hello World" example:
```cpp
#include <Thea/MatVec.hpp>

int main(int argc, char * argv[])
{
  using namespace Thea;

  Vector3 v(1.0, 2.0, 3.0);
  Matrix3 m = Math::rotationArc(Vector3(1, 0, 0), Vector3(0, 1, 0));

  // Automatically adds newline and synchronization to std::cout...
  THEA_CONSOLE << "Hello world! The product is " << toString(m * v);
}
```
For real-world samples, see the applications in the `Thea/Code/Source/Tools` folder.
