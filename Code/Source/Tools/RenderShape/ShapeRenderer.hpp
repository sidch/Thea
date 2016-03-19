#ifndef __Thea_RenderShape_ShapeRenderer__
#define __Thea_RenderShape_ShapeRenderer__

#include "../../Common.hpp"

using namespace std;

class ShapeRendererImpl;

class ShapeRenderer
{
  public:
    ShapeRenderer(int argc, char * argv[]);  // just loads plugins and initializes variables
    ~ShapeRenderer();

    int exec(string const & cmdline);
    int exec(int argc, char ** argv);

  private:
    ShapeRendererImpl * impl;

}; // class ShapeRenderer

#endif
