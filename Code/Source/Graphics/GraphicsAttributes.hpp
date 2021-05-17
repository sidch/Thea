//============================================================================
//
// This file is part of the Thea toolkit.
//
// This software is distributed under the BSD license, as detailed in the
// accompanying LICENSE.txt file. Portions are derived from other works:
// their respective licenses and copyright information are reproduced in
// LICENSE.txt and/or in the relevant source files.
//
// Author: Siddhartha Chaudhuri
// First version: 2009
//
//============================================================================

#ifndef __Thea_Graphics_GraphicsAttributes_hpp__
#define __Thea_Graphics_GraphicsAttributes_hpp__

#include "../Common.hpp"
#include "../Colors.hpp"
#include "RenderOptions.hpp"
#include "IRenderSystem.hpp"

namespace Thea {
namespace Graphics {

/** A dummy attribute that has no data members and does nothing. */
struct /* THEA_API */ NullAttribute  // for some reason VS 2008 objects to dllimport-ing this, even though it has no
                                     // header-defined static functions
{}; // class NullAttribute

/** A utility attribute class that wraps a position in space. */
template <typename PositionT>
struct PositionAttribute
{
  public:
    THEA_DECL_SMART_POINTERS(PositionAttribute)

    typedef PositionT Position;  ///< The point type.

    /** Default constructor. */
    PositionAttribute() {}

    /** Initializing constructor. */
    PositionAttribute(Position const & p_) : p(p_) {}

    /** Get the wrapped position. */
    Position const & getPosition() const { return p; }

    /** Set the wrapped position. */
    void setPosition(Position const & p_) { p = p_; }

  private:
    Position p;

}; // class PositionAttribute

/** A utility attribute class that wraps a normal. */
template <typename NormalT>
struct NormalAttribute
{
  public:
    THEA_DECL_SMART_POINTERS(NormalAttribute)

    typedef NormalT Normal;  ///< The normal type.

    /** Default constructor. */
    NormalAttribute() {}

    /** Initializing constructor. */
    NormalAttribute(Normal const & n_) : n(n_) {}

    /** Get the wrapped normal. */
    Normal const & getNormal() const { return n; }

    /** Set the wrapped normal. */
    void setNormal(Normal const & n_) { n = n_; }

  private:
    Normal n;

}; // class NormalAttribute

/** A utility attribute class that wraps a color. */
template <typename ColorT>
struct ColorAttribute
{
  public:
    THEA_DECL_SMART_POINTERS(ColorAttribute)

    typedef ColorT Color;  ///< The color type.

    /** Default constructor. */
    ColorAttribute() {}

    /** Initializing constructor. */
    ColorAttribute(Color const & c_) : c(c_) {}

    /** Get the wrapped color. */
    Color const & getColor() const { return c; }

    /** Set the wrapped color. */
    void setColor(Color const & c_) { c = c_; }

  private:
    Color c;

}; // class ColorAttribute

/** A utility attribute class that wraps a set of texture coordinates. */
template <typename TexCoordT>
struct TexCoordAttribute
{
  public:
    THEA_DECL_SMART_POINTERS(TexCoordAttribute)

    typedef TexCoordT TexCoord;  ///< The type of the texture coordinates vector.

    /** Default constructor. */
    TexCoordAttribute() {}

    /** Initializing constructor. */
    TexCoordAttribute(TexCoord const & t_) : t(t_) {}

    /** Get the wrapped normal. */
    TexCoord const & getTexCoord() const { return t; }

    /** Set the wrapped normal. */
    void setTexCoord(TexCoord const & t_) { t = t_; }

  private:
    TexCoord t;

}; // class TexCoordAttribute

/**
 * A utility attribute class that wraps both a normal and a color.
 *
 * @see NormalAttribute, ColorAttribute
 */
template <typename NormalT, typename ColorT>
struct NormalColorAttribute : public NormalAttribute<NormalT>, public ColorAttribute<ColorT>
{
  private:
    typedef NormalAttribute<NormalT>  NormalBaseT;
    typedef ColorAttribute<ColorT>    ColorBaseT;

  public:
    THEA_DECL_SMART_POINTERS(NormalColorAttribute)

    typedef typename NormalBaseT::Normal  Normal; ///< The normal type.
    typedef typename ColorBaseT::Color    Color;  ///< The color type.

    /** Default constructor. */
    NormalColorAttribute() {}

    /** Initializing constructor. */
    NormalColorAttribute(Normal const & n_, Color const & c_) : NormalBaseT(n_), ColorBaseT(c_) {}

}; // class NormalColorAttribute

/**
 * A utility attribute class that wraps both a normal and a set of texture coordinates.
 *
 * @see NormalAttribute, TexCoordAttribute
 */
template <typename NormalT, typename TexCoordT>
struct NormalTexCoordAttribute : public NormalAttribute<NormalT>, public TexCoordAttribute<TexCoordT>
{
  private:
    typedef NormalAttribute<NormalT>      NormalBaseT;
    typedef TexCoordAttribute<TexCoordT>  TexCoordBaseT;

  public:
    THEA_DECL_SMART_POINTERS(NormalTexCoordAttribute)

    typedef typename NormalBaseT  ::Normal    Normal;    ///< The normal type.
    typedef typename TexCoordBaseT::TexCoord  TexCoord;  ///< The type of the texture coordinates vector.

    /** Default constructor. */
    NormalTexCoordAttribute() {}

    /** Initializing constructor. */
    NormalTexCoordAttribute(Normal const & n_, TexCoord const & t_) : NormalBaseT(n_), TexCoordBaseT(t_) {}

}; // class NormalTexCoordAttribute

/**
 * A utility attribute class that wraps both a color and a set of texture coordinates.
 *
 * @see ColorAttribute, TexCoordAttribute
 */
template <typename ColorT, typename TexCoordT>
struct ColorTexCoordAttribute : public ColorAttribute<ColorT>, public TexCoordAttribute<TexCoordT>
{
  private:
    typedef ColorAttribute<ColorT>        ColorBaseT;
    typedef TexCoordAttribute<TexCoordT>  TexCoordBaseT;

  public:
    THEA_DECL_SMART_POINTERS(ColorTexCoordAttribute)

    typedef typename ColorBaseT  ::Color      Color;    ///< The normal type.
    typedef typename TexCoordBaseT::TexCoord  TexCoord;  ///< The type of the texture coordinates vector.

    /** Default constructor. */
    ColorTexCoordAttribute() {}

    /** Initializing constructor. */
    ColorTexCoordAttribute(Color const & c_, TexCoord const & t_) : ColorBaseT(c_), TexCoordBaseT(t_) {}

}; // class ColorTexCoordAttribute

/**
 * A utility attribute class that wraps a normal, a color and a set of texture coordinates.
 *
 * @see NormalAttribute, ColorAttribute, TexCoordAttribute
 */
template <typename NormalT, typename ColorT, typename TexCoordT>
struct NormalColorTexCoordAttribute
: public NormalAttribute<NormalT>, public ColorAttribute<ColorT>, public TexCoordAttribute<TexCoordT>
{
  private:
    typedef NormalAttribute<NormalT>      NormalBaseT;
    typedef ColorAttribute<ColorT>        ColorBaseT;
    typedef TexCoordAttribute<TexCoordT>  TexCoordBaseT;

  public:
    THEA_DECL_SMART_POINTERS(NormalColorTexCoordAttribute)

    typedef NormalT    Normal;    ///< The normal type.
    typedef ColorT     Color;     ///< The color type.
    typedef TexCoordT  TexCoord;  ///< The type of the texture coordinates vector.

    /** Default constructor. */
    NormalColorTexCoordAttribute() {}

    /** Initializing constructor. */
    NormalColorTexCoordAttribute(Normal const & n_, Color const & c_, TexCoord const & t_)
    : NormalBaseT(n_), ColorBaseT(c_), TexCoordBaseT(t_) {}

}; // class NormalColorTexCoordAttribute

// Dirty SFINAE template trick to check if the type T has get...() and set...() functions
#define THEA_HAS_GRAPHICS_ATTRIB(name, get_func_name, set_func_name) \
template <typename T> \
class name \
{ \
  public: \
    typedef char Yes; \
    typedef struct { char a[2]; } No; \
\
    struct GetBaseMixin {  void get_func_name() {}  }; \
    struct SetBaseMixin {  void set_func_name() {}  }; \
\
    struct GetBase : public T, public GetBaseMixin {}; \
    struct SetBase : public T, public SetBaseMixin {}; \
    template <typename S, S>  class Helper {}; \
\
    template <typename U> static No testGet(U *, Helper<void (GetBaseMixin::*)(), &U::get_func_name> * = 0); \
    static Yes testGet(...); \
\
    template <typename U> static No testSet(U *, Helper<void (SetBaseMixin::*)(), &U::set_func_name> * = 0); \
    static Yes testSet(...); \
\
  public: \
    static bool const has_get = (sizeof(testGet((GetBase *)(0))) == sizeof(Yes)); \
    static bool const has_set = (sizeof(testSet((SetBase *)(0))) == sizeof(Yes)); \
    static bool const value = has_get && has_set; \
};

THEA_HAS_GRAPHICS_ATTRIB(HasPositionAttrib,  getPosition,  setPosition)
THEA_HAS_GRAPHICS_ATTRIB(HasNormalAttrib,    getNormal,    setNormal)
THEA_HAS_GRAPHICS_ATTRIB(HasColorAttrib,     getColor,     setColor)
THEA_HAS_GRAPHICS_ATTRIB(HasTexCoordAttrib,  getTexCoord,  setTexCoord)

#undef THEA_HAS_GRAPHICS_ATTRIB

// Short-hand for calling HasPositionAttrib<T::Attribute>
template <typename T>
struct HasPosition
{
  static bool const value = HasPositionAttrib<typename T::Attribute>::value;
};

// Short-hand for calling HasNormalAttrib<T::Attribute>
template <typename T>
struct HasNormal
{
  static bool const value = HasNormalAttrib<typename T::Attribute>::value;
};

// Short-hand for calling HasColorAttrib<T::Attribute>
template <typename T>
struct HasColor
{
  static bool const value = HasColorAttrib<typename T::Attribute>::value;
};

// Short-hand for calling HasTexCoordAttrib<T::Attribute>
template <typename T>
struct HasTexCoord
{
  static bool const value = HasTexCoordAttrib<typename T::Attribute>::value;
};

} // namespace Graphics
} // namespace Thea

#endif
