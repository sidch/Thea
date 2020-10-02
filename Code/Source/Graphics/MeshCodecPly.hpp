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
// First version: 2016
//
//============================================================================

#ifndef __Thea_Graphics_MeshCodecPly_hpp__
#define __Thea_Graphics_MeshCodecPly_hpp__

#include "../Common.hpp"
#include "../Array.hpp"
#include "MeshGroup.hpp"
#include "MeshCodec.hpp"
#include <algorithm>

namespace Thea {

namespace CodecPlyInternal {

template <typename MeshT, typename Enable = void>
struct VertexIndexMap
{
  typedef UnorderedMap<typename MeshT::Vertex const *, intx> type;
};

template <typename MeshT>
struct VertexIndexMap<MeshT, typename std::enable_if< Graphics::IsDisplayMesh<MeshT>::value >::type>
{
  typedef UnorderedMap<std::pair<MeshT const *, intx>, intx> type;
};

} // namespace CodecPlyInternal

/** %Codec for reading and writing Stanford PLY files. @see http://paulbourke.net/dataformats/ply/ */
template <typename MeshT, typename BuilderT>
class CodecPly : public CodecPlyBase<MeshT>
{
  private:
    typedef CodecPlyBase<MeshT> BaseT;
    typedef typename CodecPlyInternal::VertexIndexMap<MeshT>::type VertexIndexMap;

    /** Types of mesh elements (enum class). */
    struct ElementType
    {
      enum Value
      {
        VERTEX,
        FACE,
        UNKNOWN,
      };

      THEA_ENUM_CLASS_BODY(ElementType)

    }; // struct ElementType

    /** Type of an element property (enum class). */
    struct PropertyType
    {
      enum Value
      {
        // Relies on encodings in NumericType:
        //   - Right-two hex digits: number of bits in type
        //   - Third digit from right: signed (0) or unsigned (1)
        //   - Fourth digit from right: int (0) or float (1)
        INVALID     =   (int)NumericType::INVALID,
        LIST        =   0xFFFF,

        // New names
        INT8        =   (int)NumericType::INT8,
        INT16       =   (int)NumericType::INT16,
        INT32       =   (int)NumericType::INT32,
        UINT8       =   (int)NumericType::UINT8,
        UINT16      =   (int)NumericType::UINT16,
        UINT32      =   (int)NumericType::UINT32,
        FLOAT32     =   (int)NumericType::FLOAT32,
        FLOAT64     =   (int)NumericType::FLOAT64,
      };

      THEA_ENUM_CLASS_BODY(PropertyType)

      THEA_ENUM_CLASS_STRINGS_BEGIN(PropertyType)
        THEA_ENUM_CLASS_STRING(INVALID,   "invalid")
        THEA_ENUM_CLASS_STRING(LIST,      "list")

        // New names
        THEA_ENUM_CLASS_STRING(INT8,      "int8")
        THEA_ENUM_CLASS_STRING(INT16,     "int16")
        THEA_ENUM_CLASS_STRING(INT32,     "int32")
        THEA_ENUM_CLASS_STRING(UINT8,     "uint8")
        THEA_ENUM_CLASS_STRING(UINT16,    "uint16")
        THEA_ENUM_CLASS_STRING(UINT32,    "uint32")
        THEA_ENUM_CLASS_STRING(FLOAT32,   "float32")
        THEA_ENUM_CLASS_STRING(FLOAT64,   "float64")

        // Old names
        THEA_ENUM_CLASS_STRING(INT8,      "char")
        THEA_ENUM_CLASS_STRING(INT16,     "short")
        THEA_ENUM_CLASS_STRING(INT32,     "int")
        THEA_ENUM_CLASS_STRING(UINT8,     "uchar")
        THEA_ENUM_CLASS_STRING(UINT16,    "ushort")
        THEA_ENUM_CLASS_STRING(UINT32,    "uint")
        THEA_ENUM_CLASS_STRING(FLOAT32,   "float")
        THEA_ENUM_CLASS_STRING(FLOAT64,   "double")
      THEA_ENUM_CLASS_STRINGS_END(PropertyType)

    }; // struct PropertyType

    /** An element property. */
    struct Property
    {
      PropertyType type;

      // Only if type == LIST
      PropertyType count_type;
      PropertyType item_type;

    }; // struct Property

    /** Specification for a block of elements. */
    struct ElementBlock
    {
      ElementType type;           ///< Type of element.
      intx num_elems;             ///< Number of elements in the block.
      Array<Property> props;  ///< Element properties

    }; // struct ElementBlock

    /** Information in the header of a PLY file. */
    struct Header
    {
      Header() : binary(false), endianness(Endianness::LITTLE) {}

      bool binary;
      Endianness endianness;
      Array<ElementBlock> elem_blocks;

    }; // struct Header

  public:
    typedef MeshT Mesh;                                   ///< The type of mesh processed by the codec.
    typedef Graphics::MeshGroup<Mesh> MeshGroup;          ///< A group of meshes.
    typedef typename MeshGroup::MeshPtr MeshPtr;          ///< A shared pointer to a mesh.
    typedef BuilderT Builder;                             ///< The mesh builder class used by the codec.
    typedef typename BaseT::ReadCallback ReadCallback;    ///< Called when a mesh element is read.
    typedef typename BaseT::WriteCallback WriteCallback;  ///< Called when a mesh element is written.
    using BaseT::getName;

    /** %Options for deserializing meshes. */
    class ReadOptions
    {
      private:
        bool skip_empty_meshes;
        bool store_vertex_indices;
        bool store_face_indices;
        bool strict;
        bool verbose;

        friend class CodecPly;

      public:
        /** Constructor. Sets default values. */
        ReadOptions()
        : skip_empty_meshes(true), store_vertex_indices(true), store_face_indices(true), strict(false), verbose(false)
        {}

        /** Skip meshes with no faces? */
        ReadOptions & setSkipEmptyMeshes(bool value) { skip_empty_meshes = value; return *this; }

        /** Store vertex indices in mesh? */
        ReadOptions & setStoreVertexIndices(bool value) { store_vertex_indices = value; return *this; }

        /** Store face indices in mesh? */
        ReadOptions & setStoreFaceIndices(bool value) { store_face_indices = value; return *this; }

        /** Treat warnings as errors */
        ReadOptions & setStrict(bool value) { strict = value; return *this; }

        /** Print debugging information? */
        ReadOptions & setVerbose(bool value) { verbose = value; return *this; }

        /**
         * The set of default options. The default options correspond to
         * ReadOptions().setSkipEmptyMeshes(true).setStoreVertexIndices(true).setStoreFaceIndices(true).setVerbose(false).
         */
        static ReadOptions const & defaults() { static ReadOptions const def; return def; }

    }; // class ReadOptions

    /** %Options for serializing meshes. */
    class WriteOptions
    {
      private:
        bool binary;
        bool verbose;

        friend class CodecPly;

      public:
        /** Constructor. Sets default values. */
        WriteOptions() : binary(false), verbose(false) {}

        /** Write in the binary format? */
        WriteOptions & setBinary(bool value) { binary = value; return *this; }

        /** Print debugging information? */
        WriteOptions & setVerbose(bool value) { verbose = value; return *this; }

        /** The set of default options. The default options correspond to WriteOptions().setBinary(false).setVerbose(false). */
        static WriteOptions const & defaults() { static WriteOptions const def; return def; }

    }; // class WriteOptions

    /** Constructor. */
    CodecPly(ReadOptions const & read_opts_ = ReadOptions::defaults(),
             WriteOptions const & write_opts_ = WriteOptions::defaults())
    : read_opts(read_opts_), write_opts(write_opts_) {}

    void readMeshGroup(MeshGroup & mesh_group, BinaryInputStream & input, Codec::BlockHeader const * block_header,
                       ReadCallback * callback) const
    {
      mesh_group.clear();

      BinaryInputStream * in = &input;
      BinaryInputStream::Ptr tmp_in;
      if (block_header)
      {
        if (block_header->data_size <= 0) { return; }
        tmp_in = std::make_shared<BinaryInputStream>(input, (int64)block_header->data_size);
        in = tmp_in.get();
      }

      Header header;
      readHeader(header, *in);

      if (header.binary)
        readBinary(mesh_group, *in, header, callback);
      else
        readAscii(mesh_group, *in, header, callback);
    }

    void writeMeshGroup(MeshGroup const & mesh_group, BinaryOutputStream & output, bool write_block_header,
                        WriteCallback * callback) const
    {
      Codec::BlockHeader bh(this->getMagic());
      if (write_block_header)
        bh.markAndSkip(output);

      { BinaryOutputStream::EndiannessScope scope(output, Endianness::LITTLE);  // default to little-endian if in binary mode

        intx num_vertices = 0, num_faces = 0;
        getStats(mesh_group, num_vertices, num_faces);

        writeDefaultHeader(output, write_opts.binary, num_vertices, num_faces);

        VertexIndexMap vertex_indices;
        writeVertices(mesh_group, output, vertex_indices, callback);

        intx next_index = 0;
        writeFaces(mesh_group, vertex_indices, output, callback, next_index);
      }

      if (write_block_header)
        bh.calcAndWrite(output);
    }

  private:
    /** Read a property specification from a line. */
    void readPropertyHeader(Property & prop, std::istream & in) const
    {
      // Assume first word "property" has already been extracted
      std::string field;
      if (!(in >> field))
        throw Error(std::string(getName()) + ": Could not read property type");

      if (!prop.type.fromString(field))
        throw Error(std::string(getName()) + ": Unknown property type '" + field + '\'');

      if (prop.type == PropertyType::LIST)
      {
        std::string count_str, item_str;
        if (!(in >> count_str >> item_str))
          throw Error(std::string(getName()) + ": Could not read list item and count types");

        if (!prop.count_type.fromString(count_str))
          throw Error(std::string(getName()) + ": Unknown list count type '" + count_str + '\'');

        if (!prop.item_type.fromString(item_str))
          throw Error(std::string(getName()) + ": Unknown list item type '" + item_str + '\'');
      }
    }

    /** Read the header of a PLY file. */
    void readHeader(Header & header, BinaryInputStream & in) const
    {
      std::string magic = trimWhitespace(in.readLine());
      if (magic != "ply")
        throw Error(std::string(getName()) + ": Invalid PLY stream (does not start with 'ply')");

      header = Header();  // reset

      bool first = true;
      std::string line, field;
      while (in.hasMore())
      {
        line = trimWhitespace(in.readLine());

        if (line.empty())
          continue;

        if (line == "end_header")
          return;

        std::istringstream line_in(line);
        if (!(line_in >> field))
          throw Error(std::string(getName()) + ": Could not read first word of line");

        if (field == "format")
        {
          if (first)
          {
            if (!(line_in >> field))
              throw Error(std::string(getName()) + ": Could not read ascii/binary specifier");

            // Don't bother looking for the version number at the end of the line
            if (field == "ascii")
              header.binary = false;
            else if (field == "binary_little_endian")
            {
              header.binary = true;
              header.endianness = Endianness::LITTLE;
            }
            else if (field == "binary_big_endian")
            {
              header.binary = true;
              header.endianness = Endianness::BIG;
            }
            else
              throw Error(std::string(getName()) + ": Invalid format specifier");

            first = false;
          }
          else
            throw Error(std::string(getName()) + ": Format specifier must be first line after magic string");
        }
        else if (field == "element")
        {
          if (!(line_in >> field))
            throw Error(std::string(getName()) + ": Could not read element type");

          ElementBlock block;
          if (!(line_in >> block.num_elems))
            throw Error(std::string(getName()) + ": Could not read number of elements of type '" + field + '\'');

          if (field == "vertex")
            block.type = ElementType::VERTEX;
          else if (field == "face")
            block.type = ElementType::FACE;
          else
            block.type = ElementType::UNKNOWN;

          header.elem_blocks.push_back(block);
        }
        else if (field == "property")
        {
          if (header.elem_blocks.empty())
            throw Error(std::string(getName()) + ": Property specification outside an element block");

          Property prop;
          readPropertyHeader(prop, line_in);
          if (prop.type == PropertyType::INVALID)
            throw Error(std::string(getName()) + ": Invalid property type");

          header.elem_blocks.back().props.push_back(prop);
        }
      }

      throw Error(std::string(getName()) + ": No end_header token, header not closed");
    }

    /** Sanity checks for an element block. */
    void checkBlock(ElementBlock const & block) const
    {
      if (block.type == ElementType::VERTEX)
      {
        if (block.props.size() < 3
          || block.props[0].type == PropertyType::LIST
          || block.props[1].type == PropertyType::LIST
          || block.props[2].type == PropertyType::LIST)
          throw Error(std::string(getName()) + ": Vertex element must start with three numerical coordinates");
      }
      else if (block.type == ElementType::FACE)
      {
        if (block.props.empty()
          || block.props[0].type != PropertyType::LIST
          || (block.props[0].count_type & 0xF000)
          || (block.props[0].item_type  & 0xF000))
          throw Error(std::string(getName()) + ": Face element must start with a list of integer indices");
      }
    }

    /** Read a mesh group in ASCII format. */
    void readAscii(MeshGroup & mesh_group, BinaryInputStream & in, Header const & header, ReadCallback * callback) const
    {
      // Create new mesh
      MeshPtr mesh(new Mesh(std::string(mesh_group.getName()) + "/Mesh"));

      // Create a builder for the mesh
      Builder builder(mesh);
      builder.begin();

      Array<typename Builder::VertexHandle> vrefs;
      Array<typename Builder::VertexHandle> face;

      std::string line;
      intx num_vertices = 0, num_faces = 0;
      for (size_t i = 0; i < header.elem_blocks.size(); ++i)
      {
        ElementBlock const & block = header.elem_blocks[i];
        checkBlock(block);

        for (intx j = 0; j < block.num_elems; ++j)
        {
          do
          {
            if (!in.hasMore())
              throw Error(std::string(getName()) + ": Unexpected end of input");

            line = trimWhitespace(in.readLine());

          } while (line.empty() || beginsWith(line, "comment"));

          switch (block.type)
          {
            case ElementType::VERTEX:
            {
              double x, y, z;
              std::istringstream vstr(line);
              if (!(vstr >> x >> y >> z))
                throw Error(std::string(getName()) + ": Could not read vertex on line '" + line + '\'');

              typename Builder::VertexHandle vref = builder.addVertex(Vector3((Real)x, (Real)y, (Real)z),
                                                                      (read_opts.store_vertex_indices ? num_vertices : -1));
              if (callback)
                callback->vertexRead(mesh.get(), num_vertices, vref);

              vrefs.push_back(vref);
              num_vertices++;

              break;
            }

            case ElementType::FACE:
            {
              intx index, num_face_vertices;
              std::istringstream vstr(line);
              if (!(vstr >> num_face_vertices))
                throw Error(std::string(getName()) + ": Could not read number of vertices in face on line '" + line + '\'');

              if (num_face_vertices > 0)
              {
                face.resize(num_face_vertices);

                bool skip = false;
                for (intx v = 0; v < num_face_vertices && !skip; ++v)
                {
                  if (!(vstr >> index))
                  {
                    if (read_opts.strict)
                      throw Error(std::string(getName()) + ": Could not read vertex index on line '" + line + '\'');
                    else
                    {
                      THEA_WARNING << getName() << ": Skipping face, could not read vertex index on line '" << line << '\'';
                      skip = true; break;
                    }
                  }

                  if (index < 0 || index >= (intx)vrefs.size())
                  {
                    if (read_opts.strict)
                      throw Error(getName() + format(": Vertex index %ld out of bounds (#vertices = %ld) on line '%s'",
                                                     index, (intx)vrefs.size(), line.c_str()));
                    else
                    {
                      THEA_WARNING << getName() << ": Vertex index " << index << " out of bounds (#vertices = "
                                   << vrefs.size() << ") on line '" << line << '\'';
                      skip = true; break;
                    }
                  }

                  face[(size_t)v] = vrefs[(size_t)index];

                  for (int w = 0; w < v; ++w)
                    if (face[w] == face[v])  // face has repeated vertices
                    {
                      if (read_opts.strict)
                        throw Error(std::string(getName()) + ": Face has repeated vertices on line '" + line + '\'');
                      else
                      {
                        THEA_WARNING << getName() << ": Skipping face with repeated vertices on line '" << line << '\'';
                        skip = true; break;
                      }
                    }
                }

                if (!skip)
                {
                  typename Builder::FaceHandle fref = builder.addFace(face.begin(), face.end(),
                                                                      (read_opts.store_face_indices ? num_faces : -1));
                  if (callback)
                    callback->faceRead(mesh.get(), num_faces, fref);

                  num_faces++;
                }
              }

              break;
            }

            default: {}
          }
        }
      }

      builder.end();

      if (builder.numVertices() > 0 || builder.numFaces() > 0)
        mesh_group.addMesh(mesh);

      THEA_CONSOLE << getName() << ": Read mesh with " << num_vertices << " vertices and " << num_faces << " faces";
    }

    /** Read a binary-encoded number and return it as a specified type. */
    template <typename T> T readBinaryNumber(BinaryInputStream & in, PropertyType const & type) const
    {
      switch (type)
      {
        case PropertyType::INT8:      return static_cast<T>(in.readInt8());
        case PropertyType::INT16:     return static_cast<T>(in.readInt16());
        case PropertyType::INT32:     return static_cast<T>(in.readInt32());
        case PropertyType::UINT8:     return static_cast<T>(in.readUInt8());
        case PropertyType::UINT16:    return static_cast<T>(in.readUInt16());
        case PropertyType::UINT32:    return static_cast<T>(in.readUInt32());
        case PropertyType::FLOAT32:   return static_cast<T>(in.readFloat32());
        case PropertyType::FLOAT64:   return static_cast<T>(in.readFloat64());
        default: throw Error(std::string(getName()) + ": Unknown property type");
      }
    }

    /** Return a binary-encoded list of elements of a specified type. */
    template <typename T> void readBinaryList(BinaryInputStream & in, Property const & prop, Array<T> & items) const
    {
      debugAssertM(prop.type == PropertyType::LIST, std::string(getName()) + ": Can't read non-list property as list");

      intx num_items = readBinaryNumber<intx>(in, prop.count_type);
      if (num_items < 0)
        throw Error(std::string(getName()) + ": List has negative size");

      items.resize((size_t)num_items);
      for (size_t i = 0; i < items.size(); ++i)
        items[i] = readBinaryNumber<T>(in, prop.item_type);
    }

    /** Skip over a binary-encoded property. */
    void skipBinaryProperty(BinaryInputStream & in, Property const & prop) const
    {
      if (prop.type == PropertyType::LIST)
      {
        intx num_items = readBinaryNumber<intx>(in, prop.count_type);
        if (num_items >= 0)
          in.skip(num_items * (prop.item_type & 0xFF));
      }
      else
        in.skip(prop.type & 0xFF);
    }

    /** Read a mesh group in binary format. */
    void readBinary(MeshGroup & mesh_group, BinaryInputStream & in, Header const & header, ReadCallback * callback) const
    {
      // Create new mesh
      MeshPtr mesh(new Mesh(std::string(mesh_group.getName()) + "/Mesh"));

      // Create a builder for the mesh
      Builder builder(mesh);
      builder.begin();

      Array<typename Builder::VertexHandle> vrefs;
      Array<typename Builder::VertexHandle> face;

      BinaryInputStream::EndiannessScope scope(in, header.endianness);

      intx num_vertices = 0, num_faces = 0;
      for (size_t i = 0; i < header.elem_blocks.size(); ++i)
      {
        ElementBlock const & block = header.elem_blocks[i];
        checkBlock(block);

        for (intx j = 0; j < block.num_elems; ++j)
        {
          switch (block.type)
          {
            case ElementType::VERTEX:
            {
              Vector3 vertex;
              vertex[0] = readBinaryNumber<Real>(in, block.props[0].type);
              vertex[1] = readBinaryNumber<Real>(in, block.props[1].type);
              vertex[2] = readBinaryNumber<Real>(in, block.props[2].type);

              for (size_t k = 3; k < block.props.size(); ++i)
                skipBinaryProperty(in, block.props[k]);

              typename Builder::VertexHandle vref = builder.addVertex(vertex,
                                                                      (read_opts.store_vertex_indices ? num_vertices : -1));
              if (callback)
                callback->vertexRead(mesh.get(), num_vertices, vref);

              vrefs.push_back(vref);
              num_vertices++;

              break;
            }

            case ElementType::FACE:
            {
              Array<intx> face_vertices;
              readBinaryList(in, block.props[0], face_vertices);

              if (!face_vertices.empty())
              {
                face.resize(face_vertices.size());

                bool skip = false;
                for (size_t v = 0; v < face_vertices.size() && !skip; ++v)
                {
                  intx index = face_vertices[v];
                  if (index < 0 || index >= (intx)vrefs.size())
                  {
                    if (read_opts.strict)
                      throw Error(getName() + format(": Vertex index %ld out of bounds (#vertices = %ld) in face %ld",
                                                     index, (intx)vrefs.size(), num_faces));
                    else
                    {
                      THEA_WARNING << getName() << ": Skipping face, vertex index " << index << " out of bounds (#vertices = "
                                                << vrefs.size() << ')';
                      skip = true; break;
                    }
                  }

                  face[v] = vrefs[(size_t)index];

                  for (size_t w = 0; w < v; ++w)
                    if (face[w] == face[v])  // face has repeated vertices
                    {
                      if (read_opts.strict)
                        throw Error(getName() + format(": Face %ld has repeated vertices", num_faces));
                      else
                      {
                        THEA_WARNING << getName() << ": Skipping face with repeated vertices";
                        skip = true; break;
                      }
                    }
                }

                if (!skip)
                {
                  typename Builder::FaceHandle fref = builder.addFace(face.begin(), face.end(),
                                                                      (read_opts.store_face_indices ? num_faces : -1));
                  if (callback)
                    callback->faceRead(mesh.get(), num_faces, fref);

                  num_faces++;
                }
              }

              break;
            }

            default: {}
          }
        }
      }

      builder.end();

      if (builder.numVertices() > 0 || builder.numFaces() > 0)
        mesh_group.addMesh(mesh);

      THEA_CONSOLE << getName() << ": Read mesh with " << num_vertices << " vertices and " << num_faces << " faces";
    }

    /** Count the number of vertices and faces in a mesh. */
    template <typename _MeshT>
    static void getStats(_MeshT const & mesh, intx & num_vertices, intx & num_faces)
    {
      num_vertices += mesh.numVertices();
      num_faces += mesh.numFaces();
    }

    /** Count the number of vertices and faces in a mesh group (increments parameters). */
    static void getStats(MeshGroup const & mesh_group, intx & num_vertices, intx & num_faces)
    {
      for (typename MeshGroup::MeshConstIterator mi = mesh_group.meshesBegin(); mi != mesh_group.meshesEnd(); ++mi)
        getStats(**mi, num_vertices, num_faces);

      for (typename MeshGroup::GroupConstIterator ci = mesh_group.childrenBegin(); ci != mesh_group.childrenEnd(); ++ci)
        getStats(**ci, num_vertices, num_faces);
    }

    /** Write the default header of a PLY file. */
    void writeDefaultHeader(BinaryOutputStream & out, bool binary, intx num_vertices, intx num_faces) const
    {
      out.printf("ply\n");

      if (binary) out.printf("format binary_little_endian 1.0\n");  // default to little-endian output
      else        out.printf("format ascii 1.0\n");

      out.printf("element vertex %ld\n", num_vertices);
      out.printf("property float x\n");  // stick to old typenames for compatibility
      out.printf("property float y\n");
      out.printf("property float z\n");

      out.printf("element face %ld\n", num_faces);
      out.printf("property list int int vertex_indices\n");

      out.printf("end_header\n");
    }

    /** Write out all the vertices from a mesh group and map them to indices. */
    void writeVertices(MeshGroup const & mesh_group, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                       WriteCallback * callback) const
    {
      for (typename MeshGroup::MeshConstIterator mi = mesh_group.meshesBegin(); mi != mesh_group.meshesEnd(); ++mi)
      {
        writeVertices(**mi, output, vertex_indices, callback);
      }

      for (typename MeshGroup::GroupConstIterator ci = mesh_group.childrenBegin(); ci != mesh_group.childrenEnd(); ++ci)
      {
        writeVertices(**ci, output, vertex_indices, callback);
      }
    }

    /** Write out all the vertices from a general or DCEL mesh and map them to indices. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsGeneralMesh<_MeshT>::value
                                                      || Graphics::IsDcelMesh<_MeshT>::value, int >::type = 0 >
    void writeVertices(_MeshT const & mesh, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                       WriteCallback * callback) const
    {
      intx vertex_index = (intx)vertex_indices.size();
      for (typename Mesh::VertexConstIterator vi = mesh.verticesBegin(); vi != mesh.verticesEnd(); ++vi, ++vertex_index)
      {
        if (write_opts.binary)
        {
          output.writeFloat32((float32)vi->getPosition().x());
          output.writeFloat32((float32)vi->getPosition().y());
          output.writeFloat32((float32)vi->getPosition().z());
        }
        else
          output.printf("%f %f %f\n", vi->getPosition().x(), vi->getPosition().y(), vi->getPosition().z());

        vertex_indices[&(*vi)] = vertex_index;
        if (callback) callback->vertexWritten(&mesh, vertex_index, &(*vi));
      }
    }

    /** Write out all the vertices from a display mesh and map them to indices. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsDisplayMesh<_MeshT>::value, int >::type = 0 >
    void writeVertices(_MeshT const & mesh, BinaryOutputStream & output, VertexIndexMap & vertex_indices,
                       WriteCallback * callback) const
    {
      typedef std::pair<_MeshT const *, intx> DisplayMeshVRef;
      typename Mesh::VertexArray const & vertices = mesh.getVertices();
      intx vertex_index = (intx)vertex_indices.size();

      for (size_t i = 0; i < vertices.size(); ++i, ++vertex_index)
      {
        Vector3 const & v = vertices[i];

        if (write_opts.binary)
        {
          output.writeFloat32((float32)v.x());
          output.writeFloat32((float32)v.y());
          output.writeFloat32((float32)v.z());
        }
        else
          output.printf("%f %f %f\n", v.x(), v.y(), v.z());

        vertex_indices[DisplayMeshVRef(&mesh, (intx)i)] = vertex_index;
        if (callback) callback->vertexWritten(&mesh, vertex_index, (intx)i);
      }
    }

    /** Write out all the faces from a mesh group. */
    void writeFaces(MeshGroup const & mesh_group, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                    WriteCallback * callback, intx & next_index) const
    {
      for (typename MeshGroup::MeshConstIterator mi = mesh_group.meshesBegin(); mi != mesh_group.meshesEnd(); ++mi)
      {
        writeFaces(**mi, vertex_indices, output, callback, next_index);
      }

      for (typename MeshGroup::GroupConstIterator ci = mesh_group.childrenBegin(); ci != mesh_group.childrenEnd(); ++ci)
      {
        writeFaces(**ci, vertex_indices, output, callback, next_index);
      }
    }

    /** Write out all the faces from a general or DCEL mesh. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsGeneralMesh<_MeshT>::value
                                                      || Graphics::IsDcelMesh<_MeshT>::value, int >::type = 0 >
    void writeFaces(_MeshT const & mesh, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                    WriteCallback * callback, intx & next_index) const
    {
      for (typename Mesh::FaceConstIterator fi = mesh.facesBegin(); fi != mesh.facesEnd(); ++fi)
      {
        typename Mesh::Face const & face = *fi;
        if (face.numVertices() < 3) continue;

        if (write_opts.binary)
        {
          output.writeInt32((int32)face.numVertices());
          for (typename Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
          {
            typename VertexIndexMap::const_iterator ii = vertex_indices.find(*vi);
            alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

            output.writeInt32((int32)ii->second);
          }
        }
        else
        {
          std::ostringstream os; os << face.numVertices();
          for (typename Mesh::Face::VertexConstIterator vi = face.verticesBegin(); vi != face.verticesEnd(); ++vi)
          {
            typename VertexIndexMap::const_iterator ii = vertex_indices.find(*vi);
            alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

            os << ' ' << ii->second;
          }

          os << '\n';
          output.writeBytes((int64)os.str().length(), os.str().data());
        }

        if (callback) callback->faceWritten(&mesh, next_index++, &face);
      }
    }

    /** Write out all the faces from a display mesh. */
    template < typename _MeshT, typename std::enable_if< Graphics::IsDisplayMesh<_MeshT>::value, int >::type = 0 >
    void writeFaces(_MeshT const & mesh, VertexIndexMap const & vertex_indices, BinaryOutputStream & output,
                    WriteCallback * callback, intx & next_index) const
    {
      typedef std::pair<_MeshT const *, intx> DisplayMeshVRef;

      for (int type = 0; type < 2; ++type)  // 0: triangles, 1: quads
      {
        typename Mesh::IndexArray indices = (type == 0 ? mesh.getTriangleIndices() : mesh.getQuadIndices());
        size_t degree = (type == 0 ? 3 : 4);

        for (size_t i = 0; i < indices.size(); i += degree)
        {
          if (write_opts.binary)
          {
            output.writeInt32((int32)degree);
            for (size_t j = 0; j < degree; ++j)
            {
              typename VertexIndexMap::const_iterator ii = vertex_indices.find(DisplayMeshVRef(&mesh, (intx)indices[i + j]));
              alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

              output.writeInt32((int32)ii->second);
            }
          }
          else
          {
            std::ostringstream os; os << degree;

            for (size_t j = 0; j < degree; ++j)
            {
              typename VertexIndexMap::const_iterator ii = vertex_indices.find(DisplayMeshVRef(&mesh, (intx)indices[i + j]));
              alwaysAssertM(ii != vertex_indices.end(), std::string(getName()) + ": Vertex index not found");

              os << ' ' << ii->second;
            }

            os << '\n';
            output.writeBytes((int64)os.str().length(), os.str().data());
          }

          if (callback)
          {
            typename Mesh::Face face(const_cast<Mesh *>(&mesh), degree, (type == 0), (intx)i, 1);
            callback->faceWritten(&mesh, next_index++, face);
          }
        }
      }
    }

    ReadOptions read_opts;
    WriteOptions write_opts;

}; // class CodecPly

} // namespace Thea

#endif
