#ifndef __Thea_SegmentSDF_Graph_hpp__
#define __Thea_SegmentSDF_Graph_hpp__

#include "../../AttributedObject.hpp"
#include <set>

/** Directed graph. */
template <typename VertexAttribute, typename EdgeAttribute>
class Graph
{
  private:
    /** An intrusive list similar to std::list that assumes the class T has <code>prev</code> and <code>next</code> pointers. */
    template <typename T> class IntrusiveList
    {
      private:
        T * first;
        T * last;
        size_t num_elems;

      public:
        template <typename P> class basic_iterator
        {
          private:
            P * p;

          public:
            basic_iterator() : p(NULL) {}

            basic_iterator(P * p_) : p(p_) {}

            basic_iterator(basic_iterator const & src) : p(src.p) {}

            basic_iterator operator++(int dummy)
            {
              if (p)
              {
                P * op = p;
                p = p->next;
                return basic_iterator(op);
              }
              else
                return *this;
            }

            basic_iterator & operator++()
            {
              if (p) p = p->next;
              return *this;
            }

            bool operator==(basic_iterator const & other) const { return p == other.p; }
            bool operator!=(basic_iterator const & other) const { return p != other.p; }

            template <typename Q> basic_iterator & operator=(basic_iterator<Q> const & src) { p = src.p; return *this; }

            P & operator*() const { return *p; }
            P * operator->() const { return p; }
        };

        typedef basic_iterator<T> iterator;
        typedef basic_iterator<T const> const_iterator;

        IntrusiveList() : first(NULL), last(NULL), num_elems(0) {}

        const_iterator begin() const { return const_iterator(first); }
        iterator begin() { return iterator(first); }
        const_iterator end() const { static const_iterator const e; return e; }
        iterator end() { static iterator const e; return e; }

        void push_back(T & t)
        {
          if (empty())
          {
            first = last = &t;
            t.prev = t.next = NULL;
          }
          else
          {
            last->next = &t;
            t.prev = last;
            t.next = NULL;
            last = &t;
          }

          num_elems++;
        }

        void push_front(T & t)
        {
          if (empty())
          {
            first = last = &t;
            t.prev = t.next = NULL;
          }
          else
          {
            first->prev = &t;
            t.next = first;
            t.prev = NULL;
            first = &t;
          }

          num_elems++;
        }

        void insert(iterator pos, T & t)
        {
          if (pos == end())
            push_back(t);
          else
          {
            t.next = pos->next;
            t.next->prev = &t;
            pos->next = &t;
            t.prev = &(*pos);
            num_elems++;
          }
        }

        void erase(iterator pos)
        {
          if (pos != end())
          {
            if (pos->prev)
              pos->prev->next = pos->next;
            else
              first = pos->next;

            if (pos->next)
              pos->next->prev = pos->prev;
            else
              last = pos->prev;

            num_elems--;
          }
        }

        size_t size() const { return num_elems; }
        bool empty() const { return num_elems == 0; }

    }; // class IntrusiveList

  public:
    class Vertex;  // forward declaration

    /** Graph edge. */
    class Edge : public Thea::AttributedObject<EdgeAttribute>
    {
      private:
        typedef Thea::AttributedObject<EdgeAttribute> BaseType;

      public:
        Vertex * origin, * end;
        Edge * prev;
        Edge * next;

        Edge(Vertex * origin_, Vertex * end_, EdgeAttribute const & attrib)
        : BaseType(attrib), origin(origin_), end(end_), prev(NULL), next(NULL) {}

        Edge() : origin(NULL), end(NULL), prev(NULL), next(NULL) {}

        Vertex const * getOrigin() const { return origin; }
        Vertex * getOrigin() { return origin; }

        Vertex const * getEnd() const { return end; }
        Vertex * getEnd() { return end; }

    }; // class Edge

    /** Graph vertex. */
    class Vertex : public Thea::AttributedObject<VertexAttribute>
    {
      private:
        typedef Thea::AttributedObject<VertexAttribute> BaseType;

      public:
        typedef std::set<Edge *> EdgeSet;

        EdgeSet incoming_edges;
        EdgeSet outgoing_edges;
        Vertex * prev;
        Vertex * next;

        typedef typename EdgeSet::iterator EdgeIterator;
        typedef typename EdgeSet::const_iterator EdgeConstIterator;

        Vertex(VertexAttribute const & attrib) : BaseType(attrib), prev(NULL), next(NULL) {}
        Vertex() : prev(NULL), next(NULL) {}

        EdgeConstIterator incomingEdgesBegin() const { return incoming_edges.begin(); }
        EdgeIterator incomingEdgesBegin() { return incoming_edges.begin(); }
        EdgeConstIterator incomingEdgesEnd() const { return incoming_edges.end(); }
        EdgeIterator incomingEdgesEnd() { return incoming_edges.end(); }

        EdgeConstIterator outgoingEdgesBegin() const { return outgoing_edges.begin(); }
        EdgeIterator outgoingEdgesBegin() { return outgoing_edges.begin(); }
        EdgeConstIterator outgoingEdgesEnd() const { return outgoing_edges.end(); }
        EdgeIterator outgoingEdgesEnd() { return outgoing_edges.end(); }

        unsigned long numIncomingEdges() const { return (unsigned long)incoming_edges.size(); }
        unsigned long numOutgoingEdges() const { return (unsigned long)outgoing_edges.size(); }
        unsigned long numEdges() const { return (unsigned long)(incoming_edges.size() + outgoing_edges.size()); }

    }; // class Vertex

  private:
    typedef IntrusiveList<Vertex> VertexList;
    typedef IntrusiveList<Edge> EdgeList;

    VertexList vertices;
    EdgeList edges;

  public:
    typedef typename VertexList::const_iterator VertexConstIterator;
    typedef typename VertexList::iterator VertexIterator;
    typedef typename EdgeList::const_iterator EdgeConstIterator;
    typedef typename EdgeList::iterator EdgeIterator;

    ~Graph() { clear(); }

    void clear()
    {
      if (vertices.empty() && edges.empty()) return;

      for (VertexIterator vi = verticesBegin(); vi != verticesEnd(); )
        delete &(*(vi++));

      for (EdgeIterator ei = edgesBegin(); ei != edgesEnd(); )
        delete &(*(ei++));

      // Reinitialize
      vertices = VertexList();
      edges = EdgeList();
    }

    VertexConstIterator verticesBegin() const { return vertices.begin(); }
    VertexIterator verticesBegin() { return vertices.begin(); }
    VertexConstIterator verticesEnd() const { return vertices.end(); }
    VertexIterator verticesEnd() { return vertices.end(); }

    EdgeConstIterator edgesBegin() const { return edges.begin(); }
    EdgeIterator edgesBegin() { return edges.begin(); }
    EdgeConstIterator edgesEnd() const { return edges.end(); }
    EdgeIterator edgesEnd() { return edges.end(); }

    unsigned long numVertices() const { return (unsigned long)vertices.size(); }
    unsigned long numEdges() const { return (unsigned long)edges.size(); }

    VertexIterator addVertex(VertexAttribute const & attrib)
    {
      Vertex * v = new Vertex(attrib);
      vertices.push_front(*v);

      return vertices.begin();
    }

    EdgeIterator addEdge(VertexIterator origin, VertexIterator end, EdgeAttribute const & attrib)
    {
      debugAssertM(origin != verticesEnd() && end != verticesEnd(), "Graph: Cannot add edge with null endpoint");

      Edge * e = new Edge(&(*origin), &(*end), attrib);
      edges.push_front(*e);

      origin->outgoing_edges.insert(e);
      end->incoming_edges.insert(e);

      return edges.begin();
    }

    void removeVertex(VertexIterator v)
    {
      if (v == verticesEnd()) return;

      for (typename Vertex::EdgeIterator ei = v->incomingEdgesBegin(); ei != v->incomingEdgesEnd(); )
        removeEdge(EdgeIterator(*(ei++)));

      for (typename Vertex::EdgeIterator ei = v->outgoingEdgesBegin(); ei != v->outgoingEdgesEnd(); )
        removeEdge(EdgeIterator(*(ei++)));

      debugAssertM(v->numEdges() == 0, "Graph: Vertex has at least one edge, where it is expected to have zero");

      vertices.erase(v);
      delete &(*v);
    }

    void removeEdge(EdgeIterator e)
    {
      e->getOrigin()->outgoing_edges.erase(&(*e));
      e->getEnd()->incoming_edges.erase(&(*e));

      edges.erase(e);
      delete &(*e);
    }

    VertexIterator collapseEdge(EdgeIterator e)
    {
      if (e == edgesEnd()) return verticesEnd();

      Vertex * vo = e->getOrigin();
      Vertex * ve = e->getEnd();

      if (vo == ve)  // self-loop
      {
        removeEdge(e);
        return VertexIterator(vo);
      }

      for (typename Vertex::EdgeIterator ei = ve->incomingEdgesBegin(); ei != ve->incomingEdgesEnd(); ++ei)
      {
        if (*ei != &(*e))
        {
          (*ei)->end = vo;
          vo->incoming_edges.insert(*ei);
        }
      }

      for (typename Vertex::EdgeIterator ei = ve->outgoingEdgesBegin(); ei != ve->outgoingEdgesEnd(); ++ei)
      {
        (*ei)->origin = vo;
        vo->outgoing_edges.insert(*ei);
      }

      vertices.erase(VertexIterator(ve));
      delete &(*ve);

      vo->outgoing_edges.erase(&(*e));
      edges.erase(e);
      delete &(*e);

      return VertexIterator(vo);
    }

    /** Merge each pair of twin edges incident at a vertex, whether parallel or opposing. */
    void mergeTwinEdges(VertexIterator v)
    {
      mergeParallelTwinEdges(v);
      mergeOpposingTwinEdges(v);
    }

    /** Merge double edges in the same direction incident at a vertex. */
    void mergeParallelTwinEdges(VertexIterator v)
    {
      if (v == verticesEnd()) return;

      // Compare incoming edges with incoming edges
      for (typename Vertex::EdgeIterator ei = v->incomingEdgesBegin(); ei != v->incomingEdgesEnd(); ++ei)
      {
        for (typename Vertex::EdgeIterator ej = v->incomingEdgesBegin(); ej != ei; ++ej)
          if ((*ei)->getOrigin() == (*ej)->getOrigin())
          {
            removeEdge(EdgeIterator(*ej));
            break;
          }
      }

      // Compare outgoing edges with outgoing edges
      for (typename Vertex::EdgeIterator ei = v->outgoingEdgesBegin(); ei != v->outgoingEdgesEnd(); ++ei)
      {
        for (typename Vertex::EdgeIterator ej = v->outgoingEdgesBegin(); ej != ei; ++ej)
          if ((*ei)->getEnd() == (*ej)->getEnd())
          {
            removeEdge(EdgeIterator(*ej));
            break;
          }
      }
    }

    /**
     * Merge double edges in opposite directions incident at a vertex. The incoming twin is always removed, so a single outgoing
     * edge remains.
     */
    void mergeOpposingTwinEdges(VertexIterator v)
    {
      if (v == verticesEnd()) return;

      for (typename Vertex::EdgeIterator ei = v->outgoingEdgesBegin(); ei != v->outgoingEdgesEnd(); ++ei)
      {
        for (typename Vertex::EdgeIterator ej = v->incomingEdgesBegin(); ej != v->incomingEdgesEnd(); ++ej)
          if ((*ei)->getEnd() == (*ej)->getOrigin())
          {
            if ((*ej)->getEnd() == (*ej)->getOrigin())
              std::cout << "WARNING: Removing self-loop!!!" << std::endl;

            removeEdge(EdgeIterator(*ej));
            break;
          }
      }
    }

}; // class Graph

template <typename MyGraph>
void
printGraph(MyGraph const & graph)
{
  std::cout << "Graph has " << graph.numEdges() << " edges and " << graph.numVertices() << " vertices" << std::endl;
  for (typename MyGraph::VertexConstIterator vi = graph.verticesBegin(); vi != graph.verticesEnd(); ++vi)
  {
    std::cout << "Vertex: " << vi->attr() << std::endl;
    std::cout << "  Outgoing: " << std::endl;
    for (typename MyGraph::Vertex::EdgeConstIterator ei = vi->outgoingEdgesBegin(); ei != vi->outgoingEdgesEnd(); ++ei)
      std::cout << "    " << (*ei)->attr() << std::endl;

    std::cout << "  Incoming: " << std::endl;
    for (typename MyGraph::Vertex::EdgeConstIterator ei = vi->incomingEdgesBegin(); ei != vi->incomingEdgesEnd(); ++ei)
      std::cout << "    " << (*ei)->attr() << std::endl;
  }

  for (typename MyGraph::EdgeConstIterator ei = graph.edgesBegin(); ei != graph.edgesEnd(); ++ei)
    std::cout << "Edge(" << ei->getOrigin()->attr() << ", " << ei->getEnd()->attr() << "): "
              << ei->attr() << std::endl;
}

#endif
