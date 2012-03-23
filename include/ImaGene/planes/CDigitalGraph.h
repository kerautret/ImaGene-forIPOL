/** @file CDigitalGraph.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : CDigitalGraph.h
//
// Creation : 2011/01/13
//
// Version : 2011/01/13
//
// Author : JOL
//
// Summary : Header file for files CDigitalGraph.ih and CDigitalGraph.cxx
//
// History :
//	2011/01/13 : ?Name? : ?What?
//
// Rcs Id : "@(#)class CDigitalGraph declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(CDigitalGraph_RECURSES)
#error Recursive header files inclusion detected in CDigitalGraph.h
#else // defined(CDigitalGraph_RECURSES)
#define CDigitalGraph_RECURSES

#if !defined CDigitalGraph_h
#define CDigitalGraph_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <set>
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  
  /////////////////////////////////////////////////////////////////////////////
  // class CDigitalGraph
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'CDigitalGraph' <p>
   * Aim: The concept of a graph whose vertices are embedded in some Z^n.
   */
  class CDigitalGraph
  {
    // ----------------------- Needed types ------------------------------
  public:
    /**
     * The type for each vertex.
     */
    typedef int Vertex;

    /**
     * An archetype for a neighborhood. Must be iterable.
     */
    struct CNeighborhood : public std::vector<Vertex>
    {
      Vertex p;
      const Vertex & center() const;
      // also has begin() and end() services.
    };

    /**
     * The type for the (proper or not) neighborhood of a vertex
     */
    typedef CNeighborhood Neighborhood;

    /**
     * The type for an arbitrary large set of vertices. Such a Set
     * should have the following operations: copy-on-write mechanism,
     * operator[] to access one element, operators +=( elem ) and -=(
     * elem ) to add or remove some element.
     */
    typedef std::set<Vertex> Set;

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~CDigitalGraph();

    /**
     * @return the empty set of vertices (subset of this graph).
     */
    Set emptySet() const;

    /**
     * @return the whole set of vertices (subset of this graph).
     */
    Set wholeSet() const;

    /**
     * Sets the neighborhood of [p] in [n].
     *
     * @param p any vertex of this graph.
     * @param n (modified) contains after the neighborhood of p.
     */
    void setNeighborhood( const Vertex & p, Neighborhood & n ) const;

    /**
     * Sets the proper neighborhood of [p] in [n].
     *
     * @param p any vertex of this graph.
     * @param n (modified) contains after the proper neighborhood of p.
     */
    void setProperNeighborhood( const Vertex & p, Neighborhood & n ) const;

    /**
     * @return the dimension of the embedding space.
     */
    int dim() const;

    /**
     * @return the diameter of the graph in its embedding space.
     */
    int diameter() const;

    /**
     * Embeds the vertex [p] in the digital space.
     *
     * @param p any vertex of this graph.
     * @param coords (modified) the integer coordinates of the embedding of [p].
     */
    void embed( const Vertex & p, int* coords ) const;

    /**
       @param v any vertex.
       @return its preferred axis when computing its disk.
    */
    unsigned int preferredAxis( const Vertex & v ) const;

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param that_stream the output stream where the object is written.
     */
    void selfDisplay( std::ostream & that_stream ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool OK() const;
  

    // ------------------------- Datas ----------------------------------------
  private:


    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    INLINE CDigitalGraph();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE CDigitalGraph( const CDigitalGraph & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE CDigitalGraph & operator=( const CDigitalGraph & other );
  
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'CDigitalGraph'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'CDigitalGraph' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const CDigitalGraph & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/planes/CDigitalGraph.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CDigitalGraph_h

#undef CDigitalGraph_RECURSES
#endif // else defined(CDigitalGraph_RECURSES)
