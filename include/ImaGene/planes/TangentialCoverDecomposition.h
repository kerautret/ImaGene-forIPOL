/** @file TangentialCoverDecomposition.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : TangentialCoverDecomposition.h
//
// Creation : 2011/02/15
//
// Version : 2011/02/15
//
// Author : JOL
//
// Summary : Header file for files TangentialCoverDecomposition.ih and TangentialCoverDecomposition.cxx
//
// History :
//	2011/02/15 : ?Name? : ?What?
//
// Rcs Id : "@(#)class TangentialCoverDecomposition declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(TangentialCoverDecomposition_RECURSES)
#error Recursive header files inclusion detected in TangentialCoverDecomposition.h
#else // defined(TangentialCoverDecomposition_RECURSES)
#define TangentialCoverDecomposition_RECURSES

#if !defined TangentialCoverDecomposition_h
#define TangentialCoverDecomposition_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <deque>
#include <queue>
#include <vector>
#include <set>
#include <map>
#include "ImaGene/base/CountedPtr.h"
#include "ImaGene/planes/NuThickDisk.h"
#include "ImaGene/planes/TangentialCoverHierarchy.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  
  /////////////////////////////////////////////////////////////////////////////
  // class TangentialCoverDecomposition
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'TangentialCoverDecomposition' <p> Aim:
   * Represents the decomposition into active vertices of the
   * tangential cover over an arbitrary digital surface graph.
   *
   * @tparam DigitalGraph the type of digital graph.
   *
   * @todo limited to 3D.
   */
  template <typename TDigitalGraph>
  class TangentialCoverDecomposition
  {
  public:
    typedef TDigitalGraph DigitalGraph;
    typedef typename DigitalGraph::Vertex Vertex;
    typedef typename DigitalGraph::Set Set;
    typedef typename DigitalGraph::Neighborhood Neighborhood;
    typedef typename Neighborhood::const_iterator NeighborhoodConstIterator;
    typedef typename Set::const_iterator SetConstIterator;
    typedef typename NuThickDisk<DigitalGraph>::Point3i Point3i;
    typedef Weighted2Cell<Vertex> WeightedVertex;
    typedef std::priority_queue< WeightedVertex > PriorityQueue;
    typedef NuThickDisk< DigitalGraph > ComputedDisk;
    typedef typename ComputedDisk::Expander Expander;
    typedef std::map< Vertex, ComputedDisk > MapVertexToDisk;
    typedef typename MapVertexToDisk::iterator MapVertexToDiskIterator;
    typedef std::set< Vertex > VertexSet;
    typedef typename VertexSet::iterator VertexSetIterator;
    typedef typename VertexSet::const_iterator VertexSetConstIterator;
    typedef std::map< Vertex, VertexSet > MapVertexToVertexSet;
    typedef TangentialCoverMethods<DigitalGraph> Methods;
    enum DecompositionType {
      Dominance, Segmentation, Accretion, MaximalPlanes, ExactMaximalPlanes 
    };

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~TangentialCoverDecomposition();

    /**
     * Constructor. 
     */
    TangentialCoverDecomposition();

    /**
     * Initialization from digital graph.
     *
     * @param dgraph the digital graph onto which the tangential cover
     * is extracted.
     *
     * @param scale the scale at which the hierarchy is computed.
     *
     * @param active_vertices the set of vertices that will be used
     * for computing the decomposition.
     */
    void init( const DigitalGraph & dgraph, 
	       double scale,
	       const Set & active_vertices );

    /**
       @return the current digital graph.
    */
    const DigitalGraph* digitalGraph() const;

    /**
       @return the scale of this decomposition.
    */
    double scale() const;

    /**
       Computes the decomposition of the tangential cover.
       @param methods useful methods for computing a tangential cover.
    */
    void compute( Methods & methods,
		  DecompositionType type );

    /**
       Computes the decomposition of the tangential cover. The method
       here is to define a dominance relation between neighboring
       vertices, then to propagate it. A vertex may be dominated by
       several vertices.

       @param methods useful methods for computing a tangential cover.
    */
    void computeDominance( Methods & methods );

    /**
       Computes the segmentation of the tangential cover. The method
       here is to define a dominance relation between neighboring
       vertices, then to propagate it. The vertex keeps its first
       dominant vertex. 

       @param methods useful methods for computing a tangential cover.
       @deprecated Do not use.
    */
    void computeSegmentation( Methods & methods );

    /**
       Computes the accretion decomposition of the tangential
       cover. The method here is to define the accretion of a vertex x
       as the set of vertices whose disk is included in the band
       defined by the disk of x.

       @param methods useful methods for computing a tangential cover.
    */
    void computeAccretion( Methods & methods );

    /**
       Computes the accretion decomposition of the tangential
       cover. The method here computes the set of maximal planes of
       the object (similar to accretion).

       @param methods useful methods for computing a tangential cover.
    */
    void computeMaximalPlanes( Methods & methods, 
			       bool exact = false );

    /**
     * Clears the decomposition as if it was just created.
     */
    void clear();

    /**
     * @param v any vertex.
     * @return the representatives of this vertex at this scale.
     */
    const VertexSet & representatives( const Vertex & v ) const;

    /**
     * @return all the representants at this scale.
     */
    const VertexSet & allRepresentatives() const;


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
    const DigitalGraph* myDGraph; /**< the digital graph. */
    double myScale;               /**< the scale for the decomposition. */
    Set* myActiveVertices;        /**< the set of active vertices. */ 
    VertexSet myAllRepresentatives; /**< the set of all representatives
				      vertices at this scale. */
    MapVertexToVertexSet myRepresentativesMap; /**< the map giving the
						 representatives for
						 each vertex. */
    MapVertexToDisk myDisks;       /**< Stores disks for each
				      representative. */
    // ------------------------- Hidden services ------------------------------
  protected:

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE TangentialCoverDecomposition( const TangentialCoverDecomposition & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE TangentialCoverDecomposition & operator=( const TangentialCoverDecomposition & other );
  
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'TangentialCoverDecomposition'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'TangentialCoverDecomposition' to write.
   * @return the output stream after the writing.
   */
  template <typename TDigitalGraph>
  std::ostream&
  operator<<( std::ostream & that_stream, 
	      const TangentialCoverDecomposition<TDigitalGraph> & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/planes/TangentialCoverDecomposition.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined TangentialCoverDecomposition_h

#undef TangentialCoverDecomposition_RECURSES
#endif // else defined(TangentialCoverDecomposition_RECURSES)
