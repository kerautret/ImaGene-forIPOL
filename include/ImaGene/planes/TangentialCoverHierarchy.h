/** @file TangentialCoverHierarchy.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : TangentialCoverHierarchy.h
//
// Creation : 2011/01/14
//
// Version : 2011/01/14
//
// Author : JOL
//
// Summary : Header file for files TangentialCoverHierarchy.ih and TangentialCoverHierarchy.cxx
//
// History :
//	2011/01/14 : ?Name? : ?What?
//
// Rcs Id : "@(#)class TangentialCoverHierarchy declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(TangentialCoverHierarchy_RECURSES)
#error Recursive header files inclusion detected in TangentialCoverHierarchy.h
#else // defined(TangentialCoverHierarchy_RECURSES)
#define TangentialCoverHierarchy_RECURSES

#if !defined TangentialCoverHierarchy_h
#define TangentialCoverHierarchy_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <vector>
#include <map>
#include "ImaGene/base/CountedPtr.h"
#include "ImaGene/planes/NuThickDisk.h"
#include "ImaGene/planes/TangentialCover.h"
#include "ImaGene/helper/ScaleProfile.h"
//////////////////////////////////////////////////////////////////////////////
#include "ImaGene/mathutils/Statistic.h"
//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{

  /////////////////////////////////////////////////////////////////////////////
  // class TangentialCoverHierarchy
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'TangentialCoverHierarchy' <p> Aim:
   * Represents the whole hierarchy of tangential cover over an
   * arbitrary digital surface graph.
   *
   * @tparam TDigitalGraph the type of digital graph.
   *
   * @todo limited to 3D.
   */
  template <typename TDigitalGraph>
  class TangentialCoverHierarchy
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
    typedef std::map< Vertex, ComputedDisk > MapVertexToDisk;
    typedef typename MapVertexToDisk::iterator MapVertexToDiskIterator;
    typedef TangentialCoverMethods<DigitalGraph> Methods;
    typedef std::map< Vertex, Vertex > MapVertexToVertex;
    typedef std::set< Vertex > VertexSet;
    typedef std::map< Vertex, VertexSet > MapVertexToVertexSet;
    typedef std::deque<WeightedVertex> WeightedVertexList;
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~TangentialCoverHierarchy();

    /**
     * Constructor. The object is not valid.
     */
    TangentialCoverHierarchy();

    /**
       Given the scale index, return its tangential cover.
       @param scale_idx the scale index (as given in init).
       @return the associated tangential cover.
    */
    const TangentialCover<DigitalGraph> & 
    tangentialCover( int scale_idx ) const;

    /**
       Given the scale index, return its tangential cover.
       @param scale_idx the scale index (as given in init).
       @return the associated tangential cover.
    */
    TangentialCover<DigitalGraph> & 
    tangentialCover( int scale_idx );

   /**
     * Initialization from digital graph.
     *
     * @param dgraph the digital graph onto which the tangential cover
     * is extracted.
     *
     * @param scales the scales at which the hierarchy is computed.
     */
    void init( const DigitalGraph & dgraph, 
	       const std::vector<double> & scales );

    /**
     * Clears the hierarchy as if it was just created.
     */
    void clear();

    /** 
     * Given the vertices sorted from greater radius to smaller
     * radius, computes the list of active vertices.
     *
     * Here, \f$ non-activity( v1 ) = \exists v2, ball_\nu( v1 ) subset
     * ball_\nu( v2 ) \quad \text{and} \quad radius_\nu( v1 ) <
     * radius_\nu( v2 ).\f$
     *
     * @param fp (modified) the priority queue of vertices (sorted from the one
     * with greater disk radius). Empty at the end.
     *
     * @param active_vtx (modified) the list of active vertices at
     * this scale.
     *
     * @param map_inactive2active (modified) the mapping that associates
     * to each inactive vertex some active vertices.
     */
    void computeActiveVertices( PriorityQueue & fp, 
				WeightedVertexList & active_vtx,
				MapVertexToVertexSet & map_inactive2active ) const;
    
    /** 
     * Given the vertices sorted from greater radius to smaller
     * radius, computes the list of active vertices.
     *
     * Here, \f$ non-activity( v1 ) = \exists v2, ball_\nu( v1 ) subset
     * ball_\nu( v2 ) \quad \text{and} \quad radius_\nu( v1 ) <
     * radius_\nu( v2 ).\f$
     *
     * @param fp (modified) the priority queue of vertices (sorted from the one
     * with greater disk radius). Empty at the end.
     *
     * @param active_vtx (modified) the list of active vertices at
     * this scale.
     */
    void computeActiveVertices( PriorityQueue & fp, 
				std::list<WeightedVertex> & active_vtx ) const;

    /** 
     * Given the vertices sorted from greater radius to smaller
     * radius, computes the list of active vertices.
     *
     * Here, \f$ non-activity( v1 ) = \exists v2, disk_\nu( v1 )
     * subset ext_\nu( v2 ) \quad \text{and} \quad radius_\nu( v1 ) <
     * radius_\nu( v2 ).\f$
     *
     * @param fp (modified) the priority queue of vertices (sorted from the one
     * with greater disk radius). Empty at the end.
     *
     * @param active_vtx (modified) the list of active vertices at
     * this scale.
     */
     void computeActiveVertices2
    ( PriorityQueue & fp, 
      std::list<WeightedVertex> & active_vtx,
      Methods & methods,
      double width ) const;

    /** 
     * Given the vertices sorted from greater radius to smaller
     * radius, computes the list of active vertices.
     *
     * Here, \f$ non-activity( v1 ) = \exists v2, v1 \in ball_\nu( v2 )
     * \quad \text{and} \quad disk_\nu( v1 ) <
     * ext_\nu( v2 ).\f$
     *
     * @param fp (modified) the priority queue of vertices (sorted from the one
     * with greater disk radius). Empty at the end.
     *
     * @param active_vtx (modified) the list of active vertices at
     * this scale.
     */
     void computeActiveVertices3
    ( PriorityQueue & fp, 
      std::list<WeightedVertex> & active_vtx,
      Methods & methods,
      double width ) const;

    /**
       Complete the new map inactive -> active vertices with the former
       map inactive -> actives (at the preceding scale).
       
       @param i2a (updated) the new map inactive -> actives, which has
       been computed at this scale, but which must be completed with
       the map at the preceding scale.

       @param former_i2a the map inactive -> actives at the preceding
       scale.
    */
    void completeActiveVertices( MapVertexToVertexSet & i2a,
				 const MapVertexToVertexSet & former_i2a ) const;
    VertexSet reduceActives( Vertex vtx, MapVertexToVertexSet & i2a );


    /**
       Computes the tangential covers at all scales. 
       
       @param methods useful methods for computing a tangential cover.
       
       @param not_in_core when 'true' computes all maximal planes
       whose center are not included in the core of a bigger maximal
       plane, otherwise retain only the maximal planes whose center
       are not included in a bigger maximal plane.
     */
    void computeTangentialCovers( Methods & methods, bool not_in_core );

    void compute( Methods & methods,
		  bool ball_inclusion = true );
    void compute2( Methods & methods );
    
    void computeNormals( Methods & methods, std::map< int,std::map< double,Statistic<float> > > & greedyNxMap,
                                                                                  std::map< int,std::map< double,Statistic<float> > > & greedyNyMap,
										  std::map< int,std::map< double,Statistic<float> > > & greedyNzMap);

    Set activeVertices( int scale_idx ) const;
    Set activeVertices( double scale ) const;

    /**
     * test if the ball centered on c1.v of radius c1.dist is included 
     * in the ball centered on C2.v of radius c2.dist
     *
     * @param c1 weighted cell representing a ball (voxel center,
     * radius of the ball)
     *
     * @param c2 weighted cell representing a ball (voxel center,
     * radius of the ball)
     *
     * @return true if the ball centered on c1 is included in the ball
     * centered on c2
     */
    bool isBallIncluded( const WeightedVertex & c1, 
			 const WeightedVertex & c2 ) const;

    /**
     * test if the vertex c1.cell is strictly included in the ball
     * centered on c2.cell of radius c2.weight2
     *
     * @param c1 weighted cell representing a ball (voxel center,
     * radius of the ball)
     *
     * @param c2 weighted cell representing a ball (voxel center,
     * radius of the ball)
     *
     * @return true if vertex c1 is included in the ball
     * centered on c2
     */
    bool isIncluded( const WeightedVertex & c1, 
		     const WeightedVertex & c2 ) const;

    // ----------------------- Noise services ---------------------------------
  public:

    /**
       @param sp (modified) is updated so as to contain the scale
       profile of [vtx].

       @param vtx any Vertex of the digital graph.

       @param mp_area use the maximal plane projected area, otherwise
       use the disk projected area.
    */
    void getScaleProfile( ScaleProfile & sp,
			  Vertex vtx,
			  bool mp_area, bool saveForMedian=false ) const;


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

    const DigitalGraph* myDGraph;  /**< the digital graph. */
    std::vector<double> myScales;  /**< the scales for the hierarchy. */
    std::vector<Set*> myActiveVtx; /**< the set of actives vertices at
				      each scale */
    /**
       The sequence of tangential covers.
    */
    std::vector< TangentialCover<DigitalGraph> > myTangentialCovers; 
    // ------------------------- Hidden services ------------------------------
  protected:

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE TangentialCoverHierarchy( const TangentialCoverHierarchy & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE TangentialCoverHierarchy & operator=( const TangentialCoverHierarchy & other );
  
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'TangentialCoverHierarchy'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'TangentialCoverHierarchy' to write.
   * @return the output stream after the writing.
   */
  template <typename TDigitalGraph>
  std::ostream&
  operator<<( std::ostream & that_stream, 
	      const TangentialCoverHierarchy<TDigitalGraph> & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/planes/TangentialCoverHierarchy.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined TangentialCoverHierarchy_h

#undef TangentialCoverHierarchy_RECURSES
#endif // else defined(TangentialCoverHierarchy_RECURSES)
