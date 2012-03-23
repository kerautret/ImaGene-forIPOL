/** @file TangentialCover.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : TangentialCover.h
//
// Creation : 2011/02/23
//
// Version : 2011/02/23
//
// Author : JOL
//
// Summary : Header file for files TangentialCover.ih and TangentialCover.cxx
//
// History :
//	2011/02/23 : ?Name? : ?What?
//
// Rcs Id : "@(#)class TangentialCover declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(TangentialCover_RECURSES)
#error Recursive header files inclusion detected in TangentialCover.h
#else // defined(TangentialCover_RECURSES)
#define TangentialCover_RECURSES

#if !defined TangentialCover_h
#define TangentialCover_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <map>
#include "ImaGene/planes/NuThickDisk.h"
#include "ImaGene/mathutils/Statistic.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  /**
   An interface specifying a set of useful methods for computing a
   tangential cover.

   @tparam TDigitalGraph the type of digital graph.
  */
  template <typename TDigitalGraph>
  struct TangentialCoverMethods {
    typedef TDigitalGraph DigitalGraph;
    typedef typename DigitalGraph::Vertex Vertex;
    virtual unsigned int preferredAxis( const Vertex & v ) = 0;
    virtual void adjustOrientation
    ( const Vertex & v, std::vector<double> & normal ) = 0;
    virtual double projectedArea( const Vertex & v, const double* normal ) = 0;
  };
  
  /////////////////////////////////////////////////////////////////////////////
  // class TangentialCover
  /////////////////////////////////////////////////////////////////////////////
  /** 
   @brief Description of class 'TangentialCover': Represents the set
   of maximal planes at a given scale.
   
   Generally computed with TangentialCoverHierarchy and
   TangentialCoverDecomposition.

   @tparam TDigitalGraph the type of digital graph.
   @todo limited to 3D.
   */
  template <typename TDigitalGraph>
  class TangentialCover
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
    typedef std::deque<WeightedVertex> WeightedVertexList;
    typedef std::priority_queue< WeightedVertex > PriorityQueue;
    typedef NuThickDisk< DigitalGraph > ComputedDisk;
    typedef std::map< Vertex, ComputedDisk > MapVertexToDisk;
    typedef typename MapVertexToDisk::iterator MapVertexToDiskIterator;
    typedef std::map< Vertex, Vertex > MapVertexToVertex;
    typedef std::set< Vertex > VertexSet;
    typedef std::map< Vertex, VertexSet > MapVertexToVertexSet;
    typedef TangentialCoverMethods<DigitalGraph> Methods;

    struct MaximalPlaneSummary {
      Vertex center;         /**< the center of the disk which defines
				the maximal plane. */
      double center_embedding[ 3 ]; /**< the position of the vertex. */
      unsigned int major_axis; /**< the major axis of this plane. */
      double radius;         /**< the maximal radius of the disk. */
      double limit_radius;   /**< the radius at which it is no more a
				nu-thick disk ( > radius). */
      double normal[ 3 ];    /**< the unit normal direction. */
      double upper;          /** upper plane offset: normal.x = upper */
      double lower;          /** lower plane offset: normal.x = lower */
      double projectedDiskArea; /**< the area of the projected disk. */
      double projectedArea;  /**< the area of the projected maximal plane. */
      void init( Methods & methods,
		 const ComputedDisk & disk, 
		 const Set & maximal_plane );
      template <typename Vector3D>
      void projectUpperPlane( Vector3D & proj, const Vector3D & p ) const;
      template <typename Vector3D>
      void projectLowerPlane( Vector3D & proj, const Vector3D & p ) const;
      template <typename Vector3D>
      void getLikelyUnitNormal( Vector3D & n ) const;
      void selfDisplay( std::ostream & out ) const;
    };
    
    typedef unsigned int MaximalPlaneIndex;
    typedef std::vector<MaximalPlaneSummary> MaximalPlaneList;
    typedef std::set<MaximalPlaneIndex> MaximalPlaneIndexSet;
    typedef std::map<Vertex,MaximalPlaneIndexSet> MapVertexToMaximalPlanes;
    typedef std::map<Vertex,MaximalPlaneIndex> MapVertexToMaximalPlane;

    enum AveragingMode { SimpleAveraging, 
			 DistanceAveraging, 
			 RadiusAndDistanceAveraging,
			 InOutAveraging,
			 MaxProjectedDisk,
			 MaxProjectedPlane };

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~TangentialCover();

    /**
     * Constructor.  Note: a tangential cover is usually filled by a
     * TangentialCoverHierarchy object.
     */
    TangentialCover();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    TangentialCover( const TangentialCover & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    TangentialCover & operator=( const TangentialCover & other );

    /**
       @return the current digital graph.
    */
    const DigitalGraph* digitalGraph() const;

    /**
       Computes the maximal planes. Note that myDGraph, myScale and
       myActiveVertices,myVtx2Active should have been initialized before.
     
       @param methods useful methods for computing a tangential cover.
       
       @param not_in_core when 'true' computes all maximal planes
       whose center are not included in the core of a bigger maximal
       plane, otherwise retain only the maximal planes whose center
       are not included in a bigger maximal plane.
     */
    void computeMaximalPlanes( Methods & methods,
			       bool not_in_core );
			       
     /**
       Computes the maximal planes. Note that myDGraph, myScale and
       myActiveVertices,myVtx2Active should have been initialized before.
     
       @param methods useful methods for computing a tangential cover.
       
       @param not_in_core when 'true' computes all maximal planes
       whose center are not included in the core of a bigger maximal
       plane, otherwise retain only the maximal planes whose center
       are not included in a bigger maximal plane.
       
       @param map_vtx2radius (modified) a map which associates each
       vertex of the surface to the radius of its centered disk. May
       be updated.
     */
    void computeMaximalPlanes( Methods & methods,
			       bool not_in_core, 
			       std::map<Vertex,double> & map_vtx2radius );

    /**
       Given a vertex, writes its maximal plane indices into the given
       collection with the specified output_iterator.
       
       @tparam OutputIterator any output_iterator (*it++ = ... is valid).
       @param vtx the vertex.
       
       @param OutputIterator the iterator in the modified collection
       that will hold the list of maximal planes of the vertex.
    */
    template <typename OutputIterator>
    void getMaximalPlanesIndices( Vertex vtx, OutputIterator it ) const;

    /**
       @param vtx the vertex.

       @return the index of its core maximal plane.
    */
    MaximalPlaneIndex getCoreMaximalPlaneIndex( Vertex vtx ) const;

    /**
       @param idx a valid index of maximal plane.
       @return a const-reference on the associated maximal plane summary.
    */
    const MaximalPlaneSummary & maximalPlane( MaximalPlaneIndex idx ) const;


    /**
       Computes the averaging coefficients for the given vertex [p] according to the choosen averaging mode.

       @param coefs (modified) contains the averaging coefficients (sum is 1.0)
       @param p any vertex.
       @param nd the mode chosen for averaging.
    */
    void getAveragingCoefficients( std::vector<double> & coefs, Vertex p, 
				   AveragingMode nd = SimpleAveraging ) const;

    /**
       Combines the normals of all maximal planes containing the vertex. 
       @tparam Vector3D any kind of double[3].
       @param n (modified) the normal estimated at [p].
       @param p any vertex.
       @param coefs the averaging coefficients.
       @see getAveragingCoefficients
    */
    template <typename Vector3D>
    void getEstimatedNormal( Vector3D & n, Vertex p,
			     const std::vector<double> & coefs ) const;

    /**
       Combines the normals of all maximal planes containing the vertex. 
       @tparam Vector3D any kind of double[3].
       @param n (modified) the normal estimated at [p].
       @param p any vertex.
       @param nd the mode chosen for computing normals.
    */
    template <typename Vector3D>
    void getEstimatedNormal( Vector3D & n, Vertex p,
			     AveragingMode nd = SimpleAveraging ) const;

    /**
       Computes the statistic of the angle of normals at [p] wrt
       vector [est_normal].

       @tparam Vector3D any kind of double[3].
       @param angle_stat (returns) the statistic object of these angles.
       @param p the vertex of interest
       @param est_normal the vector to compare with, generally obtained with a getEstimatedNormal
       
       @see getEstimatedNormal
    */
    template <typename Vector3D>
    void
    getNormalAngleStatistic( Statistic<double> & angle_stat,
			     Vertex p,
			     const Vector3D & est_normal ) const;

    /**
       Planes are associated to vertices, but their geometry are
       computed from the embedding of these vertices. Several vertices
       may thus have the same embedding. This method reassociates
       planes to vertices so that (i) two vertices having the same
       embedding have the same planes, (ii) only planes common to all
       vertices with same embedding survive.
    */
    void intersectPlanesAccordingToEmbedding();

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
  public:

    const DigitalGraph* myDGraph;  /**< the digital graph. */
    double myScale;               /**< The scale of this tangential cover. */
    WeightedVertexList myActiveVertices; /**< Active vertices ordered
					    by radius. */
    MaximalPlaneList myAllPlanes; /**< All the maximal planes. */
    MapVertexToMaximalPlanes myMapVtx2MP; /**< The mapping giving the
			     planes for each pertinent vertex. */
    MapVertexToVertexSet myVtx2Actives; /**< The mapping associating to each vertex aat least one active vertex. */
    MapVertexToMaximalPlane myMapVtx2CoreMP; /**< The mapping giving the core maximal plane of each vertex. */

    // ------------------------- Hidden services ------------------------------
  protected:

  private:


  
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'TangentialCover'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'TangentialCover' to write.
   * @return the output stream after the writing.
   */
  template <typename TDigitalGraph>
  std::ostream&
  operator<<( std::ostream & that_stream, 
	      const TangentialCover<TDigitalGraph> & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/planes/TangentialCover.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined TangentialCover_h

#undef TangentialCover_RECURSES
#endif // else defined(TangentialCover_RECURSES)
