/** @file NuThickDisk.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : NuThickDisk.h
//
// Creation : 2011/01/13
//
// Version : 2011/01/13
//
// Author : JOL
//
// Summary : Header file for files NuThickDisk.ih and NuThickDisk.cxx
//
// History :
//	2011/01/13 : ?Name? : ?What?
//
// Rcs Id : "@(#)class NuThickDisk declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(NuThickDisk_RECURSES)
#error Recursive header files inclusion detected in NuThickDisk.h
#else // defined(NuThickDisk_RECURSES)
#define NuThickDisk_RECURSES

#if !defined NuThickDisk_h
#define NuThickDisk_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <map>
#include "ImaGene/base/CowPtr.h"
#include "ImaGene/planes/CDigitalGraph.h"
#include "ImaGene/planes/IsotropicExpanderWithSet.h"
#include "ImaGene/arithmetic/COBAPlaneRecognition.h"
#include "ImaGene/arithmetic/COBAPlaneRecognitionWithSet.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{

  /////////////////////////////////////////////////////////////////////////////
  // class NuThickDisk
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'NuThickDisk' <p>
   * Aim: Represents a nu-thick disk around some vertex p of digital graph.
   * @tparam DigitalGraph the type of digital graph.
   *
   * @todo limited to 3D.
   */
  template <typename TDigitalGraph>
  class NuThickDisk
  {
  public:
    typedef TDigitalGraph DigitalGraph;
    typedef typename DigitalGraph::Vertex Vertex;
    typedef typename DigitalGraph::Set Set;
    typedef typename DigitalGraph::Neighborhood Neighborhood;
    typedef IsotropicExpanderWithSet<DigitalGraph,EuclideanDigitalDistance> Expander;
    typedef typename Expander::MarkedSet MarkedSet;
    typedef typename Neighborhood::const_iterator NeighborhoodConstIterator;
    typedef COBAPlaneRecognitionWithSet COBAAlgorithm;
    typedef COBAAlgorithm::Point3i Point3i;
    typedef WeightedCell<Vertex> WeightedVertex;
    typedef std::priority_queue< WeightedVertex > PriorityQueue;
    
    struct Computer {
      COBAAlgorithm coba;    /**< the plane recognition algorithms. */
      double radius;         /**< the maximal radius for each axis. */
      double limit_radius;   /**< the radius at which it is no more a nu-thick disk ( > radius). */
      Expander expander;
      bool finished;
    };
    // ----------------------- Standard services ------------------------------
  public:
    static double distance( const Point3i & p1, const Point3i & p2 );

    /**
     * Destructor. 
     */
    ~NuThickDisk();

    /**
     * Constructor.
     * The object is invalid.
     */
    NuThickDisk();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    NuThickDisk( const NuThickDisk & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    NuThickDisk & operator=( const NuThickDisk & other );

    /**
     * Initialization from digital graph and vertex.
     *
     * @param dgraph the digital graph onto which the nu-thick disk
     * are extracted.
     *
     * @param p the center vertex.
     *
     * @param thickness the specified thickness (at least 1.0 for
     * naive planes).
     */
    void init( const DigitalGraph & dgraph, const Vertex & p,
	        double thickness );

    /**
     * Initialization from other disk.
     *
     * @param p the center vertex.
     *
     * @param disk an other disk onto the same digital graph. The
     * nu-thick disk of this vertex should contain the other disk.
     */
    void init( const Vertex & p, 
	       const NuThickDisk<DigitalGraph> & disk );

    /**
     * Computes this nu-thick disk for the given thickness. Calls
     * computeForAxis for each possible axis.
     *
     * @param radius_upper_bound if negative, the maximal disk is
     * computed, otherwise stops at the first vertex whose distance is
     * greater than this value.
     *
     * @return 'true' if the complete nu-thick disk was extracted,
     * 'false' otherwise.
     */
    bool computeDisk( double radius_upper_bound = -1.0 );
		      
    /**
     *	Determine the thiner plane including the NuThickdisk
     */
    void refineNormal();

    /**
     * Computes the nu-thick extension from the current nu-thick disk
     * with the current thickness. Call optimizeBand.
     *
     * @param vertices (modified) the set of vertices belonging to  the
     * nu-expansion.
     */
    void computeExtension( Set & vertices );

    /**
       Computes the nu-thick maximal plane from the current nu-thick
       disk with the current thickness. The maximal plane is the
       connected set of vertices of the extension whose disk is also
       in this extension. Call optimizeBand.
       
       @param maximal_plane (modified) the union of disks included in
       the extension of this disk.

       @param core (modified) the set of vertices that is the core of
       this maximal plane, i.e. the centers of the disks included in
       the extension of this disk.
      
       @see optimizeBand
     */
    void computeMaximalPlane( Set & maximal_plane, Set & core );

    /**
       Computes the nu-thick maximal plane from the current nu-thick
       disk with the current thickness. The maximal plane is the
       connected set of vertices of the extension whose disk is also
       in this extension. Call optimizeBand.
       
       @param maximal_plane (modified) the union of disks included in
       the extension of this disk.

       @param core (modified) the set of vertices that is the core of
       this maximal plane, i.e. the centers of the disks included in
       the extension of this disk.

       @param all_cores the set of vertices that was the core of
       preceding maximal planes, which may not been in other cores.
      
       @see optimizeBand
     */
    template <typename VertexSet>
    void computeMaximalPlane( Set & maximal_plane, 
			      Set & core, const VertexSet & all_cores );
			      
    /**
       Computes the nu-thick maximal plane from the current nu-thick
       disk with the current thickness. The maximal plane is the
       connected set of vertices of the extension whose disk is also
       in this extension. Call optimizeBand.
       
       @param maximal_plane (modified) the union of disks included in
       the extension of this disk.

       @param core (modified) the set of vertices that is the core of
       this maximal plane, i.e. the centers of the disks included in
       the extension of this disk.

       @param all_cores the set of vertices that was the core of
       preceding maximal planes, which may not been in other cores.

       @param map_vtx2radius (modified) a map which associates each
       vertex of the surface to the radius of its centered disk. May
       be updated.
       
       @see optimizeBand
     */
    template <typename VertexSet>
    void computeMaximalPlane( Set & maximal_plane, 
			      Set & core, const VertexSet & all_cores, 
			      std::map< Vertex, double > & map_vtx2radius );

    /**
       Computes the nu-thick maximal plane from the current nu-thick
       disk with the current thickness. The maximal plane is the
       connected set of vertices of the extension whose disk is also
       in this extension. Call optimizeBand.
       
       @param maximal_plane (modified) the set of vertices belonging to  the
       nu-thick maximal plane.
       
       @see optimizeBand
       @deprecated use computeMaximalPlane(Set &, Set &)
     */
    void computeMaximalPlane( Set & maximal_plane );

    /**
       Computes the nu-thick maximal plane from the current nu-thick
       disk with the current thickness. The exact maximal plane is the
       connected set of vertices of the extension such that their disk
       union this disk is still a nu-thick plane.
       
       @param maximal_plane (modified) the set of vertices belonging to  the
       nu-thick maximal plane.
     */
    void computeExactMaximalPlane( Set & maximal_plane );

    /**
     * @param disk_vertices (returns) the set of vertices of the
     * nu-thick disk. Has meaning only if computeDisk was called till
     * completion.
     */
    void getDiskVertices( Set & disk_vertices ) const; 

    /**
     * @return the set of vertices of the nu-thick disk. Has meaning
     * only if computeDisk was called till completion.
     *
     * NB: MarkedSet is a shortname for std::set<Vertex>
     */
    MarkedSet diskVertices() const; 

    /**
     * @return the center vertex.
     */
    const Vertex & center() const;

    /**
     * @return the digital embedding of the center vertex.
     */
    const Point3i & centerEmbedding() const;

    /**
     * @return the major axis of the nu-thick disk (when computed).
     */
    unsigned int majorAxis() const;

    /**
     * @return the number of different vertices embedding in the
     * nu-thick disk (when computed).
     */
    unsigned int size() const;

    /**
     * @return the radius of the nu-thick disk (when computed).
     */
    double radius() const;

    /**
     * @return the limit radius of the nu-thick disk (when computed).
     */
    double limitRadius() const;

    /**
     * @return nu, ie. the thickness of the nu-thick disk.
     */
    double thickness() const;

    /**
     * @param n (modified) the normal vector to the nu-thick disk.
     */
    template <typename Vector3d>
    void getNormal( Vector3d & n ) const;

    /**
     * @param n (modified) the unit normal vector to the nu-thick disk.
     */
    template <typename Vector3d>
    void getUnitNormal( Vector3d & n ) const;

    /**
       Chooses the "best" band given the current normal, set of points
       and current thickness. Should be called before check.
       
       @see check
    */
    void optimizeBand();

    /**
       @param v any vertex.
       
       @return 'true' if v belongs to the band (intersection of two
       half-planes) spanned by the disk.
    */
    bool isInBand( const Vertex & v ) const;

    /**
       @param it any iterator.
       @param it_end any iterator.
       
       @return 'true' if all vertices pointed by the iterators belongs
       to the band (intersection of two half-planes) spanned by the
       disk.
    */
    template <typename VertexIterator>
    bool areInBand( VertexIterator it, VertexIterator it_end ) const;

    /**
     * Computes this nu-thick disk for the given axis.
     *
     * @param axis the chosen axis.
     *
     * @param radius_upper_bound if negative, the maximal disk is
     * computed, otherwise stops at the first vertex whose distance is
     * greater than this value.
     *
     * @return 'true' if the complete nu-thick disk was extracted,
     * 'false' otherwise.
     */
    bool computeForAxis( unsigned int axis,
			 double radius_upper_bound );

    /**
     * @param upper_pt (modified) a point above the disk center that
     * is lying on the upper plane.
     */
    template <typename Vector3d>
    void getUpperPoint( Vector3d & upper_pt ) const;

    /**
     * @param lower_pt (modified) a point below the disk center that
     * is lying on the lower plane.
     */
    template <typename Vector3d>
    void getLowerPoint( Vector3d & lower_pt ) const;

    /**
     * If n is the unit normal to the current plane, then n.x >= min
     * and n.x <= max are the two half-planes defining it.
     *
     * @param min the lower bound (corresponding to the unit vector).
     * @param max the upper bound (corresponding to the unit vector).
     */
    void getBounds( double & min, double & max ) const;

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
    const DigitalGraph* myDGraph;    /**< the digital graph. */
    EuclideanDigitalDistance myDist; /**< the digital distance. */
    Vertex myVtx; /**< Starting vertex. */
    Point3i myPt; /**< embedding of myVtx. */
    CowPtr<Computer> myComputer[ 3 ]; /**< One computer per axis. */
    double myThickness;
    int myMajorAxis;
    // ------------------------- Hidden services ------------------------------
  protected:

  private:

  
    // ------------------------- Internals ------------------------------------
  private:
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'NuThickDisk'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'NuThickDisk' to write.
   * @return the output stream after the writing.
   */
  template <typename TDigitalGraph>
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const NuThickDisk<TDigitalGraph> & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/planes/NuThickDisk.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined NuThickDisk_h

#undef NuThickDisk_RECURSES
#endif // else defined(NuThickDisk_RECURSES)
