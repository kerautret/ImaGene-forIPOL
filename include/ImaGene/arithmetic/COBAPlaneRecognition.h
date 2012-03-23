/** @file COBAPlaneRecognition.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : COBAPlaneRecognition.h
//
// Creation : 2011/01/12
//
// Version : 2011/01/12
//
// Author : JOL
//          Emilie Charrier
//
// Summary : Header file for files COBAPlaneRecognition.ih and COBAPlaneRecognition.cxx
//
// History :
//	2011/01/12 : ?Name? : ?What?
//
// Rcs Id : "@(#)class COBAPlaneRecognition declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(COBAPlaneRecognition_RECURSES)
#error Recursive header files inclusion detected in COBAPlaneRecognition.h
#else // defined(COBAPlaneRecognition_RECURSES)
#define COBAPlaneRecognition_RECURSES

#if !defined COBAPlaneRecognition_h
#define COBAPlaneRecognition_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "ImaGene/arithmetic/IntegerComputer.h"
#include "ImaGene/arithmetic/ConvexIntegerPolygon.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  
  /////////////////////////////////////////////////////////////////////////////
  // class COBAPlaneRecognition
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'COBAPlaneRecognition' <p> Aim: A class that
   * contains the COBA algorithm (Emilie Charrier, Lilian Buzer) for
   * recognizing pieces of digital planes of given vertical thickness.
   */
  class COBAPlaneRecognition
  {
  public:
    /**
     * The state is useful to come back to a previous step in the algorithm.
     */
    struct State {
      ConvexIntegerPolygon cip; /**< current constraint integer polygon. */
      Point3I centroid;         /**< current centroid of cip. */
      Point3I N;                /**< current normal vector. */
      Integer max;              /**< current max dot product. */
      Integer min;              /**< current min dot product. */
      int indmax;               /**< 3D point giving the max dot product. */
      int indmin;               /**< 3D point giving the min dot product. */
      Point2I grad;             /**< the current gradient. */
      int nb_used;              /**< number of meaningful 3D points in myPts. */
    };

    struct Point3i {
      int coords[ 3 ];
      inline Point3i() {}
      inline Point3i( const int * t )
      {
	for ( int i = 0; i < 3; ++i ) coords[ i ] = t[ i ];
      }
      inline Point3i( int x, int y, int z )
      {
	coords[ 0 ] = x;	
	coords[ 1 ] = y;
	coords[ 2 ] = z;
      }
      inline Point3i( const Point3i & other )
      {
	for ( int i = 0; i < 3; ++i ) 
	  coords[ i ] = other.coords[ i ];
      }
      inline Point3i & operator=( const Point3i & other )
      {
	if ( this != &other )
	  for ( int i = 0; i < 3; ++i ) 
	    coords[ i ] = other.coords[ i ];
	return *this;
      }
      inline const int* data() const
      {
	return coords;
      }
      inline int* data()
      {
	return coords;
      }
      inline int & operator[ ]( int i )
      {
	return coords[ i ];
      }
      inline int operator[ ]( int i ) const
      {
	return coords[ i ];
      }
      inline bool operator==( const Point3i & p2 ) const
      {
	return ( coords[ 0 ] == p2.coords[ 0 ] )
	  && ( coords[ 1 ] == p2.coords[ 1 ] )
	  && ( coords[ 2 ] == p2.coords[ 2 ] );
      }
      
    };

    
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~COBAPlaneRecognition();

    /**
     * Constructor.
     * The object is not valid.
     */
    COBAPlaneRecognition();

    /**
     * Clear the object, free memory. It is no more valid.
     */
    void clear();

    /**
     * All these parameters cannot be changed during the process.
     * After this call, the object is in a consistent state and can
     * accept new points for recognition.
     *
     * @param axis the main axis (0,1,2) for x, y or z.
     *
     * @param diameter the diameter for the set of points (maximum
     * distance between the given points)
     *
     * @param width the maximal axis-width (x,y,or z) for the plane. 
     *
     * @param pt the first point for initializing the plane
     * recognition algorithm.
     */
    void init( int axis, int diameter, double width, const Point3i & pt );

    /**
     * @return the i-th point stored in this object.
     */
    const Point3i & point( int i ) const;

    /**
     * @return the number of points in the current state.
     */
    unsigned int size() const;

    /**
     * Adds the point [p] and checks if we have still a digital plane
     * of specified width.
     *
     * @param p any 3D point (in the specified diameter).
     *
     * @param allowNewNormal if 'true' the normal may be updated,
     * 'false' the normal is never updated.
     *
     * @return 'true' if it is still a plane, false otherwise.
     */
    bool add( const Point3i & p, bool allowNewNormal = true );

    /**
     * Useful to come back to a previous state.
     *
     * @param state (update) gets the current state of the algorithm.
     */
    void getState( State & state ) const;

    /**
     * Useful to come back to a previous state.
     *
     * @param state the new current state for the algorithm (previous
     * state is dropped).
     */
    void setState( const State & state );

    /**
     * @tparam Vector3D any type T such that T.operator[](int i)
     * returns a reference to a double. i ranges in 0,1,2.
     *
     * @param (updates) the current normal vector 
     */
    template <typename Vector3D>
    void getNormal( Vector3D & normal ) const;

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
    int myAxis;     /**< the main axis used in all subsequent computations. */
    Integer myG;    /**< the grid step used in all subsequent computations. */
    double myWidth; /**< the plane width used in all subsequent computations. */
    State myState;  /**< the current step. */
    /**
     * the current points in the digital plane. Only myState.nb_used
     * are meaningful. */
    std::vector<Point3i> myPts; 
    Integer _v, _maxValue;
    Integer _cst1, _cst2;
    Point3I _p;
    // ------------------------- Hidden services ------------------------------
  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE COBAPlaneRecognition( const COBAPlaneRecognition & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE COBAPlaneRecognition & operator=( const COBAPlaneRecognition & other );
  
    // ------------------------- Internals ------------------------------------
  private:
    /**
     * Updates myState.centroid and myState.N.
     */
    void computeDigitalNormalVector();

    /**
     * Checks if the given point of myPts belongs to the current
     * plane. Updates myState accordingly.
     *
     * @param idx any valid index of a point of myPts.
     *
     * @return 'true' if this is the case, otherwise returns 'false' and
     * computes the gradient for finding a new valid plane.
     */
    bool oracleForPoint( int idx );

    /**
     * Checks if all the points of myPts between idx_begin (included) and
     * idx_end (excluded) belongs to the current plane. Updates myState
     * accordingly.
     *
     * @param idx_begin any valid index of a point of myPts, generally 0.
     * @param idx_end any (excluded) valid index of a point of myPts,
     * generally nb_used.
     *
     * @return 'true' if this is the case, otherwise returns 'false' and
     * computes the gradient for finding a new valid plane.
     *
     * computeDigitalNormalVector should have been called beforehand.
     */
    bool oracleForPoints( int idx_begin, int idx_end );

    /**
     * Performs the double cut in parameter space according to the current
     * gradient and width.
     */
    void doubleCut();
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'COBAPlaneRecognition'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'COBAPlaneRecognition' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const COBAPlaneRecognition & that_object_to_display );

  inline
  std::ostream &
  operator<<( std::ostream & out, const COBAPlaneRecognition::Point3i & p )
  {
    out << "(" << p[ 0 ] << "," << p[ 1 ] << "," << p[ 2 ] << ")";
    return out;
  }
  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/arithmetic/COBAPlaneRecognition.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined COBAPlaneRecognition_h

#undef COBAPlaneRecognition_RECURSES
#endif // else defined(COBAPlaneRecognition_RECURSES)
