/** @file COBAPlaneRecognitionWithSet.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : COBAPlaneRecognitionWithSet.h
//
// Creation : 2011/01/12
//
// Version : 2011/01/12
//
// Author : JOL
//          Emilie Charrier
//
// Summary : Header file for files COBAPlaneRecognitionWithSet.ih and COBAPlaneRecognitionWithSet.cxx
//
// History :
//	2011/01/12 : ?Name? : ?What?
//
// Rcs Id : "@(#)class COBAPlaneRecognitionWithSet declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(COBAPlaneRecognitionWithSet_RECURSES)
#error Recursive header files inclusion detected in COBAPlaneRecognitionWithSet.h
#else // defined(COBAPlaneRecognitionWithSet_RECURSES)
#define COBAPlaneRecognitionWithSet_RECURSES

#if !defined COBAPlaneRecognitionWithSet_h
#define COBAPlaneRecognitionWithSet_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <set>
#include "ImaGene/base/CowPtr.h"
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
  // class COBAPlaneRecognitionWithSet
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'COBAPlaneRecognitionWithSet' <p> Aim: A class that
   * contains the COBA algorithm (Emilie Charrier, Lilian Buzer) for
   * recognizing pieces of digital planes of given vertical thickness.
   */
  class COBAPlaneRecognitionWithSet
  {
  public:

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
      inline bool operator<( const Point3i & p2 ) const
      {
	return ( coords[ 0 ] < p2.coords[ 0 ] )
	  || ( ( coords[ 0 ] == p2.coords[ 0 ] )
	       && ( ( coords[ 1 ] < p2.coords[ 1 ] )
		    || ( ( coords[ 1 ] == p2.coords[ 1 ] )
			 && ( coords[ 2 ] < p2.coords[ 2 ] ) ) ) );
      }
      
    };

    typedef std::set<Point3i>::iterator Iterator;
    typedef std::set<Point3i>::const_iterator ConstIterator;

    /**
     * The state is useful to come back to a previous step in the algorithm.
     */
    struct State {
      ConvexIntegerPolygon cip; /**< current constraint integer polygon. */
      Point3I centroid;         /**< current centroid of cip. */
      Point3I N;                /**< current normal vector. */
      Integer max;              /**< current max dot product. */
      Integer min;              /**< current min dot product. */
      mutable ConstIterator indmax;  /**< 3D point giving the max dot
					product. */
      mutable ConstIterator indmin;  /**< 3D point giving the min dot
					product. */
      Point2I grad;             /**< the current gradient. */
      /** the "accepted" set of points */ 
      mutable CowPtr< std::set<Point3i> > myCorePts; 
      mutable std::set<Point3i> myPts;  /**< the "current" ring of points. */
    };

    
    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~COBAPlaneRecognitionWithSet();

    /**
     * Constructor.
     * The object is not valid.
     */
    COBAPlaneRecognitionWithSet();

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    COBAPlaneRecognitionWithSet( const COBAPlaneRecognitionWithSet & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    COBAPlaneRecognitionWithSet & operator=( const COBAPlaneRecognitionWithSet & other );

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
       Chooses the "best" band given the current normal, set of points
       and current thickness. Should be called before check.
       
       @see check
    */
    void optimizeBand();

    /**
     * Checks if  the point [p] is in the current digital plane.
     *
     * @param p any 3D point (in the specified diameter).
     *
     * @return 'true' if it is in the plane, false otherwise.
     */
    bool check( const Point3i & p ) const;

    /**
     * Useful to come back to a previous state.
     *
     * @param state (update) gets the current state of the algorithm.
     */
    void getState( State & state );

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

    /**
     * @tparam Vector3D any type T such that T.operator[](int i)
     * returns a reference to a double. i ranges in 0,1,2.
     *
     * @param (updates) the current upper point.
     */
    template <typename Vector3D>
    void getUpper( Vector3D & upper ) const;

    /**
     * @tparam Vector3D any type T such that T.operator[](int i)
     * returns a reference to a double. i ranges in 0,1,2.
     *
     * @param (updates) the current lower point.
     */
    template <typename Vector3D>
    void getLower( Vector3D & lower ) const;

    /**
     * @tparam Vector3D any type T such that T.operator[](int i)
     * returns a reference to a double. i ranges in 0,1,2.
     *
     * @param (updates) the current unit normal vector 
     */
    template <typename Vector3D>
    void getUnitNormal( Vector3D & normal ) const;

    /**
     * If n is the unit normal to the current plane, then n.x >= min
     * and n.x <= max are the two half-planes defining it.
     *
     * @param min the lower bound (corresponding to the unit vector).
     * @param max the upper bound (corresponding to the unit vector).
     */
    void getBounds( double & min, double & max ) const;
    
    /**
      * Useful when a solution vector has been found (is DP)
      * Look for the solution vector which minimizes the thickness
      *
      * cuts the search space and tests the centroid
      * returns the associated width or -1 if the cip is not bigger than on single point
    */
    double refineNormal( void );
    
    /** 
    * computes and returns the thickness of the set of points
    */
    double thickness();

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
    mutable Integer _v, _maxValue;
    mutable Integer _cst1, _cst2, tmp;
    mutable Point3I _p;
    // ------------------------- Hidden services ------------------------------
  private:

  
  public:
    /**
     * Updates myState.centroid and myState.N.
     */
    void computeDigitalNormalVector();

    // ------------------------- Internals ------------------------------------
  private:

    /**
     * Checks if the given point of myPts belongs to the current
     * plane. Updates myState accordingly.
     *
     * @param p any 3D point (in the specified diameter).
     *
     * @return 'true' if this is the case, otherwise returns 'false' and
     * computes the gradient for finding a new valid plane.
     */
    bool oracleForPoint( Iterator it );

    /**
     * Checks if the given point of myPts belongs to the current
     * plane. Updates myState accordingly.
     *
     * @param p any 3D point (in the specified diameter).
     *
     * @return 'true' if this is the case, otherwise returns 'false' and
     * computes the gradient for finding a new valid plane.
     */
    bool oracleForPointNoNewNormal( Iterator it );

    /**
     * Checks if all the points of myPts belongs to the current
     * plane. Updates myState accordingly.
     *
     * @return 'true' if this is the case, otherwise returns 'false' and
     * computes the gradient for finding a new valid plane.
     *
     * computeDigitalNormalVector should have been called beforehand.
     */
    bool oracleForPoints();
    
    /**
     * Computes the gradient to find a "better" plane (i.e. a thiner one)
     * need to determine the points defining the min and the max dot product
     *
     * @return 'true' if a not null gradient has been computed, 
     * otherwise returns 'false' (it means that the current normal vector is optimal)
     *
     */
    bool computeGradient();

    /**
     * Performs the double cut in parameter space according to the current
     * gradient and width.
     */
    void doubleCut();
    
    /**
     * Performs the simple cut in parameter space according to the current
     * gradient.
     */
    void simpleCut();
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'COBAPlaneRecognitionWithSet'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'COBAPlaneRecognitionWithSet' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const COBAPlaneRecognitionWithSet & that_object_to_display );

  inline
  std::ostream &
  operator<<( std::ostream & out, const COBAPlaneRecognitionWithSet::Point3i & p )
  {
    out << "(" << p[ 0 ] << "," << p[ 1 ] << "," << p[ 2 ] << ")";
    return out;
  }
  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/arithmetic/COBAPlaneRecognitionWithSet.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined COBAPlaneRecognitionWithSet_h

#undef COBAPlaneRecognitionWithSet_RECURSES
#endif // else defined(COBAPlaneRecognitionWithSet_RECURSES)
