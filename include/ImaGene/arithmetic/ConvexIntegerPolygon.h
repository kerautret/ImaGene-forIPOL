/** @file ConvexIntegerPolygon.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : ConvexIntegerPolygon.h
//
// Creation : 2011/01/11
//
// Version : 2011/01/11
//
// Author : JOL
//          Emilie Charrier (COBA algorithm)
//
// Summary : Header file for files ConvexIntegerPolygon.ih and ConvexIntegerPolygon.cxx
//
// History :
//	2011/01/11 : ?Name? : ?What?
//
// Rcs Id : "@(#)class ConvexIntegerPolygon declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(ConvexIntegerPolygon_RECURSES)
#error Recursive header files inclusion detected in ConvexIntegerPolygon.h
#else // defined(ConvexIntegerPolygon_RECURSES)
#define ConvexIntegerPolygon_RECURSES

#if !defined ConvexIntegerPolygon_h
#define ConvexIntegerPolygon_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <list>
#include <vector>
#include "ImaGene/arithmetic/IntegerComputer.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  
  /////////////////////////////////////////////////////////////////////////////
  // class ConvexIntegerPolygon
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'ConvexIntegerPolygon' <p>
   * Aim: Represents a convex polygon in the two-dimensional digital plane.
   */
  class ConvexIntegerPolygon : public IntegerComputer
  {

    // ----------------------- Standard types ------------------------------
  public:
    typedef std::list<Point2I> VertexList;
    typedef VertexList::iterator Iterator;
    typedef VertexList::const_iterator ConstIterator;


    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~ConvexIntegerPolygon();
    /**
     * Constructor.
     */
    ConvexIntegerPolygon();
    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    ConvexIntegerPolygon( const ConvexIntegerPolygon & other );
    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    ConvexIntegerPolygon & operator=( const ConvexIntegerPolygon & other );

    /**
     * @return the list of vertices of the polygon.
     */
    const VertexList & vertices() const;
    /**
     * @return the list of vertices of the polygon.
     */
    VertexList & vertices();
    /**
     * Reinitializes the polygon.
     */
    void clear();
    /**
     * Removes (duplicate) consecutive vertices.
     */
    void purge();

    /**
     * adds the point K to the convex polygon before position "pos".
     * @param pos any iterator
     * @param K the point to add
     * @return an iterator on the newly created element.
     */
    Iterator addBefore( const Iterator & pos, const Point2I & K );

    /**
     * adds the point K to the end of the polygon.
     * @param K the point to add
     */
    void addEnd( const Point2I & K );

    /**
     * @return 2*area of polygon.
     */
    const Integer & twiceArea();

    /**
     * if area is not 0, computes centroid, else, computes the middle
     * of the straight line segment.
     *
     * The centroid is a 2D rational point but it is represented as a
     * 3D integer point (a/d,c/d) corresponds to (a,b,d).
     *
     * @param point_test (modified) the centroid.
     *
     * The centroid is \b not in reduced form.
     */
    void centroid( Point3I & point_test);

    /**
     * cuts the convex polygon with the constraint N.(x,y) <= c
     *
     * @return 'true' if the polygon was modified, 'false' otherwise.
     */
    bool cutOpti( const Point2I & N, const Integer & c );

    // ----------------------- Helper methods ----------------------------------
    
    /**
     * initializes v and the two first vertices of resultup
     * and resultdown.
     *
     * A lying on C1, satisfying C2 and the closest to C2 B lying on
     * C1, not satisfying C2 and the closest to C2 v: valid Bezout
     * vector of the direction vector of C1 resultup and resultdown:
     * vertices of the grid convex hull satisfying and not satisfyin
     * C2 resp.
     *
     * @param resultup the first point of up intersection.
     * @param resultdown the first point of down intersection.
     *
     * @return 1 when the intersection of the supporting lines of C1
     * and C2 corresponds to an integer point, otherwise 0.
     */
    int init( const Point2I & N1, const Integer & c1, 
	      const Point2I & N2, const Integer & c2, 
	      Point2I & v,
	      std::vector<Point2I> & resultup, 
	      std::vector<Point2I> & resultdown );

    /**
     * computes the constraint of the form N.P<=c whose supporting
     * line passes through A and B such that the points refPoint1 and
     * refPoint2 satisfy the constraint
     */
    void computeConstraintFromPoints( Point2I & N, Integer & c, 
				      const Point2I & A, 
				      const Point2I & B, 
				      const Point2I & refPoint1, 
				      const Point2I & refPoint2 );

    /**
     * computes the border of the upper and of the lower convex hull
     * from the starting points resultup[0] (up) and resultdown[0]
     * down, along the constraint N2.p <= c2 while the vertices
     * satisfy the constraint N3.p <= c3. The vertices of the two
     * borders are stored at the end of resultup and resultdown.
     */
    void computeBorderPart1( const Point2I & initialBezoutVector, 
			     const Point2I & N2, const Integer & c2,
			     const Point2I & N3, const Integer & c3, 
			     std::vector<Point2I> & resultup,
			     std::vector<Point2I> & resultdown );

    /**
     * compute the convex hull of grid points satisfying the
     * constraints N1.P<=c1, N2.P<=c2 and N3.P>=c3.
     *
     * N2.P<=c2 corresponds to the cut two parts of computation: from
     * constraint 1 to constraint 3 and from constraint 3 to
     * constraint 1.
     *
     * The computed vertices are inserted at position [pos] in some list.
     *
     * @param pointRefC1 and pointRefC3 corresponds to grid point lying on
     * the supporting lines of C1 and of C3 resp.
     *
     * @param pos corresponds to an iterator in the list of vertices
     * of the convex, to add the next new vertices
     *
     * NB: the method also computes grid point satisfying N1.P<=c1 and
     * N3.P>=c3 but not satisfying N2.P<=c2. They are stored in
     * "resultdown" of size "nbverticesdown".  the algorithm uses
     * these points that's why they appear in the code.
     */
    void convexHullBorder( const Point2I & pointRefC1, 
			   const Point2I & pointRefC3,
			   const Point2I & N1, const Integer & c1,
			   const Point2I & N2, const Integer & c2, 
			   const Point2I & N3, const Integer & c3,
			   Iterator & pos );


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
  
    /**
     * Checks the class 'ConvexIntegerPolygon'.
     * @return 'true' if everything went well.
     */
    bool selfTest();

    // ------------------------- Datas ----------------------------------------
  private:
    VertexList myVertices;
    mutable Integer _area;
    mutable Integer _a, _b, _c, _den;
    mutable Point2I _p, _sum;
    mutable Point2I _directionVector;
    mutable Integer _fl, _ce;
    mutable Integer _dp1,_g;
    mutable Point2I _A, _B, _u, _v;
    mutable Point2I _A1,_B1,_A2,_B2;
    mutable std::vector<bool> _visible;
    mutable Point2I _N1, _N3;
    mutable Integer _c1, _c3;
    mutable std::vector<Point2I> _resultup;
    mutable std::vector<Point2I> _resultdown;
    // ------------------------- Hidden services ------------------------------
  protected:

  private:
 
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'ConvexIntegerPolygon'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'ConvexIntegerPolygon' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const ConvexIntegerPolygon & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/arithmetic/ConvexIntegerPolygon.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined ConvexIntegerPolygon_h

#undef ConvexIntegerPolygon_RECURSES
#endif // else defined(ConvexIntegerPolygon_RECURSES)
