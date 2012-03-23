/** @file IntegerComputer.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : IntegerComputer.h
//
// Creation : 2011/01/10
//
// Version : 2011/01/10
//
// Author : JOL
//          EC (COBA algorithm)
//
// Summary : Header file for files IntegerComputer.ih and IntegerComputer.cxx
//
// History :
//	2011/01/10 : ?Name? : ?What?
//
// Rcs Id : "@(#)class IntegerComputer declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(IntegerComputer_RECURSES)
#error Recursive header files inclusion detected in IntegerComputer.h
#else // defined(IntegerComputer_RECURSES)
#define IntegerComputer_RECURSES

#if !defined IntegerComputer_h
#define IntegerComputer_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  /////////////////////////////////////////////////////////////////////////////
  // class Integer
  /////////////////////////////////////////////////////////////////////////////

#ifdef BIG_INTEGER_IS_LONGLONG
  typedef long long Integer;
  inline int get_si( const Integer & i )
  {
    return (int) i;
  }
  
  inline double get_d( const Integer & i )
  {
    return (double) i;
  }
  const Integer I_ZERO = 0LL;
  const Integer I_ONE = 1LL;
  const Integer I_TWO = 2LL;
  const Integer I_THREE = 3LL;
  const Integer I_MINUS_ONE = -1LL;
#else
  typedef mpz_class Integer;
  inline int get_si( const Integer & i )
  {
    return i.get_si();
  }
  inline double get_d( const Integer & i )
  {
    return i.get_d();
  }
  const Integer I_ZERO = Integer( 0 );
  const Integer I_ONE = Integer( 1 );
  const Integer I_TWO = Integer( 2 );
  const Integer I_THREE = Integer( 3 );
  const Integer I_MINUS_ONE = Integer( -1 );
#endif

  Integer _abs(const Integer & a );
  Integer _max(const Integer & a,const Integer & b);
  Integer _min(const Integer & a,const Integer & b);
  //compute the floor value of na/nb
  Integer _floordiv(const Integer &na, const Integer &nb);
  //compute the ceiling value of na/nb
  Integer _ceildiv(const Integer &na, const Integer &nb);
  //compute gcd of a and b
  Integer _gcd(const Integer & a,const Integer & b);
  // std::ostream & operator<<( std::ostream & out, const Integer & a );

  /////////////////////////////////////////////////////////////////////////////
  // class Point2I
  /////////////////////////////////////////////////////////////////////////////
  /**
   * Represents a point in the two-dimensional digital plane.
   */
  class Point2I
  {
  public :
    
    Integer x;
    Integer y;

    Point2I();
    Point2I& operator = (const Point2I & A);
    Point2I(const Point2I & A);
    Point2I(const Integer & a, const Integer & b);
    
    void reduce();  // tq pgcd(composantes) = 1
    void rnd(int lg_max);    // random function
    Point2I & operator+=( const Point2I & A );
    Point2I & operator-=( const Point2I & A );
    Point2I & operator*=( const Integer & B );
    Point2I & operator/=( const Integer & B );
    void neg();
  };

  Point2I operator + (const Point2I & A , const Point2I & B);
  Point2I operator - (const Point2I & A , const Point2I & B);
  Point2I operator * (const Point2I & A , const Integer & B);
  Point2I operator * (const Integer & B , const Point2I & A);
  Integer operator ^ (const Point2I & A , const Point2I & B);
  Integer operator * (const Point2I & A , const Point2I & B);
  
  bool operator == (const Point2I & A, const Point2I & B);
  bool operator != (const Point2I & A, const Point2I & B);

  std::ostream & operator<<( std::ostream & out, const Point2I & a );

  /**
   * Represents a point in the three-dimensional digital plane.
   */
  class Point3I
  {
  public :
    Integer x;
    Integer y;
    Integer z;

    Point3I();
    Point3I( const Point3I & t );
    Point3I& operator=( const Point3I & A );
    Point3I( const Integer & _x, const Integer & _y, const Integer & _z );
    
    void reduce();  // tq pgcd(composantes) = 1
    void neg();
  };
  
  Point3I operator + (const Point3I & a, const Point3I & b);
  Point3I operator - (const Point3I & a, const Point3I & b);
  Point3I operator ^ (const Point3I & a, const Point3I & b);
  Integer operator * (const Point3I & a, const Point3I & b);
  Point3I operator * (const Integer & a, const Point3I & b);
  Point3I operator * (const Point3I & b, const Integer & a);
  Point3I operator / (const Point3I & b, const Integer & a);
  bool operator == (const Point3I & A, const Point3I & B);
  bool operator != (const Point3I & A, const Point3I & B);

  std::ostream & operator<<( std::ostream & out, const Point3I & a );

  /////////////////////////////////////////////////////////////////////////////
  // class IntegerComputer
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'IntegerComputer' <p> Aim: This class
   * gathers several types and methods to make computation with
   * integers. To be thread-safe, each thread should instantiate an
   * IntegerComputer.
   */
  class IntegerComputer
  {

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~IntegerComputer();

    /**
     * Constructor.  Each thread must have its own instance for all
     * computations. Such object stores several local variables to
     * limit the number of memory allocations.
     */
    IntegerComputer();

    // ----------------------- Integer services ------------------------------
  public:
    /**
     * computes the floor value of na/nb.
     */
    Integer floordiv( const Integer & na, const Integer & nb ) const;
    /**
     * computes the ceil value of na/nb.
     */
    Integer ceildiv( const Integer & na, const Integer & nb ) const;

    /**
     * computes the gcd of a and b (a and b may be either positive or negative).
     */
    Integer gcd( const Integer & a, const Integer & b ) const;
    /**
     * computes the gcd of a and b (a and b may be either positive or negative).
     */
    void gcd( Integer & g, const Integer & a, const Integer & b ) const;

    /**
     * computes and push_backs the simple continued fraction of a
     * (positive) and b (positive).
     * @return the gcd of a and b
     */
    Integer cfrac( std::vector<Integer> & quotients,
		const Integer & a, const Integer & b ) const;

    // ----------------------- Point2I services ------------------------------
  public:
    /**
     * makes p irreducible.
     */
    void reduce( Point2I & p ) const ;
    /**
     * computes the cross product of u and v.
     */
    Integer cross_product( const Point2I & u, const Point2I & v) const;
    /**
     * computes the cross product of u and v.
     */
    void cross_product( Integer & cp, const Point2I & u, const Point2I & v) const;
    /**
     * computes the dot product of u and v.
     */
    Integer dot_product( const Point2I & u, const Point2I & v ) const;
    /**
     * computes the dot product of u and v.
     */
    void dot_product( Integer & dp, const Point2I & u, const Point2I & v) const;
     /**
     * computes the dot product of u and v.
     */
    void dot_product( Integer & dp, const Point3I & u, const Point3I & v) const;
    /**
     * returns a solution of the Diophantine equation: a x + b y = c
     */
    Point2I ExtendedEuclid( const Integer & a, const Integer & b, 
			    const Integer & c ) const;
    /**
     * computes the floor (fl) and the ceiling (ce) value of the real
     * number k such that p + k u lies on the supporting line of the
     * constraint N.p <= c
     */
    void coefficientIntersection( const Point2I & p, 
				  const Point2I & u, const Point2I & N, 
				  const Integer & c,
				  Integer & fl, Integer & ce) const;
    /**
     * compute the valid bezout vector v of u such that B+v satifies
     * the constraints C and C2 and such that B+v+u doesn't satify the
     * constraint C2 if redef==0, v is already initialized and the
     * constraint N.p<=c is not used in the computation
     */
    void validBezout ( const Point2I & A, const Point2I & u,
		       Point2I & v, const Point2I & N, 
		       const Integer & c, const Point2I & N2, 
		       const Integer & c2, int redef ) const;

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
     * Checks the classes Integer, Point2I, Point3I, IntegerComputer.
     * @return 'true' if everything went well.
     */
    bool selfTest();

    bool selfTestExtendedEuclid( Integer a, Integer b, Integer c );

    // ------------------------- Datas ----------------------------------------
  private:
    mutable Integer _fda, _fdb;
    mutable Integer _cda, _cdb;
    mutable Integer _gaa, _gbb, _ga0, _ga1, _gr;
    mutable Integer _rt;
    mutable Integer _cp1,_cp2;
    mutable Integer _dp1,_dp2;
    mutable Integer _eu_a, _eu_b;
    mutable std::vector<Integer> _eu_tab_bezout[ 4 ];
    mutable Point2I _eu_v;
    mutable Integer _ci1, _ci2, _ci3;
    mutable Point2I _vbp1, _vbp2;
    mutable Integer _vbi1, _vbfl, _vbce;
    // ------------------------- Hidden services ------------------------------
  protected:


  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE IntegerComputer( const IntegerComputer & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE IntegerComputer & operator=( const IntegerComputer & other );
  
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'IntegerComputer'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'IntegerComputer' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const IntegerComputer & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/arithmetic/IntegerComputer.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined IntegerComputer_h

#undef IntegerComputer_RECURSES
#endif // else defined(IntegerComputer_RECURSES)
