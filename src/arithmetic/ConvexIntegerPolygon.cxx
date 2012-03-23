///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : ConvexIntegerPolygon.cxx
//
// Creation : 2011/01/11
//
// Version : 2011/01/11
//
// Author : Jacques-Olivier Lachaud
//          Emilie Charrier (COBA algorithm)
//
// email : lachaud@labri.fr
//
// Purpose : ??
//
// Distribution :
//
// Use :
//	??
//
// Todo :
//	O ??
//
// History :
//	2011/01/11 : Mr ?Name? : ?What?
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //


///////////////////////////////////////////////////////////////////////////////
#include "ImaGene/arithmetic/ConvexIntegerPolygon.h"
// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/arithmetic/ConvexIntegerPolygon.ih"
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const ConvexIntegerPolygon_RCS_ID = "@(#)class ConvexIntegerPolygon definition.";



///////////////////////////////////////////////////////////////////////////////
// class ConvexIntegerPolygon
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::ConvexIntegerPolygon::~ConvexIntegerPolygon()
{
}

/**
 * @return 2*area of polygon.
 */
const ImaGene::Integer & 
ImaGene::ConvexIntegerPolygon::twiceArea()
{
  _area = I_ZERO;
  ConstIterator it = myVertices.begin();
  ConstIterator it_end = myVertices.end();
  ConstIterator it_next = it; ++it_next;
  while ( it_next != it_end )
    {
      _area += cross_product( *it, *it_next );
      it = it_next; ++it_next;
    }
  _area += cross_product( *it, *(myVertices.begin()) );
  return _area;
  
}

/**
 * if area is not 0, computes centroid, else, computes the middle
 * of the straight line segment.
 *
 * The centroid is a 2D rational point but it is represented as a
 * 3D integer point (a/d,b/d) corresponds to (a,b,d).
 *
 * @param point_test (modified) the centroid.
 */
void
ImaGene::ConvexIntegerPolygon::centroid( Point3I & point_test)
{
  _a = I_ZERO;
  _b = I_ZERO;
  _p = Point2I( I_ZERO, I_ZERO );
  twiceArea(); // updates _area;

  if( _area > I_ZERO )
    {
      _den = I_THREE * _area;
      ConstIterator it_begin = myVertices.begin();
      ConstIterator it = it_begin;
      ConstIterator it_end = myVertices.end();
      ConstIterator it_next = it; ++it_next;
      while ( it_next != it_end )
	{
	  cross_product( _c, *it, *it_next );
	  _sum = *it + *it_next;
	  _sum *= _c;
	  _p += _sum;
	  it = it_next; ++it_next;
	}
      cross_product( _c, *it, *it_begin );
      _sum = *it + *it_begin;
      _sum *= _c;
      _p += _sum;
    }
  else
    {
      _den = I_ZERO;
      ConstIterator it = myVertices.begin();
      ConstIterator it_end = myVertices.end();
      for ( ; it != it_end; ++it )
	{
	  _p += *it;
	  ++_den;
	}
    }
  point_test.x = _p.x;
  point_test.y = _p.y;
  point_test.z = _den;
}


/**
 * cuts the convex polygon with the constraint N.(x,y) <= c
 *
 * @return 'true' if the polygon was modified, 'false' otherwise.
 */
bool 
ImaGene::ConvexIntegerPolygon::cutOpti
( const Point2I & N, const Integer & c )
{
  //cout<<"begin cut opti"<<endl;
  int n = myVertices.size(); //number of vertices
  _visible.resize( n ); //table of visibility of each vertex
  twiceArea(); // in _area: twice the area of the current convex polygon
  int index;

  // _N1, _c1, _N3, _c3 : to determine the two constraints supported
  // by vertices of the polygon

  // computation of the visibility of each vertex
  int i = 0;
  int nb_visibles = 0;
  bool is_visible;
  for ( ConstIterator it = myVertices.begin();
	it != myVertices.end();
	++it )
    {
      dot_product( _c, *it, N );
      is_visible = _c <= c;
      _visible[ i++ ] = is_visible;
      if ( is_visible ) ++nb_visibles;
    }

  // cerr << "nb=" << n << " nb_visibles=" << nb_visibles << endl;
  if ( nb_visibles == n ) { return false; }
  if ( nb_visibles == 0 ) { myVertices.clear(); return true; }

  // determines the 2 edges A1B1 and A2B2 intersected by the constaint N.P <= c
  // A1 and A2 are visible (they satisfy the constraint)
  bool A1B1Found = false;
  bool A2B2Found = false;

  i = 0;
  Iterator it_begin = myVertices.begin();
  Iterator it = it_begin;
  Iterator it_next = it; ++it_next;
  Iterator it_A1;
  while( it != myVertices.end() )
    {
      if ( ( _visible[i] ) && ( ! _visible[(i+1)%n] ) )  
	{ 
	  _A1 = *it;
	  _B1 = ( it_next == myVertices.end() ) ? *it_begin : *it_next;
	  it_A1 = it;
	  A1B1Found=true; 
	}
      else if ( ( ! _visible[i] ) && ( _visible[(i+1)%n] ) )
	{
	  _A2 = ( it_next == myVertices.end() ) ? *it_begin : *it_next;
	  _B2 = *it;
	  A2B2Found=true; 
	}
      it = it_next; ++it_next; 
      ++i;
    }

  // delete non visible vertices
  i = 0;
  it = myVertices.begin();
  while( it != myVertices.end() )
    {
      if ( ! _visible[i] )
	{
// 	  cerr << "[ removing " << *it << " ]" << endl;
	  it = myVertices.erase( it );
	}
      else
	++it;
      ++i;
    }

  // no need to determine the index of the vertex A1
  if ( *it_A1 != _A1 )
    cerr << "[ImaGene::ConvexIntegerPolygon::cutOpti] ERROR 1 !" << endl;
  else if ( ( ! A1B1Found ) || ( ! A2B2Found ) )
    cerr << "[ImaGene::ConvexIntegerPolygon::cutOpti] ERROR 2 !" << endl;

  Iterator it_after_A1 = it_A1; ++it_after_A1;

  if ( _area > I_ZERO ) //convex not reduced to a straight line segment
    {
      // cerr << "Non null area." << endl;
      //computes the constraint C whose supporting line passes through A1 and B1
      computeConstraintFromPoints( _N1, _c1, _A1, _B1, _A2, _B2);
      //computes the constraint C whose supporting line passes through A2 and B2
      computeConstraintFromPoints( _N3, _c3, _A2, _B2, _A1, _B1);

      //run the reconstruction of the convex polygon
      convexHullBorder( _A1, _A2, _N1, _c1, N, c, _N3, _c3, it_after_A1 );
      //printVertices();
    }
  else //convex reduced to a straight line segment
    {
      //compute the new extremity of the straight line segment
      _p = _B1 - _A1;
      reduce( _p );
      _a = ( c - N * _A1 ) / ( N * _p );
//       if( _a * N * _p == c - N * _A1 )
// 	{--_a; cout<<"moins moins a "<<_a<<endl;}
      _A1 += _a * _p;
      addBefore( it_after_A1, _A1 );
    }
  purge();
  return true;
}


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
void
ImaGene::ConvexIntegerPolygon::convexHullBorder
( const Point2I & pointRefC1, 
  const Point2I & pointRefC3,
  const Point2I & N1, const Integer & c1,
  const Point2I & N2, const Integer & c2, 
  const Point2I & N3, const Integer & c3,
  Iterator & pos )
{
  // _u, _v: vectors u and v to determine the next vertex
  int integerIntersection;
  // initializes A, B, u, v and the two first vertices of resultup and
  // resultdown.
  _resultup.resize( 1 );    //to store half convex hull border
  _resultdown.resize( 1 );  //to store half convex hull border
  _resultup[ 0 ] = pointRefC1;

  // integerIntersection is equal to one when
  // the intersection of the supporting lines of C1
  // and C2 corresponds to an integer point
  integerIntersection = init( N1, c1, N2, c2, 
			      _v, _resultup, _resultdown );
  if( integerIntersection != 1 ) // not integer intersection
  {
    //computation of the first part of the border
    computeBorderPart1( _v, N2, c2, N3, c3, 
			_resultup, _resultdown );
  }
  for( int i = 0; i < _resultup.size(); ++i) //fill in convexup
    {
      addBefore( pos, _resultup[ i ] );
    }  

  // second part
  // initializes A, B, u, v and the two first vertices of resultup and
  // resultdown.
  _resultup.resize( 1 );
  _resultdown.resize( 1 );
  _resultup[0] = pointRefC3;
  // integerIntersection is equal to one when the intersection of the
  // supporting lines of C3 and C2 corresponds to an integer point
  integerIntersection = init( N3, c3, N2, c2, 
			      _v, _resultup, _resultdown );

  if ( integerIntersection != 1 ) //not integer intersection
  {
    //computation of the second part of the border
    computeBorderPart1( _v, N2, c2, N1, c1, 
			_resultup, _resultdown );
  }
  
  for( int i = _resultup.size() - 1; i >= 0; --i ) //fill in convexup
    {
      addBefore( pos, _resultup[ i ] ); 
    }
}


int
ImaGene::ConvexIntegerPolygon::init
( const Point2I & N1, const Integer & c1, 
  const Point2I & N2, const Integer & c2, 
  Point2I & v, 
  vector<Point2I> & resultup, vector<Point2I> & resultdown )
{
  //cout<<"begin init"<<endl;
  int result;

  //initialize  vector directionVector (not definitive)
  _directionVector.y = N1.x;
  _directionVector.x = - N1.y;

  //compute the intersection of ray (resultup[0],directionVector) with constraint C2
  // cout<<"before coef intersect"<<endl;
  coefficientIntersection( resultup[0], _directionVector, N2, c2,
			   _fl, _ce );
  //cout<<"after coef intersect"<<endl;

  // uses the intersection to compute the first vertex of the upper convex hull
  // i.e. the grid point closest to C2 and satisfying C2
  dot_product( _dp1, resultup[0], N2 );
  if ( ( _dp1 > c2 ) && ( _fl == I_ZERO ) )
  {
    resultup[0] += _directionVector * _ce;
    _directionVector.neg();
  }
  else if ( ( ( _dp1 <= c2 ) && ( _fl >= I_ZERO ) )
	    || ( ( _dp1 > c2 ) && ( _fl <= I_ZERO ) ) )
    {
      resultup[0] += _directionVector * _fl;
    }
  else
    {
      resultup[0] += _directionVector * _ce;
      _directionVector.neg();
    }

  //compute the first vertex of the lower convex hull
  if ( _fl == _ce )  //integer intersection
    {
      resultdown[0] = resultup[0];
      result = 1;
    }
  else
    {
      resultdown[0] = resultup[0] + _directionVector;
      //initialization of v: valid Bezout vector of u
      validBezout( resultup[0], _directionVector, v, N1, c1, N2, c2, 1 );
      result = 0;
    }
  //cout<<"end init"<<endl;
  return result;
}


/**
 * computes the constraint of the form N.P<=c whose supporting
 * line passes through A and B such that the points refPoint1 and
 * refPoint2 satisfy the constraint
 */
void
ImaGene::ConvexIntegerPolygon::computeConstraintFromPoints
( Point2I & N, Integer & c, 
  const Point2I & A, const Point2I & B, 
  const Point2I & refPoint1, const Point2I & refPoint2 )
{
  //cout<<" begin compute constraint from points"<<endl;
  N.x = A.y - B.y;
  N.y = B.x - A.x;
  dot_product( c, N, A );
  dot_product( _a, N, refPoint1 );
  dot_product( _b, N, refPoint2 );
  
  if ( ( _a > c ) || ( _b > c ) )
    {
      N.neg();
      c = -c;
    }
  
  //cout<<"N=("<<N.x<<","<<N.y<<")"<<endl;
  //simplification of the constraint
  gcd( _g, N.x, N.y );
  //cout<<"gcd="<<GCD<<endl;
  N /= _g;
  c = floordiv( c, _g);
  //cout<<"compute constraint from points c="<<c<<endl;
}

/**
 * computes the border of the upper and of the lower convex hull
 * from the starting points resultup[0] (up) and resultdown[0]
 * down, along the constraint N2.p <= c2 while the vertices
 * satisfy the constraint N3.p <= c3. The vertices of the two
 * borders are stored at the end of resultup and resultdown.
 */
void
ImaGene::ConvexIntegerPolygon::computeBorderPart1
( const Point2I & BV, 
  const Point2I & N2, const Integer & c2,
  const Point2I & N3, const Integer & c3, 
  std::vector<Point2I> & resultup,
  std::vector<Point2I> & resultdown )
{
  // _A and _B represents the last computed vertex above and below the constraint resp.
  // _u, _v pair of Bezout vectors.
  _v = BV;
  //initialize A and B
  _A = resultup[0];
  _B = resultdown[0];

  // while A and B do not lie on the supporting line of C2
  // and satisfies C3 and while v is not parallel to C2
  while ( ( dot_product( _A, N2 ) != c2 ) 
	  && ( dot_product( _A, N3 ) <= c3 )
	  && (dot_product( _B, N2 ) != c2 ) 
	  && ( dot_product( _B, N3 ) <= c3 )
	  && ( dot_product( _v, N2 ) != I_ZERO ) )
    {
      if ( dot_product( _v, N2 ) < I_ZERO ) //second configuration, we find a new B
	{
	  //computation of the new vertex
	  coefficientIntersection( _B, _v, N2, c2, _fl, _ce );
	  _B += _v * _fl;
	  resultdown.push_back( _B );
	}
      else //first configuration, we find a new A
	{ 
	  //computation of the new vertex
	  coefficientIntersection( _A, _v, N2, c2, _fl, _ce );
	  _A += _v * _fl;
	  resultup.push_back( _A );
	}
      // update u and v
      _u = _B;
      _u -= _A;
      validBezout( _A, _u, _v, N2, c2, N2, c2, 0);
    }
  //when the loop finishes, we have to complete the computation
  //of each convex hull
  if ( dot_product( _A, N3 ) > c3 ) 
    // A does not satisfy C3, we remove it.
    resultup.pop_back();
  else if ( dot_product( _B, N3 ) > c3 ) 
    // B does not satisfy C3, we remove it
    resultdown.pop_back();
  else if( dot_product( _A, N2 ) == c2 ) 
    //A lies on C2 so A is also a vertex of resultdown
    resultdown.push_back( _A );
  else if( dot_product( _B, N2) == c2 )
    //B lies on C2 so B is also a vertex of resultup
    resultup.push_back( _B );
}


///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param that_stream the output stream where the object is written.
 */
void 
ImaGene::ConvexIntegerPolygon::selfDisplay( ostream& that_stream ) const
{
  that_stream << "[";
  for ( ConstIterator it = vertices().begin();
	it != vertices().end();
	++it )
    that_stream << " " << *it;
  that_stream << " ]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::ConvexIntegerPolygon::OK() const
{
  return true;
}

/**
 * Checks the class 'ConvexIntegerPolygon'.
 * @return 'true' if everything went well.
 */
bool 
ImaGene::ConvexIntegerPolygon::selfTest()
{
  unsigned int nb_ok = 0;
  unsigned int nb = 0;
  Integer t( 100 );

  addEnd( Point2I( t, t ) );
  addEnd( Point2I( -t, t ) );
  addEnd( Point2I( -t, -t ) );
  addEnd( Point2I( t, -t ) );
  
  Point3I c;
  centroid( c );
  c.reduce();
  Point2I c2( c.x, c.y );
  cerr << "2*Area = " << twiceArea()
       << " Centroid = " << c 
       << endl;
  nb_ok += twiceArea() == Integer( 80000 ) ? 1 : 0;
  nb++;
  nb_ok += c2 == Point2I( 0, 0 ) ? 1 : 0;
  nb++;

  addEnd( Point2I( 2*t, 0 ) );
  centroid( c );
  c.reduce();
  cerr << "2*Area = " << twiceArea()
       << " Centroid = " << c 
       << endl;
  nb_ok += twiceArea() == Integer( 100000 ) ? 1 : 0;
  nb++;

  cerr << *this << endl;
  cutOpti( Point2I( 1, 1 ), Integer( 50 ) );
  cerr << "[ (-100,100) (-100,-100) (100,-100) (125,-75) (-50,100) ] (should be)" << endl;
  cerr << *this << endl;
  cutOpti( Point2I( 1, 1 ), Integer( 0 ) );
  cerr << "[ (-100,-100) (100,-100) (-100,100) ] (should be)" << endl;
  cerr << *this << endl;
  cutOpti( Point2I( 5, 3 ), Integer( 7 ) );
  cerr << *this << endl;
  cutOpti( Point2I( -2, 5 ), Integer( 12 ) );
 cerr << "[ (-100,-100) (61,-100) (59,-96) (5,-6) (3,-3) (-1,1) (-3,1) (-6,0) (-96,-36) (-100,-38) ] (should be)" << endl;
  cerr << *this << endl;

  cerr << "[ConvexIntegerPolygon::selfTest] " 
       << nb_ok << "/" << nb << " passed."
       << endl;
  return nb_ok == nb;
}




///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
