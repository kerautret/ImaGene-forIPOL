///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : IntegerComputer.cxx
//
// Creation : 2011/01/10
//
// Version : 2011/01/10
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
//	2011/01/10 : Mr ?Name? : ?What?
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //


///////////////////////////////////////////////////////////////////////////////
#include "ImaGene/arithmetic/IntegerComputer.h"
// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/arithmetic/IntegerComputer.ih"
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const IntegerComputer_RCS_ID = "@(#)class IntegerComputer definition.";



///////////////////////////////////////////////////////////////////////////////
// class IntegerComputer
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::IntegerComputer::~IntegerComputer()
{
}

/**
 * Constructor.  Each thread must have its own instance for all
 * computations. Such object stores several local variables to
 * limit the number of memory allocations.
 */
ImaGene::IntegerComputer::IntegerComputer()
{
}

/**
 * returns a solution of the Diophantine equation: a x + b y = c
 */
ImaGene::Point2I 
ImaGene::IntegerComputer::ExtendedEuclid
( const Integer & a, const Integer & b, const Integer & c ) const
{
  int k;
  for ( k = 0; k < 4; ++k )
    _eu_tab_bezout[ k ].clear();

  if( a == I_ZERO )  return Point2I( I_ZERO, b * c );
  if( b == I_ZERO )  return Point2I( a * c, I_ZERO );

  _eu_a = _abs( a );
  _eu_b = _abs( b );
  
  _eu_tab_bezout[ 0 ].push_back( _eu_a );
  _eu_tab_bezout[ 0 ].push_back( _eu_b );
  _eu_tab_bezout[ 1 ].push_back( I_ZERO );
  _eu_tab_bezout[ 2 ].push_back( I_ONE );
  _eu_tab_bezout[ 2 ].push_back( I_ZERO );
  _eu_tab_bezout[ 3 ].push_back( I_ZERO );
  _eu_tab_bezout[ 3 ].push_back( I_ONE );
  
  k = 0;
  while( _eu_tab_bezout[ 0 ][ k+1 ] != I_ZERO )
    {
      _eu_tab_bezout[ 1 ].push_back( _eu_tab_bezout[ 0 ][ k ] 
				     / _eu_tab_bezout[ 0 ][ k+1 ] );
      _eu_tab_bezout[ 0 ].push_back( _eu_tab_bezout[ 0 ][ k ] 
				     % _eu_tab_bezout[ 0 ][ k+1 ] );
      _eu_tab_bezout[ 2 ].push_back( _eu_tab_bezout[ 2 ][ k ] 
				     - _eu_tab_bezout[ 1 ][ k+1 ]
				     * _eu_tab_bezout[ 2 ][ k+1 ] );
      _eu_tab_bezout[ 3 ].push_back( _eu_tab_bezout[ 3 ][ k ] 
				     - _eu_tab_bezout[ 1 ][ k+1 ]
				     *_eu_tab_bezout[ 3 ][ k+1 ] );
      k++;
    }

  _eu_v.x = _eu_tab_bezout[ 2 ][ k ];
  _eu_v.y = _eu_tab_bezout[ 3 ][ k ];

  if( a < I_ZERO )
    _eu_v.x = -_eu_v.x;
  if( b < I_ZERO )
    _eu_v.y = -_eu_v.y;
  _eu_v = _eu_v * c;
  return _eu_v;
}


///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param that_stream the output stream where the object is written.
 */
void 
ImaGene::IntegerComputer::selfDisplay( ostream& that_stream ) const
{
  that_stream << "[IntegerComputer]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::IntegerComputer::OK() const
{
  return true;
}


/**
 * Checks the classes Integer, Point2I, Point3I, IntegerComputer.
 * @return 'true' if everything went well.
 */
bool
ImaGene::IntegerComputer::selfTest()
{
  unsigned int nb_ok = 0;
  unsigned int nb = 0;
  Integer i( 12 );
  Integer j( 15 );
  Point2I p1( i, j );
  Point2I p2( 12, 15 );
  nb_ok += p1 == p2 ? 1 : 0;
  nb++;
  Point2I z( 0, 0 );
  nb_ok += ( p1 - p2 == z ) ? 1 : 0;
  nb++;

  nb_ok += selfTestExtendedEuclid( 13324243, 74328428, 4321 );
  nb++;
  nb_ok += selfTestExtendedEuclid( 8, 5, 1 );
  nb++;
  nb_ok += selfTestExtendedEuclid( 8, 5, -1 );
  nb++;

  vector<Integer> q;
  nb_ok += cfrac( q, 8, 5 ) == I_ONE;
  nb++;
  cerr << "8/5=[" << q[ 0 ] << "," << q[ 1 ] << "," << q[ 2 ] << "," 
       << q[ 3 ] << "]" << endl;
  cerr << "[IntegerComputer::selfTest] " << nb_ok << "/" << nb << " passed."
       << endl;
  return nb_ok == nb;
}

bool
ImaGene::IntegerComputer::selfTestExtendedEuclid
( Integer a, Integer b, Integer c )
{
  Point2I p = ExtendedEuclid( a, b, c );
  Integer c_check = a * p.x + b * p.y;
  cerr << a << " * " << p.x << " + " << b << " * " << p.y << " = " 
       << c_check << "( == " << c << " )" << endl;
  return c == c_check;
}

///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
