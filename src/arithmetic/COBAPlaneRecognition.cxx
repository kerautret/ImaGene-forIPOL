///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : COBAPlaneRecognition.cxx
//
// Creation : 2011/01/12
//
// Version : 2011/01/12
//
// Author : Jacques-Olivier Lachaud
//          Emilie Charrier
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
//	2011/01/12 : Mr ?Name? : ?What?
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //


///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "ImaGene/arithmetic/COBAPlaneRecognition.h"
// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/arithmetic/COBAPlaneRecognition.ih"
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const COBAPlaneRecognition_RCS_ID = "@(#)class COBAPlaneRecognition definition.";



///////////////////////////////////////////////////////////////////////////////
// class COBAPlaneRecognition
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::COBAPlaneRecognition::~COBAPlaneRecognition()
{
}

/**
 * Constructor.
 * The object is not valid.
 */
ImaGene::COBAPlaneRecognition::COBAPlaneRecognition()
{
}

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
void 
ImaGene::COBAPlaneRecognition::init
( int axis, int diameter, double width, const Point3i & pt )
{
  myAxis = axis;
  myWidth = width;
  // initialize the grid step.
  myG = 2*diameter;
  myG *= diameter;
  myG *= diameter;
  // initialize the search space as a square.
  myState.cip.clear();
  myState.cip.addEnd( Point2I( -myG, -myG ) ); 
  myState.cip.addEnd( Point2I(  myG, -myG ) ); 
  myState.cip.addEnd( Point2I(  myG,  myG ) ); 
  myState.cip.addEnd( Point2I( -myG,  myG ) ); 
  computeDigitalNormalVector();
  //   myState.cip.centroid( myState.centroid );
  //   myState.centroid.reduce();
  if ( myPts.size() < 1 ) myPts.resize( 10 );
  myPts[ 0 ] = pt;
  myState.nb_used = 1;
  _p.x = pt.coords[ 0 ];
  _p.y = pt.coords[ 1 ];
  _p.z = pt.coords[ 2 ];
  myState.cip.dot_product( _v, myState.N, _p );
  myState.min = myState.max = _v; //myState.N * pt;
  myState.indmin = myState.indmax = 0;
  myState.grad = Point2I( 0, 0 );
}


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
bool
ImaGene::COBAPlaneRecognition::add
( const Point3i & p, bool allowNewNormal )
{
  // Add point to myPts.
  int idx = myState.nb_used;
  if ( myState.nb_used < myPts.size() )
    myPts[ myState.nb_used ] = p;
  else
    myPts.push_back( p );
  ++myState.nb_used;

  // Checks if new point is still in the plane.
  if ( oracleForPoint( idx ) )
    // yes, nothing to do.
    return true;

  // Checks if we can change the normal so as to find another digital plane.
  if( ( ( myState.grad.x == I_ZERO )
	&& ( myState.grad.y == I_ZERO ) )
      || ! allowNewNormal )
    {
      // last point is removed.
      --myState.nb_used;
      return false;
    }

  // cerr << "  -- opt " << myState.cip << endl;
  // If yes, tries to optimize.
  doubleCut();


  //while at least 1 point left on the search space
  while ( myState.cip.vertices().size() >= 1 )
  {
    computeDigitalNormalVector();
    // cerr << "  -- opt " << myState.cip << endl;
    //calls oracle
    if ( oracleForPoints( 0, myState.nb_used ) )
      // found a plane, nothing more to do.
      return true;
    else
      {
	if( ( myState.grad.x == I_ZERO )
	    && ( myState.grad.y == I_ZERO ) )
	  {
	    --myState.nb_used;
	    return false;
	  }
	doubleCut();
      }
  }
  // was unable to find a correct plane.
  --myState.nb_used;
  return false;
}



/**
 * Checks if the given point of myPts belongs to the current
 * plane. Updates myState accordingly.
 *
 * @param idx any valid index of a point of myPts.
 *
 * @return 'true' if this is the case, otherwise returns 'false' and
 * computes the gradient for finding a new valid plane.
 */
bool 
ImaGene::COBAPlaneRecognition::oracleForPoint( int idx )
{
  //look for the points defining the min dot product and the max dot product
  _p.x = myPts[ idx ].coords[ 0 ];
  _p.y = myPts[ idx ].coords[ 1 ];
  _p.z = myPts[ idx ].coords[ 2 ];
  myState.cip.dot_product( _v, myState.N, _p );
  bool changed = false;
  // int axis_min, axis_max;
  // switch ( myAxis ) {
  // case 0: axis_min = 1; axis_max = 2; break;
  // case 1: axis_min = 2; axis_max = 0; break;
  // case 2: axis_min = 0; axis_max = 1; break;
  // }
  // if ( ( _v > myState.max ) || 
  //      ( ( _v == myState.max )
  // 	 && ( myPts[ idx ].coords[ axis_max ] 
  // 	      > myPts[ myState.indmax ].coords[ axis_max ] ) ) ) 
  //   { 
  //     myState.max = _v;  
  //     myState.indmax = idx;
  //   }
  // else if ( ( _v < myState.min ) || 
  // 	    ( ( _v == myState.min )
  // 	      && ( myPts[ idx ].coords[ axis_min ] 
  // 		   > myPts[ myState.indmin ].coords[ axis_min ] ) ) ) 
  //   {
  //     myState.min = _v;
  //     myState.indmin = idx;
  //   }
  
 if ( _v > myState.max ) 
    { 
      myState.max = _v;  
      myState.indmax = idx;
      changed = true;
    }
  else if ( _v < myState.min )
    {
      myState.min = _v;
      myState.indmin = idx;
      changed = true;
    }
 if ( changed )
   {
     for ( unsigned int i = 0; i < idx; ++i )
       if ( myPts[ i ] == myPts[ idx ] )
	 {
	   --myState.nb_used;
	   return true;
	 }
     switch( myAxis ){
     case 0 : _maxValue = _abs( myState.N.x ); break; //3D normal vector
     case 1 : _maxValue = _abs( myState.N.y ); break; //3D normal vector
     case 2 : _maxValue = _abs( myState.N.z ); break; //3D normal vector
     }
     
     // cerr << "[O. " << _p << " _v=" << _v 
     //      << " i=" << myState.indmin << "," << myState.indmax
     //      << " m=" << myState.min << " M=" << myState.max 
     //      << " D=" << _maxValue
     //      << " ]";
     if ( ( myState.max - myState.min ) 
	  < ( get_d( _maxValue ) * myWidth ) ) 
       return true;      // is DP
     else //is not DP
       {
	 // computation of the gradient
	 switch( myAxis ){
	 case 0 : 
	   myState.grad.x = myPts[ myState.indmin ].coords[ 1 ]
	     - myPts[ myState.indmax ].coords[ 1 ];
	   myState.grad.y = myPts[ myState.indmin ].coords[ 2 ]
	     - myPts[ myState.indmax ].coords[ 2 ];
	   break;
	 case 1:
	   myState.grad.x = myPts[ myState.indmin ].coords[ 0 ]
	     - myPts[ myState.indmax ].coords[ 0 ];
	   myState.grad.y = myPts[ myState.indmin ].coords[ 2 ]
	     - myPts[ myState.indmax ].coords[ 2 ];
	   break;
	 case 2:
	   myState.grad.x = myPts[ myState.indmin ].coords[ 0 ] 
	     - myPts[ myState.indmax ].coords[ 0 ];
	   myState.grad.y = myPts[ myState.indmin ].coords[ 1 ]
	     - myPts[ myState.indmax ].coords[ 1 ];
	   break;
	 }
	 // cerr << "[Grad=" << myState.grad << "]";
	 return false;
       }
   }
 return true;
}


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
bool
ImaGene::COBAPlaneRecognition::oracleForPoints( int idx_begin, int idx_end )
{
  _p.x = myPts[ idx_begin ].coords[ 0 ];
  _p.y = myPts[ idx_begin ].coords[ 1 ];
  _p.z = myPts[ idx_begin ].coords[ 2 ];
  myState.cip.dot_product( myState.min, myState.N, _p );
  myState.max = myState.min;
  myState.indmax = myState.indmin = idx_begin;

  //look for the points defining the min dot product and the max dot product
  // int axis_min, axis_max;
  // switch ( myAxis ) {
  // case 0: axis_min = 1; axis_max = 2; break;
  // case 1: axis_min = 2; axis_max = 0; break;
  // case 2: axis_min = 0; axis_max = 1; break;
  // }
  for ( int idx = idx_begin + 1; idx != idx_end; ++idx )
    {
      _p.x = myPts[ idx ].coords[ 0 ];
      _p.y = myPts[ idx ].coords[ 1 ];
      _p.z = myPts[ idx ].coords[ 2 ];
      myState.cip.dot_product( _v, myState.N, _p );
      if ( _v > myState.max ) 
	{ 
	  myState.max = _v;  
	  myState.indmax = idx;
	}
      else if ( _v < myState.min )
	{
	  myState.min = _v;
	  myState.indmin = idx;
	}

      // if ( ( _v > myState.max ) || 
      // 	   ( ( _v == myState.max )
      // 	     && ( myPts[ idx ].coords[ axis_max ] 
      // 		  > myPts[ myState.indmax ].coords[ axis_max ] ) ) ) 
      // 	{ 
      // 	  myState.max = _v;  
      // 	  myState.indmax = idx;
      // 	}
      // else if ( ( _v < myState.min ) || 
      // 		( ( _v == myState.min )
      // 		  && ( myPts[ idx ].coords[ axis_min ] 
      // 		       > myPts[ myState.indmin ].coords[ axis_min ] ) ) ) 
      // 	{
      // 	  myState.min = _v;
      // 	  myState.indmin = idx;
      // 	}
    }

  switch( myAxis ){
  case 0 : _maxValue = _abs( myState.N.x ); break; //3D normal vector
  case 1 : _maxValue = _abs( myState.N.y ); break; //3D normal vector
  case 2 : _maxValue = _abs( myState.N.z ); break; //3D normal vector
  }
  // cerr << "[OALL. " 
  //      << " i=" << myState.indmin << "," << myState.indmax
  //      << " m=" << myState.min << " M=" << myState.max 
  //      << " D=" << _maxValue
  //      << " ]";

  if ( ( myState.max - myState.min ) 
       < ( get_d( _maxValue ) * myWidth ) ) 
    return true;      // is DP
  else //is not DP
  {
    // computation of the gradient
    switch( myAxis ){
    case 0 : 
      myState.grad.x = myPts[ myState.indmin ].coords[ 1 ]
	- myPts[ myState.indmax ].coords[ 1 ];
      myState.grad.y = myPts[ myState.indmin ].coords[ 2 ]
	- myPts[ myState.indmax ].coords[ 2 ];
      break;
    case 1:
      myState.grad.x = myPts[ myState.indmin ].coords[ 0 ]
	- myPts[ myState.indmax ].coords[ 0 ];
      myState.grad.y = myPts[ myState.indmin ].coords[ 2 ]
	- myPts[ myState.indmax ].coords[ 2 ];
      break;
    case 2:
      myState.grad.x = myPts[ myState.indmin ].coords[ 0 ] 
	- myPts[ myState.indmax ].coords[ 0 ];
      myState.grad.y = myPts[ myState.indmin ].coords[ 1 ]
	- myPts[ myState.indmax ].coords[ 1 ];
      break;
    }
    // cerr << "[Grad=" << myState.grad << ", "
    // 	 << myPts[ myState.indmin ] << myPts[ myState.indmax ] << "]";
    return false;
  }
}

/**
 * Performs the double cut in parameter space according to the current
 * gradient and width.
 */
void
ImaGene::COBAPlaneRecognition::doubleCut()
{
  // 2 cuts on the search space:
  //  Gradient.p <= constante
  // -Gradient.p <= constantebis
 
  _cst1 = ( (int) floor( get_si( myG ) * myWidth ) ) - 1;
  _cst2 = ( (int) ceil( get_si( myG ) * myWidth ) ) - 1;

  switch( myAxis ) {
  case 0 : 
    _v = myG * ( myPts[ myState.indmin ].coords[ 0 ] 
		 - myPts[ myState.indmax ].coords[ 0 ] );
    _cst1 -= _v;
    _cst2 += _v;
    break;
  case 1 : 
    _v = myG * ( myPts[ myState.indmin ].coords[ 1 ] 
		 - myPts[ myState.indmax ].coords[ 1 ] );
    _cst1 -= _v;
    _cst2 += _v;
    break;
  case 2 : 
    _v = myG * ( myPts[ myState.indmin ].coords[ 2 ]
		 - myPts[ myState.indmax ].coords[ 2 ] );
    _cst1 -= _v;
    _cst2 += _v;
    break;
  }
  myState.cip.cutOpti( myState.grad, _cst1 );
  myState.grad.neg();
  myState.cip.cutOpti( myState.grad, _cst2 );
  myState.grad.neg();
}


void
ImaGene::COBAPlaneRecognition::computeDigitalNormalVector()
{
  myState.cip.centroid( myState.centroid );
  myState.centroid.reduce();
  switch( myAxis){
  case 0 : 
    myState.N.x = myState.centroid.z * myG;
    myState.N.y = myState.centroid.x;
    myState.N.z = myState.centroid.y;
    break;
  case 1 : 
    myState.N.x = myState.centroid.x;
    myState.N.y = myState.centroid.z * myG;
    myState.N.z = myState.centroid.y;
    break;
  case 2 : 
    myState.N.x = myState.centroid.x;
    myState.N.y = myState.centroid.y;
    myState.N.z = myState.centroid.z * myG;
    break;
  }
}


///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param that_stream the output stream where the object is written.
 */
void 
ImaGene::COBAPlaneRecognition::selfDisplay( ostream& that_stream ) const
{
  double n[ 3 ];
  getNormal( n );
  that_stream << "[ Plane" 
	      << " #pts=" << myState.nb_used 
	      << " #cvx=" << myState.cip.vertices().size()
	      << " N=(" << n[ 0 ] << "," << n[ 1 ] << "," << n[ 2 ]
	      << " ]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::COBAPlaneRecognition::OK() const
{
  return true;
}



///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
