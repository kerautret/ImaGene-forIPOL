///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : COBAPlaneRecognitionWithSet.cxx
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
#include "ImaGene/arithmetic/COBAPlaneRecognitionWithSet.h"
// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/arithmetic/COBAPlaneRecognitionWithSet.ih"
#endif
#define BETTER 1
#define WORSE 0
#define END 2
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const COBAPlaneRecognitionWithSet_RCS_ID = "@(#)class COBAPlaneRecognitionWithSet definition.";



///////////////////////////////////////////////////////////////////////////////
// class COBAPlaneRecognitionWithSet
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::COBAPlaneRecognitionWithSet::~COBAPlaneRecognitionWithSet()
{
  //delete myState.myCorePts;
}

/**
 * Constructor.
 * The object is not valid.
 */
ImaGene::COBAPlaneRecognitionWithSet::COBAPlaneRecognitionWithSet()
{
  myState.myCorePts = CowPtr< std::set< Point3i > >( new std::set< Point3i > );
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
ImaGene::COBAPlaneRecognitionWithSet::init
( int axis, int diameter, double width, const Point3i & pt )
{
  myAxis = axis;
  myWidth = width;
  // initialize the grid step.
  myG = 2*diameter;
  myG *= diameter;
  myG *= diameter;
  // std::cerr << "ImaGene::COBAPlaneRecognitionWithSet::init: myG=" 
  // 	    << myG << std::endl;
  // initialize the search space as a square.
  myState.cip.clear();
  myState.cip.addEnd( Point2I( -myG, -myG ) ); 
  myState.cip.addEnd( Point2I(  myG, -myG ) ); 
  myState.cip.addEnd( Point2I(  myG,  myG ) ); 
  myState.cip.addEnd( Point2I( -myG,  myG ) ); 
  computeDigitalNormalVector();
  //   myState.cip.centroid( myState.centroid );
  //   myState.centroid.reduce();
  myState.myCorePts->clear();
  myState.myPts.clear();
  myState.myCorePts->insert( pt );
  _p.x = pt.coords[ 0 ];
  _p.y = pt.coords[ 1 ];
  _p.z = pt.coords[ 2 ];
  myState.cip.dot_product( _v, myState.N, _p );
  myState.min = myState.max = _v; //myState.N * pt;
  // Const necessary since myCorePts is mutable. Forces to use const
  // methods and prevents a possible copy on write.
  const CowPtr< std::set<Point3i> > & constCorePts = myState.myCorePts;
  myState.indmin = constCorePts->begin();
  myState.indmax = constCorePts->begin();
  myState.grad = Point2I( 0, 0 );
  // cerr << "COBA_INIT " << pt << " " 
  //      <<  "(" <<*myState.indmin << "<" << *myState.indmax << ")"
  //      << *this << endl;
}


/**
 * Checks if  the point [p] is in the current digital plane.
 *
 * @param p any 3D point (in the specified diameter).
 *
 * @return 'true' if it is in the plane, false otherwise.
 */
bool 
ImaGene::COBAPlaneRecognitionWithSet::check( const Point3i & p ) const
{
  //look for the points defining the min dot product and the max dot product
  _p.x = p[ 0 ];
  _p.y = p[ 1 ];
  _p.z = p[ 2 ];
  myState.cip.dot_product( _v, myState.N, _p );
  return ( _v >= myState.min ) && ( _v <= myState.max );
}

/**
   Chooses the "best" band given the current normal, set of points
   and current thickness. Should be called before check.
   
   @see check
*/
void
ImaGene::COBAPlaneRecognitionWithSet::optimizeBand()
{
  computeDigitalNormalVector();
  switch( myAxis ){
  case 0 : _maxValue = _abs( myState.N.x ); break; //3D normal vector
  case 1 : _maxValue = _abs( myState.N.y ); break; //3D normal vector
  case 2 : _maxValue = _abs( myState.N.z ); break; //3D normal vector
  }
  _v = floor( get_d( _maxValue ) * myWidth );
  _v -= myState.max - myState.min;
  if ( _v > I_ZERO )
    {
      // std::cerr << "  COBA (" << myState.min << "," << myState.max << ")";
      _cst1 = _v >> 1;
      _cst2 = _v - _cst1;
      myState.min -= _cst1;
      myState.max += _cst2;
      // std::cerr << "->(" << myState.min << "," << myState.max << ")" << endl;
    }
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
ImaGene::COBAPlaneRecognitionWithSet::add
( const Point3i & p, bool allowNewNormal )
{
  // Add point to myPts.
  // Const necessary since myCorePts is mutable. Forces to use const
  // methods and prevents a possible copy on write.
  const CowPtr< std::set<Point3i> > & constCorePts = myState.myCorePts;
  ConstIterator itcore = constCorePts->lower_bound( p );
  if ( ( itcore != constCorePts->end() )
       && ( *itcore == p ) ) return true;
  Iterator it = myState.myPts.lower_bound( p );
  if ( ( it != myState.myPts.end() )
       && ( *it == p ) ) return true;
  it = myState.myPts.insert( it, p );
  // cerr << "AD3 " << *it << " " 
  //      <<  "(" <<*myState.indmin << "<" << *myState.indmax << ")"
  //      << *this << endl;
  // cerr << "Add " << p << " in " << *this << endl;
//   idx = myState.nb_used;
//   if ( myState.nb_used < myPts.size() )
//     myPts[ myState.nb_used ] = p;
//   else
//     myPts.push_back( p );

  if ( ! allowNewNormal )
    {
      bool ok = oracleForPointNoNewNormal( it );
      if ( ! ok )
	{
	  myState.myPts.erase( it );
	}
      return ok;
    }

  // Checks if new point is still in the plane.
  if ( oracleForPoint( it ) )
    // yes, nothing to do.
    return true;

  // Checks if we can change the normal so as to find another digital plane.
  if( ( ( myState.grad.x == I_ZERO )
	&& ( myState.grad.y == I_ZERO ) ) )
    {
      // last point is removed.
      myState.myPts.erase( it );
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
    if ( oracleForPoints() )
      // found a plane, nothing more to do.
      return true;
    else
      {
	if( ( myState.grad.x == I_ZERO )
	    && ( myState.grad.y == I_ZERO ) )
	  {
	    myState.myPts.erase( it );
	    return false;
	  }
	doubleCut();
      }
  }
  // was unable to find a correct plane.
  myState.myPts.erase( it );
  return false;
}



/**
 * Checks if the given point of myPts belongs to the current
 * plane. Updates myState accordingly.
 *
 * @param p any 3D point (in the specified diameter).
 *
 * @return 'true' if this is the case, otherwise returns 'false' and
 * computes the gradient for finding a new valid plane.
 */
bool 
ImaGene::COBAPlaneRecognitionWithSet::oracleForPoint( Iterator it )
{
  //look for the points defining the min dot product and the max dot product
  _p.x = it->coords[ 0 ];
  _p.y = it->coords[ 1 ];
  _p.z = it->coords[ 2 ];
  myState.cip.dot_product( _v, myState.N, _p );
  bool changed = false;

  //std::cerr << "[Oracle" << "/" << myAxis << " p=" << _p;
  if ( _v > myState.max ) 
    { 
      myState.max = _v;  
      myState.indmax = it;
      changed = true;
      //cerr <<  "(" <<*myState.indmin << "<" << *myState.indmax << "*)";
    }
  else if ( _v < myState.min )
    {
      myState.min = _v;
      myState.indmin = it;
      changed = true;
      //cerr <<  "(" <<*myState.indmin << "*<" << *myState.indmax << ")";
    }
  //cerr << "]";
  if ( changed )
   {
     switch( myAxis ){
     case 0 : _maxValue = _abs( myState.N.x ); break; //3D normal vector
     case 1 : _maxValue = _abs( myState.N.y ); break; //3D normal vector
     case 2 : _maxValue = _abs( myState.N.z ); break; //3D normal vector
     }
     
     if ( ( myState.max - myState.min ) 
	  < ( get_d( _maxValue ) * myWidth ) ) 
       return true;      // is DP
     else //is not DP
       {
	 // computation of the gradient
	 switch( myAxis ){
	 case 0 : 
	   myState.grad.x = myState.indmin->coords[ 1 ]
	     - myState.indmax->coords[ 1 ];
	   myState.grad.y = myState.indmin->coords[ 2 ]
	     - myState.indmax->coords[ 2 ];
	   break;
	 case 1:
	   myState.grad.x = myState.indmin->coords[ 0 ]
	     - myState.indmax->coords[ 0 ];
	   myState.grad.y = myState.indmin->coords[ 2 ]
	     - myState.indmax->coords[ 2 ];
	   break;
	 case 2:
	   myState.grad.x = myState.indmin->coords[ 0 ] 
	     - myState.indmax->coords[ 0 ];
	   myState.grad.y = myState.indmin->coords[ 1 ]
	     - myState.indmax->coords[ 1 ];
	   break;
	 }
	 //cerr << "[Grad=" << myState.grad << "]";
	 return false;
       }
   }
 return true;
}

/**
 * Checks if the given point of myPts belongs to the current
 * plane. Updates myState accordingly.
 *
 * @param p any 3D point (in the specified diameter).
 *
 * @return 'true' if this is the case, otherwise returns 'false' and
 * computes the gradient for finding a new valid plane.
 */
bool 
ImaGene::COBAPlaneRecognitionWithSet::oracleForPointNoNewNormal( Iterator it )
{
  //look for the points defining the min dot product and the max dot product
  _p.x = it->coords[ 0 ];
  _p.y = it->coords[ 1 ];
  _p.z = it->coords[ 2 ];
  myState.cip.dot_product( _v, myState.N, _p );
  
  return ( _v >= myState.min ) && ( _v <= myState.max );
}

/**
 * Checks if all the points of myPts belongs to the current
 * plane. Updates myState accordingly.
 *
 * @return 'true' if this is the case, otherwise returns 'false' and
 * computes the gradient for finding a new valid plane.
 *
 * computeDigitalNormalVector should have been called beforehand.
 */
bool
ImaGene::COBAPlaneRecognitionWithSet::oracleForPoints()
{
  // Const necessary since myCorePts is mutable. Forces to use const
  // methods and prevents a possible copy on write.
  const CowPtr< std::set<Point3i> > & constCorePts = myState.myCorePts;
  ConstIterator it = constCorePts->begin();
  ConstIterator it_end = constCorePts->end();
  _p.x = it->coords[ 0 ];
  _p.y = it->coords[ 1 ];
  _p.z = it->coords[ 2 ];
  myState.cip.dot_product( myState.min, myState.N, _p );
  myState.max = myState.min;
  myState.indmax = myState.indmin = it;
  ++it;

  //look for the points defining the min dot product and the max dot product
  for ( ; it != it_end; ++it )
    {
      _p.x = it->coords[ 0 ];
      _p.y = it->coords[ 1 ];
      _p.z = it->coords[ 2 ];
      myState.cip.dot_product( _v, myState.N, _p );
      if ( _v > myState.max ) 
	{ 
	  myState.max = _v;  
	  myState.indmax = it;
	}
      else if ( _v < myState.min )
	{
	  myState.min = _v;
	  myState.indmin = it;
	}
    }
  
  for ( ConstIterator it = myState.myPts.begin(), it_end = myState.myPts.end();
	it != it_end;
	++it )
    {
      _p.x = it->coords[ 0 ];
      _p.y = it->coords[ 1 ];
      _p.z = it->coords[ 2 ];
      myState.cip.dot_product( _v, myState.N, _p );
      if ( _v > myState.max ) 
	{ 
	  myState.max = _v;  
	  myState.indmax = it;
	}
      else if ( _v < myState.min )
	{
	  myState.min = _v;
	  myState.indmin = it;
	}
    }

  switch( myAxis ){
  case 0 : _maxValue = _abs( myState.N.x ); break; //3D normal vector
  case 1 : _maxValue = _abs( myState.N.y ); break; //3D normal vector
  case 2 : _maxValue = _abs( myState.N.z ); break; //3D normal vector
  }
  

  if ( ( myState.max - myState.min ) 
       < ( get_d( _maxValue ) * myWidth ) ) 
    return true;      // is DP
  else //is not DP
  {
    // computation of the gradient
    switch( myAxis ){
    case 0 : 
      myState.grad.x = myState.indmin->coords[ 1 ]
	- myState.indmax->coords[ 1 ];
      myState.grad.y = myState.indmin->coords[ 2 ]
	- myState.indmax->coords[ 2 ];
      break;
    case 1:
      myState.grad.x = myState.indmin->coords[ 0 ]
	- myState.indmax->coords[ 0 ];
      myState.grad.y = myState.indmin->coords[ 2 ]
	- myState.indmax->coords[ 2 ];
      break;
    case 2:
      myState.grad.x = myState.indmin->coords[ 0 ] 
	- myState.indmax->coords[ 0 ];
      myState.grad.y = myState.indmin->coords[ 1 ]
	- myState.indmax->coords[ 1 ];
      break;
    }
    // cerr << "[Grad=" << myState.grad << ", "
    // 	 << myPts[ myState.indmin ] << myPts[ myState.indmax ] << "]";
    return false;
  }
}

/**
 * Computes the gradient to find a "better" plane (i.e. a thiner one)
 * need to determine the points defining the min and the max dot product
 *
 * @return 'true' if a not null gradient has been computed, 
 * otherwise returns 'false' (it means that the current normal vector is optimal)
 *
 */
bool
ImaGene::COBAPlaneRecognitionWithSet::computeGradient()
{
  // JOL: for now, forbids computeGradient except when myPts is empty
  // (for instance, true after a getState).
  ASSERT_COBAPlaneRecognitionWithSet( myState.myPts.size() == 0 );
  // Const necessary since myCorePts is mutable. Forces to use const
  // methods and prevents a possible copy on write.
  const CowPtr< std::set<Point3i> > & constCorePts = myState.myCorePts;
  ConstIterator it = constCorePts->begin();
  ConstIterator it_end = constCorePts->end();
  _p.x = it->coords[ 0 ];
  _p.y = it->coords[ 1 ];
  _p.z = it->coords[ 2 ];
  myState.cip.dot_product( myState.min, myState.N, _p );
  myState.max = myState.min;
  myState.indmax = myState.indmin = it;
  ++it;
  
	              
  //look for the points defining the min dot product and the max dot product
  for ( ; it != it_end; ++it )
    {
      _p.x = it->coords[ 0 ];
      _p.y = it->coords[ 1 ];
      _p.z = it->coords[ 2 ];

      myState.cip.dot_product( _v, myState.N, _p );
      if ( _v > myState.max ) 
	{ 
	  myState.max = _v;  
	  myState.indmax = it;
	}
      else if ( _v < myState.min )
	{
	  myState.min = _v;
	  myState.indmin = it;
	}
    }

  switch( myAxis ){
  case 0 : _maxValue = _abs( myState.N.x ); break; //3D normal vector
  case 1 : _maxValue = _abs( myState.N.y ); break; //3D normal vector
  case 2 : _maxValue = _abs( myState.N.z ); break; //3D normal vector
  }

    // computation of the gradient
    switch( myAxis ){
    case 0 : 
      myState.grad.x = myState.indmin->coords[ 1 ]
	- myState.indmax->coords[ 1 ];
      myState.grad.y = myState.indmin->coords[ 2 ]
	- myState.indmax->coords[ 2 ];
      break;
    case 1:
      myState.grad.x = myState.indmin->coords[ 0 ]
	- myState.indmax->coords[ 0 ];
      myState.grad.y = myState.indmin->coords[ 2 ]
	- myState.indmax->coords[ 2 ];
      break;
    case 2:
      myState.grad.x = myState.indmin->coords[ 0 ] 
	- myState.indmax->coords[ 0 ];
      myState.grad.y = myState.indmin->coords[ 1 ]
	- myState.indmax->coords[ 1 ];
      break;
    }

    if(myState.grad.x==0 && myState.grad.y==0)
      return false;
    else
      return true;
}

/**
 * Performs the double cut in parameter space according to the current
 * gradient and width.
 */
void
ImaGene::COBAPlaneRecognitionWithSet::doubleCut()
{
  // 2 cuts on the search space:
  //  Gradient.p <= cst1
  // -Gradient.p <= cst2
 
  // _cst1 = ( (int) ceil( get_si( myG ) * myWidth ) + 1 );
  // _cst2 = ( (int) floor( get_si( myG ) * myWidth ) - 1 );
  _cst1 = Integer( ceil( get_d( myG ) * myWidth ) + 1 );
  _cst2 = Integer( floor( get_d( myG ) * myWidth ) - 1 );

  switch( myAxis ) {
  case 0 : 
    _v = myG * ( myState.indmin->coords[ 0 ] 
		 - myState.indmax->coords[ 0 ] );
    _cst1 -= _v;
    _cst2 += _v;
    break;
  case 1 : 
    _v = myG * ( myState.indmin->coords[ 1 ] 
		 - myState.indmax->coords[ 1 ] );
    _cst1 -= _v;
    _cst2 += _v;
    break;
  case 2 : 
    _v = myG * ( myState.indmin->coords[ 2 ]
		 - myState.indmax->coords[ 2 ] );
    _cst1 -= _v;
    _cst2 += _v;
    break;
  }
  myState.cip.cutOpti( myState.grad, _cst1 );
  myState.grad.neg();
  myState.cip.cutOpti( myState.grad, _cst2 );
  myState.grad.neg();
}

/**
 * Performs the simple cut in parameter space according to the current
 * gradient
 */
void
ImaGene::COBAPlaneRecognitionWithSet::simpleCut()
{
  _cst1 = (int) ceil( 
          get_d(myState.grad.x)*get_d(myState.centroid.x)/get_d(myState.centroid.z)
        + get_d(myState.grad.y)*get_d(myState.centroid.y)/get_d(myState.centroid.z)) ;
  _cst2 = (int) floor( 
          get_d(myState.grad.x)*get_d(myState.centroid.x)/get_d(myState.centroid.z)
        + get_d(myState.grad.y)*get_d(myState.centroid.y)/get_d(myState.centroid.z)) ;

  // @todo JOL: conversion (double) -> (int) not checked.
  // std::cerr << "simpleCut() :_cst1=" << _cst1 << " _cst2=" << _cst2
  // 	    << std::endl;
// according to the gradient, we decide to round off the constant in the spacific way
  if( myState.grad.y > 0 || (myState.grad.y == 0 && myState.grad.x > 0) )
  {
	if( _cst1==_cst2 )
	{
	  _cst1++;
	}
 	myState.grad.neg();
	myState.cip.cutOpti( myState.grad, -_cst1 );
  }
  else
  {
    if( myState.grad.y == 0 && myState.grad.x > 0 )
    {
	if( _cst1==_cst2 )
	{
	  _cst2--;
	}
 	myState.grad.neg();
	myState.cip.cutOpti( myState.grad, -_cst2 );
    }
    else
    {
	if( _cst1==_cst2 )
	{
	  _cst1++;
	}
 	myState.grad.neg();
	myState.cip.cutOpti( myState.grad, -_cst1 );
    }
   }
	
   myState.grad.neg();
}

void
ImaGene::COBAPlaneRecognitionWithSet::computeDigitalNormalVector()
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

/**
  * Useful when a solution vector has been found (is DP)
  * Look for the solution vector which minimizes the thickness
  *
  * cuts the search space and tests the centroid
  * returns the associated width or -1 if the cip is not bigger than on single point
 */
double 
ImaGene::COBAPlaneRecognitionWithSet::refineNormal( void )
{
  if( myState.cip.vertices().size() <= 1 ) //no pt will remain after the cut...
    { return -1.0; }
	      
  if(myState.grad.x!=0 || myState.grad.y!=0) //not null gradient
  {
    simpleCut(); //cut the search space
    if( myState.cip.vertices().size() >= 1 )
    {
       //compute the normal vector to test <-> centroid
       computeDigitalNormalVector();
       computeGradient();
       switch( myAxis ){
         case 0 : _maxValue = _abs( myState.N.x ); break;
         case 1 : _maxValue = _abs( myState.N.y ); break;
         case 2 : _maxValue = _abs( myState.N.z ); break;
       }
       return get_d(myState.max - myState.min)/get_d( _maxValue );
    }
  }
  return -1; //no pt remains after the cut or null gradient
}
///////////////////////////////////////////////////////////////////////////////
/** 
    * computes and returns the thickness of the set of points
    */
double 
ImaGene::COBAPlaneRecognitionWithSet::thickness()
    {
      switch( myAxis ){
         case 0 : _maxValue = _abs( myState.N.x ); break;
         case 1 : _maxValue = _abs( myState.N.y ); break;
         case 2 : _maxValue = _abs( myState.N.z ); break;
       }
       return get_d(myState.max - myState.min)/get_d( _maxValue );
    }
    
///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param that_stream the output stream where the object is written.
 */
void 
ImaGene::COBAPlaneRecognitionWithSet::selfDisplay( ostream& that_stream ) const
{
  double n[ 3 ];
  getNormal( n );
  that_stream << "[ Plane" 
	      << " #pts=" << size()
	      << " #cvx=" << myState.cip.vertices().size()
	      << " N=(" << n[ 0 ] << "," << n[ 1 ] << "," << n[ 2 ]
	      << ") ]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::COBAPlaneRecognitionWithSet::OK() const
{
  return true;
}



///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
