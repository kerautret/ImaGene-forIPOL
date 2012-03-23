///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : DigitalMidSurfaceGraph.cxx
//
// Creation : 2011/01/13
//
// Version : 2011/01/13
//
// Author : Jacques-Olivier Lachaud
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
//	2011/01/13 : Mr ?Name? : ?What?
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //


///////////////////////////////////////////////////////////////////////////////
#include "ImaGene/planes/DigitalMidSurfaceGraph.h"
// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/planes/DigitalMidSurfaceGraph.ih"
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const DigitalMidSurfaceGraph_RCS_ID = "@(#)class DigitalMidSurfaceGraph definition.";



///////////////////////////////////////////////////////////////////////////////
// class DigitalMidSurfaceGraph
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::DigitalMidSurfaceGraph::~DigitalMidSurfaceGraph()
{
  if ( myTracker != 0 ) delete myTracker;
}

/**
 * Constructor. 
 */
ImaGene::DigitalMidSurfaceGraph::DigitalMidSurfaceGraph
( const KnSpace & ks, 
  KnRCellSet surface, 
  const DigitalSurfaceTracker & tracker )
  : mySpace( ks ), mySurface( surface ), myDiameter( 0 )
{
  myTracker = tracker.clone();
}


///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param that_stream the output stream where the object is written.
 */
void 
ImaGene::DigitalMidSurfaceGraph::selfDisplay( ostream& that_stream ) const
{
  that_stream << "[DigitalMidSurfaceGraph]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::DigitalMidSurfaceGraph::OK() const
{
  return true;
}



///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
