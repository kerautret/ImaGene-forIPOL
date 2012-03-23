///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : CompletedDigitalSurfaceGraph.cxx
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
#include "ImaGene/planes/CompletedDigitalSurfaceGraph.h"
// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/planes/CompletedDigitalSurfaceGraph.ih"
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const CompletedDigitalSurfaceGraph_RCS_ID = "@(#)class CompletedDigitalSurfaceGraph definition.";



///////////////////////////////////////////////////////////////////////////////
// class CompletedDigitalSurfaceGraph
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::CompletedDigitalSurfaceGraph::~CompletedDigitalSurfaceGraph()
{
  if ( myTracker != 0 ) delete myTracker;
}

/**
 * Constructor. 
 */
ImaGene::CompletedDigitalSurfaceGraph::CompletedDigitalSurfaceGraph
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
ImaGene::CompletedDigitalSurfaceGraph::selfDisplay( ostream& that_stream ) const
{
  that_stream << "[CompletedDigitalSurfaceGraph]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::CompletedDigitalSurfaceGraph::OK() const
{
  return true;
}



///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
