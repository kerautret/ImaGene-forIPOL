///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : CDigitalGraph.cxx
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
#include "ImaGene/planes/CDigitalGraph.h"
// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/planes/CDigitalGraph.ih"
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const CDigitalGraph_RCS_ID = "@(#)class CDigitalGraph definition.";



///////////////////////////////////////////////////////////////////////////////
// class CDigitalGraph
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::CDigitalGraph::~CDigitalGraph()
{
}



///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param that_stream the output stream where the object is written.
 */
void 
ImaGene::CDigitalGraph::selfDisplay( ostream& that_stream ) const
{
  that_stream << "[CDigitalGraph]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::CDigitalGraph::OK() const
{
  return true;
}



///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
