///////////////////////////////////////////////////////////////////////////////
// Generates contours from pgm image
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/Proxy.h"
#include "ImaGene/base/Vector.h"
#include "ImaGene/base/Shapes.h"
#include "ImaGene/digitalnD/GridEmbedder.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/C4CIteratorOnBdry.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/helper/ContourHelper.h"
#include "ImaGene/helper/ShapeHelper.h"

using namespace std;
using namespace ImaGene;


static Arguments args;

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// M A I N
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char** argv ) 
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  args.addOption( "-threshold", "-threshold <val>: threshold value for binarizing PGM gray values (def. is 128).", "128" );
  args.addOption( "-min_size", "-min_size <m>: minimum digital length of contours for output (def. is 4).", "4" );
  args.addOption( "-badj", "-badj <0/1>: 0 is interior bel adjacency, 1 is exterior (def. is 0).", "0" );
  args.addOption("-selectContour", "-selectContour <x0> <y0> <distanceMax>: select the contours for which the first point is near (x0, y0) with a distance less than <distanceMax>","0", "0", "0" );
  args.addBooleanOption("-invertVerticalAxis", "-invertVerticalAxis used to transform the contour representation (need for DGtal), used only for the contour displayed, not for the contour selection (-selectContour). ");
  
  
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "pgm2freeman", 
			  "Extracts all 2D contours from a PGM image given on the standard input and writes them on the standard output as FreemanChain's.",
			  "" )
	   << endl;
      return 1;
    }
  
  bool yInverted = args.check("-invertVerticalAxis");
  KnSpace* ks;
  KnCharSet* voxset;
  uint threshold = (uint) args.getOption( "-threshold" )->getIntValue( 0 );
  if ( ! ShapeHelper::importFromPGM( cin, ks, voxset, threshold, 1 ) )
    {
      cerr << "Error reading PGM file." << endl;
      return 2;
    }
  
  Vector2i ptReference;
  double distanceMax=0.0;
  //Rajout (BK) 
  if(args.check("-selectContour")){
    ptReference.x()= args.getOption("-selectContour")->getIntValue(0);
    ptReference.y()= args.getOption("-selectContour")->getIntValue(1);
    distanceMax= args.getOption("-selectContour")->getIntValue(2);
  }

  bool interior = args.getOption( "-badj" )->getIntValue( 0 ) == 0;
  uint min_size = args.getOption( "-min_size" )->getIntValue( 0 );
  BelAdjacency badj( *ks, interior );
  KnRCellSet bdry = KnShapes::smakeBoundary( *ks, *voxset );
  KnRCellSet not_visited( bdry );
  uint num_contour = 0;
  for ( KnRCellSet::cell_iterator cell_it = bdry.begin();
	cell_it != bdry.end();
	++cell_it )
    {
      Kn_sid bel = *cell_it;
      uint k = *( ks->sbegin_dirs( bel ) );
      C4CIteratorOnBdry c4c_it( badj, bel, k, *voxset );
      bool is_open;
      uint nb_surfels = C4CIterator::size( c4c_it, is_open );
      if ( nb_surfels >= min_size )
	{
	  Proxy<C4CIteratorOnSurface> cp
	    ( (C4CIteratorOnSurface*) c4c_it.clone() );
	  if(!args.check("-selectContour")){
	    ContourHelper::displayFreemanChain( cout, ks, cp, 0, 1, yInverted );
	  }else{
	    //Rajout option (BK 29/07/09)
	    Frame2D frame;
	    frame.init( ks, 0, 1 );
	    Kn_sid sbel = cp->current();
	    frame.setSurfelFrame( sbel, cp->trackDir() );
	    Vector2i p1( frame.transformPoint( Vector2i( 0, 0 ) ) );
	    double distance = sqrt((p1.x() - ptReference.x())*(p1.x() - ptReference.x())+
				   (p1.y() - ptReference.y())*(p1.y() - ptReference.y()));
	    if(distance< distanceMax){
	      ContourHelper::displayFreemanChain( cout, ks, cp, 0, 1, yInverted );
	    }
	  }
	}
      // Clear contour from set of bels.
      bel = c4c_it.current();
      Kn_sid sbel = bel;
      do
	{
	  bdry[ bel ] = false;
	  if ( c4c_it.next() == 0 ) break;
	  bel = c4c_it.current();
	}
      while ( bel != sbel );

      num_contour++;
    }



  return 0;
}
  
