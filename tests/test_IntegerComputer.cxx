///////////////////////////////////////////////////////////////////////////////
// Test module for some arithmetic properties
///////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cmath>
#include <cstdlib>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/arithmetic/IntegerComputer.h"
#include "ImaGene/arithmetic/ConvexIntegerPolygon.h"
#include "ImaGene/arithmetic/COBAPlaneRecognition.h"

using namespace std;
using namespace ImaGene;

static Arguments args;

bool 
testCOBAPlaneRecognition()
{
  bool is_plane;
  COBAPlaneRecognition coba;
  typedef COBAPlaneRecognition::Point3i Point3i;
  coba.init( 2, 100, 1.6, Point3i( 0, 0, 0 ) );
  cout << "I " << coba << endl;
  is_plane = coba.add( Point3i( 1, 0, 0 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 1, 1, 0 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 2, 0, 0 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 2, 1, 0 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 3, 1, 1 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 3, 2, 2 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 1, 3, 2 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 5, 5, 4 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 5, 6, 4 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 5, 7, 4 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 5, 8, 4 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 5, 9, 4 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 5, 10, 4 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 5, 11, 4 ), true );
  cout << is_plane << " " << coba << endl;
  is_plane = coba.add( Point3i( 5, 12, 4 ), true );
  cout << is_plane << " " << coba << endl;
  return true;
}

int 
main( int argc, char** argv ) 
{
  if ( ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_IntegerComputer", 
			  "Tests the class IntegerComputer.",
			  "" )
	   << endl;
      return 1;
    }
  IntegerComputer computer;
  ConvexIntegerPolygon cvx_poly;
  bool ok = computer.selfTest()
    && cvx_poly.selfTest()
    && testCOBAPlaneRecognition();
  return 0;
}

