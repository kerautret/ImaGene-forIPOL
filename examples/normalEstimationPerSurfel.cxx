///////////////////////////////////////////////////////////////////////////////
// Example tcover3d: for a given 3D shape, extract its boundary,
// compute at different scales its tangential cover with maximal disks
// and with extended disk (greedy approach)
///////////////////////////////////////////////////////////////////////////////

#include<cmath>
#include<cstdlib>
#include<iostream>
#include<map>
#include<utility>
#include<fstream>
#include<string>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/ObjectBoundaryTracker.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/planes/DigitalSurfaceGraph.h"
#include "ImaGene/planes/NuThickDisk.h"
#include "ImaGene/planes/TangentialCoverHierarchy.h"
#include "ImaGene/mathutils/Statistic.h"
#include "ImaGene/timetools/Clock.h"

#include "ImaGene/planeRecognition/planeRecognitionOnSurfaces.h"

using namespace std;
using namespace ImaGene;

static Arguments args;

/**
* check if the normalVector has the right sign relative to the surf,i.e. the normal vector of surf and 
* normalVector have to point at the same half plane 
* change it if it's not the case
* @param ks a KnSpace
* @param surf a surfel
* @param normalVector a real vector
*/
void checkSign(const KnSpace & ks, Kn_sid surf,vector<double> & normalVector)
{
  //get the normal vector of the surfel
  Vector vAxe = ks.sorthVectorBasis(surf); 
  //check the sign of the dot product
  double dotproduct = vAxe.ro(0)*normalVector[0] + vAxe.ro(1)*normalVector[1] + vAxe.ro(2)*normalVector[2];
  if(dotproduct < 0)
  {
    normalVector[0] = -normalVector[0];
    normalVector[1] = -normalVector[1];
    normalVector[2] = -normalVector[2];
  }
  if ( ks.sdirect( surf, ks.sorthDir(surf) ) )
  {
    normalVector[0] = -normalVector[0];
    normalVector[1] = -normalVector[1];
    normalVector[2] = -normalVector[2];
  }
}

void addnoise( KnSpace & ks, KnCharSet & voxset, int noise )
{
  //DigitalSurfaceTracker initialization
  ObjectBoundaryTracker trackerInit;
//   BelAdjacency baInit(ks,false);
//   ObjectBoundary obInit(baInit,shape);
//   trackerInit.init(&obInit);
//   trackerInit.move(b);
  KnRCellSet bdryInit = KnShapes::smakeBoundary( ks, voxset );
    
  Kn_uid in, out;
  vector<double> p_in( 6 );
  vector<double> p_out( 6 );
  switch(noise){
  case 1 : cout<<"noise level 1"<<endl;
    p_in[0] = 0.0;
    p_in[1] = 0.2;
    p_in[2] = 0.1;
    p_in[3] = 0.05;
    p_in[4] = 0.025;
    p_in[5] = 0.0125;
    
    p_out[0] = 0.0;
    p_out[1] = 0.2;
    p_out[2] = 0.1;
    p_out[3] = 0.05;
    p_out[4] = 0.025;
    p_out[5] = 0.0125;
    break;
  case 2 :  cout<<"noise level 2"<<endl;
    p_in[0] = 0.0;
    p_in[1] = 0.5;
    p_in[2] = 0.15;
    p_in[3] = 0.0;
    p_in[4] = 0.0;
    p_in[5] = 0.0;
    
    p_out[0] = 0.0;
    p_out[1] = 0.3;
    p_out[2] = 0.1;
    p_out[3] = 0.05;
    p_out[4] = 0.0;
    p_out[5] = 0.0;
    break;
  default : cout<<"noise level 0"<<endl;
    p_in[0] = 0.0;
    p_in[1] = 0.0;
    p_in[2] = 0.0;
    p_in[3] = 0.0;
    p_in[4] = 0.0;
    p_in[5] = 0.0;
    
    p_out[0] = 0.0;
    p_out[1] = 0.0;
    p_out[2] = 0.0;
    p_out[3] = 0.0;
    p_out[4] = 0.0;
    p_out[5] = 0.0;
  }
  KnCharSet shapeWithNoise = KnShapes::noisifyObject
    ( ks, voxset, bdryInit, p_in, p_out, in, out); 
  cout << " --- noisy shape has " << shapeWithNoise.nbElements() << " voxels.--" << endl;
  KnCharSet main_inner_comp = KnCharSet::create (ks,true,3,0);
  Kn_sid b = ShapeHelper::findInnerObject
    (&ks, shapeWithNoise, in, main_inner_comp);
  voxset = main_inner_comp;
} //if(args.check("-noise"))
 
/**
 * s = "1.0, 1.5, 2.0, 2.5, 3.0, 3.5"
 * implies
 * t[] = { ..., 1.0, 1.5, 2.0, 2.5, 3.0, 3.5 }
 */
void string2doubles( vector<double> & t, const string & s, char sep = ',' )
{
  size_t last_pos = 0;
  size_t pos = ( s == "" ) ? string::npos : 0;
  while ( pos != string::npos )
    {
      pos = s.find_first_of( ',', last_pos );
      string v = s.substr( last_pos, pos - last_pos );
      // cerr << pos << ":" << v << endl;
      t.push_back( atof( v.c_str() ) );
      last_pos = pos + 1;
    }
}

void computeTgtCover3D( const KnSpace & ks, KnCharSet voxset, 
			vector<double> & scales )
{
  // WITH NU-THICK DISKS.
  // boundary is needed to get the set of all vertices.
  KnRCellSet bdry = KnShapes::smakeBoundary( ks, voxset );
  // the tracker is used for extracting neighbors within the graph of
  // bel-adjacency.
  ObjectBoundaryTracker tracker;
  BelAdjacency ba( ks, false );
  ObjectBoundary ob( ba, voxset );
  tracker.init( &ob );
  tracker.move( *( bdry.begin() ) );
  DigitalSurfaceGraph dsg( ks, bdry, tracker );
  NuThickDisk<DigitalSurfaceGraph> disk;
  typedef NuThickDisk<DigitalSurfaceGraph>::WeightedVertex WeightedCell;
  priority_queue< WeightedCell > fp;
  //normal vector of the recognized piece of plane
  vector<double> normalVector(3); 

  KnRCellSet exp = KnRCellSet::create( ks, 2, true, 0 );
  for ( unsigned int i = 0; i < scales.size(); ++i )
    {
      Statistic<float> Xradius;
      cerr << "---- Scale " << i << " is " << scales[ i ] << " ----" << endl;
      Clock::startClock();
      for ( KnRCellSet::cell_iterator p = bdry.begin(); p != bdry.end(); ++p )
	{
	  Kn_sid bel = *p;
	  disk.init( dsg, bel, scales[ i ] );
	  disk.computeDisk( ks.sorthDir( bel ) );
	  disk.getNormal( normalVector );
	  checkSign( ks, bel, normalVector );
	  // // favoriteAxe[itOnMap->first] = disk.majorAxis();
	  // fp.push( weightedCell( bel,  disk.radius() ) );
	  Xradius.addValue( disk.radius() );
	  // disk.computeExpansion( exp );
	  //	  cerr << "  exp=" << exp.nbElements() << endl;
// 	  cerr << "Cell=" << bel << " r=" << disk.radius() 
// 	       << " #=" << disk.size() 
// 	       << " axis=" << disk.majorAxis() 
// 	       << " N=(" << normalVector[ 0 ] 
// 	       << "," << normalVector[ 1 ] 
// 	       << "," << normalVector[ 2 ] << ")" << endl;
	}
      long t1 = Clock::stopClock();
      cerr << "   nb_computed = " << Xradius.samples() << endl;
      cerr << "   Xradius = " << Xradius.mean() 
	   << " +/- " << sqrt( Xradius.unbiasedVariance() ) << endl;
      cerr << "   in " << t1 << " ms." << endl;
    }

}


// void computeNormals
// ( const KnSpace & ks, KnCharSet voxset, vector<double> & scales, char* fileName )
// {
//   map<int,double> myMap;
//   map< int,map< double,Statistic<float> > > greedyNxMap;
//   map< int,map< double,Statistic<float> > > greedyNyMap;
//   map< int,map< double,Statistic<float> > > greedyNzMap;
//   // boundary is needed to get the set of all vertices.
//   KnRCellSet bdry = KnShapes::smakeBoundary( ks, voxset );
//   // the tracker is used for extracting neighbors within the graph of
//   // bel-adjacency.
//   ObjectBoundaryTracker tracker;
//   BelAdjacency ba( ks, false );
//   ObjectBoundary ob( ba, voxset );
//   tracker.init( &ob );
//   tracker.move( *( bdry.begin() ) );
//   DigitalSurfaceGraph dsg( ks, bdry, tracker );
//   DigitalSurfaceGraphMethods methods( dsg );
//   TangentialCoverHierarchy<DigitalSurfaceGraph> tch;
//   tch.init( dsg, scales );
//   tch.computeNormals( methods, greedyNxMap, greedyNyMap, greedyNzMap );
//   
//   //write in filename : s Nx Ny Nz
//   //with s the code of the surfel and (Nx,Ny,Nz) the associated normal vector
//   
//   for(vector<double>::iterator scaleIt= scales.begin();scaleIt!=scales.end();scaleIt++)
//     {
//       double w = *scaleIt;
//       char fullFileName[150];
//       char numero[15];
//       strcpy(fullFileName,fileName);
//       sprintf(numero,"%.0f",w);
// //       cout<<"extension "<<numero<<endl;
//       strcat(fullFileName,numero);
//       strcat(fullFileName,".txt");
// //       cout<<"ecriture ds fichier "<<fullFileName<<endl;
//       
//       ofstream f(fullFileName);
//       for(map< int,map< double,Statistic<float> > >::iterator it= greedyNxMap.begin();it!=greedyNxMap.end();it++)
//       {
//          int v = it->first; //surfel
//          f<<v<<" ";
//          f<<greedyNxMap[v][w].mean()<<" "<<greedyNyMap[v][w].mean()<<" "<<greedyNzMap[v][w].mean()<<endl;
//       }
//       f.close();
//     }
// }

void computeNormals2
( const KnSpace & ks, KnCharSet voxset, vector<double> & scales, char* fileName )
{
//   cout<<"compute Normals2"<<endl;
  map< int,map< double,Statistic<float> > > greedyNxMap;
  map< int,map< double,Statistic<float> > > greedyNyMap;
  map< int,map< double,Statistic<float> > > greedyNzMap;
  map< double,Statistic<float> > littleMap;
  Statistic<float> fstat;
  
  // boundary is needed to get the set of all vertices.
  KnRCellSet bdry = KnShapes::smakeBoundary( ks, voxset );
  // the tracker is used for extracting neighbors within the graph of
  // bel-adjacency.
  ObjectBoundaryTracker tracker;
  BelAdjacency ba( ks, false );
  ObjectBoundary ob( ba, voxset );
  tracker.init( &ob );
  tracker.move( *( bdry.begin() ) );
  
  /** test affichage**/
//   planeRecognitionOnSurfaces::printVoxelSetInFileTracker(ks,tracker,"/home2/echar/JOL-XP/imageneutils/tests/sphere.txt" );
//   return;
  /**fin test**/
  
  DigitalSurfaceGraph dsg( ks, bdry, tracker );
  DigitalSurfaceGraphMethods methods( dsg );
  TangentialCoverHierarchy<DigitalSurfaceGraph> tch;
  tch.init( dsg, scales );
  tch.compute( methods);

  //write in filename : s Nx Ny Nz
  //with s the code of the surfel and (Nx,Ny,Nz) the associated normal vector
  
  NuThickDisk<DigitalSurfaceGraph> disk;
  vector<double> normalVector(3); 
  for(vector<double>::iterator scaleIt= scales.begin();scaleIt!=scales.end();scaleIt++)
    {
      double w = *scaleIt;
      KnRCellSet * ptActVtx; //set of active vertices
      KnRCellSet actVtx = tch.activeVertices((int)w+1);
      ptActVtx = & actVtx;
      for ( KnRCellSet::iterator p = ptActVtx->begin();
	    p != ptActVtx->end(); ++p )
      {
        int vtx = *p;
	disk.init( dsg, vtx, w);
	disk.computeDisk( -1.0 );
	disk.getNormal( normalVector );
 	methods.adjustOrientation( vtx, normalVector );
	KnRCellSet * ptExp;
	KnRCellSet exp = KnRCellSet::create( ks, 2, true, 0 );
	disk.computeExtension(exp);
	ptExp = &exp;
	for ( KnRCellSet::iterator q = ptExp->begin();
	    q != ptExp->end(); ++q )
	{
	    //associer le vecteur normal au sommet
	    int v = *q;

	    greedyNxMap.insert(make_pair(v,littleMap));
	    greedyNxMap[v].insert(make_pair(w,fstat));
	    greedyNxMap[v][w].addValue(normalVector[0]);
	  
	    greedyNyMap.insert(make_pair(v,littleMap));
	    greedyNyMap[v].insert(make_pair(w,fstat));
	    greedyNyMap[v][w].addValue(normalVector[1]);
	  
	    greedyNzMap.insert(make_pair(v,littleMap));
	    greedyNzMap[v].insert(make_pair(w,fstat));
	    greedyNzMap[v][w].addValue(normalVector[2]);
	} // for ( VertexConstIterator q = expVertices->begin(); q != expVertices->end(); ++q )
      } //for ( Set<int>::iterator p = actVtx->begin(); p != actVtx->end(); ++p )

      char fullFileName[150];
      char numero[20];
      strcpy(fullFileName,fileName);
      sprintf(numero,"%.0f",w);
      strcat(fullFileName,numero);
      strcat(fullFileName,".txt");
     
  
      ofstream f(fullFileName);
      for(map< int,map< double,Statistic<float> > >::iterator it= greedyNxMap.begin();it!=greedyNxMap.end();it++)
      {
         int v = it->first; //surfel
         if(!isnan(greedyNxMap[v][w].mean()))
	 {
           f<<v<<" ";
           f<<greedyNxMap[v][w].mean()<<" "<<greedyNyMap[v][w].mean()<<" "<<greedyNzMap[v][w].mean()<<endl;
	 }
      }
      f.close();
    }
    cout<<"done"<<endl;
}



//main function
int main(int argc, char* argv[])
{
  StandardArguments::addDigitalArgs( args, 3, false, false );
  ShapeHelper::addSimple3DShapesArgs( args );
  args.addOption( "-voxset", "-voxset <filename_pgm3d> <min> <max>: the input shape is read from a pgm3d file. If I is the voxel value, the voxel is added to the shape whenever min <= I <= max.", "", "0", "128" );
  args.addOption( "-noise", "-noise <level>: adds noise to the shape.", "0");
  args.addOption( "-scales", "-scales <scales>: compute hierarchical tangential cover at the specified scales. Default is 1,2,3,4,5,6,7,8.", "1,2,3,4,5,6,7,8");
  args.addOption( "-file","-file <fileName>: writes the voxel coordinates and the estimated normal vector in the file fileName","none.txt");
  if ( ! args.readArguments( argc, argv ) ) 
    {
      cout << args.usage( "exempleWith args", 
			  "-d 3 -x <size> -y <size> -z <size> -sphere <r> -noise <level> -file <fileName>",
			  "" )
	   << endl;
      return 1;
    }
    
  // Creates a 3D space
  const uint Dim = StandardArguments::dim( args );
  if ( Dim != 3 )
    {
      cerr << "Dimension should be 3." << endl;
      return 2;
    }
  KnSpace* ks = 0;
  KnCharSet* voxset = 0;
  if ( args.check( "-voxset" ) )
    { // init from file.
      string filename = args.getOption( "-voxset" )->getValue( 0 );
      uint smin = args.getOption( "-voxset" )->getIntValue( 1 );
      uint smax = args.getOption( "-voxset" )->getIntValue( 2 );
      ifstream ifs;
      ifs.open(filename.c_str(), ifstream::in);
      if ( ! ShapeHelper::importFromPGM3d( ifs, ks, voxset, smin, smax, 1 ) )
	{
	  cerr << "[tcover3d] Error reading file <" << filename << ">" << endl;
	  return 1;
	}
      cout<<"sizes :"<<(*ks).size(0)<<" "<<(*ks).size(1)<<" "<<(*ks).size(2)<<endl;
    }
  else
    { // init from parameters.
      Kn_size sizes[ Dim ]; 
      StandardArguments::fillSizes( args, sizes );
      ks = new KnSpace( Dim, sizes );
      //create the shape  
      voxset = new KnCharSet
	( ShapeHelper::makeSimple3DShapesFromArgs(args, *ks) );
    }
  if ( ( ks == 0 ) || ( voxset == 0 ) )
    {
      cerr << "[normalEstimationPerSurfel] Unable to create space or shape." << endl;
      return 2;
    }

  //check if some noise must be added
  if ( args.check("-noise") )
    {
      int noise = args.getOption("-noise")->getIntValue( 0 );
      addnoise( *ks, *voxset, noise );
    }
  uint nb_elements = voxset->nbElements();
  cerr << "-- Shape has " << nb_elements << " voxels." << endl;
  if ( nb_elements == 0 )
    {
      cerr << "[normalEstimationPerSurfel] Shape has zero elements." << endl;
      delete voxset;
      delete ks;
      return 3;
    }
  //voxel type in the KnSpace
  Kn_uid typeVox = ks->uspel( ks->ufirst() );
  // Extracts one spel of the shape
  Kn_uid anyVox = voxset->ubegin ();
  // Extracts one bel on the shape boundary.
  Kn_sid b= KnShapes::sfindFurthestBel( *ks, anyVox, *voxset);

  vector<double> scales;
  string2doubles( scales, args.getOption( "-scales" )->getValue( 0 ), ',' );
  
  char fileName[100];
  //check if normal estimation must be computed
  if ( args.check("-file") )
  {
     strcpy(fileName,(args.getOption("-file")->getValue(0)).c_str());
     computeNormals2( *ks, *voxset, scales, fileName );
  }
  
  delete voxset;
  delete ks;
  return 0;
}
