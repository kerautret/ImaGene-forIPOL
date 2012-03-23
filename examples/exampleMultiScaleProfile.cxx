///////////////////////////////////////////////////////////////////////////////
// Example exampleMultiScaleProfile: for a given 3D shape, extract its boundary, compute at different scales 
// its tangential cover with maximal disks and with extended disk (greedy approach)
// execution : 
// ./exampleMultiScaleProfile -d 3 -x <size> -y <size> -z <size> -ellipse/cube/sphere <param> -ry <> -rz <> -noise <0,1 or 2> -file <begginningNameOfFile>
// example : ./exampleMultiScaleProfile -d 3 -x 50 -y 50 -z 50 -sphere 15 -noise 1 -file "file1"
// computes a sphere of radius 15 with noise, its tangential cover at different scales 
// and the multiscale profile of 10 points on the surface is stored in files "file1-i.txt", 1<=i<=10
// the multiscale profile is, for each voxel :
// scale nbCoveringMaxDisks minDiskArea meanDiskArea maxDiskArea varianceDiskArea nbCoveringGreedyPlanes  minGreedyPlaneArea meanGreedyPlaneArea maxGreedyPlaneArea varianceGreedyArea
// the file "file1-info.txt" contains the code of the voxels corresponding to 
// the multiscale profiles
///////////////////////////////////////////////////////////////////////////////

#include<iostream>
#include <gmp.h>  //gmp library
#include <gmpxx.h> //gmp library c++
#include <cmath>
#include<cstdlib>
#include<map>
#include<utility>
#include<fstream>
#include<string>

#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"

#include "ImaGene/planeRecognition/utils.h"
#include "ImaGene/planeRecognition/Z3.h"
#include "ImaGene/planeRecognition/COBA.h"
#include "ImaGene/planeRecognition/planeRecognitionOnSurfaces.h"

#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/KnShapes.h"

#include "ImaGene/digitalnD/DigitalSurfaceTracker.h"
#include "ImaGene/digitalnD/ObjectBoundaryTracker.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/BelAdjacency.h"

#include "ImaGene/mathutils/Statistic.h"

#include "ImaGene/helper/ScaleProfile.h"
#include "ImaGene/helper/ShapeHelper.h"
#include "ImaGene/planes/DigitalSurfaceGraph.h"
#include "ImaGene/planes/NuThickDisk.h"

#include <unistd.h>

#define MAX 100

using namespace std;
using namespace ImaGene;

static Arguments args;

/*new types*/
typedef Statistic<float> FStatistic;

typedef map<Kn_sid,FStatistic> mapVoxStat;

typedef map<float,FStatistic> mapEpStat;

typedef map<Kn_uid,mapEpStat> mapVoxEpStat;

typedef map<Kn_uid, Kn_sid> mapVoxSurf;

typedef map<Kn_uid, int> mapVoxInt;

/*classes uses in the construction of a priority queue*/
class weightedCell // voxels weighted by a Euclidean distance
{
    public:
      Kn_uid v;
      double dist;
    
    weightedCell(Kn_uid vox, double distance){v=vox; dist=distance;}
    weightedCell(const weightedCell & other){v=other.v; dist=other.dist;}
};

/*comparison between weighted cells*/
class mycomparisonCell
{
    public:
      bool operator() (const weightedCell & b1, const weightedCell & b2) const
      {
        return b2.dist > b1.dist;
      }
};  


/**
* computes and returns the area of a surface projceted on a 
* Euclidean plane of normal vector normalVector.
* @param ks a KnSpace
* @param bdry a KnRCellSet representaing the surface
* @param normalVector the normal vector of the Euclidean plane
* @return the projected area.
*/
float computeArea(const KnSpace & ks, const KnRCellSet & bdry, vector<double> & normalVector)
{
  KnRCellSet::cell_iterator it = bdry.begin();  //iterator on the boundary
  Kn_sid surf;
  Kn_sign sign;
  int simpleSign;
  uint n;
  float aire = 0.0;
  
  while(it != bdry.end()) //scan the bdry
  {
    surf = it.get();
    sign = ks.decodeSign(surf);
    if(sign==0)
      simpleSign = 1;
    else
      simpleSign = -1;
    n = ks.sorthDir(surf);
    aire += normalVector[n] * simpleSign;
    ++it;
  }
  return fabs(aire);
}

/**
* scan the surface (surfels) associated to the tracker and associates to each voxel its most 
* significant surfel (best normal vector) and its normal direction
* a surfel is the most significant if it is the only surfel associated to the voxel, or the slopes of 2
* DSS centered of the surfels on the surface are <=1
* @param ks a KnSpace
* @param tracker a tracker on a surface (surfels)
* @param correspondingSurf map associating voxels to their best surfel
* @param favoriteAxe map associating voxels to their best normal direction
* @param shape the KnCharSet corresponding to the surface
*/
void computeFavoriteAxes(const KnSpace & ks, const DigitalSurfaceTracker & tracker, mapVoxSurf & correspondingSurf, mapVoxInt & favoriteAxe, KnCharSet shape)
{
  DigitalSurfaceTracker * p = tracker.clone();
  Kn_sid b = p->current();  // current tracked bel 
  Kn_sid nsurf;   // next tracked bel 
  uint track_dir;  // current tracked direction
  queue<Kn_sid> qbels; //queue for tracking
  
  C4CSegment centeredSegment1, centeredSegment2; //maximal segments used to initialize a digital plane
  uint n; //directions of tracking and normal direction of the reference surfel b
  Kn_size xyzv[3]; //coordinates of the corresponding voxel
  Kn_uid vox; //corresponding voxel
  
  // surfels being explored.
  KnRCellSet explored = KnRCellSet::create( ks, ks.dim() - 1, true, 0 );
  
  //beginning with the initial surfel
  explored += b;
  qbels.push( b );
  n = ks.sorthDir(b); //normal of the surfel
  vox = ks.unsigns( ks.sincident( b, 
	                              n,
				      ks.sdirect( b, n )
				      ) ); //associated voxel
  
  correspondingSurf.insert(make_pair(vox,b));
  favoriteAxe.insert(make_pair(vox,n));
  
  // For all pending bels
  while( ! qbels.empty() )
    {
      b = qbels.front();
      qbels.pop();
      p->move( b );
      // Loop in all tracking directions. (Generic tracking).
      for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end(); ++q )
	{
	  track_dir = *q;
	  nsurf = p->adjacent( track_dir, false );
	  // If the surfel exists and is not already explored,
	  // checks if it is the most significant for its associated voxel.
	  if ( ( p->lastCode() != 0 ) && ( ! explored[ nsurf ] ) )
	    {
	      n = ks.sorthDir(nsurf); //normal of the surfel
	      vox = ks.unsigns( ks.sincident( nsurf, 
	                              n,
				      ks.sdirect( nsurf, n )
				      ) );
	      explored += nsurf;
	      qbels.push( nsurf );
	      
	      if(correspondingSurf.find(vox)==correspondingSurf.end()) //voxel not yet in the maps
              {      
	        correspondingSurf.insert(make_pair(vox,nsurf));
		favoriteAxe.insert(make_pair(vox,n));
              }
              else
	      {
	          planeRecognitionOnSurfaces::computeMaximalSegmentsIterator( ks, nsurf, shape, centeredSegment1, centeredSegment2); //computes 2 DSS centered on nsurf
	          if( (abs((float(centeredSegment1.a()) / (float)centeredSegment1.b())) <= 1) 
                   && (abs((float(centeredSegment2.a()) / (float)centeredSegment2.b())) <= 1) ) //slopes <= 1
                   {      
		     correspondingSurf[vox] = nsurf;
		     favoriteAxe[vox] = n;
		   }
	       }
	    }
 	  nsurf = p->adjacent( track_dir, true );
	  // If the surfel exists and is not already explored,
	  // checks if it is the most significant for its associated voxel.
 	  if ( ( p->lastCode() != 0 ) && ( ! explored[ nsurf ] ) )
 	    {
	      n = ks.sorthDir(nsurf); //normal of the surfel
	      vox = ks.unsigns( ks.sincident( nsurf, 
	                              n,
				      ks.sdirect( nsurf, n )
				      ) );
 	      explored += nsurf;
 	      qbels.push( nsurf );	
		    
	      if(correspondingSurf.find(vox)==correspondingSurf.end()) //voxel not yet in the maps
              {      
	        correspondingSurf.insert(make_pair(vox,nsurf));
		favoriteAxe.insert(make_pair(vox,n));
              }
              else
	      {
	          planeRecognitionOnSurfaces::computeMaximalSegmentsIterator( ks, nsurf, shape, centeredSegment1, centeredSegment2); //computes 2 DSS centered on nsurf
	          if( (abs((float(centeredSegment1.a()) / (float)centeredSegment1.b())) <= 1) 
                   && (abs((float(centeredSegment2.a()) / (float)centeredSegment2.b())) <= 1) )  //slopes <= 1
                   {      
		     correspondingSurf[vox] = nsurf;
		     favoriteAxe[vox] = n;
		   }
	       }
	    }
	} // for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end(); ++q )
    } // while ( ! qbels.empty() )

  delete p;
}

/**
* test if the ball centered on c1.v of radius c1.dist is included 
* in the ball centered on C2.v of radius c2.dist
* @param ks a KnSpace
* @param c1 weighted cell representing a ball (voxel center, radius of the ball)
* @param c2 weighted cell representing a ball (voxel center, radius of the ball)
* @return true if the ball centered on c1 is included in the ball centered on c2
*/
bool included(const KnSpace & ks, weightedCell c1, weightedCell c2)
{
  double c1c2;
  Kn_size coord1 [3];  //to store the coordinates of the associated voxel
  Kn_size coord2 [3];
  ks.udecodeCoords(c1.v,coord1);
  ks.udecodeCoords(c2.v,coord2);
  c1c2 = sqrt( (coord2[0]-coord1[0])*(coord2[0]-coord1[0]) + 
               (coord2[1]-coord1[1])*(coord2[1]-coord1[1]) +
	       (coord2[2]-coord1[2])*(coord2[2]-coord1[2]) );
  return (c1c2+c1.dist <= c2.dist);
}

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


//main function
int main(int argc, char* argv[])
{
  StandardArguments::addDigitalArgs( args, 3, false, false );
  ShapeHelper::addSimple3DShapesArgs( args );
  args.addOption( "-noise", "-noise <level>: adds noise to the shape.", "0");
  args.addOption("-file","-file <fileName>: writes the points of the mutisacle profile in the file fileName","none.txt");
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
  Kn_size sizes[ Dim ]; 
  StandardArguments::fillSizes( args, sizes );
  KnSpace ks( Dim, sizes );
  //voxel type in the KnSpace
  Kn_uid typeVox = ks.uspel( ks.ufirst() );
  
  //create the shape  
  KnCharSet shape = ShapeHelper::makeSimple3DShapesFromArgs(args, ks);
  // Extracts one spel of the shape
  Kn_uid anyVox = shape.ubegin ();
  // Extracts one bel on the shape boundary.
  Kn_sid b= KnShapes::sfindFurthestBel( ks, anyVox, shape);
  
  //check if some noise must be added
  if(args.check("-noise"))
  {
    //DigitalSurfaceTracker initialization
    ObjectBoundaryTracker trackerInit;
    BelAdjacency baInit(ks,false);
    ObjectBoundary obInit(baInit,shape);
    trackerInit.init(&obInit);
    trackerInit.move(b);
    KnRCellSet bdryInit = KnShapes::strack(ks, trackerInit);
    
    Kn_uid in, out;
    vector<double> p_in( 6 );
    vector<double> p_out( 6 );
    int noise = args.getOption("-noise")->getIntValue( 0 );
    switch(noise){
      case 1 : { cout<<"noise level 1"<<endl;
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
      }break;
      case 2 : { cout<<"noise level 2"<<endl;
               p_in[0] = 0.0;
               p_in[1] = 0.5;
               p_in[2] = 0.15;
	       
               p_out[0] = 0.0;
               p_out[1] = 0.3;
               p_out[2] = 0.1;
               p_out[3] = 0.05;
      }break;
      default : { cout<<"noise level 0"<<endl;
               p_in[0] = 0.0;
               p_in[1] = 0.0;
               p_in[2] = 0.0;
               p_in[3] = 0.00;
               p_in[4] = 0.0;
               p_in[5] = 0.0;

               p_out[0] = 0.0;
               p_out[1] = 0.0;
               p_out[2] = 0.0;
               p_out[3] = 0.0;
               p_out[4] = 0.0;
               p_out[5] = 0.0;
      }
    }
  KnCharSet shapeWithNoise = KnShapes::noisifyObject(ks, shape, bdryInit, p_in, p_out, in, out); 
  cout << " --- noisy shape has " << shapeWithNoise.nbElements() << " voxels.--" << endl;
  KnCharSet main_inner_comp = KnCharSet::create (ks,true,3,0);
  b = ShapeHelper::findInnerObject(&ks, shapeWithNoise, in, main_inner_comp);
  shape = main_inner_comp;
  } //if(args.check("-noise"))
  
  Kn_size coord[ Dim ]; //to store coordinates 
  Kn_uid vox; //to represent a voxel
  S set [100000];  //to store each digital plane (set of voxels of type S -> see Z3.h)
  S otherSet[100000];
  int size = 0; //to store the size of each digital plane

  vector<double> scale(10); //different scales
  scale[0] = 1; scale[1] = 2; scale[2] = 3; scale[3] = 4; scale[4] = 5;
  scale[5] = 6; scale[6] = 7; scale[7] = 8; scale[8] = 9; scale[9] = 10;
  int indWidth = 0; //index of one scale
  float width = scale[indWidth]; //current scale
  
  //different maps for the statistic
  mapVoxEpStat greedyAreaMap; // map(voxel -> map(ep ->stats(aeras of greedy planes))
  mapVoxEpStat ballAreaMap; // map(voxel -> map(ep ->stats(aeras of isotropic planes))
  mapEpStat littleMap;  //tmp map(ep ->stats)
  FStatistic fstat;
  mapVoxSurf correspondingSurf; // map(voxel -> favorite surfel)
  mapVoxInt favoriteAxe; // map(voxel -> favorite axe)
 
  //DigitalSurfaceTracker initialization
  ObjectBoundaryTracker tracker;
  BelAdjacency ba(ks,false);
  ObjectBoundary ob(ba,shape);
  tracker.init(&ob);
  tracker.move(b);  //the tracker is on a surfel on the boundary of the main inner component
  DigitalSurfaceTracker* p = tracker.clone();

  // // WITH NU-THICK DISKS.
  // KnRCellSet bdry = KnShapes::strack( ks, tracker );
  // DigitalSurfaceGraph dsg( ks, bdry, tracker );
  
  /*tracking to find the favorites axes*/
  /*we obtain a map voxel -> surfel*/
  computeFavoriteAxes(ks, tracker, correspondingSurf, favoriteAxe, shape);
  cout<<"Favorite axes computed"<<endl;

  int axe; //major axe for each recognition
  vector<double> normalVector(3); //normal vector of the recognized piece of plane
  double radius; //radius of the maximal isotropic ball
  double area; //area of the recognized piece of plane
  
  while(indWidth<3) //for each scale
  {    
    cout<<"Scale : "<<width<<endl;
    //compute the maximal isotropic ball for each voxel in the map correspondingSurf
    mapVoxSurf::iterator itOnMap;
    itOnMap=correspondingSurf.begin();
    //priority queue to store the balls, a ball : a weighted cell (Kn-uid center, radius)
    //the largest on the top
    priority_queue<weightedCell,vector<weightedCell>,mycomparisonCell> fp; 
    while(itOnMap!=correspondingSurf.end())
    {
      // // WITH NU-THICK DISKS.
      // NuThickDisk<DigitalSurfaceGraph> disk;
      // disk.init( dsg, itOnMap->second );
      // disk.compute( (double) width );
      // disk.getNormal( normalVector );
      // checkSign(ks,itOnMap->second,normalVector);
      // // favoriteAxe[itOnMap->first] = disk.majorAxis();
      // // fp.push( weightedCell(itOnMap->first,  disk.radius() ) );
      // cerr << "Cell=" << itOnMap->first << " r=" << disk.radius() 
      // 	   << " #=" << disk.size() 
      // 	   << " axis=" << disk.majorAxis() 
      // 	   << " N=(" << normalVector[ 0 ] 
      // 	   << "," << normalVector[ 1 ] 
      // 	   << "," << normalVector[ 2 ] << ")" << endl;
      // // itOnMap++;

       // get the favorite axe
       axe = favoriteAxe[itOnMap->first];
       // compute the maximal disk
       // cout << "[ " << itOnMap->second;
       KnRCellSet surface = planeRecognitionOnSurfaces::maximalIsotropicPlane3( ks,tracker, itOnMap->second, axe, width,set,size, normalVector,radius);
       // cout << " r1=" << radius << " #=" << size;
       //checks the sign of the normal vector, changes it if necessary
       checkSign(ks,itOnMap->second,normalVector);
       //try with other axes
       // an other axe can become the new favorite one if the digital disk is bigger
       for(int i=1; i<=2; i++)
       {
          int otherAxe = (axe+i)%3;
       	  double otherRadius;
       	  int otherSize;
       	  vector<double> otherNormalVector(3);
          KnRCellSet otherSurface = planeRecognitionOnSurfaces::maximalIsotropicPlane3( ks,tracker, itOnMap->second, otherAxe, width, otherSet, otherSize, otherNormalVector,otherRadius);
       	  // cout << " ro=" << otherRadius  << " #=" << otherSize;
       	  checkSign(ks,itOnMap->second,otherNormalVector);
          if(otherRadius > radius)
          {
             axe = otherAxe;
             radius = otherRadius;
       	     size = otherSize;
             for(int j=0;j<otherSize;j++)
       	       set[j] = otherSet[j];
       	     favoriteAxe[itOnMap->first] = axe;
          }
       }
      // cerr << "Cell=" << itOnMap->first << " r=" << radius
      // 	   << " #=" << size
      // 	   << " axis=" << axe
      // 	   << " N=(" << normalVector[ 0 ] 
      // 	   << "," << normalVector[ 1 ] 
      // 	   << "," << normalVector[ 2 ] << ")" << endl;

       // add the ball in the priority queue
       // cout << endl;
       fp.push( weightedCell(itOnMap->first,radius) );
       itOnMap++;
    } //while(itOnMap!=correspondingSurf.end())
    cout<<"Maximal disks are computed"<<endl;
    
    //extract the classes representatives (just keep one representative of each class)
    list<weightedCell> f; //need a simple list
    //for each disk
    while(!fp.empty())
    {
      weightedCell wv = fp.top();
      fp.pop();
      list<weightedCell>::iterator itOnList;
      itOnList=f.begin();
      // check if the extracted disk is included in a larger one
      while(itOnList!=f.end())
      {
        //the extracted disk is not maximal, removes it from the main map correspondingSurf
        if(included(ks, wv,*itOnList)) 
	{
	  correspondingSurf.erase(wv.v); //delete the corresponding voxel of the map correspondingSurf
          break;
	}
	itOnList++;
      } //while(itOnList!=f.end())
      //the extracted disk is maximal
      if(itOnList==f.end())
        f.push_back(wv); //new class
    } //while(!fp.empty())
    cout<<"Class representatives are extracted"<<endl;
    
    //re-compute maximal planes by only for the representaives of the classes
    itOnMap=correspondingSurf.begin();
    while(itOnMap!=correspondingSurf.end())
    {
       axe = favoriteAxe[itOnMap->first];
       //recompute the maximal disk
       KnRCellSet ballSurface = planeRecognitionOnSurfaces::maximalIsotropicPlane2( ks,tracker, itOnMap->second, axe, width, set, size, normalVector,radius);
       checkSign(ks,itOnMap->second,normalVector);       
       area = computeArea(ks, ballSurface, normalVector);
 
/*       cout << " --- ball of area : " << area << endl;
       cout<<"nb surfels : "<<ballSurface.nbElements()<<endl;
       cout<<"---------------------------------------------------------------------------"<<endl;    */   
       
       //for each voxel covered by a maximal disk, adds the area of the disk to the corresponding stat
       for(int i=0;i<size;i++)
       {
	 coord[0] = get_si(set[i].x); coord[1] = get_si(set[i].y); coord[2] = get_si(set[i].z);
          vox = ks.ucode( coord, typeVox );

          ballAreaMap.insert(make_pair(vox,littleMap));
          ballAreaMap[vox].insert(make_pair(width,fstat));
          ballAreaMap[vox][width].addValue(area);
       }
       
       //compute the maximal extension (greedy)
       KnRCellSet surface = planeRecognitionOnSurfaces::maximalGreedyPlane( ks,tracker, itOnMap->second, axe, width, set, size, normalVector);
       checkSign(ks,itOnMap->second,normalVector);
       area = computeArea(ks, surface, normalVector);
 
//        cout << " --- surface of area : " << area << endl;
//        cout<<"nb surfels : "<<surface.nbElements()<<endl;
//        cout<<"---------------------------------------------------------------------------"<<endl;
       
       //for each voxel covered by an extended plane, add the area of the plane to the corresponding stat
       for(int i=0;i<size;i++)
       {
	 coord[0] = get_si(set[i].x); coord[1] = get_si(set[i].y); coord[2] = get_si(set[i].z);
          vox = ks.ucode( coord, typeVox );

          greedyAreaMap.insert(make_pair(vox,littleMap));
          greedyAreaMap[vox].insert(make_pair(width,fstat));
          greedyAreaMap[vox][width].addValue(area);
       }
     
       itOnMap++;
    } //while(itOnMap!=correspondingSurf.end())
    cout<<"Greedy extensions are computed"<<endl;
    
    indWidth++;
    width = scale[indWidth];
  } //while(indWidth<10)

  //store the datas in files for multiscale profile  
  if(args.check("-file"))
  {
    char beginning[60];
    char nomFichier[80];
    int nbFichiers=10;
    char numero[20];
    
    Kn_size coord1 [3];  //to store the coordinates of the voxel
    
    strcpy(beginning,(args.getOption("-file")->getValue(0)).c_str());

    mapVoxEpStat::iterator itGlob;
    itGlob=greedyAreaMap.begin();
    
    for(int n=1;n<=nbFichiers && itGlob!=greedyAreaMap.end();n++)
    {
      strcpy(nomFichier,beginning);
      sprintf(numero,"-%d",n);
      strcat(nomFichier,numero);
      strcat(nomFichier,".txt");
    
      ofstream fichier(nomFichier);
      fichier<<"#scale nbCoveringMaxDisks minDiskArea meanDiskArea maxDiskArea varianceBallArea nbCoveringGreedyPlanes  minGreedyPlaneArea meanGreedyPlaneArea maxGreedyPlaneArea varianceGreedyArea"<<endl;
      ks.udecodeCoords(itGlob->first,coord1);
      fichier<<"# voxel "<<itGlob->first<<" correspondant à ("<<coord1[0]<<","<<coord1[1]<<","<<coord1[2]<<")"<<endl;
      
      for(vector<double>::iterator scaleIt= scale.begin();scaleIt!=scale.end();scaleIt++)
      {
        double w = *scaleIt;
        fichier<<w<<" ";
	fichier<<ballAreaMap[itGlob->first][w].samples()<<" ";
	fichier<<ballAreaMap[itGlob->first][w].min()<<" "<<ballAreaMap[itGlob->first][w].mean()<<" ";
	fichier<<ballAreaMap[itGlob->first][w].max()<<" "<<ballAreaMap[itGlob->first][w].variance()<<" ";
	
	fichier<<greedyAreaMap[itGlob->first][w].samples()<<" ";
	fichier<<greedyAreaMap[itGlob->first][w].min()<<" "<<greedyAreaMap[itGlob->first][w].mean()<<" ";
	fichier<<greedyAreaMap[itGlob->first][w].max()<<" "<<greedyAreaMap[itGlob->first][w].variance()<<endl;

      }
      
      fichier.close();

      for(int k=0;k<100;k++) itGlob++;
    } 
  }//if(args.check("-file"))

  delete p;
  return 0;
}
