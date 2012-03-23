//////////////////////////////////////////////////////////////////////////////
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
  
  args.addOption("-importPGM3D", "-importPGM3D <fileName> <threshold Min> <threshold max> import pgm3D ","pgmFile.p3D",  "1", "128");
  args.addOption("-exportNoiseLevels" , "-exportNoiseLevels <maxSlope> <minSize> <fileName.p3d", "0", "1", "noiseLevel.p3d");
  args.addOption("-file","-file <fileName>: writes the points of the mutisacle profile in the file fileName","none.txt");
  args.addOption("-maxScale", "-maxScale <maxScale>", "10");

  if ( ! args.readArguments( argc, argv ) || !args.check("-importPGM3D") ) 
    {
      cout << args.usage( "exampleExportNoiseLevels args", 
			  "-importPGM3D <fileName> <thresholdMin> <thresholdMax>",
			  "-importPGM3D -file -maxScale" )
	   << endl;
      return 1;
    }

  
  // Import 3D shape
  const uint Dim = 3; 
  Kn_size sizes[ Dim ]; 
  KnSpace *ks ;
  KnCharSet* shape;
  string fileNameImport= args.getOption("-importPGM3D")->getValue(0); 
  uint thresholdMin = args.getOption("-importPGM3D")->getIntValue(1);
  uint thresholdMax = args.getOption("-importPGM3D")->getIntValue(2);
  
  uint maxScale = args.getOption("-maxScale")->getIntValue(0) ;
  
  
  ifstream ifs;
  ifs.open(fileNameImport.c_str(), ifstream::in);
  ShapeHelper::importFromPGM3d(ifs, ks, shape, thresholdMin, thresholdMax,5);    
  ifs.close();
  Kn_uid anyVox = shape->ubegin ();
  Kn_sid b= KnShapes::sfindFurthestBel( *ks, anyVox, *shape);
  Kn_uid typeVox = ks->uspel( ks->ufirst() );

  cout << " --- noisy shape has " << shape->nbElements() << " voxels.--" << endl;


  
  Kn_size coord[ Dim ]; //to store coordinates 
  Kn_uid vox; //to represent a voxel
  S set [10000];  //to store each digital plane (set of voxels of type S -> see Z3.h)
  S otherSet[10000];
  int size = 0; //to store the size of each digital plane

  vector<double> scale(maxScale); //different scales
  
  for(int i=0; i< scale.size();i++){
    scale[i]= i+1;
  }
  

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
  BelAdjacency ba(*ks,false);
  ObjectBoundary ob(ba,*shape);
  tracker.init(&ob);
  tracker.move(b);  //the tracker is on a surfel on the boundary of the main inner component
  DigitalSurfaceTracker* p = tracker.clone();
  
  /*tracking to find the favorites axes*/
  /*we obtain a map voxel -> surfel*/
  computeFavoriteAxes(*ks, tracker, correspondingSurf, favoriteAxe, *shape);
  cout<<"Favorite axes computed"<<endl;

  int axe; //major axe for each recognition
  vector<double> normalVector(3); //normal vector of the recognized piece of plane
  double radius; //radius of the maximal isotropic ball
  double area; //area of the recognized piece of plane
  
  while(indWidth<maxScale) //for each scale
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
       // get the favorite axe
       axe = favoriteAxe[itOnMap->first];
       // compute the maximal disk
       KnRCellSet surface = planeRecognitionOnSurfaces::maximalIsotropicPlane2( *ks,tracker, itOnMap->second, axe, width,set,size, normalVector,radius);
       //checks the sign of the normal vector, changes it if necessary
       checkSign(*ks,itOnMap->second,normalVector);
       //try with other axes
       // an other axe can become the new favorite one if the digital disk is bigger
       for(int i=1; i<=2; i++)
       {
          int otherAxe = (axe+i)%3;
	  double otherRadius;
	  int otherSize;
	  vector<double> otherNormalVector(3);
          KnRCellSet otherSurface = planeRecognitionOnSurfaces::maximalIsotropicPlane2( *ks,tracker, itOnMap->second, otherAxe, width, otherSet, otherSize, otherNormalVector,otherRadius);
	  checkSign(*ks,itOnMap->second,otherNormalVector);
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
       //add the ball in the priority queue
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
        if(included(*ks, wv,*itOnList)) 
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
       KnRCellSet ballSurface = planeRecognitionOnSurfaces::maximalIsotropicPlane2( *ks,tracker, itOnMap->second, axe, width, set, size, normalVector,radius);
       checkSign(*ks,itOnMap->second,normalVector);       
       area = computeArea(*ks, ballSurface, normalVector);
 
/*       cout << " --- ball of area : " << area << endl;
       cout<<"nb surfels : "<<ballSurface.nbElements()<<endl;
       cout<<"---------------------------------------------------------------------------"<<endl;    */   
       
       //for each voxel covered by a maximal disk, adds the area of the disk to the corresponding stat
       for(int i=0;i<size;i++)
       {
	 coord[0] = get_si(set[i].x); coord[1] = get_si(set[i].y); coord[2] = get_si(set[i].z);
          vox = ks->ucode( coord, typeVox );

          ballAreaMap.insert(make_pair(vox,littleMap));
          ballAreaMap[vox].insert(make_pair(width,fstat));
          ballAreaMap[vox][width].addValue(area);
       }
       
       //compute the maximal extension (greedy)
       KnRCellSet surface = planeRecognitionOnSurfaces::maximalGreedyPlane( *ks,tracker, itOnMap->second, axe, width, set, size, normalVector);
       checkSign(*ks,itOnMap->second,normalVector);
       area = computeArea(*ks, surface, normalVector);
 
//        cout << " --- surface of area : " << area << endl;
//        cout<<"nb surfels : "<<surface.nbElements()<<endl;
//        cout<<"---------------------------------------------------------------------------"<<endl;
       
       //for each voxel covered by an extended plane, add the area of the plane to the corresponding stat
       for(int i=0;i<size;i++)
       {
	 coord[0] = get_si(set[i].x); coord[1] = get_si(set[i].y); coord[2] = get_si(set[i].z);
          vox = ks->ucode( coord, typeVox );

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
  
  //Computing the noise levels.
  
  if(args.check("-exportNoiseLevels")){
    double maxSlope = args.getOption("-exportNoiseLevels")->getFloatValue(0);
    uint minSize = args.getOption("-exportNoiseLevels")->getIntValue(1);
    string fileNameExportMax = args.getOption("-exportNoiseLevels")->getValue(2)+"Max.p3d";
    string fileNameExportMean = args.getOption("-exportNoiseLevels")->getValue(2)+"Mean.p3d";
    string fileNameExportGreedyMax = args.getOption("-exportNoiseLevels")->getValue(2)+"GreedyMax.p3d";
    string fileNameExportGreedyMean = args.getOption("-exportNoiseLevels")->getValue(2)+"GreedyMean.p3d";

    KnRUCellVector<int> noiseVoxelsToExportMax (*ks, 3);
    KnRUCellVector<int> noiseVoxelsToExportMean (*ks, 3);
    KnRUCellVector<int> noiseVoxelsToExportGreedyMax (*ks, 3);
    KnRUCellVector<int> noiseVoxelsToExportGreedyMean (*ks, 3);
    
    
    mapVoxEpStat::iterator itGlob;
    itGlob=greedyAreaMap.begin();

    while(itGlob!=greedyAreaMap.end()){
      ScaleProfile profileMean;
      ScaleProfile profileGreedyMean;
      ScaleProfile profileMax;
      ScaleProfile profileGreedyMax;
      
      profileMean.init(scale.size());
      profileMax.init(scale.size());
      profileGreedyMean.init(scale.size());
      profileGreedyMax.init(scale.size());
      
      
      uint ss=0;
      for(vector<double>::iterator scaleIt= scale.begin();scaleIt!=scale.end();scaleIt++){
	double valMean  = ballAreaMap[itGlob->first][*scaleIt].mean()/((ss+1.0)*(ss+1.0));
	profileMean.addValue(ss, valMean );
	double valMax  = ballAreaMap[itGlob->first][*scaleIt].max()/((ss+1.0)*(ss+1.0));
	profileMax.addValue(ss, valMax );
	double valMeanGreedy  = greedyAreaMap[itGlob->first][*scaleIt].mean()/((ss+1.0)*(ss+1.0));
	profileGreedyMean.addValue(ss, valMeanGreedy );
	double valMaxGreedy  = greedyAreaMap[itGlob->first][*scaleIt].max()/((ss+1.0)*(ss+1.0));
	profileGreedyMax.addValue(ss, valMaxGreedy );
	ss++;
      }
      uint noiseLevelMean= profileMean.noiseLevel(minSize, maxSlope);
      uint noiseLevelMax= profileMax.noiseLevel(minSize, maxSlope);
      uint noiseLevelGreedyMean= profileGreedyMean.noiseLevel(minSize, maxSlope);
      uint noiseLevelGreedyMax= profileGreedyMax.noiseLevel(minSize, maxSlope);
      
      noiseVoxelsToExportMax[itGlob->first]=noiseLevelMax;
      noiseVoxelsToExportMean[itGlob->first]=noiseLevelMean;
      noiseVoxelsToExportGreedyMean[itGlob->first]=noiseLevelGreedyMean;
      noiseVoxelsToExportGreedyMax[itGlob->first]=noiseLevelGreedyMax;
      
      
      itGlob++;
    } 
    ofstream ofsMax(fileNameExportMax.c_str());
    ofstream ofsMean(fileNameExportMean.c_str());
    ofstream ofsGreedyMax(fileNameExportGreedyMax.c_str());
    ofstream ofsGreedyMean(fileNameExportGreedyMean.c_str());

    ShapeHelper::exportToPGM3d(ofsMax, ks, noiseVoxelsToExportMax);
    ShapeHelper::exportToPGM3d(ofsMean, ks, noiseVoxelsToExportMean);
    ShapeHelper::exportToPGM3d(ofsGreedyMax, ks, noiseVoxelsToExportGreedyMax);
    ShapeHelper::exportToPGM3d(ofsGreedyMean, ks, noiseVoxelsToExportGreedyMean);
    ofsMax.close();
    ofsMean.close();
    ofsGreedyMax.close();
    ofsGreedyMean.close();
    
  }
  
  return 0;
}



