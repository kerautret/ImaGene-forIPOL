///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : planeRecognitionOnSurfaces.cxx
//
// Creation : 2010/11/26
//
// Version : 2010/11/26
//
// Author : EC
//
// email : Emilie.Charrier@univ-savoie.fr
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
//	2010/11/26 : Mr ?Name? : ?What?
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <queue>
// #include <priority_queue>
#include <set>
#include <list>
#include <algorithm>
#include<fstream>
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/digitalnD/KnShapes.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/digitalnD/KnSpaceCoordScanner.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/SurfelNeighborhood.h"

#include "ImaGene/digitalnD/KnTypes.h"
#include "ImaGene/dgeometry2d/C4CGeometry.h"
#include "ImaGene/digitalnD/C4CIteratorOnBdry.h"
#include "ImaGene/digitalnD/Frame2D.h"
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/VectorUtils.h"

#include "ImaGene/planeRecognition/planeRecognitionOnSurfaces.h"
#include "ImaGene/arithmetic/IntegerComputer.h"
#include "ImaGene/arithmetic/COBAPlaneRecognition.h"

using namespace std;

/**
  * Given a surfel, computes all the maximal segments (DSS) passing through b
  * selects the two most centered (one in each direction)
  * @param ks the digital space
  * @param b a surfel
  * @param obj a 3d digital object
  * @param centeredSegment1 to store the first centered DSS
  * @param centeredSegment2 to store the second centered DSS
  * @param trackdir1 first tracking direction
  * @param trackdir2 second tracking direction
  */  
void 
ImaGene::planeRecognitionOnSurfaces::computeMaximalSegmentsIterator( const KnSpace & ks,
			   const Kn_sid & b, const KnCharSet & obj, C4CSegment & centeredSegment1, C4CSegment & centeredSegment2)
{
  C4CSegment m_segments1[20];
  C4CSegment m_segments2[20];
  int gap = 100000;
  int newGap;
  double slope = 1000000;
  double newSlope;
  uint j = 0;  //beginning index
  uint k;  //ending index
  uint m = 20;  //max number of segments
  bool result;
  uint track_dir, n;
  bool oneFound = false;
  
  //to store all the maximal segments
  C4CIterator* tabIteratorFwd[20];
  C4CIterator* tabIteratorBwd[20];
  
  uint track_dir1, track_dir2; //the 2 tracking directions
  
  // Loop in all tracking directions (in practice, only two directions per surfel)
  for ( KnSpace::dir_iterator it = ks.sbegin_dirs( b ); ! it.end(); ++it )
  {
    track_dir = *it;
    if(!oneFound)
      {track_dir1 = track_dir; oneFound = true;}
    else
      track_dir2 = track_dir;
  }
  
  //compute all the maximal segments in direction track_dir1
  BelAdjacency ba(ks,true);
  C4CIteratorOnBdry it1(ba,b,track_dir1,obj);
  result = C4CGeometry::maximalSegmentsWithMemory( it1,m_segments1,
                                                   tabIteratorFwd,tabIteratorBwd, 
						   j, k, m );
  //check if the table was big enough
  if(!result)
  {
    cout<<"maximal segments : tableau trop petit!!!"<<endl;
    exit(-1);
   }
   
  // select the most centered DSS
  // prefer a DSS such that slope <= 1
   for(int i=0;i<k;i++)
   {
      newGap =  abs(abs((m_segments1[i].c_n()).x()) + abs((m_segments1[i].c_n()).y()) - abs((m_segments1[i].cp_n()).x()) - abs((m_segments1[i].cp_n()).y()));
      newSlope = abs((float(m_segments1[i].a()) / (float)m_segments1[i].b()));
      if( gap > newGap || (gap == newGap && (newSlope <= 1)) )
      {
        gap = abs((m_segments1[i].c_n()).x()) + abs((m_segments1[i].c_n()).y()) - abs((m_segments1[i].cp_n()).x()) - abs((m_segments1[i].cp_n()).y());
	gap = abs(gap);
        centeredSegment1 = m_segments1[i];
      }
    }
    gap = 100000;

    j=0;
    //compute all the maximal segments in direction track_dir2
    C4CIteratorOnBdry it2(ba,b,track_dir2,obj);
    result = C4CGeometry::maximalSegmentsWithMemory( it2, m_segments1,tabIteratorFwd,tabIteratorBwd, j, k, m );
    if(!result)
    {
      cout<<"maximal segments : tableau trop petit!!!"<<endl;
      exit(-1);
    }
    // select the most centered DSS
    // prefer a DSS such that slope <= 1
    for(int i=0;i<k;i++)
    {
      newGap =  abs(abs((m_segments1[i].c_n()).x()) + abs((m_segments1[i].c_n()).y()) - abs((m_segments1[i].cp_n()).x()) - abs((m_segments1[i].cp_n()).y()));
      newSlope = abs((float(m_segments1[i].a()) / (float)m_segments1[i].b()));
      if( gap > newGap || (gap == newGap && (newSlope <= 1)) )
      {
        gap = abs((m_segments1[i].c_n()).x()) + abs((m_segments1[i].c_n()).y()) - abs((m_segments1[i].cp_n()).x()) - abs((m_segments1[i].cp_n()).y());
	gap = abs(gap);
        centeredSegment2 = m_segments1[i];
      }
    }
    
    for(int i=0;i<k;i++)
    {
      delete tabIteratorFwd[i];
      delete tabIteratorBwd[i];
    }
}

  /**
  * Given a tracker on a surface and an initial voxel b, compute the maximal 
  * isotropic plane (disk), for a given width, centered on b
  * @param ks the digital space
  * @param tracker a tracker on the surface
  * @param b the initial surfel
  * @param axe the chosen major axis for recognition
  * @param width the maximal allowed width for the digital plane
  * @param set the resulating set of voxels (type S -> see Z3.h)
  * @param size the number of voxels of the computed disk
  * @param normalVector the normal vector of the computed disk
  * @param radius the radius of the disk
  * @return the computed disk as a surface (set of surfels)
  */ 
ImaGene::KnRCellSet 
ImaGene::planeRecognitionOnSurfaces::maximalIsotropicPlane2( const KnSpace & ks,
			   const DigitalSurfaceTracker & tracker, Kn_sid b, int axe, float width, S* set, int & size, vector<double> & normalVector, double & radius )
{
  //direction of tracking
  uint track_dir;
  //max diameter of the digital plane
  int D = max(ks.size(0),max(ks.size(1),ks.size(2)));
  //boolean table to know wether a voxel have already been added to the plane
  bool presentInSet[ks.size(0)][ks.size(1)][ks.size(2)];
  for(int i=0;i<ks.size(0);i++)
    for(int j=0;j<ks.size(1);j++)
      for(int k=0;k<ks.size(2);k++)
        presentInSet[i][j][k] = false;
  // surfels being explored
  KnRCellSet explored = KnRCellSet::create( ks, ks.dim() - 1, true, 0 );
  // plane being extracted (as a surface)
  KnRCellSet bdry = KnRCellSet::create( ks, ks.dim() - 1, true, 0 );
  // clone of the tracker
  DigitalSurfaceTracker* p = tracker.clone();
  p->move( b );
  
  Kn_sid nsurf;   // neighboor bel

  Kn_size coord [3];  //to store the coordinates of the initial voxel
  Kn_size xyzv [3];  //to store the coordinates of a neighboor voxel
  
  bool failOnce = false; //did the reco algorithm fail once?
  // priority queue to store the neighboors
  // ordered relative to the dist of the neighboor to the initial voxel
  // the closest on the top
  priority_queue<weightedBel,vector<weightedBel>,mycomparison> balls;
  balls.push( weightedBel(b,0) );
  explored += b;
  //extract the coordinates of the initial voxel b in coord
  uint  n = ks.sorthDir(b);
  Kn_sid  v = ks.sincident(b,n,ks.sdirect(b,n));
  ks.sdecodeCoords(v,coord);
  
  //inits the recognition algorithm by adding the initial voxel
  int nbAdded=0;
  int ind;
  set[0] = S(coord[0],coord[1],coord[2]);
  COBA cobaAlgo(set, D, axe);
  switch(axe) {
    case 0 : {normalVector[0] = 1;
              normalVector[1] = get_d(cobaAlgo.N.y)/get_d( cobaAlgo.N.x );
              normalVector[2] = get_d( cobaAlgo.N.z )/get_d( cobaAlgo.N.x );} break;
    case 1 : {normalVector[0] = get_d( cobaAlgo.N.x )/get_d( cobaAlgo.N.y );
              normalVector[1] = 1;
              normalVector[2] = get_d( cobaAlgo.N.z )/get_d( cobaAlgo.N.y );} break;
    case 2 : {normalVector[0] = get_d( cobaAlgo.N.x )/get_d( cobaAlgo.N.z );
              normalVector[1] = get_d( cobaAlgo.N.y )/get_d( cobaAlgo.N.z );
              normalVector[2] = 1;} break;
    }
  ind = 1;
  presentInSet[coord[0]][coord[1]][coord[2]] = true;
  cobaAlgo.cardinality = ind;
  
  Kn_sid b_init = b;
  
  //the initial disk contains one single voxel and is of radius 0
  double currentradius = 0; //radius of the current disk
  double previousradius = 0; //radius of the preceeding disk
  //currentRing of neighboors
  KnRCellSet currentRing = KnRCellSet::create( ks, ks.dim() - 1, true, 0 ); 
//   cout<<"before main loop"<<endl;
  //main loop
  while( !failOnce && !balls.empty() )
  {
    //extract a surfel
    b = (balls.top()).s;
//     if(b_init==1660638)
//       cout<<b<<" extracted (iso)"<<endl;
    balls.pop();  
    n = ks.sorthDir(b); //normal of the surfel
    v = ks.sincident(b,n,ks.sdirect(b,n));  //voxel correspondant
    ks.sdecodeCoords(v,xyzv);  //recuperation des coordonnées réelles du voxel
     
    //check the associated voxel is not alread added to the plane
    if(!presentInSet[xyzv[0]][xyzv[1]][xyzv[2]])
    {
      // add the voxel to the set
      set[ind] = S(xyzv[0],xyzv[1],xyzv[2]);;
      presentInSet[xyzv[0]][xyzv[1]][xyzv[2]] = true;
      ind++;
      cobaAlgo.cardinality = ind;
      //and test if it still corresponds to a digital plane
      if(!cobaAlgo.runFloatWidth(axe,width,failOnce, !failOnce))
      {
        //digital plane recognition failed, remove the voxel
        failOnce = true;
        ind--;
        cobaAlgo.cardinality--;
        presentInSet[xyzv[0]][xyzv[1]][xyzv[2]] = false;
// 	if(b_init==1660638)
//           cout<<"echec"<<endl;
      }
      else
      {
        //digital plane recognition succeded
        currentRing += b;
        nbAdded++; //one more voxel
      }
    } //if(!presentInSet[xyzv[0]][xyzv[1]][xyzv[2]])
    else //the associated voxel has already been added
    {
      currentRing += b;
    }
    
    //surf added -> add its neighboors to the priority queue
    if(! failOnce)
    {
      p->move( b );
      // Loop in all tracking directions. (Generic tracking).
      for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end() && !failOnce; ++q )
	{
	  track_dir = *q;
	  nsurf = p->adjacent( track_dir, false );
	  // If the surfel exists and is not already explored
	  if ( ( p->lastCode() != 0 ) && ( ! explored[ nsurf ] ) )
	  {
	      //computes its distance to the initial voxel
	      n = ks.sorthDir(nsurf);
              v = ks.sincident(nsurf,n,ks.sdirect(nsurf,n));
	      ks.sdecodeCoords(v,xyzv);
              double dist = sqrt( (xyzv[0]-coord[0])*(xyzv[0]-coord[0]) + (xyzv[1]-coord[1])*(xyzv[1]-coord[1]) + (xyzv[2]-coord[2])*(xyzv[2]-coord[2]) );
	      //and adds it to the priority queue
              balls.push(weightedBel(nsurf,dist));
	      explored += nsurf;
	  }
	  
	  nsurf = p->adjacent( track_dir, true );
	  // If the surfel exists and is not already explored,
	  if ( ( p->lastCode() != 0 ) && ( ! explored[ nsurf ] ) && !failOnce )
	  {
	      //computes its distance to the initial voxel
	      n = ks.sorthDir(nsurf);
              v = ks.sincident(nsurf,n,ks.sdirect(nsurf,n));
	      ks.sdecodeCoords(v,xyzv);
              double dist = sqrt( (xyzv[0]-coord[0])*(xyzv[0]-coord[0]) + (xyzv[1]-coord[1])*(xyzv[1]-coord[1]) + (xyzv[2]-coord[2])*(xyzv[2]-coord[2]) );
              //and adds it to the priority queue
              balls.push(weightedBel(nsurf,dist));
	      explored += nsurf;
	  }
	} //for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end() && !failOnce; ++q )
   } //if(! failOnce)
 
  // if no fail and one ring is full
  // adds the ring to the surface
  // updates the normal vector
  // change the current radius
  if( (!failOnce && !balls.empty() && (balls.top()).dist > currentradius) ||
      (!failOnce && balls.empty() ) )
  {
    bdry += currentRing;
    switch(axe) {
    case 0 : {normalVector[0] = 1;
              normalVector[1] = get_d( cobaAlgo.N.y )/get_d( cobaAlgo.N.x );
              normalVector[2] = get_d( cobaAlgo.N.z )/get_d( cobaAlgo.N.x );} break;
    case 1 : {normalVector[0] = get_d( cobaAlgo.N.x )/get_d( cobaAlgo.N.y );
              normalVector[1] = 1;
              normalVector[2] = get_d( cobaAlgo.N.z )/get_d( cobaAlgo.N.y );} break;
    case 2 : {normalVector[0] = get_d( cobaAlgo.N.x )/get_d( cobaAlgo.N.z );
              normalVector[1] = get_d( cobaAlgo.N.y )/get_d( cobaAlgo.N.z );
              normalVector[2] = 1;} break;
    }
    
    if(!balls.empty())
    {
      previousradius = currentradius;
      currentradius = (balls.top()).dist;
      currentRing -= currentRing;
      nbAdded = 0;
    }
  }

  }
  
  // consider only the voxels of full rings
  if(failOnce)
  {
    size = ind - nbAdded;
  }
  
  if(!balls.empty())
    currentradius = previousradius;
    
  radius = currentradius;
  
  /*desalloc*/
  cobaAlgo.free(true);
  delete p;
  
  return bdry;
}


  /**
  * Given a tracker on a surface and an initial voxel b, compute the extension
  * of the maximal disk, for a given width, centered on b
  * First adds isotropicly the neighboors of b, when a voxels
  * cannot be added, the normal vector of the digtal plane is fixed
  * and we try to add any other neighboor
  * @param ks the digital space
  * @param tracker a tracker on the surface
  * @param b the initial surfel
  * @param axe the chosen major axis for recognition
  * @param width the maximal allowed width for the digital plane
  * @param set the resulating set of voxels (type S -> see Z3.h)
  * @param size the number of voxels of the computed plane
  * @param normalVector the normal vector of the computed plane (not used)
  * @return the computed plane as a surface (set of surfels)
  */
ImaGene::KnRCellSet 
ImaGene::planeRecognitionOnSurfaces::maximalGreedyPlane( const KnSpace & ks,
			   const DigitalSurfaceTracker & tracker, Kn_sid b, int axe, float width, S* set, int & size, vector<double> & normalVector )
{
   //direction of tracking
  uint track_dir;
  //max diameter of the digital plane
  int D = max(ks.size(0),max(ks.size(1),ks.size(2)));
  //boolean table to know wether a voxel have already been added to the plane
  bool presentInSet[ks.size(0)][ks.size(1)][ks.size(2)];
  for(int i=0;i<ks.size(0);i++)
    for(int j=0;j<ks.size(1);j++)
      for(int k=0;k<ks.size(2);k++)
        presentInSet[i][j][k] = false;
  // surfels being explored
  KnRCellSet explored = KnRCellSet::create( ks, ks.dim() - 1, true, 0 );
  // boundary being extracted (as a surface)
  KnRCellSet bdry = KnRCellSet::create( ks, ks.dim() - 1, true, 0 );
  // clone of the tracker
  DigitalSurfaceTracker* p = tracker.clone();
  p->move( b );
  
  Kn_sid nsurf;   // next tracked bel 
  
  Kn_size coord [3];  //to store the coordinates of the initial voxel
  Kn_size xyzv [3];  //to store the coordinates of a neighboor voxel
  
  bool failPrevious = false;
  bool allowNewNormal = true;
  bool additionalSurf = false;

  // priority queue to store the neighboors
  // ordered relative to the dist of the neighboor to the initial voxel
  // the closest on the top
  priority_queue<weightedBel,vector<weightedBel>,mycomparison> balls;
  balls.push( weightedBel(b,0) );
  explored += b;

  //extract the coordinates of the initial voxel b in coord
  uint  n = ks.sorthDir(b);
  Kn_sid  v = ks.sincident(b,n,ks.sdirect(b,n));
  ks.sdecodeCoords(v,coord);
  
  //inits the recognition algorithm by adding the initial voxel
  int nbAdded=0;
  int ind;
  set[0] = S(coord[0],coord[1],coord[2]);
  COBA cobaAlgo(set, D, axe); //creates and inits the search space

  ind = 1;
  presentInSet[coord[0]][coord[1]][coord[2]] = true;
  cobaAlgo.cardinality = 1;
  
  Kn_sid b_init = b;

  //main loop
  while( !balls.empty() )
  {
    //extract a surfel
    b = (balls.top()).s;
//     if(b_init==1660638)
//       cout<<b<<" extracted (greedy)"<<endl;
    balls.pop();   
    n = ks.sorthDir(b);
    v = ks.sincident(b,n,ks.sdirect(b,n));
    ks.sdecodeCoords(v,xyzv);
      
    //check the associated voxel is not alread added to the plane
    if(!presentInSet[xyzv[0]][xyzv[1]][xyzv[2]])
    {
      // add the voxel to the set
      set[ind] = S(xyzv[0],xyzv[1],xyzv[2]);;
      presentInSet[xyzv[0]][xyzv[1]][xyzv[2]] = true;
      ind++;
      cobaAlgo.cardinality = ind;
      additionalSurf = false;
      //and test if it still corresponds to a digital plane
      if(!cobaAlgo.runFloatWidth(axe,width,failPrevious, allowNewNormal))
      {
        //digital plane recognition failed, remove the voxel
	failPrevious = true;
	allowNewNormal = false;
        ind--;
        cobaAlgo.cardinality--;
        presentInSet[xyzv[0]][xyzv[1]][xyzv[2]] = false;
      }
      else
      {
        //digital plane recognition succeded
	failPrevious = false;
        bdry += b;
        nbAdded++; //one more voxel
      }
    } //if(!presentInSet[xyzv[0]][xyzv[1]][xyzv[2]])
    else //the associated voxel has already been added
    {
      additionalSurf = true;
      bdry += b;
    }

    //voxel added to the plane -> add the neighboors of its surfel to the priority queue
    //if the surfel is just an additional surfel (the correponding voxel has been
    // added to the plane before), we add its neighboors to the priority queue
    // because the value of failPrevious has nothing to do with the surfel
    if(!failPrevious || additionalSurf)
    {
      p->move( b );
      // Loop in all tracking directions. (Generic tracking).
      for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end() ; ++q )
	{
	  track_dir = *q;
	  nsurf = p->adjacent( track_dir, false );
	  // If the surfel exists and is not already explored
	  if ( ( p->lastCode() != 0 ) && ( ! explored[ nsurf ] ) )
	  {
	      //computes its distance to the initial voxel
	      n = ks.sorthDir(nsurf); //normal of the surfel
              v = ks.sincident(nsurf,n,ks.sdirect(nsurf,n));  //voxel correspondant
	      ks.sdecodeCoords(v,xyzv);  //recuperation des coordonnées réelles du voxel
              double dist = sqrt( (xyzv[0]-coord[0])*(xyzv[0]-coord[0]) + (xyzv[1]-coord[1])*(xyzv[1]-coord[1]) + (xyzv[2]-coord[2])*(xyzv[2]-coord[2]) );
	      //and adds it to the priority queue
              balls.push(weightedBel(nsurf,dist));
	      explored += nsurf;
	  }
	  
	  nsurf = p->adjacent( track_dir, true );
	  // If the surfel exists and is not already explored,
	  if ( ( p->lastCode() != 0 ) && ( ! explored[ nsurf ] ) )
	  {
	      //computes its distance to the initial voxel
	      n = ks.sorthDir(nsurf);
              v = ks.sincident(nsurf,n,ks.sdirect(nsurf,n));
	      ks.sdecodeCoords(v,xyzv);
              double dist = sqrt( (xyzv[0]-coord[0])*(xyzv[0]-coord[0]) + (xyzv[1]-coord[1])*(xyzv[1]-coord[1]) + (xyzv[2]-coord[2])*(xyzv[2]-coord[2]) );
	      //and adds it to the priority queue
              balls.push(weightedBel(nsurf,dist));
	      explored += nsurf;
	  }
	} //for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end() ; ++q )
   } //if(!failPrevious || additionalSurf)
   
 } //while( !balls.empty() )
  
  size = ind;

  /*desalloc*/
  cobaAlgo.free(true);

  delete p;

  return bdry;
}


/**
* Given a tracker on a surface, scan the surface and stores in the file fileName 
* the coordinates of the associated voxels
* @param ks the digital space
* @param tracker a tracker on a surface
* @param fileName the name (path) of the file to write in
*/
void 
ImaGene::planeRecognitionOnSurfaces::printVoxelSetInFileTracker( const KnSpace & ks,
						  const DigitalSurfaceTracker & tracker,
						  char fileName[30] )
{
   // boundary being extracted.
  KnRCellSet bdry = KnRCellSet::create( ks, ks.dim() - 1, true, 0 );
  DigitalSurfaceTracker* p = tracker.clone();
  Kn_sid b = p->current();       // current tracked bel 
  Kn_sid nsurf;   // next tracked bel 
  uint track_dir;  // current tracked direction (one spanned by [b]).
  
  bdry += b;

  queue<Kn_sid> qbels;
  qbels.push( b );
  
  Kn_size xyzv [3];  //to store the coordinates of the associated voxel
  uint n; //normal vector
  Kn_sid v;
  
  ofstream fichier(fileName);
  if(fichier)
    cout<<"ouverture réussie!"<<endl;
  else
    cout<<"echec ouverure!"<<endl;
  Kn_uid vox;
  
  // For all pending bels
  int k=1;
  while( ! qbels.empty() )
    {
      k++;
      b = qbels.front();
      qbels.pop();

      p->move( b );
      
      n = ks.sorthDir(b); //normal of the surfel
      v = ks.sincident(b,n,ks.sdirect(b,n));  //voxel correspondant
      ks.sdecodeCoords(v,xyzv);  //recuperation des coordonnées réelles du voxel
      fichier<<xyzv[0]<<" "<<xyzv[1]<<" "<<xyzv[2]<<endl;
  
  
      // Loop in all tracking directions. (Generic tracking).
      for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end(); ++q )
	{
	  track_dir = *q;
	  nsurf = p->adjacent( track_dir, false );
	  // If the surfel exists and is not not already in boundary,
	  // add it to boundary and to the bels left for tracking.
	  if ( ( p->lastCode() != 0 ) && ( ! bdry[ nsurf ] ) )
	    {
	      bdry += nsurf;
	      qbels.push( nsurf );
	    }
 	  nsurf = p->adjacent( track_dir, true );
 	  // If the surfel exists and is not not already in boundary,
 	  // add it to boundary and to the bels left for tracking.
 	  if ( ( p->lastCode() != 0 ) && ( ! bdry[ nsurf ] ) )
 	    {
 	      bdry += nsurf;
 	      qbels.push( nsurf );
 	    }
	} // for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end(); ++q )
    } // while ( ! qbels.empty() )

  delete p;
  fichier.close();
}

/**
* Given a set of voxels (type S -> see Z3.h), stores in the file fileName 
* their coordiantes
* @param ks the digital space
* @param set the set of voxels (table of S)
* @param size the number of voxels
* @param fileName the name (path) of the file to write in
*/
void 
ImaGene::planeRecognitionOnSurfaces::printVoxelSetInFile( const KnSpace & ks,
						  S * set, int size, 
						  char fileName[30] )
{  
  ofstream fichier(fileName);
  if(fichier)
    cout<<"ouverture réussie!"<<endl;
  else
    cout<<"echec ouverure!"<<endl;
    
  for(int i=0;i<size;i++)
     fichier<<set[i].x<<" "<<set[i].y<<" "<<set[i].z<<endl;
  
  fichier.close();
}



/**
 * COBAPlaneRecognition version.
  * Given a tracker on a surface and an initial voxel b, compute the maximal 
  * isotropic plane (disk), for a given width, centered on b
  * @param ks the digital space
  * @param tracker a tracker on the surface
  * @param b the initial surfel
  * @param axe the chosen major axis for recognition
  * @param width the maximal allowed width for the digital plane
  * @param set the resulating set of voxels (type S -> see Z3.h)
  * @param size the number of voxels of the computed disk
  * @param normalVector the normal vector of the computed disk
  * @param radius the radius of the disk
  * @return the computed disk as a surface (set of surfels)
  */ 
ImaGene::KnRCellSet 
ImaGene::planeRecognitionOnSurfaces::maximalIsotropicPlane3
( const KnSpace & ks,
  const DigitalSurfaceTracker & tracker, Kn_sid b, int axe, float width, 
  S* set, int & size, vector<double> & normalVector, double & radius )
{
  //direction of tracking
  uint track_dir;
  //max diameter of the digital plane
  int D = max(ks.size(0),max(ks.size(1),ks.size(2)));
  //boolean table to know wether a voxel have already been added to the plane
  bool presentInSet[ks.size(0)][ks.size(1)][ks.size(2)];
  for(int i=0;i<ks.size(0);i++)
    for(int j=0;j<ks.size(1);j++)
      for(int k=0;k<ks.size(2);k++)
        presentInSet[i][j][k] = false;
  // surfels being explored
  KnRCellSet explored = KnRCellSet::create( ks, ks.dim() - 1, true, 0 );
  // plane being extracted (as a surface)
  KnRCellSet bdry = KnRCellSet::create( ks, ks.dim() - 1, true, 0 );
  // clone of the tracker
  DigitalSurfaceTracker* p = tracker.clone();
  p->move( b );
  
  Kn_sid nsurf;   // neighboor bel

  Kn_size coord [3];  //to store the coordinates of the initial voxel
  Kn_size xyzv [3];  //to store the coordinates of a neighboor voxel
  
  bool failOnce = false; //did the reco algorithm fail once?
  // priority queue to store the neighboors
  // ordered relative to the dist of the neighboor to the initial voxel
  // the closest on the top
  priority_queue<weightedBel,vector<weightedBel>,mycomparison> balls;
  balls.push( weightedBel(b,0) );
  explored += b;
  //extract the coordinates of the initial voxel b in coord
  uint  n = ks.sorthDir(b);
  Kn_sid  v = ks.sincident(b,n,ks.sdirect(b,n));
  ks.sdecodeCoords(v,coord);
  // cout << "[v=(" << coord[ 0 ] << "," << coord[ 1 ] << "," << coord[ 2 ] 
  //      << ") ]";
  //inits the recognition algorithm by adding the initial voxel
  typedef COBAPlaneRecognition::Point3i Point3i;
  int nbAdded=0;
  Point3i pt;
  for ( int i = 0; i < 3; ++i ) pt.coords[ i ] = coord[ i ];
  COBAPlaneRecognition coba;
  coba.init( axe, D, (double) width, pt );
  coba.getNormal( normalVector );
  presentInSet[coord[0]][coord[1]][coord[2]] = true;
  
  //the initial disk contains one single voxel and is of radius 0
  double currentradius = 0.0; //radius of the current disk
  double previousradius = 0.0; //radius of the preceeding disk
  //currentRing of neighboors
  KnRCellSet currentRing = KnRCellSet::create( ks, ks.dim() - 1, true, 0 ); 
  
  //main loop
  COBAPlaneRecognition::State saveRingState;
  coba.getState( saveRingState );
  while( ! failOnce && ! balls.empty() )
    {
      //extract a surfel
      b = ( balls.top() ).s;
      balls.pop();  
      n = ks.sorthDir(b); //normal of the surfel
      v = ks.sincident(b,n,ks.sdirect(b,n));  //voxel correspondant
      ks.sdecodeCoords(v,xyzv);  //recuperation des coordonnées réelles du voxel
     
      //check the associated voxel is not alread added to the plane
      if( ! presentInSet[xyzv[0]][xyzv[1]][xyzv[2]] )
	{
	  // add the voxel to the set
	  for ( int i = 0; i < 3; ++i ) pt.coords[ i ] = xyzv[ i ];
	  presentInSet[xyzv[0]][xyzv[1]][xyzv[2]] = true;

	  //and test if it still corresponds to a digital plane
	  COBAPlaneRecognition::State saveState;
	  coba.getState( saveState );
	  // std::cerr << "(" << pt[ 0 ] << "," << pt[ 1 ] << "," << pt[ 2 ] << ")";
	  if ( ! coba.add( pt, true ) )
	    {
	      // std::cerr << "X";
	      //digital plane recognition failed, remove the voxel
	      failOnce = true;
	      presentInSet[xyzv[0]][xyzv[1]][xyzv[2]] = false;
	      coba.setState( saveState );
	    }
	  else
	    {
	      // std::cerr << "+";
	      //digital plane recognition succeded
	      currentRing += b;
	      nbAdded++; //one more voxel
	    }
	} //if(!presentInSet[xyzv[0]][xyzv[1]][xyzv[2]])
      else //the associated voxel has already been added
	{
	  currentRing += b;
	}
    
      //surf added -> add its neighboors to the priority queue
      if(! failOnce)
	{
	  p->move( b );
	  // Loop in all tracking directions. (Generic tracking).
	  for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); 
		! q.end(); 
		++q )
	    {
	      track_dir = *q;
	      nsurf = p->adjacent( track_dir, false );
	      // If the surfel exists and is not already explored
	      if ( ( p->lastCode() != 0 ) && ( ! explored[ nsurf ] ) )
		{
		  //computes its distance to the initial voxel
		  n = ks.sorthDir(nsurf);
		  v = ks.sincident(nsurf,n,ks.sdirect(nsurf,n));
		  ks.sdecodeCoords(v,xyzv);
		  double dist = sqrt( (xyzv[0]-coord[0])*(xyzv[0]-coord[0])
				      + (xyzv[1]-coord[1])*(xyzv[1]-coord[1])
				      + (xyzv[2]-coord[2])*(xyzv[2]-coord[2]) );
		  //and adds it to the priority queue
		  balls.push(weightedBel(nsurf,dist));
		  explored += nsurf;
		}
	      nsurf = p->adjacent( track_dir, true );
	      // If the surfel exists and is not already explored,
	      if ( ( p->lastCode() != 0 ) && ( ! explored[ nsurf ] ) )
		{
		  //computes its distance to the initial voxel
		  n = ks.sorthDir(nsurf);
		  v = ks.sincident(nsurf,n,ks.sdirect(nsurf,n));
		  ks.sdecodeCoords(v,xyzv);
		  double dist = sqrt( (xyzv[0]-coord[0])*(xyzv[0]-coord[0])
				      + (xyzv[1]-coord[1])*(xyzv[1]-coord[1])
				      + (xyzv[2]-coord[2])*(xyzv[2]-coord[2]) );
		  //and adds it to the priority queue
		  balls.push(weightedBel(nsurf,dist));
		  explored += nsurf;
		}
	    } //for ( KnSpace::dir_iterator q = ks.sbegin_dirs( b ); ! q.end(); ++q )
	} //if(! failOnce)
 
      // if no fail and one ring is full
      // adds the ring to the surface
      // updates the normal vector
      // change the current radius
      if( (!failOnce && !balls.empty() && (balls.top()).dist > currentradius) ||
	  (!failOnce && balls.empty() ) )
	{
	  bdry += currentRing;
	  coba.getNormal( normalVector );
	  coba.getState( saveRingState );
	  if(!balls.empty())
	    {
	      previousradius = currentradius;
	      currentradius = (balls.top()).dist;
	      currentRing -= currentRing;
	      nbAdded = 0;
	    }
	}
      
    } // while( ! failOnce && ! balls.empty() )
  
  // consider only the voxels of full rings
  coba.setState( saveRingState );
  size = saveRingState.nb_used;
  for ( int i = 0; i < size; ++i )
    {
      set[ i ].x = coba.point( i ).coords[ 0 ];
      set[ i ].y = coba.point( i ).coords[ 1 ];
      set[ i ].z = coba.point( i ).coords[ 2 ];
    }
  if(!balls.empty())
    currentradius = previousradius;
    
  radius = currentradius;

  // std::cerr << "=>" << currentradius << std::endl;
  delete p;
  
  return bdry;
}

