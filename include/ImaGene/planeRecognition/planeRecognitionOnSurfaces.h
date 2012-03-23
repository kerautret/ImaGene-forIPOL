/** @file planeRecognitionOnSurfaces.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : planeRecognitionOnSurfaces.h
//
// Creation : 2010/11/26
//
// Version : 2010/11/26
//
// Author : EC
//
// Summary : Header file for file planeRecognitionOnSurfaces.cxx
//
// History :
//	2010/11/26
//
// Rcs Id :
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(planeRecognitionOnSurfaces_RECURSES)
#error Recursive header files inclusion detected in planeRecognitionOnSurfaces.h
#else // defined(KnShapes_RECURSES)
#define planeRecognitionOnSurfaces_RECURSES

#if !defined planeRecognitionOnSurfaces_h
#define planeRecognitionOnSurfaces_h

//////////////////////////////////////////////////////////////////////////////  
#include <iostream>
#include <vector>
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/KnCharSet.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/KnRCellVector.h"
#include "ImaGene/digitalnD/DigitalSurfaceTracker.h"  

#include <cmath>
#include <queue>
#include <set>
#include <list>
#include <algorithm>
#include "ImaGene/mathutils/Mathutils.h"
#include "ImaGene/digitalnD/KnSpaceScanner.h"
#include "ImaGene/digitalnD/KnSpaceCoordScanner.h"
#include "ImaGene/digitalnD/BelAdjacency.h"
#include "ImaGene/digitalnD/ObjectBoundary.h"
#include "ImaGene/digitalnD/SurfelNeighborhood.h"

#include "ImaGene/dgeometry2d/C4CGeometry.h"
//////////////////////////////////////////////////////////////////////////////

#include "ImaGene/planeRecognition/Z3.h"
#include "ImaGene/planeRecognition/utils.h"
#include "ImaGene/planeRecognition/COBA.h"

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

using namespace std;


namespace ImaGene {

/** 
 * Description of class 'planeRecognitionOnSurfaces' <p>
 * Aim: To provide functions that builds maximal planes on digital 3d surfaces.
 */
class planeRecognitionOnSurfaces
{

public:

  class weightedBel //a surfel weighted by its distance to the seed point
  {
    public:
      Kn_sid s;
      double dist;
    
      weightedBel(Kn_sid surf, double distance){s=surf; dist=distance;}
  };
  
  class mycomparison
  {
    public:
      bool operator() (const weightedBel & b1, const weightedBel & b2) const
      {
        return b2.dist < b1.dist;
      }
  }; 
  
  /**
  * Given a surfel, computes all the maximal segments (DSS) passing through b
  * selects the two most centered (one in each direction)
  * @param ks the digital space
  * @param b a surfel
  * @param obj a 3d digital object
  * @param centeredSegment1 to store the first centered DSS
  * @param centeredSegment2 to store the second centered DSS
  */  
  static void computeMaximalSegmentsIterator( const KnSpace & ks,
			   const Kn_sid & b, const KnCharSet & obj, C4CSegment & centeredSegment1, C4CSegment & centeredSegment2);
			   
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
  static KnRCellSet maximalIsotropicPlane2( const KnSpace & ks,
			   const DigitalSurfaceTracker & tracker, Kn_sid b, int axe, float width, S* set, int & size, vector<double> & normalVector, double & radius );

  /**
   * version with COBAPlaneRecognition.
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
  static KnRCellSet maximalIsotropicPlane3
  ( const KnSpace & ks,
    const DigitalSurfaceTracker & tracker, Kn_sid b, int axe, float width, S* set, int & size, vector<double> & normalVector, double & radius );

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
  * @param normalVector the normal vector of the computed plane
  * @return the computed plane as a surface (set of surfels)
  */ 
  static KnRCellSet maximalGreedyPlane( const KnSpace & ks,
			   const DigitalSurfaceTracker & tracker, Kn_sid b, int axe, float width, S* set, int & size, vector<double> & normalVector );
			   
  /**
   * Given a tracker on a surface, scan the surface and stores in the file fileName 
   * the coordinates of the associated voxels
   * @param ks the digital space
   * @param tracker a tracker on a surface
   * @param fileName the name (path) of the file to write in
   */
  static void printVoxelSetInFileTracker( const KnSpace & ks,
						  const DigitalSurfaceTracker & tracker,
						  char fileName[30] );

  /**
   * Given a set of voxels (type S -> see Z3.h), stores in the file fileName 
   * their coordiantes
   * @param ks the digital space
   * @param set the set of voxels (table of S)
   * @param size the number of voxels
   * @param fileName the name (path) of the file to write in
   */
  static void printVoxelSetInFile( const KnSpace & ks,
						  S * set, int size, 
						  char fileName[30] );
			   
 };
 
 }
                                                                        //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined planeRecognitionOnSurfaces_h

#undef planeRecognitionOnSurfaces_RECURSES
#endif // else defined(planeRecognitionOnSurfaces_RECURSES)
