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
#include "ImaGene/planes/DigitalMidSurfaceGraph.h"
#include "ImaGene/planes/CompletedDigitalSurfaceGraph.h"
#include "ImaGene/planes/NuThickDisk.h"
#include "ImaGene/planes/TangentialCoverHierarchy.h"
#include "ImaGene/planes/TangentialCoverDecomposition.h"
#include "ImaGene/mathutils/Statistic.h"
#include "ImaGene/timetools/Clock.h"

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
  case 3 :  cout<<"noise level 3"<<endl;
    p_in[0] = 0.0;
    p_in[1] = 0.01;
    p_in[2] = 0.001;
    p_in[3] = 0.0;
    p_in[4] = 0.0;
    p_in[5] = 0.0;
    
    p_out[0] = 0.0;
    p_out[1] = 0.01;
    p_out[2] = 0.001;
    p_out[3] = 0.0;
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




typedef float Vec3f[ 3 ];

void outputCellInColor( ostream & out, 
			const KnSpace & ks, 
			Kn_sid s, const Vec3f & color )
{
  Kn_size x[ 3 ];
  Kn_sign sign = ks.decodeSign( s );
  ks.skdecodeCoords( s, x );
  out << "Cell " << x[ 0 ] << ' ' << x[ 1 ] << ' ' << x[ 2 ] << ' '
      << sign << ' ' 
      << color[ 0 ] << ' ' << color[ 1 ] << ' ' << color[ 2 ] << endl;

}
void outputCellInColorWithNormal( ostream & out, 
				  const KnSpace & ks, 
				  Kn_sid s, const Vec3f & color,
				  const Vec3f & n )
{
  Kn_size x[ 3 ];
  Kn_sign sign = ks.decodeSign( s );
  ks.skdecodeCoords( s, x );
  out << "CellN " << x[ 0 ] << ' ' << x[ 1 ] << ' ' << x[ 2 ] << ' '
      << sign 
      << ' ' << color[ 0 ] << ' ' << color[ 1 ] << ' ' << color[ 2 ]
      << ' ' << n[ 0 ] << ' ' << n[ 1 ] << ' ' << n[ 2 ]
      << endl;

}
void outputLinelsOfSurfelInColor( ostream & out, 
				  const KnSpace & ks, 
				  Kn_sid surfel, const Vec3f & color )
{
  for ( KnSpace::dir_iterator q = ks.sbegin_dirs( surfel );
	q != 0; ++q )
    {
      outputCellInColor( out, ks, ks.sincident( surfel, *q, true ), color );
      outputCellInColor( out, ks, ks.sincident( surfel, *q, false ), color );
    }
}
template <typename Iterator>
void outputSignedSetInColor( ostream & out, 
			     const KnSpace & ks, 
			     Iterator b, Iterator e, const Vec3f & color )
{
  for ( ; b != e; ++b )
    outputCellInColor( out, ks, *b, color );
}

void
fillRandomColors( unsigned int nb, Vec3f* colors )
{
  for ( unsigned int i = 0; i < nb; ++i )
    {
      colors[ i ][ 0 ] = ((double) random()) / (double) RAND_MAX;
      colors[ i ][ 1 ] = ((double) random()) / (double) RAND_MAX;
      colors[ i ][ 2 ] = ((double) random()) / (double) RAND_MAX;
    }
}

template <typename TDigitalGraph>
void
visualizeDecomposition
( ostream & out,
  const KnSpace & ks, 
  typename TDigitalGraph::Set bdry,
  const TangentialCoverDecomposition<TDigitalGraph> & tcd )
{
  typedef typename
    TangentialCoverDecomposition<TDigitalGraph>::VertexSet
    VertexSet;
  Vec3f grey = { 0.8, 0.8, 0.8 };
  Vec3f red = { 1.0, 0.0, 0.0 };
  outputSignedSetInColor( out, ks, bdry.begin(), bdry.end(), grey );
  const VertexSet & reps = tcd.allRepresentatives();
  Vec3f* colors = new Vec3f[ reps.size() ];
  fillRandomColors( reps.size(), colors );
  unsigned int i = 0;
  std::map< Kn_sid, unsigned int > vtx2color;
  for ( typename VertexSet::const_iterator it = reps.begin(), it_end = reps.end();
	it != it_end; ++it, ++i )
    {
      outputLinelsOfSurfelInColor( out, ks, *it, colors[ i ] );
      vtx2color[ *it ] = i;
    }
  for ( KnRCellSet::cell_iterator it = bdry.begin(), it_end = bdry.end();
	it != it_end;
	++it )
    {
      const VertexSet & vtx_reps = tcd.representatives( *it );
      Vec3f mcolor;
      mcolor[ 0 ] = mcolor[ 1 ] = mcolor[ 2 ] = 0.0;
      
      for ( typename VertexSet::const_iterator 
      	      itr = vtx_reps.begin(), itr_end = vtx_reps.end();
      	    itr != itr_end; ++itr )
      	for ( unsigned int i = 0; i < 3; ++i )
      	  mcolor[ i ] += colors[ vtx2color[ *itr ] ][ i ];
      if ( vtx_reps.size() > 0 )
	{
	  for ( unsigned int i = 0; i < 3; ++i )
	    mcolor[ i ] /= vtx_reps.size();
	  outputCellInColor( out, ks, *it, mcolor );
	}
      else
	cerr << "  -- Vertex " << *it << " has no representative." << endl;
      // outputCellInColor( out, ks, *it, 
      // 			 colors[ vtx2color[ *(vtx_reps.begin()) ] ] );
    }
    //outputSignedSetInColor( out, ks, reps.begin(), reps.end(), red );

  typedef typename
    TangentialCoverDecomposition<TDigitalGraph>::ComputedDisk
    ComputedDisk;
  typedef typename TDigitalGraph::Methods Methods;
  Methods methods( *( tcd.digitalGraph() ) );
  i = 0;
  for ( typename VertexSet::const_iterator 
	  it = reps.begin(), it_end = reps.end();
	it != it_end; ++it, ++i )
    {
      ComputedDisk disk;
      disk.init( *( tcd.digitalGraph() ), *it, tcd.scale() );
      disk.computeDisk( -1.0 );
      Vec3f n;
      disk.getNormal( n );
      int x[ 3 ]; 
      tcd.digitalGraph()->embed( *it, x );
      Vec3f xup, xlow;
      xup[ 0 ] = xlow[ 0 ] = x[ 0 ] + 0.5;
      xup[ 1 ] = xlow[ 1 ] = x[ 1 ] + 0.5;
      xup[ 2 ] = xlow[ 2 ] = x[ 2 ] + 0.5;
      xup[ methods.preferredAxis( *it ) ] += 0.5 * tcd.scale();
      xlow[ methods.preferredAxis( *it ) ] -= 0.5 * tcd.scale();
      out << "Disk " << xup[ 0 ] << ' ' << xup[ 1 ] << ' ' << xup[ 2 ] << ' '
	  << disk.radius() << ' ' 
	  << colors[ i ][ 0 ] << ' ' 
	  << colors[ i ][ 1 ] << ' ' 
	  << colors[ i ][ 2 ] << ' ' 
	  << "0.5 " << n[ 0 ] << ' ' << n[ 1 ] << ' ' << n[ 2 ] << ' '
	  << endl;
      out << "Disk " << xlow[ 0 ] << ' ' << xlow[ 1 ] << ' ' << xlow[ 2 ] << ' '
	  << disk.radius() << ' ' 
	  << colors[ i ][ 0 ] << ' ' 
	  << colors[ i ][ 1 ] << ' ' 
	  << colors[ i ][ 2 ] << ' ' 
	  << "0.5 " << n[ 0 ] << ' ' << n[ 1 ] << ' ' << n[ 2 ] << ' '
	  << endl;
    }
  delete colors;
}

template <typename TDigitalGraph>
void
visualizeHierarchy
( ostream & out,
  const KnSpace & ks, 
  typename TDigitalGraph::Set bdry,
  const TangentialCoverHierarchy<TDigitalGraph> & tch )
{
  typedef 
    TangentialCoverDecomposition<CompletedDigitalSurfaceGraph>::VertexSet
    VertexSet;
  Vec3f grey = { 0.3, 0.3, 0.3 };
  Vec3f colors[ 8 ] = { 
    { 0.0, 0.0, 0.0 },
    { 0.0, 0.0, 1.0 },
    { 0.0, 1.0, 0.0 },
    { 0.0, 1.0, 1.0 },
    { 1.0, 0.0, 0.0 },
    { 1.0, 0.0, 1.0 },
    { 1.0, 1.0, 0.0 },
    { 1.0, 1.0, 1.0 },
  };
  outputSignedSetInColor( out, ks, bdry.begin(), bdry.end(), grey );
  for ( unsigned int i = 1; i < 8; ++i )
    {
      typename TDigitalGraph::Set active_set = tch.activeVertices( (double) i );
      outputSignedSetInColor( out, ks, 
			      active_set.begin(), active_set.end(), 
			      colors[ i ] );
    }
    exit(-1);
}

template <typename TDigitalGraph>
void
visualizeTangentialCover
( ostream & out,
  const KnSpace & ks, 
  typename TDigitalGraph::Set bdry,
  const TangentialCover<TDigitalGraph> & tc )
{
  typedef typename TangentialCover<TDigitalGraph>::MaximalPlaneList
    MaximalPlaneList;
  typedef typename TangentialCover<TDigitalGraph>::MaximalPlaneIndexSet
    MaximalPlaneIndexSet;
  typedef typename TangentialCover<TDigitalGraph>::MaximalPlaneSummary
    MaximalPlaneSummary;
  typedef typename MaximalPlaneList::iterator MaximalPlaneListIterator;
  Vec3f grey = { 0.8, 0.8, 0.8 };
  Vec3f red = { 1.0, 0.0, 0.0 };
  outputSignedSetInColor( out, ks, bdry.begin(), bdry.end(), grey );
  const MaximalPlaneList & all_mps = tc.myAllPlanes;
  Vec3f* colors = new Vec3f[ all_mps.size() ];
  fillRandomColors( all_mps.size(), colors );
  unsigned int i = 0;
  std::map< Kn_sid, unsigned int > vtx2color;
  for ( typename MaximalPlaneList::const_iterator 
	  it = all_mps.begin(), it_end = all_mps.end();
	it != it_end; ++it, ++i )
    {
      const MaximalPlaneSummary & mp = *it;
      outputLinelsOfSurfelInColor( out, ks, mp.center, colors[ i ] );
      vtx2color[ mp.center ] = i;
    }
  for ( KnRCellSet::cell_iterator it = bdry.begin(), it_end = bdry.end();
	it != it_end;
	++it )
    {
      typename TangentialCover<TDigitalGraph>::MapVertexToMaximalPlanes::const_iterator search_mp = tc.myMapVtx2MP.find( *it );
      if ( search_mp == tc.myMapVtx2MP.end() )
	std::cerr << "Error at vtx=" << *it << endl;
      const MaximalPlaneIndexSet & vtx_reps = search_mp->second;
      Vec3f mcolor;
      mcolor[ 0 ] = mcolor[ 1 ] = mcolor[ 2 ] = 0.0;
      
      // std::cerr << "  mp vtx=" << *it;
      unsigned int nb = 0;
      for ( typename MaximalPlaneIndexSet::const_iterator 
      	      itr = vtx_reps.begin(), itr_end = vtx_reps.end();
      	    itr != itr_end; ++itr )
	{
	  unsigned int index = *itr;
	  // std::cerr << " " << index;
	  if ( index < tc.myAllPlanes.size() )
	    {
	      ++nb;
	      const MaximalPlaneSummary & mp = tc.myAllPlanes[ index ];
	      for ( unsigned int i = 0; i < 3; ++i )
		mcolor[ i ] += colors[ vtx2color[ mp.center ] ][ i ];
	    }
	}
      //std::cerr << std::endl;
      if ( nb > 0 )
	{
	  for ( unsigned int i = 0; i < 3; ++i )
	    mcolor[ i ] /= nb;
	  // if ( view_normals )
	  //   {
	  //     Vec3f n;
	  //     tc.getEstimatedNormal( n, *it );
	  //     outputCellInColorWithNormal( out, ks, *it, mcolor, n );
	  //   }
	  outputCellInColor( out, ks, *it, mcolor );
	}
      else
	cerr << "  -- Vertex " << *it << " has no representative." << endl;
    }

  typedef typename TangentialCover<TDigitalGraph>::ComputedDisk
    ComputedDisk;
  typedef typename TDigitalGraph::Methods Methods;
  Methods methods( *( tc.digitalGraph() ) );
  i = 0;
  for ( typename MaximalPlaneList::const_iterator 
	  it = all_mps.begin(), it_end = all_mps.end();
	it != it_end; ++it, ++i )
    {
      const MaximalPlaneSummary & mp = *it;
      double xup[ 3 ];
      double xlow[ 3 ];
      mp.projectUpperPlane( xup, mp.center_embedding );
      mp.projectLowerPlane( xlow, mp.center_embedding );
      xup[ 0 ] += 0.5; xup[ 1 ] += 0.5; xup[ 2 ] += 0.5;
      xlow[ 0 ] += 0.5; xlow[ 1 ] += 0.5; xlow[ 2 ] += 0.5;
      out << "Disk " 
	  << xup[ 0 ] << ' ' << xup[ 1 ] << ' ' << xup[ 2 ] << ' '
	  << mp.radius << ' ' 
	  << colors[ i ][ 0 ] << ' ' 
	  << colors[ i ][ 1 ] << ' ' 
	  << colors[ i ][ 2 ] << ' ' 
	  << "0.5 " 
	  << mp.normal[ 0 ] << ' ' << mp.normal[ 1 ] << ' ' 
	  << mp.normal[ 2 ] << ' '
	  << endl;
      out << "Disk " 
	  << xlow[ 0 ] << ' ' << xlow[ 1 ] << ' ' << xlow[ 2 ] << ' '
	  << mp.radius << ' ' 
	  << colors[ i ][ 0 ] << ' ' 
	  << colors[ i ][ 1 ] << ' ' 
	  << colors[ i ][ 2 ] << ' ' 
	  << "0.5 "
 	  << mp.normal[ 0 ] << ' ' << mp.normal[ 1 ] << ' ' 
	  << mp.normal[ 2 ] << ' '
	  << endl;
    }
  delete colors;
}

template <typename TDigitalGraph>
void
visualizeTangentialCoverNormals
( ostream & out,
  const KnSpace & ks, 
  typename TDigitalGraph::Set bdry,
  TangentialCover<TDigitalGraph> & tc )
{
  Vec3f grey = { 0.8, 0.8, 0.8 };
  Vec3f red = { 1.0, 0.0, 0.0 };
  Vec3f n;
  //outputSignedSetInColor( out, ks, bdry.begin(), bdry.end(), grey );
  typedef 
    typename TangentialCover<TDigitalGraph>::AveragingMode AveragingMode;
  AveragingMode nd = 
    args.getOption( "-normals" )->getValue( 0 ) == "SA" ? tc.SimpleAveraging :
    args.getOption( "-normals" )->getValue( 0 ) == "DA" ? tc.DistanceAveraging :
    args.getOption( "-normals" )->getValue( 0 ) == "RDA" ? tc.RadiusAndDistanceAveraging :
    args.getOption( "-normals" )->getValue( 0 ) == "IOA" ? tc.InOutAveraging :
    args.getOption( "-normals" )->getValue( 0 ) == "MPD" ? tc.MaxProjectedDisk :
    tc.MaxProjectedPlane;

  cerr << "--- Displaying normals with "
       << ( nd == tc.SimpleAveraging ? "SimpleAveraging" :
	    nd == tc.DistanceAveraging ? "DistanceAveraging" :
	    nd == tc.RadiusAndDistanceAveraging ? "RadiusAndDistanceAveraging" :
	    nd == tc.InOutAveraging ? "InOutAveraging" :
	    nd == tc.MaxProjectedDisk ? "MaxProjectedDisk" :
	    nd == tc.MaxProjectedPlane ? "MaxProjectedPlane" :
	    "Unknown method" ) << endl;
  
  double edge_angle = args.getOption( "-edge_angle" )->getDoubleValue( 0 );
  Statistic<double> angle_stat;
  for ( KnRCellSet::cell_iterator it = bdry.begin(), it_end = bdry.end();
	it != it_end;
	++it )
    {
      std::vector<double> coefs;
      typename TangentialCover<TDigitalGraph>::Vertex p = *it;
      tc.getAveragingCoefficients( coefs, p, nd );
      tc.getEstimatedNormal( n, p, coefs );
      tc.getNormalAngleStatistic( angle_stat, p, n );
      Vec3f & color = 
	( ( angle_stat.samples() == 0 )
	  || ( sqrt( angle_stat.mean() ) > edge_angle ) ) ? red : grey;
      outputCellInColorWithNormal( out, ks, p, color, n );
    }
}

template <typename TDigitalGraph>
void
visualizeNoiseLevel
( ostream & out,
  const KnSpace & ks, 
  typename TDigitalGraph::Set bdry,
  const TangentialCoverHierarchy<TDigitalGraph> & tch,
  unsigned int min_width, double max_slope,
  bool mp_area )
{
  Vec3f noise_colors[ 10 ] = 
    { { 0.0, 1.0, 1.0 },
      { 0.8, 0.8, 0.8 },
      { 0.8, 0.8, 0.4 },
      { 1.0, 1.0, 0.0 },
      { 1.0, 0.7, 0.0 },
      { 1.0, 0.4, 0.0 },
      { 1.0, 0.0, 0.0 },
      { 1.0, 0.0, 0.4 },
      { 1.0, 0.0, 0.7 },
      { 1.0, 0.0, 1.0 }
    };
  Vec3f n;
  // outputSignedSetInColor( out, ks, bdry.begin(), bdry.end(), grey );
  for ( KnRCellSet::cell_iterator it = bdry.begin(), it_end = bdry.end();
	it != it_end;
	++it )
    {
      ScaleProfile sp;
      tch.getScaleProfile( sp, *it, mp_area );
      unsigned int nlvl = sp.noiseLevel( min_width, max_slope );
      // std::cerr << " -- " << *it << " nlvl=" << nlvl;
      std::vector<double> x, y;
      sp.getProfile( x, y );
      // for ( unsigned int j = 0; j < x.size(); ++j )
      // 	std::cerr << "(" << x[ j ] << " " << y[ j ] << ")";
      // std::cerr << std::endl;
      outputCellInColor( out, ks, *it, noise_colors[ nlvl ] );
    }
}


template <typename TDigitalGraph>
void
visualizeMaximalPlaneCores
( ostream & out,
  const KnSpace & ks, 
  typename TDigitalGraph::Set bdry,
  const TangentialCover<TDigitalGraph> & tc )
{
  typedef typename TangentialCover<TDigitalGraph>::MaximalPlaneList
    MaximalPlaneList;
  typedef typename TangentialCover<TDigitalGraph>::MaximalPlaneIndexSet
    MaximalPlaneIndexSet;
  typedef typename TangentialCover<TDigitalGraph>::MaximalPlaneSummary
    MaximalPlaneSummary;
  typedef typename MaximalPlaneList::iterator MaximalPlaneListIterator;
  Vec3f grey = { 0.8, 0.8, 0.8 };
  Vec3f red = { 1.0, 0.0, 0.0 };
  outputSignedSetInColor( out, ks, bdry.begin(), bdry.end(), grey );
  const MaximalPlaneList & all_mps = tc.myAllPlanes;
  Vec3f* colors = new Vec3f[ all_mps.size() ];
  fillRandomColors( all_mps.size(), colors );
  unsigned int i = 0;
  std::map< Kn_sid, unsigned int > vtx2color;
  for ( typename MaximalPlaneList::const_iterator 
	  it = all_mps.begin(), it_end = all_mps.end();
	it != it_end; ++it, ++i )
    {
      const MaximalPlaneSummary & mp = *it;
      outputLinelsOfSurfelInColor( out, ks, mp.center, colors[ i ] );
      vtx2color[ mp.center ] = i;
    }
  for ( KnRCellSet::cell_iterator it = bdry.begin(), it_end = bdry.end();
	it != it_end;
	++it )
    {
      unsigned int coremp_idx = tc.getCoreMaximalPlaneIndex( *it );
      const MaximalPlaneSummary & mp = tc.myAllPlanes[ coremp_idx ];
      outputCellInColor( out, ks, *it, colors[ vtx2color[ mp.center ] ] );
    }
  delete colors;
}



  
void computeTgtCoverHierarchy
( ostream & out, const KnSpace & ks, KnCharSet voxset, 
  vector<double> & scales, bool ball_inclusion, bool active_vtx )
{
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
  DigitalSurfaceGraphMethods methods( dsg );
  TangentialCoverHierarchy<DigitalSurfaceGraph> tch;
  tch.init( dsg, scales );
  tch.compute( methods, ball_inclusion );
  TangentialCoverDecomposition<DigitalSurfaceGraph> tcd;
  tcd.init( dsg, scales.back(), 
	    tch.activeVertices( (int) (scales.size() - 1) ) );
  // tcd.computeAccretion( methods );
  tcd.compute( methods, tcd.ExactMaximalPlanes );
  if ( active_vtx )
    visualizeHierarchy( out, ks, bdry, tch );
  else
    visualizeDecomposition( out, ks, bdry, tcd );
}
  

void computeTgtCoverDecomposition
( ostream & out,
  const KnSpace & ks, KnCharSet voxset, double scale )
{
  // boundary is needed to get the set of all vertices.
  KnRCellSet bdry = KnShapes::smakeBoundary( ks, voxset );
  // the tracker is used for extracting neighbors within the graph of
  // bel-adjacency.
  ObjectBoundaryTracker tracker;
  BelAdjacency ba( ks, false );
  ObjectBoundary ob( ba, voxset );
  tracker.init( &ob );
  tracker.move( *( bdry.begin() ) );
  CompletedDigitalSurfaceGraph dsg( ks, bdry, tracker );
  CompletedDigitalSurfaceGraphMethods methods( dsg );
  TangentialCoverDecomposition<CompletedDigitalSurfaceGraph> tcd;
  tcd.init( dsg, scale, bdry );
  tcd.computeAccretion( methods );
  // tcd.computeSegmentation( methods );
  // tcd.compute( methods );

  visualizeDecomposition( out, ks, bdry, tcd );
}

enum TCHDisplayMode {
  TCHNone,
  TCHDecomposition,
  TCHNormals,
  TCHDecompositionAndNormals,
  TCHActiveVertices,
  TCHCores,
  TCHNoiseLevel
};

void computeTgtCoverHierarchy2
( ostream & out, const KnSpace & ks, KnCharSet voxset, 
  vector<double> & scales, 
  TCHDisplayMode mode,
  unsigned int min_size = 1, double max_slope = -0.2,
  bool mp_area = true )
{
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
  DigitalSurfaceGraphMethods methods( dsg );
  TangentialCoverHierarchy<DigitalSurfaceGraph> tch;
  tch.init( dsg, scales );
  bool not_in_core = args.getOption( "-mp" )->getValue( 0 ) == "NIC";
  tch.computeTangentialCovers( methods, not_in_core );

  cerr << "---- Recompute common normals..." << flush;
  TangentialCover<DigitalSurfaceGraph> & tc = tch.tangentialCover( scales.size() - 1 );
  // tc.intersectPlanesAccordingToEmbedding();
  cerr << " ended." << endl;

  if ( ( mode == TCHDecomposition ) 
       || ( mode == TCHDecompositionAndNormals ) )
    visualizeTangentialCover( out, ks, bdry, 
			      tch.tangentialCover( scales.size() - 1 ) );
  if ( ( mode == TCHNormals ) 
       || ( mode == TCHDecompositionAndNormals ) )
    visualizeTangentialCoverNormals( out, ks, bdry, 
				     tch.tangentialCover( scales.size() - 1 ) );
  if ( mode == TCHNoiseLevel ) 
    visualizeNoiseLevel( out, ks, bdry, tch, min_size, max_slope, mp_area );
  if ( mode == TCHCores ) 
    visualizeMaximalPlaneCores( out, ks, bdry,
				tch.tangentialCover( scales.size() - 1 ) );
}

void computeTgtCoverHierarchy3
( ostream & out, const KnSpace & ks, KnCharSet voxset, 
  vector<double> & scales, 
  TCHDisplayMode mode,
  unsigned int min_size = 1, double max_slope = -0.2,
  bool mp_area = true )
{
  // boundary is needed to get the set of all vertices.
  KnRCellSet bdry = KnShapes::smakeBoundary( ks, voxset );
  // the tracker is used for extracting neighbors within the graph of
  // bel-adjacency.
  ObjectBoundaryTracker tracker;
  BelAdjacency ba( ks, false );
  ObjectBoundary ob( ba, voxset );
  tracker.init( &ob );
  tracker.move( *( bdry.begin() ) );
  DigitalMidSurfaceGraph dsg( ks, bdry, tracker );
  DigitalMidSurfaceGraphMethods methods( dsg );
  TangentialCoverHierarchy<DigitalMidSurfaceGraph> tch;
  tch.init( dsg, scales );
  bool not_in_core = args.getOption( "-mp" )->getValue( 0 ) == "NIC";
  tch.computeTangentialCovers( methods, not_in_core );
  
  if ( ( mode == TCHDecomposition ) 
       || ( mode == TCHDecompositionAndNormals ) )
    visualizeTangentialCover( out, ks, bdry, 
			      tch.tangentialCover( scales.size() - 1 ) );
  if ( ( mode == TCHNormals ) 
       || ( mode == TCHDecompositionAndNormals ) )
    visualizeTangentialCoverNormals( out, ks, bdry, 
				     tch.tangentialCover( scales.size() - 1 ) );
  if ( mode == TCHNoiseLevel ) 
    visualizeNoiseLevel( out, ks, bdry, tch, min_size, max_slope, mp_area );
  if ( mode == TCHCores ) 
    visualizeMaximalPlaneCores( out, ks, bdry,
				tch.tangentialCover( scales.size() - 1 ) );
}


//main function
int main(int argc, char* argv[])
{
  StandardArguments::addDigitalArgs( args, 3, false, false );
  ShapeHelper::addSimple3DShapesArgs( args );
  args.addOption( "-voxset", "-voxset <filename_pgm3d> <min> <max>: the input shape is read from a pgm3d file. If I is the voxel value, the voxel is added to the shape whenever min <= I <= max.", "", "0", "128" );
  args.addOption( "-noise", "-noise <level>: adds noise to the shape.", "0");
  args.addOption( "-scales", "-scales <scales>: compute hierarchical tangential cover at the specified scales. Default is 1,2,3,4,5,6,7,8.", "1,2,3,4,5,6,7,8");
  args.addOption( "-file","-file <fileName>: writes the points of the mutisacle profile in the file fileName","none.txt");
  args.addOption( "-hierarchy","-hierarchy <STD|BALL|EXT>: Computes the hierarchy and chooses either ball inclusion or extension inclusion for the hierarchy.","STD");
  args.addOption( "-decompose","-decompose <scale>: Decomposes the shape into maximal planes.", "1.0");
  args.addOption( "-view","-view <D|A|N|DN|L|C>: view either (D)ecomposition or (A)ctive vertices, or (N) normals, or both decomposition and normals (DN), or (L) noise level, or (C) cores of maximal planes.", "A");

  args.addOption( "-meaningfulScales", "-meaningfulScales <min_size> <max_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum slopes between consecutive samples within.", "1", "-0.2" );  
  args.addOption( "-normals", "-normals <SA|DA|RDA|IOA>: specifies the method used for generating normals starting from a set of maximal planes: SA: simple avering, DA: averaging weighted by the distance to center of the plane, RDA: mix of raidus and distance averaging, IOA: Inside radius weight is 1, outside weight decreases with the distance.", "SA" );  
  args.addOption( "-mp","-mp <NIC|NIMP>: defines maximal planes either as Not In the Core (smoothest, like maximal segments) or Not In Maximal Plane (coarsest, more greedy segmentation).", "NIC");
  args.addOption( "-surface","-surface <DSG|DMSG>: defines the kind of digital surface: DSG digital surface graph, DMSG digital mid surface graph.", "DSG");
  args.addOption( "-edge_angle","-edge_angle <angle>: indicates the angle between the normals and the estimated normals  starting from which the vertex is considered an edge.", ".78539816339744830961" );
  args.addOption( "-check_disk","-check_disk <nb>: checks disk computation <nb> times.", "10" );
  args.addOption( "-area","-area <MP|DISK>: specifies how areas of maximal planes are computed: MP means the whole maximal plane, DISK means only the disk area .", "MP" );
  if ( ! args.readArguments( argc, argv ) ) 
    {
      cout << args.usage( "tcover3d", 
			  "Computes the tangential cover of a 3D object.",
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
      cerr << "[tcover3d] Unable to create space or shape." << endl;
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
      cerr << "[tcover3d] Shape has zero elements." << endl;
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
  
  // computeTgtCover3D( *ks, *voxset, scales );
  ofstream output( args.getOption( "-file" )->getValue( 0 ).c_str(), ios::out );
  if ( args.check( "-hierarchy" ) )
    {
      bool ball_inclusion = 
	args.getOption( "-hierarchy" )->getValue( 0 ) == "BALL";
      string view = args.getOption( "-view" )->getValue( 0 );
      TCHDisplayMode mode = TCHNone;
      bool active_vtx = false;
      if ( view == "A" ) 
	{ 
	  mode = TCHActiveVertices;
	  active_vtx = true;
	}
      else if ( view == "D" )
	mode = TCHDecomposition;
      else if ( view == "DN" )
	mode = TCHDecompositionAndNormals;
      else if ( view == "N" )
	mode = TCHNormals;
      else if ( view == "L" )
	mode = TCHNoiseLevel;
      else if ( view == "C" )
	mode = TCHCores;

      uint mscales_min_size = 
	args.getOption( "-meaningfulScales" )->getIntValue( 0 );
      double mscales_max_slope = 
	args.getOption( "-meaningfulScales" )->getDoubleValue( 1 );
      bool mp_area = args.getOption( "-area" )->getValue( 0 ) == "MP";


      if ( args.getOption( "-hierarchy" )->getValue( 0 ) != "STD" )
	computeTgtCoverHierarchy( output, *ks, *voxset, scales, 
				  ball_inclusion, active_vtx );
      else if ( args.getOption( "-surface" )->getValue( 0 ) == "DSG" )
	computeTgtCoverHierarchy2( output, *ks, *voxset, scales, 
				   mode, mscales_min_size, mscales_max_slope,
                                   mp_area );
      else
	computeTgtCoverHierarchy3( output, *ks, *voxset, scales, 
				   mode, mscales_min_size, mscales_max_slope,
                                   mp_area );
    }
  if ( args.check( "-decompose" ) )
    {
      double scale = args.getOption( "-decompose" )->getDoubleValue( 0 );
      computeTgtCoverDecomposition( output, *ks, *voxset, scale );
    }
  output.close();

  if ( args.check( "-check_disk" ) )
    {
      unsigned int nb = args.getOption( "-check_disk" )->getIntValue( 0 );
      // boundary is needed to get the set of all vertices.
      KnRCellSet bdry = KnShapes::smakeBoundary( *ks, *voxset );
      // the tracker is used for extracting neighbors within the graph of
      // bel-adjacency.
      ObjectBoundaryTracker tracker;
      BelAdjacency ba( *ks, false );
      ObjectBoundary ob( ba, *voxset );
      tracker.init( &ob );
      tracker.move( *( bdry.begin() ) );
      DigitalSurfaceGraph dsg( *ks, bdry, tracker );
      DigitalSurfaceGraphMethods methods( dsg );
      NuThickDisk<DigitalSurfaceGraph> disk;
      Kn_sid vtx = *( bdry.begin()  );
      for ( unsigned int i = 0; i < nb; ++i )
	{
	  disk.init( dsg, vtx, scales[ 0 ] );
	  disk.computeDisk( -1.0 );
	  std::cerr << disk << std::endl;
	}
    }
  delete voxset;
  delete ks;
  return 0;
}
