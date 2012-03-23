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
#include "ImaGene/image/Image2D.h"
#include "ImaGene/image/PGMFilter.h"
#include "ImaGene/planes/DigitalImageGraph.h"
#include "ImaGene/planes/NuThickDisk.h"
#include "ImaGene/planes/TangentialCoverHierarchy.h"
#include "ImaGene/mathutils/Statistic.h"
#include "ImaGene/timetools/Clock.h"

using namespace std;
using namespace ImaGene;

static Arguments args;

// Types for tangential cover.
typedef Image2D<unsigned char> GreyLevelImage2D;
typedef DigitalImageGraph<GreyLevelImage2D> MyDigitalImageGraph;
typedef TangentialCoverHierarchy< MyDigitalImageGraph > MyTangentialCoverHierarchy;
typedef DigitalImageGraphMethods<GreyLevelImage2D> MyDigitalImageGraphMethods;
  


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

struct ScaleProfileParameters
{
  uint min_width;
  double max_slope;
  double min_slope;
  double lb_scale_1;
  double lb_slope;
  uint nsamples;
  double alpha;  
};

void
displayScaleProfileWithGnuplot( const MyTangentialCoverHierarchy & tch, 
				unsigned int x, unsigned int y,
				const ScaleProfileParameters & params )
{
  ofstream f_cmd, f_xy;
  f_cmd.open( "./gnuplot.tmp" );
  f_xy.open( "./gnuplot.xy" );
  if ( f_cmd.good() && f_xy.good() )
    {
      Pixel p( x, y );
      ScaleProfile sp_disk, sp_mp;
      tch.getScaleProfile( sp_disk, p, false ); // Disk area
      tch.getScaleProfile( sp_mp, p, true );   // MP area
      unsigned int nlvl = 
	sp_disk.noiseLevel( params.min_width, params.max_slope );
      unsigned int lb_nlvl = 
	sp_disk.lowerBoundedNoiseLevel( params.min_width, params.max_slope, 
					params.min_slope,
					log(params.lb_scale_1), 
					params.lb_slope );
      unsigned int ss_nlvl = 
	sp_disk.standardScale( params.nsamples, params.alpha );

      std::cout << "# " << p << " nlvl=" << nlvl
		<< " lb_nlvl=" << lb_nlvl
		<< " ss_nlvl=" << ss_nlvl << std::endl;
      std::vector<double> sx_disk, sy_disk, sx_mp, sy_mp;
      sp_disk.getProfile( sx_disk, sy_disk );
      sp_mp.getProfile( sx_mp, sy_mp );
      for ( unsigned int j = 0; j < sx_disk.size(); ++j )
       	f_xy << exp( sx_disk[ j ] ) << " " << exp( sy_disk[ j ] ) << " "
	     << exp( sx_mp[ j ] ) << " " << exp( sy_mp[ j ] ) << std::endl;
      unsigned int l = sx_disk.size() - 1;
      f_cmd << "set logscale xy" << endl;
      f_cmd << "plot \"./gnuplot.xy\" using 1:2 title 'Profile Disk' w l"
	    << ", \"./gnuplot.xy\" using 3:4 title 'Profile MP' w l" 
	    << ", (x/" << exp(sx_disk[ l ]) << ")**(-1.0)*" << exp(sy_disk[ l ]) 
	    << " title 'pente -1'" 
	    << ", (x/" << exp(sx_disk[ l ]) << ")**(-2.0)*" << exp(sy_disk[ l ])
	    << " title 'pente -2'" 
	    << ", " << params.lb_scale_1 << "*x**" 
	    << params.lb_slope 
	    << " title 'lower bound'" 
	    << endl;
      system( "gnuplot -persist gnuplot.tmp" );
    }
}

  

//main function
int main(int argc, char* argv[])
{
  // StandardArguments::addDigitalArgs( args, 3, false, false );
  // ShapeHelper::addSimple3DShapesArgs( args );
  args.addOption( "-image", "-image <filename_pgm>: the input image.", "essai.pgm" );
  args.addOption( "-outImage", "-outImage <filename_pgm>: the output image. If not given, uses the standard output stream.", "-" );
  args.addOption( "-scales", "-scales <scales>: compute hierarchical tangential cover at the specified scales. Default is 1,2,3,4,5,6,7,8.", "1,2,3,4,5,6,7,8");
  args.addOption( "-mp","-mp <NIC|NIMP>: defines maximal planes either as Not In the Core (smoothest, like maximal segments) or Not In Maximal Plane (coarsest, more greedy segmentation).", "NIC");
  args.addOption( "-meaningfulScales", "-meaningfulScales <min_size> <max_slope>: specifies parameters for defining meaningful scales: minimum size of the interval of scales and maximum slopes between consecutive samples within.", "1", "-0.2" );  
  args.addOption( "-standardScale", "-standardScale <nsamples> <alpha>: choose the standard scale as noise level estimator.", "3", "0.05" );
  args.addOption( "-lowerBound", "-lowerBound <val_scale_1> <slope>: imposes an affine lower bound in the profile when extracting the noise level. To be used in conjunction with meaningfulScales.", "1", "-3.0" );
  args.addOption( "-area","-area <MP|DISK>: specifies how areas of maximal planes are computed: MP means the whole maximal plane, DISK means only the disk area .", "MP" );
  args.addOption( "-embedding","-embedding <qd> <xf> <yf>: defines how pixels are embedding in Z3. <qd> quantification (>=1), <xf> x factor for the embedding (>=1), <yf> y factor for the embedding (>=1).", "1", "1", "1" );
  args.addBooleanOption( "-rescale","-rescale: linearly rescales the noise level to be in 0-255." );
  args.addBooleanOption( "-interactive","-interactive: launch interactive mode. The program waits for lines x y on the input stream and displays with gnuplot scale profiles." );
  if ( ! args.readArguments( argc, argv ) ) 
    {
      cout << args.usage( "tcover-image", 
			  "Reads som input image and compute its tangential cover.",
			  "" )
	   << endl;
      return 1;
    }

  PGMFilter<GreyLevelImage2D> filter;
  GreyLevelImage2D image;
  string filename = args.getOption( "-image" )->getValue( 0 );
  ifstream input( filename.c_str() );
  if ( ! filter.read( image, input ) )
    {
      cerr << "Error reading pgm file: " 
	   << filename << endl;
      return 1;
    }
  
  // Getting parameters.
  vector<double> scales;
  string2doubles( scales, args.getOption( "-scales" )->getValue( 0 ), ',' );
  uint min_width = 
    args.getOption( "-meaningfulScales" )->getIntValue( 0 );
  double max_slope = 
    args.getOption( "-meaningfulScales" )->getDoubleValue( 1 );
  bool lb_nl        = args.check( "-lowerBound" );
  double lb_scale_1 = args.getOption( "-lowerBound" )->getDoubleValue( 0 );
  double lb_slope   = args.getOption( "-lowerBound" )->getDoubleValue( 1 );

  bool mp_area = args.getOption( "-area" )->getValue( 0 ) == "MP";
  unsigned int quantification = 
    args.getOption( "-embedding" )->getIntValue( 0 );
  unsigned int xfactor = 
    args.getOption( "-embedding" )->getIntValue( 1 );
  unsigned int yfactor = 
    args.getOption( "-embedding" )->getIntValue( 2 );
  bool std_scale = args.check( "-standardScale" );
  uint nsamples = 
    args.getOption( "-standardScale" )->getIntValue( 0 );
  double alpha = 
    args.getOption( "-standardScale" )->getDoubleValue( 1 );

  MyTangentialCoverHierarchy tch;
  MyDigitalImageGraph dig( image, quantification, xfactor, yfactor );
  MyDigitalImageGraphMethods methods( dig );
  tch.init( dig, scales );
  bool not_in_core = args.getOption( "-mp" )->getValue( 0 ) == "NIC";
  tch.computeTangentialCovers( methods, not_in_core );

  GreyLevelImage2D noiseImage( image.width(), image.height() );

  unsigned int max = 0;
  for ( unsigned int y = 0; y < image.height(); ++y )
    for ( unsigned int x = 0; x < image.width(); ++x )
    {
      Pixel pixel( x, y );
      ScaleProfile sp;
      tch.getScaleProfile( sp, pixel, mp_area );

      unsigned int nlvl = std_scale
	? sp.standardScale( nsamples, alpha )
	: ( lb_nl 
	    ? sp.lowerBoundedNoiseLevel( min_width, max_slope, -1e10,
					 log(lb_scale_1), lb_slope )
	    : sp.noiseLevel( min_width, max_slope ) );
      std::cerr << " -- " << pixel << " nlvl=" << nlvl;
      std::vector<double> sx, sy;
      sp.getProfile( sx, sy );
      for ( unsigned int j = 0; j < sx.size(); ++j )
       	std::cerr << "(" << sx[ j ] << " " << sy[ j ] << ")";
      std::cerr << std::endl;
      if ( nlvl > max ) max = nlvl;
      noiseImage.set( x, y, nlvl );
    }
  if ( args.check( "-rescale" ) && ( max != 0 ) )
    {
      for ( unsigned int y = 0; y < image.height(); ++y )
	for ( unsigned int x = 0; x < image.width(); ++x )
	  {
	    unsigned int val = ( noiseImage.at( x, y ) * 255 ) / max;
	    noiseImage.set( x, y, (unsigned char) val );
	    std::cerr << " " << val;
	  }
      std::cerr << std::endl;
    }

  // outputs image.
  string fname = args.getOption( "-outImage" )->getValue( 0 );
  std::ofstream out_file;
  if ( fname != "-" )
    out_file.open( fname.c_str() );
  std::ostream & out_image_stream = ( fname == "-" ) ? std::cout : out_file;
  bool error;
  if ( error = ! filter.write( noiseImage, out_image_stream ) )
    std::cerr << "[tcover-image] Error writing image." << std::endl;

  // Checks interactive mode
  if ( args.check( "-interactive" ) )
    {
      ScaleProfileParameters params = 
	{
	  min_width, max_slope, -1e10, lb_scale_1, lb_slope, nsamples, alpha
	};
      istream & in = std::cin;
      while ( ( ! in.eof() ) && in.good() )
	{
	  unsigned int x, y;
	  in >> x >> y;
	  if ( ( x < image.width() ) && ( y < image.height() ) )
	    displayScaleProfileWithGnuplot( tch, x, y, params );
	  else 
	    std::cerr << "Invalid coordinates: " << x << " " << y << std::endl;
	}
    }
  
  return error ? 1 : 0;
}
