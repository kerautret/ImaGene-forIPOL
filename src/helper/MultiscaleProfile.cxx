///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : MultiscaleProfile.cxx
//
// Creation : 2009/07/10
//
// Version : 2009/07/10
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
//	2009/07/10 : Mr ?Name? : ?What?
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //


///////////////////////////////////////////////////////////////////////////////
#include "ImaGene/mathutils/SimpleLinearRegression.h"
#include "ImaGene/mathutils/Statistics.h"
#include "ImaGene/helper/MultiscaleProfile.h"

// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/helper/MultiscaleProfile.ih"
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const MultiscaleProfile_RCS_ID = "@(#)class MultiscaleProfile definition.";



///////////////////////////////////////////////////////////////////////////////
// class MultiscaleProfile
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::MultiscaleProfile::~MultiscaleProfile()
{
}

/**
 * Default constructor.
 */
ImaGene::MultiscaleProfile::MultiscaleProfile()
{
}

/**
 * Computes the multiresolutions of the given Freeman chain [src]
 * from (1,1) up to (h,v)=(r,r) with shifts. Starts profile
 * analysis.
 *
 * @param src the source Freeman chain.
 * @param r the maximal resolution.
 */
void
ImaGene::MultiscaleProfile::init( const FreemanChain & src, uint r )
{
  MultiscaleFreemanChain::init( src, r );
  cerr << "+--- computing length statistics " << flush;
  uint src_size = src.chain.size();
  all_stats.resize( r );
  for( int k = 0; k < r; k++ ) {
    int res = k + 1;
    cerr << "." << res << flush;
    all_stats[ k ].scale = res;
    all_stats[ k ].stats = new Statistics( src_size, true );
    bool first = true;
    for(int x0 = 0; x0 < res; x0++ ) {
      for(int y0 = 0; y0 < res; y0++ ) {	  
	MultiscaleFreemanChain::SubsampledChainKey key( res, res, x0, y0 );
	const MultiscaleFreemanChain::SubsampledChain* ptrsub
	  = get( key );
	uint size = ptrsub->subc.chain.size();
	// Computes ms length statistics for one shift. 
	FreemanChain subc_copy;
	subc_copy.x0 = ptrsub->subc.x0;
	subc_copy.y0 = ptrsub->subc.y0;
	subc_copy.chain = ptrsub->subc.chain;
	Statistics* stats1 = getStatsMaximalSegments( subc_copy );
	// Relates these statistics to surfels on the original contour.
	for ( uint i = 0; i < src_size; ++i )
	  {
	    double mean = stats1->mean( ptrsub->c2subc[ i ] );
	    double max = stats1->max( ptrsub->c2subc[ i ] );
	    all_stats[ k ].stats->addValue( i, mean );
	    if ( first )
	      {
		all_stats[ k ].longest_ms.push_back
		  ( std::make_pair( key, max ) );
		all_stats[ k ].longest_mean.push_back
		  ( std::make_pair( key, mean ) );
	      }
	    else
	      {
		if ( all_stats[ k ].longest_ms[ i ].second < max ) 
		  all_stats[ k ].longest_ms[ i ] 
		    = std::make_pair( key, max );
		if ( all_stats[ k ].longest_mean[ i ].second < mean ) 
		  all_stats[ k ].longest_mean[ i ] 
		    = std::make_pair( key, mean );
	      }
	  }
	delete stats1;
	first = false;
      }
    }
    all_stats[ k ].stats->terminate();
  }
  cerr << " ended." << endl;
}

/**
 * @param x (returns) the x-value of the profile (log(scale+1)).
 * @param y (returns) the y-value of the profile (log(length of ms)).
 * @param idx the index of the surfel.
 */
void
ImaGene::MultiscaleProfile::profile( std::vector<double> & x, 
				     std::vector<double> & y,  
				     uint idx ) const
{
  for ( uint i = 0; i < all_stats.size(); ++i )
    {
      x.push_back( log( all_stats[ i ].scale ) );
      y.push_back( log( all_stats[ i ].stats->mean( idx ) ) );
    }
}

/**
 * @param x (returns) the x-value of the profile (log(scale+1)).
 * @param y (returns) the y-value of the profile (log(length of ms)).
 * @param idx the index of the surfel.
 */
void
ImaGene::MultiscaleProfile::profileFromMedian( std::vector<double> & x, 
					       std::vector<double> & y,  
					       uint idx ) const
{
  for ( uint i = 0; i < all_stats.size(); ++i )
    {
      x.push_back( log( all_stats[ i ].scale ) );
      y.push_back( log( all_stats[ i ].stats->median( idx ) ) );
    }
}






/**
 * Rajout BK 23/09/09
 *
 * Compute the profile from several linear regressions starting from
 * maximal scales with [n] samples with the parameter [alpha]. When
 * the first linear regression stop at scale K a second linear
 * regression is applyed from the scale [K,K+n] with the same parameter.
 *  
 * 
 * @param x (returns) the x-value of the profile (log(scale+1)).
 * @param y (returns) the y-value of the profile (log(length of ms)).
 * @param idx the index of the surfel.
 * @param n the minimum number of samples (min is 3).
 * @param alpha is the proportion of rejected linear model (the ones
 * with big variance).

 */
void
ImaGene::MultiscaleProfile::profileFromLinearReg( std::vector<double> & x, 
						  std::vector<double> & y,  
						  std::vector<uint> & indiceScales, 
						  uint idx,
						  uint n, 
						  double alpha) const
{
  vector<double> xAll;
  vector<double> yAll;
  vector<uint> nb;
  this->detailedProfile( xAll, yAll, nb, idx );
  vector<double> px;
  vector<double> py;
  this->profile( px, py, idx );
  x.clear();
  y.clear();
  Vector2D intervalBeg;
  Vector2D intervalEnd;
  
    
  int posStartRec = nb.size()-1 ;
  while(posStartRec>0){
    int dec= (int) (nb.size()-1)-posStartRec;
    indiceScales.push_back(posStartRec);
    posStartRec  = this->computeRegLinearPart( xAll, yAll, px, py, intervalBeg, intervalEnd, nb, n,  alpha, dec );  
    indiceScales.push_back(posStartRec);
    x.push_back(intervalEnd.x());
    y.push_back(intervalEnd.y());    
    x.push_back(intervalBeg.x());
    y.push_back(intervalBeg.y());    
    
  }
  
}








/**
 * Computes the detailed profile for the surfel of index [idx]. All
 * values (scale, length of ms) are inserted in [x] and [y],
 * contrarily to 'profile' which inserts only the mean values. The
 * vector [nb] contains the number of samples for each scale.
 * 
 * @param x (returns) the x-value of the profile (log(scale+1)).
 * @param y (returns) the y-value of the profile (log(length of ms)).
 * @param nb (returns) the number of samples for each scale.
 * @param idx the index of the surfel.
 */
void
ImaGene::MultiscaleProfile::detailedProfile
( std::vector<double> & x, 
  std::vector<double> & y,  
  std::vector<uint> & nb,  
  uint idx ) const
{
  for ( uint i = 0; i < all_stats.size(); ++i )
    {
      double xv = log( all_stats[ i ].scale );
      nb.push_back( all_stats[ i ].stats->samples( idx ) );
      for ( uint j = 0; j < all_stats[ i ].stats->samples( idx ); ++j )
	{
	  x.push_back( xv );
	  y.push_back( log( all_stats[ i ].stats->value( idx, j ) ) );
	}
    }
}


/**
 * A meaningful scale is an interval of scales of length no
 * smaller than [min_width] and in which the profile has slopes
 * below [max_slope]. This method returns the sequence of
 * meaningful scales for surfel [idx].
 *
 * @param intervals (returns) a list of meaningful scales.
 * @param idx the surfel of interest.
 * @param min_width the minimum length for the meaningful scales.
 * @param max_slope the maximum allowed slope for length evolution.
 */
void
ImaGene::MultiscaleProfile::meaningfulScales
( std::vector< std::pair< uint, uint > > & intervals,
  uint idx,
  uint min_width,
  double max_slope ) const
{
  vector<double> x;
  vector<double> y;
  this->profile( x, y, idx );
  uint l = 0;
  for ( uint k = 1; k < x.size(); ++k )
    {
      double slope = ( y[ k ] - y[ k - 1 ] ) / ( x[ k ] - x[ k - 1 ] );
      if ( ( slope > max_slope ) || ( (k+1) == x.size() ) )
	{
	  if ( ( k - 1 - l ) >= min_width )
	    intervals.push_back( std::make_pair( l+1, k ) );
	  l = k;
	}
    }
}

/**
 * The noise level is the first scale of the first meaningful
 * scale. A meaningful scale is an interval of scales of length no
 * smaller than [min_width] and in which the profile has slopes
 * below [max_slope]. 
 *
 * @param idx the surfel of interest.
 * @param min_width the minimum length for the meaningful scales.
 * @param max_slope the maximum allowed slope for length evolution.
 * @return the noise level or zero is none was found.
 * @see meaningfulScales
 */
uint
ImaGene::MultiscaleProfile::noiseLevel
( uint idx,
  uint min_width,
  double max_slope ) const
{
  std::vector< std::pair< uint, uint > > intervals;
  this->meaningfulScales( intervals, idx, min_width, max_slope );
  return ( intervals.size() != 0 )
    ? intervals[ 0 ].first
    : 0;
}



/**
 * The standard scale is the first scale starting from which the
 * profile is straight until the maximum scale. The straightness is
 * evaluated through a statistic test based on a simple linear
 * regression model. It requires two parameters: [n] is the minimum
 * number of samples to fit a linear model, 1-[alpha] is the
 * proportion of accepted linear model of the test (99%, alpha=0.01,
 * means that 99% of all linear model with a Gaussian noise are
 * accepted).
 *
 * @param idx the surfel of interest.
 * @param n the minimum number of samples (min is 3).
 * @param alpha is the proportion of rejected linear model (the ones
 * with big variance).
 *
 * @return the standard scale at point [idx].
 */
uint
ImaGene::MultiscaleProfile::standardScale
( uint idx,
  uint n,
  double alpha ) const
{
  vector<double> x;
  vector<double> y;
  this->profile( x, y, idx );
  SimpleLinearRegression SLR;
  vector<double>::const_reverse_iterator itx = x.rbegin(); 
  vector<double>::const_reverse_iterator itxe = x.rend(); 
  vector<double>::const_reverse_iterator ity = y.rbegin(); 
  SLR.addSamples( itx, itx + n, ity );
  SLR.computeRegression();
  itx += n;
  ity += n;
  uint l = x.size() - n + 1;
  for ( ; itx != itxe; ++itx, ++ity, --l )
    {
      pair<double,double> ic;
      ic = SLR.trustIntervalForY( *itx, alpha );
      if ( ( *ity < ic.first ) || ( *ity > ic.second ) )
	break;
      SLR.addSample( *itx, *ity );
      SLR.computeRegression();
    }
  cerr << "- i=" << idx << " nl=" << l 
       << " B=[" << SLR.intercept() << " " << SLR.slope()
       <<  "] sigma2=" << SLR.estimateVariance() << endl;
  
  return l;
}




/**
 * The standard scale is the first scale starting from which the
 * profile is straight until the maximum scale. The straightness is
 * evaluated through a statistic test based on a simple linear
 * regression model. It requires two parameters: [n] is the minimum
 * number of samples to fit a linear model, 1-[alpha] is the
 * proportion of accepted linear model of the test (99%, alpha=0.01,
 * means that 99% of all linear model with a Gaussian noise are
 * accepted).
 *
 * @param idx the surfel of interest.
 * @param n the minimum number of samples (min is 3).
 * @param alpha is the proportion of rejected linear model (the ones
 * with big variance).
 *
 * @return the standard scale at point [idx].
 */
uint
ImaGene::MultiscaleProfile::detailedStandardScale
( uint idx,
  uint n,
  double alpha ) const
{
  vector<double> x;
  vector<double> y;
  vector<uint> nb;
  this->detailedProfile( x, y, nb, idx );
  vector<double> px;
  vector<double> py;
  this->profile( px, py, idx );
  SimpleLinearRegression SLR;
  int i = nb.size() - 1;
  int j = x.size() - 1;
  while ( i >= ( nb.size() - 1 - n) )
    {
      for ( uint k = 0; k < nb[ i ]; ++k )
	{
	  SLR.addSample( x[ j ], y[ j ] );
	  --j;
	}
      --i;
    }
  if ( ! SLR.computeRegression() )
    cerr << "PROBLEM in SLR" << endl;
  // cerr << "  i=" << idx 
  //      << " B=[" << SLR.intercept() << " " << SLR.slope()
  //      <<  "] sigma2=" << SLR.estimateVariance() << endl;
  while ( i >= 0 )
    {
      pair<double,double> ic;
      ic = SLR.trustIntervalForY( px[ i ], alpha );
      // cerr << "  l=" << i << " nb[l]=" << nb[ i ] 
      // 	   << " IC_(" << ((1.0-alpha)*100.0) 
      // 	   << "%) (" << px[ i ] << ", " << py[ i ]
      // 	   << ") = " << ic.first << " - " << ic.second 
      // 	   << " : " << ( ic.second - ic.first ) << endl;
      if ( ( py[ i ] < ic.first ) || ( py[ i ] > ic.second ) )
	break;
      for ( uint k = 0; k < nb[ i ]; ++k )
	{
	  // cerr << " (" << x[ j ] << " " << y[ j ] << ")";
	  SLR.addSample( x[ j ], y[ j ] );
	  --j;
	}
      if ( ! SLR.computeRegression() )
	cerr << "PROBLEM in SLR" << endl;
      --i;
    }

  ++i;
  // cerr << "  i=" << idx << " nl=" << i 
  //      << " B=[" << SLR.intercept() << " " << SLR.slope()
  //      <<  "] sigma2=" << SLR.estimateVariance()
  //      << " sx=" << SLR.m_X.size() 
  //      << " sy=" << SLR.m_Y.size() 
  //      << endl;
  return i + 1;
}



/**
 * The standard detailed scale max is the last scale obtained after
 * several linear regression with consecutive negative slopes. If no
 * negative slopes are found the resulting scale is the maximal scale.
 *
 *The straightness is
 * evaluated through a statistic test based on a simple linear
 * regression model. It requires two parameters: [n] is the minimum
 * number of samples to fit a linear model, 1-[alpha] is the
 * proportion of accepted linear model of the test (99%, alpha=0.01,
 * means that 99% of all linear model with a Gaussian noise are
 * accepted).
 *
 * @param idx the surfel of interest.
 * @param n the minimum number of samples (min is 3).
 * @param alpha is the proportion of rejected linear model (the ones
 * with big variance).
 *
 * @return the standard scale at point [idx].
 */
uint
ImaGene::MultiscaleProfile::detailedStandardScaleMax
( uint idx,
  uint n,
  double alpha ) const
{

  vector<double> x;
  vector<double> y;
  vector<uint> indiceScales;
  this->profileFromLinearReg(x, y, indiceScales, idx, n, alpha);
  
  bool penteNeg = (y.at(0)-y.at(1))/(x.at(0)-x.at(1))<0;
  int i=0;
  bool firstCandidat = false;
  while((penteNeg||!firstCandidat) && i<x.size()){
    double pente = (y.at(i)-y.at(i+1))/(x.at(i)-x.at(i+1));
    penteNeg= (pente<0); 
    if(penteNeg){
      firstCandidat=true;
    }
    i+=2;
  }
  // cas où il n'existe que des pentes négatives:
  if(!penteNeg && !firstCandidat ){
    return x.size()-1;
  }
  
  if((penteNeg||!firstCandidat)){
    i=i-2;
  }else{
    i=i-4;
  }

  return indiceScales.at(i+1)+1;

}







/**
 * Test BK 21/09/09
 *
 * The maximal valid scale is the longest interval presenting
 * slopes less than maxSlope. 
 * 
 *  @return the first scale of the longest interval at point [idx].
 */

uint
ImaGene::MultiscaleProfile::maximalValidScale
( uint idx, double maxSlope ) const
{
  uint sizeMax = all_stats.size();  
  vector< pair<uint, uint> >  vectIntervals ;
  meaningfulScales(vectIntervals, idx,1, maxSlope ); 
  if( vectIntervals.size()==0){
    return sizeMax;
  }  
  uint maxIndice= 0;
  uint tailleIntervalMax = vectIntervals.at(0).second - vectIntervals.at(0).first;
  
  for(int i=0; i< vectIntervals.size(); i++){
    pair<uint, uint> interval = vectIntervals.at(i);
    uint taille = interval.second - interval.first;
    if (taille > tailleIntervalMax  ){
      tailleIntervalMax = taille;
      maxIndice = i;
    }
  }  
  
  return vectIntervals.at(maxIndice).first;
}


/**
 * Test BK 22/09/09
 *
 * The maxMeaningfulScale is a meaningfulscale with the maximal length of maximal segments.
 *
 *  @return the first scale of the interval with longest maximal segment at point [idx].
 */


uint
ImaGene::MultiscaleProfile::maxMeaningfulScale
( uint idx, uint minSize,  double maxSlope) const
{
  vector<double> x;
  vector<double> y;
  this->profile( x, y, idx );
  
  uint sizeMax = all_stats.size();  
  vector< pair<uint, uint> >  vectIntervals ;
  meaningfulScales(vectIntervals, idx, minSize, maxSlope ); 
  if( vectIntervals.size()==0){
    return sizeMax;
  }
  
  uint  indiceIntervalMax= 0;
  uint indicePointDebIntervalMax = vectIntervals.at(indiceIntervalMax).first;
  double tailleSegmentsMax = y.at(indicePointDebIntervalMax);
  
  for(int i=1; i< vectIntervals.size(); i++){
    pair<uint, uint> interval = vectIntervals.at(i);
    uint indiceDeb = interval.first;
    double tailleSegments =  y.at(indiceDeb);
    
    if (tailleSegments > tailleSegmentsMax  ){
      tailleSegmentsMax = tailleSegments;
      indiceIntervalMax = i;
    }
  }  
  

  return vectIntervals.at(indiceIntervalMax).first;
}






/**
 * Returns the key specifying the subsampled chain where this
 * surfel has the longest maximal segment in its surrounding.
 *
 * @param idx a valid surfel index.
 * @param scale the desired scale.
 * @return the specified key.
 */
ImaGene::MultiscaleFreemanChain::SubsampledChainKey
ImaGene::MultiscaleProfile::longestMaxAtScale( uint idx, uint scale ) const
{
  ASSERT_MultiscaleProfile( ( scale > 0 ) && ( scale < all_stats.size() ) );

  return all_stats[ scale - 1 ].longest_ms[ idx ].first;
}

/**
 * Returns the key specifying the subsampled chain where this
 * surfel has the greatest mean for the lengths of maximal
 * segments in its surrounding.
 *
 * @param idx a valid surfel index.
 * @param scale the desired scale.
 * @return the specified key.
 */
ImaGene::MultiscaleFreemanChain::SubsampledChainKey
ImaGene::MultiscaleProfile::longestMeanAtScale( uint idx, uint scale ) const
{
  ASSERT_MultiscaleProfile( ( scale > 0 ) && ( scale < all_stats.size() ) );

  return all_stats[ scale - 1 ].longest_mean[ idx ].first;
}



///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param that_stream the output stream where the object is written.
 */
void 
ImaGene::MultiscaleProfile::selfDisplay( ostream& that_stream ) const
{
  that_stream << "[MultiscaleProfile]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::MultiscaleProfile::OK() const
{
  return true;
}



///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////












/**
 * Compute a linear regression starting from the maximal scale with a
 * [shift] to the initial scale. It returns the indice of the last
 * scale (indexed from 0) for which the linear is computed. The linear
 * part is represented by two points [beginInterval] and [endInterval].
 * 
 * @param xp  the x-value of the profile (log(scale+1)).
 * @param yp  the y-value (mean or median value) of the profile (log(length of ms)).
 * @param x  the x-values of the profile (log(scale+1)).
 * @param y  the y-values of the profile (log(length of ms)).
 * @param nb nb[i] contains the number of samples associated to x[i]. 
 * 
 * @param beginInterval (return) the first point of the linear part. 
 * @param endInterval (return) the last point of the linear part. 
 *
 * @param idx the index of the surfel.
 * @param n the minimum number of samples (min is 3).
 * @param alpha is the proportion of rejected linear model (the ones
 * with big variance).

 **/


uint 
ImaGene::MultiscaleProfile::computeRegLinearPart (const vector<double> & x, const vector<double> &y, 
						  const vector<double> & px, const vector<double> &py,
						  Vector2D & beginInterval, Vector2D &endInterval, 
						  const vector<uint> &nb, uint n,  double alpha, uint shift) const { 

    
  SimpleLinearRegression SLR;

  int i = nb.size() - 1;
  int j = x.size() - 1;
  
  // ajustement des valeurs de i et j par rapport au décalage
  for(uint k=0; k<shift; k++){    
    for(uint l =0;l< nb[nb.size()-1-k]; l++){
      --j;
    }
    --i;
  }
  
  int minPos = nb.size()-1-shift-n;
  // initialisation en prenant en compte les [n] points de départ
  while ( i >= minPos && i>=0)
    {
      for ( uint k = 0; k < nb[ i ]; ++k )
	{
	  SLR.addSample( x[ j ], y[ j ] );
	  --j;
	}
      --i;
    }
  if ( ! SLR.computeRegression() ){
    cerr << "PROBLEM in SLR" << endl;
  }

  // Cas particulier s'il y a pas assez d'échantillons, on retourne simplement 
  // la regression linéaire sur les échantillons les points extrémités.  

  if(i<0 ){
    beginInterval.x() = px.at(0);
    beginInterval.y() = SLR.estimateY(px.at(0));
    endInterval.x() = px.at(nb.size()-1-shift);
    endInterval.y() = SLR.estimateY(px.at(nb.size()-1-shift));    
    return 0;
  }
  
  //Calcul de la regression
  while ( i >= 0 ){
    pair<double,double> ic;
      ic = SLR.trustIntervalForY( px[ i ], alpha );
      
      if ( ( py[ i ] < ic.first ) || ( py[ i ] > ic.second ) )
	break;
      for ( uint k = 0; k < nb[ i ]; ++k )
	{
	  SLR.addSample( x[ j ], y[ j ] );
	  --j;
	}
      if ( ! SLR.computeRegression() )
	cerr << "PROBLEM in SLR" << endl;
      --i;
  }
  ++i;    
  
  beginInterval.x() = px.at(i);
  beginInterval.y() = SLR.estimateY( px.at(i) );
  endInterval.x() = px.at(nb.size()-1-shift);
  endInterval.y() = SLR.estimateY(px.at(nb.size()-1-shift));
  
  return i;
}
