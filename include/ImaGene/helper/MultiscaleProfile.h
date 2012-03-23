/** @file MultiscaleProfile.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : MultiscaleProfile.h
//
// Creation : 2009/07/10
//
// Version : 2009/07/10
//
// Author : JOL
//
// Summary : Header file for files MultiscaleProfile.ih and MultiscaleProfile.cxx
//
// History :
//	2009/07/10 : ?Name? : ?What?
//
// Rcs Id : "@(#)class MultiscaleProfile declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(MultiscaleProfile_RECURSES)
#error Recursive header files inclusion detected in MultiscaleProfile.h
#else // defined(MultiscaleProfile_RECURSES)
#define MultiscaleProfile_RECURSES

#if !defined MultiscaleProfile_h
#define MultiscaleProfile_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <utility>
#include "ImaGene/dgeometry2d/MultiscaleFreemanChain.h"
#include "ImaGene/base/Vector.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  class Statistics;

  /////////////////////////////////////////////////////////////////////////////
  // class MultiscaleProfile
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'MultiscaleProfile' <p> Aim: Useful to
   * detect noise on Freeman chains based on multiscale decomposition.
   */
  class MultiscaleProfile : public MultiscaleFreemanChain
  {

    // ----------------------- Standard services ------------------------------
  public:
    struct LengthStatsAtScale
    {
      uint scale;
      Statistics* stats;
      std::vector< std::pair<SubsampledChainKey,double> > longest_ms; 
      std::vector< std::pair<SubsampledChainKey,double> > longest_mean; 
      LengthStatsAtScale()
	: stats( 0 )
      {}
      ~LengthStatsAtScale()
      {
	if ( stats != 0 ) delete stats;
      }
    };

    std::vector<LengthStatsAtScale> all_stats;

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~MultiscaleProfile();

    /**
     * Constructor.
     */
    MultiscaleProfile();
    
    /**
     * Computes the multiresolutions of the given Freeman chain [src]
     * from (1,1) up to (h,v)=(r,r) with shifts. Starts profile
     * analysis.
     *
     * @param src the source Freeman chain.
     * @param r the maximal resolution.
     */
    void init( const FreemanChain & src, uint r );

    /**
     * @param x (returns) the x-value of the profile (log(scale+1)).
     * @param y (returns) the y-value of the profile (log(length of ms)).
     * @param idx the index of the surfel.
     */
    void profile( std::vector<double> & x, 
		  std::vector<double> & y,  
		  uint idx ) const;


    /**
     * @param x (returns) the x-value of the profile (log(scale+1)).
     * @param y (returns) the y-value of the profile (log(length of ms)).
     * @param idx the index of the surfel.
     */
    void profileFromMedian( std::vector<double> & x, 
			    std::vector<double> & y,  
			    uint idx ) const;
    
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
    void detailedProfile( std::vector<double> & x, 
			  std::vector<double> & y,  
			  std::vector<uint> & nb,  
			  uint idx ) const;



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
    profileFromLinearReg( std::vector<double> & x, 
			  std::vector<double> & y,  
			  std::vector<uint> & scales,
			  uint idx,
			  uint n, 
			  double alpha) const;
      


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
    void meaningfulScales( std::vector< std::pair< uint, uint > > & intervals,
			   uint idx,
			   uint min_width = 1,
			   double max_slope = -0.2 ) const;

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
    uint noiseLevel( uint idx,
		     uint min_width = 1,
		     double max_slope = -0.2 ) const;

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
    uint standardScale( uint idx,
			uint n = 4,
			double alpha = 0.01 ) const;

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
    uint detailedStandardScale( uint idx,
				uint n = 4,
				double alpha = 0.01 ) const;




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
    uint detailedStandardScaleMax
    ( uint idx,
      uint n,
      double alpha ) const;
    
    
    /**
     * Test BK 21/09/09
     *
     * The maximal valid scale is the longest interval presenting
     * slopes less than maxSlope. 
     * 
     *  @return the first scale of the longest interval at point [idx].
     */
    uint maximalValidScale ( uint idx, double maxSlope ) const;




    /**
     * Test BK 22/09/09
     *
     * The maxMeaningfulScale is a meaningfulscale with the maximal length of maximal segments.
     *
     *  @return the first scale of the interval with longest maximal segment at point [idx].
     */
    
    uint maxMeaningfulScale ( uint idx, uint minSize,  double maxSlope ) const;
    
    


    
    /**
     * Returns the key specifying the subsampled chain where this
     * surfel has the longest maximal segment in its surrounding.
     *
     * @param idx a valid surfel index.
     * @param scale the desired scale.
     * @return the specified key.
     */
    SubsampledChainKey longestMaxAtScale( uint idx, uint scale ) const;

    /**
     * Returns the key specifying the subsampled chain where this
     * surfel has the greatest mean for the lengths of maximal
     * segments in its surrounding.
     *
     * @param idx a valid surfel index.
     * @param scale the desired scale.
     * @return the specified key.
     */
    SubsampledChainKey longestMeanAtScale( uint idx, uint scale ) const;

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param that_stream the output stream where the object is written.
     */
    void selfDisplay( std::ostream & that_stream ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool OK() const;
  

    // ------------------------- Datas ----------------------------------------
  private:


    // ------------------------- Hidden services ------------------------------
    //protected:
  public:

    
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
    computeRegLinearPart ( const std::vector<double> & x, 
			   const std::vector<double> & y, 
			   const std::vector<double> & px, 
			   const std::vector<double> & py,
			   Vector2D & beginInterval, 
			   Vector2D & endInterval, 
			   const std::vector<uint> & nb, 
			   uint n,  
			   double alpha, 
			   uint shift ) const;

    
  private:
    
    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE MultiscaleProfile( const MultiscaleProfile & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE MultiscaleProfile & operator=( const MultiscaleProfile & other );
  
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'MultiscaleProfile'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'MultiscaleProfile' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const MultiscaleProfile & that_object_to_display );


  



  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/helper/MultiscaleProfile.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined MultiscaleProfile_h

#undef MultiscaleProfile_RECURSES
#endif // else defined(MultiscaleProfile_RECURSES)
