///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// File name : SimpleLinearRegression.cxx
//
// Creation : 2009/09/01
//
// Version : 2009/09/01
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
//	2009/09/01 : Mr ?Name? : ?What?
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//                                                                           //


///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <boost/math/distributions/students_t.hpp>
#include "ImaGene/mathutils/SimpleLinearRegression.h"
// Includes inline functions/methods if necessary.
#if !defined(INLINE)
#include "ImaGene/mathutils/SimpleLinearRegression.ih"
#endif
///////////////////////////////////////////////////////////////////////////////

using namespace std;

const char* const SimpleLinearRegression_RCS_ID = "@(#)class SimpleLinearRegression definition.";



///////////////////////////////////////////////////////////////////////////////
// class SimpleLinearRegression
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

/**
 * Destructor. 
 */
ImaGene::SimpleLinearRegression::~SimpleLinearRegression()
{
}

/**
 * Constructor.
 * The object is invalid.
 *
 * @param eps_zero the value below which the absolute value of the
 * determinant is considered null.
 */
ImaGene::SimpleLinearRegression::SimpleLinearRegression( double eps_zero )
  : epsilon_zero( eps_zero ), m_n( 0 )
{
  clear();
}

/**
 * Clear all datas.
 */
void
ImaGene::SimpleLinearRegression::clear()
{
  m_n = 0;
  m_sum_x = 0.0;
  m_sum_x2 = 0.0;
  m_sum_y = 0.0;
  m_sum_xy = 0.0;
  m_Y.clear();
  m_X.clear();
  m_B[ 0 ] = 0.0;
  m_B[ 1 ] = 0.0;
  m_d = 0.0;
}

/**
 * Computes the regression of the current parameters.
 *
 * @return 'true' if the regression was valid (non null number of
 * samples, rank of X is 2), 'false' otherwise.
 */
bool
ImaGene::SimpleLinearRegression::computeRegression()
{
  if ( m_n <= 2 ) return false;
  m_d = m_n * m_sum_x2 - ( m_sum_x * m_sum_x );
  if ( ( m_d > -epsilon_zero ) && ( m_d < epsilon_zero ) )
    {
      m_d = 0.0;
      return false;
    }
  m_B[ 0 ] = ( m_sum_x2 * m_sum_y - m_sum_x * m_sum_xy ) / m_d;
  m_B[ 1 ] = ( -m_sum_x * m_sum_y + m_n * m_sum_xy ) / m_d;

  m_U.clear();
  m_norm_U2 = 0.0;
  for ( unsigned int i = 0; i < m_n; ++i )
    {
      double u = m_Y[ i ] - m_B[ 0 ] - m_B[ 1 ] * m_X[ i ];
      m_U.push_back( u );
      m_norm_U2 += u * u;
    }
  return true;
}



/**
 * Given a test confidence value (1-[a]), return the expected interval
 * of value for Y, given a new [x], so that the model is still
 * linear. One may thus check if a new pair (y,x) is still in the
 * current linear model or not.
 *
 * @param x any x value.
 *
 * @param a the expected confidence value for the test (a=0.05
 * means 95% of confidence).
 *
 * @return the expected interval [min_y, max_y] such that any
 * value y within confirms the current linear model.
 */
std::pair<double,double>
ImaGene::SimpleLinearRegression::trustIntervalForY( double x, double a ) const
{
  double t = ( m_sum_x2 - 2.0 * x * m_sum_x + m_n * x * x ) / m_d;
  double c = sqrt( estimateVariance() * ( 1 + t ) );
  boost::math::students_t_distribution<double> T( m_n - 2 );
  double q = boost::math::quantile( T, 1.0 - a/2.0 );
  return make_pair( estimateY( x ) - q*c, estimateY( x ) + q*c );
}



///////////////////////////////////////////////////////////////////////////////
// Interface - public :

/**
 * Writes/Displays the object on an output stream.
 * @param that_stream the output stream where the object is written.
 */
void 
ImaGene::SimpleLinearRegression::selfDisplay( ostream& that_stream ) const
{
  that_stream << "[SimpleLinearRegression]";
}

/**
 * Checks the validity/consistency of the object.
 * @return 'true' if the object is valid, 'false' otherwise.
 */
bool 
ImaGene::SimpleLinearRegression::OK() const
{
  return true;
}


///////////////////////////////////////////////////////////////////////////////

/**
 * Several tests for the class.
 */
bool
ImaGene::SimpleLinearRegression::test()
{
  SimpleLinearRegression SLR;
  SLR.addSample( 1.0, 2.0 );
  SLR.addSample( 2.0, 3.0 );
  SLR.addSample( 4.0, 4.0 );
  SLR.addSample( 5.0, 5.5 );
  SLR.computeRegression();
  cout << "B=[" << SLR.intercept() << " " << SLR.slope()
       <<  "] sigma2=" << SLR.estimateVariance() << endl;
  pair<double,double> ic;
  ic = SLR.trustIntervalForY( 6.0, 0.2 );
  cout << "IC_80%(6.0)=" << ic.first << " - " << ic.second << endl;
  ic = SLR.trustIntervalForY( 6.0, 0.1 );
  cout << "IC_90%(6.0)=" << ic.first << " - " << ic.second << endl;
  ic = SLR.trustIntervalForY( 6.0, 0.05 );
  cout << "IC_95%(6.0)=" << ic.first << " - " << ic.second << endl;
  return true;
}


///////////////////////////////////////////////////////////////////////////////
// Internals - private :

//                                                                           //
///////////////////////////////////////////////////////////////////////////////
