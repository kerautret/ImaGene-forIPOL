///////////////////////////////////////////////////////////////////////////////
// Test various mathematics functions.
///////////////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <deque>
#include "ImaGene/base/Vector2i.h"
#include "ImaGene/base/Arguments.h"
#include "ImaGene/base/StandardArguments.h"
#include "ImaGene/mathutils/SimpleLinearRegression.h"


using namespace std;
using namespace ImaGene;


static Arguments args;

int
main( int argc, char** argv ) 
{
  // -------------------------------------------------------------------------
  // Prepare arguments.
  args.addBooleanOption( "-slr", "-slr: test SimpleLinearRegression." );
  args.addOption( "-analyze_log_model", "-analyze_log_model <filename> <alpha>: determines a power model (using logscale) given in <filename> as the sequence of line 'x y'. 1-<alpha> is the confidence of the test", "-", "0.05" );
  args.addOption( "-analyze_model", "-analyze_model <filename> <alpha>: determines a linear model given in <filename> as the sequence of line 'x y'. 1-<alpha> is the confidence of the test", "-", "0.05" );
  if ( ( argc <= 0 ) 
       || ! args.readArguments( argc, argv ) ) 
    {
      cerr << args.usage( "test_Math", 
			  "Test various mathematics functions.",
			  "" ) << endl;
      return 1;
    }
  if ( args.check( "-slr" ) )
    {
      SimpleLinearRegression::test();
    }
  if ( args.check( "-analyze_log_model" ) )
    {
      ifstream input( args.getOption( "-analyze_log_model" )->getValue( 0 ).c_str() );
      double alpha = args.getOption( "-analyze_log_model" )->getDoubleValue( 1 );
      double k;
      double l;
      deque<double> x;
      deque<double> y;
      input >> k >> l;
      while ( input.good() && ! input.eof() ) 
	{
	  cerr << log( k ) << " " << log( l ) << endl;
	  x.push_front( log( k ) ); 
	  y.push_front( log( l ) ); 
	  input >> k >> l;
	}
      input.close();
      SimpleLinearRegression SLR;
      SLR.addSamples( x.begin(), x.begin() + 4, y.begin() );
      SLR.computeRegression();
      int i = 4;
      for ( ; i < x.size(); ++i )
	{
	  pair<double,double> ic;
	  ic = SLR.trustIntervalForY( x[ i ], alpha );
	  cerr << "IC_(" << ((1.0-alpha)*100.0) 
	       << "%) (" << x[ i ] 
	       << ") = " << ic.first << " - " << ic.second 
	       << " : " << ( ic.second - ic.first ) << endl;
	  if ( ( y[ i ] < ic.first ) || ( y[ i ] > ic.second ) )
	    {
	      cerr << "- Linear model from 0 to " << ( i - 1 )
		   << "/" << x.size() << endl;
	      break;
	    }
	  SLR.addSample( x[ i ], y[ i ] );
	  SLR.computeRegression();
	}
      cerr << "- B=[" << SLR.intercept() << " " << SLR.slope()
	   <<  "] sigma2=" << SLR.estimateVariance() << endl;
      cout << "# Analyze multiscale profile. Resulting linear part "
	   << endl;
      int lin_idx = i;
      for ( i = 0; i < x.size(); ++i )
	{
	  cout << x[ i ] << " " << y[ i ];
	  if ( i < lin_idx )
	    cout << " " << x[ i ] << " " << SLR.estimateY( x[ i ] );
	  else
	    cout << " " << x[ i ] << " " << 0.0;
	  cout << endl;
	}
    }
  if ( args.check( "-analyze_model" ) )
    {
      ifstream input( args.getOption( "-analyze_model" )->getValue( 0 ).c_str() );
      double alpha = args.getOption( "-analyze_model" )->getDoubleValue( 1 );
      double k;
      double l;
      deque<double> x;
      deque<double> y;
      input >> k >> l;
      while ( input.good() && ! input.eof() ) 
	{
	  cerr << k << " " << l << endl;
	  x.push_front( k ); 
	  y.push_front( l ); 
	  input >> k >> l;
	}
      input.close();
      SimpleLinearRegression SLR;
      SLR.addSamples( x.begin(), x.begin() + 4, y.begin() );
      SLR.computeRegression();
      int i = 4;
      for ( ; i < x.size(); ++i )
	{
	  pair<double,double> ic;
	  ic = SLR.trustIntervalForY( x[ i ], alpha );
	  cerr << "IC_(" << ((1.0-alpha)*100.0) 
	       << "%) (" << x[ i ] 
	       << ") = " << ic.first << " - " << ic.second 
	       << " : " << ( ic.second - ic.first ) << endl;
	  if ( ( y[ i ] < ic.first ) || ( y[ i ] > ic.second ) )
	    {
	      cerr << "- Linear model from 0 to " << ( i - 1 )
		   << "/" << x.size() << endl;
	      break;
	    }
	  SLR.addSample( x[ i ], y[ i ] );
	  SLR.computeRegression();
	}
      cerr << "- B=[" << SLR.intercept() << " " << SLR.slope()
	   <<  "] sigma2=" << SLR.estimateVariance() << endl;
      cout << "# Analyze multiscale profile. Resulting linear part "
	   << endl;
      int lin_idx = i;
      for ( i = 0; i < x.size(); ++i )
	{
	  cout << x[ i ] << " " << y[ i ];
	  if ( i < lin_idx )
	    cout << " " << x[ i ] << " " << SLR.estimateY( x[ i ] );
	  else
	    cout << " " << x[ i ] << " " << 0.0;
	  cout << endl;
	}
    }
}
