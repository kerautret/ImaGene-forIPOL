// #include "stdafx.h"
#include "ImaGene/planeRecognition/utils.h"
#include<iostream>

using namespace std;



namespace ImaGene {




//return a solution of the Diophantine equation: na x + nb y = nc
P ExtendedEuclid( const I & na, const I & nb, const I & nc)
{
  int i, k, size;
  static I a,b,c;
  static I TabBezout[4][1000];
  static P v;

  size = 2;

  if( na == I_ZERO )
    return P( I_ZERO, nb * nc );
  if( nb == I_ZERO )
    {
      return P( na * nc, I_ZERO );
    }

  a=_abs(na);
  b=_abs(nb);
  c=nc;
  
  TabBezout[0][0] = a;
  TabBezout[0][1] = b;
  TabBezout[1][0] = I_ZERO;
  TabBezout[2][0] = I_ONE;
  TabBezout[3][0] = I_ZERO;
  TabBezout[2][1] = I_ZERO;
  TabBezout[3][1] = I_ONE;
  
  k=0;
  while( TabBezout[0][k+1] != I_ZERO )
    {
      size++;
      TabBezout[1][k+1] = TabBezout[0][k] / TabBezout[0][k+1];
      TabBezout[0][k+2] = TabBezout[0][k] % TabBezout[0][k+1];
      TabBezout[2][k+2] = TabBezout[2][k] - TabBezout[1][k+1]*TabBezout[2][k+1];
      TabBezout[3][k+2] = TabBezout[3][k] - TabBezout[1][k+1]*TabBezout[3][k+1];
      k++;
    }

  v.x = TabBezout[2][size-2];
  v.y = TabBezout[3][size-2];

  if( na < I_ZERO )
    v.x = -v.x;
  if( nb < I_ZERO )
    v.y = -v.y;
  v = v*c;
  return v;
}



//update the index of the last point of the set of points to consider
//we consider we consider new points for the "iteration"-th time
//D: diameter of the set (supposed odd)
bool incrementIndex(int  iteration, int & nbAddedPoints, int & indexLastPoint, int largeurPlan)
{
  nbAddedPoints = iteration*8;
  indexLastPoint += nbAddedPoints;
  if(iteration==largeurPlan/2) return true;
  else return false;
}

}
