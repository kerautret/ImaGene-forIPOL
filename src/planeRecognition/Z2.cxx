// #include "stdafx.h"
#include "ImaGene/planeRecognition/Z2.h"
#include <stdlib.h>

P::P(const P &A)               { x = A.x;  y = A.y; }
P& P::operator = (const P &A)  { x = A.x;  y = A.y;  return *this; }
P::P(const I & a,const  I & b) { x = a; y = b;}

void P::rnd(int lg_max)
{
  do
  {
    x = ( rand() % (lg_max*2+1) ) - lg_max;
    y = ( rand() % (lg_max*2+1) ) - lg_max;
  } while ( (x==I(0)) && (y==I(0)) );
}

///////////////               operators       ///////////////////////////

P operator + (const P & A ,const P & B) { P tmp;  tmp.y = A.y + B.y;    tmp.x = A.x + B.x;  return tmp; }
P operator - (const P & A ,const P & B) { P tmp;  tmp.y = A.y - B.y;    tmp.x = A.x - B.x;  return tmp; }
P operator * (const P & A ,const I & B) { P tmp;  tmp.x = A.x * B;      tmp.y = A.y * B;    return tmp; }
P operator * (const I & B ,const P & A) { P tmp;  tmp.x = A.x * B;      tmp.y = A.y * B;    return tmp; }
I operator ^ (const P & A ,const P & B) { return A.x * B.y - B.x * A.y; }
I operator * (const P & A ,const P & B) { return A.x * B.x + B.y * A.y; }

bool operator == (const P & A, const P & B)
{
  if ( ( A.x == B.x ) && ( A.y == B.y ) )   return true;
  return false;
}

bool operator != (const P & A,const  P & B) { return !(A==B); }

