#ifndef _Z3S
#define _Z3S

#include "ImaGene/planeRecognition/Z2.h"
#include "ImaGene/planeRecognition/utils.h"

// three-dimensional point

class S
{
  public :
  I x;
  I y;
  I z;

  S();
  S(const S & t);
  S & operator = (const S & A);
  S(I _x,I _y,I _z);

  void reduc();  // tq pgcd(composantes) = 1
};

S operator + (const S & a, const S & b);
S operator - (const S & a, const S & b);
S operator ^ (const S & a, const S & b);
I operator * (const S & a, const S & b);
S operator * (const I & a, const S & b);
S operator * (const S & b, const I & a);
S operator / (const S & b, const I & a);
//bool operator == (const S & b, const I & a);

#endif
