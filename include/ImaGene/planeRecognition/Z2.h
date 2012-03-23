//---------------------------------------------------------------------------
#ifndef V2_H
#define V2_H

#include <gmp.h>
#include <gmpxx.h>

#define Fori(a) for (int i = 0 ; i < (a) ; i++)
#define Forj(a) for (int j = 0 ; j < (a) ; j++)
#define Fork(a) for (int k = 0 ; k < (a) ; k++)
#define Forx(a) for (int x = 0 ; x < (a) ; x++)
#define Fory(a) for (int y = 0 ; y < (a) ; y++)

#include "ImaGene/arithmetic/IntegerComputer.h"

//the type I denotes a "big integer"
//library gmp
typedef mpz_class I;
// inline int get_si( const I & i )
// {
//   return i.get_si();
// }
// inline double get_d( const I & i )
// {
//   return i.get_d();
// }

// typedef long long I;
// inline int get_si( const I & i )
// {
//   return (int) i;
// }

// inline double get_d( const I & i )
// {
//   return (double) i;
// }

// represents a point in the plane
class P
{
  public :

  I x,y;

  P() {}
  P& operator = (const P & A);
  P(const P & A);
  P(const I & a, const I & b);

  void rnd(int lg_max);    // random function
};

P operator + (const P & A , const P & B);
P operator - (const P & A , const P & B);
P operator * (const P & A , const I & B);
P operator * (const I & B , const P & A);
I operator ^ (const P & A , const P & B);
I operator * (const P & A , const P & B);

bool operator == (const P & A, const P & B);
bool operator != (const P & A,const  P & B);


#endif
