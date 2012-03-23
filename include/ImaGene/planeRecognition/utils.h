//---------------------------------------------------------------------------
#ifndef math_Z2H
#define math_Z2H

#include <iostream>
#include "ImaGene/planeRecognition/Z2.h"

// static const I I_ZERO = I( 0 );
// static const I I_ONE = I( 1 );
// static const I I_MINUS_ONE = I( -1 );

// inline
// I _abs(const I & a )            
// { 
//   if ( a >= I_ZERO ) return a; 
//   return -a;
// }

// inline
// const I & _max(const I & a,const I & b) 
// { if ( b > a  )    return b; return a; }

// inline
// const I & _min(const I & a,const I & b) 
// { if ( a < b  )    return a; return b; }

namespace ImaGene {

//compute the floor value of na/nb
inline
I floor(const I &na, const I &nb)
{
  static I a;
  a = na;
  static I b;
  b = nb;
  if( ( a < I_ZERO ) && ( b < I_ZERO ) )
    {
      a=-a;
      b=-b;
    }
  else if( b < I_ZERO )
    {
      a=-a;
      b=-b;
    }
  if(( a > I_ZERO ) || (a % b == I_ZERO ))
    return a/b;
  else
    return a/b - I_ONE;
}

//compute the ceiling value of na/nb
inline
I ceil (const I &na, const I &nb)
{
  static I a;
  a = na;
  static I b;
  b = nb;
  if( ( a < I_ZERO ) && ( b < I_ZERO ) )
    {
      a=-a;
      b=-b;
    }
  else if( b < I_ZERO )
    {
      a=-a;
      b=-b;
    }
  
  if( ( a < I_ZERO ) || ( a % b == I_ZERO) )
    return a/b;
  else
    return ( a/b + I_ONE );
}

//compute gcd of a and b
inline
I _pgcd(const I & a,const I & b)
{
  //  std::cerr << "gcd(" << a << ", " << b << ")=";
  static I aa;
  aa = _abs(a);
  static I bb;
  bb = _abs(b);
  static I a0;
  a0 = _max(aa,bb);
  static I a1;
  a1 = _min(aa,bb);
  static I r;
  /* std::cerr << "(aa=" << aa << ")"; */
  /* std::cerr << "(bb=" << aa << ")"; */
  /* std::cerr << "(aa=" << aa << ")"; */
  /* std::cerr << "(aa=" << aa << ")"; */
  while ( a1 != I_ZERO )
  {
    r = a0 % a1;
    a0 = a1;
    a1 = r;
  }
  //  std::cerr << a0 << std::endl;
  return a0;
}

//make p irreducible
inline
void reduc(P & p)
{
  static I t;
  t = _pgcd(_abs(p.x),_abs(p.y));
  if ( ( t != I_ONE ) && ( t != I_ZERO ) )
    {
      p.x = p.x / t;
      p.y = p.y / t;
    }
}


//compute the cross product of u and v
inline
I cross_product( const P & u, const P & v)
{
  return ( u.x * v.y - u.y * v.x );
}

//compute the dot product of u and v
inline
I dot_product( const P & u, const P & v )
{
  return ( u.x * v.x + u.y * v.y );
}

//returns a solution of the Diophantine equation: a x + b y = c
P ExtendedEuclid( const I & a, const I & b, const I & c);

//computes the floor (fl) and the ceiling (ce) value of the real number k such that P+ku lies on the supporting line
//of the constraint N.p<=c
inline
void coefficientIntersection( const P & p, const P & u, const P & N, 
			      const I & c,
			      I & fl, I & ce)
{
//cout<<"begin coeff intersect"<<endl;
//cout<<"c="<<c<<"u="<<u.x<<" "<<u.y<<"p="<<p.x<<" "<<p.y<<"N="<<N.x<<" "<<N.y<<endl;
  fl = floor(c-dot_product(p,N),dot_product(u,N));
 // cout<<"inter"<<endl;
  ce = ceil(c-dot_product(p,N),dot_product(u,N));
 // cout<<"end coeff intersect"<<endl;
}


//compute the valid bezout vector v of u such that B+v satifies the
//constraints C and C2 and such that B+v+u doesn't satify the
//constraint C2 if redef==0, v is already initialized and the
//constraint N.p<=c is not used in the computation
inline
void validBezout ( const P & A, const P & u,
		   P & v, 
		   const P & N, const I & c, const P & N2, 
		   const I & c2, int redef)
{
  static I fl, ce;
  
  if(redef==1)
  {
    v = ExtendedEuclid(-u.y,u.x,1);
    if(dot_product((A+v),N)>c)
    {
      v.x = -v.x;
      v.y = -v.y;
    }
  }
  coefficientIntersection(A+v, u, N2, c2,fl,ce);
  v = v+u*(fl);
}

//updates the index of the last point of the set of points to consider
//incrementes the number of the iteration
//compute the number of added points
//largeurPlan: diameter of the set
bool incrementIndex(int iteration, int & nbAddedPoints, int & indexLastPoint, int largeurPlan);
}

#endif
