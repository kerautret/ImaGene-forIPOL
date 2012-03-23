// #include "stdafx.h"
#include "ImaGene/planeRecognition/Z3.h"

using namespace ImaGene;

S::S()                            {}
S::S(const S & t)                 { x = t.x; y = t.y; z = t.z; }
S & S::operator = (const S & A)   { x = A.x; y = A.y; z = A.z; return * this; }
S::S(I _x,I _y,I _z)              { x = _x;  y = _y;  z = _z; }

void S::reduc()
{
  I p = _pgcd(x,y);
  p = _pgcd(p,z);
  x = x / p;
  y = y / p;
  z = z / p;

}

S operator + (const S & a, const S & b)
{
  S t;
  t.x = a.x + b.x;
  t.y = a.y + b.y;
  t.z = a.z + b.z;
  return t;
}

S operator - (const S & a, const S & b)
{
  S t;
  t.x = a.x - b.x;
  t.y = a.y - b.y;
  t.z = a.z - b.z;
  return t;
}

S operator ^ (const S & a, const S & b)
{
  S t;
  t.x = a.y * b.z - a.z * b.y;
  t.y = a.z * b.x - a.x * b.z;
  t.z = a.x * b.y - a.y * b.x;
  return t;
}

I operator * (const S & a, const S & b)
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

S operator * (const I & a, const S & b)
{
  S t;
  t.x = b.x * a;
  t.y = b.y * a;
  t.z = b.z * a;
  return t;
}

S operator * (const S & b, const I & a)
{
  S t;
  t.x = b.x * a;
  t.y = b.y * a;
  t.z=  b.z * a;
  return t;
}

S operator / (const S & b, const I & a)
{
  S t;
  t.x = b.x / a;
  t.y = b.y / a;
  t.z=  b.z / a;
  return t;
}

/*bool operator == (const S & a, const S & b)
{
  if(a.x==b.x && a.y==b.y && a.z==b.z)
    return true;
  else
    return false;
}*/

