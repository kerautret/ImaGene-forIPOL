// #include "stdafx.h"
#include "ImaGene/planeRecognition/plane_generator.h"
#include <stdlib.h>

PlaneGenerator::PlaneGenerator(int DD)   
{
  D = DD;
  L  = new S[D*D];

  int R = D/2;
  L[0] = S(R,R,0);
  int pos = 1;
  int X, Y;

  for(int limite=1;limite<=R;limite++)
  {
    X=R-limite;
    for(Y=R-limite;Y<=R+limite;Y++)     L[pos++] = S(X,Y,0);
    Y=R+limite;
    for(X=R-limite+1;X<=R+limite;X++)   L[pos++] = S(X,Y,0);
    X=R+limite;
    for(Y=R+limite-1;Y>=R-limite;Y--)   L[pos++] = S(X,Y,0);
    Y=R-limite;
    for(X=R+limite-1;X>=R-limite+1;X--) L[pos++] = S(X,Y,0);
  }
}

PlaneGenerator::~PlaneGenerator()  {   delete [] L; }


// for P(x,y,?) and a normal vector N, finds z such that 0<= N.P < |N.z|
I GiveZ(const S & N, const I & x, const I & y)
{
  I k = N.x * x + N.y * y;
  I z = -k / N.z;
  I p = k + N.z * z;
  I cc = N.z;
  if ( cc < I(0) ) cc = - cc;

  if ( (p          >= I(0)) && (p       < cc) ) return z;
  if ( (p+N.z      >= I(0)) && (p+N.z   < cc) ) return z+I(1);
  if ( (p-N.z      >= I(0)) && (p-N.z   < cc) ) return z-I(1);
  if ( (p+I(2)*N.z >= I(0)) && (p+I(2)*N.z < cc) ) return z+I(2);
  if ( (p-I(2)*N.z >= I(0)) && (p-I(2)*N.z < cc) ) return z-I(2);
  throw;
}



void PlaneGenerator::AnotherNormal(S & N)
{
  N.z = -(rand()%(D*D*2)+1);
  do
  {
     N.x = I(rand())%(-N.z);
     N.y = I(rand())%(-N.z);
  } while ( _pgcd(_pgcd(N.x,N.y),N.z) != I(1) );
}


S PlaneGenerator::Random()
{
  S N;
  AnotherNormal(N);
    
  //initialize left, right, down, up
  //to recall of 4 characteristic points of the digital plane
  minX = S(0,0,GiveZ(N,0,0));
  maxX = S(D-1,0,GiveZ(N,D-1,0));
  minY = minX;
  maxY = S(0,D-1,GiveZ(N,0,D-1));

  int n=0;
  Forx(D)
    Fory(D)
      L[n++].z  = GiveZ(N,L[n].x,L[n].y);
	
  return N;
}


int PlaneGenerator::size()
{
  return D*D;
}

