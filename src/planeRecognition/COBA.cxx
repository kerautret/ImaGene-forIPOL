#pragma hdrstop

#include <cstdlib>
#include<iostream>
#include<cmath>


#include "ImaGene/planeRecognition/COBA.h"


using namespace std;
using namespace ImaGene;

bool COBA::oracleFloatWidth(int axe, float width)
{
  Centroid.reduc();  //the vector becomes irreducible
  switch(axe){
  case 0 : N = S(Centroid.z*g, Centroid.x, Centroid.y); break; //3D normal vector
  case 1 : N = S(Centroid.x, Centroid.z*g, Centroid.y); break; //3D normal vector
  case 2 : N = S(Centroid.x, Centroid.y, Centroid.z*g); break; //3D normal vector
  }
  min = N*set[0]; // dot product
  max = min;
  indmax = indmin = 0;

  //look for the points defining the min dot product and the max dot product
  static I v;
  Fori(cardinality)
  {
    v = N*set[i];
    if ( v > max ) { max = v;  indmax = i; }
    if ( v < min ) { min = v;  indmin = i; }
  }
  
  static I maxValue;
  switch(axe){
  case 0 : maxValue=_abs(N.x); break; //3D normal vector
  case 1 : maxValue=_abs(N.y); break; //3D normal vector
  case 2 : maxValue=_abs(N.z); break; //3D normal vector
  }
 // cout<<"maxvalue "<<maxValue<<" width "<<width<<" product "<<maxValue.get_d()*width<<endl;
  
  if ( max - min < get_d( maxValue )*width )
  {
    return true;      // is DP
  }
  else   //is not DP
  {
    // computation of the gradient
    switch(axe){
     case 0 : Gradient = P(set[indmin].y-set[indmax].y,set[indmin].z-set[indmax].z); break; //3D normal vector
     case 1 : Gradient = P(set[indmin].x-set[indmax].x,set[indmin].z-set[indmax].z); break; //3D normal vector
     case 2 : Gradient = P(set[indmin].x-set[indmax].x,set[indmin].y-set[indmax].y); break; //3D normal vector
    }
    //cout<<"oracle gradient : ("<<Gradient.x<<","<<Gradient.y<<")"<<endl;
   /* if(Gradient.x==0 && Gradient.y==0)
    {
      cout<<"oracle : Gradient nul!!!!!!!!!!!"<<endl;
      exit(-1);
    }*/
    return false;
  }
}

bool COBA::oracleLightFloatWidth(int axe, float width)
{
  //look for the points defining the min dot product and the max dot product
  static I v;
  v = N*set[cardinality-1];
  if ( v > max ) { max = v;  indmax = cardinality-1; }
  if ( v < min ) { min = v;  indmin = cardinality-1; }
 
  static I maxValue;
  switch(axe){
  case 0 : maxValue=_abs(N.x); break; //3D normal vector
  case 1 : maxValue=_abs(N.y); break; //3D normal vector
  case 2 : maxValue=_abs(N.z); break; //3D normal vector
  }

  if ( max - min < get_d( maxValue )*width ) return true;      // is DP
  else   //is not DP
  {
    // computation of the gradient
    switch(axe){
     case 0 : Gradient = P(set[indmin].y-set[indmax].y,set[indmin].z-set[indmax].z); break; //3D normal vector
     case 1 : Gradient = P(set[indmin].x-set[indmax].x,set[indmin].z-set[indmax].z); break; //3D normal vector
     case 2 : Gradient = P(set[indmin].x-set[indmax].x,set[indmin].y-set[indmax].y); break; //3D normal vector
    }

    return false;
  }
}


void COBA::standardInit()
{
  //initialize the search space as a square
  CC.Add(P(-g,-g));
  CC.Add(P(g,-g));
  CC.Add(P(g,g));
  CC.Add(P(-g,g));
}

COBA::COBA(S* vSet, int D, int axe)
{
  set = vSet;
  cardinality = 1;

  indmax = 0;
  indmin = 0;

  g = I(2)*D*D*D; //digitalization step

  standardInit(); //initialization of the search space
  CC.centroid(Centroid);  //Centroid.z > 0
  Centroid.reduc();  //the vector becomes irreducible
  
  CCsave.free();
  CCsave = CC;
  Centroidsave = Centroid;
  
  switch(axe){
     case 0 : N = S(Centroid.z*g, Centroid.x, Centroid.y); break; //3D normal vector
     case 1 : N = S(Centroid.x, Centroid.z*g, Centroid.y); break; //3D normal vector
     case 2 : N = S(Centroid.x, Centroid.y, Centroid.z*g); break; //3D normal vector
    }
  
  min = N*set[0]; // dot product
  max = min;
  
  Nsave = N;
  maxsave = max;
  minsave = min;
  indmaxsave = 0;
  indminsave = 0;
}


void COBA::doubleCutFloatWidth(int axe, float width)
{
  static I constante, constantebis;

  //2 cuts on the search space:
  //Gradient.p <= constante
  //-Gradient.p <= constantebis
 
  int tmp1 = (int) floor(get_si( g )*width);
  int tmp2 = (int) ceil(get_si( g )*width);

  switch(axe){
     case 0 : constante = tmp1-g*(set[indmin].x-set[indmax].x)-I(1);
              constantebis = tmp2+g*(set[indmin].x-set[indmax].x)-I(1);
              break;
     case 1 : constante =tmp1-g*(set[indmin].y-set[indmax].y)-I(1);
              constantebis = tmp2+g*(set[indmin].y-set[indmax].y)-I(1);
              break;
     case 2 : constante = tmp1-g*(set[indmin].z-set[indmax].z)-I(1);
              constantebis = tmp2+g*(set[indmin].z-set[indmax].z)-I(1);
              break;
    }
  
  CC.cutOpti(Gradient,constante);

  CC.cutOpti(P(0,0)-Gradient,constantebis);
}

bool COBA::runFloatWidth(int axe, float width, bool failPrevious, bool allowNewNormal)
{
  bool result;
  
  if(failPrevious)
  {
    if(allowNewNormal)
    {
      CC = CCsave;
      Centroid = Centroidsave;
    }
      N = Nsave;
      max = maxsave;
      min = minsave;
      indmax = indmaxsave;
      indmin = indminsave;
  }

  result = oracleLightFloatWidth(axe,width);
  if(result)
    {
      maxsave = max;
      minsave = min;
      indmaxsave = indmax;
      indminsave = indmin;
      return true;
      }
  else
  {
    if( (Gradient.x==0 && Gradient.y==0) || !allowNewNormal )
    {
       return false;
    }

    doubleCutFloatWidth(axe,width);

  }

  //while at least 1 point left on the search space
  while ( CC.nb() >= 1 )
  {
    CC.centroid(Centroid);  //Centroid.z > 0
    //calls oracle
    result = oracleFloatWidth(axe,width);
    if ( !result )
    {
       if(Gradient.x==0 && Gradient.y==0)
       {
	 return false;
       }
       doubleCutFloatWidth(axe,width);
     }
    else
    {
      CCsave.free();
      CCsave = CC;
      Centroidsave = Centroid;
      Nsave = N;
      maxsave = max;
      minsave = min;
      indmaxsave = indmax;
      indminsave = indmin;

      return true;
    }
   }
   return false;
}


void COBA::free(bool failPrevious)
{
    if(!failPrevious)
      CC.free();
    CCsave.free();
}

void COBA::reInitMinMax(void)
{
  indmax = 0;
  indmin = 0;
  min = N*set[0]; // dot product
  max = min;
}

int majorAxe(S* set, int card)
{
  bool x=true,y=true,z=true;
  for(int i=0;i<card;i++)
    for(int j=i+1;j<card;j++)
    {
      if(set[i].x==set[j].x && set[i].y==set[j].y)
        z=false;
      if(set[i].x==set[j].x && set[i].z==set[j].z)
        y=false;
      if(set[i].y==set[j].y && set[i].z==set[j].z)
        x=false;
    }
  if(!x && !y)
    return 2;
  else if(!x && !z)
    return 1;
  else if(!z && !y)
    return 0;
  else
    return -1;
}

void COBA::printSet(bool forFile)
{
  if(!forFile)
  {
    cout<<"**********ensemble de "<<cardinality<<" voxels : ***************** "<<endl;
    for (int i=0;i<cardinality;i++)
      cout<<"("<<set[i].x<<","<<set[i].y<<","<<set[i].z<<")"<<endl;
    cout<<"************************************************ "<<endl;
  }
  else
  {
    cout<<"**********ensemble de "<<cardinality<<" voxels : ***************** "<<endl;
    for (int i=0;i<cardinality;i++)
      cout<<"voxset += uvoxel( ks, "<<set[i].x<<","<<set[i].y<<","<<set[i].z<<");"<<endl;
    cout<<"************************************************ "<<endl;
  }
}

#pragma package(smart_init)
