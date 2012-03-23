// /*#include "stdafx.h"*/
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "ImaGene/planeRecognition/convexe.h"


using namespace std;
using namespace ImaGene;

void     convexe::free()          { /*cout<<"convex free"<<endl;*/  LP.free(); /*cout<<"convex free done"<<endl;*/        }
int      convexe::nb()            {  return LP.size(); }
listeP & convexe::ListeSommets()  {  return LP;        }

convexe & convexe::operator = (convexe & A)   
{ 
  LP = A.ListeSommets();
  return * this; 
}

void convexe::printVertices()
{
  cout<<"list of vertices:"<<endl;
  for(int i=0;i<nb();i++)
    cout<<"("<<LP[i].x<<","<<LP[i].y<<") ";
  cout<<endl;
}

//add the point K to the convex hull at position "pos"
//if K is not already the first or the last vertex of the hull
//increment "pos"
void convexe::AddAtIndex(int & pos, const P & K)
{
   if(pos == LP.size())
   {
      if(Add(K))
        pos++;
   }
   else if((LP[pos-1]!=K)&&(LP[pos]!=K))
   {
     LP.insert(pos,K);
     pos++;
   }
}

//add the point K at the end of the convex hull
//if K is not already the first or the last vertex of the hull, return true
bool convexe::Add(const P & K)
{
   if(LP.size()==0)
     {  LP.push(K);  return true;}
   else if((LP[LP.size()-1]!=K)&&(LP[0]!=K))
           { LP.push(K);  return true;}
        else { return false;}
}

//compute the area *2
const I & convexe::area()
{
  int i;

  myArea = 0;
  
  for(i=0;i<LP.size()-1;i++)
    myArea += LP[i]^LP[i+1];
  
  myArea += LP[i]^LP[0];
  
  return myArea;
}

//if area is not 0, computes centroid
//else, computes the middle of the straight line segment
//centroid is a 2D rational point but it is represented as a 3D integer point
// (a/d,c/d) corresponds to (a,b,d)
// store the centroid in "point_test"
void convexe::centroid(S & point_test)
{
  int i;
  static I a;
  static I b;
  static I den;
  static I aire;
  a = I_ZERO;
  b = I_ZERO;
  aire = area();

  if( aire > I_ZERO )
    {
      den = I(3)*aire;
      for(i=0;i<LP.size()-1;i++)
	{
	  a += (LP[i].x+LP[i+1].x) * ( LP[i]^LP[i+1]);
	  b += (LP[i].y+LP[i+1].y) * ( LP[i]^LP[i+1]);
	}
      a += (LP[i].x+LP[0].x) * ( LP[i]^LP[0]);
      b += (LP[i].y+LP[0].y) * ( LP[i]^LP[0]);
    }
  else
    {
      den = LP.size();
      for(i=0;i<LP.size();i++)
	{
	  a += LP[i].x;
	  b += LP[i].y;
	}
    }
  point_test = S(a,b,den);
}

// compute the convex hull of grid points satisfying the constraints N1.P<=c1, N2.P<=c2 and N3.P>=c3
// N2.P<=c2 corresponds to the cut
//two parts of computation: from constraint 1 to constraint 3 and from constraint 3 
//to constraint 1
//the computed vertices are stored in "resultup" of size "nbverticesup"
//pointRefC1 and pointRefC3 corresponds to grid point lying on the supporting lines of C1 and of C3 resp
//pos corresponds an index in the list of vertices of the convex, to add the next new vertex

//NB: the method also computes grid point satisfying N1.P<=c1 and N3.P>=c3 but not
// satisfying N2.P<=c2. They are stored in "resultdown" of size "nbverticesdown"
//the algorithm uses these points that's why they appear in the code
void convexe::convexHullBorder(P pointRefC1, P pointRefC3, P N1, I c1, P N2, I c2, P N3, I c3, int & pos)
{
  //cout<<"begin convex hull border!!!"<<endl;
  static P u, v; //vectors u and v to determine the next vertex
  //cout<<"1"<<endl;
  int integerIntersection;
//cout<<"2"<<endl;
  //to store half convex hull border
  static P resultup[1000];
 // P * resultup = new P [1000];
  //cout<<"3 "<<endl;
  static P resultdown[1000];
  //P * resultdown = new P [1000];
  //cout<<"4 "<<endl;
  int nbverticesup, nbverticesdown;
  //cout<<"5"<<endl;
  //first part
  nbverticesup = 1;
  //cout<<"6"<<endl;
  nbverticesdown = 1;
  
  //initializes A, B, u, v and the two first vertices of resultup and resultdown
  //integerIntersection is equal to one when
  //the intersection of the supporting lines of C1
  //and C2 corresponds to an integer point
  //cout<<"7"<<endl;
  //cout<<"before pointRef "<<pointRefC1.x<<" "<<pointRefC1.y<<endl;
  resultup[0] = pointRefC1;
  //cout<<"before init "<<resultup[0].x<<" "<<resultup[0].y<<endl;
  integerIntersection = init(N1,c1,N2,c2,/*u,*/v,resultup, resultdown);
  //cout<<"init done"<<endl;
  //cout<<"after init "<<resultdown[0].x<<" "<<resultdown[0].y<<endl;
  //cout<<"after init "<<v.x<<" "<<v.y<<endl;

  if(integerIntersection!=1) //not integer intersection
  {
    //cout<<"no integer intersection"<<endl;
    //computation of the first part of the border
    ComputeBorderPart1(v,N2,c2,N3,c3,resultup, resultdown,nbverticesup,nbverticesdown);
    //cout<<"border 1 done"<<endl;
  }
  /*else
    cout<<"integer intersection"<<endl;*/
    
   //cout<<"fill in begins"<<endl;

   for(int i=0; i<nbverticesup; i++) //fill in convexup
   {
     //cout<<"before add ("<<resultup[i].x<<","<<resultup[i].y<<") at index "<<pos<<endl; 
     AddAtIndex(pos,resultup[i]);
     //cout<<"done"<<endl;
    }  
   //cout<<"fill in done"<<endl;
   //printVertices();

  //second part
  nbverticesup = 1;
  nbverticesdown = 1;

  //initializes A, B, u, v and the two first vertices of resultup and resultdown
  //integerIntersection is equal to one when
  //the intersection of the supporting lines of C3
  //and C2 corresponds to an integer point
  resultup[0] = pointRefC3;
  integerIntersection = init(N3,c3,N2,c2,/*u,*/v,resultup, resultdown);

  if(integerIntersection!=1) //not integer intersection
  {
    //cout<<"no integer intersection"<<endl;
    //computation of the second part of the border
    ComputeBorderPart1(v,N2,c2,N1,c1,resultup, resultdown,nbverticesup,nbverticesdown);
    //cout<<"border 2 done"<<endl;
  }
  
  for(int i=nbverticesup-1; i>=0; i--) //fill in convexup
    {
      //cout<<"before add ("<<resultup[i].x<<","<<resultup[i].y<<") at index "<<pos<<endl; 
      AddAtIndex(pos,resultup[i]); 
      //cout<<"done"<<endl;
      //printVertices();
    }
 // cout<<"fill in done"<<endl;
 // printVertices();
}

//computes the constraint of the form N.P<=c
//whose supporting line passes through A and B
//such that the points refPoint1 and refPoint2 satisfy the constraint
void computeConstraintFromPoints(P & N, I & c, P A, P B, P refPoint1, P refPoint2)
{
  I GCD;
  //cout<<" begin compute constraint from points"<<endl;
  N.x = A.y - B.y;
  N.y = B.x - A.x;
  c = N.x*A.x + N.y*A.y;
  if((N.x*refPoint1.x + N.y*refPoint1.y > c) || (N.x*refPoint2.x + N.y*refPoint2.y > c))
  {
    N.x = -N.x;
    N.y = -N.y;
    c = -c;
  }
  //cout<<"N=("<<N.x<<","<<N.y<<")"<<endl;
  //simplification of the constraint
  GCD = _pgcd(N.x,N.y);
  //cout<<"gcd="<<GCD<<endl;
  N.x = N.x/GCD;
  N.y = N.y/GCD;
  c = floor(c,GCD);
  //cout<<"compute constraint from points c="<<c<<endl;
}

//cut the convex polygon with the constraint N.P <= c
bool convexe::cutOpti(P N, I c)
{
  //cout<<"begin cut opti"<<endl;
  static bool visible[10000]; //table of visibility of each vertex
  int n = LP.size(); //number of vertices
  //cout<<"nb vertices : "<<n<<endl;
  I aire = area(); //area of the current convex polygon
  int index;

  P N1, N3; //to determine the two constraints supported by vertices of the polygon
  I c1, c3;

  if ( n >= 10000 )
    {cout<<"n est supérieur à 1000!!!!"<<endl; throw;}

  // computation of the visibility of each vertex
  Fori(n)   visible[i] = (LP[i]*N <= c);
  int nbVisible = 0;
  Fori(n)   if ( visible[i] )  nbVisible++;
  //cout<<"visibility checked n="<<n<<" nbVisible="<<nbVisible<<endl;
  if ( nbVisible == n )  {/*cout<<"tous visibles"<<endl;*/ return false;}
  if ( nbVisible == 0 ) {/*cout<<"aucun visibles"<<endl;*/ LP.free(); /*cout<<"cutopti return true, plus de sommets"<<endl;*/ return true; }
 // cout<<"reconstruction begins"<<endl;
  // determine the 2 edges A1B1 and A2B2 intersected by the constaint N.P <= c
  //A1 and A2 are visible (they satisfy the constraint)
  P A1,B1,A2,B2;
  bool A1B1Found = false;
  bool A2B2Found = false;
  Fori(n)  if ( ( visible[i]) && (!visible[(i+1)%n]) )  { A1 = LP[i]; B1 = LP[(i+1)%n]; A1B1Found=true; }
  Fori(n)  if ( (!visible[i]) && ( visible[(i+1)%n]) )  { A2 = LP[(i+1)%n]; B2 = LP[i]; A2B2Found=true; }
 /* if(A1B1Found && A2B2Found)
  {
    cout<<"A1 and A2 are found"<<endl;
    cout<<"A2=("<<A2.x<<","<<A2.y<<")"<<endl;
    cout<<"B2=("<<B2.x<<","<<B2.y<<")"<<endl;
    cout<<"A1=("<<A1.x<<","<<A1.y<<")"<<endl;
    cout<<"B1=("<<B1.x<<","<<B1.y<<")"<<endl;
  }
  else
  {
    cout<<"A1 and A2 are not found"<<endl;
    exit(-1);
  }*/
  //cout<<"then..."<<endl;
  // delete non visible vertices
  Fori(n) 
  {
     if (! visible[i] )   LP.mark(i);
  }
  LP.pack();
  //cout<<"pack done"<<endl;
  //determine the index of the vertex A1
  n = LP.size();
  Fori(n) 
  {
     if(A1 == LP[i]) index=i;
  }
  index++;

  if(aire>I(0)) //convex not reduced to a straight line segment
  {
  //cout<<"aire non nulle"<<endl;
  //computes the constraint C whose supporting line passes through A1 and B1
  computeConstraintFromPoints(N1, c1, A1, B1, A2, B2);
  //computes the constraint C whose supporting line passes through A2 and B2
  computeConstraintFromPoints(N3, c3, A2, B2, A1, B1);
  /*cout<<"constraints computed"<<endl;
  cout<<"N3=("<<N3.x<<","<<N3.y<<")"<<endl;
  cout<<"c3="<<c3<<endl;
  cout<<"A2=("<<A2.x<<","<<A2.y<<")"<<endl;
  cout<<"B2=("<<B2.x<<","<<B2.y<<")"<<endl;
  cout<<"A1=("<<A1.x<<","<<A1.y<<")"<<endl;
  cout<<"B1=("<<B1.x<<","<<B1.y<<")"<<endl;*/
  //run the reconstruction of the convex polygon
  //printVertices();
  convexHullBorder(A1, A2, N1, c1, N, c, N3, c3, index);
 // cout<<"convex hull border done"<<endl;
  //printVertices();
  }
  else //convex reduced to a straight line segment
  {
  //cout<<"aire nulle"<<endl;
  //compute the new extremity of the straight line segment
  P periode = B1-A1;
  reduc(periode);
  I k = (c-N*A1)/(N*periode);
  if(k*N*periode==c-N*A1)
	k=k-I(1);
  P K = A1+k*periode;
  AddAtIndex(index,K);
  }
  //cout<<"reconstruction done"<<endl;
  return true;
}

  //initializes A, B, u, v and the two first vertices of resultup and resultdown
//A lying on C1, satisfying C2 and the closest to C2
//B lying on C1, not satisfying C2 and the closest to C2
//v: valid Bezout vector of the direction vector of C1
//resultup and resultdown: vertices of the CH of grid satisfying and not satisfyin C2 resp.
  //return 1 when
  //the intersection of the supporting lines of C1
  //and C2 corresponds to an integer point
int init(P N1, I c1, P N2, I c2,P & v,P* resultup, P* resultdown)
{
//cout<<"begin init"<<endl;
  I fl;
  I ce;
  int result;
  P directionVector;

  //initialize  vector directionVector (not definitive)
  directionVector.y = N1.x;
  directionVector.x = I(0)-N1.y;

  //compute the intersection of ray (resultup[0],directionVector) with constraint C2
 // cout<<"before coef intersect"<<endl;
  coefficientIntersection(resultup[0],directionVector,N2,c2,fl,ce);
  //cout<<"after coef intersect"<<endl;
  //uses the intersection to compute the first vertex of the upper convex hull
  //i.e. the grid point closest to C2 and satisfying C2
  if(dot_product(resultup[0],N2)>c2 && fl==I(0))
  {
    resultup[0] = resultup[0]+directionVector*(ce);
    directionVector = P(0,0)-(directionVector);
  }
  else if( (dot_product(resultup[0],N2)<=c2 && fl>=I(0)) || 
                                 (dot_product(resultup[0],N2)>c2 && fl<=I(0)))
    resultup[0] = resultup[0]+directionVector*(fl);
  else
  {
    resultup[0] = resultup[0]+directionVector*(ce);
    directionVector = P(0,0)-(directionVector);
  }
  //compute the first vertex of the lower convex hull
  if(fl==ce)  //integer intersection
   {
     resultdown[0] = resultup[0];
     result = 1;
   }
  else
  {
    resultdown[0]=resultup[0]+directionVector;
    //initialization of v: valid Bezout vector of u
    validBezout(resultup[0], directionVector, v, N1, c1, N2, c2, 1);
    result = 0;
  }
//cout<<"end init"<<endl;
  return result;
}

//compute the border of the upper and of the lower convex hull from the points A (up)
// and B (down),
//along the constraint N2.p<=c2 while the vertices satisfy the constraint N3.p<=c3
//the vertices of the two borders are stored in resultup and resultdown
//we also store the number of the vertices of each convex hull in nbverticesup and
//nbverticesdown
void ComputeBorderPart1(/*P & u,*/ P initialBezoutVector, P N2, I c2,P N3, I c3,
                                            P* resultup,P* resultdown,int & nbverticesup,int & nbverticesdown)
{
  P A, B; //represent the last computed vertex above and below the constraint resp.
  I fl,ce;
  P u, v=initialBezoutVector;  //u and v vector (cf paper)
  
  //initialize A and B
  A = resultup[0];
  B = resultdown[0];

  //while A and B do not lie on the supporting line of C2
  //and satisfies C3 and while v is not parallel to C2
  while(((dot_product(A,N2) != c2)&&(dot_product(A,N3) <= c3))
        &&
       ((dot_product(B,N2) != c2)&&(dot_product(B,N3) <= c3)) 
        &&
       (dot_product(v,N2)!=I(0)))
  {
    if(dot_product(v,N2)<I(0)) //second configuration, we find a new B
    {
      //computation of the new vertex
      coefficientIntersection(B, v, N2,c2,fl,ce);
      B = B+v*(fl);
      nbverticesdown++;
      resultdown[nbverticesdown-1] = B;
    }
   else //first configuration, we find a new A
   { 
      //computation of the new vertex
      coefficientIntersection(A, v, N2, c2,fl,ce);
      A = A+v*(fl);
      nbverticesup++;
      resultup[nbverticesup-1] = A;
   }
   //update u and v
   u=B-A;
   validBezout(A, u, v, N2, c2, N2, c2, 0);
  }
//when the loop finishes, we have to complete the computation
//of each convex hull

  if(dot_product(A,N3) > c3) //A does not satisfy C3, we remove it.
      nbverticesup--;
  else if(dot_product(B,N3) > c3) //B does not satisfy C3, we remove it
            nbverticesdown--;
       else if(dot_product(A,N2) == c2) //A lies on C2 so A is also a vertex of
                                        //resultdown
            {
              nbverticesdown++;
              resultdown[nbverticesdown-1] = A;
            }
            else if(dot_product(B,N2) == c2) //B lies on C2 so B is also a vertex of
                                             //resultup
                 {
                   nbverticesup++;
                   resultup[nbverticesup-1] = B;
                 }
}

#pragma package(smart_init)
