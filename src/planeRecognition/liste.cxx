// #include "stdafx.h"
// #include "stdafx.h"
#include "ImaGene/planeRecognition/liste.h"
#include <iostream>

#include <stdio.h>
#include <stdlib.h>

using namespace std;


P &  listeP::operator [] (int n) const  {  return L[n];              }
P    listeP::pop_back()                 {  nnb-- ;   return L[nnb];  }
void listeP::suppr_end()                {  nnb--;                    }
P *  listeP::adr() const                {  return L;                 }
int  listeP::size()                     {  return nnb;               }
void listeP::free()                     {  /*cout<<"listeP free"<<endl;*/  
                                           delete [] L;
                                           if ( flag_mark == 1)
                                              {delete [] f; flag_mark=0;}
					   nnb=0; 
					   reserv=0;
					   /*cout<<"free done"<<endl;*/}
                                      /*{ cout<<"listeP free"<<endl; 
                                          if(nnb!=0)
					    delete [] L; 
					  else
					    cout<<"L is already null"<<endl;
					   nnb = 0; cout<<"free done"<<endl; }*/

listeP & listeP::operator = (const listeP & ll)
{
  nnb = ll.nnb;      // number of elements
  reserv = ll.reserv;  // reserved space
  L = new P[reserv];
  for (int i = 0 ; i < nnb ; i++) L[i] = ll.L[i];                      
  flag_mark = ll.flag_mark;                        // == 1 if table of marks exists
  if(flag_mark==1)
  {
   // cout<<"flag_mark vaut 1!!!!!"<<endl;
    f = new int[nnb];
    for (int i = 0 ; i < nnb ; i++) f[i] = ll.f[i];
  }
  return * this;
}

void listeP::transport( P * dst, P * src , int nb)
{
  for (int i = 0 ; i < nb ; i++ )
     dst[i] = src[i];
}

///////////////////////    PUSH     ////////////////////////////

void listeP::push(const listeP & liste)
{
 // cout<<"begin listeP push with list"<<endl;
  int dep = nnb;
  int size = liste.nnb;
  P * adr = liste.adr();
 // cout<<"adresse recuperee"<<endl;
  resize(nnb+size);
  //cout<<"resize done"<<endl;
  for (int i = 0 ; i < size ; i++)
  {
    L[i+dep] = adr[i];
  }
  //cout<<"copy done"<<endl;
}

////////////////////////////   listeP()     ////////////////////////////

listeP::listeP()
{
  //cout<<"constructeur listeP"<<endl;
  nnb = 0;
  reserv = 5;
  L = new P[reserv];
  flag_mark = 0;
}

//////////////////////////     find      //////////////////////////////
 
bool listeP::find(const P & PP) const
{
   for ( int i = 0 ; ( i < nnb ); i++ )
   if ( *(L+i) == PP )
     return true;

   return false;
}

////////////////////////        MARK       ////////////////////////////
      // NE PAS UTILISER DES INS APPEND PUSH_BACK PENDANT MARK //

void listeP::mark(int n)
{
  if (flag_mark == 0)
  {
    f = new int[nnb];
    flag_mark = 1;
    for (int i = 0 ; i < nnb ; i++) f[i] = 0;
  }
  f[n] += 1;
}

////////////////////////        PACK       //////////////////////////////
    
void listeP::pack()
{
  if ( flag_mark == 0)
    return;
  int deb = 0;
  for (int i = 0 ; i < nnb ; i++)
  {
    if ( f[i] == 0 )
    {
      L[deb] = L[i];
      deb++;
    }
  }
  nnb = deb;
  delete [] f;
  flag_mark = 0;
}

//////////////////////////////      RESIZE      ///////////////////////////////

void listeP::resize(const int taille)
{
  //cout<<"resize!!!!"<<endl;
  //exit(-1);
  if ( taille > reserv )
  {
    P * LL = new P[taille];
    transport(LL,L,nnb);
    nnb = taille;
    reserv = taille;
    delete [] L;
    L = LL;  // a voir
    delete [] LL;
  }
  else
  {
    nnb = taille;
  }
}

////////////////////  PUSH ELEMENT   //////////////////////

void listeP::push(const P & PP)
{
  //cout<<"begin push"<<endl;
  expand();
 // cout<<"expand done"<<endl;
  L[nnb] = PP;
  //cout<<"point added"<<endl;
  nnb++;
}

//////////////////////    INSERT      //////////////////////

void listeP::insert(const int pos,const P & PP)
{
   // cout<<"begin insert"<<endl;
    expand();
    //cout<<"expand done pos = "<<pos<<" and nnb = "<<nnb<<endl;
    for ( int i = nnb ; i > pos ; i--)
      {
      //cout<<"copy ("<<L[i-1].x<<","<<L[i-1].y<<")"<<endl;
      L[i] = L[i-1]; 
      //cout<<"ok"<<endl;
      }
//       {cout<<"copy ("<<L[i-1].x<<","<<L[i-1].y<<")"<<endl; L[i].x = L[i-1].x; cout<<"entre"<<endl; L[i].y = L[i-1].y;}
    //cout<<"copy done"<<endl;
    L[pos] = PP;
    //cout<<"last done"<<endl;
    nnb++;
}

//////////////////////////         ~        ///////////////////////////


listeP::~listeP()
{
  //cout<<"DESTRUCTEUR LISTEP!!!!!!!!!!!!"<<endl;
  //exit(-1);
 /* nnb=0;
  delete [] L;
  cout<<"destruction done"<<endl;
  if ( flag_mark)
    delete [] f;*/
}

//////////////////////////       expand       /////////////////////////


void listeP::expand()
{
 // cout<<"expand nnb = "<<nnb<<" and reserv = "<<reserv<<endl;
//   if ( nnb == reserv )
  if ( nnb >= reserv )
  {
    reserv = reserv * 2;
    P * LL = new P[reserv];
    transport(LL,L,nnb);
    delete [] L;
    L=LL;
    //delete [] LL;
   //cout<<"more space"<<endl;
  }
 //cout<<"no need more space"<<endl;
}
