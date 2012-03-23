//---------------------------------------------------------------------------
#ifndef listePH
#define listePH

#include "ImaGene/planeRecognition/Z2.h"
#include<iostream>

class listeP
{
public :

       listeP();
       listeP & operator = (const listeP & ll);
      ~listeP();
       listeP(const listeP & ll);

       void push(const P & P);                       // adds an element at the end
       void push(const listeP & liste);              // adds the elements of liste at the
                                                     //end
       void insert(int pos,const P & P);             // inserts an element at location pos
       int  stock(const P & P );                     // inserts if not present, return pos
       void resize(const int taille);                // resizes
       void suppr_end();                             // deletes the last element
       void mark(const int n);                       // marks the n-th element
       void pack();                                  // deletes all marked elements

       int  size();                                  // returns the size
       P *  adr() const;                             // returns the address
       P    pop_back();                              // deletes and returns the last element
       P &  operator [] (const int n) const;         // compatibility with tables
       bool find(const P & P) const;                 // finds the firsty occurence of P
       void free();                                  // frees the list
       
       

private:
  P * L;                               // address
  int nnb;                              // number of elements
  int reserv;                           // reserved space
  void expand();                        // control memory
  void transport(P * dst, P * src, int nb); // copy
  int * f;                             // table of marks
  int flag_mark;                        // == 1 if table of marks exists
};

 
 
#endif




