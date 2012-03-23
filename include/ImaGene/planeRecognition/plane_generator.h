//---------------------------------------------------------------------------

#ifndef Plane_generatorH
#define Plane_generatorH
#include "ImaGene/planeRecognition/Z3.h"
#include "ImaGene/planeRecognition/utils.h"

// to create a piece of digital plane
 //we suppose Diameter odd
class PlaneGenerator
{
  private:
  // bool mem;
   int D;
   void InitXY();
   void AnotherNormal(S & N);

  public :
  
    ~PlaneGenerator();
    PlaneGenerator(int DD);
   
    S Random();             // creates a piece of digital plane and returns the normal vector
    void Reorganize();      //reorganize the set of points in spiral (first point in the middle)

    int size();  //return nb points in the set

    S * L;  //tab of points

    S minX, maxX, minY, maxY; /*to recall of 4 characteristic points of the digital plane 
- such that minX and MaxX have the same y-component and 
- such that minY and MaxY have the same x-component
we choose minX=minY*/

};


#endif
