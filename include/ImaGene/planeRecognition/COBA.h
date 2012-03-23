//---------------------------------------------------------------------------

#ifndef algo_optiH
#define algo_optiH

#include "ImaGene/planeRecognition/convexe.h"
#include "ImaGene/planeRecognition/Z3.h"
#include "ImaGene/planeRecognition/plane_generator.h"
#include "ImaGene/planeRecognition/utils.h"

class COBA
{
   public :
/************************************************************************************/
S * set; //set of points
int largeurPlan;  //width of the set
int cardinality;  //number of points in the current set

convexe CC; // current search space
S Centroid; // centroid of the current search space, a 2D rational point (a/d,b/d) represented as a 3D integer point (a,b,d)

I g; //digitalization step

S N;  //current normal vector
I max, min;  //current max and min dot products
int indmax, indmin; //current index of points defining the max and min dot products
P Gradient; //current gradient
//int width;

//to come back after failing
convexe CCsave;
S Centroidsave;
S Nsave;
I maxsave, minsave;
int indmaxsave, indminsave;

//int indexLastPoint;  //index of the last added point
//int nbAddedPoints; //number of added points the last time
/************************************************************************************/

// COBA(PlaneGenerator & plan, int D);  //constructor
COBA(S* vSet, int D, int axe);  //constructor

// bool run(int axe, int width, bool failPrevious, bool recordConvex);  //run the recognition algorithm, return a boolean

bool run(int axe, int width, bool failPrevious, bool allowNewNormal);  //run the recognition algorithm, return a boolean
bool runFloatWidth(int axe, float width, bool failPrevious, bool allowNewNormal);  //run the recognition algorithm, return a boolean

bool runInit(int axe, int width, bool failPrevious);  //run the recognition algorithm, return a boolean

void standardInit();  //initialisation of the search space

void free(bool failPrevious);

// void fourPointsInit(S minX, S maxX, S minY, S maxY); //initialisation of the search space using four points

void doubleCut(int axe, int width); //applies a double cut on the search space
void doubleCutFloatWidth(int axe, float width); //applies a double cut on the search space

bool oracle(int axe, int width); //evaluate the thickness function
//computes "idnmax" and "indmin", the index of the points defining the max and the min
//computes the "Gradient" relative to these two points
bool oracleFloatWidth(int axe, float width); //evaluate the thickness function

bool oracleLight(int axe, int width); //evaluate the thickness function using only the last added points and the two ones defining the max and the min dot products
//during the last call of the oracle
bool oracleLightFloatWidth(int axe, float width);

void reInitMinMax(void);

void printSet(bool forFile); //print de set of voxels
};

/*bool run_COBA(int axe);*/

int majorAxe(S* set, int card); //return the major axe (0=x,1=y,2=z) of the set of voxels or -1 if it is not possible to decide

// void printSet(S* set, int card); //print de set of voxels

#endif
