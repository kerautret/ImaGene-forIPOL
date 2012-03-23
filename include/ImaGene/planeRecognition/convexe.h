#ifndef convexeH
#define convexeH

#include "ImaGene/planeRecognition/liste.h"
#include "ImaGene/planeRecognition/Z3.h"
#include "ImaGene/planeRecognition/utils.h"

// convex polygon
// vertices in counterclockwise

class convexe
{
   public :

        void free();                     // deletes all the vertices
        int  nb() ;                      // returns the number of vertices
        listeP & ListeSommets();         // returns the list of the vertices
	
	
	convexe & operator = (convexe & A);

       //add the point K to the convex hull at position "pos"
       //if K is not already the first or the last vertex of the hull
       //increment "pos"
        void AddAtIndex(int & pos, const P & K);

        //add the point K at the end of the convex hull
        //if K is not already the first or the last vertex of the hull, return true
        bool Add(const P & K);

        // computes area*2
	const I & area();
	
        //if area is not 0, computes centroid
        //else, computes the middle of the straight line segment
        //centroid is a 2D rational point but it is represented as a 3D integer point
        // (a/d,c/d) corresponds to (a,b,d)
        // store the centroid in "point_test"
	void centroid(S & point_test);

        // cuts the polygon with the constraint N.P<=c
        bool cutOpti(P N,I c);

        // compute the convex hull of grid points satisfying the constraints N1.P<=c1, N2.P<=c2 and N3.P>=c3
        // N2.P<=c2 corresponds to the cut
        //two parts of computation: from constraint 1 to constraint 3 and from constraint 3 
        //to constraint 1
        //the computed vertices are stored in "resultup" of size "nbverticesup"

        //NB: the method also computes grid point satisfying N1.P<=c1 and N3.P>=c3 but not
        // satisfying N2.P<=c2. They are stored in "resultdown" of size "nbverticesdown"
        //the algorithm uses these points that's why they appear in the code
        void convexHullBorder(P pointRefC1, P pointRefC3, P N1, I c1, P N2, I c2, P N3, I c3, int & pos);
	
	void printVertices();

 protected :
	listeP LP;
 private:
	I myArea;
};
        //computes the constraint of the form N.P<=c
        //whose supporting line passes through A and B
        //such that the points refPoint1 and refPoint2 satisfy the constraint
        void computeConstraintFromPoints(P & N, I & c, P A, P B, P refPoint1, P refPoint2);

        //initializes A, B, u, v and the two first vertices of resultup and resultdown
        //A lying on C1, satisfying C2 and the closest to C2
        //B lying on C1, not satisfying C2 and the closest to C2
        //v: valid Bezout vector of the direction vector of C1
        //resultup and resultdown: vertices of the CH of grid satisfying and not satisfyin C2 resp.
        //return one when
        //the intersection of the supporting lines of C1
        //and C2 corresponds to an integer point
        int init(P N1, I c1, P N2, I c2,/*P & A, P & B,*//*P & u,*/P & v,P* resultup, P* resultdown);

        //compute the border of the upper and of the lower convex hull from the points A (up)
        // and B (down),
        //along the constraint N2.p<=c2 while the vertices satisfy the constraint N3.p<=c3
        //the vertices of the two borders are stored in resultup and resultdown
        //we also store the number of the vertices of each convex hull in nbverticesup and
        //nbverticesdown
        void ComputeBorderPart1(/*P & A,P & B, P & u, */P initialBezoutVector, P N2, I c2,P N3, I c3, P* resultup,P* resultdown,int & nbverticesup,int & nbverticesdown);

#endif
