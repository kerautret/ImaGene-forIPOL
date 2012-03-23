/** @file CompletedDigitalSurfaceGraph.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : CompletedDigitalSurfaceGraph.h
//
// Creation : 2011/01/13
//
// Version : 2011/01/13
//
// Author : JOL
//
// Summary : Header file for files CompletedDigitalSurfaceGraph.ih and CompletedDigitalSurfaceGraph.cxx
//
// History :
//	2011/01/13 : ?Name? : ?What?
//
// Rcs Id : "@(#)class CompletedDigitalSurfaceGraph declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(CompletedDigitalSurfaceGraph_RECURSES)
#error Recursive header files inclusion detected in CompletedDigitalSurfaceGraph.h
#else // defined(CompletedDigitalSurfaceGraph_RECURSES)
#define CompletedDigitalSurfaceGraph_RECURSES

#if !defined CompletedDigitalSurfaceGraph_h
#define CompletedDigitalSurfaceGraph_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include "ImaGene/digitalnD/KnSpace.h"
#include "ImaGene/digitalnD/KnRCellSet.h"
#include "ImaGene/digitalnD/DigitalSurfaceTracker.h"
#include "ImaGene/planes/TangentialCover.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  struct CompletedDigitalSurfaceGraphMethods;
  /////////////////////////////////////////////////////////////////////////////
  // class CompletedDigitalSurfaceGraph
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'CompletedDigitalSurfaceGraph' <p>
   * Aim: Realizes the concept CDigitalGraph on digital surfaces.
   */
  class CompletedDigitalSurfaceGraph
  {

    // ----------------------- Needed types ------------------------------
  public:
    friend class CompletedDigitalSurfaceGraphMethods;
    /**
     * The type for each vertex.
     */
    typedef Kn_sid Vertex;

    /**
     * The type for the (proper or not) neighborhood of a vertex
     */
    struct Neighborhood : public std::vector<Vertex>
    {
    private:
      Vertex myP;
    public:
      inline Neighborhood() {}
      inline Neighborhood( Vertex p ) : myP( p ) {}
      inline const Vertex & center() const { return myP; }
      inline void setCenter( Vertex p ) { myP = p; }
      // also has begin() and end() services.
    };

    /**
     * The type for an arbitrary large set of vertices.
     */
    typedef KnRCellSet Set;

    typedef CompletedDigitalSurfaceGraphMethods Methods;

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~CompletedDigitalSurfaceGraph();

    /**
     * Constructor. 
     */
    CompletedDigitalSurfaceGraph( const KnSpace & ks, 
			 KnRCellSet surface,
			 const DigitalSurfaceTracker & tracker );

    /**
     * @return the empty set of vertices (subset of this graph).
     */
    Set emptySet() const;

    /**
     * @return the whole set of vertices (subset of this graph).
     */
    Set wholeSet() const;

    /**
     * Sets the neighborhood of [p] in [n].
     *
     * @param p any vertex of this graph.
     * @param n (modified) contains after the neighborhood of p.
     */
    void setNeighborhood( const Vertex & p, Neighborhood & n ) const;

    /**
     * Sets the proper neighborhood of [p] in [n].
     *
     * @param p any vertex of this graph.
     * @param n (modified) contains after the proper neighborhood of p.
     */
    void setProperNeighborhood( const Vertex & p, Neighborhood & n ) const;

    /**
     * @return the dimension of the embedding space.
     */
    int dim() const;

    /**
     * @return the diameter of the graph in its embedding space.
     */
    int diameter() const;

    /**
     * Embeds the vertex [p] in the digital space.
     *
     * @param p any vertex of this graph.
     * @param coords (modified) the integer coordinates of the embedding of [p].
     */
    void embed( const Vertex & p, int* coords ) const;

    /**
       @param v any vertex.
       @return its preferred axis when computing its disk.
    */
    unsigned int preferredAxis( const Vertex & v ) const;

    // ----------------------- Interface --------------------------------------
  public:

    /**
     * Writes/Displays the object on an output stream.
     * @param that_stream the output stream where the object is written.
     */
    void selfDisplay( std::ostream & that_stream ) const;

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool OK() const;
  

    // ------------------------- Datas ----------------------------------------
  private:
    const KnSpace & mySpace;
    KnRCellSet mySurface;
    DigitalSurfaceTracker* myTracker;
    mutable int myDiameter;
    // ------------------------- Hidden services ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    INLINE CompletedDigitalSurfaceGraph();

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE CompletedDigitalSurfaceGraph( const CompletedDigitalSurfaceGraph & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE CompletedDigitalSurfaceGraph & operator=( const CompletedDigitalSurfaceGraph & other );
  
    // ------------------------- Internals ------------------------------------
  private:
  
  };


  struct CompletedDigitalSurfaceGraphMethods
    : public TangentialCoverMethods<CompletedDigitalSurfaceGraph>
  {
    typedef CompletedDigitalSurfaceGraph::Vertex Vertex;
    
    inline CompletedDigitalSurfaceGraphMethods( const CompletedDigitalSurfaceGraph & dsg )
      : myDSG( dsg )
    {}

    inline
    virtual unsigned int preferredAxis( const Vertex & v )
    {
      return myDSG.mySpace.sorthDir( v );
    }

    inline
    virtual void adjustOrientation
    ( const Vertex & v, std::vector<double> & normal )
    {
      //get the normal vector of the surfel
      Vector vAxe = myDSG.mySpace.sorthVectorBasis( v ); 
      //check the sign of the dot product
      double dotproduct = vAxe.ro(0)*normal[0] 
	+ vAxe.ro(1)*normal[1] + vAxe.ro(2)*normal[2];
//       std::cout<<"nor : "<<normal[0] <<" "<<normal[1] <<" "<<normal[2]<<std::endl;
      if ( dotproduct < 0 )
	{
	  normal[0] = -normal[0];
	  normal[1] = -normal[1];
	  normal[2] = -normal[2];
	}
      if ( myDSG.mySpace.sdirect( v, myDSG.mySpace.sorthDir( v ) ) )
	{
	  normal[0] = -normal[0];
	  normal[1] = -normal[1];
	  normal[2] = -normal[2];
	}
    }

    inline
    virtual double projectedArea( const Vertex & v, const double* normal )
    {
      Vector vAxe = myDSG.mySpace.sorthVectorBasis( v ); 
      double dotproduct = vAxe.ro(0)*normal[0] 
	+ vAxe.ro(1)*normal[1] + vAxe.ro(2)*normal[2];
      return dotproduct;
    }


  private:
    const CompletedDigitalSurfaceGraph & myDSG;
  };


  /**
   * Overloads 'operator<<' for displaying objects of class 'CompletedDigitalSurfaceGraph'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'CompletedDigitalSurfaceGraph' to write.
   * @return the output stream after the writing.
   */
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const CompletedDigitalSurfaceGraph & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/planes/CompletedDigitalSurfaceGraph.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CompletedDigitalSurfaceGraph_h

#undef CompletedDigitalSurfaceGraph_RECURSES
#endif // else defined(CompletedDigitalSurfaceGraph_RECURSES)
