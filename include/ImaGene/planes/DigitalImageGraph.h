/** @file DigitalImageGraph.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : DigitalImageGraph.h
//
// Creation : 2011/01/26
//
// Version : 2011/01/26
//
// Author : JOL
//
// Summary : Header file for files DigitalImageGraph.ih and DigitalImageGraph.cxx
//
// History :
//	2011/01/26 : ?Name? : ?What?
//
// Rcs Id : "@(#)class DigitalImageGraph declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(DigitalImageGraph_RECURSES)
#error Recursive header files inclusion detected in DigitalImageGraph.h
#else // defined(DigitalImageGraph_RECURSES)
#define DigitalImageGraph_RECURSES

#if !defined DigitalImageGraph_h
#define DigitalImageGraph_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <vector>
#include <set>
//#include "ImaGene/base/CowPtr.h"
#include "ImaGene/base/SmartSet.h"
#include "ImaGene/planes/CDigitalGraph.h"
#include "ImaGene/planes/TangentialCover.h"
#include "ImaGene/image/Pixel.h"
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  template <typename T>
  struct DigitalImageGraphMethods;
  
  /////////////////////////////////////////////////////////////////////////////
  // class DigitalImageGraph
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'DigitalImageGraph' <p> Aim: Realizes the
   * concept CDigitalGraph on a digital image. The value of each image
   * pixel acts as a third digital coordinate.
   *
   * @tparam TImage2D any image type. An image has the methods
   * width(), height(), valueRange() (max value - min value ) and a
   * method value( int x, int y ) that returns an int. It is only
   * referenced in this object.
   */
  template <typename TImage2D>
  class DigitalImageGraph
  {
    // ----------------------- public types ------------------------------
  public:
    typedef TImage2D Image2D;
    // struct Pixel { 
    //   int x; 
    //   int y; 
    //   inline Pixel() {}
    //   inline Pixel( const Pixel & other ) : x( other.x ), y( other.y ) {}
    //   inline Pixel( int mx, int my ) : x( mx ), y( my ) {}
    //   inline Pixel & operator=( const Pixel & other )
    //   {
    // 	x = other.x;
    // 	y = other.y;
    // 	return *this;
    //   }
    //   inline bool operator <( const Pixel & other ) const
    //   {
    // 	return ( x < other.x ) || ( ( x == other.x ) && ( y < other.y ) );
    //   }
    //   inline bool operator ==( const Pixel & other ) const
    //   {
    // 	return ( x == other.x ) && ( y == other.y );
    //   }
    //   inline bool operator !=( const Pixel & other ) const
    //   {
    // 	return ( x != other.x ) || ( y != other.y );
    //   }
    //   void selfDisplay( std::ostream & out ) const
    //   {
    // 	out << "(" << x << "," << y << ")";
    //   }
    // };
    

    /**
     * The type for each vertex.
     */
    typedef Pixel Vertex;

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

    // /**
    //  * The type for an arbitrary large set of vertices.
    //  */
    // struct Set
    // {
    //   typedef std::set<Vertex> _Set;
    //   typedef typename _Set::iterator iterator;
    //   typedef typename _Set::const_iterator const_iterator;

    //   inline Set() 
    // 	: mySet( new _Set )
    //   {}
    //   inline ~Set() 
    //   {
    // 	// automatic: delete mySet;
    //   }
    //   inline Set( const Set & other ) 
    // 	: mySet( other.mySet )
    //   {}
    //   inline Set & operator=( const Set & other ) 
    //   {
    // 	if ( this != &other )
    // 	  mySet = other.mySet;
    // 	return *this;
    //   }
    //   inline void clear()
    //   {
    // 	mySet->clear();
    //   }
    //   inline unsigned int size() const
    //   {
    // 	return mySet->size();
    //   }
    //   inline iterator begin()
    //   {
    // 	return mySet->begin();
    //   }
    //   inline iterator end()
    //   {
    // 	return mySet->end();
    //   }
    //   inline const_iterator begin() const
    //   {
    // 	return mySet->begin();
    //   }
    //   inline const_iterator end() const
    //   {
    // 	return mySet->end();
    //   }
    //   inline bool operator[]( const Vertex & vtx ) const
    //   {
    // 	return mySet->find( vtx ) != mySet->end();
    //   }
    //   inline Set& operator+=( const Vertex & vtx )
    //   {
    // 	mySet->insert( vtx );
    // 	return *this;
    //   }
    //   inline Set& operator-=( const Vertex & vtx )
    //   {
    // 	mySet->erase( vtx );
    // 	return *this;
    //   }
    // private:
    //   CowPtr< _Set > mySet;
    // };
    
  public: 
    typedef SmartSet<Vertex> Set;
    typedef DigitalImageGraphMethods<Image2D> Methods;

    // ----------------------- Standard services ------------------------------
  public:

    /**
     * Destructor. 
     */
    ~DigitalImageGraph();

    /**
     * Constructor.
     * The object is not valid.
     */
    DigitalImageGraph();

    /**
     * Constructor with image.
     * @param image any Image2D
     */
    DigitalImageGraph( const Image2D & image,
		       unsigned int quantification = 1,
		       unsigned int xfactor = 1,
		       unsigned int yfactor = 1 );

    /**
     * Init with image.
     * @param image any Image2D
     */
    void init( const Image2D & image,
	       unsigned int quantification = 1,
	       unsigned int xfactor = 1,
	       unsigned int yfactor = 1 );


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
    
    /**
     * a const-pointer to the image defining the digital graph.
     */
    const Image2D* myImage;
    int myDiameter;
    unsigned int myQuantification;
    unsigned int myXfactor;
    unsigned int myYfactor;

    // ------------------------- Hidden services ------------------------------
  protected:

  private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    INLINE DigitalImageGraph( const DigitalImageGraph & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    INLINE DigitalImageGraph & operator=( const DigitalImageGraph & other );
  
    // ------------------------- Internals ------------------------------------
  private:
  
  };


  template <typename T>
  struct DigitalImageGraphMethods
    : public TangentialCoverMethods< DigitalImageGraph<T> >
  {
    typedef DigitalImageGraph<T> _DigitalImageGraph;
    typedef typename _DigitalImageGraph::Vertex Vertex;
    
    inline DigitalImageGraphMethods( const _DigitalImageGraph & dig )
      : myDIG( dig )
    {}

    inline
    virtual unsigned int preferredAxis( const Vertex & v )
    {
      return 2;
    }

    inline
    virtual void adjustOrientation
    ( const Vertex & v, std::vector<double> & normal )
    {
      //get the normal vector of the surfel
      double vAxe[ 3 ] = { 0.0, 0.0, 1.0 };
      //check the sign of the dot product
      double dotproduct = vAxe[0]*normal[0] 
	+ vAxe[1]*normal[1] + vAxe[2]*normal[2];
      if ( dotproduct < 0.0 )
	{
	  normal[0] = -normal[0];
	  normal[1] = -normal[1];
	  normal[2] = -normal[2];
	}
    }

    inline
    virtual double projectedArea( const Vertex & v, const double* normal )
    {
      double vAxe[ 3 ] = { 0.0, 0.0, 1.0 };
      double dotproduct = vAxe[ 0 ]*normal[0] 
       	+ vAxe[ 1 ]*normal[1] + vAxe[ 2 ]*normal[2];
      return 1.0 / dotproduct;
      // return 1.0;
    }

  private:
    const _DigitalImageGraph & myDIG;
  };


  /**
   * Overloads 'operator<<' for displaying objects of class 'DigitalImageGraph'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'DigitalImageGraph' to write.
   * @return the output stream after the writing.
   */
  template <typename TImage2D>
  std::ostream&
  operator<<( std::ostream & that_stream, 
	      const DigitalImageGraph<TImage2D> & that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/planes/DigitalImageGraph.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined DigitalImageGraph_h

#undef DigitalImageGraph_RECURSES
#endif // else defined(DigitalImageGraph_RECURSES)
