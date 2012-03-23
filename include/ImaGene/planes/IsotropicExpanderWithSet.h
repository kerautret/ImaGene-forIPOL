/** @file IsotropicExpanderWithSet.h */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File name : IsotropicExpanderWithSet.h
//
// Creation : 2011/01/20
//
// Version : 2011/01/20
//
// Author : JOL
//
// Summary : Header file for files IsotropicExpanderWithSet.ih and IsotropicExpanderWithSet.cxx
//
// History :
//	2011/01/20 : ?Name? : ?What?
//
// Rcs Id : "@(#)class IsotropicExpanderWithSet declaration."
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#if defined(IsotropicExpanderWithSet_RECURSES)
#error Recursive header files inclusion detected in IsotropicExpanderWithSet.h
#else // defined(IsotropicExpanderWithSet_RECURSES)
#define IsotropicExpanderWithSet_RECURSES

#if !defined IsotropicExpanderWithSet_h
#define IsotropicExpanderWithSet_h

//////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <set>
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Forces inline if nothing is provided by the compiler.
#ifndef INLINE
#define INLINE inline
#endif

namespace ImaGene 
{
  /////////////////////////////////////////////////////////////////////////////
  // class WeightedCell
  /////////////////////////////////////////////////////////////////////////////
  /**
   * Represents an arbitray cell weighted by some double.
   */
  template <typename Cell>
  struct WeightedCell {
    Cell cell;
    double weight;
    inline WeightedCell() {}
    inline WeightedCell( const Cell & c, double w )
      : cell( c ), weight( w ) {}
    /**
     * highest weight (highest priority) is smallest value.
     */
    inline bool operator<( const WeightedCell & other ) const
    {
      return other.weight < weight;
    }
  };

  /**
   * Represents an arbitray cell weighted by some double, and storing
   * another double.
   */
  template <typename Cell>
  struct Weighted2Cell {
    Cell cell;
    double weight;
    double weight2;
    inline Weighted2Cell() {}
    inline Weighted2Cell( const Cell & c, double w, double w2 )
      : cell( c ), weight( w ), weight2( w2 ) {}
    /**
     * highest weight (highest priority) is smallest value.
     */
    inline bool operator<( const Weighted2Cell & other ) const
    {
      return other.weight < weight;
    }
  };

  /**
   * Model of concept CDigitalDistance.
   */
  struct EuclideanDigitalDistance {
    unsigned int myDim;
    inline EuclideanDigitalDistance()
      : myDim( 0 )
    {}
    inline EuclideanDigitalDistance( unsigned int dim )
      : myDim( dim )
    {}
    inline EuclideanDigitalDistance( const EuclideanDigitalDistance & other )
      : myDim( other.myDim )
    {}
    inline EuclideanDigitalDistance &
    operator=( const EuclideanDigitalDistance & other )
    {
      myDim = other.myDim;
    }
    /**
     * Computes the Euclidean distance between lattice points p1 and p2.
     *
     * @param p1 an array of integers of correct dimension.
     * @param p2 an array of integers of correct dimension.
     *
     * @return their euclidean distance.
     */
    double distance( const int* p1, const int* p2 ) const;

    /**
     * Computes the squared Euclidean distance between lattice points
     * p1 and p2. Faster than distance. Same order as distance.
     *
     * @param p1 an array of integers of correct dimension.
     * @param p2 an array of integers of correct dimension.
     *
     * @return their euclidean distance.
     */
    double squaredDistance( const int* p1, const int* p2 ) const;

  };

  /////////////////////////////////////////////////////////////////////////////
  // class IsotropicExpanderWithSet
  /////////////////////////////////////////////////////////////////////////////
  /** 
   * Description of class 'IsotropicExpanderWithSet' <p> Aim: This is
   * a helper class to model an isotropic expansion onto a Digital
   * Graph (see concept CDigitalGraph). It is used to compute nu-thick
   * disks and extensions. Compared to IsotropicExpander, this class
   * uses a std::set representation to stored the visited
   * vertices. Avoids the creation and deletion of DigitalGraph::Set.
   *
   *
   * @tparam TDigitalGraph the type of digital graph.
   *
   * @tparam TDigitalDistance the type of the object used to compute
   * distance between vertices of the digital graph.
   *
   * @see IsotropicExpander
   * @see NuThickDisk
   * @see TangentialCoverHierarchy
   */
  template <typename TDigitalGraph, 
	    typename TDigitalDistance >
  class IsotropicExpanderWithSet
  {

    // ----------------------- Standard services ------------------------------
  public:
    typedef TDigitalGraph DigitalGraph;
    typedef TDigitalDistance DigitalDistance;
    typedef typename DigitalGraph::Vertex Vertex;
    typedef typename DigitalGraph::Set Set;
    typedef typename DigitalGraph::Neighborhood Neighborhood;
    typedef typename Neighborhood::const_iterator NeighborhoodConstIterator;
    typedef WeightedCell<Vertex> WeightedVertex;
    typedef std::priority_queue< WeightedVertex > PriorityQueue;
    typedef std::set<Vertex> MarkedSet;

    /**
     * Destructor. 
     */
    ~IsotropicExpanderWithSet();

    /**
     * Constructor.
     * The object is not valid.
     */
    IsotropicExpanderWithSet();

    /**
     * Copy constructor.
     * @param other the object to clone.
     */
    IsotropicExpanderWithSet( const IsotropicExpanderWithSet & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     */
    IsotropicExpanderWithSet & operator=( const IsotropicExpanderWithSet & other );

    /**
     * Resets the object and clears memory.
     */
    void clear();

    /**
     * Initialization from digital graph and vertex.
     *
     * @param dgraph the digital graph onto which the nu-thick disk
     * are extracted.
     *
     * @param p the center vertex.
     *
     * @param distance the object used for computing distances.
     */
    void init( const DigitalGraph & dgraph, 
	       const Vertex & p,
	       const DigitalDistance & distance );

    /**
     * Reinitialization from a new point [p]. Makes sense only if this
     * new point belongs to the current expansion. After this
     * reinitialization the object starts an isotropic expansion with
     * center [p], but does not output again the points previously
     * outputed.
     */
    void reinit( const Vertex & p );

    // ----------------------- Expansion services ------------------------------
  public:
    /**
     * @return the current set of explored vertices.
     */
    const MarkedSet & exploredVertices() const;

    /**
     * @return the current set of visited vertices (a subset of
     * explored vertices; excludes the explored vertices yet to be
     * visited).
     */
    MarkedSet visitedVertices() const;

    /**
     * @return 'true' if the expander has finished its expansion.
     */
    bool finished() const;

    /**
     * @return a const reference on the next vertex that will be
     * touched by the expansion.
     *
     * NB: valid only if not finished;
     */
    const WeightedVertex & next() const; 

    /**
     * Ignore next element in the expansion. Its neighbors are no more
     * considered in the expansion.
     *
     * @param wv (modified) the next vertex.
     */
    void ignoreNext();

    /**
     * Expands the next vertex. Valid only if not finished().
     */
    void expandNext();

    /**
     * Expands the next vertex and outputs it. Valid only if not finished().
     * @param wv (modified) the next vertex.
     */
    void expandNext( WeightedVertex & wvtx );

    /**
     * Adds again \b already \b marked vertices into the queue.
     * @tparam ForwardConstIterator the type of the iterator on WeightedVertex.
     * @param itb an iterator on the first vertex.
     * @param ite an iterator on the vertex after the last.
     */
    template <typename ForwardConstIterator>
    void pushAgain( ForwardConstIterator itb, ForwardConstIterator ite );
    
    
    // void valideExpandCurrentSet();

    
    /**
     * Expands and outputs all vertices such that the last vertex is
     * the last one having a distance no greater than [d].
     *
     * @param d all output vertices have distance <= d.
     * @param it an output iterator such that *it has type WeightedVertex.
     *
     * @return 'true' iff the expander is not finished.
     */
    template <typename TOutputIterator>
    bool expandNoGreaterThan( double d, TOutputIterator it );

    /**
     * Expands and outputs all vertices such that the last vertex is
     * the last one having a distance smaller than [d].
     *
     * @param d all output vertices have distance < d.
     * @param it an output iterator such that *it has type WeightedVertex.
     *
     * @return 'true' iff the expander is not finished.
     */
    template <typename TOutputIterator>
    bool expandSmallerThan( double d, TOutputIterator it );

    // ----------------------- Distance services ------------------------------
  public:


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
    const DigitalGraph* myDGraph; /**< the digital graph. */
    const DigitalDistance* myDistance; /**< the digital distance. */
    MarkedSet myMark;              /**< used for marking explored vertices. */
    PriorityQueue myQueue;        /**< vertices in the expanding ring. */
    Vertex myVtx;                 /**< Starting vertex. */
    int* myPt;                    /**< Embedding of starting vertex. */
    int* myPt2;                   /**< Embedding of other vertices. */

    // ------------------------- Hidden services ------------------------------
  protected:


  private:


  
    // ------------------------- Internals ------------------------------------
  private:
  
  };

  /**
   * Overloads 'operator<<' for displaying objects of class 'IsotropicExpanderWithSet'.
   * @param that_stream the output stream where the object is written.
   * @param that_object_to_display the object of class 'IsotropicExpanderWithSet' to write.
   * @return the output stream after the writing.
   */
  template <typename TDigitalGraph, typename TDigitalDistance >
  INLINE std::ostream&
  operator<<( std::ostream & that_stream, 
	      const IsotropicExpanderWithSet<TDigitalGraph,TDigitalDistance> & 
	      that_object_to_display );

  
} // namespace ImaGene


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions/methods if necessary.
#if defined(INLINE)
#include "ImaGene/planes/IsotropicExpanderWithSet.ih"
#endif

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined IsotropicExpanderWithSet_h

#undef IsotropicExpanderWithSet_RECURSES
#endif // else defined(IsotropicExpanderWithSet_RECURSES)
