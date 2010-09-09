#ifdef STKADAPTIVE

#ifndef STK_ITERATOR_H
#define STK_ITERATOR_H

#include <set>
#include <vector>
#include <iterator>
#include <iostream>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/fem/FieldTraits.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyTraits.hpp>

#include "../stk_refine/stk_types.H"

namespace STK
{

  class ElementIterator // : std::iterator<std::forward_iterator_tag,stk::mesh::Element>
  {
  public:

    ElementIterator( stk::mesh::BulkData & bulk, const stk::mesh::Selector & selector );

    bool done()
    {
      return i_ == elements_.end();
    }

    void next()
    {
      ++j_;
      if ( j_ == bucket_->end() )
      {
        ++i_;
        next_bucket();
      }
    }

    ElementIterator & operator++()
    {
      next();
      return *this;
    }

    stk::mesh::Entity & element()
    {
      return *j_;
    }

    const CellTopologyData & topology()
    {
      return *celltopology_;
    }

    const stk::mesh::Selector & selector() { return selector_; }

  private:

    void next_bucket();

    stk::mesh::BulkData & bulk_;
    stk::mesh::Selector selector_;
    const std::vector<stk::mesh::Bucket*> & elements_;
    std::vector<stk::mesh::Bucket*>::const_iterator i_;
    stk::mesh::Bucket * bucket_;
    stk::mesh::Bucket::iterator j_;

    const CellTopologyData * celltopology_;
  };


  class SideIterator
  {
  public:

    SideIterator( ElementIterator & ei );

    bool done();

    void next();

    SideIterator & operator++()
    {
      next();
      return *this;
    }

    /// @name Query methods

    bool is_boundary_side() { return not hanging_side_ and side_elements_.size()==0; }

    bool is_inner_side() { return not hanging_side_ and side_elements_.size()==1; }

    bool is_big_hanging_side() { return hanging_side_ and side_elements_.size()>1; }

    bool is_small_hanging_side() { return hanging_side_ and side_elements_.size()==1; }

    //@}

    /// @name Original side

    stk::mesh::Entity & element() { return ei_.element(); }

    unsigned side_number() { return side_num_; }

    const CellTopologyData & topology() { return ei_.topology(); }

    //@}

    /// @name Other side

    stk::mesh::Entity & other_element() { return **side_elements_.begin(); }

    unsigned other_side_number() { return *other_side_num_.begin(); }

    const CellTopologyData & other_topology() { return **other_topology_.begin(); }

    //@}

    /// @name Other side hanging node version

    const std::vector<unsigned> & other_side_numbers() { return other_side_num_; }

    const EntityVector & side_elements() { return side_elements_; }

    //@}

    //const CellTopologyData & side_topology() { return *side_topology_; }

  private:

    void new_side();

    void find_side_elements( const EntitySet & side_nodes, EntitySet & side_elements );

    void find_other_element_side( stk::mesh::Entity & oe );

    void find_other_hanging_sides();

    bool contains_side( const CellTopologyData_Subcell & side,
                        const stk::mesh::PairIterRelation & nodes,
                        const EntitySet & side_nodes );

    ElementIterator & ei_;
    stk::mesh::PairIterRelation nodes_;
    unsigned side_dim_;
    unsigned side_num_;

    const CellTopologyData_Subcell * side_;
    const CellTopologyData * side_topology_;

    EntitySet side_nodes_;
    EntityVector side_elements_;

    std::vector<unsigned> other_side_num_;
    std::vector<const CellTopologyData*> other_topology_;

    bool hanging_side_;
  };


  class UniqueSideIterator : public SideIterator
  {
  public:

    UniqueSideIterator( ElementIterator & ei ) : SideIterator( ei )
    {
      skip();
    }

    void next()
    {
      SideIterator::next();
      skip();
    }

    UniqueSideIterator & operator++()
    {
      next();
      return *this;
    }

  private:

    void skip()
    {
      while ( not done() )
      {
        if ( is_boundary_side() )
          break;
        if ( is_inner_side() and element().key()<other_element().key() )
          break;
        if ( is_small_hanging_side() )
          break;
        SideIterator::next();
      }
    }
  };


  void PrintMesh( std::ostream & stream, stk::mesh::BulkData & bulk, stk::mesh::Part & active );
}

#endif
#endif
