/*!----------------------------------------------------------------------
\file drt_discret_iterator.H
\brief Declaration

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef DRT_DISCRET_ITERATOR_H
#define DRT_DISCRET_ITERATOR_H

#include "drt_discret.H"

namespace DRT
{
namespace UTILS
{

  /// abstract iterator base
  /*!
    The Discretization class does a lot of things. Among those it handles all
    elements and nodes, both the overlapping (column) and non-overlapping
    (row) lists. Sometimes we want to loop over one of those lists. And
    sometimes we have an algorithm that can be used on more that just one of
    those lists. In this case the algorithm can be formulated in terms of an
    abstract DiscretizationIterator.

    \author u.kue
    \date 08/09
   */
  class DiscretizationIterator
  {
  public:

    explicit DiscretizationIterator(const Discretization& dis) : dis_(dis) {}

    /// virtual destruction
    virtual ~DiscretizationIterator() {}

    /// Discretization we look at
    virtual const Discretization& Discret() const { return dis_; }

    /// Number of entries
    virtual int NumEntries() const = 0;

    /// Entry at LID i
    virtual ParObject* Entry(int i) const = 0;

    /// associated map
    virtual const Epetra_Map* Map() const = 0;

  protected:

    /// the actual Discretization
    const Discretization& dis_;
  };


  /// abstract node iterator base
  class DiscretizationNodeIterator : public DiscretizationIterator
  {
  public:
    explicit DiscretizationNodeIterator(const Discretization& dis) : DiscretizationIterator(dis) {}

    virtual Node* Entry(int i) const = 0;
  };


  /// abstract element iterator base
  class DiscretizationElementIterator : public DiscretizationIterator
  {
  public:
    explicit DiscretizationElementIterator(const Discretization& dis) : DiscretizationIterator(dis) {}

    virtual Element* Entry(int i) const = 0;
  };


  /// iteration over row nodes
  class RowNodeIterator : public DiscretizationNodeIterator
  {
  public:
    explicit RowNodeIterator(const Discretization& dis) : DiscretizationNodeIterator(dis) {}

    virtual int NumEntries() const { return dis_.NumMyRowNodes(); }

    virtual Node* Entry(int i) const { return dis_.lRowNode(i); }

    virtual const Epetra_Map* Map() const { return dis_.NodeRowMap(); }
  };


  /// iteration over column nodes
  class ColNodeIterator : public DiscretizationNodeIterator
  {
  public:
    explicit ColNodeIterator(const Discretization& dis) : DiscretizationNodeIterator(dis) {}

    virtual int NumEntries() const { return dis_.NumMyColNodes(); }

    virtual Node* Entry(int i) const { return dis_.lColNode(i); }

    virtual const Epetra_Map* Map() const { return dis_.NodeColMap(); }
  };


  /// iteration over row elements
  class RowElementIterator : public DiscretizationElementIterator
  {
  public:
    explicit RowElementIterator(const Discretization& dis) : DiscretizationElementIterator(dis) {}

    virtual int NumEntries() const { return dis_.NumMyRowElements(); }

    virtual Element* Entry(int i) const { return dis_.lRowElement(i); }

    virtual const Epetra_Map* Map() const { return dis_.ElementRowMap(); }
  };


  /// iteration over column elements
  class ColElementIterator : public DiscretizationElementIterator
  {
  public:
    explicit ColElementIterator(const Discretization& dis) : DiscretizationElementIterator(dis) {}

    virtual int NumEntries() const { return dis_.NumMyColElements(); }

    virtual Element* Entry(int i) const { return dis_.lColElement(i); }

    virtual const Epetra_Map* Map() const { return dis_.ElementColMap(); }
  };
}
}

#endif
