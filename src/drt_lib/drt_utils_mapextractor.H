#ifdef CCADISCRET
#ifndef DRT_UTILS_MAPEXTRACTOR_H
#define DRT_UTILS_MAPEXTRACTOR_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <algorithm>

#include <Teuchos_RCP.hpp>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "drt_discret.H"


namespace DRT
{

  /// collection of useful and important but less central stuff
  namespace UTILS
  {

    /// Split a dof row map in two and establish the communication pattern between those maps
    /*!

    The MapExtractor splits the dofs of a Discretization in two parts. Each
    part constitutes a Epetra_Map, a submap of the Discretization's dof row
    map. The communication pattern between the dof row map and the submaps is
    established.

    Examples of such splits include the velocity -- pressure split of the dof
    row map of a fluid problem or the interface -- interior split in FSI
    problems. Many more examples are possible. This is the class to use each
    time a submap needs to be build.

    \note We work on row maps. The maps we deal with are meant to be
    nonoverlapping.

    At the core there are the CondDofMap(), the map of all selected dofs, and
    OtherDofMap(), the map of all remaining dofs. This duality also exists in
    the extraction methods ExtractCondVector() and ExtractOtherVector(), that
    extract a subvector from a full one, and the insertion methods
    InsertCondVector() and InsertOtherVector(), that copy a subvector into a
    full vector. These extractions and insertions are termed communications,
    because internally an Epetra_Import class is used, even though there is no
    communication required once the Epetra_Import object is created.

    Such splits and the need to communicate between maps occur at many
    different places. However, the splitting criteria is very different in
    most cases and cannot be implemented once and for all. This is why we use
    a templated Setup() method, that is parametrized by an ExtractorCond
    class. The ExtractorCond class is supplied by the user and provides the
    decision for each dof. It is very easy to supply new ExtractorCond
    classes, all that is needed is an appropriate call operator. No
    inheritance is needed, as we use templates to get speed and polymorphic
    behavior. However, for most cases there are predefined ExtractorCond
    classes. In particular the class ExtractorCondInCondition selects all dofs
    that belong to nodes from a given DRT::Condition.

    \see ExtractorCondInCondition, ExtractorCondMaxPos, ExtractorCondAnd, ExtractorCondOr

    \author u.kue
    \date 01/08
    */
    class MapExtractor
    {
    public:

      MapExtractor();

      /** \name Setup */
      //@{

      /// generic setup based on user supplied condition class
      template<class C>
      void Setup(const Teuchos::RCP<DRT::Discretization> dis, const C& condition);

      /// setup from a known condition map, a submap of the dof row map
      void SetupMaps(Teuchos::RCP<const Epetra_Map> dofrowmap, Teuchos::RCP<Epetra_Map> condmap);

      //@}

      /** \name Maps */
      //@{

      /// get the condition dof map
      Teuchos::RCP<Epetra_Map> CondDofMap() const { return conddofmap_; }

      /// get the other dof map
      Teuchos::RCP<Epetra_Map> OtherDofMap() const { return otherdofmap_; }

      //@}

      /** \name Extract from full dof vector */
      //@{

      /// extract a condition dof vector from a full row (or column) vector
      Teuchos::RCP<Epetra_Vector> ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full) const;

      /// extract an other dof vector from a full row (or column) vector
      Teuchos::RCP<Epetra_Vector> ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full) const;

      /// extract a condition dof vector from a full row (or column) vector
      void ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> cond) const;

      /// extract an other dof vector from a full row (or column) vector
      void ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> other) const;

      //@}

      /** \name Insert from full dof vector */
      //@{

      /// Put a condition dof vector into a full row vector
      Teuchos::RCP<Epetra_Vector> InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond) const;

      /// Put an other dof vector into a full row vector
      Teuchos::RCP<Epetra_Vector> InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other) const;

      /// Put a condition dof vector into a full row vector
      void InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const;

      /// Put an other dof vector into a full row vector
      void InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const;

      //@}

    private:

      void SetupMaps(const Teuchos::RCP<DRT::Discretization> dis,
                     const std::set<int>& conddofset, const std::set<int>& otherdofset);

      /// the original map
      Teuchos::RCP<const Epetra_Map> dofrowmap_;

      /// condition dof map
      Teuchos::RCP<Epetra_Map> conddofmap_;

      /// other dof map
      Teuchos::RCP<Epetra_Map> otherdofmap_;

      /// communication between condition dof map and full row dof map
      Teuchos::RCP<Epetra_Import> condimporter_;

      /// communication between other dof map and full row dof map
      Teuchos::RCP<Epetra_Import> otherimporter_;
    };


    /// ExtractorCond that selects all dofs of the nodes from a DRT::Condition
    /*!
      This condition is meant to select all FSI interface dofs, as the
      interface nodes are marked by a DRT::Condition.
     */
    class ExtractorCondInCondition
    {
    public:
      ExtractorCondInCondition(const Teuchos::RCP<DRT::Discretization> dis,
                               std::string condname)
      {
        // we use a RCP here so we get cheap copies
        conds_ = Teuchos::rcp(new std::vector<DRT::Condition*>());
        dis->GetCondition(condname, *conds_);
      }

      /// required function call operator
      bool operator()(const Teuchos::RCP<DRT::Discretization> dis,
                      const DRT::Node* node,
                      const std::vector<int>& dof,
                      int j) const
      {
        for (unsigned i=0; i<conds_->size(); ++i)
        {
          const vector<int>* n = (*conds_)[i]->Nodes();

          // DRT::Condition nodes are ordered by design! So we can perform a
          // binary search here.
          if (std::binary_search(n->begin(), n->end(), node->Id()))
            return true;
        }
        return false;
      }

    private:
      Teuchos::RCP<std::vector<DRT::Condition*> > conds_;
    };


    /// ExtractorCond that selects all nodal dofs up to a given position
    /*!
      This condition is meant to select all dim velocity dofs, for example,
      and leave the pressure dofs to the other map.
     */
    class ExtractorCondMaxPos
    {
    public:
      ExtractorCondMaxPos(int maxpos) : maxpos_(maxpos) {}

      /// required function call operator
      bool operator()(const Teuchos::RCP<DRT::Discretization> dis,
                      const DRT::Node* node,
                      const std::vector<int>& dof,
                      int j) const
      {
        return j < maxpos_;
      }

    private:
      int maxpos_;
    };


    /// a "and" combination of two ExtractorCond classes
    /*!
      Such a combination is required, for example, to select all velocity dofs
      at the FSI interface.
     */
    template<class C1,class C2>
    class ExtractorCondAnd
    {
    public:
      ExtractorCondAnd(const C1& c1, const C2& c2) : c1_(c1), c2_(c2) {}

      /// required function call operator
      bool operator()(const Teuchos::RCP<DRT::Discretization> dis,
                      const DRT::Node* node,
                      const std::vector<int>& dof,
                      int j) const
      {
        return c1_(dis,node,dof,j) and c2_(dis,node,dof,j);
      }

    private:
      C1 c1_;
      C2 c2_;
    };

    /// helper function to create a ExtractorCondAnd
    template<class C1,class C2>
    ExtractorCondAnd<C1,C2> CondAnd(const C1& c1, const C2& c2)
    {
      return ExtractorCondAnd<C1,C2>(c1,c2);
    }


    /// a "and" combination of two ExtractorCond classes
    template<class C1,class C2>
    class ExtractorCondOr
    {
    public:
      ExtractorCondOr(C1& c1, C2& c2) : c1_(c1), c2_(c2) {}

      /// required function call operator
      bool operator()(const Teuchos::RCP<DRT::Discretization> dis,
                      const DRT::Node* node,
                      const std::vector<int>& dof,
                      int j) const
      {
        return c1_(dis,node,dof,j) or c2_(dis,node,dof,j);
      }

    private:
      C1 c1_;
      C2 c2_;
    };

    /// helper function to create a ExtractorCondOr
    template<class C1,class C2>
    ExtractorCondOr<C1,C2> CondOr(const C1& c1, const C2& c2)
    {
      return ExtractorCondOr<C1,C2>(c1,c2);
    }

  }
}


template<class C>
void DRT::UTILS::MapExtractor::Setup(const Teuchos::RCP<DRT::Discretization> dis,
                                     const C& condition)
{
  std::set<int> conddofset;
  std::set<int> otherdofset;

  int numrownodes = dis->NumMyRowNodes();
  for (int i=0; i<numrownodes; ++i)
  {
    DRT::Node* node = dis->lRowNode(i);
    std::vector<int> dof = dis->Dof(node);
    for (unsigned j=0; j<dof.size(); ++j)
    {
      if (condition(dis,node,dof,j))
      {
        conddofset.insert(dof[j]);
      }
      else
      {
        otherdofset.insert(dof[j]);
      }
    }
  }

  SetupMaps(dis, conddofset, otherdofset);
}

#endif
#endif
