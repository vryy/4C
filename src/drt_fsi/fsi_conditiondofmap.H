
#ifdef CCADISCRET

#ifndef FSI_CONDITIONDOFMAP_H
#define FSI_CONDITIONDOFMAP_H

#include <map>
#include <set>
#include <string>
#include <vector>

#include <Teuchos_RCP.hpp>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>


namespace DRT
{
class Discretization;
}


namespace FSI
{

  /// Split a dof row map in two and establish the communication pattern between those maps
  /*!

    In baci we use conditions to select a set of nodes. Oftentimes we are
    interested to work with a dof vector that contains the dofs of all nodes
    in a conditions. Sometimes we also have to deal with all dofs that are not
    selected by the condition. And we need to transfer back and forth between
    the full dof row map and one of its parts. This is exactly that this class
    is all about.

    \note We are in an FSI context here, so we actually extract only the
    velocity dofs at the interface. This makes this class special. We could
    generalize this behavior.

    \note We work on row maps. The maps we deal with are meant to be nonoverlapping.

    \author u.kue
    \date 12/07
  */
  class ConditionDofMap
  {
  public:

    /// A condition dof map is always tied to a discretization
    explicit ConditionDofMap(Teuchos::RCP<DRT::Discretization> dis);

    /// Extract the selected nodes, the nodal dofs and create the condition dof map
    /*!
      This has to be called just once during setup.
     */
    void SetupCondDofMap(std::string condname);

    /// Setup condition dof map using a given map
    /*!
      If there is no explicit condition but we already have a map that covers
      all dofs we are interessted in, it is still fine to use this class.
     */
    void SetupCondDofMap(Teuchos::RCP<Epetra_Map> condmap);

    /// Create the other dof map, the map of non-selected dofs
    /*!
      The setup of the other dof map is separated from the condition dof map
      as it is not needed in all cases. But if you need the other dof map, you
      have to call this method during setup.
     */
    void SetupOtherDofMap();

    /// Cleanup any setup data that is not needed during normal operation
    void FillComplete() { nodes_.clear(); }

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

    /// the discretization we are working on
    const DRT::Discretization& Discret() const { return *dis_; }

    /// the nodes covered by this condition
    /*!
      This is setup data. We make it available here because others might want
      to use it during their setup.
     */
    const std::vector<int>& Nodes() const { return nodes_; }

  private:

    /// find all local nodes from discretization marked with condition
    void FindCondNodes(const DRT::Discretization& dis, std::string condname, std::set<int>& nodes);

    /// our discretization
    Teuchos::RCP<DRT::Discretization> dis_;

    /// nodes selected by our condition
    std::vector<int> nodes_;

    /// condition dof map
    Teuchos::RCP<Epetra_Map> conddofmap_;

    /// other dof map
    Teuchos::RCP<Epetra_Map> otherdofmap_;

    /// communication between condition dof map and full row dof map
    Teuchos::RCP<Epetra_Import> condimporter_;

    /// communication between other dof map and full row dof map
    Teuchos::RCP<Epetra_Import> otherimporter_;
  };
}

#endif
#endif
