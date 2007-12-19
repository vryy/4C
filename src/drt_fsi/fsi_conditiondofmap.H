
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

  ///
  class ConditionDofMap
  {
  public:

    explicit ConditionDofMap(Teuchos::RCP<DRT::Discretization> dis);

    ///
    void SetupCondDofMap(std::string condname);

    ///
    void SetupCondDofMap(Teuchos::RCP<Epetra_Map> condmap);

    ///
    void SetupOtherDofMap();

    ///
    void FillComplete() { nodes_.clear(); }

    /** \name Maps */
    //@{

    ///
    Teuchos::RCP<Epetra_Map> CondDofMap() const { return conddofmap_; }

    ///
    Teuchos::RCP<Epetra_Map> OtherDofMap() const { return otherdofmap_; }

    //@}

    /** \name Extract from full dof vector */
    //@{

    ///
    Teuchos::RCP<Epetra_Vector> ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full) const;

    ///
    Teuchos::RCP<Epetra_Vector> ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full) const;

    ///
    void ExtractCondVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> cond) const;

    ///
    void ExtractOtherVector(Teuchos::RCP<const Epetra_Vector> full, Teuchos::RCP<Epetra_Vector> other) const;

    //@}

    /** \name Insert from full dof vector */
    //@{

    ///
    Teuchos::RCP<Epetra_Vector> InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond) const;

    ///
    Teuchos::RCP<Epetra_Vector> InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other) const;

    ///
    void InsertCondVector(Teuchos::RCP<const Epetra_Vector> cond, Teuchos::RCP<Epetra_Vector> full) const;

    ///
    void InsertOtherVector(Teuchos::RCP<const Epetra_Vector> other, Teuchos::RCP<Epetra_Vector> full) const;

    //@}

    ///
    const DRT::Discretization& Discret() const { return *dis_; }

    ///
    const std::vector<int>& Nodes() const { return nodes_; }

  private:

    /// find all local nodes from discretization marked with condition
    void FindCondNodes(const DRT::Discretization& dis, std::string condname, std::set<int>& nodes);

    ///
    Teuchos::RCP<DRT::Discretization> dis_;

    ///
    std::vector<int> nodes_;

    ///
    //std::map<int, std::vector<int> > dofs_;

    ///
    Teuchos::RCP<Epetra_Map> conddofmap_;

    ///
    Teuchos::RCP<Epetra_Map> otherdofmap_;

    ///
    Teuchos::RCP<Epetra_Import> condimporter_;

    ///
    Teuchos::RCP<Epetra_Import> otherimporter_;
  };
}

#endif
#endif
