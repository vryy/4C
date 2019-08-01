/*!----------------------------------------------------------------------

\brief adapter for the volmortar framework

\level 2

\maintainer Matthias Mayr

*----------------------------------------------------------------------*/
#ifndef ADAPTER_COUPLING_VOLMORTAR_H_
#define ADAPTER_COUPLING_VOLMORTAR_H_

/*---------------------------------------------------------------------*
 | headers                                                 farah 10/13 |
 *---------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include "../linalg/linalg_fixedsizematrix.H"

#include "adapter_coupling.H"

/*---------------------------------------------------------------------*
 | forward declarations                                    farah 10/13 |
 *---------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT

namespace LINALG
{
  class SparseMatrix;
}

namespace VOLMORTAR
{
  namespace UTILS
  {
    class DefaultMaterialStrategy;
  }
}  // namespace VOLMORTAR

namespace ADAPTER
{
  /// Class for calling volmortar coupling and proper parallel redistr.
  class MortarVolCoupl : public CouplingBase
  {
   public:
    /*!
    \brief Empty constructor

    */
    MortarVolCoupl();

    /*!
    \brief Call parallel redistr. and evaluate volmortar coupl.

    */
    void Init(Teuchos::RCP<DRT::Discretization> dis1, Teuchos::RCP<DRT::Discretization> dis2,
        std::vector<int>* coupleddof12 = NULL, std::vector<int>* coupleddof21 = NULL,
        std::pair<int, int>* dofsets12 = NULL, std::pair<int, int>* dofsets21 = NULL,
        Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy = Teuchos::null,
        bool createauxdofs = true);

    /*!
    \brief Setup this class

    */
    void Setup();

    /*!
    \brief Redistribute discretizations to meet needs of volmortar coupling

    \note Call this method in your global control algorithm inbetween \ref Init()
          and \ref Setup(), in case you need parallel redistribution


    \date   09/16
    \author rauch

    */
    void Redistribute();


    /*!
    \brief Get coupling matrices for field 1 and 2

    */
    Teuchos::RCP<const LINALG::SparseMatrix> GetPMatrix12() const { return P12_; };
    Teuchos::RCP<const LINALG::SparseMatrix> GetPMatrix21() const { return P21_; };

    /*!
    \brief Mortar mapping for 1 to 2 and 2 to 1 - for vectors

    */
    Teuchos::RCP<const Epetra_Vector> ApplyVectorMapping12(
        Teuchos::RCP<const Epetra_Vector> vec) const;
    Teuchos::RCP<const Epetra_Vector> ApplyVectorMapping21(
        Teuchos::RCP<const Epetra_Vector> vec) const;

    /*!
    \brief Mortar mapping for 1 to 2 and 2 to 1 - for matrices

    */
    Teuchos::RCP<LINALG::SparseMatrix> ApplyMatrixMapping12(
        Teuchos::RCP<const LINALG::SparseMatrix> mat) const;
    Teuchos::RCP<LINALG::SparseMatrix> ApplyMatrixMapping21(
        Teuchos::RCP<const LINALG::SparseMatrix> mat) const;

    //@}

    //! assign materials
    void AssignMaterials(Teuchos::RCP<DRT::Discretization> dis1,
        Teuchos::RCP<DRT::Discretization> dis2,
        Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy = Teuchos::null);


    /** \name Conversion between master and slave */
    //@{
    /// There are different versions to satisfy all needs. The basic
    /// idea is the same for all of them.

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave(Teuchos::RCP<Epetra_Vector> mv) const
    {
      return MasterToSlave(Teuchos::rcp_static_cast<const Epetra_Vector>(mv));
    }

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(Teuchos::RCP<Epetra_Vector> sv) const
    {
      return SlaveToMaster(Teuchos::rcp_static_cast<const Epetra_Vector>(sv));
    }

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave(Teuchos::RCP<Epetra_MultiVector> mv) const
    {
      return MasterToSlave(Teuchos::rcp_static_cast<const Epetra_MultiVector>(mv));
    }

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(Teuchos::RCP<Epetra_MultiVector> sv) const
    {
      return SlaveToMaster(Teuchos::rcp_static_cast<const Epetra_MultiVector>(sv));
    }

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_Vector> MasterToSlave(Teuchos::RCP<const Epetra_Vector> mv) const;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_Vector> SlaveToMaster(Teuchos::RCP<const Epetra_Vector> sv) const;

    /// transfer a dof vector from master to slave
    virtual Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv) const;

    /// transfer a dof vector from slave to master
    virtual Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv) const;

    /// transfer a dof vector from master to slave
    virtual void MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv, Teuchos::RCP<Epetra_MultiVector> sv) const;

    /// transfer a dof vector from slave to master
    virtual void SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv, Teuchos::RCP<Epetra_MultiVector> mv) const;

    //@}

    /** \name Coupled maps */
    //@{

    /// the interface dof map of the master side
    virtual Teuchos::RCP<const Epetra_Map> MasterDofMap() const;

    /// the interface dof map of the slave side
    virtual Teuchos::RCP<const Epetra_Map> SlaveDofMap() const;

    //@}

   private:
    /*!
    \brief Create auxiliary dofsets for multiphysics if necessary

    */
    void CreateAuxDofsets(Teuchos::RCP<DRT::Discretization> dis1,
        Teuchos::RCP<DRT::Discretization> dis2, std::vector<int>* coupleddof12,
        std::vector<int>* coupleddof21);

    /// check setup call
    const bool& IsSetup() const { return issetup_; };

    /// check setup call
    const bool& IsInit() const { return isinit_; };

    /// check init and setup call
    void CheckSetup() const
    {
      if (not IsSetup()) dserror("ERROR: Call Setup() first!");
    }

    /// check init and setup call
    void CheckInit() const
    {
      if (not IsInit()) dserror("ERROR: Call Init() first!");
    }

   private:
    bool issetup_;  ///< check for setup
    bool isinit_;   ///< check for init

    // mortar matrices and projector
    // s1 = P12 * s2
    // s2 = P21 * s1
    Teuchos::RCP<LINALG::SparseMatrix>
        P12_;  ///< global Mortar projection matrix P Omega_2 -> Omega_1
    Teuchos::RCP<LINALG::SparseMatrix>
        P21_;  ///< global Mortar projection matrix P Omega_1 -> Omega_2

    Teuchos::RCP<DRT::Discretization> masterdis_;
    Teuchos::RCP<DRT::Discretization> slavedis_;

    std::vector<int>* coupleddof12_;
    std::vector<int>* coupleddof21_;
    std::pair<int, int>* dofsets12_;
    std::pair<int, int>* dofsets21_;
    Teuchos::RCP<VOLMORTAR::UTILS::DefaultMaterialStrategy> materialstrategy_;
  };
}  // namespace ADAPTER

#endif /* ADAPTER_COUPLING_VOLMORTAR_H_ */
