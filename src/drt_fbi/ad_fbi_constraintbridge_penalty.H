/*----------------------------------------------------------------------*/
/*! \file

\brief Abstract class to be overloaded by different adapter implementations connecting the penalty
constraint enforcement technique with a discretization approach for Fluid-beam interaction.

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/

#ifndef AD_FBI_CONSTRAINTBRIDGE_PENALTY_H_
#define AD_FBI_CONSTRAINTBRIDGE_PENALTY_H_

#include <Teuchos_RCP.hpp>
#include "../linalg/linalg_utils.H"
#include <Epetra_FEVector.h>
#include "ad_fbi_constraintbridge.H"

namespace ADAPTER
{
  class ConstraintEnforcerFactory;
  class FBIConstraintBridgePenalty : public FBIConstraintBridge  // todo fix this
  {
    friend ConstraintEnforcerFactory;

   public:
    /// empty destructor
    virtual ~FBIConstraintBridgePenalty(){};

    /**
     * \brief Initializes all members of the class     *
     */
    virtual void Setup(const Epetra_Map* beam_map, const Epetra_Map* fluid_map);

    /**
     * \brief Computes the coupling matrices
     *
     * This is where the magic happens. The global meshtying contributions are integrated using
     * information of the beam elements, the fluid elements and their position relative to each
     * other. Depending on
     *
     */

    virtual void Evaluate(
        Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>> discretizations)
        override;

    /// Clears the pair vector and builds the pairs anew
    virtual void ResetBridge() override;

    /// Matrix containing only slave side contributions \f$C_ss\f$
    Teuchos::RCP<LINALG::SparseMatrix const> const GetCss() { return Css_; };

    /// Matrix containing only master side contributions \f$C_mm\f$
    Teuchos::RCP<LINALG::SparseMatrix const> const GetCmm() { return Cmm_; };

    /// Matrix containing mixed slave side contributions \f$C_ms\f$
    Teuchos::RCP<LINALG::SparseMatrix const> const GetCms() { return Cms_; };

    /// Matrix containing mixed master side contributions \f$C_sm\f$
    Teuchos::RCP<LINALG::SparseMatrix const> const GetCsm() { return Csm_; };

    /// Force vector acting on the master side \f$f_m\f$
    Teuchos::RCP<Epetra_FEVector const> const GetFm() { return fm_; };

    /// Force vector acting on the slave side \f$f_s\f$
    Teuchos::RCP<Epetra_FEVector const> const GetFs() { return fs_; };

   protected:
    /// Standard constructor
    FBIConstraintBridgePenalty()
        : Css_(Teuchos::null),
          Cmm_(Teuchos::null),
          Cms_(Teuchos::null),
          Csm_(Teuchos::null),
          fm_(Teuchos::null),
          fs_(Teuchos::null){};

   private:
    /// Coupling matrix containing only slave side contributions \f$C_ss\f$
    Teuchos::RCP<LINALG::SparseMatrix> Css_;

    /// Coupling matrix containing only master side contributions \f$C_mm\f$
    Teuchos::RCP<LINALG::SparseMatrix> Cmm_;

    /// Coupling matrix containing mixed slave side contributions \f$C_ms\f$
    Teuchos::RCP<LINALG::SparseMatrix> Cms_;

    /// Coupling matrix containing mixed master side contributions \f$C_sm\f$
    Teuchos::RCP<LINALG::SparseMatrix> Csm_;

    /// Force vector acting on the master side \f$f_m\f$
    Teuchos::RCP<Epetra_FEVector> fm_;

    /// Force vector acting on the slave side \f$f_s\f$
    Teuchos::RCP<Epetra_FEVector> fs_;
  };
}  // namespace ADAPTER

#endif /* AD_FBI_CONSTRAINTBRIDGE_PENALTY_H_ */
