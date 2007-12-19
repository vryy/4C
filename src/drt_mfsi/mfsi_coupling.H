/*----------------------------------------------------------------------*/
/*!
\file mfsi_coupling.H

\brief Enhancement of the FSI coupling for monolithic methods

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef MFSI_COUPLING_H
#define MFSI_COUPLING_H

#include "../drt_fsi/fsi_coupling.H"

namespace MFSI
{

  class CouplingOperator;

  /// monolithic FSI coupling
  /*!

  A FSI coupling class that additionally sets up and stores coupling
  matrices to be used in a monolithic scheme. Here these coupling
  matrices are just (rectangular) identity matrices. But any set of
  mortar matrices could be handled in the same fassion.

  \note We build rectangular matrices here that match the block maps
  exactly. This is important. This way these matrices not only
  transform numbers but move vectors between different maps.

  The condition block equation is always formulated in terms of the
  master interface map. (That is the structural interface map for FSI
  coupling.)

  \author u.kue
  \date 11/07
   */
  class Coupling : public FSI::Coupling
  {
  public:

    /// create coupling matrices for Lagrangian coupling conditions
    void SetupCouplingMatrices(const Epetra_Map& masterdomainmap, const Epetra_Map& slavedomainmap);

    Teuchos::RCP<Epetra_CrsMatrix> MasterToMasterMat() const { return matmm_; }
    Teuchos::RCP<Epetra_CrsMatrix> SlaveToMasterMat() const { return matsm_; }
    Teuchos::RCP<Epetra_CrsMatrix> MasterToMasterMatTrans() const { return matmm_trans_; }
    Teuchos::RCP<Epetra_CrsMatrix> SlaveToMasterMatTrans() const { return matsm_trans_; }

  private:

    /// coupling matrices for Lagrangian multiplier coupling
    Teuchos::RCP<Epetra_CrsMatrix> matmm_;
    Teuchos::RCP<Epetra_CrsMatrix> matsm_;
    Teuchos::RCP<Epetra_CrsMatrix> matmm_trans_;
    Teuchos::RCP<Epetra_CrsMatrix> matsm_trans_;
  };
}

#endif
#endif
