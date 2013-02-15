/*----------------------------------------------------------------------*/
/*!
\file adapter_coupling_mortar.H

\brief A class providing coupling capabilities based on mortar methods

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_COUPLING_MORTAR_H
#define FSI_COUPLING_MORTAR_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include "../drt_lib/drt_discret.H"

// forward declarations

namespace LINALG
{
  class SparseMatrix;
}

namespace MORTAR
{
  class MortarInterface;
}

namespace ADAPTER
{
  /// Couple nonmatching FSI interfaces using mortar method
  /*!
    This is a special class used FSI coupling algorithms to interpolate
    between fluid and structure fields at the interface. The
    hard work is done by the MORTAR::MortarInterface.
    (Thus we use the mortar method here.)

    The major part of this code is the Setup() method that feeds the FSI
    interface and extracts the coupling matrices D and M in
    return.

    The actual coupling methods MasterToSlave() and SlaveToMaster() just
    evaluate one simple equation each.

    \author tk
    \date 06/09
   */
  class CouplingMortar
  {
  public:

    /// empty constructor
    CouplingMortar();

    /// setup the machinery
    void Setup( DRT::Discretization& masterdis,
                DRT::Discretization& slavedis,
                DRT::Discretization& aledis, ///> we need this for mesh initialization
                const Epetra_Comm& comm,
                bool structslave = false);

    /// setup the machinery
    bool Setup( DRT::Discretization& dis,
                const Epetra_Comm& comm,
                int meshtyingoption,
                bool structslave = false);

    /// setup the machinery
    void Setup( DRT::Discretization& dis );

    /// perform mesh initialization
    void MeshInit( DRT::Discretization& masterdis,
                   DRT::Discretization& slavedis,
                   DRT::Discretization& aledis,
                   Teuchos::RCP<Epetra_Map> masterdofrowmap,
                   Teuchos::RCP<Epetra_Map> slavedofrowmap,
                   const Epetra_Comm& comm,
                   bool structslave);

    ///filling matrices (standard case)
    void Evaluate(Teuchos::RCP<Epetra_Vector> idisp);

    ///filling matrices (case of transferring same dofs on two different meshes)
    void Evaluate(Teuchos::RCP<Epetra_Vector> idispma, Teuchos::RCP<Epetra_Vector> idispsl);

    Teuchos::RCP<Epetra_Vector> MasterToSlave
    (
        Teuchos::RCP<Epetra_Vector> mv ///< structural interface displacements
    ) const;

    /// const version
    Teuchos::RCP<Epetra_Vector> MasterToSlave
    (
        Teuchos::RCP<const Epetra_Vector> mv ///< structural interface displacements
    ) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster
    (
        Teuchos::RCP<Epetra_Vector> sv   ///< fluid interface loads
    ) const;

    /// const version
    Teuchos::RCP<Epetra_Vector> SlaveToMaster
    (
        Teuchos::RCP<const Epetra_Vector> sv   ///< fluid interface loads
    ) const;

    /// the mortar interface
    Teuchos::RCP<MORTAR::MortarInterface> Interface() const { return interface_;}

    /// the interface dof map of the master side
    Teuchos::RCP<Epetra_Map> MasterDofRowMap() const { return masterdofrowmap_; }

    /// the interface dof map of the slave side
    Teuchos::RCP<Epetra_Map> SlaveDofRowMap()  const { return slavedofrowmap_; }

    /// the interface dof map of the lagrainge multiplier
    Teuchos::RCP<Epetra_Map> LmDofRowMap()  const { return glmdofrowmap_; }

    /// the mortar projection operator
    Teuchos::RCP<LINALG::SparseMatrix> GetMortarTrafo() const { return DinvM_;}

    /// return constraint matrix
    Teuchos::RCP<LINALG::SparseMatrix> GetConMatrix() const { return conmatrix_;}

    /// return mortar M matrix
    Teuchos::RCP<LINALG::SparseMatrix> GetMMatrix() const { return M_;}

    /// return mortar D matrix
    Teuchos::RCP<LINALG::SparseMatrix> GetDMatrix() const { return D_;}

    /// return inverted mortar matrix \f$D^{-1}\f$
    Teuchos::RCP<LINALG::SparseMatrix> GetDinvMatrix() const { return Dinv_;}


  private:

    ///filling matrices by using mortar interface
    void Evaluate();

    Teuchos::RCP<MORTAR::MortarInterface> interface_;       ///< interface
    Teuchos::RCP<Epetra_Map>              masterdofrowmap_; ///< map of master row dofs
    Teuchos::RCP<Epetra_Map>              slavedofrowmap_;  ///< map of slave row dofs
    Teuchos::RCP<Epetra_Map>              glmdofrowmap_;    ///< map of slave row dofs
    Teuchos::RCP<LINALG::SparseMatrix>    conmatrix_ ;

    Teuchos::RCP<LINALG::SparseMatrix> D_;            ///< mortar matrix
    Teuchos::RCP<LINALG::SparseMatrix> Dinv_;         ///< inverse of mortar matrix #D_
    Teuchos::RCP<LINALG::SparseMatrix> M_;            ///< mortar matrix
    Teuchos::RCP<LINALG::SparseMatrix> DinvM_;        ///< product of Dinv_ \times M_
  };
}

#endif // FSI_COUPLING_MORTAR_CONTACT_H
