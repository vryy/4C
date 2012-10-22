/*----------------------------------------------------------------------*/
/*!
\file adapter_coupling_mortar.H

\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef FSI_COUPLING_MORTAR_H
#define FSI_COUPLING_MORTAR_H

#include <Teuchos_RefCountPtr.hpp>
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
                   RCP<Epetra_Map> masterdofrowmap,
                   RCP<Epetra_Map> slavedofrowmap,
                   const Epetra_Comm& comm,
                   bool structslave);

    ///filling matrices (standard case)
    void Evaluate(RCP<Epetra_Vector> idisp);

    ///filling matrices (case of transferring same dofs on two different meshes)
    void Evaluate(RCP<Epetra_Vector> idispma, RCP<Epetra_Vector> idispsl);

    RefCountPtr<Epetra_Vector> MasterToSlave
    (
        RefCountPtr<Epetra_Vector> mv ///< structural interface displacements
    ) const;

    /// const version
    RefCountPtr<Epetra_Vector> MasterToSlave
    (
        RefCountPtr<const Epetra_Vector> mv ///< structural interface displacements
    ) const;

    /// transfer a dof vector from slave to master
    RefCountPtr<Epetra_Vector> SlaveToMaster
    (
        RefCountPtr<Epetra_Vector> sv   ///< fluid interface loads
    ) const;

    /// const version
    RefCountPtr<Epetra_Vector> SlaveToMaster
    (
        RefCountPtr<const Epetra_Vector> sv   ///< fluid interface loads
    ) const;

    /// the mortar interface
    RCP<MORTAR::MortarInterface> Interface() const { return interface_;}

    /// the interface dof map of the master side
    RefCountPtr<Epetra_Map> MasterDofRowMap() const { return masterdofrowmap_; }

    /// the interface dof map of the slave side
    RefCountPtr<Epetra_Map> SlaveDofRowMap()  const { return slavedofrowmap_; }

    /// the interface dof map of the lagrainge multiplier
    RefCountPtr<Epetra_Map> LmDofRowMap()  const { return glmdofrowmap_; }

    /// the mortar projection operator
    RCP<LINALG::SparseMatrix> GetMortarTrafo() const { return DinvM_;}

    /// return constraint matrix
    RCP<LINALG::SparseMatrix> GetConMatrix() const { return conmatrix_;}

    /// return mortar M matrix
    RCP<LINALG::SparseMatrix> GetMMatrix() const { return M_;}

    /// return mortar D matrix
    RCP<LINALG::SparseMatrix> GetDMatrix() const { return D_;}

    /// return inverted mortar matrix \f$D^{-1}\f$
    RCP<LINALG::SparseMatrix> GetDinvMatrix() const { return Dinv_;}


  private:

    ///filling matrices by using mortar interface
    void Evaluate();

    RCP<MORTAR::MortarInterface>  interface_;        ///< interface
    RefCountPtr<Epetra_Map>       masterdofrowmap_;  ///< map of master row dofs
    RefCountPtr<Epetra_Map>       slavedofrowmap_;   ///< map of slave row dofs
    RefCountPtr<Epetra_Map>       glmdofrowmap_;   ///< map of slave row dofs
    RCP<LINALG::SparseMatrix>     conmatrix_ ;

    RefCountPtr<LINALG::SparseMatrix> D_;            ///< mortar matrix
    RefCountPtr<LINALG::SparseMatrix> Dinv_;         ///< inverse of mortar matrix #D_
    RefCountPtr<LINALG::SparseMatrix> M_;            ///< mortar matrix
    RefCountPtr<LINALG::SparseMatrix> DinvM_;        ///< product of Dinv_ \times M_
  };
}

#endif // FSI_COUPLING_MORTAR_CONTACT_H
