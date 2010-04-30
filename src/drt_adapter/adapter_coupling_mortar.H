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
#ifdef CCADISCRET

#ifndef FSI_COUPLING_MORTAR_H
#define FSI_COUPLING_MORTAR_H

#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Map.h>
#include "../drt_mortar/mortar_interface.H"
#include "../drt_lib/drt_discret.H"

namespace ADAPTER
{
  /// Couple nonmatching FSI interfaces using mortar method
  /*!
    This is a special class used by the Dirichlet-Neumann coupling algorithm
    to interpolate between fluid and structure fields at the interface. The
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
    void Setup( const DRT::Discretization& masterdis,
                const DRT::Discretization& slavedis,
                Epetra_Comm& comm );


    ///filling matrices
    void Evaluate(RCP<Epetra_Vector> idisp);

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

    /// the interface dof map of the master side
    RefCountPtr<Epetra_Map> MasterDofRowMap() const { return masterdofrowmap_; }

    /// the interface dof map of the slave side
    RefCountPtr<Epetra_Map> SlaveDofRowMap()  const { return slavedofrowmap_; }
    
    /// the interface dof map of the slave side
    RefCountPtr<Epetra_Map> SlaveDofColMap()  const { return slavedofcolmap_; }
    
    RCP<MORTAR::MortarInterface> Interface() const { return interface_;}
    
    RCP<LINALG::SparseMatrix> GetMortarTrafo();

    

  private:
	
    RCP<MORTAR::MortarInterface>  interface_;       ///< interface
    RefCountPtr<Epetra_Map>       masterdofrowmap_; ///< map of master row dofs
    RefCountPtr<Epetra_Map>       slavedofrowmap_;  ///< map of slave row dofs without pressure
    RefCountPtr<Epetra_Map>       slavedofcolmap_;  ///< map of slave col dofs without pressure

    RefCountPtr<LINALG::SparseMatrix> D_;            ///< mortar matrix
    RefCountPtr<LINALG::SparseMatrix> Dinv_;         ///< inverse of mortar matrix #D_
    RefCountPtr<LINALG::SparseMatrix> M_;            ///< mortar matrix
    RefCountPtr<LINALG::SparseMatrix> DinvM_;        ///< product of Dinv_ \times M_
  };
}

#endif // FSI_COUPLING_MORTAR_CONTACT_H
#endif // CCADISCRET
