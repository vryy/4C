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

    \note This class is specific to FSI coupling with structure being the
    master and fluid the slave.

    \note The force coupling in SlaveToMaster() expects interface loads on the
    fluid side and returns consistent interface forces on the structural side.

    \author tk (Georg Hammerl)
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

    /// transfer a dof vector from master to slave
    /*!
      Expects the structural interface displacements and returns interface
      displacements on the fluid side.
      Compute sv = D^{-1}M(mv)
     */
    RefCountPtr<Epetra_Vector> MasterToSlave
    (
        RefCountPtr<Epetra_Vector> mv ///< structural interface displacements
    );

    /// transfer a dof vector from slave to master
    /*!
      Expects the fluid interface loads and returns consistent interface
      forces at the structural side. So here we have a subtle change of the
      meaning these numbers.
      Computes mv = M^{T}(sv).
     */
    RefCountPtr<Epetra_Vector> SlaveToMaster
    (
        RefCountPtr<Epetra_Vector> sv   ///< fluid interface loads
    );

    /// the interface dof map of the master side
    RefCountPtr<Epetra_Map> MasterDofMap() const { return masterdofmap_; }

    /// the interface dof map of the slave side
    RefCountPtr<Epetra_Map> SlaveDofMap()  const { return slavedofmap_; }
    
    /// the interface dof map of the slave side
    RefCountPtr<Epetra_Map> SlaveColMap()  const { return slavecolmap_; }
    
    RCP<MORTAR::MortarInterface> Interface() const { return interface_;}

    

  private:
	
    RCP<MORTAR::MortarInterface>  interface_;    ///< interface
    RefCountPtr<Epetra_Map>       masterdofmap_; ///< map of master dofs
    RefCountPtr<Epetra_Map>       slavedofmap_;  ///< map of slave dofs without pressure
    RefCountPtr<Epetra_Map> slavecolmap_;  ///< map of slave dofs without pressure, ghosted nodes included

    RefCountPtr<Epetra_CrsMatrix> D_;            ///< mortar matrix
    RefCountPtr<Epetra_CrsMatrix> Dinv_;         ///< inverse of mortar matrix #D_
    RefCountPtr<Epetra_CrsMatrix> M_;            ///< mortar matrix
  };
}

#endif // FSI_COUPLING_MORTAR_CONTACT_H
#endif // CCADISCRET
