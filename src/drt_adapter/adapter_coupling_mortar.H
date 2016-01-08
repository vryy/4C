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

#ifndef ADAPTER_COUPLING_MORTAR_H
#define ADAPTER_COUPLING_MORTAR_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include "../drt_lib/drt_dserror.H"

// forward declarations
namespace LINALG
{
  class SparseMatrix;
}

namespace DRT
{
  class Discretization;
  class Element;
  class Node;
}

namespace MORTAR
{
  class MortarInterface;
  class IntCell;
}

namespace ADAPTER
{
  /// Couple non-matching interface meshes using mortar method
  /*!
    This is a generic class used to couple any non-matching meshes
    (or more general: discretizations) at interfaces. The current
    applications in BACI encompass FSI coupling algorithms (i.e. to
    interpolate between fluid and structure fields at the interface)
    and fluid mesh tying algorithms (i.e. to couple non-matching
    Eulerian fluid meshes). All the hard work is actually done by
    the MORTAR::MortarInterface class (thus we use the mortar method).

    The major part of this code is the Setup() method that gets the
    non-matching interface meshes on input, initializes the mortar
    interface and computes the so-called coupling matrices D and M.

    The actual coupling methods MasterToSlave() and SlaveToMaster()
    just evaluate one simple equation each, i.e. primal variables
    are projected from master to slave side via \f$D^{-1} M\f$ when
    calling MasterToSlaver(), and dual variables are projected from
    slave to master side via \f$M^T D^{-T}\f$ when calling SlaveToMaster().

    Whenever you want to add a new problem class, check whether you
    can re-use one of the already existing Setup() methods. If not,
    feel free to write your own tailored Setup() method.
   */
  class CouplingMortar
  {
  public:

    /// empty constructor
    CouplingMortar();

    //! setup the machinery (generalized version)
    /// master and slave discretizations are identical in case of sliding ALE or fluid/scatra meshtying
    /// ALE discretization is Teuchos::null in case of sliding ALE or fluid/scatra meshtying
    void Setup(
        const Teuchos::RCP<DRT::Discretization>&   masterdis,            ///< master discretization
        const Teuchos::RCP<DRT::Discretization>&   slavedis,             ///< slave discretization
        const Teuchos::RCP<DRT::Discretization>&   aledis,               ///< ALE discretization
        const std::vector<int>&                    coupleddof,           ///< vector defining coupled degrees of freedom
        const std::string&                         couplingcond,         ///< string for coupling condition
        const Epetra_Comm&                         comm,                 ///< communicator
        const bool                                 slavewithale=false,   ///< flag defining if slave is ALE
        const bool                                 slidingale=false,     ///< flag indicating sliding ALE case
        const int                                  nds_master=0,         ///< master dofset number
        const int                                  nds_slave=0           ///< slave dofset number
        );

    //! setup the machinery (generalized version)
    /// master and slave discretizations are identical in case of sliding ALE or fluid/scatra meshtying
    /// ALE discretization is Teuchos::null in case of sliding ALE or fluid/scatra meshtying
    void Setup(
        const Teuchos::RCP<DRT::Discretization>&           masterdis,            ///< master discretization
        const Teuchos::RCP<DRT::Discretization>&           slavedis,             ///< slave discretization
        const Teuchos::RCP<DRT::Discretization>&           aledis,               ///< ALE discretization
        const std::vector<int>&                            coupleddof,           ///< vector defining coupled degrees of freedom
        const std::map<int,DRT::Node*>&                    mastergnodes,         ///< master nodes, including ghosted nodes
        const std::map<int,DRT::Node*>&                    slavegnodes,          ///< slave nodes, including ghosted nodes
        const std::map<int,Teuchos::RCP<DRT::Element> >&   masterelements,       ///< master elements
        const std::map<int,Teuchos::RCP<DRT::Element> >&   slaveelements,        ///< slave elements
        const Epetra_Comm&                                 comm,                 ///< communicator
        const bool                                         slavewithale=false,   ///< flag defining if slave is ALE
        const bool                                         slidingale=false,     ///< flag indicating sliding ALE case
        const int                                          nds_master=0,         ///< master dofset number
        const int                                          nds_slave=0           ///< slave dofset number
        );

    //! setup the machinery (hack version (used in UQ to obtain surface normals of slave interface ONLY))
    void SetupForUQAbuseNormalCalculation(
        Teuchos::RCP<DRT::Discretization>   slavedis,   ///< slave discretization
        const Epetra_Comm&                  comm        ///< communicator
        );

    /// perform mesh relocation
    void MeshRelocation(
        Teuchos::RCP<DRT::Discretization> slavedis,
        Teuchos::RCP<DRT::Discretization> aledis,
        Teuchos::RCP<Epetra_Map>          masterdofrowmap,
        Teuchos::RCP<Epetra_Map>          slavedofrowmap,
        Teuchos::RCP<Epetra_Vector>&      idisp,
        const Epetra_Comm&                comm,
        bool                              slavewithale
        );

    /// create integration cells
    std::vector<Teuchos::RCP<MORTAR::IntCell> > EvaluateGeometry();

    /// filling matrices
    void Evaluate(Teuchos::RCP<Epetra_Vector> idisp);

    /// filling matrices (case of transferring same dofs on two different meshes)
    void Evaluate(Teuchos::RCP<Epetra_Vector> idispma, Teuchos::RCP<Epetra_Vector> idispsl);

    //! filling matrices, including mesh correction
    void EvaluateWithMeshRelocation(
        Teuchos::RCP<DRT::Discretization> slavedis,      ///< slave discretization
        Teuchos::RCP<DRT::Discretization> aledis,        ///< ALE discretization
        Teuchos::RCP<Epetra_Vector>&      idisp,         ///< ALE displacements
        const Epetra_Comm&                comm,          ///< communicator
        bool                              slavewithale   ///< flag defining if slave is ALE
        );

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave
    (
        Teuchos::RCP<Epetra_Vector> mv ///< e.g. for FSI: structural interface displacements
    ) const;

    /// transfer a dof vector from master to slave (const version)
    Teuchos::RCP<Epetra_Vector> MasterToSlave
    (
        Teuchos::RCP<const Epetra_Vector> mv ///< e.g. for FSI: structural interface displacements
    ) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster
    (
        Teuchos::RCP<Epetra_Vector> sv ///< e.g. for FSI: fluid interface loads
    ) const;

    /// transfer a dof vector from slave to master (const version)
    Teuchos::RCP<Epetra_Vector> SlaveToMaster
    (
        Teuchos::RCP<const Epetra_Vector> sv ///< e.g. for FSI: fluid interface loads
    ) const;

    /// the mortar interface itself
    Teuchos::RCP<MORTAR::MortarInterface> Interface() const { return interface_;}

    /// the interface dof map of the master side
    Teuchos::RCP<Epetra_Map> MasterDofRowMap() const { return masterdofrowmap_; }

    /// the interface dof map of the slave side
    Teuchos::RCP<Epetra_Map> SlaveDofRowMap()  const { return slavedofrowmap_; }

    /// return mortar D matrix
    Teuchos::RCP<LINALG::SparseMatrix> GetDMatrix() const { return D_;}

    /// return inverted mortar matrix \f$D^{-1}\f$
    Teuchos::RCP<LINALG::SparseMatrix> GetDinvMatrix() const { return Dinv_;}

    /// return mortar M matrix
    Teuchos::RCP<LINALG::SparseMatrix> GetMMatrix() const { return M_;}

    /// return the mortar projection operator \f$P=D^{-1} M\f$
    Teuchos::RCP<LINALG::SparseMatrix> GetMortarTrafo() const { return DinvM_;}

  protected:

    /// check setup call
    const bool& IsSetup() const { return issetup_; };

    /// check init and setup call
    void CheckSetup() const
    {
      if ( !IsSetup())
        dserror("ERROR: Call Setup() first!");
    }

  private:

    /// filling matrices by using mortar interface
    void Evaluate();

    bool                                  issetup_;         ///< check for setup

    Teuchos::RCP<MORTAR::MortarInterface> interface_;       ///< interface
    Teuchos::RCP<Epetra_Map>              masterdofrowmap_; ///< map of master row dofs
    Teuchos::RCP<Epetra_Map>              slavedofrowmap_;  ///< map of slave row dofs

    Teuchos::RCP<LINALG::SparseMatrix>    D_;               ///< mortar matrix
    Teuchos::RCP<LINALG::SparseMatrix>    Dinv_;            ///< inverse of mortar matrix D_
    Teuchos::RCP<LINALG::SparseMatrix>    M_;               ///< mortar matrix
    Teuchos::RCP<LINALG::SparseMatrix>    DinvM_;           ///< product of Dinv_ \times M_
  };
}

#endif /* ADAPTER_COUPLING_MORTAR_H */
