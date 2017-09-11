/*----------------------------------------------------------------------*/
/*!
\file adapter_coupling_mortar.H

\brief A class providing coupling capabilities based on mortar methods

\level 2

<pre>
\maintainer Rui Fang
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef ADAPTER_COUPLING_MORTAR_H
#define ADAPTER_COUPLING_MORTAR_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include "../drt_lib/drt_dserror.H"
#include "adapter_coupling.H"

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
  class CouplingMortar : public CouplingBase
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
    void SetupInterface(
        const Teuchos::RCP<DRT::Discretization>&           masterdis,            ///< master discretization
        const Teuchos::RCP<DRT::Discretization>&           slavedis,             ///< slave discretization
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

    /// create integration cells
    void EvaluateGeometry(
        std::vector<Teuchos::RCP<MORTAR::IntCell> >&   intcells   //!< vector of mortar integration cells
        );

    /// filling matrices by using mortar interface using reference configuration
    void Evaluate();

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

    /// the mortar interface itself
    Teuchos::RCP<MORTAR::MortarInterface> Interface() const { return interface_;}

    virtual Teuchos::RCP<LINALG::SparseMatrix> GetDMatrix() const
    {
      if(D_==Teuchos::null)
        dserror("ERROR: D Matrix is null pointer!");
      return D_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> GetDinvMatrix() const
    {
      if(Dinv_==Teuchos::null)
        dserror("ERROR: DInv Matrix is null pointer!");
      return Dinv_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMMatrix() const
    {
      if(M_==Teuchos::null)
        dserror("ERROR: M Matrix is null pointer!");
      return M_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> GetPMatrix() const
    {
      if(P_==Teuchos::null)
        dserror("ERROR: P Matrix is null pointer!");
      return P_;
    };

    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMortarTrafo() const
    {
      if(P_==Teuchos::null)
        dserror("ERROR: P Matrix is null pointer!");
      return P_;
    };
    /// @name Conversion between master and slave
    //@{

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_Vector> MasterToSlave
    (
        Teuchos::RCP<Epetra_Vector> mv ///< e.g. for FSI: structural interface displacements
    ) const;

    /// transfer a dof vector from master to slave
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave
    (
        Teuchos::RCP<Epetra_MultiVector> mv ///< e.g. for FSI: structural interface displacements
    ) const;

    /// transfer a dof vector from master to slave (const version)
    Teuchos::RCP<Epetra_Vector> MasterToSlave
    (
        Teuchos::RCP<const Epetra_Vector> mv ///< e.g. for FSI: structural interface displacements
    ) const;

    /// transfer a dof vector from master to slave (const version)
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave
    (
        Teuchos::RCP<const Epetra_MultiVector> mv ///< e.g. for FSI: structural interface displacements
    ) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_Vector> SlaveToMaster
    (
        Teuchos::RCP<Epetra_Vector> sv ///< e.g. for FSI: fluid interface loads
    ) const;

    /// transfer a dof vector from slave to master
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster
    (
        Teuchos::RCP<Epetra_MultiVector> sv ///< e.g. for FSI: fluid interface loads
    ) const;

    /// transfer a dof vector from slave to master (const version)
    Teuchos::RCP<Epetra_Vector> SlaveToMaster
    (
        Teuchos::RCP<const Epetra_Vector> sv ///< e.g. for FSI: fluid interface loads
    ) const;

    /// transfer a dof vector from slave to master (const version)
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster
    (
        Teuchos::RCP<const Epetra_MultiVector> sv ///< e.g. for FSI: fluid interface loads
    ) const;

    /// transfer a dof vector from master to slave
    void MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv, ///< master vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> sv ///< slave vector (containing result)
        ) const;

    /// transfer a dof vector from slave to master
    void SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv, ///< slave vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> mv ///< master vector (containing result)
        ) const;

    //@}

    /** \name Coupled maps */
    //@{

    /// the interface dof row map of the master side
    Teuchos::RCP<const Epetra_Map> MasterDofMap() const { return pmasterdofrowmap_; }

    /// the interface dof row map of the slave side
    Teuchos::RCP<const Epetra_Map> SlaveDofMap()  const { return pslavedofrowmap_; }

    //@}

    /** \name Condensation methods */
    //@{

    /// do condensation of Lagrange multiplier and slave-sided dofs
    void MortarCondensation(
        Teuchos::RCP<LINALG::SparseMatrix>& k,    ///< in:  tangent matrix w/o condensation
                                                  ///< out: tangent matrix w/  condensation
        Teuchos::RCP<Epetra_Vector>& rhs          ///< in:  rhs vector     w/o condensation
                                                  ///< out: rhs vector     w/  condensation
        ) const ;

    /// recover slave-sided dofs
    void MortarRecover(
        Teuchos::RCP<LINALG::SparseMatrix>& k,    ///< in: tangent matrix
        Teuchos::RCP<Epetra_Vector>& inc          ///< in:  solution vector     w/o condensation
                                                  ///< out: solution vector     w/  condensation
        ) const ;

    //@}

  protected:

    // create projection operator Dinv*M
    virtual void CreateP();

    // check is slave dofs have dirichlet constraints
    // this is not allows --> over constrained
    void CheckSlaveDirichletOverlap(
        const Teuchos::RCP<DRT::Discretization>& slavedis,
        const Epetra_Comm& comm);

    // back transformation to initial parallel distribution
    void MatrixRowColTransform();

    /// check setup call
    const bool& IsSetup() const { return issetup_; };

    /// check init and setup call
    virtual void CheckSetup() const
    {
      if ( !IsSetup())
        dserror("ERROR: Call Setup() first!");
    }

  private:

    /// perform mesh relocation
    void MeshRelocation(
        Teuchos::RCP<DRT::Discretization> slavedis,
        Teuchos::RCP<DRT::Discretization> aledis,
        Teuchos::RCP<const Epetra_Map>    masterdofrowmap,
        Teuchos::RCP<const Epetra_Map>    slavedofrowmap,
        Teuchos::RCP<Epetra_Vector>&      idisp,
        const Epetra_Comm&                comm,
        bool                              slavewithale
        );

  protected:
    bool                                  issetup_;         ///< check for setup

    Teuchos::RCP<MORTAR::MortarInterface> interface_;       ///< interface
    Teuchos::RCP<const Epetra_Map>        masterdofrowmap_; ///< map of master row dofs (after parallel redist.)
    Teuchos::RCP<const Epetra_Map>        slavedofrowmap_;  ///< map of slave row dofs  (after parallel redist.)
    Teuchos::RCP<const Epetra_Map>        pmasterdofrowmap_;///< map of master row dofs (before parallel redist.)
    Teuchos::RCP<const Epetra_Map>        pslavedofrowmap_; ///< map of slave row dofs  (before parallel redist.)

    Teuchos::RCP<LINALG::SparseMatrix>    D_;               ///< mortar matrix
    Teuchos::RCP<LINALG::SparseMatrix>    Dinv_;            ///< inverse of mortar matrix D_
    Teuchos::RCP<LINALG::SparseMatrix>    M_;               ///< mortar matrix
    Teuchos::RCP<LINALG::SparseMatrix>    P_;               ///< product of Dinv_ \times M_
  };
}

#endif /* ADAPTER_COUPLING_MORTAR_H */
