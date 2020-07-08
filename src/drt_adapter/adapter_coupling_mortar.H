/*----------------------------------------------------------------------*/
/*! \file

\brief A class providing coupling capabilities based on mortar methods

\level 2

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
}  // namespace DRT

namespace MORTAR
{
  class MortarInterface;
  class IntCell;
}  // namespace MORTAR

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
    interface and computes the so-called coupling matrices \f$D\f$ and \f$M\f$.

    The actual coupling methods MasterToSlave() and SlaveToMaster()
    just evaluate one simple equation each, i.e. primal variables
    are projected from master to slave side via \f$D^{-1} M\f$ when
    calling MasterToSlave(), and dual variables are projected from
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

    /*! setup the machinery (generalized version)
     *
     *  \note
     *  - Master and slave discretizations are identical in case of sliding ALE or fluid/scatra
     * meshtying
     *  - ALE discretization is Teuchos::null in case of sliding ALE or fluid/scatra meshtying
     */
    void Setup(const Teuchos::RCP<DRT::Discretization>& masterdis,  ///< master discretization
        const Teuchos::RCP<DRT::Discretization>& slavedis,          ///< slave discretization
        const Teuchos::RCP<DRT::Discretization>& aledis,            ///< ALE discretization
        const std::vector<int>& coupleddof,  ///< vector defining coupled degrees of freedom
        const std::string& couplingcond,     ///< string for coupling condition
        const Epetra_Comm& comm,             ///< communicator
        const bool slavewithale = false,     ///< flag defining if slave is ALE
        const bool slidingale = false,       ///< flag indicating sliding ALE case
        const int nds_master = 0,            ///< master dofset number
        const int nds_slave = 0              ///< slave dofset number
    );

    /*! setup the machinery (generalized version)
     *
     *  \note
     *  - Master and slave discretizations are identical in case of sliding ALE or fluid/scatra
     * meshtying
     *  - ALE discretization is Teuchos::null in case of sliding ALE or fluid/scatra meshtying
     */
    void SetupInterface(
        const Teuchos::RCP<DRT::Discretization>& masterdis,  ///< master discretization
        const Teuchos::RCP<DRT::Discretization>& slavedis,   ///< slave discretization
        const std::vector<int>& coupleddof,  ///< vector defining coupled degrees of freedom
        const std::map<int, DRT::Node*>& mastergnodes,  ///< master nodes, including ghosted nodes
        const std::map<int, DRT::Node*>& slavegnodes,   ///< slave nodes, including ghosted nodes
        const std::map<int, Teuchos::RCP<DRT::Element>>& masterelements,  ///< master elements
        const std::map<int, Teuchos::RCP<DRT::Element>>& slaveelements,   ///< slave elements
        const Epetra_Comm& comm,                                          ///< communicator
        const bool slavewithale = false,  ///< flag defining if slave is ALE
        const bool slidingale = false,    ///< flag indicating sliding ALE case
        const int nds_master = 0,         ///< master dofset number
        const int nds_slave = 0           ///< slave dofset number
    );

    //! setup the machinery (hack version (used in UQ to obtain surface normals of slave interface
    //! ONLY))
    virtual void SetupForUQAbuseNormalCalculation(
        Teuchos::RCP<DRT::Discretization> slavedis,  ///< slave discretization
        const Epetra_Comm& comm                      ///< communicator
    );

    /// create integration cells
    virtual void EvaluateGeometry(std::vector<Teuchos::RCP<MORTAR::IntCell>>&
            intcells  //!< vector of mortar integration cells
    );

    /// Compute mortar matrices by using mortar interface using reference configuration
    virtual void Evaluate();

    /// Compute mortar matrices
    virtual void Evaluate(Teuchos::RCP<Epetra_Vector> idisp  ///< [in] ??
    );

    /// Compute mortar matrices (case of transferring same dofs on two different meshes)
    virtual void Evaluate(Teuchos::RCP<Epetra_Vector> idispma,  ///< [in] ??
        Teuchos::RCP<Epetra_Vector> idispsl                     ///< [in] ??
    );

    //! Compute mortar matrices after performing a mesh correction step
    virtual void EvaluateWithMeshRelocation(
        Teuchos::RCP<DRT::Discretization> slavedis,  ///< slave discretization
        Teuchos::RCP<DRT::Discretization> aledis,    ///< ALE discretization
        Teuchos::RCP<Epetra_Vector>& idisp,          ///< ALE displacements
        const Epetra_Comm& comm,                     ///< communicator
        bool slavewithale                            ///< flag defining if slave is ALE
    );

    //! Get the mortar interface itself
    Teuchos::RCP<MORTAR::MortarInterface> Interface() const { return interface_; }

    //! Access to slave side mortar matrix \f$D\f$
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetDMatrix() const
    {
      if (D_ == Teuchos::null) dserror("ERROR: D Matrix is null pointer!");
      return D_;
    };

    //! Access to inverse of slave side mortar matrix \f$D^{-1}\f$
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetDinvMatrix() const
    {
      if (Dinv_ == Teuchos::null) dserror("ERROR: DInv Matrix is null pointer!");
      return Dinv_;
    };

    //! Access to master side mortar matrix \f$M\f$
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMMatrix() const
    {
      if (M_ == Teuchos::null) dserror("ERROR: M Matrix is null pointer!");
      return M_;
    };

    //! Access to mortar projection operator \f$P\f$
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetPMatrix() const
    {
      if (P_ == Teuchos::null) dserror("ERROR: P Matrix is null pointer!");
      return P_;
    };

    /*! \brief Access to mortar projection operator \f$P\f$
     *
     *  \warning Redundant with GetPMatrix().
     *
     *  ToDo (mayr.mt) Remove this routine and replace its calls by GetPMatrix().
     */
    virtual Teuchos::RCP<LINALG::SparseMatrix> GetMortarTrafo() const
    {
      if (P_ == Teuchos::null) dserror("ERROR: P Matrix is null pointer!");
      return P_;
    };

    /// @name Conversion between master and slave
    //@{

    /*! \brief Transfer a dof vector from master to slave
     *
     *  \return Slave vector
     */
    Teuchos::RCP<Epetra_Vector> MasterToSlave(
        Teuchos::RCP<Epetra_Vector> mv  ///< [in] master vector (to be transferred)
        ) const;

    /*! \brief Transfer a dof vector from master to slave
     *
     *  \return Slave vector
     */
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<Epetra_MultiVector> mv  ///< [in] master vector (to be transferred)
        ) const;

    /*! \brief Transfer a dof vector from master to slave (const version)
     *
     *  \return Slave vector
     */
    Teuchos::RCP<Epetra_Vector> MasterToSlave(
        Teuchos::RCP<const Epetra_Vector> mv  ///< [in] master vector (to be transferred)
        ) const;

    /*! \brief Transfer a dof vector from master to slave (const version)
     *
     *  \return Slave vector
     */
    Teuchos::RCP<Epetra_MultiVector> MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv  ///< [in] master vector (to be transferred)
        ) const;

    /*! \brief Transfer a dof vector from slave to master
     *
     *  \return Master vector
     */
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<Epetra_Vector> sv  ///< [in] slave vector (to be transferred)
        ) const;

    /*! \brief Transfer a dof vector from slave to master
     *
     *  \return Master vector
     */
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<Epetra_MultiVector> sv  ///< [in] slave vector (to be transferred)
        ) const;

    /*! \brief Transfer a dof vector from slave to master (const version)
     *
     *  \return Master vector
     */
    Teuchos::RCP<Epetra_Vector> SlaveToMaster(
        Teuchos::RCP<const Epetra_Vector> sv  ///< [in] slave vector (to be transferred)
        ) const;

    /*! \brief Transfer a dof vector from slave to master (const version)
     *
     *  \return Master vector
     */
    Teuchos::RCP<Epetra_MultiVector> SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv  ///< [in] slave vector (to be transferred)
        ) const;

    /// transfer a dof vector from master to slave
    void MasterToSlave(
        Teuchos::RCP<const Epetra_MultiVector> mv,  ///< [in] master vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> sv         ///< [out] slave vector (containing result)
        ) const;

    /// transfer a dof vector from slave to master
    void SlaveToMaster(
        Teuchos::RCP<const Epetra_MultiVector> sv,  ///< [in] slave vector (to be transferred)
        Teuchos::RCP<Epetra_MultiVector> mv         ///< [out] master vector (containing result)
        ) const;

    //@}

    /** \name Coupled maps */
    //@{

    /// Get the interface dof row map of the master side
    Teuchos::RCP<const Epetra_Map> MasterDofMap() const { return pmasterdofrowmap_; }

    /// Get the interface dof row map of the slave side
    Teuchos::RCP<const Epetra_Map> SlaveDofMap() const { return pslavedofrowmap_; }

    //@}

    /** \name Condensation methods */
    //@{

    /// do condensation of Lagrange multiplier and slave-sided dofs
    void MortarCondensation(
        Teuchos::RCP<LINALG::SparseMatrix>& k,  ///< in:  tangent matrix w/o condensation
                                                ///< out: tangent matrix w/  condensation
        Teuchos::RCP<Epetra_Vector>& rhs        ///< in:  rhs vector     w/o condensation
                                                ///< out: rhs vector     w/  condensation
        ) const;

    /// recover slave-sided dofs
    void MortarRecover(Teuchos::RCP<LINALG::SparseMatrix>& k,  ///< in: tangent matrix
        Teuchos::RCP<Epetra_Vector>& inc  ///< in:  solution vector     w/o condensation
                                          ///< out: solution vector     w/  condensation
        ) const;

    //@}

   protected:
    /// Create mortar projection operator \f$P=D{^1}M\f$
    virtual void CreateP();

    /*! \brief Check if slave dofs have Dirichlet constraints
     *
     *  Slave DOFs are not allowed to carry Dirichlet boundary conditions to
     *  avoid over-constraining the problem, cf. [1]
     *
     *  <h3> References </h3>
     *  [1] Puso, M and Laursen, TA: Mesh tying on curved interfaces in 3D,
     *      Engineering Computation, 20:305-319 (2003)
     */
    void CheckSlaveDirichletOverlap(
        const Teuchos::RCP<DRT::Discretization>& slavedis,  ///< [in] Slave discretization
        const Epetra_Comm& comm                             ///< [in] Communicator
    );

    /// back transformation to initial parallel distribution
    void MatrixRowColTransform();

    /// check setup call
    const bool& IsSetup() const { return issetup_; };

    /// check init and setup call
    virtual void CheckSetup() const
    {
      if (!IsSetup()) dserror("ERROR: Call Setup() first!");
    }

   private:
    /// perform mesh relocation
    void MeshRelocation(Teuchos::RCP<DRT::Discretization> slavedis,  ///< [in] Slave discretization
        Teuchos::RCP<DRT::Discretization> aledis,                    ///< [in] ALE discretization
        Teuchos::RCP<const Epetra_Map>
            masterdofrowmap,  ///< [in] DOF row map of master discretization
        Teuchos::RCP<const Epetra_Map>
            slavedofrowmap,                  ///< [in] DOF row map of slave discretization
        Teuchos::RCP<Epetra_Vector>& idisp,  ///< [in] ALE displacements
        const Epetra_Comm& comm,             ///< [in] Communicator
        bool slavewithale                    ///< [in] Flag defining if slave is ALE
    );

   protected:
    /// Check for setup
    bool issetup_;

    /// Interface
    Teuchos::RCP<MORTAR::MortarInterface> interface_;

    /// Map of master row dofs (after parallel redist.)
    Teuchos::RCP<const Epetra_Map> masterdofrowmap_;

    /// Map of slave row dofs  (after parallel redist.)
    Teuchos::RCP<const Epetra_Map> slavedofrowmap_;

    /// Map of master row dofs (before parallel redist.)
    Teuchos::RCP<const Epetra_Map> pmasterdofrowmap_;

    /// Map of slave row dofs  (before parallel redist.)
    Teuchos::RCP<const Epetra_Map> pslavedofrowmap_;

    /// Slave side mortar matrix \f$D\f$
    Teuchos::RCP<LINALG::SparseMatrix> D_;

    /// Inverse \f$D^{-1}\f$ of slave side mortar matrix \f$D\f$
    Teuchos::RCP<LINALG::SparseMatrix> Dinv_;

    /// Master side mortar matrix \f$M\f$
    Teuchos::RCP<LINALG::SparseMatrix> M_;

    /// Mortar projection operator \f$P=D^{-1}M\f$
    Teuchos::RCP<LINALG::SparseMatrix> P_;
  };
}  // namespace ADAPTER

#endif /* ADAPTER_COUPLING_MORTAR_H */
