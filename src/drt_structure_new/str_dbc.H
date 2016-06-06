/*-----------------------------------------------------------*/
/*!
\file str_dbc.H

\brief Wrapper for all Dirichlet boundary condition tasks.

\maintainer Michael Hiermeier

\date Nov 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_DBC_H_
#define SRC_DRT_STRUCTURE_NEW_STR_DBC_H_

#include <Teuchos_RCP.hpp>
#include "../solver_nonlin_nox/nox_nln_abstract_prepostoperator.H"

// forward declaration ...
class Epetra_Vector;
namespace DRT {
  class Discretization;
namespace UTILS {
  class LocsysManager;
} // namespace UTILS
} // namespace DRT
namespace LINALG {
  class SparseOperator;
  class SparseMatrix;
  class MapExtractor;
} // namespace LINALG
namespace STR {
namespace TIMINT {
  class Base;
} // namespace TIMINT
  class Dbc
  {
  public:
    //! Constructor
    Dbc();

    //! Destructor
    virtual ~Dbc() {};

    //! Initialize class variables
    virtual void Init(const Teuchos::RCP<DRT::Discretization>& discret,
        const Teuchos::RCP<Epetra_Vector>& freact,
        const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr);

    //! Setup class variables
    virtual void Setup();

    /*! \brief Apply the DBC to system of equations
     *
     *  Note: Stay in the local coordinate system and do not rotate back (if locSys is defined).*/
    void ApplyDirichletToLocalSystem(Teuchos::RCP<LINALG::SparseOperator> A,
        Teuchos::RCP<Epetra_Vector>& b) const;

    /*! \brief Apply the DBC to a vector
     *
     *  Note: Stay in the global coordinate system (Rotation: global-->local-->global).*/
    void ApplyDirichletToVector(Teuchos::RCP<Epetra_Vector>& vec) const;

    /*! \brief Apply the DBC to the rhs vector and calculate and save the reaction forces
     *
     *  Note: Stay in the global coordinate system (Rotation: global-->local-->global).*/
    void ApplyDirichletToRhs(Teuchos::RCP<Epetra_Vector>& b);

    //! Update the locsys manager
    void UpdateLocSysManager();

    //! Calculate the dirichlet increment of the current (time) step
    Teuchos::RCP<Epetra_Vector> GetDirichletIncrement();

    /*! \brief Evaluate and apply the DBC
     *
     * Note: Stay in the global coordinate system (Rotation: global-->local-->global).*/
    void ApplyDirichletBC(
        const double& time,
        Teuchos::RCP<Epetra_Vector> dis,
        Teuchos::RCP<Epetra_Vector> vel,
        Teuchos::RCP<Epetra_Vector> acc,
        bool recreatemap);

    //! Insert non-dbc dof values of the source vector into the non-dbc dofs of the target vector
    void InsertVectorInNonDbcDofs(Teuchos::RCP<const Epetra_Vector> source_ptr,
        Teuchos::RCP<Epetra_Vector> target_ptr) const;

    //! @name Access functions
    //!@{
    //! Get the Dirichlet Boundary Condition map extractor
    Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const;

    //! Get a pointer to the local system manager
    Teuchos::RCP<DRT::UTILS::LocsysManager> LocSysManagerPtr();

    //! Get the zeros vector
    const Epetra_Vector& GetZeros() const;
    Teuchos::RCP<const Epetra_Vector> GetZerosPtr() const;
    //!@}

  protected:
    //! Returns the initialization status
    const bool& IsInit() const
    { return isinit_; };

    //! Returns the setup status
    const bool& IsSetup() const
    { return issetup_; };

    //! Checks the initialization status
    void CheckInit() const;

    //! Checks the initialization and setup status
    void CheckInitSetup() const;

    /*! Rotate the system matrix from a global to a local coordinate system,
     *  if the locsysman_ is defined.
     *
     *  Works only for LINALG::SparseMatrices! */
    bool RotateGlobalToLocal(Teuchos::RCP<LINALG::SparseOperator> A) const;

    /*! Rotate the rhs vector from the global to the local coordinate system, if
     *  the locsysman_ is defined. */
    bool RotateGlobalToLocal(Teuchos::RCP<Epetra_Vector> v) const;
    bool RotateGlobalToLocal(Teuchos::RCP<Epetra_Vector> v, bool offset) const;

    /*! Rotate a vector from the local to the global coordinate system, if
     *  the locsysman_ is defined. */
    bool RotateLocalToGlobal(Teuchos::RCP<Epetra_Vector> v) const;
    bool RotateLocalToGlobal(Teuchos::RCP<Epetra_Vector> v, bool offset) const;

    //! Get the locsys transformation matrix
    Teuchos::RCP<const LINALG::SparseMatrix> GetLocSysTrafo() const;

    const bool& IsLocSys() const
    { return islocsys_; };

    //! Extract the reaction forces
    void ExtractFreact(Teuchos::RCP<Epetra_Vector>& b) const;

    /*! Apply the DBC to the right hand side in the local coordinate system and
     *  do not rotate it back to the global coordinate system. */
    void ApplyDirichletToLocalRhs(Teuchos::RCP<Epetra_Vector>& b) const;

    /*! Apply the DBC to the jacobian in the local coordinate system and
     *  do not rotate it back to the global coordinate system. */
    void ApplyDirichletToLocalJacobian(Teuchos::RCP<LINALG::SparseOperator> A) const;

  protected:
    //! Flag indicating the initialization status.
    bool isinit_;

    //! Flag indicating the setup status.
    bool issetup_;

    //! Flag indicating if a locsysman_ was defined.
    bool islocsys_;

    //! Discretization pointer
    Teuchos::RCP<DRT::Discretization> discret_ptr_;

    //! Reaction force
    Teuchos::RCP<Epetra_Vector> freact_ptr_;

    //! pointer to the overlying time integrator (read-only)
    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;

    //! Pointer to the local coordinate system manager
    Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_ptr_;

    //! Some vector with system size and filled with zeros.
    Teuchos::RCP<Epetra_Vector> zeros_ptr_;

    //! Dirichlet boundary condition map extractor.
    Teuchos::RCP<LINALG::MapExtractor> dbcmap_ptr_;
  };  // class  Dbc
} // namespace STR

namespace NOX {
namespace NLN {
namespace LinSystem {
namespace PrePostOp {
  /*! \brief PrePostOperator class to modify the linear system before the linear system is going to be solved.
   *
   * We use this pre/post operator to apply the DBC on the linear system of equations before the
   * linear system is going to be solved. This gives us the opportunity to rotate the matrix only
   * once, if locSys is defined and to apply possible modifications to the linear system at different
   * places without the need to re-apply the DBC (see PTC for an example).
   *
   * \author Hiermeier */
  class Dbc : public NOX::NLN::Abstract::PrePostOperator
  {
  public:
    //! constructor
    Dbc(const Teuchos::RCP<const ::STR::Dbc>& dbc_ptr);

    //! destructor
    virtual ~Dbc() {};

    //! \brief Apply the DBC and rotate the system of equations if necessary (derived)
    virtual void runPreApplyJacobianInverse(
        NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys);

  private:
    //! pointer to the underlying class, which provides the whole functionality
    Teuchos::RCP<const ::STR::Dbc> dbc_ptr_;
  };  // class Dbc (pre/post operator)
} // namespace PrePostOp
} // namespace LinSystem
} // namespace NLN
} // namespace NOX


#endif /* SRC_DRT_STRUCTURE_NEW_STR_DBC_H_ */
