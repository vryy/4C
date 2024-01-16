/*----------------------------------------------------------------------*/
/*! \file

\brief Utils routines for xfluid state class

\level 0

 */
/*----------------------------------------------------------------------*/


#ifndef BACI_XFEM_XFIELD_STATE_UTILS_H
#define BACI_XFEM_XFIELD_STATE_UTILS_H


#include "baci_config.H"

#include "baci_linalg_blocksparsematrix.H"
#include "baci_solver_nonlin_nox_aux.H"

#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


namespace XFEM
{
  /** \brief Destroy the CORE::LINALG::SparseOperator object and it's date
   *
   *  \author hiermeier
   *  \date 07/16 */
  inline void DestroyMatrix(
      Teuchos::RCP<CORE::LINALG::SparseOperator>& mat, bool throw_exception = true)
  {
    // reference-counted object can be deleted by setting RCP = Teuchos::null when strong_count() ==
    // 1 given a weak RCP we do not have the permission to delete the reference-counted object given
    // a strong RCP with strong_count() > 1 we only can decrement the strong reference counter

    if (mat.strength() == Teuchos::RCP_STRONG)  // strong RCP
    {
      if (mat.strong_count() == 1)
      {
        // which operator type do we have?
        NOX::NLN::LinSystem::OperatorType optype = NOX::NLN::AUX::GetOperatorType(*mat);
        // destroy underlying Epetra objects of the reference-counted object
        switch (optype)
        {
          case NOX::NLN::LinSystem::LinalgSparseMatrix:
          {
            Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(mat)->Destroy();
            break;
          }
          case NOX::NLN::LinSystem::LinalgBlockSparseMatrix:
          {
            Teuchos::RCP<CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>>
                block_mat = Teuchos::rcp_dynamic_cast<
                    CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>>(
                    mat, true);
            block_mat->Destroy(false);
            break;
          }
          default:
          {
            std::stringstream msg;
            msg << "The given CORE::LINALG::SparseOperator type is not supported! ( "
                << NOX::NLN::LinSystem::OperatorType2String(optype) << " )";
            dserror(msg.str());
            break;
          }
        }
        mat = Teuchos::null;  // destroy the rcp itself and delete the reference-counted object
      }
      else if (mat.strong_count() > 1)
      {
        if (throw_exception)
        {
          std::stringstream msg;
          msg << "Could not destroy matrix object: " << mat.strong_count() << "!=1 pointers";
          dserror(msg.str());
        }
        else
          mat = Teuchos::null;  // decrement the strong reference counter
      }
    }
    else if (mat.strength() == Teuchos::RCP_WEAK)  // weak RCP
    {
      mat = Teuchos::null;  // invalidate the RCP, reference-counted object won't be deleted by this
                            // weak pointer
    }
    else
    {
      std::stringstream msg;
      msg << "invalid strength of RCP: "
          << Teuchos::ToStringTraits<Teuchos::ERCPStrength>::toString(mat.strength());
      dserror(msg.str());
    }
  }

  /** \brief Destroy the CORE::LINALG::SparseMatrix object and it's date
   *
   *  \author schott
   *  \date 01/15 */
  inline void DestroyMatrix(
      Teuchos::RCP<CORE::LINALG::SparseMatrix>& mat, bool throw_exception = true)
  {
    // reference-counted object can be deleted by setting RCP = Teuchos::null when strong_count() ==
    // 1 given a weak RCP we do not have the permission to delete the reference-counted object given
    // a strong RCP with strong_count() > 1 we only can decrement the strong reference counter

    if (mat.strength() == Teuchos::RCP_STRONG)  // strong RCP
    {
      if (mat.strong_count() == 1)
      {
        mat->Destroy();       // destroy underlying Epetra objects of the reference-counted object
        mat = Teuchos::null;  // destroy the rcp itself and delete the reference-counted object
      }
      else if (mat.strong_count() > 1)
      {
        if (throw_exception)
        {
          std::ostringstream msg;
          msg << "Could not destroy matrix object: " << mat.strong_count() << "!=1 pointers";
          dserror(msg.str());
        }
        else
          mat = Teuchos::null;  // decrement the strong reference counter
      }
    }
    else if (mat.strength() == Teuchos::RCP_WEAK)  // weak RCP
    {
      mat = Teuchos::null;  // invalidate the RCP, reference-counted object won't be deleted by this
                            // weak pointer
    }
    else
      dserror("invalid strength of RCP");
  }


  /** \brief Destroy the reference counted object and the reference counter
   *
   *  \author schott
   *  \date 01/15 */
  template <class OBJECT>
  inline void DestroyRCPObject(Teuchos::RCP<OBJECT>& obj_rcp, bool throw_exception = true)
  {
    // reference-counted object can be deleted by setting RCP = Teuchos::null when strong_count() ==
    // 1 given a weak RCP we do not have the permission to delete the reference-counted object given
    // a strong RCP with strong_count() > 1 we only can decrement the strong reference counter

    if (obj_rcp == Teuchos::null) return;

    if (obj_rcp.strength() == Teuchos::RCP_STRONG)  // strong RCP
    {
      if (obj_rcp.strong_count() == 1)
      {
        obj_rcp = Teuchos::null;  // destroy the rcp itself and delete the reference-counted object
      }
      else if (obj_rcp.strong_count() > 1)
      {
        if (throw_exception)
        {
          std::ostringstream msg;
          msg << "Could not destroy reference-counted object! strong_count() = "
              << obj_rcp.strong_count() << " (!=1 pointers)";
          dserror(msg.str());
        }
        else
          obj_rcp = Teuchos::null;  // decrement the strong reference counter
      }
    }
    else if (obj_rcp.strength() == Teuchos::RCP_WEAK)  // weak RCP
    {
      obj_rcp = Teuchos::null;  // invalidate the RCP, reference-counted object won't be deleted by
                                // this weak pointer
    }
    else
      dserror("invalid strength of RCP");
  }


  /** \brief More efficient and memory safe Zero routine for system matrix
   *
   *  \author schott
   *  \date 01/15 */
  inline void ZeroMatrix(const Teuchos::RCP<CORE::LINALG::SparseMatrix>& mat)
  {
    if (mat->ExplicitDirichlet())
    {
      mat->Zero();  // matrix could have been changed due to Dirichlet conditions, go back to
                    // original Graph if savegraph == true
    }
    else
    {
      // do not create a new matrix via Zero() but zero entries
      mat->PutScalar(0.0);
    }
  }

}  // namespace XFEM

BACI_NAMESPACE_CLOSE

#endif  // XFEM_XFIELD_STATE_UTILS_H
