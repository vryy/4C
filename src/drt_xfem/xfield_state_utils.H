/*----------------------------------------------------------------------*/
/*!
\file xfield_state_utils.H

\brief Utils routines for xfluid state class

\level 0

<pre>
\maintainer  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>
 */
/*----------------------------------------------------------------------*/


#ifndef XFIELD_STATE_UTILS_H_
#define XFIELD_STATE_UTILS_H_


#include "../linalg/linalg_blocksparsematrix.H"
#include "../solver_nonlin_nox/nox_nln_aux.H"
#include <Teuchos_RCP.hpp>


namespace XFEM
{

/** \brief Destroy the LINALG::SparseOperator object and it's date
 *
 *  \author hiermeier
 *  \date 07/16 */
inline void DestroyMatrix(
    Teuchos::RCP<LINALG::SparseOperator> & mat,
    bool throw_exception = true
)
{
  // reference-counted object can be deleted by setting RCP = Teuchos::null when strong_count() == 1
  // given a weak RCP we do not have the permission to delete the reference-counted object
  // given a strong RCP with strong_count() > 1 we only can decrement the strong reference counter

  if (mat.strength() == Teuchos::RCP_STRONG) // strong RCP
  {
    if(mat.strong_count() == 1)
    {
      // which operator type do we have?
      NOX::NLN::LinSystem::OperatorType optype =
          NOX::NLN::AUX::GetOperatorType(*mat);
      // destroy underlying Epetra objects of the reference-counted object
      switch (optype)
      {
        case NOX::NLN::LinSystem::LinalgSparseMatrix:
        {
          Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(mat)->Destroy();
          break;
        }
        case NOX::NLN::LinSystem::LinalgBlockSparseMatrix:
        {
          Teuchos::rcp_dynamic_cast
              <LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> >(mat)->Destroy();
          break;
        }
        default:
        {
          dserror("The given LINALG::SparseOperator type is not supported! (%s)",
              NOX::NLN::LinSystem::OperatorType2String(optype).c_str());
          break;
        }
      }
      mat = Teuchos::null;    // destroy the rcp itself and delete the reference-counted object
    }
    else if (mat.strong_count() > 1)
    {
      if(throw_exception)
        dserror("could not destroy matrix object: %i!=1 pointers", mat.strong_count());
      else
        mat = Teuchos::null; // decrement the strong reference counter
    }
  }
  else if(mat.strength() == Teuchos::RCP_WEAK) // weak RCP
  {
    mat = Teuchos::null; // invalidate the RCP, reference-counted object won't be deleted by this weak pointer
  }
  else dserror("invalid strength of RCP");

}

/** \brief Destroy the LINALG::SparseMatrix object and it's date
 *
 *  \author schott
 *  \date 01/15 */
inline void DestroyMatrix(
    Teuchos::RCP<LINALG::SparseMatrix> & mat,
    bool throw_exception = true
)
{
  // reference-counted object can be deleted by setting RCP = Teuchos::null when strong_count() == 1
  // given a weak RCP we do not have the permission to delete the reference-counted object
  // given a strong RCP with strong_count() > 1 we only can decrement the strong reference counter

  if (mat.strength() == Teuchos::RCP_STRONG) // strong RCP
  {
    if(mat.strong_count() == 1)
    {
      mat->Destroy();         // destroy underlying Epetra objects of the reference-counted object
      mat = Teuchos::null;    // destroy the rcp itself and delete the reference-counted object
    }
    else if (mat.strong_count() > 1)
    {
      if(throw_exception)
        dserror("could not destroy matrix object: %i!=1 pointers", mat.strong_count());
      else
        mat = Teuchos::null; // decrement the strong reference counter
    }
  }
  else if(mat.strength() == Teuchos::RCP_WEAK) // weak RCP
  {
    mat = Teuchos::null; // invalidate the RCP, reference-counted object won't be deleted by this weak pointer
  }
  else dserror("invalid strength of RCP");

}


/** \brief Destroy the reference counted object and the reference counter
 *
 *  \author schott
 *  \date 01/15 */
template<class OBJECT>
inline void DestroyRCPObject(
    Teuchos::RCP<OBJECT> & obj_rcp,
    bool throw_exception = true
)
{
  // reference-counted object can be deleted by setting RCP = Teuchos::null when strong_count() == 1
  // given a weak RCP we do not have the permission to delete the reference-counted object
  // given a strong RCP with strong_count() > 1 we only can decrement the strong reference counter

  if(obj_rcp == Teuchos::null) return;

  if (obj_rcp.strength() == Teuchos::RCP_STRONG) // strong RCP
  {
    if(obj_rcp.strong_count() == 1)
    {
      obj_rcp = Teuchos::null;    // destroy the rcp itself and delete the reference-counted object
    }
    else if (obj_rcp.strong_count() > 1)
    {
      if(throw_exception)
        dserror("could not destroy reference-counted object: %i!=1 strong pointers", obj_rcp.strong_count());
      else
        obj_rcp = Teuchos::null; // decrement the strong reference counter
    }
  }
  else if(obj_rcp.strength() == Teuchos::RCP_WEAK) // weak RCP
  {
    obj_rcp = Teuchos::null; // invalidate the RCP, reference-counted object won't be deleted by this weak pointer
  }
  else dserror("invalid strength of RCP");

}


/** \brief More efficient and memory safe Zero routine for system matrix
 *
 *  \author schott
 *  \date 01/15 */
inline void ZeroMatrix( const Teuchos::RCP<LINALG::SparseMatrix> & mat )
{
  if( mat->ExplicitDirichlet() )
  {
    mat->Zero(); // matrix could have been changed due to Dirichlet conditions, go back to original Graph if savegraph == true
  }
  else
  {
    // do not create a new matrix via Zero() but zero entries
    mat->PutScalar(0.0);
  }
}

} // namespace XFEM

#endif /* XFIELD_STATE_UTILS_H_ */

