/*----------------------------------------------------------------------*/
/*!
\file xfluid_utils.H

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


#ifndef XFLUID_UTILS_H_
#define XFLUID_UTILS_H_


#include "../linalg/linalg_sparsematrix.H"
#include <Teuchos_RCP.hpp>


namespace FLD
{

/*----------------------------------------------------------------------*
 |  destroy the LINALG::SparseMatrix object and it's date  schott 01/15 |
 *----------------------------------------------------------------------*/
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


/*----------------------------------------------------------------------*
 |  destroy the reference counted object and the reference counter      |
 |                                                         schott 01/15 |
 *----------------------------------------------------------------------*/
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


/*----------------------------------------------------------------------*
 |  more efficient and memory safe Zero routine for system matrix       |
 |                                                         schott 01/15 |
 *----------------------------------------------------------------------*/
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

} // namespace FLD

#endif /* XFLUID_UTILS_H_ */

