/*----------------------------------------------------------------------*/
/*!
\file crackUtils.cpp

\brief Utility functions for crack propagation problem

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#include "crackUtils.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_utils.H"

/*----------------------------------------------------------------------------------------*
 * When a new node is introduced in the discret, DOFs corresponding to the new node are zeros.
 * We copy the field values for this node from the old node that is already         sudhakar 03/14
 * existing at the same position
 *----------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::UpdateThisEpetraVectorCrack( const Teuchos::RCP<DRT::Discretization>& discret,
                                                     Teuchos::RCP<Epetra_Vector>& vec,
                                                     const std::map<int,int>& oldnewIds )
{
  Teuchos::RCP<Epetra_Vector> old = vec;
  vec = LINALG::CreateVector(*discret->DofRowMap(),true);

  LINALG::Export( *old, *vec );

  for( std::map<int,int>::const_iterator it = oldnewIds.begin(); it != oldnewIds.end(); it++ )
  {
    int oldid = it->first;
    int newid = it->second;

    DRT::UTILS::EquateValuesAtTheseNodes( *vec, discret, oldid, newid );
  }
}
