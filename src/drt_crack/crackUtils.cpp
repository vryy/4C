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
#include "crack_tolerance.H"

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

/*--------------------------------------------------------------------------------------*
 * Find angle of a vector pointing from point a to b and horizontal plane       sudhakar 04/14
 * The calculated angle is in the range [0,2pi]
 *--------------------------------------------------------------------------------------*/
double DRT::CRACK::UTILS::FindAngle( const double a[3], const double b[3] )
{
  double normal[3];
  normal[0] = b[0] - a[0];
  normal[1] = b[1] - a[1];
  normal[2] = 0.0;

  double angle = atan2( normal[1], normal[0] );

  convertAngleTo_02PI_range( angle );

  return angle;
}

/*-----------------------------------------------------------------------------------------*
 * convert the angle in range [-pi,pi] to [0,2pi]                                 sudhakar 04/14
 * This function is usually called after using atan2() to find angle
 *-----------------------------------------------------------------------------------------*/
void DRT::CRACK::UTILS::convertAngleTo_02PI_range( double & ang )
{
  double PI = 22.0/7.0;
  if( (fabs(ang) < ANGLE_TOL_ZERO) or (fabs(ang - 2.0*PI) < ANGLE_TOL_ZERO) )
    ang = 0.0;
  else if( ang < 0.0 )
    ang = 2.0 * PI + ang;
}

/*-----------------------------------------------------------------------------------------*
 * Returns true if the given element has given nodeid                                 sudhakar 04/14
 *-----------------------------------------------------------------------------------------*/
bool DRT::CRACK::UTILS::ElementHasThisNodeId( const DRT::Element* ele, int nodeid )
{
  bool found = false;

  const int* nodes = ele->NodeIds();
  const int numNodes = ele->NumNode();

  for( int iNode = 0; iNode < numNodes; iNode++ )
  {
    int surnodeid = nodes[iNode];

    if( surnodeid == nodeid )
    {
      found = true;
      break;
    }
  }
  return found;
}

