/*-----------------------------------------------------------*/
/*!
\file biopolynet_calc_utils.H

\brief utils for biopolymer network business

\maintainer Jonas Eichinger

\date August, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BIOPOLYNET_CALC_UTILS_H_
#define BIOPOLYNET_CALC_UTILS_H_


#include "../drt_lib/drt_globalproblem.H"
#include "../drt_beamcontact/beam3contact_utils.H"

namespace STATMECH
{
/*!
\brief namespace of utility functions for statmech methods

*/
namespace UTILS
{

/*----------------------------------------------------------------------------*
 | Shift dis according to periodic boundary conditions         eichinger 08/16|
 *----------------------------------------------------------------------------*/
inline void PeriodicBoundaryConsistentDis(Teuchos::RCP<Epetra_Vector> dis,
    const Teuchos::RCP<const std::vector<double> >             periodlength,
    const Teuchos::RCP<const DRT::Discretization>              discret)
{
  /* In statmech simulations with periodic boundary conditions we need to
   * ensure that during the whole simulation, including each iteration state,
   * no current node position (x = X + dis) lies outside the periodic bounding
   * box. This is necessary for crosslinking, dirichlet, beamcontact and
   * potential to work correctly. For element evaluation, the unshifted
   * configuration is calculated on (beam)element level.
   *
   *  \author Jonas Eichinger                                         08/16  */


  // get problem dimension
  const int ndim = DRT::Problem::Instance()->NDim();

  for (int i=0; i<discret->NumMyRowNodes(); i++)
  {
    //get a pointer at i-th row node
    DRT::Node* node = discret->lRowNode(i);

    /* Hermite Interpolation: Check whether node is a beam node which is NOT
     * used for centerline interpolation if so, we simply skip it because
     * it does not have position DoFs */
    if (BEAMCONTACT::BeamNode(*node) and not BEAMCONTACT::BeamCenterlineNode(*node))
      continue;

    //get GIDs of this node's degrees of freedom
    std::vector<int> dofnode = discret->Dof(node);

    for (int j=ndim-1; j>-1; j--)
    {
      int doflid = dis->Map().LID(dofnode[j]);
      //current coordinate value
      double xcurr = node->X()[j] + (*dis)[doflid];

      /*if node currently has coordinate value greater than periodlength, it is
       *shifted by -periodlength sufficiently often to lie again in the domain*/
      if (xcurr > (*periodlength)[j])
      {
        (*dis)[doflid] -= (*periodlength)[j]*floor(xcurr/(*periodlength)[j]);
      }
      /*if node currently has coordinate value smaller than zero, it is shifted
       * by periodlength sufficiently often to lie again in the domain*/
      else if (xcurr < 0.0)
      {
        (*dis)[doflid] -= (*periodlength)[j]*floor(xcurr/(*periodlength)[j]);
      }
    }
  }

  // that's it
  return;

}; // PeriodicBoundaryConsistentDis()

/*----------------------------------------------------------------------*
 | generates a vector with a random permutation of the numbers between 0|
 | and N - 1                                        (public) cyron 06/10|
 *----------------------------------------------------------------------*/
inline std::vector<int> Permutation(const int& number)
{
  //auxiliary variable
  int j = 0;

  DRT::Problem::Instance()->Random()->SetRandRange(0.0,1.0);

  //result vector initialized with ordered numbers from 0 to N-1
  std::vector<int> randorder(number, 0);
  for (int i=0; i<(int)randorder.size(); i++)
    randorder[i] = i;

  for (int i=0; i<number; ++i)
  {
    //generate random number between 0 and i
    j = (int)floor((i + 1.0)*DRT::Problem::Instance()->Random()->Uni());

    /*exchange values at positions i and j (note: value at position i is i due to above initialization
     *and because so far only positions <=i have been changed*/
    randorder[i] = randorder[j];
    randorder[j] = i;
  }

  return randorder;
}; // Permutation()


}//namespace UTILS
}//namespace STATMECH

#endif /* BIOPOLYNET_CALC_UTILS_H_ */
