/*-----------------------------------------------------------*/
/*!
\file biopolynet_calc_utils.H

\brief utils for biopolymer network business

\maintainer Jonas Eichinger, Maximilian Grill

\date August, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BIOPOLYNET_CALC_UTILS_H_
#define BIOPOLYNET_CALC_UTILS_H_


#include "../drt_lib/drt_globalproblem.H"
#include "../drt_beamcontact/beam3contact_utils.H"

namespace BIOPOLYNET
{
/*!
\brief namespace of utility functions for statmech methods

*/
namespace UTILS
{

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void PeriodicBoundaryShift1D(
    double&       val,
    const double& periodlength,
    const double& initval = 0.0)
{
  // in case we want to manipulate the displacement only, we add the
  // constant reference value here and leave it untouched
  const double abspos = val + initval;

  if (periodlength>0.0)
  {
    /* add or subtract the period length sufficiently often
     * to lie again in the domain [0,periodlength] */
    if (abspos < 0.0)
      val -= periodlength * std::floor(abspos/periodlength);  // minus sign because abspos is negative
    else if (abspos > periodlength)
      val -= periodlength * std::floor(abspos/periodlength);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void PeriodicBoundaryUnShift1D(
    double&       val,
    const double& abspos_refpoint,
    const double& periodlength,
    const double& initval = 0.0)
{
  // in case we want to manipulate the displacement only, we add the
  // constant reference value here and leave it untouched
  const double abspos = val + initval;

  /* Check the distance to a reference point position (e.g. node of the
   * same element). If the distance is larger than half of the period
   * length, the point position has been shifted before.
   *
   * Warning: This assumes that the distance between point and reference
   *          point is not larger than half of the period length unless we
   *          shift it. For beam elements, this restricts the element length
   *          to be smaller than this value throughout the entire simulation.
   *          So far, we only check this once in the beginning.
   *
   * Note: this should be equivalent to the previously applied criterion
   *       that the distance between given point and reference point
   *       decreases by either adding or subtracting the period length. */
  if (periodlength>0.0)
  {
    if (abspos - abspos_refpoint < -0.5 * periodlength)
      val += periodlength;

    else if(abspos - abspos_refpoint > 0.5 * periodlength)
      val -= periodlength;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void PeriodicBoundaryShift(
    LINALG::Matrix<3,1>&       pos,
    const std::vector<double>& periodlength)
{
  for (unsigned int idim=0; idim<3; ++idim)
    PeriodicBoundaryShift1D(pos(idim),periodlength[idim]);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void PeriodicBoundaryUnShift(
    LINALG::Matrix<3,1>&       pos,
    const LINALG::Matrix<3,1>& pos_refpoint,
    const std::vector<double>& periodlength)
{
  for (unsigned int idim=0; idim<3; ++idim)
    PeriodicBoundaryUnShift1D(pos(idim),pos_refpoint(idim),periodlength[idim]);
}

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
      PeriodicBoundaryShift1D((*dis)[doflid],(*periodlength)[j],node->X()[j]);
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
}//namespace BIOPOLYNET

#endif /* BIOPOLYNET_CALC_UTILS_H_ */
