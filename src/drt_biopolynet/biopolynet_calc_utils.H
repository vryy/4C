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
#include "../drt_beamcontact/beam_to_beam_linkage.H"
#include "../drt_beam3/beam3_base.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_biopolynet/periodic_boundingbox.H"

#include <Epetra_FEVector.h>

namespace BIOPOLYNET
{
namespace UTILS
{

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
inline void PeriodicBoundaryConsistentDisVector(
    Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    Teuchos::RCP<DRT::Discretization> const&        discret)
{
  for (int i=0; i<discret->NumMyRowNodes(); ++i)
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

    for (int dim=0; dim<3; ++dim)
    {
      int doflid = dis->Map().LID(dofnode[dim]);
      pbb->Shift1D( dim, (*dis)[doflid], node->X()[dim] );
    }
  }
}

/*----------------------------------------------------------------------*
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
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
inline void GetCurrentElementDis(
    DRT::Discretization const& discret,
    DRT::Element const* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    std::vector<double>& eledisp)
{
  // clear
  eledisp.clear();

  std::vector<int> lm, lmowner, lmstride;

  ele->LocationVector(discret,lm,lmowner,lmstride);
  DRT::UTILS::ExtractMyValues(*ia_discolnp,eledisp,lm);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
inline void GetPosAndTriadOfBindingSpot(
    DRT::Element* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    int const locbspotnum,
    LINALG::Matrix<3,1>& bspotpos,
    LINALG::Matrix<3,3>& bspottriad,
    std::vector<double>& eledisp)
{
  // cast to beambase element
  DRT::ELEMENTS::Beam3Base* beamele =
      dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele);

#ifdef DEBUG
      if(beamele == NULL)
        dserror("Dynamic cast to beam3base failed");
#endif

  // get current position at binding spot xi
  beamele->GetPosOfBindingSpot(bspotpos,eledisp,locbspotnum,pbb);

  // get current triad at binding spot xi
  beamele->GetTriadOfBindingSpot(bspottriad,eledisp,locbspotnum);

}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
inline void GetPosAndTriadOfBindingSpot(
    DRT::Discretization const& discret,
    DRT::Element* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    int const locbspotnum,
    LINALG::Matrix<3,1>& bspotpos,
    LINALG::Matrix<3,3>& bspottriad)
{
  std::vector<double> eledisp;
  GetCurrentElementDis(discret, ele, ia_discolnp, eledisp);

  GetPosAndTriadOfBindingSpot(ele, ia_discolnp, pbb, locbspotnum,
      bspotpos, bspottriad, eledisp);

}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
inline void FEAssembleEleForceStiffIntoSystemVectorMatrix(
  const DRT::Discretization&         discret,
  const int                          elegid1,
  const int                          elegid2,
  const Epetra_SerialDenseVector&    elevec1,
  const Epetra_SerialDenseVector&    elevec2,
  const Epetra_SerialDenseMatrix&    elemat11,
  const Epetra_SerialDenseMatrix&    elemat12,
  const Epetra_SerialDenseMatrix&    elemat21,
  const Epetra_SerialDenseMatrix&    elemat22,
  Teuchos::RCP<Epetra_FEVector>      fe_sysvec,
  Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat
  )
{
  // the entries of elevecX  belong to the Dofs of the element with GID elegidX
  // the rows    of elematXY belong to the Dofs of the element with GID elegidX
  // the columns of elematXY belong to the Dofs of the element with GID elegidY
  const DRT::Element* ele1 = discret.gElement(elegid1);
  const DRT::Element* ele2 = discret.gElement(elegid2);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->LocationVector(discret,lmrow1,lmrowowner1,lmstride);
  ele2->LocationVector(discret,lmrow2,lmrowowner2,lmstride);

  // assemble both element vectors into global system vector
  if( fe_sysvec != Teuchos::null )
  {
    fe_sysvec->SumIntoGlobalValues(elevec1.Length(),&lmrow1[0],elevec1.Values());
    fe_sysvec->SumIntoGlobalValues(elevec2.Length(),&lmrow2[0],elevec2.Values());
  }

  // and finally also assemble stiffness contributions
  if( fe_sysmat != Teuchos::null )
  {
    fe_sysmat->FEAssemble(elemat11,lmrow1,lmrow1);
    fe_sysmat->FEAssemble(elemat12,lmrow1,lmrow2);
    fe_sysmat->FEAssemble(elemat21,lmrow2,lmrow1);
    fe_sysmat->FEAssemble(elemat22,lmrow2,lmrow2);
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
inline void AssembleCenterlineDofForceStiffIntoElementForceStiff(
    const DRT::Discretization&      discret,
    const int                       elegid1,
    const int                       elegid2,
    const Epetra_SerialDenseVector& ele1force_centerlineDOFs,
    const Epetra_SerialDenseVector& ele2force_centerlineDOFs,
    const Epetra_SerialDenseMatrix& ele11stiff_centerlineDOFs,
    const Epetra_SerialDenseMatrix& ele12stiff_centerlineDOFs,
    const Epetra_SerialDenseMatrix& ele21stiff_centerlineDOFs,
    const Epetra_SerialDenseMatrix& ele22stiff_centerlineDOFs,
    Epetra_SerialDenseVector*       ele1force,
    Epetra_SerialDenseVector*       ele2force,
    Epetra_SerialDenseMatrix*       ele11stiff,
    Epetra_SerialDenseMatrix*       ele12stiff,
    Epetra_SerialDenseMatrix*       ele21stiff,
    Epetra_SerialDenseMatrix*       ele22stiff
  )
{
  // Todo generalize this method and make the ifs nicer

  // assemble centerline DOF values correctly into element DOFvec vectors/matrices
  for (unsigned int idof=0; idof<3; ++idof)
  {
    if (ele1force != NULL)
    {
      (*ele1force)(idof) = ele1force_centerlineDOFs(idof);
      (*ele1force)(6+idof) = ele1force_centerlineDOFs(3+idof);
    }

    if (ele2force != NULL)
    {
      (*ele2force)(idof) = ele2force_centerlineDOFs(idof);
      (*ele2force)(6+idof) = ele2force_centerlineDOFs(3+idof);
    }

    for (unsigned int jdof=0; jdof<3; ++jdof)
    {
      if (ele11stiff != NULL)
      {
        (*ele11stiff)(idof,jdof) = ele11stiff_centerlineDOFs(idof,jdof);
        (*ele11stiff)(6+idof,jdof) = ele11stiff_centerlineDOFs(3+idof,jdof);
        (*ele11stiff)(idof,6+jdof) = ele11stiff_centerlineDOFs(idof,3+jdof);
        (*ele11stiff)(6+idof,6+jdof) = ele11stiff_centerlineDOFs(3+idof,3+jdof);
      }

      if (ele11stiff != NULL)
      {
        (*ele12stiff)(idof,jdof) = ele12stiff_centerlineDOFs(idof,jdof);
        (*ele12stiff)(6+idof,jdof) = ele12stiff_centerlineDOFs(3+idof,jdof);
        (*ele12stiff)(idof,6+jdof) = ele12stiff_centerlineDOFs(idof,3+jdof);
        (*ele12stiff)(6+idof,6+jdof) = ele12stiff_centerlineDOFs(3+idof,3+jdof);
      }

      if (ele11stiff != NULL)
      {
        (*ele21stiff)(idof,jdof) = ele21stiff_centerlineDOFs(idof,jdof);
        (*ele21stiff)(6+idof,jdof) = ele21stiff_centerlineDOFs(3+idof,jdof);
        (*ele21stiff)(idof,6+jdof) = ele21stiff_centerlineDOFs(idof,3+jdof);
        (*ele21stiff)(6+idof,6+jdof) = ele21stiff_centerlineDOFs(3+idof,3+jdof);
      }

      if (ele11stiff != NULL)
      {
        (*ele22stiff)(idof,jdof) = ele22stiff_centerlineDOFs(idof,jdof);
        (*ele22stiff)(6+idof,jdof) = ele22stiff_centerlineDOFs(3+idof,jdof);
        (*ele22stiff)(idof,6+jdof) = ele22stiff_centerlineDOFs(idof,3+jdof);
        (*ele22stiff)(6+idof,6+jdof) = ele22stiff_centerlineDOFs(3+idof,3+jdof);
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
inline void ApplyBpotForceStiffToParentElements(
    DRT::Discretization const&                             discret,
    const Teuchos::RCP<Epetra_Vector>                      disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr,
    const LINALG::SerialDenseVector& bspotforce1,
    const LINALG::SerialDenseVector& bspotforce2,
    const LINALG::SerialDenseMatrix& bspotstiff11,
    const LINALG::SerialDenseMatrix& bspotstiff12,
    const LINALG::SerialDenseMatrix& bspotstiff21,
    const LINALG::SerialDenseMatrix& bspotstiff22,
    LINALG::SerialDenseVector*       ele1force,
    LINALG::SerialDenseVector*       ele2force,
    LINALG::SerialDenseMatrix*       ele11stiff,
    LINALG::SerialDenseMatrix*       ele12stiff,
    LINALG::SerialDenseMatrix*       ele21stiff,
    LINALG::SerialDenseMatrix*       ele22stiff)
{
  // apply force on binding spots and corresponding linearizations to parent elements
  // and get discrete element force vectors and stiffness matrices

  /* Of course, the underlying idea is much more general:
   * Consider some kind of abstract interaction between two points \xi_1,\xi_2 \in [-1,1]
   * on two beam elements (or the same beam element).
   * We computed a discrete force/moment exerted on each of the two points; in addition,
   * we computed the linearizations with respect to the position/rotation vector of
   * this point and also the linearizations with respect to the position/rotation vector
   * of the other point (in general, the forces/moments depend on position/rotation
   * vector of both points).
   * Starting from here, we calculate the generalized interpolation matrices for
   * variations and increments of the position/rotation vectors at the two points
   * by expressing their dependency on the primary nodal DoFs of the corresponding
   * element (this is element specific, of course). Finally, we can 'transform' the
   * force vectors and stiffness matrices to discrete element force vectors and
   * stiffness matrices */

  DRT::Element* ele1 = discret.gElement(elepairptr->GetEleGid(0));
  DRT::Element* ele2 = discret.gElement(elepairptr->GetEleGid(1));

  DRT::ELEMENTS::Beam3Base* beamele1 =
      dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele1);
  DRT::ELEMENTS::Beam3Base* beamele2 =
      dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele2);

  // get current element displacements
  std::vector<double> ele1disp;
  GetCurrentElementDis(discret, ele1,disp_np_col,ele1disp);
  const int numdof_ele1 = ele1disp.size();

  std::vector<double> ele2disp;
  GetCurrentElementDis(discret,ele2,disp_np_col,ele2disp);
  const int numdof_ele2 = ele2disp.size();


  // transformation matrix, will be resized and reused for various needs
  LINALG::SerialDenseMatrix trafomatrix;

  // auxiliary matrices required to store intermediate results after first of two
  // consecutive matrix-matrix products for each of four stiffness matrices
  LINALG::SerialDenseMatrix auxmat11;
  LINALG::SerialDenseMatrix auxmat12;
  LINALG::SerialDenseMatrix auxmat21;
  LINALG::SerialDenseMatrix auxmat22;

  // zero out and set correct size of transformation matrix
  trafomatrix.Shape(6,numdof_ele1);


  // element 1: I_variations
  beamele1->GetGeneralizedInterpolationMatrixVariationsAtXi(
      trafomatrix,
      beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)));

  if (ele1force != NULL)
  {
    ele1force->Size(numdof_ele1);
    ele1force->Multiply('T','N',1.0,trafomatrix,bspotforce1,0.0);
  }

  if (ele11stiff != NULL)
  {
    auxmat11.Shape(numdof_ele1,6);
    auxmat11.Multiply('T','N',1.0,trafomatrix,bspotstiff11,0.0);
  }

  if (ele12stiff != NULL)
  {
    auxmat12.Shape(numdof_ele1,6);
    auxmat12.Multiply('T','N',1.0,trafomatrix,bspotstiff12,0.0);
  }


  // element 1: I_increments
  if (ele11stiff != NULL or ele21stiff != NULL)
  {
    trafomatrix.Shape(6,numdof_ele1);

    beamele1->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomatrix,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)),
        ele1disp);


    if (ele11stiff != NULL)
    {
      ele11stiff->Shape(numdof_ele1,numdof_ele1);
      ele11stiff->Multiply('N','N',1.0,auxmat11,trafomatrix,0.0);
    }

    if (ele21stiff != NULL)
    {
      auxmat21.Shape(6,numdof_ele1);
      auxmat21.Multiply('N','N',1.0,bspotstiff21,trafomatrix,0.0);
    }
  }


  // element 2: I_variations
  trafomatrix.Shape(6,numdof_ele2);

  beamele2->GetGeneralizedInterpolationMatrixVariationsAtXi(
      trafomatrix,
      beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)));

  if (ele2force != NULL)
  {
    ele2force->Size(numdof_ele2);
    ele2force->Multiply('T','N',1.0,trafomatrix,bspotforce2,0.0);
  }

  if (ele21stiff != NULL)
  {
    ele21stiff->Shape(numdof_ele2,numdof_ele1);
    ele21stiff->Multiply('T','N',1.0,trafomatrix,auxmat21,0.0);
  }

  if (ele22stiff != NULL)
  {
    auxmat22.Shape(numdof_ele2,6);
    auxmat22.Multiply('T','N',1.0,trafomatrix,bspotstiff22,0.0);
  }

  // element 2: I_increments
  if (ele11stiff != NULL or ele21stiff != NULL)
  {
    trafomatrix.Shape(6,numdof_ele2);

    beamele2->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomatrix,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)),
        ele2disp);

    if (ele12stiff != NULL)
    {
      ele12stiff->Shape(numdof_ele1,numdof_ele2);
      ele12stiff->Multiply('N','N',1.0,auxmat12,trafomatrix,0.0);
    }

    if (ele22stiff != NULL)
    {
      ele22stiff->Shape(numdof_ele2,numdof_ele2);
      ele22stiff->Multiply('N','N',1.0,auxmat22,trafomatrix,0.0);
    }
  }

}

}
}

#endif
