/*-----------------------------------------------------------*/
/*!
\file biopolynet_calc_utils.cpp

\brief utils for biopolymer network business

\maintainer Jonas Eichinger, Maximilian Grill

\level 3

*/
/*-----------------------------------------------------------*/

#include "../drt_beaminteraction/biopolynet_calc_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_beaminteraction/beam3contact_utils.H"
#include "../drt_beaminteraction/beam_to_beam_linkage.H"
#include <Epetra_FEVector.h>
#include "../drt_beaminteraction/periodic_boundingbox.H"

#include "../drt_beam3/beam3_base.H"
#include "../drt_rigidsphere/rigidsphere.H"

namespace BIOPOLYNET
{
namespace UTILS
{

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void PeriodicBoundaryConsistentDisVector(
    Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    Teuchos::RCP<DRT::Discretization> const& discret)
{
  for ( int i = 0; i < discret->NumMyRowNodes(); ++i )
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

    for ( int dim = 0; dim < 3; ++dim )
    {
      int doflid = dis->Map().LID(dofnode[dim]);
      pbb->Shift1D( dim, (*dis)[doflid], node->X()[dim] );
    }
  }
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void UpdateCellsPositionRandomly(
    Teuchos::RCP<Epetra_Vector> dis,
    Teuchos::RCP<DRT::Discretization> const& discret,
    int const step)
{
//  DRT::Problem::Instance()->Random()->SetRandRange( 0.0, 1.0 );
//  const int updateevery = 200;
//
//  //todo: this is of course not nice, this needs to be done somewhere else
//  for( int i = 0; i < discret->ElementRowMap()->NumMyElements(); ++i )
//  {
//    static std::vector< std::vector< double > > Xnew(discret->ElementRowMap()->NumMyElements(), std::vector< double >(3) );
//    Xnew.resize(discret->ElementRowMap()->NumMyElements());
//
//    DRT::Element* eleptr = discret->gElement( discret->ElementRowMap()->GID(i) );
//
//    if ( dynamic_cast< const DRT::ELEMENTS::Rigidsphere* >(eleptr) != NULL and step >= updateevery )
//    {
//
//      std::vector<int> dofnode  = discret->Dof(eleptr->Nodes()[0]);
//
//      // random position of cell inside box
//      if( step % updateevery == 0 )
//      {
//        for ( int dim = 0; dim < 2; ++dim )
//        {
//          double edgelength = 5.6;
//          double min = 0.4;
//          Xnew[i][dim] = min + ( edgelength * DRT::Problem::Instance()->Random()->Uni() );
//        }
//      }
//
//      // loop over all dofs
//      for( int dim = 0; dim < 2; ++dim )
//      {
//        int doflid = dis->Map().LID(dofnode[dim]);
//        (*dis )[doflid] = Xnew[i][dim] - eleptr->Nodes()[0]->X()[dim];
//      }
//    }
//  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::vector<int> Permutation( int const number )
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
void GetCurrentElementDis(
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
void GetCurrentUnshiftedElementDis(
    DRT::Discretization const& discret,
    DRT::Element const* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    std::vector<double>& eledisp)
{
  GetCurrentElementDis(discret, ele, ia_discolnp, eledisp);

  // cast to beambase element
  DRT::ELEMENTS::Beam3Base const* beamele =
      dynamic_cast<DRT::ELEMENTS::Beam3Base const*>(ele);

#ifdef DEBUG
  if( beamele == NULL )
    dserror("Dynamic cast to beam3base failed");
#endif

  beamele->UnShiftNodePosition(eledisp, pbb);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GetPosAndTriadOfBindingSpot(
    DRT::Element* ele,
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
  if( beamele == NULL )
    dserror("Dynamic cast to beam3base failed");
#endif

  // get current position at binding spot xi
  beamele->GetPosOfBindingSpot( bspotpos, eledisp, locbspotnum, pbb );

  // get current triad at binding spot xi
  beamele->GetTriadOfBindingSpot( bspottriad, eledisp, locbspotnum );

}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void GetPosAndTriadOfBindingSpot(
    DRT::Discretization const& discret,
    DRT::Element* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const& pbb,
    int const locbspotnum,
    LINALG::Matrix<3,1>& bspotpos,
    LINALG::Matrix<3,3>& bspottriad)
{
  std::vector<double> eledisp;
  GetCurrentUnshiftedElementDis(discret, ele, ia_discolnp, pbb, eledisp);

  GetPosAndTriadOfBindingSpot(ele, pbb, locbspotnum,
      bspotpos, bspottriad, eledisp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void FEAssembleEleForceStiffIntoSystemVectorMatrix(
    const DRT::Discretization&         discret,
    std::vector<int> const&            elegid,
    std::vector< LINALG::SerialDenseVector > const& elevec,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& elemat,
    Teuchos::RCP<Epetra_FEVector>      fe_sysvec,
    Teuchos::RCP<LINALG::SparseMatrix> fe_sysmat)
{
  // the entries of elevecX  belong to the Dofs of the element with GID elegidX
  // the rows    of elematXY belong to the Dofs of the element with GID elegidX
  // the columns of elematXY belong to the Dofs of the element with GID elegidY
  const DRT::Element* ele1 = discret.gElement(elegid[0]);
  const DRT::Element* ele2 = discret.gElement(elegid[1]);

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
    fe_sysvec->SumIntoGlobalValues(elevec[0].Length(),&lmrow1[0],elevec[0].Values());
    fe_sysvec->SumIntoGlobalValues(elevec[1].Length(),&lmrow2[0],elevec[1].Values());
  }

  // and finally also assemble stiffness contributions
  if( fe_sysmat != Teuchos::null )
  {
    fe_sysmat->FEAssemble( elemat[0][0], lmrow1, lmrow1 );
    fe_sysmat->FEAssemble( elemat[0][1], lmrow1, lmrow2 );
    fe_sysmat->FEAssemble( elemat[1][0], lmrow2, lmrow1 );
    fe_sysmat->FEAssemble( elemat[1][1], lmrow2, lmrow2 );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void AssembleCenterlineDofForceStiffIntoElementForceStiff(
    DRT::Discretization const&                                     discret,
    std::vector<int> const&                                        elegid,
    std::vector< LINALG::SerialDenseVector > const&                eleforce_centerlineDOFs,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& elestiff_centerlineDOFs,
    std::vector< LINALG::SerialDenseVector >*                      eleforce,
    std::vector< std::vector< LINALG::SerialDenseMatrix > >*       elestiff)
{
  std::vector<unsigned int> numdof_ele(2);
  std::vector< std::vector<unsigned int> > ele_centerlinedofindices(2);

  for (unsigned int iele=0; iele<2; ++iele)
  {
    DRT::Element* ele = discret.gElement(elegid[iele]);

    // Todo implement method in DRT::Element or find alternative way of doing this
    // find out the elements' number of Dofs (=dimension of element vector/matrices)
    std::vector<int> lmrow;
    std::vector<int> dummy1, dummy2;

    ele->LocationVector(discret,lmrow,dummy1,dummy2);
    numdof_ele[iele] = lmrow.size();

    DRT::ELEMENTS::Beam3Base* beamele =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele);

    if( beamele != NULL)
    {
      beamele->CenterlineDofIndicesOfElement(ele_centerlinedofindices[iele]);
    }
    else
    {
      ele_centerlinedofindices[iele].resize( numdof_ele[iele] );
      for(unsigned int i = 0; i < numdof_ele[iele]; ++i )
        ele_centerlinedofindices[iele][i] = i;
    }
  }


  // assemble centerline DOF values correctly into element DOFvec vectors/matrices
  if (eleforce != NULL)
  {
    for (unsigned int iele=0; iele<2; ++iele)
    {
      // resize and clear variable
      ((*eleforce)[iele]).Size(numdof_ele[iele]);

      // safety check: dimensions
      if ((unsigned int) eleforce_centerlineDOFs[iele].RowDim() != ele_centerlinedofindices[iele].size())
        dserror("size mismatch! need to assemble %d values of centerline-Dof based "
            "force vector into element vector but only got %d element-local Dof indices",
            eleforce_centerlineDOFs[iele].RowDim(), ele_centerlinedofindices[iele].size());

      // Todo maybe use a more general 'SerialDenseAssemble' method here
      for (unsigned int idof=0; idof<ele_centerlinedofindices[iele].size(); ++idof)
        ((*eleforce)[iele])(ele_centerlinedofindices[iele][idof]) = eleforce_centerlineDOFs[iele](idof);
    }
  }

  if (elestiff != NULL)
  {
    for (unsigned int iele=0; iele<2; ++iele)
    {
      for (unsigned int jele=0; jele<2; ++jele)
      {
        // resize and clear variable
        ((*elestiff)[iele][jele]).Shape(numdof_ele[iele],numdof_ele[jele]);

        // safety check: dimensions
        if ((unsigned int) elestiff_centerlineDOFs[iele][jele].RowDim() != ele_centerlinedofindices[iele].size())
          dserror("size mismatch! need to assemble %d row values of centerline-Dof based "
              "stiffness matrix into element matrix but only got %d element-local Dof indices",
              elestiff_centerlineDOFs[iele][jele].RowDim(), ele_centerlinedofindices[iele].size());

        if ((unsigned int) elestiff_centerlineDOFs[iele][jele].ColDim() != ele_centerlinedofindices[jele].size())
          dserror("size mismatch! need to assemble %d column values of centerline-Dof based "
              "stiffness matrix into element matrix but only got %d element-local Dof indices",
              elestiff_centerlineDOFs[iele][jele].ColDim(), ele_centerlinedofindices[jele].size());

        for (unsigned int idof=0; idof<ele_centerlinedofindices[iele].size(); ++idof)
          for (unsigned int jdof=0; jdof<ele_centerlinedofindices[jele].size(); ++jdof)
            ((*elestiff)[iele][jele])(ele_centerlinedofindices[iele][idof],ele_centerlinedofindices[jele][jdof]) =
                elestiff_centerlineDOFs[iele][jele](idof,jdof);
      }
    }
  }

}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ExtractPosDofVecAbsoluteValues(
    DRT::Discretization const& discret,
    DRT::Element const* ele,
    Teuchos::RCP<Epetra_Vector> const& ia_discolnp,
    std::vector<double>& element_posdofvec_absolutevalues)
{
  std::vector<double> eledispvec;

  // extract the Dof values of this element from displacement vector
  GetCurrentElementDis( discret, ele, ia_discolnp, eledispvec );

   DRT::ELEMENTS::Beam3Base const* beam_element_ptr =
      dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

   if( beam_element_ptr != NULL )
   {
    // get the current absolute values for those Dofs relevant for centerline interpolation
    // initial values are added by element itself
    beam_element_ptr->ExtractCenterlineDofValuesFromElementStateVector(
        eledispvec,
        element_posdofvec_absolutevalues,
        true);
   }
   else
   {
     element_posdofvec_absolutevalues = eledispvec;
     for ( unsigned int dim = 0; dim < 3; ++dim )
       for ( int node = 0; node < ele->NumNode(); ++node )
         element_posdofvec_absolutevalues[ 3*node + dim] += ele->Nodes()[node]->X()[dim];
   }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBindingSpotForceToParentElements(
    DRT::Discretization const&                             discret,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const&        pbb,
    const Teuchos::RCP<Epetra_Vector>                      disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage> elepairptr,
    std::vector< LINALG::SerialDenseVector > const&        bspotforce,
    std::vector< LINALG::SerialDenseVector >&              eleforce)
{
  // auxiliary transformation matrix, will be resized and reused
  LINALG::SerialDenseMatrix trafomatrix;

  for( int elei = 0; elei < 2; ++elei )
  {
    DRT::Element* ele = discret.gElement( elepairptr->GetEleGid(elei) );

    DRT::ELEMENTS::Beam3Base* beamele =
        dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele);

    // get current element displacements
    std::vector<double> eledisp;
    GetCurrentUnshiftedElementDis( discret, ele, disp_np_col, pbb, eledisp);
    const int numdof_ele = eledisp.size();

    // zero out and set correct size of transformation matrix
    trafomatrix.Shape( 6, numdof_ele );

    // I_variations
    beamele->GetGeneralizedInterpolationMatrixVariationsAtXi( trafomatrix,
        beamele->GetBindingSpotXi( elepairptr->GetLocBSpotNum(elei) ), eledisp );

    eleforce[elei].Size(numdof_ele);
    eleforce[elei].Multiply('T','N',1.0,trafomatrix,bspotforce[elei],0.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBindingSpotStiffToParentElements(
    DRT::Discretization const&                                     discret,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const&                pbb,
    const Teuchos::RCP<Epetra_Vector>                              disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage>         elepairptr,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& bspotstiff,
    std::vector< std::vector< LINALG::SerialDenseMatrix > >&       elestiff)
{
  // Todo grill 02/17
  dserror("we can not evaluate the tangent stiffness matrix without evaluating "
      "the residual in case of Kirchhoff beam Beam3k. This is because we have a coupled "
      "term and the residual vector is required to compute the full linearization. "
      "Use the combined evaluation method instead!");

  // todo: put this in a loop
  DRT::Element* ele1 = discret.gElement(elepairptr->GetEleGid(0));
  DRT::Element* ele2 = discret.gElement(elepairptr->GetEleGid(1));

  DRT::ELEMENTS::Beam3Base* beamele1 =
      dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele1);
  DRT::ELEMENTS::Beam3Base* beamele2 =
      dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele2);

  // get current element displacements
  std::vector<double> ele1disp;
  GetCurrentUnshiftedElementDis(discret, ele1,disp_np_col,pbb,ele1disp);
  const int numdof_ele1 = ele1disp.size();

  std::vector<double> ele2disp;
  GetCurrentUnshiftedElementDis(discret,ele2,disp_np_col,pbb,ele2disp);
  const int numdof_ele2 = ele2disp.size();

  // transformation matrix, will be resized and reused for various needs
  LINALG::SerialDenseMatrix trafomatrix;

  // auxiliary matrices required to store intermediate results after first of two
  // consecutive matrix-matrix products for each of four stiffness matrices
  std::vector< std::vector< LINALG::SerialDenseMatrix > > auxmat( 2,
      std::vector< LINALG::SerialDenseMatrix >(2) );

  // element 1:
  {
    // zero out and set correct size of transformation matrix
     trafomatrix.Shape(6,numdof_ele1);

    // i) I_variations
    beamele1->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomatrix,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)), ele1disp);

    auxmat[0][0].Shape(numdof_ele1,6);
    auxmat[0][0].Multiply('T','N',1.0,trafomatrix,bspotstiff[0][0],0.0);

    auxmat[0][1].Shape(numdof_ele1,6);
    auxmat[0][1].Multiply('T','N',1.0,trafomatrix,bspotstiff[0][1],0.0);

    // ii) I_increments
    trafomatrix.Shape(6,numdof_ele1);

    beamele1->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomatrix,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)),
        ele1disp);

    elestiff[0][0].Shape( numdof_ele1, numdof_ele1 );
    elestiff[0][0].Multiply('N','N',1.0,auxmat[0][0],trafomatrix,0.0);

    auxmat[1][0].Shape(6,numdof_ele1);
    auxmat[1][0].Multiply('N','N',1.0,bspotstiff[1][0],trafomatrix,0.0);
  }

  // element 2:
  {
    // i) I_variations
    trafomatrix.Shape(6,numdof_ele2);

    beamele2->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomatrix,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)), ele2disp);

    elestiff[1][0].Shape(numdof_ele2,numdof_ele1);
    elestiff[1][0].Multiply('T','N',1.0,trafomatrix,auxmat[1][0],0.0);

    auxmat[1][1].Shape(numdof_ele2,6);
    auxmat[1][1].Multiply('T','N',1.0,trafomatrix,bspotstiff[1][1],0.0);

  // ii) I_increments
    trafomatrix.Shape(6,numdof_ele2);

    beamele2->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomatrix,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)),
        ele2disp);

    elestiff[0][1].Shape(numdof_ele1,numdof_ele2);
    elestiff[0][1].Multiply('N','N',1.0,auxmat[1][0],trafomatrix,0.0);

    elestiff[1][1].Shape(numdof_ele2,numdof_ele2);
    elestiff[1][1].Multiply('N','N',1.0,auxmat[1][1],trafomatrix,0.0);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ApplyBindingSpotForceStiffToParentElements(
    DRT::Discretization const&                                     discret,
    Teuchos::RCP<GEO::MESHFREE::BoundingBox> const&                pbb,
    const Teuchos::RCP<Epetra_Vector>                              disp_np_col,
    const Teuchos::RCP<BEAMINTERACTION::BeamToBeamLinkage>         elepairptr,
    std::vector< LINALG::SerialDenseVector > const&                bspotforce,
    std::vector< std::vector< LINALG::SerialDenseMatrix > > const& bspotstiff,
    std::vector< LINALG::SerialDenseVector >&                      eleforce,
    std::vector< std::vector< LINALG::SerialDenseMatrix > >&       elestiff)
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

  // todo @grill: put this in a loop
  DRT::Element* ele1 = discret.gElement(elepairptr->GetEleGid(0));
  DRT::Element* ele2 = discret.gElement(elepairptr->GetEleGid(1));

  DRT::ELEMENTS::Beam3Base* beamele1 =
      dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele1);
  DRT::ELEMENTS::Beam3Base* beamele2 =
      dynamic_cast<DRT::ELEMENTS::Beam3Base*>(ele2);

  // get current element displacements
  std::vector<double> ele1disp;
  GetCurrentUnshiftedElementDis(discret, ele1,disp_np_col,pbb,ele1disp);
  const int numdof_ele1 = ele1disp.size();

  std::vector<double> ele2disp;
  GetCurrentUnshiftedElementDis(discret,ele2,disp_np_col,pbb,ele2disp);
  const int numdof_ele2 = ele2disp.size();

  // transformation matrix, will be resized and reused for various needs
  LINALG::SerialDenseMatrix trafomatrix;

  // auxiliary matrices required to store intermediate results after first of two
  // consecutive matrix-matrix products for each of four stiffness matrices
  std::vector< std::vector< LINALG::SerialDenseMatrix > > auxmat( 2,
      std::vector< LINALG::SerialDenseMatrix >(2) );

  // contribution to stiffmat from linearization of generalized interpolation matrix for variations
  LINALG::SerialDenseMatrix stiffmat_lin_Ivar;

  // zero out and set correct size of transformation matrix
  trafomatrix.Shape(6,numdof_ele1);

  // todo:
  // element 1:
  {
    // i) I_variations
    beamele1->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomatrix,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)),
        ele1disp);

    eleforce[0].Size(numdof_ele1);
    eleforce[0].Multiply('T','N',1.0,trafomatrix,bspotforce[0],0.0);

    auxmat[0][0].Shape(numdof_ele1,6);
    auxmat[0][0].Multiply('T','N',1.0,trafomatrix,bspotstiff[0][0],0.0);

    auxmat[0][1].Shape(numdof_ele1,6);
    auxmat[0][1].Multiply('T','N',1.0,trafomatrix,bspotstiff[0][1],0.0);

    // ii) I_increments
    trafomatrix.Shape(6,numdof_ele1);

    beamele1->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomatrix,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)),
        ele1disp);

    elestiff[0][0].Shape(numdof_ele1,numdof_ele1);
    elestiff[0][0].Multiply('N','N',1.0,auxmat[0][0],trafomatrix,0.0);

    auxmat[1][0].Shape(6,numdof_ele1);
    auxmat[1][0].Multiply('N','N',1.0,bspotstiff[1][0],trafomatrix,0.0);


    // additional contribution from linearization of generalized interpolation matrix for variations
    stiffmat_lin_Ivar.Shape(numdof_ele1,numdof_ele1);

    beamele1->GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(
        stiffmat_lin_Ivar,
        beamele1->GetBindingSpotXi(elepairptr->GetLocBSpotNum(0)),
        ele1disp,
        bspotforce[0]);

    elestiff[0][0].Update(1.0, stiffmat_lin_Ivar, 1.0);

  }

  // element 2
  {
    // i) I_variations
    trafomatrix.Shape(6,numdof_ele2);

    beamele2->GetGeneralizedInterpolationMatrixVariationsAtXi(
        trafomatrix,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)),
        ele2disp);

    eleforce[1].Size(numdof_ele2);
    eleforce[1].Multiply('T','N',1.0,trafomatrix,bspotforce[1],0.0);

    elestiff[1][0].Shape(numdof_ele2,numdof_ele1);
    elestiff[1][0].Multiply('T','N',1.0,trafomatrix,auxmat[1][0],0.0);

    auxmat[1][1].Shape(numdof_ele2,6);
    auxmat[1][1].Multiply('T','N',1.0,trafomatrix,bspotstiff[1][1],0.0);

    // ii) I_increments
    trafomatrix.Shape(6,numdof_ele2);

    beamele2->GetGeneralizedInterpolationMatrixIncrementsAtXi(
        trafomatrix,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)),
        ele2disp);

    elestiff[0][1].Shape(numdof_ele1,numdof_ele2);
    elestiff[0][1].Multiply('N','N',1.0,auxmat[0][1],trafomatrix,0.0);

    elestiff[1][1].Shape(numdof_ele2,numdof_ele2);
    elestiff[1][1].Multiply('N','N',1.0,auxmat[1][1],trafomatrix,0.0);


    // additional contribution from linearization of generalized interpolation matrix for variations
    stiffmat_lin_Ivar.Shape(numdof_ele2,numdof_ele2);

    beamele2->GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(
        stiffmat_lin_Ivar,
        beamele2->GetBindingSpotXi(elepairptr->GetLocBSpotNum(1)),
        ele2disp,
        bspotforce[1]);

    elestiff[1][1].Update(1.0, stiffmat_lin_Ivar, 1.0);

  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void SetupEleTypeMapExtractor(
    Teuchos::RCP<DRT::Discretization> const& discret,
    Teuchos::RCP<LINALG::MultiMapExtractor>& eletypeextractor)
{
  std::vector< std::set<int> > eletypeset(2);

  for ( int i = 0; i < discret->NumMyRowElements(); ++i )
  {
    // get ele pointer
    DRT::Element* eleptr = discret->lRowElement(i);

    if( dynamic_cast< const DRT::ELEMENTS::Beam3Base* >(eleptr) != NULL )
    {
      eletypeset[0].insert( eleptr->Id() );
    }
    else if ( dynamic_cast< const DRT::ELEMENTS::Rigidsphere* >(eleptr) != NULL)
    {
      eletypeset[1].insert( eleptr->Id() );
    }
    else
    {
      dserror("eletype multi map extractor cannot yet handle current element type.");
    }
  }

  std::vector< Teuchos::RCP< const Epetra_Map > > maps( eletypeset.size() );
  for( int i = 0; i < static_cast<int>( eletypeset.size() ); ++i )
  {
    std::vector<int> mapvec( eletypeset[i].begin(), eletypeset[i].end() );
    eletypeset[0].clear();
    maps[i] = Teuchos::rcp( new Epetra_Map(-1, mapvec.size(),
        &mapvec[0], 0, discret->Comm()));
  }

  eletypeextractor->Setup( *discret()->ElementRowMap(), maps );
}

}
}
