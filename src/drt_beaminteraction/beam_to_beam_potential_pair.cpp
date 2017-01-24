/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_beam_potential_pair.cpp

\brief One beam-to-beam potential-based interacting pair (two beam elements)

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#include "beam_to_beam_potential_pair.H"

#include "beam_potential_params.H"

// Todo get rid of outdated header inclusions
#include "beam3contact_utils.H"
#include "../drt_inpar/inpar_beampotential.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"
#include "../headers/FAD_utils.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include <Teuchos_RCP.hpp>

#include "Teuchos_TimeMonitor.hpp"
#include "beam_potential_params.H"
#include "beam3contact_defines.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::BeamToBeamPotentialPair():
    BeamPotentialPair(),
    beam_element1_(NULL),
    beam_element2_(NULL),
    k_(0.0),
    m_(0.0),
    ele1length_(0.0),
    ele2length_(0.0),
    radius1_(0.0),
    radius2_(0.0)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::Setup()
{
  CheckInit();

  // call setup of base class first
  BeamPotentialPair::Setup();


  ele1pos_.Clear();
  ele2pos_.Clear();

  fpot1_.Clear();
  fpot2_.Clear();
  stiffpot1_.Clear();
  stiffpot2_.Clear();


  // Calculate initial length of beam elements
  beam_element1_ = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element1());
  ele1length_ = beam_element1_->RefLength();
  beam_element2_ = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element2());
  ele2length_ = beam_element2_->RefLength();

  radius1_ = BEAMCONTACT::CalcEleRadius(Element1());
  radius2_ = BEAMCONTACT::CalcEleRadius(Element2());

  if (Element1()->ElementType() != Element2()->ElementType())
    dserror("The class BeamToBeamPotentialPair currently only supports element "
        "pairs of the same beam element type!");

  // Todo check this prerequisite, necessary?
  if (Element1()->Id() >= Element2()->Id())
  {
    dserror("Element 1 has to have the smaller element-ID. Ele1GID: %d, Ele2GID: %d. Adapt your contact search!",
        Element1()->Id(),Element2()->Id());
  }

  // initialize line charge conditions applied to element1 and element2
  linechargeconds_.clear();

  issetup_ = true;

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
bool BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::Evaluate(
    LINALG::SerialDenseVector* forcevec1,
    LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11,
    LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21,
    LINALG::SerialDenseMatrix* stiffmat22,
    const std::vector<DRT::Condition*> linechargeconds,
    const double k,
    const double m)
{
  // nothing to do in case of k==0.0
  if (k==0.0) return false;

   // reset fpot and stiffpot class variables
   fpot1_.Clear();
   fpot2_.Clear();
   stiffpot1_.Clear();
   stiffpot2_.Clear();

  const unsigned int dim1 = 3*numnodes*numnodalvalues;
  const unsigned int dim2 = 3*numnodes*numnodalvalues;

  // set class variables
  if (linechargeconds.size() == 2)
  {
    for (unsigned int i=0; i<2; ++i)
    {
      if (linechargeconds[i]->Type() == DRT::Condition::BeamPotential_LineChargeDensity)
        linechargeconds_.push_back(linechargeconds[i]);
      else
        dserror("Provided line charge condition is not of correct type BeamPotential_LineChargeDensity!");
    }
  }
  else
    dserror("Expected TWO dline charge conditions!");

  k_=k;
  m_=m;

  // prepare FAD
#ifdef AUTOMATICDIFF
  // The 2*3*numnodes*numnodalvalues primary DoFs are the components of the nodal positions / tangents.
  for (unsigned int i=0;i<3*numnodes*numnodalvalues;i++)
    ele1pos_(i).diff(i,2*3*numnodes*numnodalvalues);

  for (unsigned int i=0;i<3*numnodes*numnodalvalues;i++)
    ele2pos_(i).diff(3*numnodes*numnodalvalues+i,2*3*numnodes*numnodalvalues);
#endif

  // compute the values for fpot1, fpot2, stiffpot1, stiffpot2
  // TODO specify approximation type in input file: FullInt, LargeSepApprox, SmallSepApprox (Langbein)
  // should be chosen depending on potential law
  EvaluateFpotandStiffpot_LargeSepApprox();


  // resize variables and fill with pre-computed values
  if (forcevec1 != NULL)
  {
    forcevec1->Size(dim1);
    for (unsigned int i=0; i<dim1; ++i)
      (*forcevec1)(i) = FADUTILS::CastToDouble(fpot1_(i));
  }
  if (forcevec2 != NULL)
  {
    forcevec2->Size(dim2);
    for (unsigned int i=0; i<dim2; ++i)
      (*forcevec2)(i) = FADUTILS::CastToDouble(fpot2_(i));
  }

  if (stiffmat11 != NULL)
  {
    stiffmat11->Shape(dim1,dim1);
    for (unsigned int irow=0; irow<dim1; ++irow)
      for (unsigned int icol=0; icol<dim1; ++icol)
        (*stiffmat11)(irow,icol) =  FADUTILS::CastToDouble(stiffpot1_(irow,icol));
  }
  if (stiffmat12 != NULL)
  {
    stiffmat12->Shape(dim1,dim2);
    for (unsigned int irow=0; irow<dim1; ++irow)
      for (unsigned int icol=0; icol<dim2; ++icol)
        (*stiffmat12)(irow,icol) =  FADUTILS::CastToDouble(stiffpot1_(irow,dim1+icol));
  }
  if (stiffmat21 != NULL)
  {
    stiffmat21->Shape(dim2,dim1);
    for (unsigned int irow=0; irow<dim2; ++irow)
      for (unsigned int icol=0; icol<dim1; ++icol)
        (*stiffmat21)(irow,icol) =  FADUTILS::CastToDouble(stiffpot2_(irow,icol));
  }
  if (stiffmat22 != NULL)
  {
    stiffmat22->Shape(dim2,dim2);
    for (unsigned int irow=0; irow<dim2; ++irow)
      for (unsigned int icol=0; icol<dim2; ++icol)
        (*stiffmat22)(irow,icol) =  FADUTILS::CastToDouble(stiffpot2_(irow,dim1+icol));
  }

  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::EvaluateFpotandStiffpot_LargeSepApprox()
{
  // Set gauss integration rule
  DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule_line_10point;
//  DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule_line_32point;

  // Get gauss points (gp) for integration
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);
  // number of gps
  const int numgp=gausspoints.nquad;

  // vectors for shape functions and their derivatives
  // Attention: these are individual shape function values, NOT shape function matrices
  // values at all gauss points are stored in advance (more efficient due to double integral)
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N1_i(numgp);        // = N1_i
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N2_i(numgp);        // = N2_i
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N1_i_xi(numgp);     // = N1_i,xi
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N2_i_xi(numgp);     // = N2_i,eta

  // coords and derivatives of the two gauss points
  LINALG::TMatrix<TYPE, 3, 1> r1(true);                               // = r1
  LINALG::TMatrix<TYPE, 3, 1> r2(true);                               // = r2
  LINALG::TMatrix<TYPE, 3, 1> dist(true);                             // = r1-r2
  TYPE norm_dist= 0.0;                                                // = |r1-r2|

  // Evaluate shape functions at gauss points and store values
  GetShapeFunctions(N1_i,N2_i,N1_i_xi,N2_i_xi,gausspoints);

  // evaluate charge densities from DLINE charge condition specified in input file
  double q1 = linechargeconds_[0]->GetDouble("val");
  double q2 = linechargeconds_[1]->GetDouble("val");

  // TODO evaluate given functions in line charge conditions! for now: dserror
  if (linechargeconds_[0]->GetInt("funct") != -1 or linechargeconds_[1]->GetInt("funct") != -1)
    dserror("DLINE beam potential charge condition: No functions allowed yet! Set 'funct' to '-1' -> off");

  // auxiliary variable
  LINALG::TMatrix<TYPE, 3, 1> fpot_tmp(true);

  // determine prefactor of the integral (depends on whether surface or volume potential is applied)
  double prefactor = k_ * m_;

  switch ( Params()->PotentialType() )
  {
  case INPAR::BEAMPOTENTIAL::beampot_surf:
    prefactor *= 4 * radius1_ * radius2_ * std::pow(M_PI,2);
    break;
  case INPAR::BEAMPOTENTIAL::beampot_vol:
    prefactor *= std::pow(radius1_,2) * std::pow(radius2_,2) * std::pow(M_PI,2);
    break;
  default:
    dserror("No valid BEAMPOTENTIAL_TYPE specified. Choose either Surface or Volume in input file!");
  }

  // loop over gauss points on ele1
  for (int gp1=0; gp1 < numgp; ++gp1)
  {
    ComputeCoords(r1,N1_i[gp1],ele1pos_);

    double jacobifac1 = BeamElement1()->GetJacobiFacAtXi(gausspoints.qxg[gp1][0]);

    // loop over gauss points on ele2
    for (int gp2=0; gp2 < numgp; ++gp2)
    {

      ComputeCoords(r2,N2_i[gp2],ele2pos_);

      dist = FADUTILS::DiffVector(r1,r2);

      norm_dist = FADUTILS::VectorNorm<3>(dist);

      // auxiliary variables to store pre-calculated common terms
      TYPE norm_dist_exp1 = 0.0;
      if(norm_dist !=0.0)
      {
        norm_dist_exp1 = std::pow(norm_dist,-m_-2);
      }
      else
      {
        dserror("\n|r1-r2|=0 ! Interacting points are identical! Potential law not defined in this case!"
            " Think about shifting nodes in unconverged state?!");
      }

      double q1q2_JacFac_GaussWeights =
          q1 * q2 * jacobifac1 * BeamElement2()->GetJacobiFacAtXi(gausspoints.qxg[gp2][0]) * gausspoints.qwgt[gp1] * gausspoints.qwgt[gp2];

      // compute fpot_tmp here, same for both element forces
      for (unsigned int i=0; i<3; ++i)
        fpot_tmp(i) = q1q2_JacFac_GaussWeights * norm_dist_exp1 * dist(i);

      //********************************************************************
      // calculate fpot1: force on element 1
      //********************************************************************
      // sum up the contributions of all nodes (in all dimensions)
      for (unsigned int i=0; i<(numnodes*numnodalvalues); ++i)
      {
        // loop over dimensions
        for (unsigned int j=0; j<3; ++j)
        {
          fpot1_(3*i+j) -= N1_i[gp1](i)*fpot_tmp(j);
        }
      }

      //********************************************************************
      // calculate fpot2: force on element 2
      //********************************************************************
      // sum up the contributions of all nodes (in all dimensions)
      for (unsigned int i=0; i<(numnodes*numnodalvalues); ++i)
      {
        // loop over dimensions
        for (unsigned int j=0; j<3; ++j)
        {
          fpot2_(3*i+j) += N2_i[gp2](i)*fpot_tmp(j);
        }
      }

      //********************************************************************
      // calculate stiffpot1
      //********************************************************************
      // auxiliary variables (same for both elements)
      TYPE norm_dist_exp2 = (m_+2) * std::pow(norm_dist,-m_-4);

      LINALG::TMatrix<TYPE, 3, 3> dist_dist_T(true);

      for (unsigned int i=0; i<3; ++i)
      {
        for (unsigned int j=0; j<=i; ++j)
        {
          dist_dist_T(i,j) = dist(i) * dist(j);
          if(i!=j) dist_dist_T(j,i) = dist_dist_T(i,j);
        }
      }

      for (unsigned int i=0; i<(numnodes*numnodalvalues); ++i)
      {

        // d (Res_1) / d (d_1)
        for (unsigned int j=0; j<(numnodes*numnodalvalues); ++j)
        {

          for (unsigned int idim=0; idim<3; ++idim)
          {
            stiffpot1_(3*i+idim,3*j+idim) -=
                norm_dist_exp1 * N1_i[gp1](i)*N1_i[gp1](j) *  q1q2_JacFac_GaussWeights;

            for (unsigned int jdim=0; jdim<3; ++jdim)
            {
              stiffpot1_(3*i+idim,3*j+jdim) +=
                  norm_dist_exp2 * N1_i[gp1](i) * dist_dist_T(idim,jdim) * N1_i[gp1](j) *  q1q2_JacFac_GaussWeights;
            }
          }
        }

        // d (Res_1) / d (d_2)
        for (unsigned int j=0; j<(numnodes*numnodalvalues); ++j)
        {

          for (unsigned int idim=0; idim<3; ++idim)
          {
            stiffpot1_(3*i+idim,3*(numnodes*numnodalvalues+j)+idim) +=
                norm_dist_exp1 * N1_i[gp1](i)*N2_i[gp2](j) *  q1q2_JacFac_GaussWeights;

            for (unsigned int jdim=0; jdim<3; ++jdim)
            {
              stiffpot1_(3*i+idim,3*(numnodes*numnodalvalues+j)+jdim) -=
                  norm_dist_exp2 * N1_i[gp1](i) * dist_dist_T(idim,jdim) * N2_i[gp2](j) *  q1q2_JacFac_GaussWeights;
            }
          }
        }

      }

      //********************************************************************
      // calculate stiffpot2
      //********************************************************************
      for (unsigned int i=0; i<(numnodes*numnodalvalues); ++i)
      {

        // d (Res_2) / d (d_1)
        for (unsigned int j=0; j<(numnodes*numnodalvalues); ++j)
        {

          for (unsigned int idim=0; idim<3; ++idim)
          {
            stiffpot2_(3*i+idim,3*j+idim) +=
                norm_dist_exp1 * N2_i[gp2](i)*N1_i[gp1](j) *  q1q2_JacFac_GaussWeights;

            for (unsigned int jdim=0; jdim<3; ++jdim)
            {
              stiffpot2_(3*i+idim,3*j+jdim) -=
                  norm_dist_exp2 * N2_i[gp2](i) * dist_dist_T(idim,jdim) * N1_i[gp1](j) *  q1q2_JacFac_GaussWeights;
            }
          }
        }

        // d (Res_2) / d (d_2)
        for (unsigned int j=0; j<(numnodes*numnodalvalues); ++j)
        {

          for (unsigned int idim=0; idim<3; ++idim)
          {
            stiffpot2_(3*i+idim,3*(numnodes*numnodalvalues+j)+idim) -=
                norm_dist_exp1 * N2_i[gp2](i)*N2_i[gp2](j) *  q1q2_JacFac_GaussWeights;

            for (unsigned int jdim=0; jdim<3; ++jdim)
            {
              stiffpot2_(3*i+idim,3*(numnodes*numnodalvalues+j)+jdim) +=
                  norm_dist_exp2 * N2_i[gp2](i) * dist_dist_T(idim,jdim) * N2_i[gp2](j) *  q1q2_JacFac_GaussWeights;
            }
          }
        }

      }

    } // end gauss quadrature loop (element 2)
  } // end gauss quadrature loop (element 1)


  // apply constant prefactor
  for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i) fpot1_(i)*=prefactor;
  for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i) fpot2_(i)*=prefactor;

  stiffpot1_.Scale(prefactor);
  stiffpot2_.Scale(prefactor);

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::Print(std::ostream& out) const
{
  CheckInitSetup();

  out << "\nInstance of BeamToBeamPotentialPair (EleGIDs " << Element1()->Id()
      << " & " << Element2()->Id() << "):";
  out << "\nele1 dofvec: " << ele1pos_;
  out << "\nele2 dofvec: " << ele2pos_;

  out << "\n";
  // Todo add more relevant information here
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::
PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

  // Todo difficulty here is that the same element pair is evaluated more than once
  //      to be more precise, once for every common potlaw;
  //      contribution of previous evaluations is overwritten if multiple potlaws are applied

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::GetShapeFunctions( std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i,
                                                                                 std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N2_i,
                                                                                 std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i_xi,
                                                                                 std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N2_i_xi,
                                                                                 DRT::UTILS::IntegrationPoints1D& gausspoints)
{
  // get both discretization types
  const DRT::Element::DiscretizationType distype1 = BeamElement1()->Shape();
  const DRT::Element::DiscretizationType distype2 = BeamElement2()->Shape();

  if (numnodalvalues==1)
  {
    for (int gp=0; gp<gausspoints.nquad; ++gp)
    {
      // get values and derivatives of shape functions
      DRT::UTILS::shape_function_1D(N1_i[gp], gausspoints.qxg[gp][0], distype1);
      DRT::UTILS::shape_function_1D(N2_i[gp], gausspoints.qxg[gp][0], distype2);
      DRT::UTILS::shape_function_1D_deriv1(N1_i_xi[gp], gausspoints.qxg[gp][0], distype1);
      DRT::UTILS::shape_function_1D_deriv1(N2_i_xi[gp], gausspoints.qxg[gp][0], distype2);
    }
  }
  else if (numnodalvalues==2)
  {
    for (int gp=0; gp<gausspoints.nquad; ++gp)
    {
      // get values and derivatives of shape functions
      DRT::UTILS::shape_function_hermite_1D(N1_i[gp], gausspoints.qxg[gp][0], ele1length_, distype1);
      DRT::UTILS::shape_function_hermite_1D(N2_i[gp], gausspoints.qxg[gp][0], ele2length_, distype2);
      DRT::UTILS::shape_function_hermite_1D_deriv1(N1_i_xi[gp], gausspoints.qxg[gp][0], ele1length_, distype1);
      DRT::UTILS::shape_function_hermite_1D_deriv1(N2_i_xi[gp], gausspoints.qxg[gp][0], ele2length_, distype2);
    }
  }
  else
    dserror("Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) values are valid!");

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::ComputeCoords(LINALG::TMatrix<TYPE,3,1>& r,
                                                                            const LINALG::Matrix<1,numnodes*numnodalvalues>& N_i,
                                                                            const LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> elepos)
{
  r.Clear();

  // compute output variable
  for (unsigned int i=0;i<3;i++)
  {
    for (unsigned int j=0;j<numnodes*numnodalvalues;j++)
    {
      r(i)+=N_i(j)*elepos(3*j+i);
    }
  }

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues>::ResetState(
    const std::vector<double>& centerline_dofvec_ele1,
    const std::vector<double>& centerline_dofvec_ele2)
{
  for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i)
  {
      ele1pos_(i) = centerline_dofvec_ele1[i];
      ele2pos_(i) = centerline_dofvec_ele2[i];
  }

}


// explicit template instantiations
template class BEAMINTERACTION::BeamToBeamPotentialPair<2,1>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<3,1>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<4,1>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<5,1>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<2,2>;
