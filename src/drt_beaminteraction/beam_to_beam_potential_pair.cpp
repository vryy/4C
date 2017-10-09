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
#include "../drt_fem_general/drt_utils_integration.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"
#include "../headers/FAD_utils.H"

#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include <Teuchos_RCP.hpp>
#include <Sacado.hpp>

#include "Teuchos_TimeMonitor.hpp"
#include "beam_potential_params.H"
#include "beam3contact_defines.H"

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::BeamToBeamPotentialPair():
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
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::Setup()
{
  CheckInit();

  // call setup of base class first
  BeamPotentialPair::Setup();


  ele1pos_.Clear();
  ele2pos_.Clear();


  // get initial length of beam elements
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
    dserror("Element 1 has to have the smaller element-ID. Ele1GID: %d, Ele2GID: %d."
        " Adapt your search algorithm!", Element1()->Id(),Element2()->Id());
  }

  // initialize line charge conditions applied to element1 and element2
  linechargeconds_.clear();

  issetup_ = true;

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
bool BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::Evaluate(
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


  // set class variables
  if (linechargeconds.size() == 2)
  {
    for (unsigned int i=0; i<2; ++i)
    {
      if (linechargeconds[i]->Type() == DRT::Condition::BeamPotential_LineChargeDensity)
        linechargeconds_.push_back(linechargeconds[i]);
      else
        dserror("Provided line charge condition is not of correct type"
            "BeamPotential_LineChargeDensity!");
    }
  }
  else
    dserror("Expected TWO dline charge conditions!");

  k_=k;
  m_=m;



  const unsigned int dim1 = 3*numnodes*numnodalvalues;
  const unsigned int dim2 = 3*numnodes*numnodalvalues;

  LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1> force_pot1(true);
  LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1> force_pot2(true);


  if (stiffmat11 != NULL)
    stiffmat11->Shape(dim1,dim1);
  if (stiffmat12 != NULL)
    stiffmat12->Shape(dim1,dim2);
  if (stiffmat21 != NULL)
    stiffmat21->Shape(dim2,dim1);
  if (stiffmat22 != NULL)
    stiffmat22->Shape(dim2,dim2);


  // prepare differentiation via FAD if desired
  SetAutomaticDifferentiationVariablesIfRequired( ele1pos_, ele2pos_ );


  // compute the values for element residual vectors ('force') and linearizations ('stiff')
  switch ( Params()->Strategy() )
  {
    case INPAR::BEAMPOTENTIAL::strategy_doublelengthspec_largesepapprox:
    {
      EvaluateFpotandStiffpot_LargeSepApprox( force_pot1, force_pot2,
          stiffmat11, stiffmat12, stiffmat21, stiffmat22);
      break;
    }

    case INPAR::BEAMPOTENTIAL::strategy_doublelengthspec_smallsepapprox:
    {
      EvaluateFpotandStiffpot_DoubleLengthSpecific_SmallSepApprox( force_pot1, force_pot2,
          stiffmat11, stiffmat12, stiffmat21, stiffmat22);
      break;
    }

    default:
      dserror("Invalid strategy to evaluate beam interaction potential!");
  }


  // resize variables and fill with pre-computed values
  if (forcevec1 != NULL)
  {
    forcevec1->Size(dim1);
    for (unsigned int i=0; i<dim1; ++i)
      (*forcevec1)(i) = FADUTILS::CastToDouble(force_pot1(i));
  }
  if (forcevec2 != NULL)
  {
    forcevec2->Size(dim2);
    for (unsigned int i=0; i<dim2; ++i)
      (*forcevec2)(i) = FADUTILS::CastToDouble(force_pot2(i));
  }

  if ( stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL )
  {
    CalcStiffmatAutomaticDifferentiationIfRequired(
        force_pot1,
        force_pot2,
        *stiffmat11,
        *stiffmat12,
        *stiffmat21,
        *stiffmat22);


  }

  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
EvaluateFpotandStiffpot_LargeSepApprox(
    LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1>& force_pot1,
    LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1>& force_pot2,
    LINALG::SerialDenseMatrix* stiffmat11,
    LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21,
    LINALG::SerialDenseMatrix* stiffmat22)
{
  // Set gauss integration rule
  DRT::UTILS::GaussRule1D gaussrule = GetGaussRule();

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
  LINALG::TMatrix<T, 3, 1> r1(true);                               // = r1
  LINALG::TMatrix<T, 3, 1> r2(true);                               // = r2
  LINALG::TMatrix<T, 3, 1> dist(true);                             // = r1-r2
  T norm_dist= 0.0;                                                // = |r1-r2|

  // Evaluate shape functions at gauss points and store values
  GetShapeFunctions(N1_i,N2_i,N1_i_xi,N2_i_xi,gausspoints);

  // evaluate charge densities from DLINE charge condition specified in input file
  double q1 = linechargeconds_[0]->GetDouble("val");
  double q2 = linechargeconds_[1]->GetDouble("val");

  // TODO evaluate given functions in line charge conditions! for now: dserror
  if (linechargeconds_[0]->GetInt("funct") != -1 or linechargeconds_[1]->GetInt("funct") != -1)
    dserror("DLINE beam potential charge condition: No functions allowed yet! Set 'funct' to '-1' -> off");

  // auxiliary variable
  LINALG::TMatrix<T, 3, 1> fpot_tmp(true);

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
      T norm_dist_exp1 = 0.0;
      if ( norm_dist!=0.0 )
      {
        norm_dist_exp1 = std::pow(norm_dist,-m_-2);
      }
      else
      {
        dserror("\n|r1-r2|=0 ! Interacting points are identical! Potential law not defined in this"
            " case! Think about shifting nodes in unconverged state?!");
      }

      double q1q2_JacFac_GaussWeights =
          q1 * q2 * jacobifac1 * BeamElement2()->GetJacobiFacAtXi(gausspoints.qxg[gp2][0])
          * gausspoints.qwgt[gp1] * gausspoints.qwgt[gp2];

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
          force_pot1(3*i+j) -= N1_i[gp1](i)*fpot_tmp(j);
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
          force_pot2(3*i+j) += N2_i[gp2](i)*fpot_tmp(j);
        }
      }

      // evaluate analytic contributions to linearization
      if ( stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL )
      {
        EvaluateStiffpotAnalyticContributions_LargeSepApprox(
            dist,
            norm_dist,
            norm_dist_exp1,
            q1q2_JacFac_GaussWeights,
            N1_i[gp1],
            N2_i[gp2],
            *stiffmat11,
            *stiffmat12,
            *stiffmat21,
            *stiffmat22);
      }

    } // end gauss quadrature loop (element 2)
  } // end gauss quadrature loop (element 1)


  // apply constant prefactor
  force_pot1.Scale(prefactor);
  force_pot2.Scale(prefactor);

  if ( stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL )
  {
    ScaleStiffpotAnalyticContributions(
        prefactor,
        *stiffmat11,
        *stiffmat12,
        *stiffmat21,
        *stiffmat22);
  }

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
EvaluateStiffpotAnalyticContributions_LargeSepApprox(
    LINALG::TMatrix<double, 3, 1> const& dist,
    double const& norm_dist,
    double const& norm_dist_exp1,
    double q1q2_JacFac_GaussWeights,
    LINALG::Matrix<1, numnodes*numnodalvalues> const& N1_i_GP1,
    LINALG::Matrix<1, numnodes*numnodalvalues> const& N2_i_GP2,
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22) const
{

  //********************************************************************
  // calculate stiffpot1
  //********************************************************************
  // auxiliary variables (same for both elements)
  double norm_dist_exp2 = (m_+2) * std::pow(norm_dist,-m_-4);

  LINALG::TMatrix<double, 3, 3> dist_dist_T(true);

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
        stiffmat11(3*i+idim,3*j+idim) -=
            norm_dist_exp1 * N1_i_GP1(i)*N1_i_GP1(j) *  q1q2_JacFac_GaussWeights;

        for (unsigned int jdim=0; jdim<3; ++jdim)
        {
          stiffmat11(3*i+idim,3*j+jdim) +=
              norm_dist_exp2 * N1_i_GP1(i) * dist_dist_T(idim,jdim) * N1_i_GP1(j) *  q1q2_JacFac_GaussWeights;
        }
      }
    }

    // d (Res_1) / d (d_2)
    for (unsigned int j=0; j<(numnodes*numnodalvalues); ++j)
    {

      for (unsigned int idim=0; idim<3; ++idim)
      {
        stiffmat12(3*i+idim,3*j+idim) +=
            norm_dist_exp1 * N1_i_GP1(i)*N2_i_GP2(j) *  q1q2_JacFac_GaussWeights;

        for (unsigned int jdim=0; jdim<3; ++jdim)
        {
          stiffmat12(3*i+idim,3*j+jdim) -=
              norm_dist_exp2 * N1_i_GP1(i) * dist_dist_T(idim,jdim) * N2_i_GP2(j) *  q1q2_JacFac_GaussWeights;
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
        stiffmat21(3*i+idim,3*j+idim) +=
            norm_dist_exp1 * N2_i_GP2(i)*N1_i_GP1(j) *  q1q2_JacFac_GaussWeights;

        for (unsigned int jdim=0; jdim<3; ++jdim)
        {
          stiffmat21(3*i+idim,3*j+jdim) -=
              norm_dist_exp2 * N2_i_GP2(i) * dist_dist_T(idim,jdim) * N1_i_GP1(j) *  q1q2_JacFac_GaussWeights;
        }
      }
    }

    // d (Res_2) / d (d_2)
    for (unsigned int j=0; j<(numnodes*numnodalvalues); ++j)
    {

      for (unsigned int idim=0; idim<3; ++idim)
      {
        stiffmat22(3*i+idim,3*j+idim) -=
            norm_dist_exp1 * N2_i_GP2(i)*N2_i_GP2(j) *  q1q2_JacFac_GaussWeights;

        for (unsigned int jdim=0; jdim<3; ++jdim)
        {
          stiffmat22(3*i+idim,3*j+jdim) +=
              norm_dist_exp2 * N2_i_GP2(i) * dist_dist_T(idim,jdim) * N2_i_GP2(j) *  q1q2_JacFac_GaussWeights;
        }
      }
    }

  }


}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
EvaluateFpotandStiffpot_DoubleLengthSpecific_SmallSepApprox(
    LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1>& force_pot1,
    LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1>& force_pot2,
    LINALG::SerialDenseMatrix* stiffmat11,
    LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21,
    LINALG::SerialDenseMatrix* stiffmat22)
{
  // safety check
  if ( m_<3.5 )
    dserror("This strategy to evaluate the interaction potential is not applicable for exponents "
        "of the point potential law smaller than 3.5!");

  // Set gauss integration rule
  DRT::UTILS::GaussRule1D gaussrule = GetGaussRule();

  // Get gauss points (gp) for integration
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);
  // number of gps
  int numgp=gausspoints.nquad;

  // vectors for shape functions and their derivatives
  // Attention: these are individual shape function values, NOT shape function matrices
  // values at all gauss points are stored in advance (more efficient due to double integral)
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N1_i(numgp);        // = N1_i
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N2_i(numgp);        // = N2_i
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N1_i_xi(numgp);     // = N1_i,xi
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N2_i_xi(numgp);     // = N2_i,eta

  // Evaluate shape functions at gauss points and store values
  GetShapeFunctions(N1_i, N2_i, N1_i_xi, N2_i_xi, gausspoints);

  // coords and derivatives of the two gauss points
  LINALG::TMatrix<T, 3, 1> r1(true);                               // = r1
  LINALG::TMatrix<T, 3, 1> r2(true);                               // = r2
  LINALG::TMatrix<T, 3, 1> dist(true);                             // = r1-r2
  T norm_dist = 0.0;                                                // = |r1-r2|
  T gap = 0.0;                                                   // = |r1-r2|-R1-R2

  // evaluate charge/particle densities from DLINE charge condition specified in input file
  double q1 = linechargeconds_[0]->GetDouble("val");
  double q2 = linechargeconds_[1]->GetDouble("val");

  // Evaluation of the Gamma-Function term:
  // gamma(nue-3.5)*gamma(0.5*(nue-1))/gamma(nue-2)/gamma(0.5*nue-1)
  double C = 0.0;

  // safety check via split of exponent m in integer and fractional part
  double integerpart;
  if ( std::modf(m_, &integerpart) != 0.0 )
    dserror("You specified a non-integer exponent of the point potential law!");

  switch ( (int) m_ )
  {
    case 4: C=1.570796326794897; break;
    case 5: C=0.5; break;
    case 6: C=0.294524311274043; break;
    case 7: C=0.208333333333333; break;
    case 8: C=0.161067982727992; break;
    case 9: C=0.13125; break;
    case 10: C=0.110734238125495; break;
    case 11: C=0.095758928571429; break;
    case 12: C=0.084348345447154; break;
    default: dserror("Gamma-Function values not known for this exponent m of potential law"); break;
  }


  // auxiliary variables
  LINALG::TMatrix<T, 3, 1> fpot_tmp(true);

  double prefactor = k_ * 2* M_PI* (m_-3.5) * pow(m_-2,-2) *
      std::sqrt( 2*radius1_*radius2_/(radius1_+radius2_) ) * C;

  // loop over gauss points of element 1
  for (int gp1 = 0; gp1 < gausspoints.nquad; ++gp1)
  {
    // compute coord vector
    ComputeCoords(r1,N1_i[gp1],ele1pos_);

    double jacobifac1 = BeamElement1()->GetJacobiFacAtXi(gausspoints.qxg[gp1][0]);

    // loop over gauss points of element 2
    for (int gp2 = 0; gp2 < gausspoints.nquad; ++gp2)
    {
      // compute coord vector
      ComputeCoords(r2,N2_i[gp2],ele2pos_);


      dist = FADUTILS::DiffVector(r1,r2);

      norm_dist = FADUTILS::VectorNorm<3>(dist);

      gap = norm_dist - radius1_ - radius2_;


      // auxiliary variables to store pre-calculated common terms
      T gap_exp1 = 0.0;

      if ( gap >= 0.0 )
        gap_exp1 = std::pow(gap,-m_+2.5);
      else
        dserror("gap<=0!");

      double q1q2_JacFac_GaussWeights =
          q1 * q2 * jacobifac1 * BeamElement2()->GetJacobiFacAtXi(gausspoints.qxg[gp2][0])
          * gausspoints.qwgt[gp1] * gausspoints.qwgt[gp2];

      // auxiliary term, same for both element forces
      for (unsigned int i=0; i<3; i++)
        fpot_tmp(i) = q1q2_JacFac_GaussWeights * dist(i) / norm_dist * gap_exp1;


      //********************************************************************
      // calculate fpot1: force on element 1
      //********************************************************************
      // sum up the contributions of all nodes (in all dimensions)
      for (unsigned int i=0; i<(numnodes*numnodalvalues); ++i)
      {
        // loop over dimensions
        for (unsigned int j=0; j<3; ++j)
        {
          force_pot1(3*i+j) -= N1_i[gp1](i)*fpot_tmp(j);
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
          force_pot2(3*i+j) += N2_i[gp2](i)*fpot_tmp(j);
        }
      }

      // evaluate analytic contributions to linearization
      if ( stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL )
      {
        EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
            dist,
            norm_dist,
            gap,
            gap_exp1,
            q1q2_JacFac_GaussWeights,
            N1_i[gp1],
            N2_i[gp2],
            *stiffmat11,
            *stiffmat12,
            *stiffmat21,
            *stiffmat22);
      }

    }// end: loop over gauss points of element 2
  }// end: loop over gauss points of element 1

  // apply constant prefactor
  force_pot1.Scale(prefactor);
  force_pot2.Scale(prefactor);

  if ( stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL )
  {
    ScaleStiffpotAnalyticContributions(
        prefactor,
        *stiffmat11,
        *stiffmat12,
        *stiffmat21,
        *stiffmat22);
  }

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
      LINALG::TMatrix<double, 3, 1> const& dist,
      double const& norm_dist,
      double const& gap,
      double const& gap_exp1,
      double q1q2_JacFac_GaussWeights,
      LINALG::Matrix<1, numnodes*numnodalvalues> const& N1_i_GP1,
      LINALG::Matrix<1, numnodes*numnodalvalues> const& N2_i_GP2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const
{
  //********************************************************************
  // calculate stiffpot1
  //********************************************************************
  // auxiliary variables (same for both elements)
  double gap_exp2 = std::pow(gap,-m_+1.5);

  double aux_fac1 = gap_exp1 / norm_dist * q1q2_JacFac_GaussWeights;
  double aux_fac2 = ( gap_exp1 * std::pow(norm_dist,-3)
                     + (m_-2.5) * gap_exp2 * std::pow(norm_dist,-2) )* q1q2_JacFac_GaussWeights;

  LINALG::TMatrix<double, 3, 3> dist_dist_T(true);

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
        stiffmat11(3*i+idim,3*j+idim) -=
            aux_fac1 * N1_i_GP1(i)*N1_i_GP1(j);

        for (unsigned int jdim=0; jdim<3; ++jdim)
        {
          stiffmat11(3*i+idim,3*j+jdim) +=
              aux_fac2 * N1_i_GP1(i) * dist_dist_T(idim,jdim) * N1_i_GP1(j);
        }
      }
    }

    // d (Res_1) / d (d_2)
    for (unsigned int j=0; j<(numnodes*numnodalvalues); ++j)
    {

      for (unsigned int idim=0; idim<3; ++idim)
      {
        stiffmat12(3*i+idim,3*j+idim) +=
            aux_fac1 * N1_i_GP1(i)*N2_i_GP2(j);

        for (unsigned int jdim=0; jdim<3; ++jdim)
        {
          stiffmat12(3*i+idim,3*j+jdim) -=
              aux_fac2 * N1_i_GP1(i) * dist_dist_T(idim,jdim) * N2_i_GP2(j);
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
        stiffmat21(3*i+idim,3*j+idim) +=
            aux_fac1 * N2_i_GP2(i)*N1_i_GP1(j);

        for (unsigned int jdim=0; jdim<3; ++jdim)
        {
          stiffmat21(3*i+idim,3*j+jdim) -=
              aux_fac2 * N2_i_GP2(i) * dist_dist_T(idim,jdim) * N1_i_GP1(j);
        }
      }
    }

    // d (Res_2) / d (d_2)
    for (unsigned int j=0; j<(numnodes*numnodalvalues); ++j)
    {

      for (unsigned int idim=0; idim<3; ++idim)
      {
        stiffmat22(3*i+idim,3*j+idim) -=
            aux_fac1 * N2_i_GP2(i)*N2_i_GP2(j);

        for (unsigned int jdim=0; jdim<3; ++jdim)
        {
          stiffmat22(3*i+idim,3*j+jdim) +=
              aux_fac2 * N2_i_GP2(i) * dist_dist_T(idim,jdim) * N2_i_GP2(j);
        }
      }
    }

  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::ScaleStiffpotAnalyticContributions(
    double const& scalefactor,
    LINALG::SerialDenseMatrix& stiffmat11,
    LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21,
    LINALG::SerialDenseMatrix& stiffmat22) const
{
  stiffmat11.Scale(scalefactor);
  stiffmat12.Scale(scalefactor);
  stiffmat21.Scale(scalefactor);
  stiffmat22.Scale(scalefactor);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
CalcStiffmatAutomaticDifferentiationIfRequired(
      LINALG::TMatrix<Sacado::Fad::DFad<double>, 3*numnodes*numnodalvalues, 1>& force_pot1,
      LINALG::TMatrix<Sacado::Fad::DFad<double>, 3*numnodes*numnodalvalues, 1>& force_pot2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const
{
  for (unsigned int idof=0; idof<3*numnodes*numnodalvalues; ++idof)
  {
    for (unsigned int jdof=0; jdof<3*numnodes*numnodalvalues; ++jdof)
    {
      // d (Res_1) / d (d_1)
      stiffmat11(idof,jdof) = force_pot1(idof).dx(jdof);

      // d (Res_1) / d (d_2)
      stiffmat12(idof,jdof) = force_pot1(idof).dx(3*numnodes*numnodalvalues+jdof);

      // d (Res_2) / d (d_1)
      stiffmat21(idof,jdof) = force_pot2(idof).dx(jdof);

      // d (Res_2) / d (d_2)
      stiffmat22(idof,jdof) = force_pot2(idof).dx(3*numnodes*numnodalvalues+jdof);
    }
  }

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::Print(std::ostream& out) const
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
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

  // Todo difficulty here is that the same element pair is evaluated more than once
  //      to be more precise, once for every common potlaw;
  //      contribution of previous evaluations is overwritten if multiple potlaws are applied

}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::GetShapeFunctions(
    std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i,
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
    /* TODO hard set distype to line2 in case of numnodalvalues_=2 because
     *  only 3rd order Hermite interpolation is used (always 2 nodes) */
    const DRT::Element::DiscretizationType distypeherm = DRT::Element::line2;

    for (int gp=0; gp<gausspoints.nquad; ++gp)
    {
      // get values and derivatives of shape functions
      DRT::UTILS::shape_function_hermite_1D(N1_i[gp], gausspoints.qxg[gp][0], ele1length_, distypeherm);
      DRT::UTILS::shape_function_hermite_1D(N2_i[gp], gausspoints.qxg[gp][0], ele2length_, distypeherm);
      DRT::UTILS::shape_function_hermite_1D_deriv1(N1_i_xi[gp], gausspoints.qxg[gp][0], ele1length_, distypeherm);
      DRT::UTILS::shape_function_hermite_1D_deriv1(N2_i_xi[gp], gausspoints.qxg[gp][0], ele2length_, distypeherm);
    }
  }
  else
    dserror("Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) values are valid!");

  return;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::ComputeCoords(
    LINALG::TMatrix<T,3,1>& r,
    const LINALG::Matrix<1,numnodes*numnodalvalues>& N_i,
    const LINALG::TMatrix<T,3*numnodes*numnodalvalues,1> elepos)
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
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::ResetState(
    const std::vector<double>& centerline_dofvec_ele1,
    const std::vector<double>& centerline_dofvec_ele2)
{
  for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i)
  {
    ele1pos_(i) = centerline_dofvec_ele1[i];
    ele2pos_(i) = centerline_dofvec_ele2[i];
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::SetAutomaticDifferentiationVariablesIfRequired(
    LINALG::TMatrix<Sacado::Fad::DFad<double>,3*numnodes*numnodalvalues,1> & ele1centerlinedofvec,
    LINALG::TMatrix<Sacado::Fad::DFad<double>,3*numnodes*numnodalvalues,1> & ele2centerlinedofvec)
{
  // The 2*3*numnodes*numnodalvalues primary DoFs consist of all nodal positions and tangents
  for (unsigned int i=0;i<3*numnodes*numnodalvalues;i++)
    ele1centerlinedofvec(i).diff(i,2*3*numnodes*numnodalvalues);

  for (unsigned int i=0;i<3*numnodes*numnodalvalues;i++)
    ele2centerlinedofvec(i).diff(3*numnodes*numnodalvalues+i,2*3*numnodes*numnodalvalues);
}


// explicit template instantiations
template class BEAMINTERACTION::BeamToBeamPotentialPair<2,1,double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<2,1,Sacado::Fad::DFad<double> >;
template class BEAMINTERACTION::BeamToBeamPotentialPair<3,1,double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<3,1,Sacado::Fad::DFad<double> >;
template class BEAMINTERACTION::BeamToBeamPotentialPair<4,1,double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<4,1,Sacado::Fad::DFad<double> >;
template class BEAMINTERACTION::BeamToBeamPotentialPair<5,1,double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<5,1,Sacado::Fad::DFad<double> >;
template class BEAMINTERACTION::BeamToBeamPotentialPair<2,2,double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<2,2,Sacado::Fad::DFad<double> >;
