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
#include "beaminteraction_geometry_utils.H"

#include "../drt_beam3/beam3_base.H"
#include "../drt_beam3/beam_spatial_discretization_utils.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/largerotations.H"

#include "../drt_inpar/inpar_beampotential.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"

#include "../headers/FAD_utils.H"

#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include <Sacado.hpp>


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::BeamToBeamPotentialPair()
    : BeamPotentialPair(),
      beam_element1_(NULL),
      beam_element2_(NULL),
      time_(0.0),
      k_(0.0),
      m_(0.0),
      ele1length_(0.0),
      ele2length_(0.0),
      radius1_(0.0),
      radius2_(0.0),
      interaction_potential_(0.0)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::Setup()
{
  CheckInit();

  // call setup of base class first
  BeamPotentialPair::Setup();


  ele1pos_.Clear();
  ele2pos_.Clear();


  // get initial length of beam elements
  beam_element1_ = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element1());
  ele1length_ = BeamElement1()->RefLength();
  beam_element2_ = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element2());
  ele2length_ = BeamElement2()->RefLength();

  radius1_ = BeamElement1()->GetCircularCrossSectionRadiusForInteractions();
  radius2_ = BeamElement2()->GetCircularCrossSectionRadiusForInteractions();

  if (Element1()->ElementType() != Element2()->ElementType())
    dserror(
        "The class BeamToBeamPotentialPair currently only supports element "
        "pairs of the same beam element type!");

  // Todo check this prerequisite, still necessary?
  if (Element1()->Id() >= Element2()->Id())
  {
    dserror(
        "Element 1 has to have the smaller element-ID. Ele1GID: %d, Ele2GID: %d."
        " Adapt your search algorithm!",
        Element1()->Id(), Element2()->Id());
  }

  // initialize line charge conditions applied to element1 and element2
  linechargeconds_.resize(2);

  issetup_ = true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
bool BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22,
    const std::vector<DRT::Condition*> linechargeconds, const double k, const double m)
{
  // nothing to do in case of k==0.0
  if (k == 0.0) return false;


  // set class variables
  if (linechargeconds.size() == 2)
  {
    for (unsigned int i = 0; i < 2; ++i)
    {
      if (linechargeconds[i]->Type() == DRT::Condition::BeamPotential_LineChargeDensity)
        linechargeconds_[i] = linechargeconds[i];
      else
        dserror(
            "Provided line charge condition is not of correct type"
            "BeamPotential_LineChargeDensity!");
    }
  }
  else
    dserror("Expected TWO dline charge conditions!");

  k_ = k;
  m_ = m;



  const unsigned int dim1 = 3 * numnodes * numnodalvalues;
  const unsigned int dim2 = 3 * numnodes * numnodalvalues;

  LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> force_pot1(true);
  LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> force_pot2(true);


  if (stiffmat11 != NULL) stiffmat11->Shape(dim1, dim1);
  if (stiffmat12 != NULL) stiffmat12->Shape(dim1, dim2);
  if (stiffmat21 != NULL) stiffmat21->Shape(dim2, dim1);
  if (stiffmat22 != NULL) stiffmat22->Shape(dim2, dim2);


  // compute the values for element residual vectors ('force') and linearizations ('stiff')
  switch (Params()->Strategy())
  {
    case INPAR::BEAMPOTENTIAL::strategy_doublelengthspec_largesepapprox:
    {
      EvaluateFpotandStiffpot_LargeSepApprox(
          force_pot1, force_pot2, stiffmat11, stiffmat12, stiffmat21, stiffmat22);
      break;
    }

    case INPAR::BEAMPOTENTIAL::strategy_doublelengthspec_smallsepapprox:
    {
      EvaluateFpotandStiffpot_DoubleLengthSpecific_SmallSepApprox(
          force_pot1, force_pot2, stiffmat11, stiffmat12, stiffmat21, stiffmat22);
      break;
    }

    case INPAR::BEAMPOTENTIAL::strategy_singlelengthspec_smallsepapprox:
    {
      EvaluateFpotandStiffpot_SingleLengthSpecific_SmallSepApprox(
          force_pot1, force_pot2, stiffmat11, stiffmat12, stiffmat21, stiffmat22);
      break;
    }

    default:
      dserror("Invalid strategy to evaluate beam interaction potential!");
  }


  // resize variables and fill with pre-computed values
  if (forcevec1 != NULL)
  {
    forcevec1->Size(dim1);
    for (unsigned int i = 0; i < dim1; ++i) (*forcevec1)(i) = FADUTILS::CastToDouble(force_pot1(i));
  }
  if (forcevec2 != NULL)
  {
    forcevec2->Size(dim2);
    for (unsigned int i = 0; i < dim2; ++i) (*forcevec2)(i) = FADUTILS::CastToDouble(force_pot2(i));
  }

  return true;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues,
    T>::EvaluateFpotandStiffpot_LargeSepApprox(LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>&
                                                   force_pot1,
    LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // prepare differentiation via FAD if desired
  SetAutomaticDifferentiationVariablesIfRequired(ele1pos_, ele2pos_);

  // get cutoff radius
  const double cutoff_radius = Params()->CutoffRadius();

  // number of integration segments per element
  const unsigned int num_integration_segments = Params()->NumberIntegrationSegments();

  // Set Gauss integration rule applied in each integration segment
  DRT::UTILS::GaussRule1D gaussrule = GetGaussRule();

  // Get Gauss points (gp) for integration
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);
  // number of Gauss points per integration segment and in total per element
  int numgp_persegment = gausspoints.nquad;
  int numgp_perelement = num_integration_segments * numgp_persegment;

  // vectors for shape functions
  // Attention: these are individual shape function values, NOT shape function matrices
  std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>> N1_i(numgp_persegment);
  std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>> N2_i(numgp_persegment);

  // Evaluate shape functions at gauss points and store values
  // Todo think about pre-computing and storing values for inner Gauss point loops here

  // coords of the two gauss points
  LINALG::TMatrix<T, 3, 1> r1(true);    // = r1
  LINALG::TMatrix<T, 3, 1> r2(true);    // = r2
  LINALG::TMatrix<T, 3, 1> dist(true);  // = r1-r2
  T norm_dist = 0.0;                    // = |r1-r2|

  // evaluate charge densities from DLINE charge condition specified in input file
  double q1 = linechargeconds_[0]->GetDouble("val");
  double q2 = linechargeconds_[1]->GetDouble("val");

  // evaluate function in time if specified in line charge conditions
  // TODO allow for functions in space, i.e. varying charge along beam centerline
  int function_number = linechargeconds_[0]->GetInt("funct");

  if (function_number != -1)
    q1 *= DRT::Problem::Instance()->Funct(function_number - 1).EvaluateTime(time_);

  function_number = linechargeconds_[1]->GetInt("funct");

  if (function_number != -1)
    q2 *= DRT::Problem::Instance()->Funct(function_number - 1).EvaluateTime(time_);


  // auxiliary variable
  LINALG::TMatrix<T, 3, 1> fpot_tmp(true);

  // determine prefactor of the integral (depends on whether surface or volume potential is applied)
  double prefactor = k_ * m_;

  switch (Params()->PotentialType())
  {
    case INPAR::BEAMPOTENTIAL::beampot_surf:
      prefactor *= 4 * radius1_ * radius2_ * M_PI * M_PI;
      break;
    case INPAR::BEAMPOTENTIAL::beampot_vol:
      prefactor *= radius1_ * radius1_ * radius2_ * radius2_ * M_PI * M_PI;
      break;
    default:
      dserror(
          "No valid BEAMPOTENTIAL_TYPE specified. Choose either Surface or Volume in input file!");
  }

  // reset interaction potential of this pair
  interaction_potential_ = 0.0;

  // prepare data storage for vtk visualization
  centerline_coords_GP1_.resize(numgp_perelement);
  centerline_coords_GP2_.resize(numgp_perelement);
  forces_pot_GP1_.resize(numgp_perelement, LINALG::TMatrix<double, 3, 1>(true));
  forces_pot_GP2_.resize(numgp_perelement, LINALG::TMatrix<double, 3, 1>(true));
  moments_pot_GP1_.resize(numgp_perelement, LINALG::TMatrix<double, 3, 1>(true));
  moments_pot_GP2_.resize(numgp_perelement, LINALG::TMatrix<double, 3, 1>(true));

  for (unsigned int isegment1 = 0; isegment1 < num_integration_segments; ++isegment1)
  {
    // compute element parameter coordinate for lower and upper limit of current integration segment
    double integration_segment1_lower_limit = -1.0 + isegment1 * 2.0 / num_integration_segments;
    double integration_segment1_upper_limit =
        -1.0 + (isegment1 + 1) * 2.0 / num_integration_segments;

    double jacobifactor_segment1 =
        0.5 * (integration_segment1_upper_limit - integration_segment1_lower_limit);

    DRT::UTILS::BEAM::EvaluateShapeFunctionsAllGPs<numnodes, numnodalvalues>(gausspoints, N1_i,
        BeamElement1()->Shape(), ele1length_, integration_segment1_lower_limit,
        integration_segment1_upper_limit);

    for (unsigned int isegment2 = 0; isegment2 < num_integration_segments; ++isegment2)
    {
      // compute element parameter coordinate for lower and upper limit of current integration
      // segment
      double integration_segment2_lower_limit = -1.0 + isegment2 * 2.0 / num_integration_segments;
      double integration_segment2_upper_limit =
          -1.0 + (isegment2 + 1) * 2.0 / num_integration_segments;

      double jacobifactor_segment2 =
          0.5 * (integration_segment2_upper_limit - integration_segment2_lower_limit);

      DRT::UTILS::BEAM::EvaluateShapeFunctionsAllGPs<numnodes, numnodalvalues>(gausspoints, N2_i,
          BeamElement2()->Shape(), ele2length_, integration_segment2_lower_limit,
          integration_segment2_upper_limit);


      // loop over Gauss points in current segment on ele1
      for (int igp1 = 0; igp1 < numgp_persegment; ++igp1)
      {
        int igp1_total = isegment1 * numgp_persegment + igp1;

        // Get location of GP in element parameter space xi \in [-1;1]
        const double xi_GP1_tilde = gausspoints.qxg[igp1][0];

        /* do a mapping into integration segment, i.e. coordinate transformation
         * note: this has no effect if integration interval is [-1;1] */
        const double xi_GP1 = 0.5 * ((1.0 - xi_GP1_tilde) * integration_segment1_lower_limit +
                                        (1.0 + xi_GP1_tilde) * integration_segment1_upper_limit);

        ComputeCenterlinePosition(r1, N1_i[igp1], ele1pos_);

        // store for vtk visualization
        centerline_coords_GP1_[igp1_total] = FADUTILS::CastToDouble<T, 3, 1>(r1);

        double jacobifac1 = BeamElement1()->GetJacobiFacAtXi(xi_GP1);

        // loop over Gauss points in current segment on ele2
        for (int igp2 = 0; igp2 < numgp_persegment; ++igp2)
        {
          int igp2_total = isegment2 * numgp_persegment + igp2;

          // Get location of GP in element parameter space xi \in [-1;1]
          const double xi_GP2_tilde = gausspoints.qxg[igp2][0];

          /* do a mapping into integration segment, i.e. coordinate transformation
           * note: this has no effect if integration interval is [-1;1] */
          const double xi_GP2 = 0.5 * ((1.0 - xi_GP2_tilde) * integration_segment2_lower_limit +
                                          (1.0 + xi_GP2_tilde) * integration_segment2_upper_limit);

          // compute coord vector
          ComputeCenterlinePosition(r2, N2_i[igp2], ele2pos_);

          // store for vtk visualization
          centerline_coords_GP2_[igp2_total] = FADUTILS::CastToDouble<T, 3, 1>(r2);

          double jacobifac2 = BeamElement2()->GetJacobiFacAtXi(xi_GP2);

          dist = FADUTILS::DiffVector(r1, r2);

          norm_dist = FADUTILS::VectorNorm(dist);

          // check cutoff criterion: if specified, contributions are neglected at larger separation
          if (cutoff_radius != -1.0 and FADUTILS::CastToDouble(norm_dist) > cutoff_radius) continue;

          // auxiliary variables to store pre-calculated common terms
          T norm_dist_exp1 = 0.0;
          if (norm_dist != 0.0)
          {
            norm_dist_exp1 = std::pow(norm_dist, -m_ - 2);
          }
          else
          {
            dserror(
                "\n|r1-r2|=0 ! Interacting points are identical! Potential law not defined in this"
                " case! Think about shifting nodes in unconverged state?!");
          }

          double q1q2_JacFac_GaussWeights = q1 * q2 * jacobifac1 * jacobifactor_segment1 *
                                            jacobifac2 * jacobifactor_segment2 *
                                            gausspoints.qwgt[igp1] * gausspoints.qwgt[igp2];

          // compute fpot_tmp here, same for both element forces
          for (unsigned int i = 0; i < 3; ++i)
            fpot_tmp(i) = q1q2_JacFac_GaussWeights * norm_dist_exp1 * dist(i);

          //********************************************************************
          // calculate fpot1: force on element 1
          //********************************************************************
          // sum up the contributions of all nodes (in all dimensions)
          for (unsigned int i = 0; i < (numnodes * numnodalvalues); ++i)
          {
            // loop over dimensions
            for (unsigned int j = 0; j < 3; ++j)
            {
              force_pot1(3 * i + j) -= N1_i[igp1](i) * fpot_tmp(j);
            }
          }

          //********************************************************************
          // calculate fpot2: force on element 2
          //********************************************************************
          // sum up the contributions of all nodes (in all dimensions)
          for (unsigned int i = 0; i < (numnodes * numnodalvalues); ++i)
          {
            // loop over dimensions
            for (unsigned int j = 0; j < 3; ++j)
            {
              force_pot2(3 * i + j) += N2_i[igp2](i) * fpot_tmp(j);
            }
          }

          // evaluate analytic contributions to linearization
          if (stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and
              stiffmat22 != NULL)
          {
            EvaluateStiffpotAnalyticContributions_LargeSepApprox(dist, norm_dist, norm_dist_exp1,
                q1q2_JacFac_GaussWeights, N1_i[igp1], N2_i[igp2], *stiffmat11, *stiffmat12,
                *stiffmat21, *stiffmat22);
          }

          // store for vtk visualization
          forces_pot_GP1_[igp1_total].Update(
              1.0 * prefactor * q1 * q2 * FADUTILS::CastToDouble(norm_dist_exp1) * jacobifac2 *
                  jacobifactor_segment2 * gausspoints.qwgt[igp2],
              FADUTILS::CastToDouble<T, 3, 1>(dist), 1.0);
          forces_pot_GP2_[igp2_total].Update(
              -1.0 * prefactor * q1 * q2 * FADUTILS::CastToDouble(norm_dist_exp1) * jacobifac1 *
                  jacobifactor_segment1 * gausspoints.qwgt[igp1],
              FADUTILS::CastToDouble<T, 3, 1>(dist), 1.0);

          // store for energy output
          interaction_potential_ += prefactor / m_ * q1q2_JacFac_GaussWeights *
                                    std::pow(FADUTILS::CastToDouble(norm_dist), -m_);

        }  // end gauss quadrature loop (element 2)
      }    // end gauss quadrature loop (element 1)

    }  // end: loop over integration segments of element 2
  }    // end: loop over integration segments of element 1

  // apply constant prefactor
  force_pot1.Scale(prefactor);
  force_pot2.Scale(prefactor);

  if (stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL)
  {
    ScaleStiffpotAnalyticContributionsIfRequired(
        prefactor, *stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);

    CalcStiffmatAutomaticDifferentiationIfRequired(
        force_pot1, force_pot2, *stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues,
    T>::EvaluateStiffpotAnalyticContributions_LargeSepApprox(LINALG::TMatrix<double, 3, 1> const&
                                                                 dist,
    double const& norm_dist, double const& norm_dist_exp1, double q1q2_JacFac_GaussWeights,
    LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N1_i_GP1,
    LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N2_i_GP2,
    LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
{
  //********************************************************************
  // calculate stiffpot1
  //********************************************************************
  // auxiliary variables (same for both elements)
  double norm_dist_exp2 = (m_ + 2) * std::pow(norm_dist, -m_ - 4);

  LINALG::TMatrix<double, 3, 3> dist_dist_T(true);

  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
    {
      dist_dist_T(i, j) = dist(i) * dist(j);
      if (i != j) dist_dist_T(j, i) = dist_dist_T(i, j);
    }
  }

  for (unsigned int i = 0; i < (numnodes * numnodalvalues); ++i)
  {
    // d (Res_1) / d (d_1)
    for (unsigned int j = 0; j < (numnodes * numnodalvalues); ++j)
    {
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        stiffmat11(3 * i + idim, 3 * j + idim) -=
            norm_dist_exp1 * N1_i_GP1(i) * N1_i_GP1(j) * q1q2_JacFac_GaussWeights;

        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          stiffmat11(3 * i + idim, 3 * j + jdim) += norm_dist_exp2 * N1_i_GP1(i) *
                                                    dist_dist_T(idim, jdim) * N1_i_GP1(j) *
                                                    q1q2_JacFac_GaussWeights;
        }
      }
    }

    // d (Res_1) / d (d_2)
    for (unsigned int j = 0; j < (numnodes * numnodalvalues); ++j)
    {
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        stiffmat12(3 * i + idim, 3 * j + idim) +=
            norm_dist_exp1 * N1_i_GP1(i) * N2_i_GP2(j) * q1q2_JacFac_GaussWeights;

        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          stiffmat12(3 * i + idim, 3 * j + jdim) -= norm_dist_exp2 * N1_i_GP1(i) *
                                                    dist_dist_T(idim, jdim) * N2_i_GP2(j) *
                                                    q1q2_JacFac_GaussWeights;
        }
      }
    }
  }

  //********************************************************************
  // calculate stiffpot2
  //********************************************************************
  for (unsigned int i = 0; i < (numnodes * numnodalvalues); ++i)
  {
    // d (Res_2) / d (d_1)
    for (unsigned int j = 0; j < (numnodes * numnodalvalues); ++j)
    {
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        stiffmat21(3 * i + idim, 3 * j + idim) +=
            norm_dist_exp1 * N2_i_GP2(i) * N1_i_GP1(j) * q1q2_JacFac_GaussWeights;

        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          stiffmat21(3 * i + idim, 3 * j + jdim) -= norm_dist_exp2 * N2_i_GP2(i) *
                                                    dist_dist_T(idim, jdim) * N1_i_GP1(j) *
                                                    q1q2_JacFac_GaussWeights;
        }
      }
    }

    // d (Res_2) / d (d_2)
    for (unsigned int j = 0; j < (numnodes * numnodalvalues); ++j)
    {
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        stiffmat22(3 * i + idim, 3 * j + idim) -=
            norm_dist_exp1 * N2_i_GP2(i) * N2_i_GP2(j) * q1q2_JacFac_GaussWeights;

        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          stiffmat22(3 * i + idim, 3 * j + jdim) += norm_dist_exp2 * N2_i_GP2(i) *
                                                    dist_dist_T(idim, jdim) * N2_i_GP2(j) *
                                                    q1q2_JacFac_GaussWeights;
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
    EvaluateFpotandStiffpot_DoubleLengthSpecific_SmallSepApprox(
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // safety check
  if (m_ < 3.5)
    dserror(
        "This strategy to evaluate the interaction potential is not applicable for exponents "
        "of the point potential law smaller than 3.5!");

  // prepare differentiation via FAD if desired
  SetAutomaticDifferentiationVariablesIfRequired(ele1pos_, ele2pos_);

  // get cutoff radius
  const double cutoff_radius = Params()->CutoffRadius();

  // get regularization type and separation
  const INPAR::BEAMPOTENTIAL::BeamPotentialRegularizationType regularization_type =
      Params()->RegularizationType();

  const double regularization_separation = Params()->RegularizationSeparation();

  // number of integration segments per element
  const unsigned int num_integration_segments = Params()->NumberIntegrationSegments();

  // Set Gauss integration rule applied in each integration segment
  DRT::UTILS::GaussRule1D gaussrule = GetGaussRule();

  // Get Gauss points (gp) for integration
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);
  // number of Gauss points per integration segment and in total per element
  int numgp_persegment = gausspoints.nquad;
  int numgp_perelement = num_integration_segments * numgp_persegment;

  // vectors for shape function values
  // Attention: these are individual shape function values, NOT shape function matrices
  std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>> N1_i(numgp_persegment);
  std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>> N2_i(numgp_persegment);

  // Evaluate shape functions at gauss points and store values
  // Todo think about pre-computing and storing values for inner Gauss point loops here

  // coords of the two gauss points
  LINALG::TMatrix<T, 3, 1> r1(true);    // = r1
  LINALG::TMatrix<T, 3, 1> r2(true);    // = r2
  LINALG::TMatrix<T, 3, 1> dist(true);  // = r1-r2
  T norm_dist = 0.0;                    // = |r1-r2|
  T gap = 0.0;                          // = |r1-r2|-R1-R2
  T gap_regularized = 0.0;  // modified gap if a regularization of the force law is applied

  // evaluate charge/particle densities from DLINE charge condition specified in input file
  double q1 = linechargeconds_[0]->GetDouble("val");
  double q2 = linechargeconds_[1]->GetDouble("val");

  // evaluate function in time if specified in line charge conditions
  // TODO allow for functions in space, i.e. varying charge along beam centerline
  int function_number = linechargeconds_[0]->GetInt("funct");

  if (function_number != -1)
    q1 *= DRT::Problem::Instance()->Funct(function_number - 1).EvaluateTime(time_);

  function_number = linechargeconds_[1]->GetInt("funct");

  if (function_number != -1)
    q2 *= DRT::Problem::Instance()->Funct(function_number - 1).EvaluateTime(time_);


  // Evaluation of the Gamma-Function term:
  // gamma(nue-3.5)*gamma(0.5*(nue-1))/gamma(nue-2)/gamma(0.5*nue-1)
  double C = 0.0;

  // safety check via split of exponent m in integer and fractional part
  double integerpart;
  if (std::modf(m_, &integerpart) != 0.0)
    dserror("You specified a non-integer exponent of the point potential law!");

  switch ((int)m_)
  {
    case 4:
      C = 1.570796326794897;
      break;
    case 5:
      C = 0.5;
      break;
    case 6:
      C = 0.294524311274043;
      break;
    case 7:
      C = 0.208333333333333;
      break;
    case 8:
      C = 0.161067982727992;
      break;
    case 9:
      C = 0.13125;
      break;
    case 10:
      C = 0.110734238125495;
      break;
    case 11:
      C = 0.095758928571429;
      break;
    case 12:
      C = 0.084348345447154;
      break;
    default:
      dserror("Gamma-Function values not known for this exponent m of potential law");
      break;
  }


  // reset interaction potential of this pair
  interaction_potential_ = 0.0;

  // prepare data storage for vtk visualization
  centerline_coords_GP1_.resize(numgp_perelement);
  centerline_coords_GP2_.resize(numgp_perelement);
  forces_pot_GP1_.resize(numgp_perelement, LINALG::TMatrix<double, 3, 1>(true));
  forces_pot_GP2_.resize(numgp_perelement, LINALG::TMatrix<double, 3, 1>(true));
  moments_pot_GP1_.resize(numgp_perelement, LINALG::TMatrix<double, 3, 1>(true));
  moments_pot_GP2_.resize(numgp_perelement, LINALG::TMatrix<double, 3, 1>(true));

  // auxiliary variables
  LINALG::TMatrix<T, 3, 1> fpot_tmp(true);

  double prefactor = k_ * 2 * M_PI * (m_ - 3.5) / (m_ - 2) / (m_ - 2) *
                     std::sqrt(2 * radius1_ * radius2_ / (radius1_ + radius2_)) * C;

  for (unsigned int isegment1 = 0; isegment1 < num_integration_segments; ++isegment1)
  {
    // compute element parameter coordinate for lower and upper limit of current integration segment
    double integration_segment1_lower_limit = -1.0 + isegment1 * 2.0 / num_integration_segments;
    double integration_segment1_upper_limit =
        -1.0 + (isegment1 + 1) * 2.0 / num_integration_segments;

    double jacobifactor_segment1 =
        0.5 * (integration_segment1_upper_limit - integration_segment1_lower_limit);

    DRT::UTILS::BEAM::EvaluateShapeFunctionsAllGPs<numnodes, numnodalvalues>(gausspoints, N1_i,
        BeamElement1()->Shape(), ele1length_, integration_segment1_lower_limit,
        integration_segment1_upper_limit);

    for (unsigned int isegment2 = 0; isegment2 < num_integration_segments; ++isegment2)
    {
      // compute element parameter coordinate for lower and upper limit of current integration
      // segment
      double integration_segment2_lower_limit = -1.0 + isegment2 * 2.0 / num_integration_segments;
      double integration_segment2_upper_limit =
          -1.0 + (isegment2 + 1) * 2.0 / num_integration_segments;

      double jacobifactor_segment2 =
          0.5 * (integration_segment2_upper_limit - integration_segment2_lower_limit);

      DRT::UTILS::BEAM::EvaluateShapeFunctionsAllGPs<numnodes, numnodalvalues>(gausspoints, N2_i,
          BeamElement2()->Shape(), ele2length_, integration_segment2_lower_limit,
          integration_segment2_upper_limit);

      // loop over gauss points of current segment on element 1
      for (int igp1 = 0; igp1 < numgp_persegment; ++igp1)
      {
        int igp1_total = isegment1 * numgp_persegment + igp1;

        // Get location of GP in element parameter space xi \in [-1;1]
        const double xi_GP1_tilde = gausspoints.qxg[igp1][0];

        /* do a mapping into integration segment, i.e. coordinate transformation
         * note: this has no effect if integration interval is [-1;1] */
        const double xi_GP1 = 0.5 * ((1.0 - xi_GP1_tilde) * integration_segment1_lower_limit +
                                        (1.0 + xi_GP1_tilde) * integration_segment1_upper_limit);

        // compute coord vector
        ComputeCenterlinePosition(r1, N1_i[igp1], ele1pos_);

        // store for vtk visualization
        centerline_coords_GP1_[igp1_total] = FADUTILS::CastToDouble<T, 3, 1>(r1);

        double jacobifac1 = BeamElement1()->GetJacobiFacAtXi(xi_GP1);

        // loop over gauss points of current segment on element 2
        for (int igp2 = 0; igp2 < numgp_persegment; ++igp2)
        {
          int igp2_total = isegment2 * numgp_persegment + igp2;

          // Get location of GP in element parameter space xi \in [-1;1]
          const double xi_GP2_tilde = gausspoints.qxg[igp2][0];

          /* do a mapping into integration segment, i.e. coordinate transformation
           * note: this has no effect if integration interval is [-1;1] */
          const double xi_GP2 = 0.5 * ((1.0 - xi_GP2_tilde) * integration_segment2_lower_limit +
                                          (1.0 + xi_GP2_tilde) * integration_segment2_upper_limit);

          // compute coord vector
          ComputeCenterlinePosition(r2, N2_i[igp2], ele2pos_);

          // store for vtk visualization
          centerline_coords_GP2_[igp2_total] = FADUTILS::CastToDouble<T, 3, 1>(r2);

          double jacobifac2 = BeamElement2()->GetJacobiFacAtXi(xi_GP2);

          dist = FADUTILS::DiffVector(r1, r2);

          norm_dist = FADUTILS::VectorNorm(dist);

          if (norm_dist == 0.0)
          {
            this->Print(std::cout);
            std::cout << "\nGP pair: igp1_total=" << igp1_total << " & igp2_total=" << igp2_total
                      << ": |r1-r2|=" << norm_dist;

            dserror("centerline separation |r1-r2|=0! Fatal error.");
          }

          // check cutoff criterion: if specified, contributions are neglected at larger separation
          if (cutoff_radius != -1.0 and FADUTILS::CastToDouble(norm_dist) > cutoff_radius) continue;

          gap = norm_dist - radius1_ - radius2_;



          if (regularization_type == INPAR::BEAMPOTENTIAL::regularization_none and gap <= 0.0)
          {
            this->Print(std::cout);
            std::cout << "\nGP pair: igp1_total=" << igp1_total << " & igp2_total=" << igp2_total
                      << ": gap=" << gap;

            dserror(
                "gap<=0! Force law resulting from specified interaction potential law is "
                "not defined for zero/negative gaps! Use/implement a regularization!");
          }

          gap_regularized = gap;

          if (regularization_type == INPAR::BEAMPOTENTIAL::regularization_constant and
              gap < regularization_separation)
          {
            gap_regularized = regularization_separation;
          }

          if (gap_regularized <= 0)
            dserror(
                "regularized gap <= 0! Fatal error since force law is not defined for "
                "zero/negative gaps! Use positive regularization separation!");


          // auxiliary variables to store pre-calculated common terms
          T gap_exp1 = std::pow(gap_regularized, -m_ + 2.5);

          double q1q2_JacFac_GaussWeights = q1 * q2 * jacobifac1 * jacobifactor_segment1 *
                                            jacobifac2 * jacobifactor_segment2 *
                                            gausspoints.qwgt[igp1] * gausspoints.qwgt[igp2];

          // store for energy output
          interaction_potential_ += prefactor / (m_ - 3.5) * q1q2_JacFac_GaussWeights *
                                    std::pow(FADUTILS::CastToDouble(gap_regularized), -m_ + 3.5);

          if (regularization_type == INPAR::BEAMPOTENTIAL::regularization_constant and
              gap < regularization_separation)
          {
            // potential law is linear in the regime of constant extrapolation of force law
            // add the contribution from this part of the force law
            interaction_potential_ += prefactor * q1q2_JacFac_GaussWeights *
                                      FADUTILS::CastToDouble(gap_exp1) *
                                      (regularization_separation - FADUTILS::CastToDouble(gap));
          }

          // auxiliary term, same for both element forces
          for (unsigned int i = 0; i < 3; i++)
            fpot_tmp(i) = q1q2_JacFac_GaussWeights * dist(i) / norm_dist * gap_exp1;


          //********************************************************************
          // calculate fpot1: force on element 1
          //********************************************************************
          // sum up the contributions of all nodes (in all dimensions)
          for (unsigned int i = 0; i < (numnodes * numnodalvalues); ++i)
          {
            // loop over dimensions
            for (unsigned int j = 0; j < 3; ++j)
            {
              force_pot1(3 * i + j) -= N1_i[igp1](i) * fpot_tmp(j);
            }
          }

          //********************************************************************
          // calculate fpot2: force on element 2
          //********************************************************************
          // sum up the contributions of all nodes (in all dimensions)
          for (unsigned int i = 0; i < (numnodes * numnodalvalues); ++i)
          {
            // loop over dimensions
            for (unsigned int j = 0; j < 3; ++j)
            {
              force_pot2(3 * i + j) += N2_i[igp2](i) * fpot_tmp(j);
            }
          }

          // evaluate analytic contributions to linearization
          if (stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and
              stiffmat22 != NULL)
          {
            EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(dist,
                norm_dist, gap, gap_regularized, gap_exp1, q1q2_JacFac_GaussWeights, N1_i[igp1],
                N2_i[igp2], *stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);
          }

          // store for vtk visualization
          forces_pot_GP1_[igp1_total].Update(1.0 * prefactor * q1 * q2 *
                                                 FADUTILS::CastToDouble(gap_exp1) /
                                                 FADUTILS::CastToDouble(norm_dist) * jacobifac2 *
                                                 jacobifactor_segment2 * gausspoints.qwgt[igp2],
              FADUTILS::CastToDouble<T, 3, 1>(dist), 1.0);
          forces_pot_GP2_[igp2_total].Update(-1.0 * prefactor * q1 * q2 *
                                                 FADUTILS::CastToDouble(gap_exp1) /
                                                 FADUTILS::CastToDouble(norm_dist) * jacobifac1 *
                                                 jacobifactor_segment1 * gausspoints.qwgt[igp1],
              FADUTILS::CastToDouble<T, 3, 1>(dist), 1.0);

        }  // end: loop over gauss points of element 2
      }    // end: loop over gauss points of element 1

    }  // end: loop over integration segments of element 2
  }    // end: loop over integration segments of element 1

  // apply constant prefactor
  force_pot1.Scale(prefactor);
  force_pot2.Scale(prefactor);

  if (stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL)
  {
    ScaleStiffpotAnalyticContributionsIfRequired(
        prefactor, *stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);

    CalcStiffmatAutomaticDifferentiationIfRequired(
        force_pot1, force_pot2, *stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
    EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
        LINALG::TMatrix<double, 3, 1> const& dist, double const& norm_dist, double const& gap,
        double const& gap_regularized, double const& gap_exp1, double q1q2_JacFac_GaussWeights,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N1_i_GP1,
        LINALG::TMatrix<double, 1, numnodes * numnodalvalues> const& N2_i_GP2,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
{
  //********************************************************************
  // calculate stiffpot1
  //********************************************************************
  // auxiliary variables (same for both elements)
  double gap_exp2 = std::pow(gap_regularized, -m_ + 1.5);

  if (Params()->RegularizationType() == INPAR::BEAMPOTENTIAL::regularization_constant and
      gap < Params()->RegularizationSeparation())
  {
    /* in case of constant extrapolation of force law, the derivative of the force is zero
     * and this contribution to the stiffness matrix vanishes */
    gap_exp2 = 0.0;
  }

  double aux_fac1 = gap_exp1 / norm_dist * q1q2_JacFac_GaussWeights;
  double aux_fac2 =
      (gap_exp1 * std::pow(norm_dist, -3) + (m_ - 2.5) * gap_exp2 * std::pow(norm_dist, -2)) *
      q1q2_JacFac_GaussWeights;

  LINALG::TMatrix<double, 3, 3> dist_dist_T(true);

  for (unsigned int i = 0; i < 3; ++i)
  {
    for (unsigned int j = 0; j <= i; ++j)
    {
      dist_dist_T(i, j) = dist(i) * dist(j);
      if (i != j) dist_dist_T(j, i) = dist_dist_T(i, j);
    }
  }

  for (unsigned int i = 0; i < (numnodes * numnodalvalues); ++i)
  {
    // d (Res_1) / d (d_1)
    for (unsigned int j = 0; j < (numnodes * numnodalvalues); ++j)
    {
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        stiffmat11(3 * i + idim, 3 * j + idim) -= aux_fac1 * N1_i_GP1(i) * N1_i_GP1(j);

        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          stiffmat11(3 * i + idim, 3 * j + jdim) +=
              aux_fac2 * N1_i_GP1(i) * dist_dist_T(idim, jdim) * N1_i_GP1(j);
        }
      }
    }

    // d (Res_1) / d (d_2)
    for (unsigned int j = 0; j < (numnodes * numnodalvalues); ++j)
    {
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        stiffmat12(3 * i + idim, 3 * j + idim) += aux_fac1 * N1_i_GP1(i) * N2_i_GP2(j);

        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          stiffmat12(3 * i + idim, 3 * j + jdim) -=
              aux_fac2 * N1_i_GP1(i) * dist_dist_T(idim, jdim) * N2_i_GP2(j);
        }
      }
    }
  }

  //********************************************************************
  // calculate stiffpot2
  //********************************************************************
  for (unsigned int i = 0; i < (numnodes * numnodalvalues); ++i)
  {
    // d (Res_2) / d (d_1)
    for (unsigned int j = 0; j < (numnodes * numnodalvalues); ++j)
    {
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        stiffmat21(3 * i + idim, 3 * j + idim) += aux_fac1 * N2_i_GP2(i) * N1_i_GP1(j);

        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          stiffmat21(3 * i + idim, 3 * j + jdim) -=
              aux_fac2 * N2_i_GP2(i) * dist_dist_T(idim, jdim) * N1_i_GP1(j);
        }
      }
    }

    // d (Res_2) / d (d_2)
    for (unsigned int j = 0; j < (numnodes * numnodalvalues); ++j)
    {
      for (unsigned int idim = 0; idim < 3; ++idim)
      {
        stiffmat22(3 * i + idim, 3 * j + idim) -= aux_fac1 * N2_i_GP2(i) * N2_i_GP2(j);

        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          stiffmat22(3 * i + idim, 3 * j + jdim) +=
              aux_fac2 * N2_i_GP2(i) * dist_dist_T(idim, jdim) * N2_i_GP2(j);
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
    EvaluateFpotandStiffpot_SingleLengthSpecific_SmallSepApprox(
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
        LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // safety checks
  if (m_ < 6.0)
    dserror(
        "Invalid exponent m=%f. The strategy 'SingleLengthSpecific_SmallSepApprox' to evaluate the "
        "interaction potential is only applicable for exponents m>=6 of the point potential law, "
        "e.g. van der Waals (m=6) or the repulsive part of Lennard-Jones (m=12)!",
        m_);

  if (not Params()->UseFAD())
    dserror(
        "The strategy 'SingleLengthSpecific_SmallSepApprox' to evaluate the interaction "
        "potential requires automatic differentiation via FAD!");

  if (radius1_ != radius2_)
    dserror(
        "The strategy 'SingleLengthSpecific_SmallSepApprox' to evaluate the interaction "
        "potential requires the beam radii to be identical!");


  // get cutoff radius
  const double cutoff_radius = Params()->CutoffRadius();

  /* parameter coordinate of the closest point on the master beam,
   * determined via point-to-curve projection */
  T xi_master = 0.0;

  // prepare differentiation via FAD if desired
  /* xi_master is used as additional, auxiliary primary Dof
   * since there is no closed-form expression for how xi_master depends on the 'real' primary Dofs.
   * It is determined iteratively via point-to-curve projection */
  SetAutomaticDifferentiationVariablesIfRequired(ele1pos_, ele2pos_, xi_master);


  // number of integration segments per element
  const unsigned int num_integration_segments = Params()->NumberIntegrationSegments();

  // Set Gauss integration rule applied in each integration segment
  DRT::UTILS::GaussRule1D gaussrule = GetGaussRule();

  // Get Gauss points (gp) for integration
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);

  // number of Gauss points per integration segment and in total
  int numgp_persegment = gausspoints.nquad;
  int numgp_total = num_integration_segments * numgp_persegment;


  // vectors for shape functions and their derivatives
  // Attention: these are individual shape function values, NOT shape function matrices
  // values at all gauss points are stored in advance (more efficient due to double integral)
  std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>> N_i_slave(
      numgp_persegment);  // = N1_i
  std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>> N_i_xi_slave(
      numgp_persegment);  // = N1_i,xi

  LINALG::TMatrix<T, 1, numnodes * numnodalvalues> N_i_master(true);
  LINALG::TMatrix<T, 1, numnodes * numnodalvalues> N_i_xi_master(true);
  LINALG::TMatrix<T, 1, numnodes * numnodalvalues> N_i_xixi_master(true);

  // assembled shape function matrices: Todo maybe avoid these matrices
  LINALG::TMatrix<double, 3, 3 * numnodes * numnodalvalues> N_slave(true);

  LINALG::TMatrix<T, 3, 3 * numnodes * numnodalvalues> N_master(true);
  LINALG::TMatrix<T, 3, 3 * numnodes * numnodalvalues> N_xi_master(true);
  LINALG::TMatrix<T, 3, 3 * numnodes * numnodalvalues> N_xixi_master(true);


  // coords and derivatives of the slave Gauss point and projected point on master element
  LINALG::TMatrix<T, 3, 1> r_slave(true);        // centerline position vector on slave
  LINALG::TMatrix<T, 3, 1> r_xi_slave(true);     // centerline tangent vector on slave
  LINALG::TMatrix<T, 3, 1> g1_slave(true);       // unit centerline tangent vector on slave
  LINALG::TMatrix<T, 3, 1> r_master(true);       // centerline position vector on master
  LINALG::TMatrix<T, 3, 1> r_xi_master(true);    // centerline tangent vector on master
  LINALG::TMatrix<T, 3, 1> r_xixi_master(true);  // 2nd deriv of master curve
  LINALG::TMatrix<T, 3, 1> g1_master(true);      // unit centerline tangent vector on master
  LINALG::TMatrix<T, 3, 1> dist_ul(true);        // = r_slave-r_master
  LINALG::TMatrix<T, 3, 1> normal_bl(true);      // normal vector at bilateral closest point
  T norm_normal_bl_tilde = 0.0;                  // norm of vector defining bilateral normal vector
  T gap_bl = 0.0;                                // gap of bilateral closest point
  T alpha = 0.0;                                 // mutual angle of tangent vectors
  T cos_alpha = 0.0;                             // cosine of mutual angle of tangent vectors
  T sin_alpha = 0.0;                             // sine of mutual angle of tangent vectors
  T sin_2alpha = 0.0;                            // sine of 2*mutual angle of tangent vectors
  T x = 0.0;  // distance between Gauss point and bilateral closest point on slave
  LINALG::TMatrix<T, 3, 1> aux_plane_normal(true);  // normal vector of auxiliary plane n*

  T beta = 0.0;       // auxiliary quantity
  T beta_exp2 = 0.0;  // beta^2
  T beta_exp3 = 0.0;  // beta^3
  T beta_exp4 = 0.0;  // beta^3
  T a = 0.0;          // auxiliary quantity
  T Delta = 0.0;      // auxiliary quantity

  T beta_partial_x = 0.0;
  T beta_partial_gap_bl = 0.0;
  T beta_partial_cos_alpha = 0.0;

  T a_partial_beta = 0.0;
  T a_partial_x = 0.0;
  T a_partial_gap_bl = 0.0;
  T a_partial_cos_alpha = 0.0;

  T Delta_partial_beta = 0.0;
  T Delta_partial_x = 0.0;
  T Delta_partial_gap_bl = 0.0;
  T Delta_partial_cos_alpha = 0.0;

  T interaction_potential_GP = 0.0;

  T pot_ia_partial_beta = 0.0;
  T pot_ia_partial_a = 0.0;
  T pot_ia_partial_Delta = 0.0;
  T pot_ia_partial_x = 0.0;
  T pot_ia_partial_gap_bl = 0.0;
  T pot_ia_partial_cos_alpha = 0.0;

  // components from variation of bilateral gap
  LINALG::TMatrix<T, 3, 1> gap_bl_partial_r_slave(true);
  LINALG::TMatrix<T, 3, 1> gap_bl_partial_r_master(true);
  LINALG::TMatrix<T, 3, 1> gap_bl_partial_r_xi_slave(true);
  LINALG::TMatrix<T, 3, 1> gap_bl_partial_r_xi_master(true);
  T gap_bl_partial_xi_master = 0.0;

  // components from variation of cosine of enclosed angle
  LINALG::TMatrix<T, 3, 1> cos_alpha_partial_r_xi_slave(true);
  LINALG::TMatrix<T, 3, 1> cos_alpha_partial_r_xi_master(true);
  T cos_alpha_partial_xi_master = 0.0;

  // components from variation of distance from bilateral closest point on slave
  LINALG::TMatrix<T, 3, 1> x_partial_r_slave(true);
  LINALG::TMatrix<T, 3, 1> x_partial_r_master(true);
  LINALG::TMatrix<T, 3, 1> x_partial_r_xi_slave(true);
  LINALG::TMatrix<T, 3, 1> x_partial_r_xi_master(true);

  LINALG::TMatrix<T, 3, 1> x_partial_aux_plane_normal(true);
  T x_partial_xi_master = 0.0;

  // components from variation of parameter coordinate on master beam
  LINALG::TMatrix<T, 1, 3> xi_master_partial_r_slave(true);
  LINALG::TMatrix<T, 1, 3> xi_master_partial_r_master(true);
  LINALG::TMatrix<T, 1, 3> xi_master_partial_r_xi_master(true);


  // linearization of parameter coordinate on master resulting from point-to-curve projection
  LINALG::TMatrix<T, 1, 3 * numnodes * numnodalvalues> lin_xi_master_slaveDofs(true);
  LINALG::TMatrix<T, 1, 3 * numnodes * numnodalvalues> lin_xi_master_masterDofs(true);

  /* contribution of one Gauss point (required for automatic differentiation with contributions
   * from xi_master) */
  LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> force_pot_slave_GP(true);
  LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> force_pot_master_GP(true);

  // Todo remove the following two variables, which are required for verification via FAD
  LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1> force_pot_slave_GP_calc_via_FAD(true);
  LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1> force_pot_master_GP_calc_via_FAD(true);

  // evaluate charge/particle densities from DLINE charge condition specified in input file
  double rho1 = linechargeconds_[0]->GetDouble("val");
  double rho2 = linechargeconds_[1]->GetDouble("val");

  // evaluate function in time if specified in line charge conditions
  // TODO allow for functions in space, i.e. varying charge along beam centerline
  int function_number = linechargeconds_[0]->GetInt("funct");

  if (function_number != -1)
    rho1 *= DRT::Problem::Instance()->Funct(function_number - 1).EvaluateTime(time_);

  function_number = linechargeconds_[1]->GetInt("funct");

  if (function_number != -1)
    rho2 *= DRT::Problem::Instance()->Funct(function_number - 1).EvaluateTime(time_);


  // auxiliary variables
  LINALG::TMatrix<T, 3, 1> fpot_tmp(true);
  LINALG::TMatrix<T, 3, 3> v_mat_tmp(true);

  // constant prefactor of the disk-cylinder interaction potential
  double prefactor = k_ * M_PI * M_PI * M_PI;

  switch ((int)m_)
  {
    case 6:
      prefactor *= 0.25;
      break;
    case 12:
      prefactor *= 3003 / 128;
      break;
    default:
      dserror(
          "Please implement the prefactor of the analytical disk-cylinder potential law for "
          "exponent m=%f. So far, only exponent 6 and 12 is supported.",
          m_);
      break;
  }

  //************************** DEBUG ******************************************
  //  this->Print( std::cout );
  //  std::cout << "\n\n\nStart evaluation via Gauss integration..." << std::endl;
  //*********************** END DEBUG *****************************************

  // reset interaction potential of this pair
  interaction_potential_ = 0.0;

  // prepare data storage for vtk visualization
  centerline_coords_GP1_.resize(numgp_total);
  centerline_coords_GP2_.resize(numgp_total);
  forces_pot_GP1_.resize(numgp_total, LINALG::TMatrix<double, 3, 1>(true));
  forces_pot_GP2_.resize(numgp_total, LINALG::TMatrix<double, 3, 1>(true));
  moments_pot_GP1_.resize(numgp_total, LINALG::TMatrix<double, 3, 1>(true));
  moments_pot_GP2_.resize(numgp_total, LINALG::TMatrix<double, 3, 1>(true));


  for (unsigned int isegment = 0; isegment < num_integration_segments; ++isegment)
  {
    // compute element parameter coordinate for lower and upper limit of current integration segment
    double integration_segment_lower_limit = -1.0 + isegment * 2.0 / num_integration_segments;
    double integration_segment_upper_limit = -1.0 + (isegment + 1) * 2.0 / num_integration_segments;

    double jacobifactor_segment =
        0.5 * (integration_segment_upper_limit - integration_segment_lower_limit);

    // Evaluate shape functions at Gauss points of slave element and store values
    DRT::UTILS::BEAM::EvaluateShapeFunctionsAndDerivsAllGPs<numnodes, numnodalvalues>(
        gausspoints, N_i_slave, N_i_xi_slave, BeamElement1()->Shape(), ele1length_);


    // loop over gauss points of element 1
    // so far, element 1 is always treated as the slave element!
    for (int igp = 0; igp < numgp_persegment; ++igp)
    {
      const int igp_total = isegment * numgp_persegment + igp;

      // Get location of GP in element parameter space xi \in [-1;1]
      const double xi_GP_tilde = gausspoints.qxg[igp][0];

      /* do a mapping into integration segment, i.e. coordinate transformation
       * note: this has no effect if integration interval is [-1;1] */
      const double xi_GP = 0.5 * ((1.0 - xi_GP_tilde) * integration_segment_lower_limit +
                                     (1.0 + xi_GP_tilde) * integration_segment_upper_limit);


      // compute coord vector and tangent vector on slave side
      ComputeCenterlinePosition(r_slave, N_i_slave[igp], ele1pos_);
      ComputeCenterlineTangent(r_xi_slave, N_i_xi_slave[igp], ele1pos_);

      T norm_r_xi_slave = FADUTILS::VectorNorm(r_xi_slave);

      g1_slave.Update(1.0 / norm_r_xi_slave, r_xi_slave);

      // store for vtk visualization
      centerline_coords_GP1_[igp_total] = FADUTILS::CastToDouble<T, 3, 1>(r_slave);

      double rho1rho2_JacFac_GaussWeight = rho1 * rho2 * jacobifactor_segment *
                                           BeamElement1()->GetJacobiFacAtXi(xi_GP) *
                                           gausspoints.qwgt[igp];

      //************************** DEBUG ******************************************
      //    std::cout << "\n\nGP " << igp_total << ":";
      //    std::cout << "\nr_slave: " << FADUTILS::CastToDouble<T,3,1>( r_slave );
      //    std::cout << "\nr_xi_slave: " << FADUTILS::CastToDouble<T,3,1>( r_xi_slave );
      //    std::cout << "\n|r_xi_slave|: " << FADUTILS::VectorNorm<3>(
      //    FADUTILS::CastToDouble<T,3,1>( r_xi_slave ) ); std::cout << "\ng1_slave: " <<
      //    FADUTILS::CastToDouble<T,3,1>( g1_slave ); std::cout << "\n|g1_slave|: " <<
      //    FADUTILS::VectorNorm<3>( FADUTILS::CastToDouble<T,3,1>( g1_slave ) );
      //*********************** END DEBUG *****************************************

      /* point-to-curve projection, i.e. 'unilateral' closest-point projection
       * to determine point on master beam (i.e. parameter coordinate xi_master) */
      const unsigned int num_initial_values = 9;
      double xi_master_initial_guess_values[num_initial_values] = {
          0.0, -1.0, 1.0, -0.5, 0.5, -0.75, -0.25, 0.25, 0.75};

      unsigned int i = 0;

      while (i < num_initial_values and
             (not BEAMINTERACTION::GEO::PointToCurveProjection<numnodes, numnodalvalues, T>(r_slave,
                 xi_master, xi_master_initial_guess_values[i], ele2pos_, Element2()->Shape(),
                 ele2length_)))
      {
        i++;
      }

      if (i == num_initial_values)
      {
        // Fixme
        //      dserror("Point-to-curve projection ultimately failed!");
        std::cout << "\n\nWARNING: Point-to-Curve Projection ultimately failed! "
                     "Assuming invalid projection and continue with next GP..."
                  << std::endl;
        continue;
      }

      //************************** DEBUG ******************************************
      //    if ( i!=0 )
      //    {
      //      std::cout << "\n\nINFO: Point-to-Curve Projection succeeded with initial guess "
      //          << xi_master_initial_guess_values[i] << ": xi_master="
      //          << FADUTILS::CastToDouble(xi_master) << std::endl;
      //    }
      //*********************** END DEBUG *****************************************

      //************************** DEBUG ******************************************
      //    std::cout << "\nxi_master: " << FADUTILS::CastToDouble( xi_master );
      //*********************** END DEBUG *****************************************

      // Todo: specify tolerance value in a more central place
      if (FADUTILS::Norm(xi_master) > 1.0 + 1.0e-10)
      {
        //************************** DEBUG ******************************************
        //      std::cout << "\nxi_master not in valid range ... proceed to next GP\n";
        //*********************** END DEBUG *****************************************
        continue;
      }
      else if (FADUTILS::Norm(xi_master) >= 1.0 - 1.0e-10)
      {
        dserror(
            "Point-to-curve projection yields xi_master= %f. This is a critical case "
            "since it is very close to the element boundary!",
            FADUTILS::CastToDouble(xi_master));
      }

      DRT::UTILS::BEAM::EvaluateShapeFunctionsAndDerivsAnd2ndDerivsAtXi<numnodes, numnodalvalues>(
          xi_master, N_i_master, N_i_xi_master, N_i_xixi_master, BeamElement2()->Shape(),
          ele2length_);

      // compute coord vector and tangent vector on master side
      ComputeCenterlinePosition(r_master, N_i_master, ele2pos_);
      ComputeCenterlineTangent(r_xi_master, N_i_xi_master, ele2pos_);
      ComputeCenterlineTangent(r_xixi_master, N_i_xixi_master, ele2pos_);


      T norm_r_xi_master = FADUTILS::VectorNorm(r_xi_master);

      g1_master.Update(1.0 / norm_r_xi_master, r_xi_master);


      //************************** DEBUG ******************************************
      //    std::cout << "\nr_master: " << FADUTILS::CastToDouble<T,3,1>( r_master );
      //    std::cout << "\nr_xi_master: " << FADUTILS::CastToDouble<T,3,1>( r_xi_master );
      //    std::cout << "\n|r_xi_master|: " << FADUTILS::VectorNorm<3>(
      //    FADUTILS::CastToDouble<T,3,1>( r_xi_master ) ); std::cout << "\ng1_master: " <<
      //    FADUTILS::CastToDouble<T,3,1>( g1_master ); std::cout << "\n|g1_master|: " <<
      //    FADUTILS::VectorNorm<3>( FADUTILS::CastToDouble<T,3,1>( g1_master ) );
      //*********************** END DEBUG *****************************************

      // store for vtk visualization
      centerline_coords_GP2_[igp_total] = FADUTILS::CastToDouble<T, 3, 1>(r_master);

      // distance vector between unilateral closest points
      dist_ul = FADUTILS::DiffVector(r_slave, r_master);

      //************************** DEBUG ******************************************
      //    std::cout << "\ndist_ul: " << FADUTILS::CastToDouble<T,3,1>( dist_ul );
      //*********************** END DEBUG *****************************************

      // check cutoff criterion: if specified, contributions are neglected at larger separation
      if (cutoff_radius != -1.0 and FADUTILS::CastToDouble(FADUTILS::Norm(dist_ul)) > cutoff_radius)
        continue;

      // mutual angle of tangent vectors at unilateral closest points
      BEAMINTERACTION::GEO::CalcEnclosedAngle(alpha, cos_alpha, r_xi_slave, r_xi_master);
      sin_alpha = std::sin(alpha);
      sin_2alpha = std::sin(2 * alpha);


      //************************** DEBUG ******************************************
      //    std::cout << "\nalpha: " << FADUTILS::CastToDouble( alpha*180/M_PI ) << "°";
      //    std::cout << "\ncos(alpha): " << FADUTILS::CastToDouble( cos_alpha );
      //    std::cout << "\nsin(alpha): " << FADUTILS::CastToDouble( sin_alpha );
      //*********************** END DEBUG *****************************************


      //************************** DEBUG ******************************************
      if (alpha < 0.0 or alpha > M_PI_2)
        dserror("alpha=%f, should be in [0,pi/2]", FADUTILS::CastToDouble(alpha));
      //*********************** END DEBUG *****************************************

      const double BEAMSCOLINEARANGLETHRESHOLD = 5.0 / 180.0 * M_PI;  // 5 works best so far

      // Fixme what about the FAD linearization in this case?
      if (FADUTILS::Norm(alpha) < BEAMSCOLINEARANGLETHRESHOLD)
      {
        //************************** DEBUG ******************************************
        //      std::cout << "\n\nINFO: Enclosed angle is close to zero: alpha="
        //          << FADUTILS::CastToDouble(alpha)*180/M_PI << "°\n"
        //          << std::endl;
        //*********************** END DEBUG *****************************************

        // there is no unique bilateral closest point pair in case of alpha=0
        // hence, we can use the current (unilateral) closest point pair
        normal_bl = dist_ul;
        norm_normal_bl_tilde = FADUTILS::VectorNorm(normal_bl);
        normal_bl.Scale(1.0 / norm_normal_bl_tilde);

        aux_plane_normal.Clear();
        x = 0.0;
      }
      else
      {
        // normal vector at bilateral closest point Fixme
        normal_bl.CrossProduct(r_xi_slave, r_xi_master);
        norm_normal_bl_tilde = FADUTILS::VectorNorm(normal_bl);
        normal_bl.Scale(1.0 / norm_normal_bl_tilde);

        // distance between Gauss point and bilateral closest point on slave
        aux_plane_normal.Update(r_xi_master.Dot(r_xi_master), r_xi_slave,
            -1.0 * r_xi_master.Dot(r_xi_slave), r_xi_master);

        x = FADUTILS::VectorNorm(r_xi_slave) *
            (r_master.Dot(aux_plane_normal) - r_slave.Dot(aux_plane_normal)) /
            r_xi_slave.Dot(aux_plane_normal);
      }

      // gap of bilateral closest point (also valid for special case alpha=0)
      gap_bl = FADUTILS::Norm(dist_ul.Dot(normal_bl)) - radius1_ - radius2_;

      if (gap_bl < 1e-14)
        dserror("gap=%f is negative or very close to zero! Take care of this case!",
            FADUTILS::CastToDouble(gap_bl));

      //************************** DEBUG ******************************************
      //    std::cout << "\nnormal_bl: " << FADUTILS::CastToDouble<T,3,1>( normal_bl );
      //    std::cout << "\ngap_bl: " << FADUTILS::CastToDouble( gap_bl );
      //    std::cout << "\naux_plane_normal: " << FADUTILS::CastToDouble<T,3,1>( aux_plane_normal
      //    ); std::cout << "\nx: " << FADUTILS::CastToDouble( x ) << std::endl;
      //*********************** END DEBUG *****************************************

      beta = FADUTILS::sqrt<T>(
          (gap_bl + radius2_) * (gap_bl + radius2_) + x * x * sin_alpha * sin_alpha);
      beta_exp2 = beta * beta;
      beta_exp3 = beta_exp2 * beta;
      beta_exp4 = beta_exp2 * beta_exp2;

      a = 0.5 / beta *
          ((gap_bl + radius2_) / radius2_ + cos_alpha * cos_alpha -
              x * x * sin_2alpha * sin_2alpha / (4.0 * beta_exp2));

      Delta = 4 * a * (beta - radius2_) - x * x * sin_2alpha * sin_2alpha / (4 * beta_exp2);

      //************************** DEBUG ******************************************
      //    std::cout << "\nbeta: " << FADUTILS::CastToDouble( beta );
      //    std::cout << "\nbeta^2: " << FADUTILS::CastToDouble( beta*beta );
      //    std::cout << "\nbeta^3: " << FADUTILS::CastToDouble( beta*beta*beta );
      //    std::cout << "\nDelta: " << FADUTILS::CastToDouble( Delta );
      //*********************** END DEBUG *****************************************


      // interaction potential
      interaction_potential_GP =
          beta / (gap_bl + radius2_) * std::pow(a, m_ - 5) * std::pow(Delta, -m_ + 4.5);


      // partial derivatives of beta
      beta_partial_x = x * sin_alpha * sin_alpha / beta;
      beta_partial_gap_bl = (gap_bl + radius2_) / beta;
      beta_partial_cos_alpha = -x * x * cos_alpha / beta;

      // partial derivatives of a
      a_partial_beta = -a / beta + x * x * sin_2alpha * sin_2alpha / (4 * beta_exp4);

      a_partial_x = a_partial_beta * beta_partial_x - x * sin_2alpha * sin_2alpha / (4 * beta_exp3);

      a_partial_gap_bl = a_partial_beta * beta_partial_gap_bl + 0.5 / (radius2_ * beta);

      a_partial_cos_alpha = a_partial_beta * beta_partial_cos_alpha + cos_alpha / beta -
                            x * x * cos_alpha * (1 - 2 * cos_alpha * cos_alpha) / beta_exp3;

      // partial derivatives of Delta
      Delta_partial_beta =
          2.0 * (((-(gap_bl + radius2_) / radius2_ - cos_alpha * cos_alpha) / beta_exp2 +
                     3.0 * x * x * sin_2alpha * sin_2alpha / 4.0 / beta_exp4) *
                        (beta - radius2_) +
                    (gap_bl + radius2_) / beta / radius2_ + cos_alpha * cos_alpha / beta);

      Delta_partial_x = Delta_partial_beta * beta_partial_x -
                        x * sin_2alpha * sin_2alpha / beta_exp3 * (1.5 * beta - radius2_);

      Delta_partial_gap_bl =
          Delta_partial_beta * beta_partial_gap_bl + 2 / beta / radius2_ * (beta - radius2_);

      Delta_partial_cos_alpha =
          Delta_partial_beta * beta_partial_cos_alpha +
          4.0 *
              (cos_alpha / beta -
                  x * x * (cos_alpha - 2.0 * cos_alpha * cos_alpha * cos_alpha) / beta_exp3) *
              (beta - radius2_) -
          2.0 * x * x * (cos_alpha - 2.0 * cos_alpha * cos_alpha * cos_alpha) / beta_exp2;


      //************************** DEBUG ******************************************
      //    std::cout << "\nDelta_partial_beta: " << FADUTILS::CastToDouble( Delta_partial_beta );
      //    std::cout << "\nbeta_partial_cos_alpha: " << FADUTILS::CastToDouble(
      //    beta_partial_cos_alpha ); std::cout << "\nDelta_partial_cos_alpha: " <<
      //    FADUTILS::CastToDouble( Delta_partial_cos_alpha );
      //*********************** END DEBUG *****************************************


      // partial derivatives of single length specific interaction potential

      pot_ia_partial_beta = interaction_potential_GP / beta;

      pot_ia_partial_a = (m_ - 5.0) * interaction_potential_GP / a;

      pot_ia_partial_Delta = (-m_ + 4.5) * interaction_potential_GP / Delta;


      pot_ia_partial_x =
          pot_ia_partial_beta * beta_partial_x + pot_ia_partial_a * a_partial_x +
          pot_ia_partial_Delta * Delta_partial_x;  // + interaction_potential_GP / beta;

      pot_ia_partial_gap_bl = pot_ia_partial_beta * beta_partial_gap_bl +
                              pot_ia_partial_a * a_partial_gap_bl +
                              pot_ia_partial_Delta * Delta_partial_gap_bl -
                              interaction_potential_GP / (gap_bl + radius2_);

      pot_ia_partial_cos_alpha = pot_ia_partial_beta * beta_partial_cos_alpha +
                                 pot_ia_partial_a * a_partial_cos_alpha +
                                 pot_ia_partial_Delta * Delta_partial_cos_alpha;

      //************************** DEBUG ******************************************
      //    std::cout << "\npot_ia_partial_Delta: " << FADUTILS::CastToDouble( pot_ia_partial_Delta
      //    ); std::cout << "\npot_ia_partial_cos_alpha: " << FADUTILS::CastToDouble(
      //    pot_ia_partial_cos_alpha );
      //*********************** END DEBUG *****************************************


      // components from variation of bilateral gap
      T signum_dist_bl_tilde = FADUTILS::Signum(dist_ul.Dot(normal_bl));

      gap_bl_partial_r_slave.Update(signum_dist_bl_tilde, normal_bl);
      gap_bl_partial_r_master.Update(-1.0 * signum_dist_bl_tilde, normal_bl);

      // Todo: check and remove: the following term should vanish since r_xi_master \perp normal_bl
      gap_bl_partial_xi_master = -1.0 * signum_dist_bl_tilde * r_xi_master.Dot(normal_bl);

      v_mat_tmp.Clear();
      for (unsigned int i = 0; i < 3; ++i)
      {
        v_mat_tmp(i, i) += 1.0;
        for (unsigned int j = 0; j < 3; ++j) v_mat_tmp(i, j) -= normal_bl(i) * normal_bl(j);
      }

      LINALG::TMatrix<T, 3, 1> vec_tmp(true);
      vec_tmp.Multiply(v_mat_tmp, dist_ul);


      if (alpha < BEAMSCOLINEARANGLETHRESHOLD)
      {
        // Todo: check and remove: signum_dist_bl_tilde should always be +1.0 in this case!

        gap_bl_partial_r_xi_slave.Clear();

        gap_bl_partial_r_slave.Update(signum_dist_bl_tilde / norm_normal_bl_tilde, vec_tmp, 1.0);

        gap_bl_partial_r_xi_master.Clear();

        gap_bl_partial_r_master.Update(
            -1.0 * signum_dist_bl_tilde / norm_normal_bl_tilde, vec_tmp, 1.0);

        // Todo: check and remove: the following term should vanish
        gap_bl_partial_xi_master -=
            signum_dist_bl_tilde / norm_normal_bl_tilde * r_xi_master.Dot(vec_tmp);
      }
      else
      {
        gap_bl_partial_r_xi_slave.Update(signum_dist_bl_tilde / norm_normal_bl_tilde,
            FADUTILS::VectorProduct(r_xi_master, vec_tmp));

        gap_bl_partial_r_xi_master.Update(-1.0 * signum_dist_bl_tilde / norm_normal_bl_tilde,
            FADUTILS::VectorProduct(r_xi_slave, vec_tmp));

        gap_bl_partial_xi_master += r_xixi_master.Dot(gap_bl_partial_r_xi_master);
      }


      // components from variation of cosine of enclosed angle
      T signum_tangentsscalarproduct = FADUTILS::Signum(r_xi_slave.Dot(r_xi_master));

      v_mat_tmp.Clear();
      for (unsigned int i = 0; i < 3; ++i)
      {
        v_mat_tmp(i, i) += 1.0;
        for (unsigned int j = 0; j < 3; ++j) v_mat_tmp(i, j) -= g1_slave(i) * g1_slave(j);
      }

      cos_alpha_partial_r_xi_slave.Multiply(
          signum_tangentsscalarproduct / FADUTILS::VectorNorm(r_xi_slave), v_mat_tmp, g1_master);

      v_mat_tmp.Clear();
      for (unsigned int i = 0; i < 3; ++i)
      {
        v_mat_tmp(i, i) += 1.0;
        for (unsigned int j = 0; j < 3; ++j) v_mat_tmp(i, j) -= g1_master(i) * g1_master(j);
      }

      cos_alpha_partial_r_xi_master.Multiply(
          signum_tangentsscalarproduct / FADUTILS::VectorNorm(r_xi_master), v_mat_tmp, g1_slave);

      cos_alpha_partial_xi_master = r_xixi_master.Dot(cos_alpha_partial_r_xi_master);


      // components from variation of distance from bilateral closest point on slave
      if (FADUTILS::CastToDouble(alpha) < BEAMSCOLINEARANGLETHRESHOLD)
      {
        /* set the following quantities to zero since they are not required in this case
         * because pot_ia_partial_x is zero anyway */
        x_partial_r_slave.Clear();
        x_partial_r_master.Clear();
        x_partial_r_xi_slave.Clear();
        x_partial_r_xi_master.Clear();
        x_partial_xi_master = 0.0;
      }
      else
      {
        x_partial_r_slave.Update(-1.0 / g1_slave.Dot(aux_plane_normal), aux_plane_normal);
        x_partial_r_master.Update(1.0 / g1_slave.Dot(aux_plane_normal), aux_plane_normal);

        v_mat_tmp.Clear();
        for (unsigned int i = 0; i < 3; ++i)
        {
          v_mat_tmp(i, i) += 1.0;
          for (unsigned int j = 0; j < 3; ++j) v_mat_tmp(i, j) -= g1_slave(i) * g1_slave(j);
        }

        x_partial_r_xi_slave.Multiply(1.0 / FADUTILS::VectorNorm(r_xi_slave) *
                                          dist_ul.Dot(aux_plane_normal) /
                                          std::pow(g1_slave.Dot(aux_plane_normal), 2),
            v_mat_tmp, aux_plane_normal);

        x_partial_aux_plane_normal.Update(-1.0 / g1_slave.Dot(aux_plane_normal), dist_ul,
            dist_ul.Dot(aux_plane_normal) / std::pow(g1_slave.Dot(aux_plane_normal), 2), g1_slave);

        v_mat_tmp.Clear();
        for (unsigned int i = 0; i < 3; ++i)
          for (unsigned int j = 0; j < 3; ++j) v_mat_tmp(i, j) += r_xi_master(i) * r_xi_master(j);

        x_partial_r_xi_slave.Update(r_xi_master.Dot(r_xi_master), x_partial_aux_plane_normal, 1.0);

        x_partial_r_xi_slave.Multiply(-1.0, v_mat_tmp, x_partial_aux_plane_normal, 1.0);


        x_partial_r_xi_master.Update(
            -1.0 * r_xi_master.Dot(r_xi_slave), x_partial_aux_plane_normal);

        v_mat_tmp.Clear();
        for (unsigned int i = 0; i < 3; ++i)
          for (unsigned int j = 0; j < 3; ++j)
          {
            v_mat_tmp(i, j) += 2.0 * r_xi_master(i) * r_xi_slave(j);
            v_mat_tmp(i, j) -= 1.0 * r_xi_slave(i) * r_xi_master(j);
          }

        x_partial_r_xi_master.Multiply(1.0, v_mat_tmp, x_partial_aux_plane_normal, 1.0);

        x_partial_xi_master = r_xixi_master.Dot(x_partial_r_xi_master);
      }


      // Todo: maybe avoid this assembly of shape fcns into matrices
      DRT::UTILS::BEAM::AssembleShapeFunctions<numnodes, numnodalvalues>(N_i_slave[igp], N_slave);

      DRT::UTILS::BEAM::AssembleShapeFunctionsAndDerivsAnd2ndDerivs<numnodes, numnodalvalues>(
          N_i_master, N_i_xi_master, N_i_xixi_master, N_master, N_xi_master, N_xixi_master);

      BEAMINTERACTION::GEO::CalcLinearizationPointToCurveProjectionParameterCoordMaster<numnodes,
          numnodalvalues>(lin_xi_master_slaveDofs, lin_xi_master_masterDofs, dist_ul, r_xi_master,
          r_xixi_master, N_slave, N_master, N_xixi_master);


      BEAMINTERACTION::GEO::CalcPointToCurveProjectionParameterCoordMasterPartialDerivs(
          xi_master_partial_r_slave, xi_master_partial_r_master, xi_master_partial_r_xi_master,
          dist_ul, r_xi_master, r_xixi_master);

      // store for vtk visualization
      forces_pot_GP1_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_gap_bl),
          FADUTILS::CastToDouble<T, 3, 1>(gap_bl_partial_r_slave), 1.0);

      forces_pot_GP1_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_gap_bl * gap_bl_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_slave), 1.0);

      forces_pot_GP1_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_cos_alpha * cos_alpha_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_slave), 1.0);

      forces_pot_GP1_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_x),
          FADUTILS::CastToDouble<T, 3, 1>(x_partial_r_slave), 1.0);

      forces_pot_GP1_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_x * x_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_slave), 1.0);

      forces_pot_GP1_[igp_total].Scale(-1.0 * prefactor * rho1 * rho2);


      moments_pot_GP1_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_gap_bl),
          FADUTILS::CastToDouble<T, 3, 1>(gap_bl_partial_r_xi_slave), 1.0);

      moments_pot_GP1_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_cos_alpha),
          FADUTILS::CastToDouble<T, 3, 1>(cos_alpha_partial_r_xi_slave), 1.0);

      moments_pot_GP1_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_x),
          FADUTILS::CastToDouble<T, 3, 1>(x_partial_r_xi_slave), 1.0);

      /* note: relation between variation of tangent vector r_xi and variation of (transversal part
       *       of) rotation vector theta_perp describing cross-section orientation can be used to
       *       identify (distributed) moments as follows: m_pot = 1/|r_xi| * ( m_pot_pseudo x g1 )
       */
      LINALG::TMatrix<double, 3, 3> spin_pseudo_moment_tmp(true);
      LARGEROTATIONS::computespin(spin_pseudo_moment_tmp, moments_pot_GP1_[igp_total]);

      moments_pot_GP1_[igp_total].Multiply(1.0 / FADUTILS::CastToDouble(norm_r_xi_slave),
          spin_pseudo_moment_tmp, FADUTILS::CastToDouble<T, 3, 1>(g1_slave));

      moments_pot_GP1_[igp_total].Scale(-1.0 * prefactor * rho1 * rho2);


      forces_pot_GP2_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_gap_bl),
          FADUTILS::CastToDouble<T, 3, 1>(gap_bl_partial_r_master), 1.0);

      forces_pot_GP2_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_gap_bl * gap_bl_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_master), 1.0);

      forces_pot_GP2_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_cos_alpha * cos_alpha_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_master), 1.0);

      forces_pot_GP2_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_x),
          FADUTILS::CastToDouble<T, 3, 1>(x_partial_r_master), 1.0);

      forces_pot_GP2_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_x * x_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_master), 1.0);

      forces_pot_GP2_[igp_total].Scale(-1.0 * prefactor * rho1 * rho2);


      moments_pot_GP2_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_gap_bl),
          FADUTILS::CastToDouble<T, 3, 1>(gap_bl_partial_r_xi_master), 1.0);

      moments_pot_GP2_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_gap_bl * gap_bl_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_xi_master), 1.0);

      moments_pot_GP2_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_cos_alpha),
          FADUTILS::CastToDouble<T, 3, 1>(cos_alpha_partial_r_xi_master), 1.0);

      moments_pot_GP2_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_cos_alpha * cos_alpha_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_xi_master), 1.0);

      moments_pot_GP2_[igp_total].Update(FADUTILS::CastToDouble(pot_ia_partial_x),
          FADUTILS::CastToDouble<T, 3, 1>(x_partial_r_xi_master), 1.0);

      moments_pot_GP2_[igp_total].UpdateT(
          FADUTILS::CastToDouble(pot_ia_partial_x * x_partial_xi_master),
          FADUTILS::CastToDouble<T, 1, 3>(xi_master_partial_r_xi_master), 1.0);

      LARGEROTATIONS::computespin(spin_pseudo_moment_tmp, moments_pot_GP2_[igp_total]);

      moments_pot_GP2_[igp_total].Multiply(1.0 / FADUTILS::CastToDouble(norm_r_xi_master),
          spin_pseudo_moment_tmp, FADUTILS::CastToDouble<T, 3, 1>(g1_master));

      moments_pot_GP2_[igp_total].Scale(-1.0 * prefactor * rho1 * rho2);


      // integration factor
      pot_ia_partial_gap_bl *= rho1rho2_JacFac_GaussWeight;
      pot_ia_partial_cos_alpha *= rho1rho2_JacFac_GaussWeight;
      pot_ia_partial_x *= rho1rho2_JacFac_GaussWeight;

      //************************** DEBUG ******************************************
      //    std::cout << "\nprefactor: " << FADUTILS::CastToDouble( prefactor );
      //    std::cout << "\nrho1rho2_JacFac_GaussWeight: " << FADUTILS::CastToDouble(
      //    rho1rho2_JacFac_GaussWeight ); std::cout << "\npot_ia_partial_cos_alpha: " <<
      //    FADUTILS::CastToDouble( pot_ia_partial_cos_alpha );
      //*********************** END DEBUG *****************************************


      //********************************************************************
      // calculate fpot1: force/residual vector on slave element
      //********************************************************************
      force_pot_slave_GP.Clear();
      // sum up the contributions of all nodes (in all dimensions)
      for (unsigned int idofperdim = 0; idofperdim < (numnodes * numnodalvalues); ++idofperdim)
      {
        // loop over dimensions
        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          force_pot_slave_GP(3 * idofperdim + jdim) +=
              N_i_slave[igp](idofperdim) * pot_ia_partial_gap_bl *
              (gap_bl_partial_r_slave(jdim) +
                  gap_bl_partial_xi_master * xi_master_partial_r_slave(jdim));

          force_pot_slave_GP(3 * idofperdim + jdim) += N_i_xi_slave[igp](idofperdim) *
                                                       pot_ia_partial_gap_bl *
                                                       gap_bl_partial_r_xi_slave(jdim);


          force_pot_slave_GP(3 * idofperdim + jdim) +=
              N_i_slave[igp](idofperdim) * pot_ia_partial_cos_alpha * cos_alpha_partial_xi_master *
              xi_master_partial_r_slave(jdim);

          force_pot_slave_GP(3 * idofperdim + jdim) += N_i_xi_slave[igp](idofperdim) *
                                                       pot_ia_partial_cos_alpha *
                                                       cos_alpha_partial_r_xi_slave(jdim);


          force_pot_slave_GP(3 * idofperdim + jdim) +=
              N_i_slave[igp](idofperdim) * pot_ia_partial_x *
              (x_partial_r_slave(jdim) + x_partial_xi_master * xi_master_partial_r_slave(jdim));

          force_pot_slave_GP(3 * idofperdim + jdim) +=
              N_i_xi_slave[igp](idofperdim) * pot_ia_partial_x * x_partial_r_xi_slave(jdim);
        }
      }

      //********************************************************************
      // calculate fpot2: force/residual vector on master element
      //********************************************************************
      force_pot_master_GP.Clear();
      // sum up the contributions of all nodes (in all dimensions)
      for (unsigned int idofperdim = 0; idofperdim < (numnodes * numnodalvalues); ++idofperdim)
      {
        // loop over dimensions
        for (unsigned int jdim = 0; jdim < 3; ++jdim)
        {
          force_pot_master_GP(3 * idofperdim + jdim) +=
              N_i_master(idofperdim) * pot_ia_partial_gap_bl *
              (gap_bl_partial_r_master(jdim) +
                  gap_bl_partial_xi_master * xi_master_partial_r_master(jdim));

          force_pot_master_GP(3 * idofperdim + jdim) +=
              N_i_xi_master(idofperdim) * pot_ia_partial_gap_bl *
              (gap_bl_partial_r_xi_master(jdim) +
                  gap_bl_partial_xi_master * xi_master_partial_r_xi_master(jdim));


          force_pot_master_GP(3 * idofperdim + jdim) +=
              N_i_master(idofperdim) * pot_ia_partial_cos_alpha * cos_alpha_partial_xi_master *
              xi_master_partial_r_master(jdim);

          force_pot_master_GP(3 * idofperdim + jdim) +=
              N_i_xi_master(idofperdim) * pot_ia_partial_cos_alpha *
              (cos_alpha_partial_r_xi_master(jdim) +
                  cos_alpha_partial_xi_master * xi_master_partial_r_xi_master(jdim));


          force_pot_master_GP(3 * idofperdim + jdim) +=
              N_i_master(idofperdim) * pot_ia_partial_x *
              (x_partial_r_master(jdim) + x_partial_xi_master * xi_master_partial_r_master(jdim));

          force_pot_master_GP(3 * idofperdim + jdim) +=
              N_i_xi_master(idofperdim) * pot_ia_partial_x *
              (x_partial_r_xi_master(jdim) +
                  x_partial_xi_master * xi_master_partial_r_xi_master(jdim));
        }
      }

      // apply constant prefactor
      force_pot_slave_GP.Scale(prefactor);
      force_pot_master_GP.Scale(prefactor);

      // sum contributions from all Gauss points
      force_pot1.Update(1.0, force_pot_slave_GP, 1.0);
      force_pot2.Update(1.0, force_pot_master_GP, 1.0);

      if (stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL)
      {
        AddStiffmatContributionsXiMasterAutomaticDifferentiationIfRequired(force_pot_slave_GP,
            force_pot_master_GP, lin_xi_master_slaveDofs, lin_xi_master_masterDofs, *stiffmat11,
            *stiffmat12, *stiffmat21, *stiffmat22);
      }

      // do this scaling down here, because the value without the prefactors is meant to be used
      // in the force calculation above!
      interaction_potential_GP *= rho1rho2_JacFac_GaussWeight * prefactor;

      // store for energy output
      interaction_potential_ += FADUTILS::CastToDouble(interaction_potential_GP);

      // ************************************ DEBUG
      // *************************************************
      CalcFpotGausspointAutomaticDifferentiationIfRequired(force_pot_slave_GP_calc_via_FAD,
          force_pot_master_GP_calc_via_FAD, interaction_potential_GP, lin_xi_master_slaveDofs,
          lin_xi_master_masterDofs);

      // hard-set values below double precision to zero to ease comparison
      //    for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i)
      //    {
      //      if ( FADUTILS::Norm( force_pot_slave_GP(i) ) < 1e-15 )
      //        force_pot_slave_GP(i) = 0.0;
      //      if ( FADUTILS::Norm( force_pot_master_GP(i) ) < 1e-15 )
      //        force_pot_master_GP(i) = 0.0;
      //
      //      if ( FADUTILS::Norm( force_pot_slave_GP_calc_via_FAD(i) ) < 1e-15 )
      //        force_pot_slave_GP_calc_via_FAD(i) = 0.0;
      //      if ( FADUTILS::Norm( force_pot_master_GP_calc_via_FAD(i) ) < 1e-15 )
      //        force_pot_master_GP_calc_via_FAD(i) = 0.0;
      //    }

      // print FAD based residuum vectors
      //    std::cout << std::scientific << std::setprecision(10)
      //        << "\n\nforce_pot_slave_GP: " << force_pot_slave_GP_calc_via_FAD << std::endl;
      //    std::cout << std::scientific << std::setprecision(10)
      //        << "\n\nforce_pot_master_GP: " << force_pot_master_GP_calc_via_FAD << std::endl;

      // print analytic residuum vectors
      //    std::cout << "\n\nforce_pot_slave_GP: "
      //        << std::scientific << std::setprecision(10)
      //        << FADUTILS::CastToDouble<T,3*numnodes*numnodalvalues,1>( force_pot_slave_GP ) <<
      //        std::endl;
      //    std::cout << "\n\nforce_pot_master_GP: "
      //        << std::scientific << std::setprecision(10)
      //        << FADUTILS::CastToDouble<T,3*numnodes*numnodalvalues,1>( force_pot_master_GP ) <<
      //        std::endl;



      // comparison of residuum vectors and print results

      const double ABSDIFFTOL = 1e-12;
      const double RELDIFFTOL = 1e-10;

      for (unsigned int i = 0; i < force_pot_slave_GP.M(); ++i)
      {
        if (std::abs(FADUTILS::CastToDouble(force_pot_slave_GP(i)) -
                     force_pot_slave_GP_calc_via_FAD(i)) > ABSDIFFTOL and
            (std::abs(force_pot_slave_GP_calc_via_FAD(i)) > 1e-15 and
                std::abs((FADUTILS::CastToDouble(force_pot_slave_GP(i)) -
                             force_pot_slave_GP_calc_via_FAD(i)) /
                         force_pot_slave_GP_calc_via_FAD(i)) > RELDIFFTOL))
        {
          std::cout << "\n\ndetected difference in force_pot_slave_GP for i=" << i
                    << ": abs_diff=" << std::scientific << std::setprecision(10)
                    << FADUTILS::CastToDouble(force_pot_slave_GP(i)) -
                           force_pot_slave_GP_calc_via_FAD(i)
                    << ", rel_diff="
                    << (FADUTILS::CastToDouble(force_pot_slave_GP(i)) -
                           force_pot_slave_GP_calc_via_FAD(i)) /
                           force_pot_slave_GP_calc_via_FAD(i)
                    << std::endl;
        }
      }



      //    std::cout << "\n\nforce_pot_slave_GP (dbl): "
      //        << std::scientific << std::setprecision(10)
      //        << FADUTILS::CastToDouble<T,3*numnodes*numnodalvalues,1>( force_pot_slave_GP ) <<
      //        std::endl;
      //
      //    std::cout << "force_pot_slave_GP (FAD): "
      //        << std::scientific << std::setprecision(10)
      //        << force_pot_slave_GP_calc_via_FAD << std::endl;
      //
      //
      //    std::cout << "force_pot_slave_GP (ratio dbl/FAD): ";
      //
      //    for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i)
      //    {
      //      if ( std::abs( FADUTILS::CastToDouble( force_pot_slave_GP(i) )) > 1e-15)
      //        std::cout << std::scientific << std::setprecision(10)
      //                  << force_pot_slave_GP_calc_via_FAD(i)/FADUTILS::CastToDouble(
      //                  force_pot_slave_GP(i) )<< " ";
      //      else
      //        std::cout << "nan ";
      //    }
      //
      //    std::cout << std::endl;


      for (unsigned int i = 0; i < force_pot_master_GP.M(); ++i)
      {
        if (std::abs(FADUTILS::CastToDouble(force_pot_master_GP(i)) -
                     force_pot_master_GP_calc_via_FAD(i)) > ABSDIFFTOL and
            (std::abs(force_pot_master_GP_calc_via_FAD(i)) > 1e-15 and
                std::abs((FADUTILS::CastToDouble(force_pot_master_GP(i)) -
                             force_pot_master_GP_calc_via_FAD(i)) /
                         force_pot_master_GP_calc_via_FAD(i)) > RELDIFFTOL))
        {
          std::cout << "\n\ndetected difference in force_pot_master_GP for i=" << i
                    << ": abs_diff=" << std::scientific << std::setprecision(10)
                    << FADUTILS::CastToDouble(force_pot_master_GP(i)) -
                           force_pot_master_GP_calc_via_FAD(i)
                    << ", rel_diff="
                    << (FADUTILS::CastToDouble(force_pot_master_GP(i)) -
                           force_pot_master_GP_calc_via_FAD(i)) /
                           force_pot_master_GP_calc_via_FAD(i)
                    << std::endl;
        }
      }


      //    std::cout << "\n\nforce_pot_master_GP (dbl): "
      //        << std::scientific << std::setprecision(10)
      //        << FADUTILS::CastToDouble<T,3*numnodes*numnodalvalues,1>( force_pot_master_GP ) <<
      //        std::endl;
      //
      //    std::cout << "force_pot_master_GP (FAD): "
      //        << std::scientific << std::setprecision(10)
      //        << force_pot_master_GP_calc_via_FAD << std::endl;
      //
      //
      //    std::cout << "force_pot_master_GP (ratio dbl/FAD): ";
      //
      //
      //    for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i)
      //    {
      //      if (std::abs( FADUTILS::CastToDouble( force_pot_master_GP(i)) ) > 1e-15)
      //        std::cout << std::scientific << std::setprecision(10)
      //                  << force_pot_master_GP_calc_via_FAD(i)/FADUTILS::CastToDouble(
      //                  force_pot_master_GP(i) ) << " ";
      //      else
      //        std::cout << "nan ";
      //    }
      //
      //    std::cout << std::endl;
      // ******************************* END DEBUG **********************************************

    }  // end loop over Gauss points per segment
  }    // end loop over integration segments

  if (stiffmat11 != NULL and stiffmat12 != NULL and stiffmat21 != NULL and stiffmat22 != NULL)
  {
    ScaleStiffpotAnalyticContributionsIfRequired(
        prefactor, *stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);


    CalcStiffmatAutomaticDifferentiationIfRequired(
        force_pot1, force_pot2, *stiffmat11, *stiffmat12, *stiffmat21, *stiffmat22);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues,
    T>::ScaleStiffpotAnalyticContributionsIfRequired(double const& scalefactor,
    LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
{
  stiffmat11.Scale(scalefactor);
  stiffmat12.Scale(scalefactor);
  stiffmat21.Scale(scalefactor);
  stiffmat22.Scale(scalefactor);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues,
    T>::CalcStiffmatAutomaticDifferentiationIfRequired(LINALG::TMatrix<Sacado::Fad::DFad<double>,
                                                           3 * numnodes * numnodalvalues, 1> const&
                                                           force_pot1,
    LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1> const& force_pot2,
    LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
    LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
{
  for (unsigned int idof = 0; idof < 3 * numnodes * numnodalvalues; ++idof)
  {
    for (unsigned int jdof = 0; jdof < 3 * numnodes * numnodalvalues; ++jdof)
    {
      // d (Res_1) / d (d_1)
      stiffmat11(idof, jdof) += force_pot1(idof).dx(jdof);

      // d (Res_1) / d (d_2)
      stiffmat12(idof, jdof) += force_pot1(idof).dx(3 * numnodes * numnodalvalues + jdof);

      // d (Res_2) / d (d_1)
      stiffmat21(idof, jdof) += force_pot2(idof).dx(jdof);

      // d (Res_2) / d (d_2)
      stiffmat22(idof, jdof) += force_pot2(idof).dx(3 * numnodes * numnodalvalues + jdof);
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
    AddStiffmatContributionsXiMasterAutomaticDifferentiationIfRequired(
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1> const&
            force_pot1,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1> const&
            force_pot2,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_slaveDofs,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_masterDofs,
        LINALG::SerialDenseMatrix& stiffmat11, LINALG::SerialDenseMatrix& stiffmat12,
        LINALG::SerialDenseMatrix& stiffmat21, LINALG::SerialDenseMatrix& stiffmat22) const
{
  const unsigned int dim = 3 * numnodes * numnodalvalues;

  for (unsigned int idof = 0; idof < dim; ++idof)
  {
    for (unsigned int jdof = 0; jdof < dim; ++jdof)
    {
      // d (Res_1) / d (d_1)
      stiffmat11(idof, jdof) +=
          force_pot1(idof).dx(2 * dim) * FADUTILS::CastToDouble(lin_xi_master_slaveDofs(jdof));

      // d (Res_1) / d (d_2)
      stiffmat12(idof, jdof) +=
          force_pot1(idof).dx(2 * dim) * FADUTILS::CastToDouble(lin_xi_master_masterDofs(jdof));

      // d (Res_2) / d (d_1)
      stiffmat21(idof, jdof) +=
          force_pot2(idof).dx(2 * dim) * FADUTILS::CastToDouble(lin_xi_master_slaveDofs(jdof));

      // d (Res_2) / d (d_2)
      stiffmat22(idof, jdof) +=
          force_pot2(idof).dx(2 * dim) * FADUTILS::CastToDouble(lin_xi_master_masterDofs(jdof));
    }
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::
    CalcFpotGausspointAutomaticDifferentiationIfRequired(
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& force_pot1,
        LINALG::TMatrix<double, 3 * numnodes * numnodalvalues, 1>& force_pot2,
        Sacado::Fad::DFad<double> const& interaction_potential,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_slaveDofs,
        LINALG::TMatrix<Sacado::Fad::DFad<double>, 1, 3 * numnodes * numnodalvalues> const&
            lin_xi_master_masterDofs) const
{
  const unsigned int dim = 3 * numnodalvalues * numnodes;

  for (unsigned int idof = 0; idof < dim; ++idof)
  {
    force_pot1(idof) =
        interaction_potential.dx(idof) +
        interaction_potential.dx(2 * dim) * FADUTILS::CastToDouble(lin_xi_master_slaveDofs(idof));

    force_pot2(idof) =
        interaction_potential.dx(dim + idof) +
        interaction_potential.dx(2 * dim) * FADUTILS::CastToDouble(lin_xi_master_masterDofs(idof));
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::Print(
    std::ostream& out) const
{
  CheckInitSetup();

  out << "\nInstance of BeamToBeamPotentialPair (EleGIDs " << Element1()->Id() << " & "
      << Element2()->Id() << "):";
  out << "\nele1 dofvec: " << FADUTILS::CastToDouble<T, 3 * numnodes * numnodalvalues, 1>(ele1pos_);
  out << "\nele2 dofvec: " << FADUTILS::CastToDouble<T, 3 * numnodes * numnodalvalues, 1>(ele2pos_);

  out << "\n";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues,
    T>::PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const
{
  CheckInitSetup();

  // Todo difficulty here is that the same element pair is evaluated more than once
  //      to be more precise, once for every common potlaw;
  //      contribution of previous evaluations is overwritten if multiple potlaws are applied
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::GetShapeFunctions(
    std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>>& N1_i,
    std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>>& N2_i,
    std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>>& N1_i_xi,
    std::vector<LINALG::TMatrix<double, 1, numnodes * numnodalvalues>>& N2_i_xi,
    DRT::UTILS::IntegrationPoints1D& gausspoints) const
{
  DRT::UTILS::BEAM::EvaluateShapeFunctionsAndDerivsAllGPs<numnodes, numnodalvalues>(
      gausspoints, N1_i, N1_i_xi, BeamElement1()->Shape(), ele1length_);

  DRT::UTILS::BEAM::EvaluateShapeFunctionsAndDerivsAllGPs<numnodes, numnodalvalues>(
      gausspoints, N2_i, N2_i_xi, BeamElement2()->Shape(), ele2length_);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/

template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
template <typename T2>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues,
    T>::ComputeCenterlinePosition(LINALG::TMatrix<T, 3, 1>& r,
    const LINALG::TMatrix<T2, 1, numnodes * numnodalvalues>& N_i,
    const LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> eledofvec) const
{
  DRT::UTILS::BEAM::CalcInterpolation<numnodes, numnodalvalues, 3, T, T2>(eledofvec, N_i, r);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
template <typename T2>
void
    BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::ComputeCenterlineTangent(
        LINALG::TMatrix<T, 3, 1>& r_xi,
        const LINALG::TMatrix<T2, 1, numnodes * numnodalvalues>& N_i_xi,
        const LINALG::TMatrix<T, 3 * numnodes * numnodalvalues, 1> eledofvec) const
{
  DRT::UTILS::BEAM::CalcInterpolation<numnodes, numnodalvalues, 3, T, T2>(eledofvec, N_i_xi, r_xi);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues, T>::ResetState(double time,
    const std::vector<double>& centerline_dofvec_ele1,
    const std::vector<double>& centerline_dofvec_ele2)
{
  time_ = time;

  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
  {
    ele1pos_(i) = centerline_dofvec_ele1[i];
    ele2pos_(i) = centerline_dofvec_ele2[i];
  }


  forces_pot_GP1_.clear();
  forces_pot_GP2_.clear();
  moments_pot_GP1_.clear();
  moments_pot_GP2_.clear();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues,
    T>::SetAutomaticDifferentiationVariablesIfRequired(LINALG::TMatrix<Sacado::Fad::DFad<double>,
                                                           3 * numnodes * numnodalvalues, 1>&
                                                           ele1centerlinedofvec,
    LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1>&
        ele2centerlinedofvec)
{
  // The 2*3*numnodes*numnodalvalues primary DoFs consist of all nodal positions and tangents
  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
    ele1centerlinedofvec(i).diff(i, 2 * 3 * numnodes * numnodalvalues);

  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
    ele2centerlinedofvec(i).diff(
        3 * numnodes * numnodalvalues + i, 2 * 3 * numnodes * numnodalvalues);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <unsigned int numnodes, unsigned int numnodalvalues, typename T>
void BEAMINTERACTION::BeamToBeamPotentialPair<numnodes, numnodalvalues,
    T>::SetAutomaticDifferentiationVariablesIfRequired(LINALG::TMatrix<Sacado::Fad::DFad<double>,
                                                           3 * numnodes * numnodalvalues, 1>&
                                                           ele1centerlinedofvec,
    LINALG::TMatrix<Sacado::Fad::DFad<double>, 3 * numnodes * numnodalvalues, 1>&
        ele2centerlinedofvec,
    Sacado::Fad::DFad<double>& xi_master)
{
  // The 2*3*numnodes*numnodalvalues primary DoFs consist of all nodal positions and tangents
  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
    ele1centerlinedofvec(i).diff(i, 2 * 3 * numnodes * numnodalvalues + 1);

  for (unsigned int i = 0; i < 3 * numnodes * numnodalvalues; ++i)
    ele2centerlinedofvec(i).diff(
        3 * numnodes * numnodalvalues + i, 2 * 3 * numnodes * numnodalvalues + 1);

  // Additionally, we set the parameter coordinate on the master side xi_master as primary Dof
  xi_master.diff(2 * 3 * numnodes * numnodalvalues, 2 * 3 * numnodes * numnodalvalues + 1);
}

// explicit template instantiations
template class BEAMINTERACTION::BeamToBeamPotentialPair<2, 1, double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<2, 1, Sacado::Fad::DFad<double>>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<3, 1, double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<3, 1, Sacado::Fad::DFad<double>>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<4, 1, double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<4, 1, Sacado::Fad::DFad<double>>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<5, 1, double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<5, 1, Sacado::Fad::DFad<double>>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<2, 2, double>;
template class BEAMINTERACTION::BeamToBeamPotentialPair<2, 2, Sacado::Fad::DFad<double>>;
