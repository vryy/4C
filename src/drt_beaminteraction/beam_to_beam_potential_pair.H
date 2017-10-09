/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_beam_potential_pair.H

\brief One beam-to-beam potential-based interacting pair (two beam elements)

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BEAMTOBEAMPOTENTIALPAIR_H
#define BEAMTOBEAMPOTENTIALPAIR_H

#include "../drt_beaminteraction/beam_potential_pair.H"

// Todo get rid of header inclusions
#include "beam3contact_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "Sacado.hpp"


#include "../drt_lib/standardtypes_cpp.H"
#include "beam3contact_defines.H"
#include "beam3contact_tangentsmoothing.H"

// forward declarations
namespace LINALG {
  class SerialDenseVector;
  class SerialDenseMatrix;
}
namespace DRT {
namespace ELEMENTS {
  class Beam3Base;
}
}

namespace BEAMINTERACTION
{

/*!
 \brief class for potential-based interaction between two 3D beam elements
 */
template<unsigned int numnodes, unsigned int numnodalvalues, typename T>
class BeamToBeamPotentialPair: public BeamPotentialPair
{
public:

  //! @name Friends
  //no friend classes defined
  //@}

  //! @name Constructors and destructors and related methods
  /*!
  \brief Standard Constructor
  */
  BeamToBeamPotentialPair();

  /*!
  \brief Destructor
  */
  virtual ~BeamToBeamPotentialPair() {};

  //! Setup
  virtual void Setup();
  //@}

  //! @name Derived methods from base class
  /*!
  \brief Evaluate this contact element pair, return value indicates whether pair is active,
         i.e. non-zero values for force and stiffmat are returned
  */
  virtual bool Evaluate(
      LINALG::SerialDenseVector* forcevec1,
      LINALG::SerialDenseVector* forcevec2,
      LINALG::SerialDenseMatrix* stiffmat11,
      LINALG::SerialDenseMatrix* stiffmat12,
      LINALG::SerialDenseMatrix* stiffmat21,
      LINALG::SerialDenseMatrix* stiffmat22,
      const std::vector<DRT::Condition*> linechargeconds,
      const double k,
      const double m);

  /*
  \brief Update state of translational nodal DoFs (absolute positions and tangents) of both elements
  */
  virtual void ResetState(
      const std::vector<double>& centerline_dofvec_ele1,
      const std::vector<double>& centerline_dofvec_ele2);

  /** \brief print this beam potential-based element pair to screen
   *
   *  \author grill */
  virtual void Print(std::ostream& out) const;

  /** \brief print this beam potential element pair to screen
   *
   *  \author grill */
  virtual void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const;
  //@}

  //! @name Access methods
  /*!
  \brief Get ptr to first beam element
  */
  inline const DRT::ELEMENTS::Beam3Base* BeamElement1() { return beam_element1_;};

  /*!
  \brief Get ptr to second beam element
  */
  inline const DRT::ELEMENTS::Beam3Base* BeamElement2() { return beam_element2_;};
  //@}

private:

  //! @name member variables

  //! first element of pair
  DRT::ELEMENTS::Beam3Base const* beam_element1_;

  //! second element of pair
  DRT::ELEMENTS::Beam3Base const* beam_element2_;

  //! line charge conditions
  std::vector<DRT::Condition*> linechargeconds_;

  //! current node coordinates of the two elements
  LINALG::TMatrix<T,3*numnodes*numnodalvalues,1> ele1pos_;
  LINALG::TMatrix<T,3*numnodes*numnodalvalues,1> ele2pos_;

  //! parameters of the applied (point-point) potential law Phi(r)=k_ * r^(-m_)
  double k_;
  double m_;

  //! initial element lengths
  double ele1length_;
  double ele2length_;

  //! Cross-section radius of beam 1
  double radius1_;

  //! Cross-section radius of beam 2
  double radius2_;

  //@}

  //! @name Private evaluation methods

  /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
   *         using double length specific potential with approximation for large separations
   *
   *  \author grill
   *  \date 10/17 */
  void EvaluateFpotandStiffpot_LargeSepApprox(
      LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1>& force_pot1,
      LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1>& force_pot2,
      LINALG::SerialDenseMatrix* stiffmat11,
      LINALG::SerialDenseMatrix* stiffmat12,
      LINALG::SerialDenseMatrix* stiffmat21,
      LINALG::SerialDenseMatrix* stiffmat22);

  /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
   *         Gauss point pair
   *
   *  \author grill
   *  \date 10/17 */
  void EvaluateStiffpotAnalyticContributions_LargeSepApprox(
      LINALG::TMatrix<double, 3, 1> const& dist,
      double const& norm_dist,
      double const& norm_dist_exp1,
      double q1q2_JacFac_GaussWeights,
      LINALG::Matrix<1, numnodes*numnodalvalues> const& N1_i_GP1,
      LINALG::Matrix<1, numnodes*numnodalvalues> const& N2_i_GP2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const;

  /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
   *         Gauss point pair
   *
   *  \author grill
   *  \date 10/17 */
  void EvaluateStiffpotAnalyticContributions_LargeSepApprox(
      LINALG::TMatrix<Sacado::Fad::DFad<double>, 3, 1> const& dist,
      Sacado::Fad::DFad<double> const& norm_dist,
      Sacado::Fad::DFad<double> const& norm_dist_exp1,
      double q1q2_JacFac_GaussWeights,
      LINALG::Matrix<1, numnodes*numnodalvalues> const& N1_i_GP1,
      LINALG::Matrix<1, numnodes*numnodalvalues> const& N2_i_GP2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const
  {
    // this is a dummy since for type FAD, no analytic linearization is required
  }

  /** \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
   *         using double length specific potential with approximation for small separations
   *
   *  \author grill
   *  \date 10/17 */
  void EvaluateFpotandStiffpot_DoubleLengthSpecific_SmallSepApprox(
      LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1>& force_pot1,
      LINALG::TMatrix<T, 3*numnodes*numnodalvalues, 1>& force_pot2,
      LINALG::SerialDenseMatrix* stiffmat11,
      LINALG::SerialDenseMatrix* stiffmat12,
      LINALG::SerialDenseMatrix* stiffmat21,
      LINALG::SerialDenseMatrix* stiffmat22);


  /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
   *         Gauss point pair
   *
   *  \author grill
   *  \date 10/17 */
  void EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
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
      LINALG::SerialDenseMatrix& stiffmat22) const;

  /** \brief compute contributions to analytic linearization (i.e. stiffness matrices) at current
   *         Gauss point pair
   *
   *  \author grill
   *  \date 10/17 */
  void EvaluateStiffpotAnalyticContributions_DoubleLengthSpecific_SmallSepApprox(
      LINALG::TMatrix<Sacado::Fad::DFad<double>, 3, 1> const& dist,
      Sacado::Fad::DFad<double> const& norm_dist,
      Sacado::Fad::DFad<double> const& gap,
      Sacado::Fad::DFad<double> const& gap_exp1,
      double q1q2_JacFac_GaussWeights,
      LINALG::Matrix<1, numnodes*numnodalvalues> const& N1_i_GP1,
      LINALG::Matrix<1, numnodes*numnodalvalues> const& N2_i_GP2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const
  {
    // this is a dummy since for type FAD, no analytic linearization is required
  }

  /** \brief scale given stiffness matrices by given scalar factor
   *
   *  \author grill
   *  \date 10/17 */
  void ScaleStiffpotAnalyticContributions(
      double const& scalefactor,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const;

  /** \brief scale given stiffness matrices by given scalar factor
   *
   *  \author grill
   *  \date 10/17 */
  void ScaleStiffpotAnalyticContributions(
      Sacado::Fad::DFad<double> const& scalefactor,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const
  {
    // this is a dummy since for type FAD, no analytic linearization is required
  }

  /** \brief compute linearization (i.e. stiffness matrices) of given force vectors
   *         using automatic differentiation based on Sacado::Fad package
   *
   *  \author grill
   *  \date 10/17 */
  void CalcStiffmatAutomaticDifferentiationIfRequired(
      LINALG::TMatrix<double, 3*numnodes*numnodalvalues, 1>& force_pot1,
      LINALG::TMatrix<double, 3*numnodes*numnodalvalues, 1>& force_pot2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const
  {
    // this is a dummy since for type double, no automatic differentiation is required
  }

  /** \brief compute linearization (i.e. stiffness matrices) of given force vectors
   *         using automatic differentiation based on Sacado::Fad package
   *
   *  \author grill
   *  \date 10/17 */
  void CalcStiffmatAutomaticDifferentiationIfRequired(
      LINALG::TMatrix<Sacado::Fad::DFad<double>, 3*numnodes*numnodalvalues, 1>& force_pot1,
      LINALG::TMatrix<Sacado::Fad::DFad<double>, 3*numnodes*numnodalvalues, 1>& force_pot2,
      LINALG::SerialDenseMatrix& stiffmat11,
      LINALG::SerialDenseMatrix& stiffmat12,
      LINALG::SerialDenseMatrix& stiffmat21,
      LINALG::SerialDenseMatrix& stiffmat22) const;

  /*!
  \brief Calculate shape function values for given parameter values
  Todo call more general utils method
  */
  void GetShapeFunctions( std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i,
                          std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N2_i,
                          std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i_xi,
                          std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N2_i_xi,
                          DRT::UTILS::IntegrationPoints1D& gausspoints);

  /*!
  \brief Compute coordinates of centreline points from the discretization
    Todo call more general utils method
  */
  void ComputeCoords(LINALG::TMatrix<T,3,1>& r,
                     const LINALG::Matrix<1,numnodes*numnodalvalues>& N_i,
                     const LINALG::TMatrix<T,3*numnodes*numnodalvalues,1> elepos);

  /** \brief set primary variables for FAD if required
   *
   *  \author grill
   *  \date 10/17 */
  void SetAutomaticDifferentiationVariablesIfRequired(
      LINALG::TMatrix<double,3*numnodes*numnodalvalues,1> & ele1centerlinedofvec,
      LINALG::TMatrix<double,3*numnodes*numnodalvalues,1> & ele2centerlinedofvec)
  {
    // do nothing in case of type double (analytic differentiation)
  }

  /** \brief set primary variables for FAD if required
   *
   *  \author grill
   *  \date 10/17 */
  void SetAutomaticDifferentiationVariablesIfRequired(
      LINALG::TMatrix<Sacado::Fad::DFad<double>,3*numnodes*numnodalvalues,1> & ele1centerlinedofvec,
      LINALG::TMatrix<Sacado::Fad::DFad<double>,3*numnodes*numnodalvalues,1> & ele2centerlinedofvec);

  //@}

};
}

#endif
