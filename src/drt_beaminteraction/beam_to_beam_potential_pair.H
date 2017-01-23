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

// Todo get rid of header inclusions
#include "beam3contact_utils.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "Sacado.hpp"

#include "../drt_beaminteraction/beam_potential_pair.H"
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
template<unsigned int numnodes, unsigned int numnodalvalues>
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

  //! first element of contact pair
  DRT::ELEMENTS::Beam3Base const* beam_element1_;

  //! second element of contact pair
  DRT::ELEMENTS::Beam3Base const* beam_element2_;

  //! beam contact parameter list Todo change to correct data container!
//  Teuchos::ParameterList bpotparams_;

  //! line charge conditions
  std::vector<DRT::Condition*> linechargeconds_;

  //! current node coordinates of the two elements
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele1pos_;
  LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> ele2pos_;

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

  //! resulting forces on element 1 and 2
  LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1> fpot1_;
  LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 1> fpot2_;

  //! stiffness contributions
  LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 2*3*numnodes*numnodalvalues> stiffpot1_;
  LINALG::TMatrix<TYPE, 3*numnodes*numnodalvalues, 2*3*numnodes*numnodalvalues> stiffpot2_;


  //@}

  //! @name Private evaluation methods

  /*!
  \brief Evaluate forces and stiffness contribution resulting from potential-based interaction
  */
  void EvaluateFpotandStiffpot_LargeSepApprox();

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
  void ComputeCoords(LINALG::TMatrix<TYPE,3,1>& r,
                     const LINALG::Matrix<1,numnodes*numnodalvalues>& N_i,
                     const LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> elepos);

  //@}

};
}

#endif
