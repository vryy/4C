/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief One beam-to-sphere potential-based interacting pair

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAMINTERACTION_BEAM_TO_SPHERE_POTENTIAL_PAIR_H
#define BEAMINTERACTION_BEAM_TO_SPHERE_POTENTIAL_PAIR_H

#include "beaminteraction_potential_pair.H"

// Todo get rid of header inclusions
#include "linalg_fixedsizematrix.H"
#include "linalg_sparsematrix.H"
#include "discretization_fem_general_utils_local_connectivity_matrices.H"
#include "discretization_fem_general_utils_integration.H"
#include <Sacado.hpp>

#include "beaminteraction_beam3contact_defines.H"
#include "beaminteraction_beam3contact_tangentsmoothing.H"

// forward declarations
namespace CORE::LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace CORE::LINALG

namespace DRT
{
  namespace ELEMENTS
  {
    class Rigidsphere;
    class Beam3Base;
  }  // namespace ELEMENTS
}  // namespace DRT

namespace BEAMINTERACTION
{
  /*!
   \brief class for potential-based interaction of a 3D beam element and a rigid sphere
   */

  template <unsigned int numnodes, unsigned int numnodalvalues>
  class BeamToSpherePotentialPair : public BEAMINTERACTION::BeamPotentialPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods
    /*!
    \brief Standard Constructor
    */
    BeamToSpherePotentialPair();

    /*!
    \brief Destructor
    */
    virtual ~BeamToSpherePotentialPair(){};

    //! Setup
    virtual void Setup();
    //@}

    //! @name Derived methods from base class
    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    virtual bool Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22,
        const std::vector<DRT::Condition*> linechargeconds, const double k, const double m);

    /*
    \brief Update state of translational nodal DoFs (absolute positions and tangents) of both
    elements
    */
    virtual void ResetState(double time, const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2);

    /*!
    \brief Get coordinates of all interacting points on element1 and element2
    */
    virtual void GetAllInteractingPointCoordsElement1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const
    {
      dserror("not implemented yet");
    }

    virtual void GetAllInteractingPointCoordsElement2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const
    {
      dserror("not implemented yet");
    }

    /*!
    \brief Get forces at all interacting points on element1 and element2
    */
    virtual void GetForcesAtAllInteractingPointsElement1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& forces) const
    {
      dserror("not implemented yet");
    }

    virtual void GetForcesAtAllInteractingPointsElement2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& forces) const
    {
      dserror("not implemented yet");
    }

    /*!
    \brief Get moments at all interacting points on element1 and element2
    */
    virtual void GetMomentsAtAllInteractingPointsElement1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& moments) const
    {
      dserror("not implemented yet");
    }

    virtual void GetMomentsAtAllInteractingPointsElement2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& moments) const
    {
      dserror("not implemented yet");
    }

    /*!
    \brief Get interaction free energy / potential
    */
    double GetEnergy() const override { return interaction_potential_; }

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
    \brief Get first element (beam)
    */
    inline const DRT::ELEMENTS::Beam3Base* BeamElement() { return beam_element_; };

    /*!
    \brief Get second element (sphere)
    */
    inline const DRT::ELEMENTS::Rigidsphere* SphereElement() { return sphere_element_; };
    //@}

   private:
    //! @name member variables

    //! first element of pair
    DRT::ELEMENTS::Beam3Base const* beam_element_;

    //! second element of pair
    DRT::ELEMENTS::Rigidsphere const* sphere_element_;

    //! line and point charge condition
    std::vector<DRT::Condition*> chargeconds_;

    //! current time
    double time_;

    //! current absolute Dof values of the two elements
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele1pos_;
    CORE::LINALG::Matrix<3, 1, TYPE> ele2pos_;

    //! parameters of the applied (point-point) potential law Phi(r)=k_ * r^(-m_)
    double k_;
    double m_;

    //! beam element arc-length in stress-free reference configuration
    double beamele_reflength_;

    //! Cross-section radius of beam
    double radius1_;

    //! Cross-section radius of sphere
    double radius2_;

    //! resulting forces on element 1 and 2
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> fpot1_;
    CORE::LINALG::Matrix<3, 1, TYPE> fpot2_;

    //! stiffness contributions
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 3 * numnodes * numnodalvalues + 3, TYPE>
        stiffpot1_;
    CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3, TYPE> stiffpot2_;

    //! total interaction potential of this pair
    double interaction_potential_;

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
    void GetShapeFunctions(std::vector<CORE::LINALG::Matrix<1, numnodes * numnodalvalues>>& N1_i,
        std::vector<CORE::LINALG::Matrix<1, numnodes * numnodalvalues>>& N1_i_xi,
        CORE::DRT::UTILS::IntegrationPoints1D& gausspoints);

    /*!
    \brief Compute coordinates of centreline points from the discretization
        Todo call more general utils method
    */
    void ComputeCoords(CORE::LINALG::Matrix<3, 1, TYPE>& r,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues>& N_i,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> elepos);

    //@}
  };
}  // namespace BEAMINTERACTION

#endif
