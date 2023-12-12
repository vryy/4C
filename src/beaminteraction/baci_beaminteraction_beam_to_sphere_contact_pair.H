/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief class to handle contact between a 3D beam element and a rigid sphere

\level 3

*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BACI_BEAMINTERACTION_BEAM_TO_SPHERE_CONTACT_PAIR_H
#define BACI_BEAMINTERACTION_BEAM_TO_SPHERE_CONTACT_PAIR_H

#include "baci_beaminteraction_beam_to_beam_contact_defines.H"
#include "baci_beaminteraction_contact_pair.H"
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_lib_condition.H"
#include "baci_lib_element.H"
#include "baci_lib_node.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include "baci_linalg_sparsematrix.H"
#include "baci_utils_fad.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declaration ...
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
  template <unsigned int numnodes, unsigned int numnodalvalues>
  class BeamToSphereContactPair : public BeamContactPair
  {
   public:
    //! @name Friends
    // no friend classes defined
    //@}

    //! @name Constructors and destructors and related methods

    /*!
    \brief Standard Constructor
    */
    BeamToSphereContactPair();



    //! Setup
    void Setup() override;

    //@}

    //! @name Derived methods from base class
    /*!
    \brief things that need to be done in a separate loop before the actual evaluation loop
           over all contact pairs
    */
    void PreEvaluate() override;

    /*!
    \brief Evaluate this contact element pair, return value indicates whether pair is active,
           i.e. non-zero values for force and stiffmat are returned
    */
    bool Evaluate(CORE::LINALG::SerialDenseVector* forcevec1,
        CORE::LINALG::SerialDenseVector* forcevec2, CORE::LINALG::SerialDenseMatrix* stiffmat11,
        CORE::LINALG::SerialDenseMatrix* stiffmat12, CORE::LINALG::SerialDenseMatrix* stiffmat21,
        CORE::LINALG::SerialDenseMatrix* stiffmat22) override;

    /*
    \brief Update state of translational nodal DoFs (absolute positions (and tangents)) of both
    elements
    */
    void ResetState(const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2) override;

    /** \brief print information about this beam contact element pair to screen
     *
     *  \author grill
     *  \date 05/16 */
    void Print(std::ostream& out) const override;


    /** \brief print this beam contact element pair to screen
     *
     *  \author grill
     *  \date 12/16 */
    void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const override;
    //@}

    //! @name Access methods

    /*!
    \brief Get beam element
    */
    inline DRT::ELEMENTS::Beam3Base const* BeamElement() { return beam_element_; };

    /*!
    \brief Get sphere element
    */
    inline DRT::ELEMENTS::Rigidsphere const* SphereElement() { return sphere_element_; };

    /*!
    \brief Get flag indicating whether contact is active (true) or inactive (false)
    */
    inline bool GetContactFlag() const override
    {
      return (contactflag_ or nodalcontactflag_[0] or nodalcontactflag_[1]);
    };

    /*!
    \brief Get number of active contact point pairs on this element pair
    */
    unsigned int GetNumAllActiveContactPointPairs() const override { return 1; }


    /*!
    \brief Get coordinates of all active contact points on element1 and element2
    */
    inline void GetAllActiveContactPointCoordsElement1(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const override
    {
      dserror("not implemented yet!");
    }

    inline void GetAllActiveContactPointCoordsElement2(
        std::vector<CORE::LINALG::Matrix<3, 1, double>>& coords) const override
    {
      dserror("not implemented yet!");
    }

    /*!
    \brief Get all (scalar) contact forces of this contact pair
    */
    inline void GetAllActiveContactForces(std::vector<double>& forces) const override
    {
      dserror("not implemented yet!");
    }

    /*!
    \brief Get all (scalar) gap values of this contact pair
    */
    void GetAllActiveContactGaps(std::vector<double>& gaps) const override
    {
      dserror("not implemented yet!");
    }

    /*!
    \brief Get energy of penalty contact.
    */
    double GetEnergy() const override
    {
      dserror("not implemented yet!");
      return 0.0;
    }
    //@}


   private:
    //! @name Private evaluation methods

    /*!
    \brief Find contact point via closest point projection
    */
    void ClosestPointProjection();

    /*!
    \brief Utility method for CPP (evaluate nonlinear function f)
    */
    void EvaluateOrthogonalityCondition(TYPE& f, const CORE::LINALG::Matrix<3, 1, TYPE>& delta_x,
        const double norm_delta_x, const CORE::LINALG::Matrix<3, 1, TYPE>& dx1);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void EvaluateLinOrthogonalityCondition(TYPE& df, CORE::LINALG::Matrix<3, 1, TYPE>& delta_x,
        const double norm_delta_x, const CORE::LINALG::Matrix<3, 1, TYPE>& dx1,
        const CORE::LINALG::Matrix<3, 1, TYPE>& ddx1);

    /*!
    \brief Evaluate and assemble contact forces
    */
    void EvaluateFcContact(CORE::LINALG::SerialDenseVector& forcevec1,
        CORE::LINALG::SerialDenseVector& forcevec2, const double& pp, const TYPE& gap,
        const CORE::LINALG::Matrix<3, 1, TYPE>& normal,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        const bool contactactive);

    /*!
    \brief Evaluate and assemble contact stiffness
    */
    void EvaluateStiffcContact(CORE::LINALG::SerialDenseMatrix& stiffmat11,
        CORE::LINALG::SerialDenseMatrix& stiffmat12, CORE::LINALG::SerialDenseMatrix& stiffmat21,
        CORE::LINALG::SerialDenseMatrix& stiffmat22, const double& pp, const TYPE& gap,
        const CORE::LINALG::Matrix<3, 1, TYPE>& normal, const TYPE& norm,
        const CORE::LINALG::Matrix<3, 1, TYPE>& x1, const CORE::LINALG::Matrix<3, 1, TYPE>& x2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& dx1, const CORE::LINALG::Matrix<3, 1, TYPE>& ddx1,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xi,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xixi,
        bool activecontact, bool linxi = true);

    /*!
    \brief Compute normal vector in contact point
    */
    void ComputeNormal(CORE::LINALG::Matrix<3, 1, TYPE>& normal, TYPE& gap, TYPE& norm,
        const CORE::LINALG::Matrix<3, 1, TYPE>& x1, const CORE::LINALG::Matrix<3, 1, TYPE>& x2);

    /*!
    \brief Evaluate gap function
    */
    void ComputeGap(TYPE& gap, const TYPE& norm);

    /*!
    \brief Compute radius of cross section based on moment of inertia
    */
    void ComputeEleRadius(double& radius, const double& moi);

    /*!
    \brief Compute coordinates and their derivatives from the discretization
    */
    void ComputeCoordsAndDerivs(CORE::LINALG::Matrix<3, 1, TYPE>& x1,
        CORE::LINALG::Matrix<3, 1, TYPE>& x2, CORE::LINALG::Matrix<3, 1, TYPE>& dx1,
        CORE::LINALG::Matrix<3, 1, TYPE>& ddx1,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xi,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xixi);

    /*!
    \brief Get shape functions and their derivatives at eta
    */
    void GetShapeFunctions(CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xi,
        CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xixi, const TYPE& eta);

    /*!
    \brief Compute linearizations of contact point
    */
    void ComputeLinXi(CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& x1, const CORE::LINALG::Matrix<3, 1, TYPE>& x2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& dx1, const CORE::LINALG::Matrix<3, 1, TYPE>& ddx1,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i_xi);

    /*!
    \brief Compute linearization of gap
    */
    void ComputeLinGap(CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3, 1, TYPE>& delta_gap,
        CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& x1, const CORE::LINALG::Matrix<3, 1, TYPE>& x2,
        const CORE::LINALG::Matrix<3, 1, TYPE>& dx1,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i, const TYPE& normdist,
        const CORE::LINALG::Matrix<3, 1, TYPE>& normal, const TYPE& norm, const TYPE& gap,
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3, TYPE>& delta_x1_minus_x2);

    /*!
    \brief Compute linearization of normal
    */

    void ComputeLinNormal(
        CORE::LINALG::Matrix<3, 3 * numnodes * numnodalvalues + 3, TYPE>& delta_normal,
        const CORE::LINALG::Matrix<3 * numnodes * numnodalvalues + 3, 1, TYPE>& delta_xi,
        const CORE::LINALG::Matrix<3, 1, TYPE>& normal, const TYPE& norm_delta_x,
        const CORE::LINALG::Matrix<3, 1, TYPE>& x1_xi,
        const CORE::LINALG::Matrix<1, numnodes * numnodalvalues, TYPE>& N1_i);

    /*!
    \brief Compute normal contact disctance
    */
    void ComputeDistance(CORE::LINALG::Matrix<3, 1, TYPE>& distance, TYPE& normdist,
        const CORE::LINALG::Matrix<3, 1, TYPE>& normal, const TYPE& norm);

    /*!
    \brief Check if contact is active and set flag accordingly
    */
    void CheckAndSetContactStatus();

    //@}

   private:
    //! @name member variables

    //! first element of contact pair
    DRT::ELEMENTS::Beam3Base const* beam_element_;

    //! second element of contact pair
    DRT::ELEMENTS::Rigidsphere const* sphere_element_;

    //! current node coordinates of the two elements
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> ele1pos_;
    CORE::LINALG::Matrix<3, 1, TYPE> ele2pos_;

    //! beam element arc-length in stress-free reference configuration
    double beamele_reflength_;

    //! Cross-section radius of beam
    double radius1_;

    //! Cross-section radius of sphere
    double radius2_;

    //! gap function
    TYPE gap_;

    //! flag indicating contact (active/inactive)
    bool contactflag_;

    //! flag indicating contact of beam end points (nodes) (active/inactive)
    std::vector<bool> nodalcontactflag_;

    //! resulting nodal contact forces on ele 1/2
    CORE::LINALG::Matrix<3 * numnodes * numnodalvalues, 1, TYPE> fc1_;
    CORE::LINALG::Matrix<3, 1, TYPE> fc2_;

    //! coordinates of contact point on center lines of beams
    CORE::LINALG::Matrix<3, 1, TYPE> x1_;
    CORE::LINALG::Matrix<3, 1, TYPE> x2_;

    //! parameter value of contact point on beam element
    TYPE xicontact_;

    //! normal vector of current time step
    CORE::LINALG::Matrix<3, 1, TYPE> normal_;

    //@}
  };
}  // namespace BEAMINTERACTION

BACI_NAMESPACE_CLOSE

#endif
