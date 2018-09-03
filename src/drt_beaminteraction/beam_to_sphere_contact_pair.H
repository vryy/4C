/*-----------------------------------------------------------------------------------------------*/
/*!
\file beam_to_sphere_contact_pair.H

\brief class to handle contact between a 3D beam element and a rigid sphere

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef BEAMTOSPHERECONTACTPAIR_H
#define BEAMTOSPHERECONTACTPAIR_H

#include "beam_contact_pair.H"

// Todo get rid of header inclusions
#include "Teuchos_RCP.hpp"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"

#include "beam3contact_defines.H"
#include "../drt_lib/drt_condition.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../headers/FAD_utils.H"

// forward declaration ...
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
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

    /*!
    \brief Destructor
    */
    virtual ~BeamToSphereContactPair(){};

    //! Setup
    virtual void Setup();

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
    virtual bool Evaluate(LINALG::SerialDenseVector* forcevec1,
        LINALG::SerialDenseVector* forcevec2, LINALG::SerialDenseMatrix* stiffmat11,
        LINALG::SerialDenseMatrix* stiffmat12, LINALG::SerialDenseMatrix* stiffmat21,
        LINALG::SerialDenseMatrix* stiffmat22);

    /*
    \brief Update state of translational nodal DoFs (absolute positions (and tangents)) of both
    elements
    */
    virtual void ResetState(const std::vector<double>& centerline_dofvec_ele1,
        const std::vector<double>& centerline_dofvec_ele2);

    /** \brief print information about this beam contact element pair to screen
     *
     *  \author grill
     *  \date 05/16 */
    virtual void Print(std::ostream& out) const;


    /** \brief print this beam contact element pair to screen
     *
     *  \author grill
     *  \date 12/16 */
    virtual void PrintSummaryOneLinePerActiveSegmentPair(std::ostream& out) const;
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
    virtual inline bool GetContactFlag() const
    {
      return (contactflag_ or nodalcontactflag_[0] or nodalcontactflag_[1]);
    };

    /*!
    \brief Get number of active contact point pairs on this element pair
    */
    virtual unsigned int GetNumAllActiveContactPointPairs() const { return 1; }


    /*!
    \brief Get coordinates of all active contact points on element1 and element2
    */
    virtual inline void GetAllActiveContactPointCoordsElement1(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const
    {
      dserror("not implemented yet!");
    }

    virtual inline void GetAllActiveContactPointCoordsElement2(
        std::vector<LINALG::TMatrix<double, 3, 1>>& coords) const
    {
      dserror("not implemented yet!");
    }

    /*!
    \brief Get all (scalar) contact forces of this contact pair
    */
    virtual inline void GetAllActiveContactForces(std::vector<double>& forces) const
    {
      dserror("not implemented yet!");
    }

    /*!
    \brief Get all (scalar) gap values of this contact pair
    */
    virtual void GetAllActiveContactGaps(std::vector<double>& gaps) const
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
    void EvaluateOrthogonalityCondition(TYPE& f, const LINALG::TMatrix<TYPE, 3, 1>& delta_x,
        const double norm_delta_x, const LINALG::TMatrix<TYPE, 3, 1>& dx1);

    /*!
    \brief Utility method for CPP (evaluate Jacobian of nonlinear function f)
    */
    void EvaluateLinOrthogonalityCondition(TYPE& df, LINALG::TMatrix<TYPE, 3, 1>& delta_x,
        const double norm_delta_x, const LINALG::TMatrix<TYPE, 3, 1>& dx1,
        const LINALG::TMatrix<TYPE, 3, 1>& ddx1);

    /*!
    \brief Evaluate and assemble contact forces
    */
    void EvaluateFcContact(LINALG::SerialDenseVector& forcevec1,
        LINALG::SerialDenseVector& forcevec2, const double& pp, const TYPE& gap,
        const LINALG::TMatrix<TYPE, 3, 1>& normal,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i, const bool contactactive);

    /*!
    \brief Evaluate and assemble contact stiffness
    */
    void EvaluateStiffcContact(LINALG::SerialDenseMatrix& stiffmat11,
        LINALG::SerialDenseMatrix& stiffmat12, LINALG::SerialDenseMatrix& stiffmat21,
        LINALG::SerialDenseMatrix& stiffmat22, const double& pp, const TYPE& gap,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm,
        const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2,
        const LINALG::TMatrix<TYPE, 3, 1>& dx1, const LINALG::TMatrix<TYPE, 3, 1>& ddx1,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi, bool activecontact,
        bool linxi = true);

    /*!
    \brief Compute normal vector in contact point
    */
    void ComputeNormal(LINALG::TMatrix<TYPE, 3, 1>& normal, TYPE& gap, TYPE& norm,
        const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2);

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
    void ComputeCoordsAndDerivs(LINALG::TMatrix<TYPE, 3, 1>& x1, LINALG::TMatrix<TYPE, 3, 1>& x2,
        LINALG::TMatrix<TYPE, 3, 1>& dx1, LINALG::TMatrix<TYPE, 3, 1>& ddx1,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi);

    /*!
    \brief Get shape functions and their derivatives at eta
    */
    void GetShapeFunctions(LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
        LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi,
        LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xixi, const TYPE& xi);

    /*!
    \brief Compute linearizations of contact point
    */
    void ComputeLinXi(LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
        const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2,
        const LINALG::TMatrix<TYPE, 3, 1>& dx1, const LINALG::TMatrix<TYPE, 3, 1>& ddx1,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i_xi);

    /*!
    \brief Compute linearization of gap
    */
    void ComputeLinGap(LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_gap,
        LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
        const LINALG::TMatrix<TYPE, 3, 1>& x1, const LINALG::TMatrix<TYPE, 3, 1>& x2,
        const LINALG::TMatrix<TYPE, 3, 1>& dx1,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i, const TYPE& normdist,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm, const TYPE& gap,
        LINALG::TMatrix<TYPE, 3, 3 * numnodes * numnodalvalues + 3>& delta_x1_minus_x2);

    /*!
    \brief Compute linearization of normal
    */

    void ComputeLinNormal(LINALG::TMatrix<TYPE, 3, 3 * numnodes * numnodalvalues + 3>& delta_normal,
        const LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues + 3, 1>& delta_xi,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm_delta_x,
        const LINALG::TMatrix<TYPE, 3, 1>& x1_xi,
        const LINALG::TMatrix<TYPE, 1, numnodes * numnodalvalues>& N1_i);

    /*!
    \brief Compute normal contact disctance
    */
    void ComputeDistance(LINALG::TMatrix<TYPE, 3, 1>& distance, TYPE& normdist,
        const LINALG::TMatrix<TYPE, 3, 1>& normal, const TYPE& norm);

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
    LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues, 1> ele1pos_;
    LINALG::TMatrix<TYPE, 3, 1> ele2pos_;

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
    LINALG::TMatrix<TYPE, 3 * numnodes * numnodalvalues, 1> fc1_;
    LINALG::TMatrix<TYPE, 3, 1> fc2_;

    //! coordinates of contact point on center lines of beams
    LINALG::TMatrix<TYPE, 3, 1> x1_;
    LINALG::TMatrix<TYPE, 3, 1> x2_;

    //! parameter value of contact point on beam element
    TYPE xicontact_;

    //! normal vector of current time step
    LINALG::TMatrix<TYPE, 3, 1> normal_;

    //@}
  };
}  // namespace BEAMINTERACTION

#endif
