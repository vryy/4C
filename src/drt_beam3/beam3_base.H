/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief base class for all beam elements

\level 2

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BEAM3_BASE_H_
#define BEAM3_BASE_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "Epetra_Vector.h"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_beaminteraction.H"
#include "../drt_lib/drt_element.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include "beam_spatial_discretization_utils.H"

/* number of linear subsegments used for visualization of beam centerline
 * in vtu format (see post_drt_vtu_writer.cpp and beam_discretization_runtime_vtu_writer.cpp) */
#define BEAMSVTUVISUALSUBSEGMENTS 5

// forward declaration ...
namespace LINALG
{
  class SerialDenseVector;
  class SerialDenseMatrix;
}  // namespace LINALG
namespace GEO
{
  namespace MESHFREE
  {
    class BoundingBox;
  }
}  // namespace GEO
namespace STR
{
  namespace ELEMENTS
  {
    class ParamsInterface;
  }
}  // namespace STR

namespace BROWNIANDYN
{
  class ParamsInterface;
}

namespace MAT
{
  class BeamElastHyperMaterial;
}

namespace DRT
{
  namespace ELEMENTS
  {
    //! base class for all beam elements
    class Beam3Base : public DRT::Element
    {
     public:
      /*!
      \brief Standard Constructor

      \param id    (in): A globally unique element id
      \param etype (in): Type of element
      \param owner (in): owner processor of the element
      */
      Beam3Base(int id, int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a Element
      */
      Beam3Base(const Beam3Base& old);

      //! Destructor
      virtual ~Beam3Base(){};

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this element

      */
      virtual void Unpack(const std::vector<char>& data);

      /** \brief set the parameter interface ptr for the solid elements
       *
       *  \param p (in): Parameter list coming from the time integrator.
       *
       *  \author hiermeier
       *  \date 04/16 */
      virtual void SetParamsInterfacePtr(const Teuchos::ParameterList& p);

      virtual void SetBrownianDynParamsInterfacePtr();

      /** \brief returns true if the parameter interface is defined and initialized, otherwise false
       *
       *  \author hiermeier
       *  \date 04/16 */
      virtual inline bool IsParamsInterface() const { return (not interface_ptr_.is_null()); }

      /** \brief get access to the parameter interface pointer
       *
       *  \author hiermeier
       *  \date 04/16 */
      virtual Teuchos::RCP<DRT::ELEMENTS::ParamsInterface> ParamsInterfacePtr();
      virtual Teuchos::RCP<BROWNIANDYN::ParamsInterface> BrownianDynParamsInterfacePtr() const;

      //! computes the number of different random numbers required in each time step for generation
      //! of stochastic forces
      virtual int HowManyRandomNumbersINeed() const = 0;

      /** \brief get access to the element reference length
       *        (i.e. arc-length in stress-free configuration)
       *
       *  \author grill
       *  \date 05/16 */
      virtual double RefLength() const = 0;

      /** \brief get the radius of the element which is used for interactions (contact, viscous,
       *         potential-based, ...)
       *         - if needed, extend this to other than circular cross-section shapes and dimensions
       *           to be specified via input file
       *         - allow for different assumed shapes for different interaction types if needed
       *
       *  \author grill
       *  \date 02/17 */
      double GetCircularCrossSectionRadiusForInteractions() const;

      /** \brief get number of nodes used for centerline interpolation
       *
       *  \author grill
       *  \date 05/16 */
      virtual int NumCenterlineNodes() const = 0;

      /** \brief find out whether given node is used for centerline interpolation
       *
       *  \author grill
       *  \date 10/16 */
      virtual bool IsCenterlineNode(const DRT::Node& node) const = 0;

      /** \brief return flag whether Hermite polynomials are applied for centerline interpolation
       *
       *  \author grill
       *  \date 05/16 */
      virtual bool HermiteCenterlineInterpolation() const = 0;

      /** \brief return GIDs of all additive DoFs for a given node
       *
       *  \author grill
       *  \date 07/16 */
      std::vector<int> GetAdditiveDofGIDs(
          const DRT::Discretization& discret, const DRT::Node& node) const;

      /** \brief return GIDs of all non-additive, i.e. rotation pseudo vector DoFs for a given node
       *
       *  \author grill
       *  \date 07/16 */
      std::vector<int> GetRotVecDofGIDs(
          const DRT::Discretization& discret, const DRT::Node& node) const;

      /** \brief add indices of those DOFs of a given node that are positions
       *
       *  \author grill
       *  \date 07/16 */
      virtual void PositionDofIndices(std::vector<int>& posdofs, const DRT::Node& node) const = 0;

      /** \brief add indices of those DOFs of a given node that are tangents (in the case of Hermite
       * interpolation)
       *
       *  \author grill
       *  \date 07/16 */
      virtual void TangentDofIndices(std::vector<int>& tangdofs, const DRT::Node& node) const = 0;

      /** \brief add indices of those DOFs of a given node that are rotation DOFs (non-additive
       * rotation vectors)
       *
       *  \author grill
       *  \date 07/16 */
      virtual void RotationVecDofIndices(
          std::vector<int>& rotvecdofs, const DRT::Node& node) const = 0;

      /** \brief add indices of those DOFs of a given node that are 1D rotation DOFs
       *         (planar rotations are additive, e.g. in case of relative twist DOF of beam3k with
       * rotvec=false)
       *
       *  \author grill
       *  \date 07/16 */
      virtual void Rotation1DDofIndices(
          std::vector<int>& twistdofs, const DRT::Node& node) const = 0;

      /** \brief add indices of those DOFs of a given node that represent norm of tangent vector
       *         (additive, e.g. in case of beam3k with rotvec=true)
       *
       *  \author grill
       *  \date 07/16 */
      virtual void TangentLengthDofIndices(
          std::vector<int>& tangnormdofs, const DRT::Node& node) const = 0;

      /** \brief get element local indices of those Dofs that are used for centerline interpolation
       *
       *  \author grill
       *  \date 12/16 */
      virtual void CenterlineDofIndicesOfElement(
          std::vector<unsigned int>& centerlinedofindices) const = 0;

      /** \brief get Jacobi factor ds/dxi(xi) at xi \in [-1;1]
       *
       *  \author grill
       *  \date 06/16 */
      virtual double GetJacobiFacAtXi(const double& xi) const = 0;

      /** \brief Get material cross-section deformation measures, i.e. strain resultants
       *
       *  \author grill
       *  \date 04/17 */
      virtual inline void GetMaterialStrainResultantsAtAllGPs(std::vector<double>& axial_strain_GPs,
          std::vector<double>& shear_strain_2_GPs, std::vector<double>& shear_strain_3_GPs,
          std::vector<double>& twist_GPs, std::vector<double>& curvature_2_GPs,
          std::vector<double>& curvature_3_GPs) const
      {
        dserror("not implemented");
      }

      /** \brief Get spatial cross-section stress resultants
       *
       *  \author eichinger
       *  \date 05/17 */
      virtual inline void GetSpatialStressResultantsAtAllGPs(
          std::vector<double>& spatial_axial_force_GPs,
          std::vector<double>& spatial_shear_force_2_GPs,
          std::vector<double>& spatial_shear_force_3_GPs, std::vector<double>& spatial_torque_GPs,
          std::vector<double>& spatial_bending_moment_2_GPs,
          std::vector<double>& spatial_bending_moment_3_GPs) const
      {
        dserror("not implemented");
      }

      /** \brief Get spatial cross-section stress resultants
       *
       *  \author eichinger
       *  \date 05/17 */
      virtual inline void GetSpatialForcesAtAllGPs(std::vector<double>& spatial_axial_force_GPs,
          std::vector<double>& spatial_shear_force_2_GPs,
          std::vector<double>& spatial_shear_force_3_GPs) const
      {
        dserror("not implemented");
      }

      /** \brief Get spatial cross-section stress resultants
       *
       *  \author eichinger
       *  \date 05/17 */
      virtual inline void GetSpatialMomentsAtAllGPs(std::vector<double>& spatial_torque_GPs,
          std::vector<double>& spatial_bending_moment_2_GPs,
          std::vector<double>& spatial_bending_moment_3_GPs) const
      {
        dserror("not implemented");
      }

      /** \brief Get material cross-section stress resultants
       *
       *  \author grill
       *  \date 04/17 */
      virtual inline void GetMaterialStressResultantsAtAllGPs(
          std::vector<double>& material_axial_force_GPs,
          std::vector<double>& material_shear_force_2_GPs,
          std::vector<double>& material_shear_force_3_GPs, std::vector<double>& material_torque_GPs,
          std::vector<double>& material_bending_moment_2_GPs,
          std::vector<double>& material_bending_moment_3_GPs) const
      {
        dserror("not implemented");
      }

      /** \brief Get number of degrees of freedom of a single node
       *
       *  \author eichinger
       *  \date 08/16 */
      virtual int NumDofPerNode(const DRT::Node& node) const
      {
        dserror("not implemented");
        return -1;
      }

      /** \brief get centerline position at xi \in [-1,1] (element parameter space) in stress-free
       * reference configuration
       *
       *  \author grill
       *  \date 06/16 */
      void GetRefPosAtXi(LINALG::Matrix<3, 1>& refpos, const double& xi) const;

      /** \brief get unit tangent vector in reference configuration at i-th node of beam element
       * (element-internal numbering)
       *
       *  \author grill
       *  \date 06/16 */
      virtual void GetRefTangentAtNode(LINALG::Matrix<3, 1>& Tref_i, const int& i) const = 0;

      /** \brief get centerline position at xi \in [-1,1] (element parameter space) from
       * displacement state vector
       *
       *  \author grill
       *  \date 06/16 */
      virtual void GetPosAtXi(
          LINALG::Matrix<3, 1>& pos, const double& xi, const std::vector<double>& disp) const = 0;

      /** \brief get triad at xi \in [-1,1] (element parameter space)
       *
       *  \author grill
       *  \date 07/16 */
      virtual void GetTriadAtXi(
          LINALG::Matrix<3, 3>& triad, const double& xi, const std::vector<double>& disp) const
      {
        // ToDo make pure virtual and add/generalize implementations in beam eles
        dserror("not implemented");
      }

      /** \brief get generalized interpolation matrix which yields the variation of the position and
       *         orientation at xi \in [-1,1] if multiplied with the vector of primary DoF
       * variations
       *
       *  \author grill
       *  \date 11/16 */
      virtual void GetGeneralizedInterpolationMatrixVariationsAtXi(
          LINALG::SerialDenseMatrix& Ivar, const double& xi, const std::vector<double>& disp) const
      {
        dserror("not implemented");
      }

      /** \brief get linearization of the product of (generalized interpolation matrix for
       * variations (see above) and applied force vector) with respect to the primary DoFs of this
       * element
       *
       *  \author grill
       *  \date 01/17 */
      virtual void GetStiffmatResultingFromGeneralizedInterpolationMatrixAtXi(
          LINALG::SerialDenseMatrix& stiffmat, const double& xi, const std::vector<double>& disp,
          const LINALG::SerialDenseVector& force) const
      {
        dserror("not implemented");
      }

      /** \brief get generalized interpolation matrix which yields the increments of the position
       * and orientation at xi \in [-1,1] if multiplied with the vector of primary DoF increments
       *
       *  \author grill
       *  \date 11/16 */
      virtual void GetGeneralizedInterpolationMatrixIncrementsAtXi(
          LINALG::SerialDenseMatrix& Iinc, const double& xi, const std::vector<double>& disp) const
      {
        dserror("not implemented");
      }

      //! get internal (elastic) energy of element
      virtual double GetInternalEnergy() const = 0;

      //! get kinetic energy of element
      virtual double GetKineticEnergy() const = 0;

      //! shifts nodes so that proper evaluation is possible even in case of periodic boundary
      //! conditions
      virtual void UnShiftNodePosition(std::vector<double>& disp,  //!< element disp vector
          GEO::MESHFREE::BoundingBox const& periodic_boundingbox) const;

      //! get directions in which element might be cut by a periodic boundary
      virtual void GetDirectionsOfShifts(std::vector<double>& disp,
          GEO::MESHFREE::BoundingBox const& periodic_boundingbox,
          std::vector<bool>& shift_in_dim) const;

      /** \brief extract values for those Dofs relevant for centerline-interpolation from total
       * state vector
       *
       *  \author grill
       *  \date 11/16 */
      virtual void ExtractCenterlineDofValuesFromElementStateVector(
          const std::vector<double>& dofvec, std::vector<double>& dofvec_centerline,
          bool add_reference_values = false) const = 0;

     protected:
      /** \brief get access to the interface
       *
       *  \author hiermeier
       *  \date 04/16 */
      inline STR::ELEMENTS::ParamsInterface& ParamsInterface() const
      {
        if (not IsParamsInterface()) dserror("The interface ptr is not set!");
        return *interface_ptr_;
      }

      inline BROWNIANDYN::ParamsInterface& BrownianDynParamsInterface() const
      {
        return *browndyn_interface_ptr_;
      }

      /** \brief compute beam centerline position vector at position \xi in element parameter space
       * [-1,1] via interpolation of nodal DoFs based on given shape function values
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void Calc_r(const LINALG::Matrix<3 * vpernode * nnode, 1, T>& disp_totlag_centerline,
          const LINALG::Matrix<1, vpernode * nnode, double>& funct,
          LINALG::Matrix<3, 1, T>& r) const
      {
        DRT::UTILS::BEAM::CalcInterpolation<nnode, vpernode, 3, T>(
            disp_totlag_centerline, funct, r);
      }

      /** \brief compute derivative of beam centerline (i.e. tangent vector) at position \xi in
       *         element parameter space [-1,1] with respect to \xi via interpolation of nodal DoFs
       *         based on given shape function derivative values
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void Calc_r_xi(const LINALG::Matrix<3 * vpernode * nnode, 1, T>& disp_totlag_centerline,
          const LINALG::Matrix<1, vpernode * nnode, double>& deriv,
          LINALG::Matrix<3, 1, T>& r_xi) const
      {
        DRT::UTILS::BEAM::CalcInterpolation<nnode, vpernode, 3, T>(
            disp_totlag_centerline, deriv, r_xi);
      }

      /** \brief compute derivative of beam centerline (i.e. tangent vector) at position \xi in
       *         element parameter space [-1,1] with respect to arc-length parameter s in reference
       *         configuration via interpolation of nodal DoFs based on given shape function
       * derivative values
       *
       *  \author grill
       *  \date 03/16 */
      template <unsigned int nnode, unsigned int vpernode, typename T>
      void Calc_r_s(const LINALG::Matrix<3 * vpernode * nnode, 1, T>& disp_totlag_centerline,
          const LINALG::Matrix<1, vpernode * nnode, double>& deriv, const double& jacobi,
          LINALG::Matrix<3, 1, T>& r_s) const
      {
        Calc_r_xi<nnode, vpernode, T>(disp_totlag_centerline, deriv, r_s);

        /* at this point we have computed derivative with respect to the element parameter \xi \in
         * [-1;1]; as we want derivative with respect to the reference arc-length parameter s, we
         * have to divide it by the Jacobi determinant at the respective point*/
        r_s.Scale(1.0 / jacobi);
      }

      /** \brief get applied beam material law object
       *
       *  \author grill
       *  \date 02/17 */
      const MAT::BeamElastHyperMaterial& GetBeamMaterial() const;

      /** \brief setup constitutive matrices from material law
       *
       *  \author grill
       *  \date 03/16 */
      template <typename T>
      void GetConstitutiveMatrices(LINALG::Matrix<3, 3, T>& CN, LINALG::Matrix<3, 3, T>& CM) const;

      /** \brief setup mass inertia tensors from material law
       *
       *  \author grill
       *  \date 03/16 */
      template <typename T>
      void GetTranslationalAndRotationalMassInertiaTensor(
          double& mass_inertia_translational, LINALG::Matrix<3, 3, T>& J) const;

      /** \brief setup only translational mass inertia factor from material law
       *      this method is called by reduced beam formulation which don't inlcude
       *      rotational mass inertia
       *
       *  \author grill
       *  \date 03/17 */
      void GetTranslationalMassInertiaFactor(double& mass_inertia_translational) const;

      //! @name Methods and variables for Brownian dynamics or beaminteraction simulations
      //! @{
      //! computes damping coefficients
      void GetDampingCoefficients(LINALG::Matrix<3, 1>& gamma) const;

      //! computes velocity of background fluid and gradient of that velocity at a certain
      //! evaluation point in the physical space and adds respective terms to internal forces and
      //! damping matrix
      template <unsigned int ndim, typename T>  // number of dimensions of embedding space
      void GetBackgroundVelocity(Teuchos::ParameterList& params,  //!< parameter list
          const LINALG::Matrix<ndim, 1, T>&
              evaluationpoint,  //!< point at which background velocity and its gradient has to be
                                //!< computed
          LINALG::Matrix<ndim, 1, T>& velbackground,  //!< velocity of background fluid
          LINALG::Matrix<ndim, ndim, T>& velbackgroundgrad)
          const;  //!< gradient of velocity of background fluid

     public:
      //! get centerline pos at binding spot with locn x stored in element parameter space
      //! coordinates \in [-1,1] from displacement state vector
      void GetPosOfBindingSpot(LINALG::Matrix<3, 1>& pos, std::vector<double>& disp,
          INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int bspotlocn,
          GEO::MESHFREE::BoundingBox const& periodic_boundingbox) const;

      //! get triad at binding spot with locn x stored in element parameter space coordinates \in
      //! [-1,1] from displacement state vector
      void GetTriadOfBindingSpot(LINALG::Matrix<3, 3>& triad, std::vector<double>& disp,
          INPAR::BEAMINTERACTION::CrosslinkerType linkertype, int bspotlocn) const;

      /** \brief get entire binding spot information of element
       *
       *  \author eichinger
       *  \date 06/17 */
      std::map<INPAR::BEAMINTERACTION::CrosslinkerType, std::vector<double>> const&
      GetBindingSpots() const
      {
        return bspotposxi_;
      }

      /** \brief get number of binding spot types on this element
       *
       *  \author eichinger
       *  \date 06/17 */
      unsigned int GetNumberOfBindingSpotTypes() const { return bspotposxi_.size(); }

      /** \brief get number of binding spots of certain binding spot type on this element
       *
       *  \author eichinger
       *  \date 06/17 */
      unsigned int GetNumberOfBindingSpots(INPAR::BEAMINTERACTION::CrosslinkerType linkertype) const
      {
        return bspotposxi_.at(linkertype).size();
      }

      /** \brief get binding spot positions xi
       *
       *  \author eichinger
       *  \date 03/17 */
      double GetBindingSpotXi(
          INPAR::BEAMINTERACTION::CrosslinkerType linkertype, unsigned int bspotlocn) const
      {
        if (bspotlocn > bspotposxi_.at(linkertype).size())
          dserror("number of requested binding spot exceeds total number of binding spots");

        return bspotposxi_.at(linkertype)[bspotlocn];
      }

      /** \brief set binding spot positions and status in crosslinker simulation
       *
       *  \author eichinger
       *  \date 03/17 */
      void SetBindingSpots(
          std::map<INPAR::BEAMINTERACTION::CrosslinkerType, std::vector<double>> bspotposxi)
      {
        bspotposxi_.clear();
        bspotposxi_ = bspotposxi;
      }

      /** \brief set binding spot positions and status in crosslinker simulation
       *
       *  \author eichinger
       *  \date 03/17 */
      void SetPositionsOfBindingSpotType(
          INPAR::BEAMINTERACTION::CrosslinkerType linkertype, std::vector<double> const& bspotposxi)
      {
        bspotposxi_[linkertype] = bspotposxi;
      }

      /** \brief set/get type of filament the element is part of
       *
       *  \author eichinger
       *  \date 03/17 */
      void SetFilamentType(INPAR::BEAMINTERACTION::FilamentType filamenttype)
      {
        filamenttype_ = filamenttype;
      }

      INPAR::BEAMINTERACTION::FilamentType GetFilamentType() const { return filamenttype_; }


     private:
      //! position of binding spots on beam element in local coordinate system
      //! size of vector equals number of binding spots on this element
      std::map<INPAR::BEAMINTERACTION::CrosslinkerType, std::vector<double>> bspotposxi_;

      //! type of filament element belongs to
      INPAR::BEAMINTERACTION::FilamentType filamenttype_;

      //! @}

     private:
      /** \brief interface ptr
       *
       *  data exchange between the element and the time integrator. */
      Teuchos::RCP<STR::ELEMENTS::ParamsInterface> interface_ptr_;

      Teuchos::RCP<BROWNIANDYN::ParamsInterface> browndyn_interface_ptr_;

      /*!
      \brief Default Constructor must not be called
      */
      Beam3Base();
    };

  }  // namespace ELEMENTS
}  // namespace DRT

#endif
