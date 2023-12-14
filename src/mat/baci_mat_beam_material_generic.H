/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief interface for constitutive relations for beam cross-section resultants


\level 1
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_MAT_BEAM_MATERIAL_GENERIC_H
#define BACI_MAT_BEAM_MATERIAL_GENERIC_H

#include "baci_comm_parobjectfactory.H"
#include "baci_inpar_material.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_mat_material.H"

#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


// forward declaration
namespace DRT
{
  class ParObject;
}

namespace MAT
{
  // forward declaration
  namespace PAR
  {
    class BeamElastHyperMaterialParameterGeneric;
  }

  /*---------------------------------------------------------------------------------------------*/
  /// constitutive relations for beam cross-section resultants (hyperelastic stored energy function)
  class BeamMaterial : public Material
  {
   public:
    /**
     * \brief Initialize and setup element specific variables
     *
     */
    virtual void Setup(int numgp_force, int numgp_moment) = 0;

    /** \brief get the radius of a circular cross-section that is ONLY to be used for evaluation of
     *         any kinds of beam interactions (contact, potentials, viscous drag forces ...)
     *
     */
    virtual double GetInteractionRadius() const = 0;

    /** \brief get mass inertia factor with respect to translational accelerations
     *         (usually: density * cross-section area)
     *
     */
    virtual double GetTranslationalMassInertiaFactor() const = 0;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     *
     */
    virtual void GetMassMomentOfInertiaTensorMaterialFrame(CORE::LINALG::Matrix<3, 3>& J) const = 0;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     *
     */
    virtual void GetMassMomentOfInertiaTensorMaterialFrame(
        CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>& J) const = 0;

    /** \brief Update all material related variables at the end of a time step
     *
     */
    virtual void Update() = 0;

    /** \brief Resets all material related variables i.e. in case of adaptive time stepping
     *
     */
    virtual void Reset() = 0;
  };
}  // namespace MAT

BACI_NAMESPACE_CLOSE

#endif
