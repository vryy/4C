/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief interface for constitutive relations for beam cross-section resultants


\level 1
*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_BEAM_MATERIAL_GENERIC_HPP
#define FOUR_C_MAT_BEAM_MATERIAL_GENERIC_HPP

#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_inpar_material.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"

#include <Sacado.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


// forward declaration
namespace Discret
{
  class ParObject;
}

namespace Mat
{
  // forward declaration
  namespace PAR
  {
    class BeamElastHyperMaterialParameterGeneric;
  }

  /*---------------------------------------------------------------------------------------------*/
  /// constitutive relations for beam cross-section resultants (hyperelastic stored energy function)
  class BeamMaterial : public Core::Mat::Material
  {
   public:
    /**
     * \brief Initialize and setup element specific variables
     *
     */
    virtual void setup(int numgp_force, int numgp_moment) = 0;

    /** \brief get the radius of a circular cross-section that is ONLY to be used for evaluation of
     *         any kinds of beam interactions (contact, potentials, viscous drag forces ...)
     *
     */
    virtual double get_interaction_radius() const = 0;

    /** \brief get mass inertia factor with respect to translational accelerations
     *         (usually: density * cross-section area)
     *
     */
    virtual double get_translational_mass_inertia_factor() const = 0;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     *
     */
    virtual void get_mass_moment_of_inertia_tensor_material_frame(
        Core::LinAlg::Matrix<3, 3>& J) const = 0;

    /** \brief get mass moment of inertia tensor, expressed w.r.t. material frame
     *
     */
    virtual void get_mass_moment_of_inertia_tensor_material_frame(
        Core::LinAlg::Matrix<3, 3, Sacado::Fad::DFad<double>>& J) const = 0;

    /** \brief Update all material related variables at the end of a time step
     *
     */
    virtual void update() = 0;

    /** \brief Resets all material related variables i.e. in case of adaptive time stepping
     *
     */
    virtual void reset() = 0;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
