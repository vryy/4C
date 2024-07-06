/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for materials at Gauss points

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MATERIAL_BASE_HPP
#define FOUR_C_MATERIAL_BASE_HPP

#include "4C_config.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_legacy_enum_definitions_materials.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


/// MAT: materials
namespace Core::Mat
{
  namespace PAR
  {
    class Parameter;
  }

  /*!
   * @brief Interface class for materials in elements
   * The Material class provides a general way to store material history
   * variables in elements. Each element owns one object derived from
   * this class. Each derived class implements a specific material and
   * might even store history data.
   *
   * @note There is one material object to each element, so one material
   * object must store the history data of all Gauss points.
   *
   * Different material classes will (in general) not share a common set
   * of functions, so the element needs to cast to the right
   * material implementation. The element code will typically perform a
   * test on material_type() and do the appropriate cast afterwards.
   *
   * <h3>Storage</h3>
   *
   * In order to support storage of material data (for restart, post
   * processing, ...) this class does implement Core::Communication::ParObject. The
   * elements pack and unpack its materials along with its other data.
   */
  class Material : public Core::Communication::ParObject
  {
   public:
    /// return type of this material
    virtual Core::Materials::MaterialType material_type() const = 0;

    /// return copy of this material object
    virtual Teuchos::RCP<Core::Mat::Material> clone() const = 0;

    /// return quick accessible material parameter data
    virtual Core::Mat::PAR::Parameter* parameter() const = 0;

    /// return material density at gauss point (if provided by the specific material)
    virtual double density(int gp) const { return density(); }

    /// return material density (if provided by the specific material)
    virtual double density() const
    {
      FOUR_C_THROW("The material you are using does not provide a density");
    }

    /// return internal state
    virtual double get_internal_state(int k) const
    {
      FOUR_C_THROW("The material you are using does not provide an internal state");
    }
  };

}  // namespace Core::Mat

FOUR_C_NAMESPACE_CLOSE

#endif
