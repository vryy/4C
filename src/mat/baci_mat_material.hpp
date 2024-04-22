/*----------------------------------------------------------------------*/
/*! \file
\brief Interface class for materials at Gauss points

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_MATERIAL_HPP
#define FOUR_C_MAT_MATERIAL_HPP



#include "baci_config.hpp"

#include "baci_comm_parobject.hpp"
#include "baci_inpar_material.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


/// MAT: materials
namespace MAT
{
  namespace PAR
  {
    class Parameter;
  }

  const int NUM_STRESS_3D = 6;  ///< 6 stresses for 3D

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
   * test on MaterialType() and do the appropriate cast afterwards.
   *
   * <h3>Storage</h3>
   *
   * In order to support storage of material data (for restart, post
   * processing, ...) this class does implement CORE::COMM::ParObject. The
   * elements pack and unpack its materials along with its other data.
   */
  class Material : public CORE::COMM::ParObject
  {
   public:
    /// return type of this material
    virtual INPAR::MAT::MaterialType MaterialType() const = 0;

    /// return copy of this material object
    virtual Teuchos::RCP<Material> Clone() const = 0;

    /// return quick accessible material parameter data
    virtual MAT::PAR::Parameter* Parameter() const = 0;

    /// create element material object given the number of a material definition
    static Teuchos::RCP<Material> Factory(int matnum  ///< material ID
    );

    /// return material density at gauss point (if provided by the specific material)
    virtual double Density(int gp) const { return Density(); }

    /// return material density (if provided by the specific material)
    virtual double Density() const
    {
      FOUR_C_THROW("The material you are using does not provide a density");
      return 0.0;
    }

    /// return internal state
    virtual double GetInternalState(int k) const
    {
      FOUR_C_THROW("The material you are using does not provide an internal state");
      return 0.0;
    }
  };

}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
