/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a cylinder coordinate system manager

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_MANAGER_HPP
#define FOUR_C_MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_MANAGER_HPP

#include "baci_config.hpp"

#include "baci_mat_anisotropy_cylinder_coordinate_system_provider.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace CORE::COMM
{
  class PackBuffer;
}
namespace INPUT
{
  class LineDefinition;
}

namespace MAT
{
  /*!
   * \brief A manager that handles reading of a cylinder coordinate manager, distribution to
   * multiple processors and getter/setter methods
   */
  class CylinderCoordinateSystemManager : public CylinderCoordinateSystemProvider
  {
   public:
    /*!
     * \brief Constructor of the cylinder coordinate system manager
     */
    explicit CylinderCoordinateSystemManager();

    ///@name Packing and Unpacking
    ///@{
    /*!
     * Pack all data for parallel distribution
     *
     * @param data (in/out) : data object
     */
    void Pack(CORE::COMM::PackBuffer& data) const;

    /*!
     * Unpack all data from another processor
     *
     * @param data (in) : data object
     * @param position (in/out) : current position in the data
     */
    void Unpack(const std::vector<char>& data, std::vector<char>::size_type& position);

    ///@}


    /*!
     * Reads the line definition of an element to get the coordinate system defined on the element
     *
     * @param linedef (in) : Input line of the corresponding element
     */
    void ReadFromElementLineDefinition(INPUT::LineDefinition* linedef);

    /*!
     * \brief Flag
     *
     * \return true
     * \return false
     */
    bool IsDefined() const { return is_defined_; }

    const CORE::LINALG::Matrix<3, 1>& GetRad() const override
    {
      if (!is_defined_)
      {
        FOUR_C_THROW("The coordinate system is not yet defined.");
      }
      return radial_;
    };

    const CORE::LINALG::Matrix<3, 1>& GetAxi() const override
    {
      if (!is_defined_)
      {
        FOUR_C_THROW("The coordinate system is not yet defined.");
      }
      return axial_;
    }

    const CORE::LINALG::Matrix<3, 1>& GetCir() const override
    {
      if (!is_defined_)
      {
        FOUR_C_THROW("The coordinate system is not yet defined.");
      }
      return circumferential_;
    };

    /*!
     * \brief Evaluation
     *
     * \param cosy
     */
    void EvaluateLocalCoordinateSystem(CORE::LINALG::Matrix<3, 3>& cosy) const;

   private:
    /// Flag whether coordinate system is already set
    bool is_defined_ = false;

    /*!
     * \brief Unit vector in radial direction
     */
    CORE::LINALG::Matrix<3, 1> radial_;

    /*!
     * \brief unit vector in axial direction
     */
    CORE::LINALG::Matrix<3, 1> axial_;


    /*!
     * \brief unit vector in circumferential direction
     */
    CORE::LINALG::Matrix<3, 1> circumferential_;
  };
}  // namespace MAT

FOUR_C_NAMESPACE_CLOSE

#endif
