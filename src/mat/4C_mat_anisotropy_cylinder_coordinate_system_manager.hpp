/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a cylinder coordinate system manager

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_MANAGER_HPP
#define FOUR_C_MAT_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_mat_anisotropy_cylinder_coordinate_system_provider.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::Communication
{
  class PackBuffer;
}
namespace Input
{
  class LineDefinition;
}

namespace Mat
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
    void pack(Core::Communication::PackBuffer& data) const;

    /*!
     * Unpack all data from another processor
     *
     * @param data (in) : data object
     * @param position (in/out) : current position in the data
     */
    void unpack(const std::vector<char>& data, std::vector<char>::size_type& position);

    ///@}


    /*!
     * Reads the line definition of an element to get the coordinate system defined on the element
     *
     * @param linedef (in) : Input line of the corresponding element
     */
    void read_from_element_line_definition(Input::LineDefinition* linedef);

    /*!
     * \brief Flag
     *
     * \return true
     * \return false
     */
    bool is_defined() const { return is_defined_; }

    const Core::LinAlg::Matrix<3, 1>& get_rad() const override
    {
      if (!is_defined_)
      {
        FOUR_C_THROW("The coordinate system is not yet defined.");
      }
      return radial_;
    };

    const Core::LinAlg::Matrix<3, 1>& get_axi() const override
    {
      if (!is_defined_)
      {
        FOUR_C_THROW("The coordinate system is not yet defined.");
      }
      return axial_;
    }

    const Core::LinAlg::Matrix<3, 1>& get_cir() const override
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
    void evaluate_local_coordinate_system(Core::LinAlg::Matrix<3, 3>& cosy) const;

   private:
    /// Flag whether coordinate system is already set
    bool is_defined_ = false;

    /*!
     * \brief Unit vector in radial direction
     */
    Core::LinAlg::Matrix<3, 1> radial_;

    /*!
     * \brief unit vector in axial direction
     */
    Core::LinAlg::Matrix<3, 1> axial_;


    /*!
     * \brief unit vector in circumferential direction
     */
    Core::LinAlg::Matrix<3, 1> circumferential_;
  };
}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif
