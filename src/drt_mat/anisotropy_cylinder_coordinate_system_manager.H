/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a cylinder coordinate system manager

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_MANAGER_H
#define BACI_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_MANAGER_H

#include "anisotropy_cylinder_coordinate_system_provider.H"

// forward declarations
namespace DRT
{
  class PackBuffer;

  namespace INPUT
  {
    class LineDefinition;
  }
}  // namespace DRT

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
    void Pack(DRT::PackBuffer& data) const;

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
    void ReadFromElementLineDefinition(DRT::INPUT::LineDefinition* linedef);

    /*!
     * \brief Flag
     *
     * \return true
     * \return false
     */
    bool IsDefined() const { return isDefined_; }

    /*!
     * \brief Returns the unit vector pointing in radial direction
     *
     * \return const LINALG::Matrix<3, 1>& Unit vector pointing in radial direction
     */
    const LINALG::Matrix<3, 1>& GetRad() const override
    {
      if (!isDefined_)
      {
        dserror("The coordinate system is not yet defined.");
      }
      return radial_;
    };

    /*!
     * \brief Returns the unit vector pointing in axial direction
     *
     * \return const LINALG::Matrix<3, 1>& Unit vector pointing in axial direction
     */
    const LINALG::Matrix<3, 1>& GetAxi() const override
    {
      if (!isDefined_)
      {
        dserror("The coordinate system is not yet defined.");
      }
      return axial_;
    }

    /*!
     * \brief Returns the unit vector pointing in circumferential direction
     *
     * \return const LINALG::Matrix<3, 1>& Unit vector pointing in circumferential direction
     */
    const LINALG::Matrix<3, 1>& GetCir() const override
    {
      if (!isDefined_)
      {
        dserror("The coordinate system is not yet defined.");
      }
      return circumferential_;
    };

    /*!
     * \brief Evaluation
     *
     * \param cosy
     */
    void EvaluateLocalCoordinateSystem(LINALG::Matrix<3, 3>& cosy) const;

   private:
    /// Flag whether coordinate system is already set
    bool isDefined_ = false;

    /*!
     * \brief List of unit vectors in radial direction (normally 1 for an element coordinate system,
     * otherwise one per GP)
     */
    LINALG::Matrix<3, 1> radial_;

    /*!
     * \brief List of unit vectors in axial direction (normally 1 for an element coordinate system,
     * otherwise one per GP)
     */
    LINALG::Matrix<3, 1> axial_;


    /*!
     * \brief List of unit vectors in circumferential direction (normally 1 for an element
     * coordinate system, otherwise one per GP)
     */
    LINALG::Matrix<3, 1> circumferential_;
  };
}  // namespace MAT

#endif  // BACI_ANISOTROPY_CYLINDER_COORDINATE_SYSTEM_MANAGER_H
