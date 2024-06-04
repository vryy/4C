/*----------------------------------------------------------------------*/
/*! \file

\brief This file defines a class that holds different nodal fibers

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_DISCRETIZATION_FEM_GENERAL_FIBER_NODE_HOLDER_HPP
#define FOUR_C_DISCRETIZATION_FEM_GENERAL_FIBER_NODE_HOLDER_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_fiber_node.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace CORE::COMM
{
  class PackBuffer;
}

namespace CORE::Nodes
{
  /*!
   * \brief This class holds the different fibers on the Gauss point level.
   *
   * It is internally a std::unordered_map<> that adds the << operator so that it can be stored
   * in a Teuchos::ParameterList.
   */
  class NodalFiberHolder
  {
   public:
    /*!
     * \brief Sets a coordinate system direction of a specific point at every Gauss point
     *
     * @param type Type of the fiber
     * @param fiber Vector of the coordinate system direction at the Gauss-points
     */
    void set_coordinate_system_direction(
        CoordinateSystemDirection type, const std::vector<CORE::LINALG::Matrix<3, 1>>& fiber);

    /*!
     * \brief Returns a reference to coordinate system directions at every Gauss point belonging
     * to the specific type
     *
     * @param type Type of the coordinate system direction
     * @return Reference to the coordinate system direction
     */
    const std::vector<CORE::LINALG::Matrix<3, 1>>& get_coordinate_system_direction(
        CoordinateSystemDirection type) const;

    /*!
     * \brief Returns a reference to coordinate system directions at every Gauss point belonging
     * to the specific type
     *
     * @param type Type of the coordinate system direction
     * @return Reference to the coordinate system direction
     */
    std::vector<CORE::LINALG::Matrix<3, 1>>& get_coordinate_system_direction_mutual(
        CoordinateSystemDirection type);
    /*!
     * \brief Inserts a fiber
     *
     * @param fiber Vector of the fiber at the Gauss-points
     */
    void AddFiber(const std::vector<CORE::LINALG::Matrix<3, 1>>& fiber);

    /*!
     * @brief Return the list of all fibers at all Gauss points
     *
     * @return Reference to the list of all fibers at all Gauss points
     */
    const std::vector<std::vector<CORE::LINALG::Matrix<3, 1>>>& GetFibers() const
    {
      return fibers_;
    }

    /*!
     * \brief Returns a reference to the specific fiber at every Gauss point
     *
     * @param fiberid Id of the fiber (starting from 0)
     * @return Reference to the coordinate system direction
     */
    const std::vector<CORE::LINALG::Matrix<3, 1>>& GetFiber(std::size_t fiberid) const;

    /*!
     * \brief Returns a reference to the specific fiber at every Gauss point
     *
     * @param fiberid Id of the fiber (starting from 0)
     * @return Reference to the fiber
     */
    std::vector<CORE::LINALG::Matrix<3, 1>>& GetFiberMutual(std::size_t fiberid);

    /*!
     * \brief Sets a type of angle
     *
     * @param type Type of the angle
     * @param fiber angle value
     */
    void SetAngle(AngleType type, const std::vector<double>& angle);

    /*!
     * \brief Returns the value of the angle of the type for each Gauss point
     *
     * @param type Type of the angle
     * @return value of the angle
     */
    const std::vector<double>& GetAngle(AngleType type) const;

    /*!
     * \brief Checks, whether the holder contains fibers of the type @type
     *
     * \param type Type of the fiber
     * \return true type is contained in the holder
     * \return false type is not contained in the holder
     */
    bool contains_coordinate_system_direction(CoordinateSystemDirection type) const;

    /**
     * Returns the number of coordinate system directions stored in the fiber holder
     * @return number of coordinate system directions
     */
    std::size_t coordinate_system_size() const;

    /**
     * Returns the number of fibers stored in the fiber holder
     * @return number of fibers
     */
    std::size_t FibersSize() const;

    /**
     * Returns the number of angles stored in the fiber holder
     * @return number of fibers
     */
    std::size_t AnglesSize() const;

    /*!
     * \brief Print a short summary of the Map into an output stream.
     *
     * This operator is needed to be able to store this map in a Teuchos::ParameterList
     *
     * @param os The output stream
     * @param holder The holder that should be printed
     * @return The output stream
     */
    friend std::ostream& operator<<(std::ostream& os, const NodalFiberHolder& holder)
    {
      os << "NodalFiberHolder {Cosy count=" << holder.coordinate_system_size()
         << ", Fibers count=" << holder.FibersSize() << " Angles count = " << holder.AnglesSize()
         << "}";
      return os;
    }

    /*!
     * \brief Returns true of the map of both NodalFiberHolders are the same
     * @param holder1 reference to the first holder
     * @param holder2 reference to the second holder
     * @return
     */
    friend bool operator==(const NodalFiberHolder& holder1, const NodalFiberHolder& holder2)
    {
      return holder1.fibers_ == holder2.fibers_ && holder1.angles_ == holder2.angles_ &&
             holder1.coordinate_system_directions_ == holder2.coordinate_system_directions_;
    }

   private:
    /// Map holding GP coordinate system directions
    std::map<CoordinateSystemDirection, std::vector<CORE::LINALG::Matrix<3, 1>>>
        coordinate_system_directions_;

    /// List of numbered fibers
    std::vector<std::vector<CORE::LINALG::Matrix<3, 1>>> fibers_;

    /// Map holding GP angles of different names
    std::map<AngleType, std::vector<double>> angles_;
  };
}  // namespace CORE::Nodes
FOUR_C_NAMESPACE_CLOSE

#endif
