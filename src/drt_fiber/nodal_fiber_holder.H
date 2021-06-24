/*----------------------------------------------------------------------*/
/*! \file

\brief This file defines a class that holds different nodal fibers

\level 3


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_NODAL_FIBER_HOLDER_H
#define BACI_NODAL_FIBER_HOLDER_H

#include <vector>
#include <unordered_map>
#include "../linalg/linalg_fixedsizematrix.H"
#include "drt_fiber_node.H"

namespace DRT
{
  class PackBuffer;
  namespace FIBER
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
       * \brief Sets a type of fiber
       *
       * @param type Type of the fiber
       * @param fiber Vector of the fibers at the Gauss-points
       */
      void SetFiber(FiberType type, const std::vector<LINALG::Matrix<3, 1>>& fiber);

      /*!
       * \brief Returns a reference to a list of fibers belonging to the specific type
       *
       * @param type Type of the fiber
       * @return Reference to the fiber
       */
      const std::vector<LINALG::Matrix<3, 1>>& GetFiber(FiberType type) const;

      /*!
       * \brief Returns a reference to a list of fibers belonging to the specific type for each
       * Gauss point
       *
       * @param type Type of the fiber
       * @return Reference to the fiber
       */
      std::vector<LINALG::Matrix<3, 1>>& GetFiberMutual(FiberType type);

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
      bool ContainsFiber(FiberType type) const;

      /**
       * Returns the number of fibers stored in the fiber holder
       * @return number of fibers
       */
      unsigned int FibersSize() const;

      /**
       * Returns the number of angles stored in the fiber holder
       * @return number of fibers
       */
      unsigned int AnglesSize() const;

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
        os << "NodalFiberHolder {Fibers count=" << holder.FibersSize()
           << ", Angles count=" << holder.AnglesSize() << "}";
        return os;
      }

      /*!
       * \brief Returns true of the map of both NodalFiberHolders are the same
       * @param holder1 reference to the first holder
       * @param holder2 reference to the second holder
       * @return
       */
      friend bool operator==(
          const DRT::FIBER::NodalFiberHolder& holder1, const DRT::FIBER::NodalFiberHolder& holder2)
      {
        return holder1.fibers_ == holder2.fibers_ && holder1.angles_ == holder2.angles_;
      }

     private:
      /// Map holding GP fibers of different names
      std::map<FiberType, std::vector<LINALG::Matrix<3, 1>>> fibers_;

      /// Map holding GP angles of different names
      std::map<AngleType, std::vector<double>> angles_;
    };
  }  // namespace FIBER
}  // namespace DRT
#endif  // BACI_NODAL_FIBER_HOLDER_H
