/*----------------------------------------------------------------------*/
/*! \file

\brief This file defines a class that holds different nodal fibers

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_NODAL_FIBER_HOLDER_H
#define BACI_NODAL_FIBER_HOLDER_H

#include <vector>
#include <unordered_map>
#include "../linalg/linalg_fixedsizematrix.H"

namespace DRT
{
  namespace FIBER
  {
    /*!
     * \brief Type of the fiber
     *
     * In future, the FiberType should differentiate between different kinds of fiber definition,
     * like the "FIBERX" notation, "RAD AXI CIR" notation, ...
     * The way nodal fibers are implemented now however does not allow this differentiation.
     */
    enum class FiberType
    {
      Fiber1,
      Fiber2
    };

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

      /**
       * Returns the number of fibers stored in the fiber holder
       * @return number of fibers
       */
      unsigned int Size() const;

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
        os << "NodalFiberHolder {Count=" << holder.Size() << "}";
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
        return holder1.map_ == holder2.map_;
      }

     private:
      /// Map holding GP fibers of different names
      std::unordered_map<FiberType, const std::vector<LINALG::Matrix<3, 1>>> map_;
    };
  }  // namespace FIBER
}  // namespace DRT
#endif  // BACI_NODAL_FIBER_HOLDER_H
