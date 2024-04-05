/*----------------------------------------------------------------------*/
/*! \file

\brief   This is basically a (3d-) node with an additional fiber direction.

\level 2
*----------------------------------------------------------------------*/
#ifndef FOUR_C_FIBER_NODE_HPP
#define FOUR_C_FIBER_NODE_HPP

#include "baci_config.hpp"

#include "baci_comm_parobjectfactory.hpp"
#include "baci_lib_node.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace FIBER
  {
    /*!
     * \brief Type of the coordinate system direction
     */
    enum class CoordinateSystemDirection
    {
      Circular,
      Tangential,
      Radial
    };

    /*!
     * \brief Type of the Angle
     */
    enum class AngleType
    {
      Helix,
      Transverse
    };

    class FiberNodeType : public CORE::COMM::ParObjectType
    {
     public:
      std::string Name() const override { return "FiberNodeType"; }

      static FiberNodeType& Instance() { return instance_; };

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

     private:
      static FiberNodeType instance_;
    };

    /*!
    \brief Node with additional fiber information

    */
    class FiberNode : public DRT::Node
    {
     public:
      //! @name Enums and Friends

      /*!
      \brief The Discretization is a friend of the fiber node
      */
      friend class Discretization;

      //@}

      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id     (in): A globally unique fiber node
      \param coords (in): vector of nodal coordinates, length 3
      \param coordinateSystemDirections (in): map of coordinate system directions
      \param fibers (in): list of fiber vectors
      \param angles (in): map of angles
      \param owner  (in): Owner of this node.
      */
      FiberNode(int id, const std::vector<double>& coords,
          std::map<FIBER::CoordinateSystemDirection, std::array<double, 3>>
              coordinateSystemDirections,
          std::vector<std::array<double, 3>> fibers, std::map<FIBER::AngleType, double> angles,
          const int owner);

      FiberNode(const DRT::FIBER::FiberNode& old) = default;

      /*!
      \brief Deep copy the derived class and return
             pointer to it
      */
      DRT::FIBER::FiberNode* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique
      id defined at the top of parobject.H.

      \return the parobject id
      */
      int UniqueParObjectId() const override
      {
        return FiberNodeType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this node

      \param data (in/out): a char vector to pack the data into

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this node

      \param data (in): a char vector to unpack the data from

      */
      void Unpack(const std::vector<char>& data) override;

      //@}

      inline const std::map<CoordinateSystemDirection, std::array<double, 3>>&
      CoordinateSystemDirections() const
      {
        return coordinateSystemDirections_;
      }

      inline const std::vector<std::array<double, 3>>& Fibers() const { return fibers_; }

      inline const std::map<AngleType, double>& Angles() const { return angles_; }

      /*!
      \brief Print this node

      \param os ofstrem
      */
      void Print(std::ostream& os) const override;

      //@}

     protected:
      /// map of coordinate system directions
      std::map<CoordinateSystemDirection, std::array<double, 3>> coordinateSystemDirections_;

      /// list of fibers
      std::vector<std::array<double, 3>> fibers_;

      /// map of angles
      std::map<AngleType, double> angles_;

    };  // class FiberNode

  }  // namespace FIBER

}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // FIBER_NODE_H
