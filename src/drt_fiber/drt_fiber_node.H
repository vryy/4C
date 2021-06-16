/*----------------------------------------------------------------------*/
/*! \file

\brief   This is basically a (3d-) node with an additional fiber direction.

\level 2
*----------------------------------------------------------------------*/
#ifndef FIBER_NODE_H
#define FIBER_NODE_H

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_parobjectfactory.H"

namespace DRT
{
  namespace FIBER
  {
    /*!
     * \brief Type of the fiber
     */
    enum FiberType
    {
      Fiber1,
      Fiber2,
      Fiber3,
      Circular,
      Tangential,
      Radial
    };

    /*!
     * \brief Type of the Angle
     */
    enum AngleType
    {
      Helix,
      Transverse
    };
    class FiberNodeType : public DRT::ParObjectType
    {
     public:
      std::string Name() const override { return "FiberNodeType"; }

      static FiberNodeType& Instance() { return instance_; };

      DRT::ParObject* Create(const std::vector<char>& data) override;

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
      \param fibers (in): unordered map of fiber vectors
      \param angles (in): unordered map of angles
      \param owner  (in): Owner of this node.
      */
      FiberNode(int id, std::array<double, 3> coords,
          std::unordered_map<FIBER::FiberType, std::array<double, 3>> fibers,
          std::unordered_map<FIBER::AngleType, double> angles, const int owner);

      FiberNode(const DRT::FIBER::FiberNode& old) = default;

      /*!
      \brief Deep copy the derived class and return
             pointer to it
      */
      DRT::FIBER::FiberNode* Clone() const override;

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique
      id defined at the top of drt_parobject.H.

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
      void Pack(DRT::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this node

      \param data (in): a char vector to unpack the data from

      */
      void Unpack(const std::vector<char>& data) override;

      //@}

      inline const std::unordered_map<FiberType, std::array<double, 3>>& Fibers() const
      {
        return fibers_;
      }

      inline const std::unordered_map<AngleType, double>& Angles() const { return angles_; }

      /*!
      \brief Print this node

      \param os ofstrem
      */
      void Print(std::ostream& os) const override;

      //@}

     protected:
      /// Unordered map of fibers
      std::unordered_map<FiberType, std::array<double, 3>> fibers_;

      /// Unordered map of angles
      std::unordered_map<AngleType, double> angles_;

    };  // class FiberNode

  }  // namespace FIBER

}  // namespace DRT

#endif  // #ifndef FIBER_NODE_H
