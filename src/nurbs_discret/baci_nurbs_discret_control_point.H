/*----------------------------------------------------------------------*/
/*! \file

\brief   This is basically a (3d-) node with an additional weight.
   The weight is required for the evaluation of the nurbs
   basis functions.

   note that X() is not the coordinate of some grid point
   anymore, it's just the control point position


\level 2


*----------------------------------------------------------------------*/
#ifndef BACI_NURBS_DISCRET_CONTROL_POINT_H
#define BACI_NURBS_DISCRET_CONTROL_POINT_H

#include "baci_config.H"

#include "baci_comm_parobjectfactory.H"
#include "baci_lib_node.H"

BACI_NAMESPACE_OPEN

namespace DRT
{
  namespace NURBS
  {
    class ControlPointType : public CORE::COMM::ParObjectType
    {
     public:
      std::string Name() const override { return "ControlPointType"; }

      static ControlPointType& Instance() { return instance_; };

      CORE::COMM::ParObject* Create(const std::vector<char>& data) override;

     private:
      static ControlPointType instance_;
    };

    /*!
    \brief Control points for nurbs surfaces/volumes

    Control points are derived from nodes. For isogeometric
    analysis, they replace the nodes in the (nurbs)
    discretisation.
    They are connected to elements, have degrees of freedom
    etc just like normal nodes (and hence the discretisation
    can handle them correctly), the only difference is that
    their coordinate is just a control point coordinate and
    not a mesh point coordinate.

    */
    class ControlPoint : public DRT::Node
    {
     public:
      //! @name Enums and Friends

      /*!
      \brief The Discretization is a friend of the control point
      */
      friend class Discretization;

      //@}

      //! @name Constructors and destructors and related methods

      /*!
      \brief Standard Constructor

      \param id     (in): A globally unique control point id
      \param coords (in): vector of nodal coordinates, length 3
      \param weight (in): nurbs weight
      \param owner  (in): Owner of this node.
      */
      ControlPoint(int id, const std::vector<double>& coords, const double weight, const int owner);

      /*!
      \brief Copy Constructor

      Makes a deep copy of a control point

      \param old (in): The control point to copy

      */
      ControlPoint(const DRT::NURBS::ControlPoint& old);

      /*!
      \brief Deep copy the derived class and return
             pointer to it

      */
      DRT::NURBS::ControlPoint* Clone() const override;


      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique
      id defined at the top of parobject.H.

      \return the parobject id
      */
      int UniqueParObjectId() const override
      {
        return ControlPointType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this CP

      \param data (in/out): a char vector to pack the data into

      */
      void Pack(CORE::COMM::PackBuffer& data) const override;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this CP

      \param data (in): a char vector to unpack the data from

      */
      void Unpack(const std::vector<char>& data) override;

      //@}

      //! @name Additional access methods

      /*!
      \brief Return weight

      \return weight

      */
      virtual inline double W() const { return w_; }

      /*!
      \brief Print this node

      \param os ofstrem
      */
      void Print(std::ostream& os) const override;

      //@}

     protected:
      //! nurbs weight
      double w_;

    };  // class ControlPoint

  }  // namespace NURBS

}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif  // BACI_NURBS_DISCRET_CONTROL_POINT_H
