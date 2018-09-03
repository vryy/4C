/*!----------------------------------------------------------------------
\file drt_control_point.H

\brief   This is basically a (3d-) node with an additional weight.
   The weight is required for the evaluation of the nurbs
   basis functions.

   note that X() is not the coordinate of some grid point
   anymore, it's just the control point position

\maintainer Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271

\level 2


*----------------------------------------------------------------------*/
#ifndef CONTROL_POINT_H
#define CONTROL_POINT_H

#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_parobjectfactory.H"

/*!
\brief DRT: namespace of the ccarat discretization module

*/
namespace DRT
{
  /*!
  \brief NURBS: namespace of all isogeometric analysis tools

  */
  namespace NURBS
  {
    class ControlPointType : public DRT::ParObjectType
    {
     public:
      std::string Name() const { return "ControlPointType"; }

      static ControlPointType& Instance() { return instance_; };

      virtual DRT::ParObject* Create(const std::vector<char>& data);

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

    \author gammi (gamnitzer@lnm.mw.tum.de)
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
      ControlPoint(int id, const double* coords, const double weight, const int owner);

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
      virtual DRT::NURBS::ControlPoint* Clone() const;

      /*!
      \brief Destructor

      */
      virtual ~ControlPoint();

      /*!
      \brief Return unique ParObject id

      every class implementing ParObject needs a unique
      id defined at the top of drt_parobject.H.

      \return the parobject id
      */
      virtual int UniqueParObjectId() const
      {
        return ControlPointType::Instance().UniqueParObjectId();
      }

      /*!
      \brief Pack this class so it can be communicated

      \ref Pack and \ref Unpack are used to communicate this CP

      \param data (in/out): a char vector to pack the data into

      */
      virtual void Pack(DRT::PackBuffer& data) const;

      /*!
      \brief Unpack data from a char vector into this class

      \ref Pack and \ref Unpack are used to communicate this CP

      \param data (in): a char vector to unpack the data from

      */
      virtual void Unpack(const std::vector<char>& data);

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
      virtual void Print(std::ostream& os) const;

      //@}

     protected:
      //! nurbs weight
      double w_;

    };  // class ControlPoint

  }  // namespace NURBS

}  // namespace DRT

#endif  // #ifndef CONTROL_POINT_H
