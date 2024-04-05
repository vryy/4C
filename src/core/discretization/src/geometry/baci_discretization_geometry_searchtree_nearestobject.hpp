/*----------------------------------------------------------------------*/
/*! \file

\brief stores data about nearest object in search tree

\level 3


*----------------------------------------------------------------------*/
#ifndef FOUR_C_DISCRETIZATION_GEOMETRY_SEARCHTREE_NEARESTOBJECT_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRY_SEARCHTREE_NEARESTOBJECT_HPP


#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"

BACI_NAMESPACE_OPEN


namespace CORE::GEO
{
  //! possible positions of a point with respect to an element
  enum ObjectType
  {
    NOTYPE_OBJECT,   ///< closest object not defined
    SURFACE_OBJECT,  ///< closest object is a point
    LINE_OBJECT,     ///< closest object is a line
    NODE_OBJECT      ///< closest object is a surface
  };


  /*!
  \brief  NearestObject stores and delivers all data , which is important
          during a nearest object in tree node search
  */
  class NearestObject
  {
   public:
    /*!
    \brief constructor
    */
    NearestObject();

    /*!
    \brief copy constructor
    */
    NearestObject(const CORE::GEO::NearestObject& old);

    /*!
    \brief assignment operator
    */
    CORE::GEO::NearestObject& operator=(const CORE::GEO::NearestObject& old);

    /*!
    \brief clear nearest object
    */
    void clear();

    /*!
    \brief Set node object type
    \param nodeId       (in)        : node gid
    \param label        (in)        : label
    \param physcoord    (in)        : physical coordinates of point on object
    */
    void setNodeObjectType(
        const int nodeId, const int label, const CORE::LINALG::Matrix<3, 1>& physcoord);

    /*!
    \brief Set line object type
    \param lineId       (in)        : line gid
    \param surfId       (in)        : surf gid
    \param label        (in)        : label
    \param physcoord    (in)        : physical coordinates of point on object
    */
    void setLineObjectType(const int lineId, const int surfId, const int label,
        const CORE::LINALG::Matrix<3, 1>& physcoord);

    /*!
    \brief Set surface object type
    \param surfId       (in)        : surf gid
    \param label        (in)        : label
    \param physcoord    (in)        : physical coordinates of point on object
    */
    void setSurfaceObjectType(
        const int surfId, const int label, const CORE::LINALG::Matrix<3, 1>& physcoord);

    /*!
    \brief Return object type
     */
    inline ObjectType getObjectType() const { return objectType_; }

    /*!
    \brief Return label
     */
    inline int getLabel() const { return label_; }

    /*!
    \brief Return vector of physical coordinates
     */
    inline CORE::LINALG::Matrix<3, 1> getPhysCoord() const
    {
      if (objectType_ == NOTYPE_OBJECT) dserror("no object type and physical coordinates are set");
      return physcoord_;
    }


   private:
    //! ObjectType NOTYPE SURFACE LINE NODE
    ObjectType objectType_;

    //! id of node
    int nodeId_;

    //! id of line
    int lineId_;

    //! id of surface
    int surfId_;

    //! label of object
    int label_;

    //! physical coordinates of point on nearest object
    CORE::LINALG::Matrix<3, 1> physcoord_;
  };

}  // namespace CORE::GEO

BACI_NAMESPACE_CLOSE

#endif  // DISCRETIZATION_GEOMETRY_SEARCHTREE_NEARESTOBJECT_H
