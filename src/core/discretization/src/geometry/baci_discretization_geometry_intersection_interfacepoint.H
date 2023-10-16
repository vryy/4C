/*----------------------------------------------------------------------*/
/*! \file

\brief      InterfacePoint stores and delivers all data a point lying on the
            intersection interface has to know
\level 2

*/


#ifndef BACI_DISCRETIZATION_GEOMETRY_INTERSECTION_INTERFACEPOINT_H
#define BACI_DISCRETIZATION_GEOMETRY_INTERSECTION_INTERFACEPOINT_H


#include "baci_linalg_fixedsizematrix.H"


namespace CORE::GEO
{
  //! possible positions of a point with respect to an element
  enum pointType
  {
    NOTYPE,
    INTERNAL,
    SURFACE,
    LINE,
    NODE
  };


  /*!
  \brief      InterfacePoint stores and delivers all data a point lying on the
              intersection interface has to know
  */
  class InterfacePoint
  {
   public:
    /*!
    \brief constructor
    */
    InterfacePoint();


    /*!
    \brief constructor
    \param pType                (in)        : point type
    \param nodeId               (in)        : xfem node Ids the interface point is lying on
    \param lineId               (in)        : xfem line Ids ids the interface point is lying on
    \param surfId               (in)        : xfem surface Ids ids the interface point is lying on
    \param coordinates          (in)        : coordinates of the interface points
    */
    InterfacePoint(CORE::GEO::pointType& pType, int nodeId, std::vector<int>& lineId,
        std::vector<int>& surfId, CORE::LINALG::Matrix<3, 1>& coordinates);


    /*!
    \brief copy constructor
    */
    InterfacePoint(const CORE::GEO::InterfacePoint& old);


    //! destructor
    ~InterfacePoint();


    //! assignemnent operator
    InterfacePoint& operator=(const InterfacePoint& point);

    /*!
    \brief  set number of xfem nodes lines and surfaces the interface point is lying

    \param pType            (in)        : point Type of interface point
    */
    void setNodeLineSurfNumbers(const CORE::GEO::pointType pType);

    /*!
      \brief Return coordinates the interface points
     */
    inline CORE::LINALG::Matrix<3, 1> getCoord() const { return coord_; }


   private:
    //! pointType NOTYPE INTERNAL SURFACE LINE NODE
    CORE::GEO::pointType pType_;

    //! number of xfem nodes the point is lying on
    int nnode_;

    //! number of xfem lines the point is lying on
    int nline_;

    //! number of xfem surfaces the point is lying on
    int nsurf_;

    //! ids of xfem nodes the point is lying on
    int nodeId_;

    //! ids of xfem lines the point is lying on
    std::vector<int> lineId_;

    //! ids of xfem surfaces the point is lying on
    std::vector<int> surfId_;

    //! coordinates of interface points
    CORE::LINALG::Matrix<3, 1> coord_;
  };


}  // namespace CORE::GEO


#endif  // DISCRETIZATION_GEOMETRY_INTERSECTION_INTERFACEPOINT_H
