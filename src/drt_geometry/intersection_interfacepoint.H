/*!
\file intersection_interfacepoint.H

\brief      InterfacePoint stores and delivers all data a point lying on the
            intersection interface has to know
\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>
*/


#ifndef INTERSECTION_INTERFACEPOINT_H_
#define INTERSECTION_INTERFACEPOINT_H_


#include "../linalg/linalg_fixedsizematrix.H"


namespace GEO
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
    InterfacePoint(GEO::pointType& pType, int nodeId, std::vector<int>& lineId,
        std::vector<int>& surfId, LINALG::Matrix<3, 1>& coordinates);


    /*!
    \brief copy constructor
    */
    InterfacePoint(const GEO::InterfacePoint& old);


    //! destructor
    ~InterfacePoint();


    //! assignemnent operator
    InterfacePoint& operator=(const InterfacePoint& point);


    /*!
    \brief  set number of xfem nodes lines and surfaces the interface point is lying

    \param pType            (in)        : point Type of interface point
    */
    void setNodeLineSurfNumbers(const GEO::pointType pType);


    /*!
    \brief  set xfem node Ids the interface points is lying on

    \param pType            (in)        : point Type of interface point
    */
    void setPointType(const GEO::pointType pType);

    /*!
    \brief  set xfem node Ids the interface points is lying on

    \param nodeId            (in)        : node Id the interface point is lying on
    */
    void setNodeId(const int nodeId);

    /*!
    \brief  set xfem line Ids the interface points is lying on

    \param lineId            (in)        : vector of xfem line Ids the interface point is lying on
    */
    void setLineId(const std::vector<int>& lineId);

    /*!
    \brief  set xfem surface Ids the interface points is lying on

    \param surfId            (in)        : vector of xfem surface Ids the interface point is lying
    on
    */
    void setSurfaceId(const std::vector<int>& surfId);

    /*!
    \brief  set coordinates of the interface point

    \param coordinates            (in)        : coordinates of interface point
    */
    void setCoord(const LINALG::Matrix<3, 1>& coordinates);

    /*!
    \brief  set X-coordinate of the interface point

    \param coordX            (in)        : X-coordinate of interface point
    */
    void setCoordX(const double coordX);

    /*!
    \brief  set Y-coordinate of the interface point

    \param coordY            (in)        : Y-coordinate of interface point
    */
    void setCoordY(const double coordY);

    /*!
    \brief  set Z-coordinate of the interface point

    \param coordZ            (in)        : Z-coordinate of interface point
    */
    void setCoordZ(const double coordZ);

    /*!
    \brief  set single coordinate of the interface point

    \param index            (in)        : index in coordinate vector
    \param coord            (in)        : coordinate of interface point
    */
    void setSingleCoord(const int index, const double coord);

    /*!
      \brief Return number of xfem nodes the interface point is lying on
    */
    inline int getNumNode() const { return nnode_; }

    /*!
      \brief Return number of xfem lines the interface point is lying on
    */
    inline int getNumLine() const { return nline_; }

    /*!
      \brief Return number of xfem surface the interface point is lying on
    */
    inline int getNumSurface() const { return nsurf_; }

    /*!
    \brief Return Id of xfem nodes the interface point is lying on
    */
    inline int getNodeId() const { return nodeId_; }

    /*!
    \brief Return Ids of xfem lines the interface point is lying on
     */
    inline std::vector<int> getLineId() const { return lineId_; }

    /*!
    \brief Return Ids of xfem surface the interface point is lying on
     */
    inline std::vector<int> getSurfId() const { return surfId_; }

    /*!
      \brief Return coordinates the interface points
     */
    inline LINALG::Matrix<3, 1> getCoord() const { return coord_; }

    /*!
    \brief Return point type the interface point ( NOTYPE INTERNAL SURFACE LINE NODE  )
     */
    inline GEO::pointType getPointType() const { return pType_; }



   private:
    //! pointType NOTYPE INTERNAL SURFACE LINE NODE
    GEO::pointType pType_;

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
    LINALG::Matrix<3, 1> coord_;
  };


}  // namespace GEO


#endif /*INTERSECTION_INTERFACEPOINT_H_*/
