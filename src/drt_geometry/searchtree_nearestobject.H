/*!
\file searchtree_nearestobject.H

\brief stores data about nearest object in search tree

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions and of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>
 */


#ifndef SEARCHTREE_NEARESTOBJECT_H_
#define SEARCHTREE_NEARESTOBJECT_H_


#include "../linalg/linalg_fixedsizematrix.H"


namespace GEO
{

//! possible positions of a point with respect to an element
enum ObjectType {
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
    NearestObject(
        const GEO::NearestObject& old
        );

    /*!
    \brief assignment operator
    */
    GEO::NearestObject& operator=(const GEO::NearestObject& old);


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
          const int nodeId,
          const int label,
          const LINALG::Matrix<3,1>& physcoord
          );

    /*!
    \brief Set line object type
    \param lineId       (in)        : line gid
    \param surfId       (in)        : surf gid
    \param label        (in)        : label
    \param physcoord    (in)        : physical coordinates of point on object
    */
    void setLineObjectType(
          const int lineId,
          const int surfId,
          const int label,
          const LINALG::Matrix<3,1>& physcoord
          );

    /*!
    \brief Set surface object type
    \param surfId       (in)        : surf gid
    \param label        (in)        : label
    \param physcoord    (in)        : physical coordinates of point on object
    */
    void setSurfaceObjectType(
          const int surfId,
          const int label,
          const LINALG::Matrix<3,1>& physcoord
          );

    /*!
    \brief Return object type
     */
    inline ObjectType getObjectType() const { return objectType_;}

    /*!
    \brief Return node Id of object
     */
    inline int getNodeId() const {
      if(objectType_ != NODE_OBJECT || nodeId_ == -1)
        dserror("objecttype is not correct (nodeId) or nodeId not set");
      return nodeId_;}

    /*!
    \brief Return line Id of object
     */
    inline int getLineId() const {
      if(objectType_ != LINE_OBJECT || lineId_ == -1)
        dserror("objecttype is not correct (lineId) or line id not set");
      return lineId_;}

    /*!
    \brief Return surface Id of object
     */
    inline int getSurfaceId() const {
      if( (objectType_ != LINE_OBJECT && objectType_ != SURFACE_OBJECT) || ( lineId_ == -1 && surfId_ == -1)  )
        dserror("objecttype is not correct (surfId) or ids are not set");
      return surfId_;}

    /*!
    \brief Return label
     */
    inline int getLabel() const { return label_;}

    /*!
    \brief Return vector of physical coordinates
     */
    inline LINALG::Matrix<3,1> getPhysCoord() const {
      if(objectType_ == NOTYPE_OBJECT)
        dserror("no object type and physical coordinates are set");
      return physcoord_;}



  private:

    //! ObjectType NOTYPE SURFACE LINE NODE
    ObjectType                        objectType_;

    //! id of node
    int                               nodeId_;

    //! id of line
    int                               lineId_;

    //! id of surface
    int                               surfId_;

    //! label of object
    int                               label_;

    //! physical coordinates of point on nearest object
    LINALG::Matrix<3,1>               physcoord_;

};

}

#endif /*OCTTREE_NEARESTOBJECT_H_*/

