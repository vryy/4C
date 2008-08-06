/*!
\file octtree_nearestobject.H

\brief stores data about nearest object in oct tree

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
 */

#ifdef CCADISCRET
#ifndef OCTTREE_NEARESTOBJECT_H_
#define OCTTREE_NEARESTOBJECT_H_

#include "../drt_lib/drt_dserror.H"
#include "../drt_geometry/vector_definitions.H"


namespace GEO
{

//! possible positions of a point with respect to an element
enum ObjectType {NOTYPE_OBJECT, SURFACE_OBJECT, LINE_OBJECT, NODE_OBJECT};


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
    \brief Set node object type
    \param nodeId       (in)        : node gid
    \param label        (in)        : label
    \param physcoord    (in)        : physical coordinates of point on object
    */
    void setNodeObjectType(
          const int nodeId, 
          const int label,
          const BlitzVec3& physcoord
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
          const BlitzVec3& physcoord
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
          const BlitzVec3& physcoord
          );  
    
    /*!
    \brief Return object type 
     */
    inline GEO::ObjectType getObjectType() const { return objectType_;}
           
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
    inline BlitzVec3 getPhysCoord() const {
      if(objectType_ == NOTYPE_OBJECT)
        dserror("no object type and physical coordinates are set");
      return physcoord_;}
    
    

  private:
  
    //! ObjectType NOTYPE SURFACE LINE NODE 
    GEO::ObjectType       objectType_;
  
    //! id of node
    int                   nodeId_;
  
    //! id of line 
    int                   lineId_;
    
    //! id of surface
    int                   surfId_;
    
    //! label of object
    int                   label_;
    
    //! physical coordinates of point on nearest object
    BlitzVec3             physcoord_;

};

}

#endif /*OCTTREE_NEARESTOBJECT_H_*/
#endif /*CCADISCRET*/ 


