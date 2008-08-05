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


namespace GEO
{

//! possible positions of a point with respect to an element
enum ObjectType {NOTYPE, SURFACE, LINE, NODE};


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
    \brief  set node object type
    \param nodeId       (in)        : node gid
    \param label        (in)        : label
    */
    void setNodeObjectType(
          const int nodeId, 
          const int label
          );
   
    /*!
    \brief  set line object type 
    \param lineId       (in)        : line gid
    \param surfId       (in)        : surf gid
    \param label        (in)        : label
    */
    void setLineObjectType(
          const int lineId,
          const int surfId, 
          const int label
          );
    
    /*!
    \brief  set surface object type 
    \param surfId       (in)        : surf gid
    \param label        (in)        : label
    */
    void setSurfaceObjectType(
          const int surfId, 
          const int label
          );  
    
    /*!
    \brief  set node gId of object
    \param nodeId            (in)        : node gId of object
    */
    void setNodeId(
      const int    nodeId
      ); 
       
    /*!
    \brief  set line Id of object 
            (note:: lineId is a local id with respect to a surface element, 
            whose id has to be provided in addition)
    \param lineId            (in)        : line lId of object of a surface element
    \param surfId            (in)        : surf gId of object 
    */
    void setLineId(
      const int    lineId,
      const int    surfId
      );  
         
    /*!
    \brief  set surf gId of object
    \param surfId            (in)        : surf gId of object
    */
    void setSurfaceId(
      const int    surfId
      );      
    
    /*!
    \brief  set label
    \param label           (in)        : label
    */
    void setLabel(
      const int    label
      );      
    
    /*!
    \brief Return object type 
     */
    inline GEO::ObjectType getObjectType() const { return objectType_;}
           
    /*!
    \brief Return node Id of object
     */
    inline int getNodeId() const { 
      if(objectType_ != NODE || nodeId_ == -1)
        dserror("objecttype is not correct (nodeId) or nodeId not set");
      return nodeId_;}   
    
    /*!
    \brief Return line Id of object
     */
    inline int getLineId() const { 
      if(objectType_ != LINE || lineId_ == -1)
        dserror("objecttype is not correct (lineId) or line id not set");
      return lineId_;}
       
    /*!
    \brief Return surface Id of object
     */
    inline int getSurfaceId() const { 
      if(objectType_ != LINE || objectType_ != SURFACE || lineId_ == -1 || surfId_ == -1 )
        dserror("objecttype is not correct (surfeId) or ids are not set");
      return surfId_;}
       
    /*!
    \brief Return label
     */
    inline int getLabel() const { return label_;}
    
    

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

};

}

#endif /*OCTTREE_NEARESTOBJECT_H_*/
#endif /*CCADISCRET*/ 


