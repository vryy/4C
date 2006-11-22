/*!----------------------------------------------------------------------
\file utils.cpp
\brief A collection of helper methods for namespace CCADISCRETIZATION

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "utils.H"
#include "node.H"
#include "designnode.H"
#include "designelement.H"
#include "shell8.H"
#include "dserror.H"


/*----------------------------------------------------------------------*
 |  create an instance of ParObject  (public)                mwgee 11/06|
 *----------------------------------------------------------------------*/
CCADISCRETIZATION::ParObject* CCADISCRETIZATION::Factory(const char* data)
{
  // mv ptr behind the size record
  const int* ptr = (const int*)data;
  ptr++;
  
  // get the type
  const int type = *ptr;
  
  switch(type)
  {
    case ParObject_Container:
    {
      CCADISCRETIZATION::Container* object = new CCADISCRETIZATION::Container();
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Condition:
    {
      CCADISCRETIZATION::Condition* object = 
        new CCADISCRETIZATION::Condition(CCADISCRETIZATION::Condition::condition_none);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Node:
    {
      double dummycoord[3] = {999.,999.,999.};
      CCADISCRETIZATION::Node* object = new CCADISCRETIZATION::Node(-1,dummycoord,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_DesignNode:
    {
      double dummycoord[3] = {999.,999.,999.};
      CCADISCRETIZATION::DesignNode* object = new CCADISCRETIZATION::DesignNode(-1,dummycoord,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Element:
    {
      dserror("CCADISCRETIZATION::Element is pure virtual, cannot create instance");
    }
    break;
    case ParObject_DesignElement:
    {
      CCADISCRETIZATION::DesignElement* object = 
        new CCADISCRETIZATION::DesignElement(-1,CCADISCRETIZATION::Element::element_none,-1);
      object->Unpack(data);
      return object;
    }
    break;
    case ParObject_Shell8:
    {
      CCADISCRETIZATION::Shell8* object = new CCADISCRETIZATION::Shell8(-1,-1);
      object->Unpack(data);
      return object;
    }
    break;
    default:
      dserror("Unknown type of ParObject instance: %d",type);
    break;
  }
  
  return NULL;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
