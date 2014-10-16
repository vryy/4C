/*!----------------------------------------------------------------------
\file beam3contact_utils.cpp
\brief A set of utility functions for beam contact

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

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
Maintainer: Christoh Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#include "beam3contact_utils.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "beam3contact_manager.H"


//Cast of FAD to double
double BEAMCONTACT::CastToDouble(FAD a)
{
  return a.val();
}

//Cast of double to double
double BEAMCONTACT::CastToDouble(double a)
{
  return a;
}

//Calculate Norm of a scalar FAD or double quantity
double BEAMCONTACT::Norm(double a)
{
  return sqrt(a*a);
}

//Calculate Norm of a scalar FAD or double quantity
FAD BEAMCONTACT::Norm(FAD a)
{
  return pow(a*a,0.5);
}
/*----------------------------------------------------------------------*
 |  Check, if current node belongs to a beam element         meier 05/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::BeamNode(DRT::Node& node)
{
  bool beameles = false;
  bool othereles = false;

  //TODO: actually we would have to check all elements of all processors!!! Gather?
  for (int i=0; i< (int)(node.NumElement()); i++)
  {
    if(BeamElement(*(node.Elements())[i]))
      beameles = true;
    else
      othereles = true;
  }

  if (beameles and othereles)
    dserror("Beam elements and other (solid, rigid sphere) elements sharing the same node is currently not allowed in BACI!");

  return beameles;
}

/*----------------------------------------------------------------------*
 |  Check, if current node belongs to a rigid sphere element   grill 09/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::RigidsphereNode(DRT::Node& node)
{
  bool sphereeles = false;
  bool othereles = false;

  //TODO: actually we would have to check all elements of all processors!!! Gather?
  for (int i=0; i< (int)(node.NumElement()); i++)
  {
    if(RigidsphereElement(*(node.Elements())[i]))
      sphereeles = true;
    else
      othereles = true;
  }

  if (sphereeles and othereles)
    dserror("Rigid sphere elements and other (solid, beam) elements sharing the same node is currently not allowed in BACI!");

  return sphereeles;
}

/*----------------------------------------------------------------------*
 |  Check, if current element is a beam element         meier 05/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::BeamElement(DRT::Element& element)
{
  const DRT::ElementType& ele_type = element.ElementType();

  if (ele_type == DRT::ELEMENTS::Beam3ebType::Instance() or ele_type ==DRT::ELEMENTS::Beam3Type::Instance() or ele_type ==DRT::ELEMENTS::Beam3iiType::Instance())
    return true; //TODO: Print Warning, that only these three types of beam elements are supported!!!
  else
    return false;
}

/*----------------------------------------------------------------------*
 |  Check, if current element is a rigid sphere element       grill 09/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::RigidsphereElement(DRT::Element& element)
{
  const DRT::ElementType& ele_type = element.ElementType();

  if (ele_type == DRT::ELEMENTS::RigidsphereType::Instance())
    return true;
  else
    return false;
}

/*----------------------------------------------------------------------*
 |  Check, if two elements share a node -> neighbor elements meier 05/14|
 *----------------------------------------------------------------------*/
bool BEAMCONTACT::ElementsShareNode(DRT::Element& element1,DRT::Element& element2)
{
  bool sharenode = false;

  for (int i=0;i<element1.NumNode();i++)
  {
    int id = element1.NodeIds()[i];

    for (int j=0;j<element2.NumNode();j++)
    {
      if(id==element2.NodeIds()[j])
        sharenode = true;
    }
  }

    return sharenode;
}
