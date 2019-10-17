/*----------------------------------------------------------------------*/
/*! \file

\brief base class for two or more elements that interact geometrically, e.g. contact or meshtying
problems.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair.H"


/**
 *
 */
GEOMETRYPAIR::GeometryPair::GeometryPair()
    : isinit_(false), issetup_(false), element1_(nullptr), element2_(nullptr)
{
  // Empty constructor.
}

/**
 *
 */
void GEOMETRYPAIR::GeometryPair::Init(const DRT::Element* element1, const DRT::Element* element2)
{
  element1_ = element1;
  element2_ = element2;

  isinit_ = true;
}

/**
 *
 */
void GEOMETRYPAIR::GeometryPair::Setup()
{
  CheckInit();

  issetup_ = true;
}
