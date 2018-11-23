/*!
\file geometry_pair.cpp

\brief base class for two or more elements that interact geometrically, e.g. contact or meshtying
problems.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


#include "geometry_pair.H"

#include "../drt_lib/drt_dserror.H"


/**
 *
 */
GEOMETRYPAIR::GeometryPair::GeometryPair()
    : isinit_(false),
      issetup_(false),
      evaluation_data_(Teuchos::null),
      element1_(NULL),
      element2_(NULL)
{
  // Empty constructor.
}

/**
 *
 */
void GEOMETRYPAIR::GeometryPair::Init(
    Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> evaluation_data_ptr,
    const DRT::Element* element1, const DRT::Element* element2)
{
  evaluation_data_ = evaluation_data_ptr;
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

/**
 *
 */
void GEOMETRYPAIR::GeometryPair::CheckInit() const
{
  if (!isinit_) dserror("Init() has not been called yet!");
}

/**
 *
 */
void GEOMETRYPAIR::GeometryPair::CheckInitSetup() const
{
  if (!isinit_ || !issetup_) dserror("Init() and Setup() have not been called yet!");
}
