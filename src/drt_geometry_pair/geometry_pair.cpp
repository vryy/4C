/*!

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
