/*----------------------------------------------------------------------*/
/*! \file

\brief base class for geometry pair evaluation data subcontainers.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_evaluation_data_base.H"

#include "../drt_lib/drt_dserror.H"


/**
 *
 */
GEOMETRYPAIR::GeometryEvaluationDataBase::GeometryEvaluationDataBase()
    : isinit_(false), issetup_(false)
{
  // Empty constructor.
}

/**
 *
 */
void GEOMETRYPAIR::GeometryEvaluationDataBase::Init()
{
  // The flag isinit_ will be set in the derived method.
}

/**
 *
 */
void GEOMETRYPAIR::GeometryEvaluationDataBase::Setup()
{
  CheckInit();

  // The flag issetup_ will be set in the derived method.
}

/**
 *
 */
void GEOMETRYPAIR::GeometryEvaluationDataBase::CheckInit() const
{
  if (!isinit_) dserror("Init() has not been called yet!");
}

/**
 *
 */
void GEOMETRYPAIR::GeometryEvaluationDataBase::CheckInitSetup() const
{
  if (!isinit_ || !issetup_) dserror("Init() and Setup() have not been called yet!");
}
