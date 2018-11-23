/*!
\file geometry_pair_evaluation_data_base.cpp

\brief Data container holding pointers to all subcontainers that in turn hold all evaluation data
specific to their problem type, that can not be stored pairwise.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


#include "geometry_pair_evaluation_data_global.H"
#include "geometry_pair_line_to_volume_evaluation_data.H"


/**
 *
 */
GEOMETRYPAIR::GeometryEvaluationDataGlobal::GeometryEvaluationDataGlobal()
    : line_to_volume_evaluation_data_(Teuchos::null), sub_container_list_(0)
{
  // Empty Constructor
}


/**
 *
 */
void GEOMETRYPAIR::GeometryEvaluationDataGlobal::ResetGeometryEvaluationDataGlobal()
{
  for (auto& sub_container : sub_container_list_)
  {
    sub_container->Init();
    sub_container->Setup();
  }
}


/**
 *
 */
void GEOMETRYPAIR::GeometryEvaluationDataGlobal::BuildLineToVolumeEvaluationData()
{
  line_to_volume_evaluation_data_ = Teuchos::rcp(new GEOMETRYPAIR::LineToVolumeEvaluationData());
  line_to_volume_evaluation_data_->Init();
  line_to_volume_evaluation_data_->Setup();
  sub_container_list_.push_back(line_to_volume_evaluation_data_);
}
