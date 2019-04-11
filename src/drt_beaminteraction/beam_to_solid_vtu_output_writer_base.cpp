/*!
\file beam_to_solid_vtu_output_writer_base.cpp

\brief Base object that stores all relevant data for beam to solid output

\level 3

\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_vtu_output_writer_base.H"

#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidVtuOutputWriterBase::BeamToSolidVtuOutputWriterBase(
    const std::string& base_output_name,
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params)
    : vtk_params_(vtk_params)
{
  base_output_name_ =
      DRT::Problem::Instance()->OutputControlFile()->FileNameOnlyPrefix() + "-" + base_output_name;
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>
BEAMINTERACTION::BeamToSolidVtuOutputWriterBase::AddVisualizationWriter(
    const std::string& writer_name)
{
  const auto& it = visualization_writers_.find(writer_name);
  if (it != visualization_writers_.end())
  {
    dserror("The output writer '%s' you want to add already exists.", writer_name.c_str());
    return Teuchos::null;
  }
  else
  {
    std::string writer_full_name = base_output_name_ + "-" + writer_name;
    Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization> new_writer =
        Teuchos::rcp<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>(
            new BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization(
                writer_full_name, vtk_params_));
    visualization_writers_[writer_name] = new_writer;
    return new_writer;
  }
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>
BEAMINTERACTION::BeamToSolidVtuOutputWriterBase::GetVisualizationWriter(
    const std::string& writer_name)
{
  const auto& it = visualization_writers_.find(writer_name);
  if (it != visualization_writers_.end())
    return it->second;
  else
    return Teuchos::null;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidVtuOutputWriterBase::Write(
    const unsigned int timestep_number, const double time)
{
  for (auto& it : visualization_writers_) it.second->Write(timestep_number, time);
  // TODO steinbrecher delete data here.
}
