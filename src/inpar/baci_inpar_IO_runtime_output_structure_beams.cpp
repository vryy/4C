/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of beam discretization at runtime

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_inpar_IO_runtime_output_structure_beams.hpp"

#include "baci_inpar.hpp"
#include "baci_inpar_parameterlist_utils.hpp"
#include "baci_inpar_validparameters.hpp"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN


namespace INPAR
{
  namespace IO_RUNTIME_OUTPUT
  {
    namespace BEAMS
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
      {
        using namespace INPUT;
        using Teuchos::setStringToIntegralParameter;
        using Teuchos::tuple;

        // related sublist
        Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
        Teuchos::ParameterList& sublist_IO_VTK_structure =
            sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");
        Teuchos::ParameterList& sublist_IO_output_beams =
            sublist_IO_VTK_structure.sublist("BEAMS", false, "");

        // whether to write special output for beam elements
        BoolParameter("OUTPUT_BEAMS", "No", "write special output for beam elements",
            &sublist_IO_output_beams);

        // whether to write displacement state
        BoolParameter("DISPLACEMENT", "No", "write displacement output", &sublist_IO_output_beams);

        // use absolute positions or initial positions for vtu geometry (i.e. point coordinates)
        // 'absolute positions' requires writing geometry in every output step (default for now)
        BoolParameter("USE_ABSOLUTE_POSITIONS", "Yes",
            "use absolute positions or initial positions for vtu geometry (i.e. point coordinates)",
            &sublist_IO_output_beams);

        // write internal (elastic) energy of element
        BoolParameter("INTERNAL_ENERGY_ELEMENT", "No",
            "write internal (elastic) energy for each element", &sublist_IO_output_beams);

        // write kinetic energy of element
        BoolParameter("KINETIC_ENERGY_ELEMENT", "No", "write kinetic energy for each element",
            &sublist_IO_output_beams);

        // write triads as three orthonormal base vectors at every visualization point
        BoolParameter("TRIAD_VISUALIZATIONPOINT", "No", "write triads at every visualization point",
            &sublist_IO_output_beams);

        // write material cross-section strains at the Gauss points:
        // axial & shear strains, twist & curvatures
        BoolParameter("STRAINS_GAUSSPOINT", "No",
            "write material cross-section strains at the Gauss points", &sublist_IO_output_beams);

        // write material cross-section strains at the visualization points:
        // axial & shear strains, twist & curvatures
        BoolParameter("STRAINS_CONTINUOUS", "No",
            "write material cross-section strains at the visualization points",
            &sublist_IO_output_beams);

        // write material cross-section stresses at the Gauss points:
        // axial and shear forces, torque and bending moments
        BoolParameter("MATERIAL_FORCES_GAUSSPOINT", "No",
            "write material cross-section stresses at the Gauss points", &sublist_IO_output_beams);

        // write material cross-section stresses at the visualization points:
        // axial and shear forces, torque and bending moments
        BoolParameter("MATERIAL_FORCES_CONTINUOUS", "No",
            "write material cross-section stresses at the visualization points",
            &sublist_IO_output_beams);

        // write spatial cross-section stresses at the Gauss points:
        // axial and shear forces, torque and bending moments
        BoolParameter("SPATIAL_FORCES_GAUSSPOINT", "No",
            "write material cross-section stresses at the Gauss points", &sublist_IO_output_beams);

        // write element filament numbers and type
        BoolParameter("BEAMFILAMENTCONDITION", "No", "write element filament numbers",
            &sublist_IO_output_beams);

        // write element and network orientation parameter
        BoolParameter("ORIENTATION_PARAMETER", "No", "write element filament numbers",
            &sublist_IO_output_beams);

        // write crossection forces of periodic RVE
        BoolParameter("RVE_CROSSSECTION_FORCES", "No", " get sum of all internal forces of  ",
            &sublist_IO_output_beams);

        // write reference length of beams
        BoolParameter(
            "REF_LENGTH", "No", "write reference length of all beams", &sublist_IO_output_beams);

        // write element GIDs
        BoolParameter(
            "ELEMENT_GID", "No", "write the BACI internal element GIDs", &sublist_IO_output_beams);

        // write element ghosting information
        BoolParameter("ELEMENT_GHOSTING", "No", "write which processors ghost the elements",
            &sublist_IO_output_beams);

        // number of subsegments along a single beam element for visualization
        IntParameter("NUMBER_SUBSEGMENTS", 5,
            "Number of subsegments along a single beam element for visualization",
            &sublist_IO_output_beams);
      }
    }  // namespace BEAMS
  }    // namespace IO_RUNTIME_OUTPUT
}  // namespace INPAR

BACI_NAMESPACE_CLOSE
