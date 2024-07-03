/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_inpar_IO_runtime_vtk_output_structure.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Inpar
{
  namespace IORuntimeOutput
  {
    namespace Solid
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
      {
        using namespace Input;
        using Teuchos::setStringToIntegralParameter;
        using Teuchos::tuple;

        // related sublist
        Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
        Teuchos::ParameterList& sublist_IO_VTK =
            sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");
        Teuchos::ParameterList& sublist_IO_VTK_structure =
            sublist_IO_VTK.sublist("STRUCTURE", false, "");

        // whether to write output for structure
        Core::UTILS::BoolParameter(
            "OUTPUT_STRUCTURE", "No", "write structure output", &sublist_IO_VTK_structure);

        // whether to write displacement state
        Core::UTILS::BoolParameter(
            "DISPLACEMENT", "No", "write displacement output", &sublist_IO_VTK_structure);

        // whether to write velocity state
        Core::UTILS::BoolParameter(
            "VELOCITY", "No", "write velocity output", &sublist_IO_VTK_structure);

        // whether to write element owner
        Core::UTILS::BoolParameter(
            "ELEMENT_OWNER", "No", "write element owner", &sublist_IO_VTK_structure);

        // whether to write element GIDs
        Core::UTILS::BoolParameter(
            "ELEMENT_GID", "No", "write 4C internal element GIDs", &sublist_IO_VTK_structure);

        // write element ghosting information
        Core::UTILS::BoolParameter("ELEMENT_GHOSTING", "No",
            "write which processors ghost the elements", &sublist_IO_VTK_structure);

        // whether to write node GIDs
        Core::UTILS::BoolParameter(
            "NODE_GID", "No", "write 4C internal node GIDs", &sublist_IO_VTK_structure);

        // whether to write stress and / or strain data
        Core::UTILS::BoolParameter("STRESS_STRAIN", "No",
            "Write element stress and / or strain  data. The type of stress / strain has to be "
            "selected in the --IO input section",
            &sublist_IO_VTK_structure);

        // mode to write gauss point data
        setStringToIntegralParameter<Inpar::Solid::GaussPointDataOutputType>(
            "GAUSS_POINT_DATA_OUTPUT_TYPE", "none",
            "Where to write gauss point data. (none, projected to nodes, projected to element "
            "center, raw at gauss points)",
            tuple<std::string>("none", "nodes", "element_center", "gauss_points"),
            tuple<Inpar::Solid::GaussPointDataOutputType>(
                Inpar::Solid::GaussPointDataOutputType::none,
                Inpar::Solid::GaussPointDataOutputType::nodes,
                Inpar::Solid::GaussPointDataOutputType::element_center,
                Inpar::Solid::GaussPointDataOutputType::gauss_points),
            &sublist_IO_VTK_structure);
      }


    }  // namespace Solid
  }    // namespace IORuntimeOutput
}  // namespace Inpar

FOUR_C_NAMESPACE_CLOSE
