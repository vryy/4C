/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "baci_inpar_IO_runtime_vtk_output_structure.hpp"

#include "baci_inpar.hpp"
#include "baci_inpar_structure.hpp"
#include "baci_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace INPAR
{
  namespace IO_RUNTIME_OUTPUT
  {
    namespace STRUCTURE
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
        Teuchos::ParameterList& sublist_IO_VTK =
            sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");
        Teuchos::ParameterList& sublist_IO_VTK_structure =
            sublist_IO_VTK.sublist("STRUCTURE", false, "");

        // whether to write output for structure
        CORE::UTILS::BoolParameter(
            "OUTPUT_STRUCTURE", "No", "write structure output", &sublist_IO_VTK_structure);

        // whether to write displacement state
        CORE::UTILS::BoolParameter(
            "DISPLACEMENT", "No", "write displacement output", &sublist_IO_VTK_structure);

        // whether to write velocity state
        CORE::UTILS::BoolParameter(
            "VELOCITY", "No", "write velocity output", &sublist_IO_VTK_structure);

        // whether to write element owner
        CORE::UTILS::BoolParameter(
            "ELEMENT_OWNER", "No", "write element owner", &sublist_IO_VTK_structure);

        // whether to write element GIDs
        CORE::UTILS::BoolParameter(
            "ELEMENT_GID", "No", "write baci internal element GIDs", &sublist_IO_VTK_structure);

        // write element ghosting information
        CORE::UTILS::BoolParameter("ELEMENT_GHOSTING", "No",
            "write which processors ghost the elements", &sublist_IO_VTK_structure);

        // whether to write node GIDs
        CORE::UTILS::BoolParameter(
            "NODE_GID", "No", "write baci internal node GIDs", &sublist_IO_VTK_structure);

        // whether to write stress and / or strain data
        CORE::UTILS::BoolParameter("STRESS_STRAIN", "No",
            "Write element stress and / or strain  data. The type of stress / strain has to be "
            "selected in the --IO input section",
            &sublist_IO_VTK_structure);

        // mode to write gauss point data
        setStringToIntegralParameter<INPAR::STR::GaussPointDataOutputType>(
            "GAUSS_POINT_DATA_OUTPUT_TYPE", "none",
            "Where to write gauss point data. (none, projected to nodes, projected to element "
            "center, raw at gauss points)",
            tuple<std::string>("none", "nodes", "element_center", "gauss_points"),
            tuple<INPAR::STR::GaussPointDataOutputType>(INPAR::STR::GaussPointDataOutputType::none,
                INPAR::STR::GaussPointDataOutputType::nodes,
                INPAR::STR::GaussPointDataOutputType::element_center,
                INPAR::STR::GaussPointDataOutputType::gauss_points),
            &sublist_IO_VTK_structure);
      }


    }  // namespace STRUCTURE
  }    // namespace IO_RUNTIME_OUTPUT
}  // namespace INPAR

FOUR_C_NAMESPACE_CLOSE
