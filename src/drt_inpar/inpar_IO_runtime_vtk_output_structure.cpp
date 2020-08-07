/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for VTK output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "inpar_IO_runtime_vtk_output_structure.H"

#include "drt_validparameters.H"
#include "inpar.H"
#include "inpar_parameterlist_utils.H"

#include <Teuchos_ParameterList.hpp>

namespace INPAR
{
  namespace IO_RUNTIME_VTK
  {
    namespace STRUCTURE
    {
      /*----------------------------------------------------------------------*
       *----------------------------------------------------------------------*/
      void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
      {
        using namespace DRT::INPUT;
        using Teuchos::setStringToIntegralParameter;
        using Teuchos::tuple;

        Teuchos::Array<std::string> yesnotuple =
            tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
        Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

        // related sublist
        Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
        Teuchos::ParameterList& sublist_IO_VTK =
            sublist_IO.sublist("RUNTIME VTK OUTPUT", false, "");
        Teuchos::ParameterList& sublist_IO_VTK_structure =
            sublist_IO_VTK.sublist("STRUCTURE", false, "");

        // whether to write output for structure
        setStringToIntegralParameter<int>("OUTPUT_STRUCTURE", "No", "write structure output",
            yesnotuple, yesnovalue, &sublist_IO_VTK_structure);

        // whether to write displacement state
        setStringToIntegralParameter<int>("DISPLACEMENT", "No", "write displacement output",
            yesnotuple, yesnovalue, &sublist_IO_VTK_structure);

        // whether to write element owner
        setStringToIntegralParameter<int>("ELEMENT_OWNER", "No", "write element owner", yesnotuple,
            yesnovalue, &sublist_IO_VTK_structure);

        // whether to write element GID
        setStringToIntegralParameter<int>("ELEMENT_GID", "No", "write baci internal element GIDs",
            yesnotuple, yesnovalue, &sublist_IO_VTK_structure);

        // whether to write stress and / or strain data
        setStringToIntegralParameter<int>("STRESS_STRAIN", "No",
            "Write element stress and / or strain  data. The type of stress / strain has to be "
            "selected in the --IO input section",
            yesnotuple, yesnovalue, &sublist_IO_VTK_structure);
      }


    }  // namespace STRUCTURE
  }    // namespace IO_RUNTIME_VTK
}  // namespace INPAR
