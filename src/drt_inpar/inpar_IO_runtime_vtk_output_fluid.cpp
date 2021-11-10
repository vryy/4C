/*----------------------------------------------------------------------*/
/*! \file

\brief input parameters for VTK output of a fluid field at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "inpar_IO_runtime_vtk_output_fluid.H"

#include "drt_validparameters.H"
#include "inpar.H"
#include "inpar_parameterlist_utils.H"
#include "inpar_fluid.H"

#include <Teuchos_ParameterList.hpp>

namespace INPAR
{
  namespace IO_RUNTIME_VTK
  {
    namespace FLUID
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
        Teuchos::ParameterList& sublist_IO_VTK_fluid = sublist_IO_VTK.sublist("FLUID", false, "");

        // whether to write output for fluid
        setStringToIntegralParameter<int>("OUTPUT_FLUID", "No", "write fluid output", yesnotuple,
            yesnovalue, &sublist_IO_VTK_fluid);

        // whether to write velocity state
        setStringToIntegralParameter<int>("VELOCITY", "No", "write velocity output", yesnotuple,
            yesnovalue, &sublist_IO_VTK_fluid);

        // whether to write pressure state
        setStringToIntegralParameter<int>("PRESSURE", "No", "write pressure output", yesnotuple,
            yesnovalue, &sublist_IO_VTK_fluid);

        // whether to write acceleration state
        setStringToIntegralParameter<int>("ACCELERATION", "No", "write acceleration output",
            yesnotuple, yesnovalue, &sublist_IO_VTK_fluid);

        // whether to write displacement state
        setStringToIntegralParameter<int>("DISPLACEMENT", "No", "write displacement output",
            yesnotuple, yesnovalue, &sublist_IO_VTK_fluid);

        // whether to write displacement state
        setStringToIntegralParameter<int>("GRIDVELOCITY", "No", "write grid velocity output",
            yesnotuple, yesnovalue, &sublist_IO_VTK_fluid);

        // whether to write element owner
        setStringToIntegralParameter<int>("ELEMENT_OWNER", "No", "write element owner", yesnotuple,
            yesnovalue, &sublist_IO_VTK_fluid);

        // whether to write element GIDs
        setStringToIntegralParameter<int>("ELEMENT_GID", "No", "write baci internal element GIDs",
            yesnotuple, yesnovalue, &sublist_IO_VTK_fluid);

        // whether to write node GIDs
        setStringToIntegralParameter<int>("NODE_GID", "No", "write baci internal node GIDs",
            yesnotuple, yesnovalue, &sublist_IO_VTK_fluid);
      }
    }  // namespace FLUID
  }    // namespace IO_RUNTIME_VTK
}  // namespace INPAR
