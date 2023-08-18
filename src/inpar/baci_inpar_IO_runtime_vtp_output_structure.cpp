/*----------------------------------------------------------------------*/
/*! \file
\brief input parameters for VTP output of structural problem at runtime

\level 2

*/
/*----------------------------------------------------------------------*/

#include "baci_inpar_IO_runtime_vtp_output_structure.H"

#include "baci_inpar.H"
#include "baci_inpar_parameterlist_utils.H"
#include "baci_inpar_validparameters.H"

#include <Teuchos_ParameterList.hpp>

namespace INPAR
{
  namespace IO_RUNTIME_VTP_STRUCTURE
  {
    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
    {
      using namespace DRT::INPUT;
      using Teuchos::setStringToIntegralParameter;
      using Teuchos::tuple;

      // related sublist
      Teuchos::ParameterList& sublist_IO = list->sublist("IO", false, "");
      Teuchos::ParameterList& sublist_IO_VTP_structure =
          sublist_IO.sublist("RUNTIME VTP OUTPUT STRUCTURE", false, "");


      // output interval regarding steps: write output every INTERVAL_STEPS steps
      IntParameter("INTERVAL_STEPS", -1, "write VTP output at runtime every INTERVAL_STEPS steps",
          &sublist_IO_VTP_structure);

      IntParameter("STEP_OFFSET", 0,
          "An offset added to the current step to shift the steps to be written.",
          &sublist_IO_VTP_structure);

      // whether to write output in every iteration of the nonlinear solver
      BoolParameter("EVERY_ITERATION", "No",
          "write output in every iteration of the nonlinear solver", &sublist_IO_VTP_structure);

      // write owner at every visualization point
      BoolParameter("OWNER", "No", "write owner of every point", &sublist_IO_VTP_structure);

      // write orientation at every visualization point
      BoolParameter("ORIENTATIONANDLENGTH", "No", "write orientation at every point",
          &sublist_IO_VTP_structure);

      // write number of bonds at every visualization point
      BoolParameter(
          "NUMBEROFBONDS", "No", "write number of bonds of every point", &sublist_IO_VTP_structure);

      // write force actin in linker
      BoolParameter(
          "LINKINGFORCE", "No", "write force acting in linker", &sublist_IO_VTP_structure);
    }


  }  // namespace IO_RUNTIME_VTP_STRUCTURE
}  // namespace INPAR
