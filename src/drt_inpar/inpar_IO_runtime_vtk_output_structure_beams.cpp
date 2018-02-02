/*----------------------------------------------------------------------*/
/*!
\file inpar_IO_runtime_vtk_output_structure_beams.cpp

\brief input parameters for VTK output of beam discretization at runtime

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/

#include "inpar_IO_runtime_vtk_output_structure_beams.H"

#include "drt_validparameters.H"
#include "inpar.H"
#include "inpar_parameterlist_utils.H"

#include <Teuchos_ParameterList.hpp>


namespace INPAR
{
namespace IO_RUNTIME_VTK
{
namespace BEAMS
{

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
  void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
  {
    using namespace DRT::INPUT;
    using Teuchos::tuple;
    using Teuchos::setStringToIntegralParameter;

    Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
    Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

    // related sublist
    Teuchos::ParameterList& sublist_IO = list->sublist("IO",false,"");
    Teuchos::ParameterList& sublist_IO_VTK_structure =
        sublist_IO.sublist("RUNTIME VTK OUTPUT",false,"");
    Teuchos::ParameterList& sublist_IO_VTK_beams =
        sublist_IO_VTK_structure.sublist("BEAMS",false,"");

    // whether to write special output for beam elements
    setStringToIntegralParameter<int>("OUTPUT_BEAMS","No",
                                 "write special output for beam elements",
                                 yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // whether to write displacement state
    setStringToIntegralParameter<int>("DISPLACEMENT","No",
                                 "write displacement output",
                                 yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // use absolute positions or initial positions for vtu geometry (i.e. point coordinates)
    // 'absolute positions' requires writing geometry in every output step (default for now)
    setStringToIntegralParameter<int>("USE_ABSOLUTE_POSITIONS","Yes",
        "use absolute positions or initial positions for vtu geometry (i.e. point coordinates)",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // write internal (elastic) energy of element
    setStringToIntegralParameter<int>("INTERNAL_ENERGY_ELEMENT","No",
        "write internal (elastic) energy for each element",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // write kinetic energy of element
    setStringToIntegralParameter<int>("KINETIC_ENERGY_ELEMENT","No",
        "write kinetic energy for each element",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // write triads as three orthonormal base vectors at every visualization point
    setStringToIntegralParameter<int>("TRIAD_VISUALIZATIONPOINT","No",
        "write triads at every visualization point",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // write material cross-section strains at the Gauss points:
    // axial & shear strains, twist & curvatures
    setStringToIntegralParameter<int>("STRAINS_GAUSSPOINT","No",
        "write material cross-section strains at the Gauss points",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // write material cross-section stresses at the Gauss points:
    // axial and shear forces, torque and bending moments
    setStringToIntegralParameter<int>("MATERIAL_FORCES_GAUSSPOINT","No",
        "write material cross-section stresses at the Gauss points",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // write spatial cross-section stresses at the Gauss points:
    // axial and shear forces, torque and bending moments
    setStringToIntegralParameter<int>("SPATIAL_FORCES_GAUSSPOINT","No",
        "write material cross-section stresses at the Gauss points",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // write element filament numbers and type
    setStringToIntegralParameter<int>("BEAMFILAMENTCONDITION","No",
        "write element filament numbers",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);

    // write element and network orientation parameter
    setStringToIntegralParameter<int>("ORIENTATION_PARAMETER","No",
        "write element filament numbers",
        yesnotuple, yesnovalue, &sublist_IO_VTK_beams);
  }

}
}
}

