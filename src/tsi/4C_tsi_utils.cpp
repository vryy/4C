/*----------------------------------------------------------------------*/
/*! \file

\brief utility functions for tsi problems

\level 2


*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#include "4C_tsi_utils.hpp"

#include "4C_coupling_volmortar_utils.hpp"
#include "4C_discretization_condition_periodic.hpp"
#include "4C_discretization_dofset.hpp"
#include "4C_discretization_dofset_predefineddofnumber.hpp"
#include "4C_discretization_fem_general_element_center.hpp"
#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_utils_createdis.hpp"
#include "4C_so3_plast_ssn.hpp"
#include "4C_so3_thermo.hpp"
#include "4C_thermo_ele_impl_utils.hpp"
#include "4C_thermo_element.hpp"

#include <Epetra_MpiComm.h>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | remove flag thermo from condition                         dano 12/11 |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string> TSI::UTILS::ThermoStructureCloneStrategy::ConditionsToCopy()
    const
{
  return {{"ThermoDirichlet", "Dirichlet"}, {"ThermoPointNeumann", "PointNeumann"},
      {"ThermoLineNeumann", "LineNeumann"}, {"ThermoSurfaceNeumann", "SurfaceNeumann"},
      {"ThermoVolumeNeumann", "VolumeNeumann"}, {"ThermoConvections", "ThermoConvections"},
      {"LinePeriodic", "LinePeriodic"}, {"SurfacePeriodic", "SurfacePeriodic"},
      {"ThermoInitfield", "Initfield"}, {"MortarMulti", "MortarMulti"}};
}


/*----------------------------------------------------------------------*
 | check material of cloned element                          dano 12/11 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::ThermoStructureCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  //  CORE::Materials::MaterialType mtype =
  //  GLOBAL::Problem::Instance()->Materials()->ById(matid)->Type(); if ((mtype !=
  //  CORE::Materials::m_th_fourier_iso))
  //    FOUR_C_THROW("Material with ID %d is not admissible for thermo elements",matid);

}  // CheckMaterialType()


/*----------------------------------------------------------------------*
 | set element data for cloned element                       dano 12/11 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::ThermoStructureCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // initialise kinematic type to geo_linear.
  // kintype is passed to the cloned thermo element
  INPAR::STR::KinemType kintype = INPAR::STR::KinemType::linear;
  // if oldele is a so3_base element or a so3_Plast element
  DRT::ELEMENTS::SoBase* so_base = dynamic_cast<DRT::ELEMENTS::SoBase*>(oldele);
  if (so_base != nullptr)
    kintype = so_base->KinematicType();
  else
    FOUR_C_THROW("oldele is neither a So_base element!");

  // note: SetMaterial() was reimplemented by the thermo element!

  Teuchos::RCP<DRT::ELEMENTS::Thermo> therm =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Thermo>(newele);
  if (therm != Teuchos::null)
  {
    // cloning to same material id -> use the same material instance
    if (so_base->Material()->Parameter()->Id() == matid)
      therm->SetMaterial(0, so_base->Material());
    else
      therm->SetMaterial(0, MAT::Factory(matid));
    therm->SetDisType(oldele->Shape());  // set distype as well!
    therm->SetKinematicType(kintype);    // set kintype in cloned thermal element
  }
  else
  {
    FOUR_C_THROW("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}  // SetElementData()


/*----------------------------------------------------------------------*
 | cloned element has to be a THERMO element                 dano 12/11 |
 *----------------------------------------------------------------------*/
bool TSI::UTILS::ThermoStructureCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // we only support thermo elements here
  eletype.push_back("THERMO");

  return true;  // yes, we copy EVERY element (no submeshes)
}  // DetermineEleType()


/*----------------------------------------------------------------------*
 | setup TSI                                                 dano 12/11 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::SetupTSI(const Epetra_Comm& comm)
{
  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis;
  structdis = GLOBAL::Problem::Instance()->GetDis("structure");
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs())
  {
    structdis->fill_complete();
    Epetra_Map nc = *(structdis->NodeColMap());
    Epetra_Map nr = *(structdis->NodeRowMap());
    structdis->Redistribute(nr, nc);
  }

  // access the thermo discretization
  Teuchos::RCP<DRT::Discretization> thermdis;
  thermdis = GLOBAL::Problem::Instance()->GetDis("thermo");
  if (!thermdis->Filled()) thermdis->fill_complete();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& tsidyn = GLOBAL::Problem::Instance()->TSIDynamicParams();

  bool matchinggrid = CORE::UTILS::IntegralValue<bool>(tsidyn, "MATCHINGGRID");

  // we use the structure discretization as layout for the temperature discretization
  if (structdis->NumGlobalNodes() == 0) FOUR_C_THROW("Structure discretization is empty!");

  // create thermo elements if the temperature discretization is empty
  if (thermdis->NumGlobalNodes() == 0)
  {
    if (!matchinggrid)
      FOUR_C_THROW(
          "MATCHINGGRID is set to 'no' in TSI DYNAMIC section, but thermo discretization is "
          "empty!");

    DRT::UTILS::CloneDiscretization<TSI::UTILS::ThermoStructureCloneStrategy>(structdis, thermdis);
    thermdis->fill_complete();

    // connect degrees of freedom for periodic boundary conditions
    {
      CORE::Conditions::PeriodicBoundaryConditions pbc_struct(structdis);

      if (pbc_struct.HasPBC())
      {
        pbc_struct.update_dofs_for_periodic_boundary_conditions();
      }
    }

    // connect degrees of freedom for periodic boundary conditions
    {
      CORE::Conditions::PeriodicBoundaryConditions pbc(thermdis);

      if (pbc.HasPBC())
      {
        pbc.update_dofs_for_periodic_boundary_conditions();
      }
    }

    // TSI must know the other discretization
    // build a proxy of the structure discretization for the temperature field
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> structdofset = structdis->GetDofSetProxy();
    // build a proxy of the temperature discretization for the structure field
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> thermodofset = thermdis->GetDofSetProxy();

    // check if ThermoField has 2 discretizations, so that coupling is possible
    if (thermdis->AddDofSet(structdofset) != 1) FOUR_C_THROW("unexpected dof sets in thermo field");
    if (structdis->AddDofSet(thermodofset) != 1)
      FOUR_C_THROW("unexpected dof sets in structure field");

    structdis->fill_complete(true, true, true);
    thermdis->fill_complete(true, true, true);

    TSI::UTILS::SetMaterialPointersMatchingGrid(structdis, thermdis);
  }
  else
  {
    if (matchinggrid)
      FOUR_C_THROW(
          "MATCHINGGRID is set to 'yes' in TSI DYNAMIC section, but thermo discretization is not "
          "empty!");

    // first call fill_complete for single discretizations.
    // This way the physical dofs are numbered successively
    structdis->fill_complete();
    thermdis->fill_complete();

    // build auxiliary dofsets, i.e. pseudo dofs on each discretization
    const int ndofpernode_thermo = 1;
    const int ndofperelement_thermo = 0;
    const int ndofpernode_struct = GLOBAL::Problem::Instance()->NDim();
    const int ndofperelement_struct = 0;
    Teuchos::RCP<CORE::Dofsets::DofSetInterface> dofsetaux;
    dofsetaux = Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(
        ndofpernode_thermo, ndofperelement_thermo, 0, true));
    if (structdis->AddDofSet(dofsetaux) != 1)
      FOUR_C_THROW("unexpected dof sets in structure field");
    dofsetaux = Teuchos::rcp(new CORE::Dofsets::DofSetPredefinedDoFNumber(
        ndofpernode_struct, ndofperelement_struct, 0, true));
    if (thermdis->AddDofSet(dofsetaux) != 1) FOUR_C_THROW("unexpected dof sets in thermo field");

    // call assign_degrees_of_freedom also for auxiliary dofsets
    // note: the order of fill_complete() calls determines the gid numbering!
    // 1. structure dofs
    // 2. thermo dofs
    // 3. structure auxiliary dofs
    // 4. thermo auxiliary dofs
    structdis->fill_complete(true, false, false);
    thermdis->fill_complete(true, false, false);
  }

}  // SetupTSI()


/*----------------------------------------------------------------------*
 | print TSI-logo                                            dano 03/10 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::SetMaterialPointersMatchingGrid(Teuchos::RCP<const DRT::Discretization> sourcedis,
    Teuchos::RCP<const DRT::Discretization> targetdis)
{
  const int numelements = targetdis->NumMyColElements();

  for (int i = 0; i < numelements; ++i)
  {
    DRT::Element* targetele = targetdis->lColElement(i);
    const int gid = targetele->Id();

    DRT::Element* sourceele = sourcedis->gElement(gid);

    // for coupling we add the source material to the target element and vice versa
    targetele->AddMaterial(sourceele->Material());
    sourceele->AddMaterial(targetele->Material());
  }
}

/*----------------------------------------------------------------------*
 |  assign material to discretization A                       vuong 09/14|
 *----------------------------------------------------------------------*/
void TSI::UTILS::TSIMaterialStrategy::AssignMaterial2To1(
    const CORE::VOLMORTAR::VolMortarCoupl* volmortar, DRT::Element* ele1,
    const std::vector<int>& ids_2, Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  // call default assignment
  CORE::VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial2To1(
      volmortar, ele1, ids_2, dis1, dis2);

  // done
  return;
};


/*----------------------------------------------------------------------*
|  assign material to discretization B                       vuong 09/14|
 *----------------------------------------------------------------------*/
void TSI::UTILS::TSIMaterialStrategy::AssignMaterial1To2(
    const CORE::VOLMORTAR::VolMortarCoupl* volmortar, DRT::Element* ele2,
    const std::vector<int>& ids_1, Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  // if no corresponding element found -> leave
  if (ids_1.empty()) return;

  // call default assignment
  CORE::VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial1To2(
      volmortar, ele2, ids_1, dis1, dis2);

  // initialise kinematic type to geo_linear.
  // kintype is passed to the corresponding thermo element
  INPAR::STR::KinemType kintype = INPAR::STR::KinemType::linear;

  // default strategy: take material of element with closest center in reference coordinates
  DRT::Element* ele1 = nullptr;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords2 = CORE::FE::element_center_refe_coords(*ele2);

    for (unsigned i = 0; i < ids_1.size(); ++i)
    {
      DRT::Element* actele1 = dis1->gElement(ids_1[i]);
      std::vector<double> centercoords1 = CORE::FE::element_center_refe_coords(*actele1);

      CORE::LINALG::Matrix<3, 1> diffcoords(true);

      for (int j = 0; j < 3; ++j) diffcoords(j, 0) = centercoords1[j] - centercoords2[j];

      if (diffcoords.Norm2() - mindistance < 1e-16)
      {
        mindistance = diffcoords.Norm2();
        ele1 = actele1;
      }
    }
  }

  // if Aele is a so3_base element
  DRT::ELEMENTS::SoBase* so_base = dynamic_cast<DRT::ELEMENTS::SoBase*>(ele1);
  if (so_base != nullptr)
    kintype = so_base->KinematicType();
  else
    FOUR_C_THROW("ele1 is not a so3_thermo element!");

  DRT::ELEMENTS::Thermo* therm = dynamic_cast<DRT::ELEMENTS::Thermo*>(ele2);
  if (therm != nullptr)
  {
    therm->SetKinematicType(kintype);  // set kintype in cloned thermal element
  }

  // done
  return;
}


/*----------------------------------------------------------------------*
 | print TSI-logo                                            dano 03/10 |
 *----------------------------------------------------------------------*/
void TSI::printlogo()
{
  // more at http://www.ascii-art.de under entry "rockets"
  std::cout << "Welcome to Thermo-Structure-Interaction " << std::endl;
  std::cout << "         !\n"
            << "         !\n"
            << "         ^\n"
            << "        / \\\n"
            << "       /___\\\n"
            << "      |=   =|\n"
            << "      |     |\n"
            << "      |     |\n"
            << "      |     |\n"
            << "      |     |\n"
            << "      |     |\n"
            << "      | TSI |\n"
            << "      |     |\n"
            << "     /|##!##|\\\n"
            << "    / |##!##| \\\n"
            << "   /  |##!##|  \\\n"
            << "  |  / ^ | ^ \\  |\n"
            << "  | /  ( | )  \\ |\n"
            << "  |/   ( | )   \\|\n"
            << "      ((   ))\n"
            << "     ((  :  ))\n"
            << "     ((  :  ))\n"
            << "      ((   ))\n"
            << "       (( ))\n"
            << "        ( )\n"
            << "         .\n"
            << "         .\n"
            << "         .\n"
            << "\n"
            << std::endl;
}  // printlogo()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
