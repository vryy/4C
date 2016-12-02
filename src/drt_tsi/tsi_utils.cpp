/*----------------------------------------------------------------------*/
/*!
\file tsi_utils.cpp

\brief utility functions for tsi problems

\level 2

\maintainer Alexander Seitz

*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif


/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include "tsi_utils.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_thermo/thermo_element.H"
#include "../drt_so3/so3_thermo.H"
#include "../drt_so3/so3_ssn_plast.H"

#include "../drt_lib/drt_dofset.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_periodicbc.H"

#include "../drt_volmortar/volmortar_utils.H"

/*----------------------------------------------------------------------*
 | remove flag thermo from condition                         dano 12/11 |
 *----------------------------------------------------------------------*/
std::map<std::string,std::string> TSI::UTILS::ThermoStructureCloneStrategy::ConditionsToCopy()
{
  std::map<std::string,std::string> conditions_to_copy;

  // special Thermo conditions
  conditions_to_copy.insert(std::pair<std::string,std::string>("ThermoDirichlet","Dirichlet"));
  conditions_to_copy.insert(std::pair<std::string,std::string>("ThermoPointNeumann","PointNeumann"));
  conditions_to_copy.insert(std::pair<std::string,std::string>("ThermoLineNeumann","LineNeumann"));
  conditions_to_copy.insert(std::pair<std::string,std::string>("ThermoSurfaceNeumann","SurfaceNeumann"));
  conditions_to_copy.insert(std::pair<std::string,std::string>("ThermoVolumeNeumann","VolumeNeumann"));

  // special Thermo convective heat transfer conditions (Newton's law of heat
  // transfer)
  conditions_to_copy.insert(std::pair<std::string,std::string>("ThermoConvections","ThermoConvections"));

  // conditions for periodic boundary conditions
  conditions_to_copy.insert(std::pair<std::string,std::string>("LinePeriodic","LinePeriodic"));
  conditions_to_copy.insert(std::pair<std::string,std::string>("SurfacePeriodic","SurfacePeriodic"));

  /// initial field
  conditions_to_copy.insert(std::pair<std::string,std::string>("ThermoInitfield","Initfield"));

  return conditions_to_copy;
}  // ConditionsToCopy()


/*----------------------------------------------------------------------*
 | check material of cloned element                          dano 12/11 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::ThermoStructureCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
//  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
//  if ((mtype != INPAR::MAT::m_th_fourier_iso))
//    dserror("Material with ID %d is not admissible for thermo elements",matid);

}  // CheckMaterialType()


/*----------------------------------------------------------------------*
 | set element data for cloned element                       dano 12/11 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::ThermoStructureCloneStrategy::SetElementData(
  Teuchos::RCP<DRT::Element> newele,
  DRT::Element* oldele,
  const int matid,
  const bool isnurbs
  )
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // initialise kinematic type to geo_linear.
  // kintype is passed to the cloned thermo element
  INPAR::STR::KinemType kintype = INPAR::STR::kinem_linear;
  // if oldele is a so3_base element or a so3_Plast element
  DRT::ELEMENTS::So_base*  so_base  = dynamic_cast<DRT::ELEMENTS::So_base*>(oldele);
  if (so_base != NULL)
    kintype = so_base->KinematicType();
  else
    dserror("oldele is neither a So_base element!");

  // note: SetMaterial() was reimplemented by the thermo element!
#if defined(D_THERMO)
  Teuchos::RCP<DRT::ELEMENTS::Thermo> therm = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Thermo>(newele);
  if (therm != Teuchos::null)
  {
    therm->SetMaterial(matid);
    therm->SetDisType(oldele->Shape());  // set distype as well!
    therm->SetKinematicType(kintype);  // set kintype in cloned thermal element
  }
  else
#endif
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}  // SetElementData()


/*----------------------------------------------------------------------*
 | cloned element has to be a THERMO element                 dano 12/11 |
 *----------------------------------------------------------------------*/
bool TSI::UTILS::ThermoStructureCloneStrategy::DetermineEleType(
  DRT::Element* actele,
  const bool ismyele,
  std::vector<std::string>& eletype
  )
{
  // we only support thermo elements here
  eletype.push_back("THERMO");

  return true; // yes, we copy EVERY element (no submeshes)
}  // DetermineEleType()


/*----------------------------------------------------------------------*
 | setup TSI                                                 dano 12/11 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::SetupTSI(const Epetra_Comm& comm)
{
  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs())
  {
    structdis->FillComplete();
    Epetra_Map nc=*(structdis->NodeColMap());
    Epetra_Map nr=*(structdis->NodeRowMap());
    structdis->Redistribute(nr,nc);
  }

  // access the thermo discretization
  Teuchos::RCP<DRT::Discretization> thermdis;
  thermdis = DRT::Problem::Instance()->GetDis("thermo");
  if (!thermdis->Filled()) thermdis->FillComplete();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& tsidyn
    = DRT::Problem::Instance()->TSIDynamicParams();

  bool matchinggrid = DRT::INPUT::IntegralValue<bool>(tsidyn,"MATCHINGGRID");

   // we use the structure discretization as layout for the temperature discretization
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");

  // create thermo elements if the temperature discretization is empty
  if (thermdis->NumGlobalNodes()==0)
  {
    if(!matchinggrid)
      dserror("MATCHINGGRID is set to 'no' in TSI DYNAMIC section, but thermo discretization is empty!");

    DRT::UTILS::CloneDiscretization<TSI::UTILS::ThermoStructureCloneStrategy>(structdis,thermdis);
    thermdis->FillComplete();

    // connect degrees of freedom for periodic boundary conditions
    {
      PeriodicBoundaryConditions pbc_struct(structdis);

      if (pbc_struct.HasPBC())
      {
        pbc_struct.UpdateDofsForPeriodicBoundaryConditions();
      }
    }

    // connect degrees of freedom for periodic boundary conditions
    {
      PeriodicBoundaryConditions pbc(thermdis);

      if (pbc.HasPBC())
      {
        pbc.UpdateDofsForPeriodicBoundaryConditions();
      }
    }

    // TSI must know the other discretization
    // build a proxy of the structure discretization for the temperature field
    Teuchos::RCP<DRT::DofSetInterface> structdofset = structdis->GetDofSetProxy();
    // build a proxy of the temperature discretization for the structure field
    Teuchos::RCP<DRT::DofSetInterface> thermodofset = thermdis->GetDofSetProxy();

    // check if ThermoField has 2 discretizations, so that coupling is possible
    if (thermdis->AddDofSet(structdofset)!=1)
      dserror("unexpected dof sets in thermo field");
    if (structdis->AddDofSet(thermodofset)!=1)
      dserror("unexpected dof sets in structure field");

    structdis->FillComplete(true,true,true);
    thermdis ->FillComplete(true,true,true);

    TSI::UTILS::SetMaterialPointersMatchingGrid(structdis,thermdis);
  }
  else
  {
    if(matchinggrid)
      dserror("MATCHINGGRID is set to 'yes' in TSI DYNAMIC section, but thermo discretization is not empty!");

    //first call FillComplete for single discretizations.
    //This way the physical dofs are numbered successively
    structdis->FillComplete();
    thermdis->FillComplete();

    //build auxiliary dofsets, i.e. pseudo dofs on each discretization
    const int ndofpernode_thermo = 1;
    const int ndofperelement_thermo  = 0;
    const int ndofpernode_struct = DRT::Problem::Instance()->NDim();
    const int ndofperelement_struct = 0;
    Teuchos::RCP<DRT::DofSetInterface> dofsetaux;
    dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(ndofpernode_thermo, ndofperelement_thermo, 0, true ));
    if (structdis->AddDofSet(dofsetaux) != 1)
      dserror("unexpected dof sets in structure field");
    dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(ndofpernode_struct, ndofperelement_struct, 0, true));
    if (thermdis->AddDofSet(dofsetaux) != 1)
      dserror("unexpected dof sets in thermo field");

    //call AssignDegreesOfFreedom also for auxiliary dofsets
    //note: the order of FillComplete() calls determines the gid numbering!
    // 1. structure dofs
    // 2. thermo dofs
    // 3. structure auxiliary dofs
    // 4. thermo auxiliary dofs
    structdis->FillComplete(true, false,false);
    thermdis->FillComplete(true, false,false);
  }

}  // SetupTSI()


/*----------------------------------------------------------------------*
 | print TSI-logo                                            dano 03/10 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::SetMaterialPointersMatchingGrid(
    Teuchos::RCP<const DRT::Discretization> sourcedis,
    Teuchos::RCP<const DRT::Discretization> targetdis)
{
  const int numelements = targetdis->NumMyColElements();

  for (int i=0; i<numelements; ++i)
  {
    DRT::Element* targetele = targetdis->lColElement(i);
    const int gid = targetele->Id();

    DRT::Element* sourceele = sourcedis->gElement(gid);

    //for coupling we add the source material to the target element and vice versa
    targetele->AddMaterial(sourceele->Material());
    sourceele->AddMaterial(targetele->Material());
  }
}

/*----------------------------------------------------------------------*
 |  assign material to discretization A                       vuong 09/14|
 *----------------------------------------------------------------------*/
void TSI::UTILS::TSIMaterialStrategy::AssignMaterial2To1(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* ele1,
    const std::vector<int>& ids_2,
    Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  //call default assignment
  VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial2To1(volmortar,ele1,ids_2,dis1,dis2);

  //done
  return;
};


/*----------------------------------------------------------------------*
|  assign material to discretization B                       vuong 09/14|
 *----------------------------------------------------------------------*/
void TSI::UTILS::TSIMaterialStrategy::AssignMaterial1To2(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* ele2,
    const std::vector<int>& ids_1,
    Teuchos::RCP<DRT::Discretization> dis1,
    Teuchos::RCP<DRT::Discretization> dis2)
{
  //if no corresponding element found -> leave
  if(ids_1.empty())
    return;

  //call default assignment
  VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterial1To2(volmortar,ele2,ids_1,dis1,dis2);

  // initialise kinematic type to geo_linear.
  // kintype is passed to the corresponding thermo element
  INPAR::STR::KinemType kintype = INPAR::STR::kinem_linear;

  //default strategy: take material of element with closest center in reference coordinates
  DRT::Element* ele1 = NULL;
  double mindistance = 1e10;
  {
    std::vector<double> centercoords2 = DRT::UTILS::ElementCenterRefeCoords(ele2);

    for (unsigned i=0; i<ids_1.size(); ++i)
    {
      DRT::Element* actele1= dis1->gElement(ids_1[i]);
      std::vector<double> centercoords1 = DRT::UTILS::ElementCenterRefeCoords(actele1);

      LINALG::Matrix<3,1> diffcoords(true);

      for (int j=0; j<3; ++j)
        diffcoords(j,0)=centercoords1[j]-centercoords2[j];

      if(diffcoords.Norm2()-mindistance<1e-16)
      {
        mindistance=diffcoords.Norm2();
        ele1 = actele1;
      }
    }
  }

  // if Aele is a so3_base element
  DRT::ELEMENTS::So_base* so_base =
      dynamic_cast<DRT::ELEMENTS::So_base*>(ele1);
  if (so_base != NULL)
    kintype = so_base->KinematicType();
  else
    dserror("ele1 is not a so3_thermo element!");

  DRT::ELEMENTS::Thermo* therm = dynamic_cast<DRT::ELEMENTS::Thermo*>(ele2);
  if (therm != NULL)
  {
    therm->SetKinematicType(kintype); // set kintype in cloned thermal element
  }

  //done
  return;
}


/*----------------------------------------------------------------------*
 | print TSI-logo                                            dano 03/10 |
 *----------------------------------------------------------------------*/
void TSI::printlogo()
{
  // more at http://www.ascii-art.de under entry "rockets"
  std::cout << "Welcome to Thermo-Structure-Interaction " << std::endl;
  std::cout<<"         !\n"
  <<"         !\n"
  <<"         ^\n"
  <<"        / \\\n"
  <<"       /___\\\n"
  <<"      |=   =|\n"
  <<"      |     |\n"
  <<"      |     |\n"
  <<"      |     |\n"
  <<"      |     |\n"
  <<"      |     |\n"
  <<"      | TSI |\n"
  <<"      |     |\n"
  <<"     /|##!##|\\\n"
  <<"    / |##!##| \\\n"
  <<"   /  |##!##|  \\\n"
  <<"  |  / ^ | ^ \\  |\n"
  <<"  | /  ( | )  \\ |\n"
  <<"  |/   ( | )   \\|\n"
  <<"      ((   ))\n"
  <<"     ((  :  ))\n"
  <<"     ((  :  ))\n"
  <<"      ((   ))\n"
  <<"       (( ))\n"
  <<"        ( )\n"
  <<"         .\n"
  <<"         .\n"
  <<"         .\n"
  <<"\n"<<std::endl;
}  // printlogo()


/*----------------------------------------------------------------------*/
