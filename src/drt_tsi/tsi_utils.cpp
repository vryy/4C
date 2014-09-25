/*----------------------------------------------------------------------*/
/*!
\file tsi_utils.cpp

\brief utility functions for tsi problems

<pre>
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>
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

#include "../drt_lib/drt_discret.H"

#include "../drt_fluid/drt_periodicbc.H"

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
  GenKinematicType kintype = geo_linear;
  // if oldele is a so3_base element or a so3_Plast element
  DRT::ELEMENTS::So3_Base*  so3_base  = dynamic_cast<DRT::ELEMENTS::So3_Base*>(oldele);
  DRT::ELEMENTS::So3_Plast<DRT::Element::hex8>* so3_plast_h8 =
      dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::Element::hex8> * >(oldele);
  DRT::ELEMENTS::So3_Plast<DRT::Element::hex27>* so3_plast_h27 =
      dynamic_cast<DRT::ELEMENTS::So3_Plast<DRT::Element::hex27> * >(oldele);
  if (so3_base != NULL)
    kintype = so3_base->GetKinematicType();
  else if (so3_plast_h8 != NULL || so3_plast_h27!=NULL)
    kintype = geo_nonlinear; // So3_Plast only with nonlinear kinematics
  else
    dserror("oldele is neither a so3_thermo nor a So3_Plast element!");

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
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
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
  Teuchos::RCP<DRT::Discretization> thermdis = Teuchos::null;
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
    Teuchos::RCP<DRT::DofSet> structdofset
      = structdis->GetDofSetProxy();
    // build a proxy of the temperature discretization for the structure field
    Teuchos::RCP<DRT::DofSet> thermodofset
      = thermdis->GetDofSetProxy();

    // check if ThermoField has 2 discretizations, so that coupling is possible
    if (thermdis->AddDofSet(structdofset)!=1)
      dserror("unexpected dof sets in thermo field");
    if (structdis->AddDofSet(thermodofset)!=1)
      dserror("unexpected dof sets in structure field");

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
    if (structdis->BuildDofSetAuxProxy(ndofpernode_thermo, ndofperelement_thermo, 0, true ) != 1)
      dserror("unexpected dof sets in structure field");
    if (thermdis->BuildDofSetAuxProxy(ndofpernode_struct, ndofperelement_struct, 0, true) != 1)
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
 *----------------------------------------------------------------------*/
void TSI::UTILS::TSIMaterialStrategy::AssignMaterialBToA(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* Aele,
    const std::vector<int>& Bids,
    Teuchos::RCP<DRT::Discretization> disA,
    Teuchos::RCP<DRT::Discretization> disB)
{
  //call default assignment
  VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterialBToA(volmortar,Aele,Bids,disA,disB);

  //done
  return;
};


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TSI::UTILS::TSIMaterialStrategy::AssignMaterialAToB(
    const VOLMORTAR::VolMortarCoupl* volmortar,
    DRT::Element* Bele,
    const std::vector<int>& Aids,
    Teuchos::RCP<DRT::Discretization> disA,
    Teuchos::RCP<DRT::Discretization> disB)
{
  //call default assignment
  VOLMORTAR::UTILS::DefaultMaterialStrategy::AssignMaterialAToB(volmortar,Bele,Aids,disA,disB);

  // initialise kinematic type to geo_linear.
  // kintype is passed to the corresponding thermo element
  GenKinematicType kintype = geo_linear;

  //default strategy: take only material of first element found
  DRT::Element* Aele = disA->gElement(Aids[0]);

  // if Aele is a so3_base element
  DRT::ELEMENTS::So3_Base* so3_base =
      dynamic_cast<DRT::ELEMENTS::So3_Base*>(Aele);
  if (so3_base != NULL)
    kintype = so3_base->GetKinematicType();
  else
    dserror("Aele is not a so3_thermo element!");

  DRT::ELEMENTS::Thermo* therm = dynamic_cast<DRT::ELEMENTS::Thermo*>(Bele);
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
