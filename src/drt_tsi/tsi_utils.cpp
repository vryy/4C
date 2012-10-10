/*----------------------------------------------------------------------*/
/*!
\file tsi_utils.cpp

\brief utility functions for tsi problems

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
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
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_thermo/thermo_element.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_condition_utils.H"
#include <Epetra_Time.h>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string,std::string> TSI::UTILS::ThermoStructureCloneStrategy::ConditionsToCopy()
{
  std::map<std::string,std::string> conditions_to_copy;

  // special Thermo conditions
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoDirichlet","Dirichlet"));
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoPointNeumann","PointNeumann"));
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoLineNeumann","LineNeumann"));
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoSurfaceNeumann","SurfaceNeumann"));
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoVolumeNeumann","VolumeNeumann"));

  // special Thermo convective heat transfer conditions (Newton's law of heat
  // transfer)
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoConvections","ThermoConvections"));

  return conditions_to_copy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TSI::UTILS::ThermoStructureCloneStrategy::CheckMaterialType(const int matid)
{
//  //// We take the material with the ID specified by the user
//  //// Here we check first, whether this material is of admissible type
//  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
//  if ((mtype != INPAR::MAT::m_th_fourier_iso))
//  dserror("Material with ID %d is not admissible for thermo elements",matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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

  // note: SetMaterial() was reimplemented by the thermo element!
#if defined(D_THERMO)
      DRT::ELEMENTS::Thermo* therm = dynamic_cast<DRT::ELEMENTS::Thermo*>(newele.get());
      if (therm!=NULL)
      {
        therm->SetMaterial(matid);
        therm->SetDisType(oldele->Shape()); // set distype as well!
      }
      else
#endif
    {
      dserror("unsupported element type '%s'", typeid(*newele).name());
    }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool TSI::UTILS::ThermoStructureCloneStrategy::DetermineEleType(
  DRT::Element* actele,
  const bool ismyele,
  vector<string>& eletype
  )
{
  // we only support thermo elements here
  eletype.push_back("THERMO");

  return true; // yes, we copy EVERY element (no submeshes)
}


/*----------------------------------------------------------------------*
 | setup TSI                                                 dano 12/11 |
 *----------------------------------------------------------------------*/
void TSI::UTILS::SetupTSI(const Epetra_Comm& comm)
{
  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs()) structdis->FillComplete();

  // access the thermo discretization
  Teuchos::RCP<DRT::Discretization> thermdis = Teuchos::null;
  thermdis = DRT::Problem::Instance()->GetDis("thermo");
  if (!thermdis->Filled()) thermdis->FillComplete();

   // we use the structure discretization as layout for the temperature discretization
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");

  // create thermo elements if the temperature discretization is empty
  if (thermdis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<TSI::UTILS::ThermoStructureCloneStrategy>(structdis,thermdis);
  }
  else
      dserror("Structure AND Thermo discretization present. This is not supported.");
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
}


/*----------------------------------------------------------------------*/
