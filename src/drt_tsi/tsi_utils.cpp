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
#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
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
#include "../drt_thermo/thermo_element.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string,std::string> TSI::UTILS::ThermoStructureCloneStrategy::ConditionsToCopy()
{
  std::map<std::string,std::string> conditions_to_copy;

  //the special Thermo conditions
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoDirichlet","Dirichlet"));
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoPointNeumann","PointNeumann"));
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoLineNeumann","LineNeumann"));
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoSurfaceNeumann","SurfaceNeumann"));
  conditions_to_copy.insert(pair<std::string,std::string>("ThermoVolumeNeumann","VolumeNeumann"));

  return conditions_to_copy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TSI::UTILS::ThermoStructureCloneStrategy::CheckMaterialType(
  const int matid
  )
{
  //// We take the material with the ID specified by the user
  //// Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if (
       (mtype != INPAR::MAT::m_th_fourier_iso)&&
       (mtype != INPAR::MAT::m_matlist)
      )
  dserror("Material with ID %d is not admissible for thermo elements",matid);
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
//  // We must not add a new material type here because that might move
//  // the internal material vector. And each element material might
//  // have a pointer to that vector. Too bad.
//  // So we search for a Fourier material and take the first one we find.
//  // => matid from outside remains unused!
//
//  const int matnr =
//    DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_th_fourier_iso);
//  if (matnr==-1)
//    dserror("No isotropic Fourier material defined. Cannot generate thermo mesh.");
//
//  // We need to set material and possibly other things to complete element setup.
//  // This is again really ugly as we have to extract the actual
//  // element type in order to access the material property
//
//  // note: SetMaterial() was reimplemented by the thermo element!
//#if defined(D_THERMO)
//      DRT::ELEMENTS::Thermo* therm = dynamic_cast<DRT::ELEMENTS::Thermo*>(newele.get());
//      if (therm!=NULL)
//      {
//        therm->SetMaterial(matnr);
//        therm->SetDisType(oldele->Shape()); // set distype as well!
//      }
//      else
//#endif
//    {
//      dserror("unsupported element type '%s'", typeid(*newele).name());
//    }
//  return;


  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // like in ScatraFluidCloneStrategy
  // note: SetMaterial() was reimplemented by the thermo element!
#if defined(D_THERMO)
      DRT::ELEMENTS::Thermo* therm = dynamic_cast<DRT::ELEMENTS::Thermo*>(newele.get());
      if (therm!=NULL)
      {
        // now set the material
        therm->SetMaterial(matid);

     // TODO:
        // we want to have a copy of the structure material
        // --> in case of coupling: copy and set the structure material to the
        // thermo element

        // get Material of oldele (solid material)

        // get the 2nd material of the matlist (solid material in thermo field)
        // (= structmaterial but without correct values for parameters)


        // set Material of oldele on place of thrstvenantk in newele (thermo)

        // get the first thermostvenant material in thermo element

        // copy at this place the structure material
    // End of TODO

        // Set distype to thermo element (same shape like structure)
        therm->SetDisType(oldele->Shape());
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
#endif // CCADISCRET
