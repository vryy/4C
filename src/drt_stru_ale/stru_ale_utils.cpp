/*----------------------------------------------------------------------*/
/*!
\file stru_ale_utils.cpp

\brief utility functions for structure with ale problems

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*
 | definitions                                               mgit 04/11 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

/*----------------------------------------------------------------------*
 | headers                                                   mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "stru_ale_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_ale2/ale2.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string,std::string> STRU_ALE::UTILS::AleStructureCloneStrategy::ConditionsToCopy()
{
  std::map<std::string,std::string> conditions_to_copy;

  // special Thermo conditions
  conditions_to_copy.insert(pair<std::string,std::string>("ALEDirichlet","Dirichlet"));
  conditions_to_copy.insert(pair<std::string,std::string>("AleWear","AleWear"));
  
  return conditions_to_copy;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STRU_ALE::UTILS::AleStructureCloneStrategy::CheckMaterialType(const int matid)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STRU_ALE::UTILS::AleStructureCloneStrategy::SetElementData(
  Teuchos::RCP<DRT::Element> newele,
  DRT::Element* oldele,
  const int matid,
  const bool isnurbs
  )
{
  // We must not add a new material type here because that might move
  // the internal material vector. And each element material might
  // have a pointer to that vector. Too bad.
  // So we search for a Fourier material and take the first one we find.
  // => matid from outside remains unused!
  const int matnr = DRT::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_stvenant);
  if (matnr==-1)
    dserror("No StVenantKirchhoff material defined. Cannot generate ale mesh.");

  // note: SetMaterial() was reimplemented by the thermo element!
#ifdef D_ALE
      DRT::ELEMENTS::Ale2* ale2 = dynamic_cast<DRT::ELEMENTS::Ale2*>(newele.get());
      //DRT::ELEMENTS::Ale2* ale2 = dynamic_cast<DRT::ELEMENTS::Ale2*>(newele.get());
      if (ale2!=NULL)
      {
        ale2->SetMaterial(matnr);
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
bool STRU_ALE::UTILS::AleStructureCloneStrategy::DetermineEleType(
  DRT::Element* actele,
  const bool ismyele,
  vector<string>& eletype
  )
{
  // we only support ale elements here
  eletype.push_back("ALE2");
  
  return true; // yes, we copy EVERY element (no submeshes)
}

/*----------------------------------------------------------------------*/
#endif // CCADISCRET
