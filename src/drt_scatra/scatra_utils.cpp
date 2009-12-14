/*----------------------------------------------------------------------*/
/*!
\file scatra_utils.cpp

\brief utility functions for scalar transport problems

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "scatra_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_scatra/scatra_element.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<string,string> SCATRA::ScatraFluidCloneStrategy::ConditionsToCopy()
{
  std::map<string,string> conditions_to_copy;

  conditions_to_copy.insert(pair<string,string>("TransportDirichlet","Dirichlet"));
  conditions_to_copy.insert(pair<string,string>("TransportPointNeumann","PointNeumann"));
  conditions_to_copy.insert(pair<string,string>("TransportLineNeumann","LineNeumann"));
  conditions_to_copy.insert(pair<string,string>("TransportSurfaceNeumann","SurfaceNeumann"));
  conditions_to_copy.insert(pair<string,string>("TransportVolumeNeumann","VolumeNeumann"));
  conditions_to_copy.insert(pair<string,string>("TransportNeumannInflow","TransportNeumannInflow"));
  // when the fluid problem is periodic we also expect the mass transport to be so:
  conditions_to_copy.insert(pair<string,string>("LinePeriodic","LinePeriodic"));
  conditions_to_copy.insert(pair<string,string>("SurfacePeriodic","SurfacePeriodic"));

  conditions_to_copy.insert(pair<string,string>("LineNeumann","FluidLineNeumann"));
  conditions_to_copy.insert(pair<string,string>("SurfaceNeumann","FluidSurfaceNeumann"));
  conditions_to_copy.insert(pair<string,string>("VolumeNeumann","FluidVolumeNeumann"));
  conditions_to_copy.insert(pair<string,string>("KrylovSpaceProjection","KrylovSpaceProjection"));
  conditions_to_copy.insert(pair<string,string>("ElectrodeKinetics","ElectrodeKinetics"));

  // a hack:
  conditions_to_copy.insert(pair<string,string>("FluidStressCalc","FluxCalculation"));

  return conditions_to_copy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScatraFluidCloneStrategy::CheckMaterialType(const int matid)
{
// We take the material with the ID specified by the user
// Here we check first, whether this material is of admissible type
INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
if ((mtype != INPAR::MAT::m_scatra) &&
    (mtype != INPAR::MAT::m_mixfrac) &&
    (mtype != INPAR::MAT::m_sutherland) &&
    (mtype != INPAR::MAT::m_arrhenius_pv) &&
    (mtype != INPAR::MAT::m_matlist))
  dserror("Material with ID %d is not admissible for scalar transport elements",matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScatraFluidCloneStrategy::SetElementData(
    RCP<DRT::Element> newele,
    DRT::Element* oldele,
    const int matid,
    const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
#if defined(D_FLUID2) || defined(D_FLUID3)
      DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
      if (trans!=NULL)
      {
        trans->SetMaterial(matid);
        trans->SetDisType(oldele->Shape()); // set distype as well!
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
bool SCATRA::ScatraFluidCloneStrategy::DetermineEleType(
    DRT::Element* actele,
    const bool ismyele,
    vector<string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.push_back("TRANSP");

  return true; // yes, we copy EVERY element (no submeshes)
}




#endif  // CCADISCRET
