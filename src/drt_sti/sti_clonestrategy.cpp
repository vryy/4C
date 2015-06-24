/*----------------------------------------------------------------------*/
/*!
\file sti_clonestrategy.cpp

\brief strategy for cloning thermo discretization from scatra discretization

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_inpar/inpar_material.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"

#include "../drt_scatra_ele/scatra_ele.H"

#include "sti_clonestrategy.H"

/*----------------------------------------------------------------------*
 | check material of cloned element                          fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::ScatraThermoCloneStrategy::CheckMaterialType(
    const int matid   //! material of cloned element
    )
{
  // check whether material with specified ID is compatible with cloned element or not
  if(DRT::Problem::Instance()->Materials()->ById(matid)->Type() != INPAR::MAT::m_scatra)
    dserror("Material with ID %d is not compatible with cloned transport element!",matid);

  return;
}


/*--------------------------------------------------------------------------*
 | return map with original names of conditions to be cloned as key values, |
 | and final names of cloned conditions as mapped values         fang 04/15 |
 *--------------------------------------------------------------------------*/
std::map<std::string,std::string> STI::ScatraThermoCloneStrategy::ConditionsToCopy()
{
  // initialize map
  std::map<std::string,std::string> conditions;

  // insert thermo conditions
  conditions.insert(std::pair<std::string,std::string>("ThermoDirichlet","Dirichlet"));
  conditions.insert(std::pair<std::string,std::string>("ThermoPointNeumann","PointNeumann"));
  conditions.insert(std::pair<std::string,std::string>("ThermoLineNeumann","LineNeumann"));
  conditions.insert(std::pair<std::string,std::string>("ThermoSurfaceNeumann","SurfaceNeumann"));
  conditions.insert(std::pair<std::string,std::string>("ThermoVolumeNeumann","VolumeNeumann"));
  conditions.insert(std::pair<std::string,std::string>("ThermoInitfield","Initfield"));

  // return map
  return conditions;
} // STI::ScatraThermoCloneStrategy::ConditionsToCopy()


/*----------------------------------------------------------------------------------------------------------*
 | decide whether element should be cloned or not, and if so, determine type of cloned element   fang 04/15 |
 *----------------------------------------------------------------------------------------------------------*/
bool STI::ScatraThermoCloneStrategy::DetermineEleType(
    DRT::Element*               actele,    //! current element on source discretization
    const bool                  ismyele,   //! ownership flag
    std::vector<std::string>&   eletype    //! vector storing types of cloned elements
    )
{
  // set type of cloned element to transport type
  eletype.push_back("TRANSP");

  // element should always be cloned
  return true;
} // STI::ScatraThermoCloneStrategy::DetermineEleType


/*----------------------------------------------------------------------*
 | provide cloned element with element specific data         fang 04/15 |
 *----------------------------------------------------------------------*/
void STI::ScatraThermoCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element>   newele,   //! current cloned element on target discretization
    DRT::Element*                oldele,   //! current element on source discretization
    const int                    matid,    //! material of cloned element
    const bool                   isnurbs   //! nurbs flag
    )
{
  // cast pointer to current cloned element on target discretization
  Teuchos::RCP<DRT::ELEMENTS::Transport> newele_transport = Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Transport>(newele);

  // safety check
  if(newele_transport == Teuchos::null)
    dserror("Expected transport element, but received element of type '%s'!",typeid(*newele).name());

  // provide cloned element with material
  newele_transport->SetMaterial(matid,oldele);

  // provide cloned element with discretization type
  newele_transport->SetDisType(oldele->Shape());

  // provide cloned element with physical implementation type
  newele_transport->SetImplType(INPAR::SCATRA::impltype_thermo);

  return;
}
