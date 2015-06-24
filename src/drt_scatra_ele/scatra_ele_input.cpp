/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_input.cpp
\brief

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_linedefinition.H"

#include "../drt_mat/myocard.H"

#include "scatra_ele.H"


/*----------------------------------------------------------------------*
 | read element input                                        fang 02/15 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Transport::ReadElement(
    const std::string& eletype,
    const std::string& distype,
    DRT::INPUT::LineDefinition* linedef
    )
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  // read implementation type
  std::string impltype;
  linedef->ExtractString("TYPE",impltype);
  if(impltype == "Std")
    impltype_ = INPAR::SCATRA::impltype_std;
  else if(impltype == "AdvReac")
    impltype_ = INPAR::SCATRA::impltype_advreac;
  else if(impltype == "Aniso")
    impltype_ = INPAR::SCATRA::impltype_aniso;
  else if(impltype == "CardMono")
    impltype_ = INPAR::SCATRA::impltype_cardiac_monodomain;
  else if(impltype == "ElchDiffCond")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond;
  else if(impltype == "ElchDiffCondThermo")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond_thermo;
  else if(impltype == "ElchElectrode")
    impltype_ = INPAR::SCATRA::impltype_elch_electrode;
  else if(impltype == "ElchElectrodeThermo")
    impltype_ = INPAR::SCATRA::impltype_elch_electrode_thermo;
  else if(impltype == "ElchNP")
    impltype_ = INPAR::SCATRA::impltype_elch_NP;
  else if(impltype == "Loma")
    impltype_ = INPAR::SCATRA::impltype_loma;
  else if(impltype == "Ls")
    impltype_ = INPAR::SCATRA::impltype_levelset;
  else if(impltype == "LsReinit")
    impltype_ = INPAR::SCATRA::impltype_lsreinit;
  else if(impltype == "Poro")
    impltype_ = INPAR::SCATRA::impltype_poro;
  else if(impltype == "PoroReac")
    impltype_ = INPAR::SCATRA::impltype_pororeac;
  else
    dserror("Transport element received invalid implementation type!");

  // set discretization type
  SetDisType(DRT::StringToDistype(distype));

  if (Material()->MaterialType() == INPAR::MAT::m_myocard)
  {
    Teuchos::RCP<MAT::Myocard> myocard = Teuchos::rcp_dynamic_cast<MAT::Myocard>(Material());
    myocard->Setup(linedef);
  }

  return true;
}
