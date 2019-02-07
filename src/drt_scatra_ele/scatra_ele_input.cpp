/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_input.cpp
\brief Read scalar transport element

<pre>
\level 1

\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15237
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
    const std::string& eletype, const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  // read implementation type
  std::string impltype;
  linedef->ExtractString("TYPE", impltype);
  if (impltype == "Std")
    impltype_ = INPAR::SCATRA::impltype_std;
  else if (impltype == "AdvReac")
    impltype_ = INPAR::SCATRA::impltype_advreac;
  else if (impltype == "RefConcReac")
    impltype_ = INPAR::SCATRA::impltype_refconcreac;
  else if (impltype == "Chemo")
    impltype_ = INPAR::SCATRA::impltype_chemo;
  else if (impltype == "ChemoReac")
    impltype_ = INPAR::SCATRA::impltype_chemoreac;
  else if (impltype == "Aniso")
    impltype_ = INPAR::SCATRA::impltype_aniso;
  else if (impltype == "CardMono")
    impltype_ = INPAR::SCATRA::impltype_cardiac_monodomain;
  else if (impltype == "ElchDiffCond")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond;
  else if (impltype == "ElchDiffCondMultiScale")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond_multiscale;
  else if (impltype == "ElchDiffCondThermo")
    impltype_ = INPAR::SCATRA::impltype_elch_diffcond_thermo;
  else if (impltype == "ElchElectrode")
    impltype_ = INPAR::SCATRA::impltype_elch_electrode;
  else if (impltype == "ElchElectrodeGrowth")
    impltype_ = INPAR::SCATRA::impltype_elch_electrode_growth;
  else if (impltype == "ElchElectrodeThermo")
    impltype_ = INPAR::SCATRA::impltype_elch_electrode_thermo;
  else if (impltype == "ElchNP")
    impltype_ = INPAR::SCATRA::impltype_elch_NP;
  else if (impltype == "Loma")
    impltype_ = INPAR::SCATRA::impltype_loma;
  else if (impltype == "Ls")
    impltype_ = INPAR::SCATRA::impltype_levelset;
  else if (impltype == "LsReinit")
    impltype_ = INPAR::SCATRA::impltype_lsreinit;
  else if (impltype == "Poro")
    impltype_ = INPAR::SCATRA::impltype_poro;
  else if (impltype == "PoroReac")
    impltype_ = INPAR::SCATRA::impltype_pororeac;
  else if (impltype == "PoroReacECM")
    impltype_ = INPAR::SCATRA::impltype_pororeacECM;
  else if (impltype == "Hdg")
    impltype_ = INPAR::SCATRA::impltype_std_hdg;
  else if (impltype == "HdgCardMono")
    impltype_ = INPAR::SCATRA::impltype_cardiac_monodomain_hdg;
  else if (impltype == "BondReac")
    impltype_ = INPAR::SCATRA::impltype_bondreac;
  else if (impltype == "VariationalDiffusion")
    impltype_ = INPAR::SCATRA::impltype_variational_diffusion;
  else
    dserror("Transport element received invalid implementation type!");

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(material);

  // set discretization type
  SetDisType(DRT::StringToDistype(distype));

  if (Material()->MaterialType() == INPAR::MAT::m_myocard)
  {
    Teuchos::RCP<MAT::Myocard> myocard = Teuchos::rcp_dynamic_cast<MAT::Myocard>(Material());
    myocard->Setup(linedef);
  }

  return true;
}
