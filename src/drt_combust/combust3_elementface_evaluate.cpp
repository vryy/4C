/*!
\file combust3_elementface_evaluate.cpp
\brief

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/

#include <Teuchos_TimeMonitor.hpp>

#include "combust3.H"
#include "combust3_sysmat.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/enrichment_utils.H"
#include "../drt_inpar/inpar_xfem.H"
#include "combust_flamefront.H"
#include "combust_interface.H"


/*----------------------------------------------------------------------*
 | converts a string into an Action for this element    rasthofer 02/13 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Combust3IntFace::ActionType DRT::ELEMENTS::Combust3IntFace::convertStringToActionType(
              const std::string& action) const
{
  DRT::ELEMENTS::Combust3IntFace::ActionType act = Combust3IntFace::none;
  if (action == "calc_edge_based_stab_terms")
    act = Combust3IntFace::calc_edge_based_stab_terms;
  else
    dserror("Unknown type of action for Combust3IntFace");
  return act;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                       rasthofer 02/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3IntFace::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  dserror("Don't use it.");
  return 1;
}


/*----------------------------------------------------------------------*
|  evaluate the element (public)                       rasthofer 02/13 |
*----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3IntFace::Evaluate(
    Teuchos::ParameterList&            params,               ///< parameter list
    DRT::Discretization&               discretization,       ///< discretization
    std::vector<int>&                  lm_master,            ///< master local map
    std::vector<int>&                  lm_slave,             ///< slave local map
    std::map<XFEM::PHYSICS::Field,std::vector<int> >&   lm_masterDofPerFieldToPatch, ///< local map between master nodes and nodes in patch
    std::map<XFEM::PHYSICS::Field,std::vector<int> >&   lm_slaveDofPerFieldToPatch,  ///< local map between slave nodes and nodes in patch
    std::map<std::pair<XFEM::PHYSICS::Field,XFEM::PHYSICS::Field>,Epetra_SerialDenseMatrix>&  elemat_blocks,   ///< element matrix blocks
    std::map<XFEM::PHYSICS::Field,Epetra_SerialDenseVector>&  elevec_blocks)    ///< element vector blocks
{
  // get the action required
  const std::string action(params.get<std::string>("action","none"));
  const DRT::ELEMENTS::Combust3IntFace::ActionType act = convertStringToActionType(action);

  // get the list of materials
  // note: it is assumed that both master and slave have the same material
  // we can not ask the face here!
  const Teuchos::RCP<MAT::Material> material = ParentMasterElement()->Material();

  switch (act)
  {
    case calc_edge_based_stab_terms:
    {
      // get combust type
      const INPAR::COMBUST::CombustionType combusttype   = DRT::INPUT::get<INPAR::COMBUST::CombustionType>(params, "combusttype");
      const INPAR::COMBUST::SelectedEnrichment selected_enr = DRT::INPUT::get<INPAR::COMBUST::SelectedEnrichment>(params, "selectedenrichment");
      if (combusttype == INPAR::COMBUST::combusttype_twophaseflow or (combusttype == INPAR::COMBUST::combusttype_twophaseflow_surf
          and selected_enr != INPAR::COMBUST::selectedenrichment_pressure))
        dserror("Face-based terms not supported for this problem");

      // get interface handler
      const COMBUST::InterfaceHandleCombust* ih = (params.get<Teuchos::RCP<COMBUST::InterfaceHandleCombust> >("interface handle")).get();

      // get time integration specific parameters
      const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");

      const double time   = params.get<double>("time");
      const double dt     = params.get<double>("dt");
      const double theta  = params.get<double>("theta");
      const double ga_gamma  = params.get<double>("gamma");
      const double ga_alphaF = params.get<double>("alphaF");
      const double ga_alphaM = params.get<double>("alphaM");

      // generalized alpha time integration scheme
      bool genalpha = false;
      if (timealgo == INPAR::FLUID::timeint_afgenalpha) genalpha = true;

      // get stabilization specific parameter
      const INPAR::FLUID::EOS_Pres pres_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Pres>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_PRES");
      const INPAR::FLUID::EOS_Conv_Stream conv_stream_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Conv_Stream>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_CONV_STREAM");
      const INPAR::FLUID::EOS_Conv_Cross conv_cross_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Conv_Cross>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_CONV_CROSS");
      const INPAR::FLUID::EOS_Div conti_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Div>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_DIV");
      if (conv_cross_stab != INPAR::FLUID::EOS_CONV_CROSS_none)
          dserror("EOS_CONV_CROSS not yet supported!");
      if (pres_stab == INPAR::FLUID::EOS_PRES_xfem_gp or conv_stream_stab == INPAR::FLUID::EOS_CONV_STREAM_xfem_gp
          or conti_stab == INPAR::FLUID::EOS_DIV_div_jump_xfem_gp or conti_stab == INPAR::FLUID::EOS_DIV_vel_jump_xfem_gp)
          dserror("EOS_xfem_gp not yet supported by combustion module!");
      // flag to indicate inclusion of ghost penalties
      const bool xfem_stab = params.get<bool>("xfemstab");

      const INPAR::FLUID::EOS_ElementLength hk_def = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_ElementLength>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_H_DEFINITION");
      const INPAR::FLUID::EOS_TauType tau_def = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_TauType>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_DEFINITION_TAU");

      const INPAR::FLUID::EOS_GP_Pattern pattern =params.sublist("EDGE-BASED STABILIZATION").get<INPAR::FLUID::EOS_GP_Pattern>("EOS_PATTERN");

      // get face type: defines whether face has to be integrated twice
      const INPAR::XFEM::FaceType facetype = params.get<INPAR::XFEM::FaceType>("facetype");

      // extract nodal valus from state vectors
      std::vector<double> m_velnp;     ///< velocity/pressure unknowns at n+1 (master)
      std::vector<double> m_velaf;     ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (master)
      std::vector<double> m_phinp;     ///< phi at n+1 (master)
      // get phinp
      Teuchos::RCP<Epetra_Vector> phinp = params.get<Teuchos::RCP<Epetra_Vector> >("phinp",Teuchos::null);
      if (phinp==Teuchos::null) dserror("Could not get phinp!");
      if (!genalpha)
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),m_velnp,lm_master);
      else
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velaf") ,m_velaf ,lm_master);

      DRT::UTILS::ExtractMyNodeBasedValues(ParentMasterElement(), m_phinp,*phinp);

      std::vector<double> s_velnp;     ///< velocity/pressure unknowns at n+1 (slave)
      std::vector<double> s_velaf;     ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (slave)
      std::vector<double> s_phinp;     ///< phi at n+1 (slave)
      if (!genalpha)
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),s_velnp,lm_slave);
      else
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velaf") ,s_velaf ,lm_slave);

      DRT::UTILS::ExtractMyNodeBasedValues(ParentSlaveElement(), s_phinp,*phinp);

      // get assembly type for master and slave element
      XFEM::AssemblyType m_assembly_type=XFEM::standard_assembly;
      // suppress enrichment dofs
      if (params.get<int>("selectedenrichment") != INPAR::COMBUST::selectedenrichment_none)
        m_assembly_type = XFEM::ComputeAssemblyType(*(ParentMasterElement()->GetEleDofManager()),
                                                    ParentMasterElement()->NumNode(),
                                                    ParentMasterElement()->NodeIds());
      XFEM::AssemblyType s_assembly_type=XFEM::standard_assembly;
      // suppress enrichment dofs
      if (params.get<int>("selectedenrichment") != INPAR::COMBUST::selectedenrichment_none)
        s_assembly_type = XFEM::ComputeAssemblyType(*(ParentSlaveElement()->GetEleDofManager()),
                                                      ParentSlaveElement()->NumNode(),
                                                      ParentSlaveElement()->NodeIds());

      // call calculation of face mat
      COMBUST::callFacemat(this, ih, elemat_blocks, elevec_blocks, combusttype,
                           m_velnp, m_velaf, m_phinp,
                           s_velnp, s_velaf, s_phinp,
                           material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma,
                           pres_stab, conv_stream_stab, conv_cross_stab, conti_stab, hk_def, tau_def, pattern, xfem_stab, facetype,
                           m_assembly_type, s_assembly_type, lm_masterDofPerFieldToPatch, lm_slaveDofPerFieldToPatch,
                           *(ParentMasterElement()->GetEleDofManager()), *(ParentSlaveElement()->GetEleDofManager()));
    }
    break;
    default:
      dserror("Unknown action type for combust face element!");
      break;
  }

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition rasthofer 02/13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Combust3IntFace::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  dserror("not available");

  return 0;
}
