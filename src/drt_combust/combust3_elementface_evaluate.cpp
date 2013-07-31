/*!
\file combust3_elementface_evaluate.cpp
\brief

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/

#include <Teuchos_TimeMonitor.hpp>

#include "combust3.H"
#include "combust3_sysmat.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/enrichment_utils.H"


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
    std::vector<int>&                  lm_patch,             ///< patch local map
    std::vector<int>&                  lm_face,              ///< face local map
    std::vector<int>&                  lm_master,            ///< master local map
    std::vector<int>&                  lm_slave,             ///< slave local map
    std::vector<int>&                  lm_masterToPatch,     ///< local map between master dofs and patchlm
    std::vector<int>&                  lm_slaveToPatch,      ///< local map between slave dofs and patchlm
    std::vector<int>&                  lm_faceToPatch,       ///< local map between face dofs and patchlm
    std::vector<int>&                  lm_masterNodeToPatch, ///< local map between master nodes and nodes in patch
    std::vector<int>&                  lm_slaveNodeToPatch,  ///< local map between slave nodes and nodes in patch
    std::vector<Epetra_SerialDenseMatrix>&  elemat_blocks,   ///< element matrix blocks
    std::vector<Epetra_SerialDenseVector>&  elevec_blocks)    ///< element vector blocks
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
        std::cout << "Combust3IntFace Evaluate()" << std::endl;
      // im Standardfall werden an dieser Stelle
      // - alle parmeter ausgelesen
      // - die element state-Groessen geholt
      // - eine assembly strategie definiert
      // - call sysmat gerufen

      // get combust type
      const INPAR::COMBUST::CombustionType combusttype   = DRT::INPUT::get<INPAR::COMBUST::CombustionType>(params, "combusttype");

      // get time integration specific parameters
      const INPAR::FLUID::TimeIntegrationScheme timealgo = DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(params, "timealgo");

      const double time   = params.get<double>("time");
      const double dt     = params.get<double>("dt");
      const double theta  = params.get<double>("theta");
      const double ga_gamma  = params.get<double>("gamma");
      const double ga_alphaF = params.get<double>("alphaF");
      const double ga_alphaM = params.get<double>("alphaM");

      // instationary formulation
      bool instationary = true;
      if (timealgo == INPAR::FLUID::timeint_stationary) instationary = false;
      // generalized alpha time integration scheme
      bool genalpha = false;
      if (timealgo == INPAR::FLUID::timeint_afgenalpha) genalpha = true;

      // get stabilization specific parameter
      const INPAR::FLUID::EOS_Pres pres_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Pres>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_PRES");
      const INPAR::FLUID::EOS_Conv_Stream conv_stream_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Conv_Stream>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_CONV_STREAM");
      const INPAR::FLUID::EOS_Conv_Cross conv_cross_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Conv_Cross>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_CONV_CROSS");
      const INPAR::FLUID::EOS_Div conti_stab = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_Div>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_DIV");
      if (conti_stab == INPAR::FLUID::EOS_DIV_vel_jump_std_eos or conti_stab == INPAR::FLUID::EOS_DIV_vel_jump_xfem_gp)
          dserror("EOS_DIV_vel_jump not yet supported");
      if (conv_cross_stab != INPAR::FLUID::EOS_CONV_CROSS_none)
          dserror("EOS_CONV_CROSS not yet supported");
      // TODO: was ist fuer xfem interessant: brauche ich conti stab ueberhaupt

      const INPAR::FLUID::EOS_ElementLength hk_def = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_ElementLength>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_H_DEFINITION");
      const INPAR::FLUID::EOS_TauType tau_def = DRT::INPUT::IntegralValue<INPAR::FLUID::EOS_TauType>(params.sublist("EDGE-BASED STABILIZATION"),"EOS_DEFINITION_TAU");

      const INPAR::FLUID::EOS_GP_Pattern pattern =params.sublist("EDGE-BASED STABILIZATION").get<INPAR::FLUID::EOS_GP_Pattern>("EOS_PATTERN");
//      // extract local (element level) vectors from global state vectors
//      DRT::ELEMENTS::Combust3::MyState mystate(
//          discretization,
//          lm,
//          instationary,
//          genalpha,
//          gradphi,
//          this,
//          epetra_phinp_,
//          gradphi_,
//          curvature_);

      // TODO: welche vectoren brauche ich wirklich
      std::vector<double> m_velnp;     ///< velocity/pressure unknowns at n+1 (master)
      std::vector<double> m_veln;      ///< velocity/pressure unknowns at n (master)
      std::vector<double> m_velnm;     ///< velocity/pressure unknowns at n-1 (master)
      std::vector<double> m_velaf;     ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (master)
      std::vector<double> m_accn;      ///< time derivative of unknowns at n (master)
      std::vector<double> m_accam;     ///< time derivative of unknowns at alpha+m for generalized alpha scheme (master)
      if (!genalpha)
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),m_velnp,lm_master);
      if (instationary)
      {
        DRT::UTILS::ExtractMyValues(*discretization.GetState("veln") ,m_veln ,lm_master);
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velnm"),m_velnm,lm_master);
        DRT::UTILS::ExtractMyValues(*discretization.GetState("accn") ,m_accn ,lm_master);
        if(genalpha)
        {
          DRT::UTILS::ExtractMyValues(*discretization.GetState("velaf") ,m_velaf ,lm_master);
          DRT::UTILS::ExtractMyValues(*discretization.GetState("accam") ,m_accam ,lm_master);
        }
      }

      std::vector<double> s_velnp;     ///< velocity/pressure unknowns at n+1 (slave)
      std::vector<double> s_veln;      ///< velocity/pressure unknowns at n (slave)
      std::vector<double> s_velnm;     ///< velocity/pressure unknowns at n-1 (slave)
      std::vector<double> s_velaf;     ///< velocity/pressure unknowns at alpha+f for generalized alpha scheme (slave)
      std::vector<double> s_accn;      ///< time derivative of unknowns at n (slave)
      std::vector<double> s_accam;     ///< time derivative of unknowns at alpha+m for generalized alpha scheme (slave)
      if (!genalpha)
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velnp"),s_velnp,lm_slave);
      if (instationary)
      {
        DRT::UTILS::ExtractMyValues(*discretization.GetState("veln") ,s_veln ,lm_slave);
        DRT::UTILS::ExtractMyValues(*discretization.GetState("velnm"),s_velnm,lm_slave);
        DRT::UTILS::ExtractMyValues(*discretization.GetState("accn") ,s_accn ,lm_slave);
        if(genalpha)
        {
          DRT::UTILS::ExtractMyValues(*discretization.GetState("velaf") ,s_velaf ,lm_slave);
          DRT::UTILS::ExtractMyValues(*discretization.GetState("accam") ,s_accam ,lm_slave);
        }
      }
//      std::vector<double> phinp_;     ///< G-function unknowns at n+1
//      std::vector<double> gradphinp_; ///< smoothed G-function gradient at n+1
//      std::vector<double> curv_;      ///< curvature based on G-function gradient at n+1

      // get assmebly type for master and slave element
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
      COMBUST::callFacemat(this, elemat_blocks, elevec_blocks, combusttype,
                           m_velnp, m_veln, m_velnm, m_velaf, m_accn, m_accam,
                           s_velnp, s_veln, s_velnm, s_velaf, s_accn, s_accam,
                           material, timealgo, time, dt, theta, ga_alphaF, ga_alphaM, ga_gamma,
                           pres_stab, conv_stream_stab, conv_cross_stab, conti_stab, hk_def, tau_def, pattern,
                           m_assembly_type, s_assembly_type, lm_masterNodeToPatch, lm_slaveNodeToPatch,
                           *(ParentMasterElement()->GetEleDofManager()), *(ParentSlaveElement()->GetEleDofManager()));

    }
    break;
    default:
      dserror("Unknown action type for combust face element!");
      break;
  }

  //dserror("not available");

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
