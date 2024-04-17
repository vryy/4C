/*----------------------------------------------------------------------*/
/*! \file

\brief Communicates between xfluid and NIT contact --> for XFSCI and XFPSCI

\level 3

*/
/*----------------------------------------------------------------------*/

#include "baci_xfem_xfluid_contact_communicator.hpp"

#include "baci_contact_element.hpp"
#include "baci_contact_nitsche_strategy_fpi.hpp"
#include "baci_contact_nitsche_strategy_fsi.hpp"
#include "baci_cut_boundingbox.hpp"
#include "baci_cut_cutwizard.hpp"
#include "baci_cut_element.hpp"
#include "baci_cut_elementhandle.hpp"
#include "baci_cut_facet.hpp"
#include "baci_cut_output.hpp"
#include "baci_cut_position.hpp"
#include "baci_cut_sidehandle.hpp"
#include "baci_cut_volumecell.hpp"
#include "baci_discretization_fem_general_extract_values.hpp"
#include "baci_fluid_ele.hpp"
#include "baci_mat_newtonianfluid.hpp"
#include "baci_mortar_element.hpp"
#include "baci_so3_hex8.hpp"
#include "baci_so3_surface.hpp"
#include "baci_xfem_condition_manager.hpp"

FOUR_C_NAMESPACE_OPEN

void XFEM::XFluid_Contact_Comm::InitializeFluidState(Teuchos::RCP<CORE::GEO::CutWizard> cutwizard,
    Teuchos::RCP<DRT::Discretization> fluiddis,
    Teuchos::RCP<XFEM::ConditionManager> condition_manager,
    Teuchos::RCP<Teuchos::ParameterList> fluidparams)
{
  if (fluiddis != Teuchos::null && condition_manager != Teuchos::null) fluid_init_ = true;
  cutwizard_ = cutwizard;
  fluiddis_ = fluiddis;
  condition_manager_ = condition_manager;

  parallel_ = (fluiddis_->Comm().NumProc() > 1);

  Teuchos::ParameterList& params_xf_stab = fluidparams->sublist("XFLUID DYNAMIC/STABILIZATION");

  visc_stab_trace_estimate_ = CORE::UTILS::IntegralValue<INPAR::XFEM::ViscStab_TraceEstimate>(
      params_xf_stab, "VISC_STAB_TRACE_ESTIMATE");
  visc_stab_hk_ =
      CORE::UTILS::IntegralValue<INPAR::XFEM::ViscStab_hk>(params_xf_stab, "VISC_STAB_HK");
  nit_stab_gamma_ = params_xf_stab.get<double>("NIT_STAB_FAC");
  is_pseudo_2D_ = (bool)CORE::UTILS::IntegralValue<int>(params_xf_stab, "IS_PSEUDO_2D");
  mass_conservation_scaling_ = CORE::UTILS::IntegralValue<INPAR::XFEM::MassConservationScaling>(
      params_xf_stab, "MASS_CONSERVATION_SCALING");
  mass_conservation_combination_ =
      CORE::UTILS::IntegralValue<INPAR::XFEM::MassConservationCombination>(
          params_xf_stab, "MASS_CONSERVATION_COMBO");


  INPAR::XFEM::ConvStabScaling ConvStabScaling =
      CORE::UTILS::IntegralValue<INPAR::XFEM::ConvStabScaling>(params_xf_stab, "CONV_STAB_SCALING");
  if (ConvStabScaling != INPAR::XFEM::ConvStabScaling_none)
    dserror("ConvStabScaling not handled correctly!");

  extrapolate_to_zero_ = (bool)CORE::UTILS::IntegralValue<int>(
      GLOBAL::Problem::Instance()->XFluidDynamicParams().sublist("XFPSI MONOLITHIC"),
      "EXTRAPOLATE_TO_ZERO");

  if (extrapolate_to_zero_)
    std::cout << "==| The Fluid Stress Extrapolation is relaxed to zero! |==" << std::endl;

  dt_ = fluidparams->get<double>("time step size");
  theta_ = fluidparams->get<double>("theta");

  mc_.clear();
  Teuchos::RCP<XFEM::MeshCouplingFSI> mcfsi = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSI>(
      condition_manager->GetMeshCoupling("XFEMSurfFSIMono"));
  if (mcfsi != Teuchos::null) mc_.push_back(mcfsi);

  Teuchos::RCP<XFEM::MeshCouplingFPI> mcfpi = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(
      condition_manager->GetMeshCoupling("XFEMSurfFPIMono_ps_ps"));
  if (mcfpi != Teuchos::null)
  {
    mc_.push_back(mcfpi);
    mcfpi_ps_pf_ = Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(
        condition_manager->GetMeshCoupling("XFEMSurfFPIMono_ps_pf"));
    if (mcfpi_ps_pf_ == Teuchos::null) dserror("Couldn't find mcfpi_ps_pf_ object!");
  }
  else
  {
    mcfpi_ps_pf_ = Teuchos::null;
  }

  if (mc_.size())
  {
    if (!fluiddis_->Comm().MyPID())
      std::cout << "==| XFluid_Contact_Comm: Loaded " << mc_.size()
                << " Mesh Coupling Objects! |==" << std::endl;
  }
  else
    dserror("Didn't find any mesh coupling object!");

  for (std::size_t mc = 0; mc < mc_.size(); ++mc)
  {
    if (mc_[mc]->GetAveragingStrategy() != INPAR::XFEM::Xfluid_Sided &&
        mass_conservation_scaling_ != INPAR::XFEM::MassConservationScaling_only_visc)
      dserror("The implementation does not what you expect!");
    else if (mc_[mc]->GetAveragingStrategy() != INPAR::XFEM::Xfluid_Sided &&
             visc_stab_trace_estimate_ != INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue)
      dserror("The implementation does not what you expect!");
  }

  my_sele_ids_.clear();

  PrepareIterationStep();

  last_physical_sides_ =
      std::pair<int, std::set<CORE::GEO::CUT::Side*>>(-2, std::set<CORE::GEO::CUT::Side*>());

  last_ele_h_ = std::pair<int, double>(-1, -1);

  Create_New_Gmsh_files();
}

double XFEM::XFluid_Contact_Comm::Get_FSI_Traction(MORTAR::Element* ele,
    const CORE::LINALG::Matrix<3, 1>& xsi_parent, const CORE::LINALG::Matrix<2, 1>& xsi_boundary,
    const CORE::LINALG::Matrix<3, 1>& normal, bool& FSI_integrated,
    bool& gp_on_this_proc,  // for serial run
    double* poropressure)
{
  gp_on_this_proc = true;
  if (!fluid_init_) dserror("Fluid not initialized!");
  if (ele == nullptr) dserror("Contact Element not set!");

  if (parallel_)
  {
    if (contact_ele_rowmap_fluidownerbased_->LID(ele->Id()) == -1)
    {
      gp_on_this_proc = false;
      return 0.0;
    }
  }

  DRT::ELEMENTS::StructuralSurface* sele = GetSurfEle(ele->Id());
  mcidx_ = GetSurfMc(ele->Id());

  if (Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(mc_[mcidx_]) != Teuchos::null)
    isporo_ = true;
  else
    isporo_ = false;

  // 1 // get CutSide--> BC --> VolumeCell --> FluidElement --> Check for Dofs
  CORE::GEO::CUT::SideHandle* sidehandle = cutwizard_->GetCutSide(GetSurfSid(ele->Id()));
  if (sidehandle == nullptr) dserror("Coundn't find Sidehandle for this structural surface!");

  std::vector<int> nds;
  int eleid;

  CORE::GEO::CUT::VolumeCell* volumecell = nullptr;
  static CORE::LINALG::Matrix<3, 1> elenormal(true);
  static CORE::LINALG::Matrix<3, 1> x(false);
  CORE::LINALG::Matrix<2, 1> new_xsi(xsi_boundary.A(), false);
  double distance = 0.0;
  if (!GetVolumecell(sele, new_xsi, sidehandle, nds, eleid, volumecell, elenormal, x,
          FSI_integrated, distance))
  {
    gp_on_this_proc = false;
    return 0.0;
  }

  if (!volumecell)
  {
    std::cout << "==| You have a mail from XFluid_Contact_Comm: As I couldn't find an appropriate "
                 "fluid solution at your gausspoint I decided that the solution is 0.0! |=="
              << std::endl;
    return 0.0;
  }

  static std::vector<double> velpres;
  static std::vector<double> disp;
  static std::vector<double> ivel;
  DRT::Element* fluidele = nullptr;
  CORE::LINALG::SerialDenseMatrix ele_xyze;
  double pres_m;
  static CORE::LINALG::Matrix<3, 1> vel_m;
  static CORE::LINALG::Matrix<3, 1> vel_s;
  static CORE::LINALG::Matrix<3, 1> velpf_s;
  static CORE::LINALG::Matrix<3, 3> vderxy_m;
  Get_States(eleid, nds, sele, new_xsi, x, fluidele, ele_xyze, velpres, disp, ivel, pres_m, vel_m,
      vel_s, vderxy_m, velpf_s);

  double penalty_fac = 0.0;

  // To get the actual element normal we take the normal from the Contact, as this is the actual
  // one!
  elenormal.Update(-1, normal, 0.0);

  if (mc_[mcidx_]->GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided)
    Get_Penalty_Param(fluidele, volumecell, ele_xyze, elenormal, penalty_fac, vel_m);
  else if (mc_[mcidx_]->GetAveragingStrategy() == INPAR::XFEM::Embedded_Sided)
    Get_Penalty_Param(sele, penalty_fac);
  else
    dserror(
        "Your interface stress averaging strategy is not yet implemented for XFSCI and XFPSCI!");

  double visc_m;
  mc_[mcidx_]->GetViscosityMaster(fluidele, visc_m);

  static double porosity = -1;
  if (isporo_)
  {
    CORE::LINALG::Matrix<3, 1> xsi3(new_xsi.A(), true);
    double J = 0;
    porosity =
        Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFPI>(mc_[mcidx_])->CalcPorosity(sele, xsi3, J);
  }

  if (mc_[mcidx_]->GetAveragingStrategy() == INPAR::XFEM::Xfluid_Sided)
  {
    if (poropressure && distance > 1e-10)
    {
      double scaling = (distance) / (mc_[mcidx_]->Get_h());
      if (scaling < 1.0)
      {
        return -(*poropressure) * (scaling) +
               (1 - scaling) * XFEM::UTILS::Evaluate_Full_Traction(pres_m, vderxy_m, visc_m,
                                   penalty_fac, vel_m, vel_s, elenormal, elenormal, velpf_s,
                                   porosity);
      }
      else
      {
        return -(*poropressure);
      }
    }
    else
    {
      if (!extrapolate_to_zero_)
        return XFEM::UTILS::Evaluate_Full_Traction(pres_m, vderxy_m, visc_m, penalty_fac, vel_m,
            vel_s, elenormal, elenormal, velpf_s, porosity);
      else
      {
        double scaling = (distance) / (mc_[mcidx_]->Get_h());
        if (scaling > 1.0) scaling = 1;
        return XFEM::UTILS::Evaluate_Full_Traction(pres_m, vderxy_m, visc_m, penalty_fac, vel_m,
                   vel_s, elenormal, elenormal, velpf_s, porosity) *
               (1 - scaling);
      }
    }
  }
  else
  {
    dserror("Other Nitsche FSI variant than fluid-sided weighting not implemented yet!");
    return 0.0;
  }
}

bool XFEM::XFluid_Contact_Comm::CheckNitscheContactState(CONTACT::Element* cele,
    const CORE::LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
    const double& full_fsi_traction,        // stressfluid + penalty
    double& gap                             // gap
)
{
  if (contact_strategy_fsi_)
    return contact_strategy_fsi_->CheckNitscheContactState(cele, xsi, full_fsi_traction, gap);
  else if (contact_strategy_fpi_)
    return contact_strategy_fpi_->CheckNitscheContactState(cele, xsi, full_fsi_traction, gap);
  else
    dserror("CheckNitscheContactState: Not adequate contact strategy is assigned!");
  return false;  // dummy to make compiler happy :)
}

bool XFEM::XFluid_Contact_Comm::Get_Contact_State(int sid,      // Solid Surface Element
    std::string mcname, const CORE::LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
    const double& full_fsi_traction,                            // stressfluid + penalty ...
    double& gap)
{
  if (!ele_ptrs_already_setup_) return true;

  int startgid =
      condition_manager_->GetMeshCouplingStartGID(condition_manager_->GetCouplingIndex(mcname));
  CONTACT::Element* cele = GetContactEle(sid + startgid);
  if (cele)
    return CheckNitscheContactState(cele, xsi, full_fsi_traction, gap);
  else
  {
    gap = 1e12;
    return true;
  }
}

void XFEM::XFluid_Contact_Comm::Get_States(const int fluidele_id, const std::vector<int>& fluid_nds,
    const DRT::ELEMENTS::StructuralSurface* sele, const CORE::LINALG::Matrix<2, 1>& selexsi,
    const CORE::LINALG::Matrix<3, 1>& x, DRT::Element*& fluidele,
    CORE::LINALG::SerialDenseMatrix& ele_xyze, std::vector<double>& velpres,
    std::vector<double>& disp, std::vector<double>& ivel, double& pres_m,
    CORE::LINALG::Matrix<3, 1>& vel_m, CORE::LINALG::Matrix<3, 1>& vel_s,
    CORE::LINALG::Matrix<3, 3>& vderxy_m, CORE::LINALG::Matrix<3, 1>& velpf_s)
{
  fluidele = fluiddis_->gElement(fluidele_id);
  DRT::ELEMENTS::Fluid* ffluidele = dynamic_cast<DRT::ELEMENTS::Fluid*>(fluidele);
  // 1 // get element states
  {
    DRT::Element::LocationArray laf(1);
    fluidele->LocationVector(*fluiddis_, fluid_nds, laf, false);
    Teuchos::RCP<const Epetra_Vector> matrix_state = fluiddis_->GetState("velaf");
    CORE::FE::ExtractMyValues(*matrix_state, velpres, laf[0].lm_);

    std::vector<int> lmdisp;
    lmdisp.resize(fluid_nds.size() * 3);
    if (!ffluidele) dserror("Cast to Fluidelement failed");
    if (ffluidele->IsAle())
    {
      for (std::size_t n = 0; n < fluid_nds.size(); ++n)
        for (int dof = 0; dof < 3; ++dof) lmdisp[n * 3 + dof] = laf[0].lm_[n * 4 + dof];
      Teuchos::RCP<const Epetra_Vector> matrix_state_disp = fluiddis_->GetState("dispnp");
      CORE::FE::ExtractMyValues(*matrix_state_disp, disp, lmdisp);
    }
  }
  {
    DRT::Element::LocationArray las(1);
    sele->LocationVector(*mc_[mcidx_]->GetCutterDis(), las, false);
    Teuchos::RCP<const Epetra_Vector> matrix_state =
        mc_[mcidx_]->GetCutterDis()->GetState("ivelnp");
    CORE::FE::ExtractMyValues(*matrix_state, ivel, las[0].lm_);
  }
  static std::vector<double> ipfvel;
  if (isporo_)
  {
    DRT::Element::LocationArray las(1);
    sele->LocationVector(*mcfpi_ps_pf_->GetCutterDis(), las, false);
    Teuchos::RCP<const Epetra_Vector> matrix_state =
        mcfpi_ps_pf_->GetCutterDis()->GetState("ivelnp");
    CORE::FE::ExtractMyValues(*matrix_state, ipfvel, las[0].lm_);
  }

  // 2 // get element xyze
  /// element coordinates in EpetraMatrix
  ele_xyze.shape(3, fluidele->NumNode());
  for (int i = 0; i < fluidele->NumNode(); ++i)
  {
    for (int j = 0; j < 3; j++)
    {
      if (ffluidele->IsAle())
        ele_xyze(j, i) = fluidele->Nodes()[i]->X()[j] + disp[i * 3 + j];
      else
        ele_xyze(j, i) = fluidele->Nodes()[i]->X()[j];
    }
  }

  // 3 // get quantities in gp
  {
    CORE::LINALG::Matrix<3, 1> fluidele_xsi(true);
    if (fluidele->Shape() == CORE::FE::CellType::hex8)
    {
      CORE::LINALG::Matrix<3, 8> xyze(ele_xyze.values(), true);
      // find element local position of gauss point
      Teuchos::RCP<CORE::GEO::CUT::Position> pos =
          CORE::GEO::CUT::PositionFactory::BuildPosition<3, CORE::FE::CellType::hex8>(xyze, x);
      if (!pos->Compute(1e-1))  // if we are a litte bit outside of the element we don't care ...
      {
        pos->LocalCoordinates(fluidele_xsi);
        std::cout << "fluidele_xsi: " << fluidele_xsi << std::endl;
        std::ofstream file("DEBUG_OUT_D007.pos");
        CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The point");
        CORE::GEO::CUT::OUTPUT::GmshCoordDump(file, x, -1);
        CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The Element", true);
        file << "SH(";
        for (int i = 0; i < 7; ++i)
          for (int j = 0; j < 3; ++j) file << xyze(j, i) << ", ";
        file << xyze(0, 7) << ", " << xyze(1, 7) << ", " << xyze(2, 7) << "){1,2,3,4,5,6,7,8};";
        CORE::GEO::CUT::OUTPUT::GmshEndSection(file, true);
        dserror("Couldn'd compute local coordinate for fluid element (DEBUG_OUT_D007.pos)!");
      }
      pos->LocalCoordinates(fluidele_xsi);

      static CORE::LINALG::Matrix<3, 3> xji;
      static CORE::LINALG::Matrix<3, 3> xjm;
      static CORE::LINALG::Matrix<3, 8> deriv;
      static CORE::LINALG::Matrix<8, 1> funct;
      static CORE::LINALG::Matrix<3, 8> derxy;

      // evaluate shape functions
      CORE::FE::shape_function<CORE::FE::CellType::hex8>(fluidele_xsi, funct);

      // evaluate the derivatives of shape functions
      CORE::FE::shape_function_deriv1<CORE::FE::CellType::hex8>(fluidele_xsi, deriv);
      xjm.MultiplyNT(deriv, xyze);
      // double det = xji.Invert(xjm); //if we need this at some point
      xji.Invert(xjm);

      // compute global first derivates
      derxy.Multiply(xji, deriv);

      static CORE::LINALG::Matrix<3, 8> vel;
      static CORE::LINALG::Matrix<8, 1> pres;
      for (int n = 0; n < fluidele->NumNode(); ++n)
      {
        pres(n, 0) = velpres[n * 4 + 3];
        for (int dof = 0; dof < 3; ++dof)
        {
          vel(dof, n) = velpres[n * 4 + dof];
        }
      }
      pres_m = pres.Dot(funct);
      vel_m.Multiply(1., vel, funct, 0.);
      vderxy_m.MultiplyNT(vel, derxy);
    }
    else
      dserror("fluidele is not hex8!");
  }

  // 4 // evaluate slave velocity at guasspoint
  if (sele->Shape() == CORE::FE::CellType::quad4)
  {
    static CORE::LINALG::Matrix<3, 4> vels;
    static CORE::LINALG::Matrix<3, 4> velpfs;
    for (int n = 0; n < sele->NumNode(); ++n)
    {
      for (int dof = 0; dof < 3; ++dof)
      {
        vels(dof, n) = ivel[n * 3 + dof];
        if (isporo_) velpfs(dof, n) = ipfvel[n * 3 + dof];
      }
    }

    const int numnodes = CORE::FE::num_nodes<CORE::FE::CellType::quad4>;
    static CORE::LINALG::Matrix<numnodes, 1> funct(false);
    CORE::FE::shape_function_2D(funct, selexsi(0), selexsi(1), CORE::FE::CellType::quad4);
    vel_s.Multiply(vels, funct);
    if (isporo_) velpf_s.Multiply(velpfs, funct);
  }
  else
    dserror("Your Slave Element is not a quad4?!");

  return;
}

void XFEM::XFluid_Contact_Comm::Get_Penalty_Param(DRT::Element* fluidele,
    CORE::GEO::CUT::VolumeCell* volumecell, CORE::LINALG::SerialDenseMatrix& ele_xyze,
    const CORE::LINALG::Matrix<3, 1>& elenormal, double& penalty_fac,
    const CORE::LINALG::Matrix<3, 1>& vel_m)
{
  double h_k;
  double inv_h_k;
  if (last_ele_h_.first != fluidele->Id() ||
      visc_stab_hk_ != INPAR::XFEM::ViscStab_hk_ele_vol_div_by_max_ele_surf)
  {
    // 1 // Get Boundary Cells and Gausspoints of this Boundarycells for this fluid element!
    std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>> bcells;
    std::map<int, std::vector<CORE::FE::GaussIntegration>> bintpoints;

    CORE::GEO::CUT::ElementHandle* cele = cutwizard_->GetElement(fluidele);
    if (cele == nullptr) dserror("Couldn't find cut element for ele %d", fluidele->Id());

    std::vector<CORE::GEO::CUT::plain_volumecell_set> cell_sets;
    {
      std::vector<std::vector<int>> nds_sets;
      std::vector<std::vector<CORE::FE::GaussIntegration>> intpoints_sets;
      cele->GetCellSets_DofSets_GaussPoints(cell_sets, nds_sets, intpoints_sets, false);
    }

    // find the right cell_set
    std::vector<CORE::GEO::CUT::plain_volumecell_set>::iterator sit = cell_sets.end();
    // if we have just one dof in this element this isn't a loop
    for (std::vector<CORE::GEO::CUT::plain_volumecell_set>::iterator s = cell_sets.begin();
         s != cell_sets.end(); s++)
    {
      CORE::GEO::CUT::plain_volumecell_set& cells = *s;
      for (CORE::GEO::CUT::plain_volumecell_set::iterator i = cells.begin(); i != cells.end(); ++i)
      {
        if (volumecell == (*i))
        {
          sit = s;
          break;
        }
        if (sit != cell_sets.end()) break;
      }
    }

    if (sit == cell_sets.end()) dserror("Couldn't identify a cell set!");

    CORE::GEO::CUT::plain_volumecell_set& cells = *sit;
    std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>> element_bcells;
    for (CORE::GEO::CUT::plain_volumecell_set::iterator i = cells.begin(); i != cells.end(); ++i)
    {
      CORE::GEO::CUT::VolumeCell* vc = *i;
      vc->GetBoundaryCellsToBeIntegrated(element_bcells);
    }

    for (std::map<int, std::vector<CORE::GEO::CUT::BoundaryCell*>>::const_iterator bc =
             element_bcells.begin();
         bc != element_bcells.end(); ++bc)
    {
      std::vector<CORE::GEO::CUT::BoundaryCell*>& bc_new = bcells[bc->first];
      bc_new.clear();
      std::copy(bc->second.begin(), bc->second.end(), std::inserter(bc_new, bc_new.end()));
    }

    if (bcells.size() > 0)
    {
      // get boundary cell Gaussian points
      cele->BoundaryCellGaussPointsLin(bcells, bintpoints, cutwizard_->Get_BC_Cubaturedegree());
    }
    else
    {
      std::cout << "I didn't identify any boundary cells, this happens if I'm outside of fluid "
                   "solution elements! --> penalty ~ nue/h = nue/V*O = 0"
                << std::endl;
      penalty_fac = 0.0;
      return;
    }
    if (fluidele->Shape() != CORE::FE::CellType::hex8) dserror("Add hex8 shapes here!");

    h_k = XFEM::UTILS::ComputeCharEleLength<CORE::FE::CellType::hex8>(
        fluidele, ele_xyze, condition_manager_, cells, bcells, bintpoints, visc_stab_hk_);

    inv_h_k = 1.0 / h_k;
    last_ele_h_ = std::pair<int, double>(fluidele->Id(), h_k);
  }
  else
  {
    h_k = last_ele_h_.second;
    inv_h_k = 1. / h_k;
  }

  // 3 // Compute Penalty Param
  double dummy;
  double kappa_m = 1.0;
  double kappa_s = 1.0 - kappa_m;
  double visc_stab_fac = 0.0;
  mc_[mcidx_]->Get_ViscPenalty_Stabfac(fluidele, nullptr, kappa_m, kappa_s, inv_h_k, visc_stab_fac,
      dummy, nit_stab_gamma_, nit_stab_gamma_, is_pseudo_2D_, visc_stab_trace_estimate_);

  Teuchos::RCP<MAT::Material> mat;
  XFEM::UTILS::GetVolumeCellMaterial(fluidele, mat);
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
  if (actmat == nullptr) dserror("Cast of Fluidmat failed!");

  XFEM::UTILS::NIT_Compute_FullPenalty_Stabfac(
      penalty_fac,  ///< to be filled: full Nitsche's penalty term scaling (viscous+convective part)
      elenormal, h_k,
      kappa_m,  // weights (only existing for Nitsche currently!!)
      kappa_s,  // weights (only existing for Nitsche currently!!)
      vel_m, vel_m,
      visc_stab_fac,  ///< Nitsche's viscous scaling part of penalty term
      theta_ * dt_, false, actmat->Density(), actmat->Density(), mass_conservation_scaling_,
      mass_conservation_combination_, nit_stab_gamma_, INPAR::XFEM::ConvStabScaling_none,
      INPAR::XFEM::XFF_ConvStabScaling_none, false, false);
  return;
}

void XFEM::XFluid_Contact_Comm::Get_Penalty_Param(
    DRT::ELEMENTS::StructuralSurface* sele, double& penalty_fac)
{
  penalty_fac = nit_stab_gamma_ *
                Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSI>(mc_[mcidx_])->GetTimeFac() *
                Teuchos::rcp_dynamic_cast<XFEM::MeshCouplingFSI>(mc_[mcidx_])
                    ->Get_EstimateNitscheTraceMaxEigenvalue(sele);
  return;
}

void XFEM::XFluid_Contact_Comm::SetupSurfElePtrs(DRT::Discretization& contact_interface_dis)
{
  if (ele_ptrs_already_setup_) return;

  if (!fluid_init_) dserror("fluid not initialized");

  if (!mc_.size()) dserror("SetupSurfElePtrs: Don't have any mesh coupling objects ....!");

  int startgid = condition_manager_->GetMeshCouplingStartGID(
      condition_manager_->GetCouplingIndex(mc_[0]->GetName()));
  min_surf_id_ = mc_[0]->GetCutterDis()->ElementColMap()->MinAllGID() + startgid;
  int max_surf_gid = mc_[0]->GetCutterDis()->ElementColMap()->MaxAllGID() + startgid;
  for (std::size_t mc = 1; mc < mc_.size(); ++mc)
  {
    int startgid = condition_manager_->GetMeshCouplingStartGID(
        condition_manager_->GetCouplingIndex(mc_[mc]->GetName()));
    if (min_surf_id_ > mc_[mc]->GetCutterDis()->ElementColMap()->MinAllGID() + startgid)
      min_surf_id_ = mc_[mc]->GetCutterDis()->ElementColMap()->MinAllGID() + startgid;
    if (max_surf_gid < mc_[mc]->GetCutterDis()->ElementColMap()->MaxAllGID() + startgid)
      max_surf_gid = mc_[mc]->GetCutterDis()->ElementColMap()->MaxAllGID() + startgid;
  }
  if (mcfpi_ps_pf_ != Teuchos::null)
  {
    int startgid = condition_manager_->GetMeshCouplingStartGID(
        condition_manager_->GetCouplingIndex(mcfpi_ps_pf_->GetName()));
    if (min_surf_id_ > mcfpi_ps_pf_->GetCutterDis()->ElementColMap()->MinAllGID() + startgid)
      min_surf_id_ = mcfpi_ps_pf_->GetCutterDis()->ElementColMap()->MinAllGID() + startgid;
    if (max_surf_gid < mcfpi_ps_pf_->GetCutterDis()->ElementColMap()->MaxAllGID() + startgid)
      max_surf_gid = mcfpi_ps_pf_->GetCutterDis()->ElementColMap()->MaxAllGID() + startgid;
  }
  soSurfId_to_mortar_ele_.resize(max_surf_gid - min_surf_id_ + 1, nullptr);
  min_mortar_id_ = contact_interface_dis.ElementColMap()->MinAllGID();
  const int max_mortar_gid = contact_interface_dis.ElementColMap()->MaxAllGID();
  mortarId_to_soSurf_ele_.resize(max_mortar_gid - min_mortar_id_ + 1, nullptr);
  mortarId_to_somc_.resize(max_mortar_gid - min_mortar_id_ + 1, -1);
  mortarId_to_sosid_.resize(max_mortar_gid - min_mortar_id_ + 1, -1);

  for (int i = 0; i < contact_interface_dis.ElementColMap()->NumMyElements(); ++i)
  {
    CONTACT::Element* cele = dynamic_cast<CONTACT::Element*>(contact_interface_dis.lColElement(i));
    if (!cele) dserror("no contact element or no element at all");
    const int c_parent_id = cele->ParentElementId();
    const int c_parent_surf = cele->FaceParentNumber();

    for (std::size_t mc = 0; mc < mc_.size(); ++mc)
    {
      for (int j = 0; j < mc_[mc]->GetCutterDis()->NumMyColElements(); ++j)
      {
        int startgid = condition_manager_->GetMeshCouplingStartGID(
            condition_manager_->GetCouplingIndex(mc_[mc]->GetName()));
        DRT::ELEMENTS::StructuralSurface* fele = dynamic_cast<DRT::ELEMENTS::StructuralSurface*>(
            mc_[mc]->GetCutterDis()->lColElement(j));
        if (!fele) dserror("no face element or no element at all");
        const int f_parent_id = fele->ParentElementId();
        const int f_parent_surf = fele->FaceParentNumber();
        if (c_parent_id == f_parent_id && c_parent_surf == f_parent_surf)
        {
          mortarId_to_soSurf_ele_[cele->Id() - min_mortar_id_] = fele;
          mortarId_to_somc_[cele->Id() - min_mortar_id_] = mc;
          mortarId_to_sosid_[cele->Id() - min_mortar_id_] = fele->Id() + startgid;

          soSurfId_to_mortar_ele_[fele->Id() - min_surf_id_ + startgid] = cele;
        }
      }
    }
    if (mcfpi_ps_pf_ != Teuchos::null)
    {
      for (int j = 0; j < mcfpi_ps_pf_->GetCutterDis()->NumMyColElements(); ++j)
      {
        int startgid = condition_manager_->GetMeshCouplingStartGID(
            condition_manager_->GetCouplingIndex(mcfpi_ps_pf_->GetName()));
        DRT::ELEMENTS::StructuralSurface* fele = dynamic_cast<DRT::ELEMENTS::StructuralSurface*>(
            mcfpi_ps_pf_->GetCutterDis()->lColElement(j));
        if (!fele) dserror("no face element or no element at all");
        const int f_parent_id = fele->ParentElementId();
        const int f_parent_surf = fele->FaceParentNumber();

        if (c_parent_id == f_parent_id && c_parent_surf == f_parent_surf)
        {
          // mortarId_to_soSurf_ele_[cele->Id()-min_mortar_id_]=fele; //we dont need this connection
          // as this is already available from ps_ps
          soSurfId_to_mortar_ele_[fele->Id() - min_surf_id_ + startgid] = cele;
        }
      }
    }
  }
  ele_ptrs_already_setup_ = true;

  contact_strategy_fsi_ =
      dynamic_cast<CONTACT::NitscheStrategyFsi*>(&contact_strategy_);  // might be nullptr
  contact_strategy_fpi_ =
      dynamic_cast<CONTACT::NitscheStrategyFpi*>(&contact_strategy_);  // might be nullptr
}

bool XFEM::XFluid_Contact_Comm::GetVolumecell(DRT::ELEMENTS::StructuralSurface*& sele,
    CORE::LINALG::Matrix<2, 1>& xsi, CORE::GEO::CUT::SideHandle*& sidehandle, std::vector<int>& nds,
    int& eleid, CORE::GEO::CUT::VolumeCell*& volumecell, CORE::LINALG::Matrix<3, 1>& elenormal,
    CORE::LINALG::Matrix<3, 1>& x, bool& FSI_integrated, double& distance)
{
  distance = 0.0;
  FSI_integrated = true;
  // 1 // Compute global coord x
  volumecell = nullptr;
  if (sele->Shape() == CORE::FE::CellType::quad4)
  {
    const int numnodes = CORE::FE::num_nodes<CORE::FE::CellType::quad4>;

    CORE::LINALG::SerialDenseMatrix xyze_m;
    CORE::LINALG::Matrix<numnodes, 1> funct(false);

    sidehandle->Coordinates(xyze_m);
    CORE::LINALG::Matrix<3, numnodes> xyze(xyze_m.values(), true);
    CORE::FE::shape_function_2D(funct, xsi(0), xsi(1), CORE::FE::CellType::quad4);
    x.Multiply(xyze, funct);
  }
  else
    dserror("GetFacet: Your solid face is not a quad4, please add your element type here!");

  // 2 //Identify Subside
  CORE::GEO::CUT::Facet* facet = nullptr;

  CORE::GEO::CUT::plain_side_set subsides;
  sidehandle->CollectSides(subsides);
  int found_side = -1;
  static CORE::LINALG::Matrix<3, 1> tmpxsi(true);
  static CORE::LINALG::Matrix<3, 1> tmpxsi_tmp(true);
  for (std::size_t ss = 0; ss < subsides.size(); ++ss)
  {
    CORE::GEO::CUT::Side* s = subsides[ss];
    if (s->LocalCoordinates(x, tmpxsi_tmp, true, 1e-10))
    {
      found_side = ss;
      tmpxsi = tmpxsi_tmp;
      CORE::LINALG::Matrix<2, 1> tmp2xsi(tmpxsi.A(), true);
      s->Normal(tmp2xsi, elenormal, true);
      elenormal.Scale(-1.0);          // flip direction
      if (fabs(tmpxsi(2, 0)) < 1e-1)  // do not search for better canidates anymore
        break;
    }
  }
  if (found_side == -1)  // do the same thing with reduced tolerance
  {
    for (std::size_t ss = 0; ss < subsides.size(); ++ss)
    {
      CORE::GEO::CUT::Side* s = subsides[ss];
      if (s->LocalCoordinates(x, tmpxsi_tmp, true, 1e-6))
      {
        found_side = ss;
        tmpxsi = tmpxsi_tmp;
        CORE::LINALG::Matrix<2, 1> tmp2xsi(tmpxsi.A(), true);
        s->Normal(tmp2xsi, elenormal, true);
        elenormal.Scale(-1.0);          // flip direction
        if (fabs(tmpxsi(2, 0)) < 1e-1)  // do not search for better canidates anymore
          break;
      }
    }
  }
  if (found_side == -1)
  {
    std::ofstream file("DEBUG_OUT_D001.pos");
    CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The point to idenity");
    CORE::GEO::CUT::OUTPUT::GmshCoordDump(file, x, -1);
    CORE::GEO::CUT::OUTPUT::GmshEndSection(file);
    CORE::GEO::CUT::OUTPUT::GmshWriteSection(file, "All subsides", subsides, true);
    dserror("Coundn't identify side (DEBUG_OUT_D001.pos)!");
  }
  else
  {
    std::vector<CORE::GEO::CUT::Facet*> facets = subsides[found_side]->Facets();

    CORE::GEO::CUT::Side* side = subsides[found_side];
    // Handle here unphysical sides ...
    if (sidehandle->IsunphysicalSubSide(subsides[found_side]))
    {
      FSI_integrated = false;

      // Find Closest Point on physical boundary ...
      side = FindnextPhysicalSide(x, subsides[found_side], sidehandle, xsi, distance);
      side->Normal(xsi, elenormal, true);
      elenormal.Scale(-1.0);  // flip direction
      sele =
          dynamic_cast<DRT::ELEMENTS::StructuralSurface*>(condition_manager_->GetSide(side->Id()));
      if (!sele) dserror("Couldn't Identify new sele %d", side->Id());
      facets = side->Facets();
    }

    if (facets.size() == 1 && !parallel_)
    {
      facet = facets[0];
    }
    else if (facets.size() > 1 || parallel_)
    {
      // a //simple case --> still all facets belong to the same volumecells
      bool SameVCs = true;
      for (std::vector<CORE::GEO::CUT::Facet*>::iterator fit = facets.begin(); fit != facets.end();
           ++fit)
      {
        if (facets[0]->Cells().size() != (*fit)->Cells().size())  // simplest check
        {
          SameVCs = false;
          break;
        }
        else
        {
          for (std::size_t checkcell = 0; checkcell < (*fit)->Cells().size(); ++checkcell)
          {
            if (facets[0]->Cells()[checkcell] != (*fit)->Cells()[checkcell])
            {
              SameVCs = false;
              break;
            }
          }
        }
      }

      if (SameVCs && !parallel_)  // we can take any facet ...
      {
        facet = facets[0];
      }
      else  // b // more complex need to really identify the facet ...
      {
        for (std::vector<CORE::GEO::CUT::Facet*>::iterator fit = facets.begin();
             fit != facets.end(); ++fit)
        {
          CORE::GEO::CUT::Facet* afacet = *fit;
          std::vector<std::vector<CORE::GEO::CUT::Point*>> triangulation;
          if (afacet->Triangulation().size())
            triangulation = afacet->Triangulation();
          else if (afacet->GetSplitCells().size())
            triangulation = afacet->GetSplitCells();
          else if (afacet->NumPoints() == 3)
            triangulation.push_back(afacet->Points());
          else
          {
            std::ofstream file("DEBUG_OUT_D003.pos");
            CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The point to idenity");
            CORE::GEO::CUT::OUTPUT::GmshCoordDump(file, x, -1);
            CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "facet", true);
            CORE::GEO::CUT::OUTPUT::GmshFacetDump(file, afacet, "sides", true);
            CORE::GEO::CUT::OUTPUT::GmshEndSection(file);
            CORE::GEO::CUT::OUTPUT::GmshWriteSection(file, "ALLFACETS", facets, true);
            std::cout << "==| Warning from of your friendly XFluid_Contact_Comm: I have untriagled "
                         "faces (DEBUG_OUT_D003.pos)! |=="
                      << std::endl;
          }

          for (std::size_t tri = 0; tri < triangulation.size(); ++tri)
          {
            if (triangulation[tri].size() != 3)
              dserror("Triangulation with another number of points than 3 (%d)?",
                  triangulation[tri].size());

            // Compute local coords and take first possible facet ...
            CORE::LINALG::Matrix<3, 3> xyzf;
            triangulation[tri][0]->Coordinates(xyzf.A());
            triangulation[tri][1]->Coordinates(xyzf.A() + 3);
            triangulation[tri][2]->Coordinates(xyzf.A() + 6);

            Teuchos::RCP<CORE::GEO::CUT::Position> pos =
                CORE::GEO::CUT::PositionFactory::BuildPosition<3, CORE::FE::CellType::tri3>(
                    xyzf, x);
            bool success = pos->Compute(1e-6, true);
            if (success)
            {
              pos->LocalCoordinates(tmpxsi);
              if (fabs(tmpxsi(2, 0)) > 1e-3)
                dserror("To far away from this facet %f!", tmpxsi(2, 0));
              facet = afacet;
              break;
            }
          }
          if (facet) break;
        }

        if (!facet && parallel_)
        {
          return false;  // in parallel this is ok
        }
        else if (!facet)
        {
          std::ofstream file("DEBUG_OUT_D004.pos");
          CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The point to idenity");
          CORE::GEO::CUT::OUTPUT::GmshCoordDump(file, x, -1);
          CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The selected subside", true);
          CORE::GEO::CUT::OUTPUT::GmshSideDump(file, subsides[found_side]);
          CORE::GEO::CUT::OUTPUT::GmshEndSection(file);
          CORE::GEO::CUT::OUTPUT::GmshWriteSection(file, "All facets", facets);
          for (std::size_t f = 0; f < facets.size(); ++f)
          {
            std::stringstream strf;
            strf << "Facet (" << f << ")";
            std::stringstream strv;
            strv << "Volumecell of facet (" << f << ")";
            CORE::GEO::CUT::OUTPUT::GmshWriteSection(file, strf.str(), facets[f]);
            CORE::GEO::CUT::OUTPUT::GmshWriteSection(file, strv.str(), facets[f]->Cells());
          }
          file.close();

          dserror("Couldn't identify facet (DEBUG_OUT_D004.pos)!");
        }
      }
    }
    else
    {
      if (parallel_)
      {
        return false;  // in parallel this is ok
      }

      std::ofstream file("DEBUG_OUT_D005.pos");
      CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The point to idenity");
      CORE::GEO::CUT::OUTPUT::GmshCoordDump(file, x, -1);
      CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The selected subside", true);
      CORE::GEO::CUT::OUTPUT::GmshSideDump(file, side);
      CORE::GEO::CUT::OUTPUT::GmshEndSection(file, true);
      dserror("Side has no facets but is physical (DEBUG_OUT_D005.pos)!");
    }
  }

  // 3 // Get Volumecells
  if (!volumecell)
  {
    for (std::size_t vc = 0; vc < facet->Cells().size(); ++vc)
    {
      if (facet->Cells()[vc]->Position() == CORE::GEO::CUT::Point::outside)
      {
        if (volumecell)
        {
          std::ofstream file("DEBUG_OUT_D006.pos");
          CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The point to idenity");
          CORE::GEO::CUT::OUTPUT::GmshCoordDump(file, x, -1);
          CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "The selected subside", true);
          CORE::GEO::CUT::OUTPUT::GmshSideDump(file, subsides[found_side]);
          CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "VC1", true);
          CORE::GEO::CUT::OUTPUT::GmshVolumecellDump(file, volumecell, "sides", true);
          CORE::GEO::CUT::OUTPUT::GmshNewSection(file, "VC2", true);
          CORE::GEO::CUT::OUTPUT::GmshVolumecellDump(file, facet->Cells()[vc], "sides", true);
          file.close();
          dserror("Facet has at least two volumecells which are outside (DEBUG_OUT_D006.pos)!");
        }
        volumecell = facet->Cells()[vc];
        if (parallel_)
        {
          if (fluiddis_->ElementRowMap()->LID(volumecell->ParentElement()->Id()) ==
              -1)  // fluidele not owned by this proc
          {
            return false;
          }
        }
      }
    }
    if (!volumecell) dserror("Facet has no volumecell which is outside?");
  }

  nds = volumecell->NodalDofSet();
  eleid = volumecell->GetParentElementId();

  return true;
}

CORE::GEO::CUT::Side* XFEM::XFluid_Contact_Comm::FindnextPhysicalSide(CORE::LINALG::Matrix<3, 1>& x,
    CORE::GEO::CUT::Side* initSide, CORE::GEO::CUT::SideHandle*& sidehandle,
    CORE::LINALG::Matrix<2, 1>& newxsi, double& distance)
{
  std::set<CORE::GEO::CUT::Side*> performed_sides;
  std::set<CORE::GEO::CUT::Side*> physical_sides;
  if (last_physical_sides_.first != initSide->Id())
  {
    Update_physical_sides(initSide, performed_sides, physical_sides);
    last_physical_sides_ =
        std::pair<int, std::set<CORE::GEO::CUT::Side*>>(initSide->Id(), physical_sides);
  }
  else
    physical_sides = last_physical_sides_.second;

  distance = 1e200;
  CORE::LINALG::Matrix<3, 1> newx(true);
  CORE::GEO::CUT::Side* newSide = nullptr;

  for (std::set<CORE::GEO::CUT::Side*>::iterator psit = physical_sides.begin();
       psit != physical_sides.end(); ++psit)
  {
    static CORE::LINALG::Matrix<3, 1> tmpx(true);
    double tmpdistance = DistancetoSide(x, *psit, tmpx);
    if (distance > tmpdistance)
    {
      distance = tmpdistance;
      newx.Update(1.0, tmpx, 0.0);
      newSide = *psit;
    }
  }

  if (!newSide)
  {
    std::stringstream str;
    str << "CINS_" << fluiddis_->Comm().MyPID() << ".pos";
    std::ofstream file(str.str().c_str());
    CORE::GEO::CUT::OUTPUT::GmshWriteSection(file, "InitSide", initSide);
    CORE::GEO::CUT::OUTPUT::GmshWriteSection(file, "PerFormedSides", performed_sides, true);
    dserror("Couldn't identify a new side (number of identified physical sides: %d)!",
        physical_sides.size());
  }

  // Update global position...
  x.Update(1.0, newx, 0.0);

  // compute again surface local coords ...
  sidehandle = cutwizard_->GetCutSide(newSide->Id());
  if (!sidehandle)
  {
    std::cout << "The Side pointer is " << newSide << std::endl;
    std::cout << "Couldn't get Sidehandle for side " << newSide->Id() << std::endl;
    // newSide->Print();
    dserror("Couldn't get Sidehandle for side %f", newSide->Id());
  }

  CORE::LINALG::SerialDenseMatrix xyzs;
  sidehandle->Coordinates(xyzs);
  if (sidehandle->Shape() == CORE::FE::CellType::quad4)
  {
    CORE::LINALG::Matrix<3, 4> xyze(xyzs.values(), true);
    Teuchos::RCP<CORE::GEO::CUT::Position> pos =
        CORE::GEO::CUT::PositionFactory::BuildPosition<3, CORE::FE::CellType::quad4>(xyze, newx);
    pos->Compute(1e-15, true);
    pos->LocalCoordinates(newxsi);
  }
  else
    dserror("Not a quad4!");
  return newSide;
}

double XFEM::XFluid_Contact_Comm::DistancetoSide(CORE::LINALG::Matrix<3, 1>& x,
    CORE::GEO::CUT::Side* side, CORE::LINALG::Matrix<3, 1>& closest_x)
{
  double distance = 1e200;
  for (std::vector<CORE::GEO::CUT::Edge*>::const_iterator eit = side->Edges().begin();
       eit != side->Edges().end(); ++eit)
  {
    CORE::GEO::CUT::Edge* e = *eit;
    CORE::LINALG::Matrix<3, 2> xyzl;
    e->Coordinates(xyzl);
    Teuchos::RCP<CORE::GEO::CUT::Position> pos =
        CORE::GEO::CUT::PositionFactory::BuildPosition<3, CORE::FE::CellType::line2>(xyzl, x);
    pos->Compute(true);
    CORE::LINALG::Matrix<1, 1> rst;
    pos->LocalCoordinates(rst);
    if (fabs(rst(0.0)) < 1)
    {
      static CORE::LINALG::Matrix<2, 1> funct;
      // evaluate shape functions
      CORE::FE::shape_function<CORE::FE::CellType::line2>(rst, funct);
      static CORE::LINALG::Matrix<3, 1> posx;
      posx.Multiply(xyzl, funct);
      posx.Update(-1, x, 1);
      double tmpdistance = posx.Norm2();
      if (distance > tmpdistance)
      {
        distance = tmpdistance;
        closest_x.Multiply(xyzl, funct);
      }
    }
  }
  for (std::vector<CORE::GEO::CUT::Node*>::const_iterator nit = side->Nodes().begin();
       nit != side->Nodes().end(); ++nit)
  {
    CORE::GEO::CUT::Node* n = *nit;
    CORE::LINALG::Matrix<3, 1> xyzn;
    n->Coordinates(xyzn.A());
    xyzn.Update(-1, x, 1);
    double tmpdistance = xyzn.Norm2();
    if (distance > tmpdistance)
    {
      distance = tmpdistance;
      n->Coordinates(closest_x.A());
    }
  }
  return distance;
}

double XFEM::XFluid_Contact_Comm::Get_h() { return mc_[0]->Get_h(); }

void XFEM::XFluid_Contact_Comm::Update_physical_sides(CORE::GEO::CUT::Side* side,
    std::set<CORE::GEO::CUT::Side*>& performed_sides,
    std::set<CORE::GEO::CUT::Side*>& physical_sides)
{
  std::vector<CORE::GEO::CUT::Side*> neibs = GetNewNeighboringSides(side, performed_sides);
  for (std::size_t sid = 0; sid < neibs.size(); ++sid)
  {
    performed_sides.insert(neibs[sid]);
    CORE::GEO::CUT::SideHandle* sh = cutwizard_->GetCutSide(neibs[sid]->Id());
    if (!sh) dserror("Couldn't Get Sidehandle %d!", neibs[sid]->Id());
    if (sh->IsunphysicalSubSide(neibs[sid]))
      Update_physical_sides(neibs[sid], performed_sides, physical_sides);
    else
    {
      CORE::LINALG::Matrix<3, 1> normal;
      CORE::LINALG::Matrix<2, 1> center(true);
      neibs[sid]->Normal(center, normal, false);
      double norm = normal.Norm2();
      if (norm > 1e-10)
        physical_sides.insert(neibs[sid]);
      else
        Update_physical_sides(neibs[sid], performed_sides, physical_sides);
    }
  }
}

std::vector<CORE::GEO::CUT::Side*> XFEM::XFluid_Contact_Comm::GetNewNeighboringSides(
    CORE::GEO::CUT::Side* side, std::set<CORE::GEO::CUT::Side*>& performed_sides)
{
  std::vector<CORE::GEO::CUT::Side*> neighbors;
  for (std::vector<CORE::GEO::CUT::Node*>::const_iterator nit = side->Nodes().begin();
       nit != side->Nodes().end(); ++nit)
  {
    CORE::GEO::CUT::Node* n = *nit;
    for (CORE::GEO::CUT::plain_side_set::const_iterator sit = n->Sides().begin();
         sit != n->Sides().end(); ++sit)
    {
      CORE::GEO::CUT::Side* s = *sit;
      if (s == side) continue;
      if (s->Id() < 0) continue;
      if (performed_sides.find(s) != performed_sides.end()) continue;
      if (!IsRegisteredSurface(s->Id())) continue;
      if (s->Id() == side->Id())  // on the same contact element
        neighbors.push_back(s);
      else  // do the contact elements share a common edge?
      {
        int common_nodes = 0;
        for (std::size_t neighbor = 0; neighbor < neighbors.size(); ++neighbor)
        {
          if (neighbors[neighbor]->Id() == s->Id())
          {
            neighbors.push_back(s);
            common_nodes = -1;
            break;
          }
        }
        if (common_nodes == -1) continue;  // we assigned already a side with the same id

        CORE::LINALG::SerialDenseMatrix xyze1;
        CORE::LINALG::SerialDenseMatrix xyze2;
        CORE::GEO::CUT::SideHandle* sh1 = cutwizard_->GetCutSide(side->Id());
        CORE::GEO::CUT::SideHandle* sh2 = cutwizard_->GetCutSide(s->Id());

        for (std::size_t nidx1 = 0; nidx1 < sh1->GetNodes().size(); ++nidx1)
          for (std::size_t nidx2 = 0; nidx2 < sh2->GetNodes().size(); ++nidx2)
          {
            if (sh1->GetNodes()[nidx1]->Id() == sh2->GetNodes()[nidx2]->Id())
            {
              ++common_nodes;
              break;
            }
          }
        if (common_nodes == 2)
          neighbors.push_back(s);
        else if (common_nodes > 2)
        {
          std::cout << "==| Rejected side " << s->Id() << "( " << side->Id()
                    << ") as there are only " << common_nodes << " common nodes! |==" << std::endl;
          s->Print();
          side->Print();
        }
      }
    }
  }
  return neighbors;
}

CORE::GEO::CUT::Element* XFEM::XFluid_Contact_Comm::GetNextElement(CORE::GEO::CUT::Element* ele,
    std::set<CORE::GEO::CUT::Element*>& performed_elements, int& lastid)
{
  CORE::GEO::CUT::Element* newele = nullptr;
  if (lastid == -1 && ele != nullptr)
  {
    for (std::size_t s = 0; s < ele->Sides().size(); ++s)
    {
      for (std::size_t e = 0; e < ele->Sides()[s]->Elements().size(); ++e)
      {
        CORE::GEO::CUT::Element* element = ele->Sides()[s]->Elements()[e];
        if (element == ele)
          break;
        else if (performed_elements.find(element) != performed_elements.end())
          break;
        else
        {
          newele = element;
          break;
        }
      }
      if (newele) break;
    }
  }
  if (!newele)  // loop structure over all elements
  {
    while (!newele)
    {
      ++lastid;
      if (lastid > 2 * fluiddis_->ElementColMap()->NumGlobalElements())
      {
        return nullptr;
      }
      std::cout << "==| Doing the expensive Version of finding an element! |==" << lastid
                << std::endl;

      CORE::GEO::CUT::ElementHandle* elementh = cutwizard_->GetElement(lastid);
      if (elementh == nullptr) continue;
      CORE::GEO::CUT::plain_element_set pes;
      elementh->CollectElements(pes);

      if (pes.size() != 1) dserror("Collected Elements != 1");

      CORE::GEO::CUT::Element* element = pes[0];

      if (performed_elements.find(element) != performed_elements.end())
        continue;
      else
      {
        newele = element;
        break;
      }
    }
  }

  return newele;
}

void XFEM::XFluid_Contact_Comm::RegisterSideProc(int sid)
{
  if (!parallel_) return;
  if (GetContactEle(sid)) my_sele_ids_.insert(GetContactEle(sid)->Id());
}


void XFEM::XFluid_Contact_Comm::GetCutSideIntegrationPoints(
    int sid, CORE::LINALG::SerialDenseMatrix& coords, std::vector<double>& weights, int& npg)
{
  CORE::GEO::CUT::SideHandle* sh = cutwizard_->GetCutSide(GetSurfSid(sid));
  if (!sh) dserror("Couldn't get SideHandle!");
  if (sh->Shape() != CORE::FE::CellType::quad4) dserror("Not a quad4!");
  const int numnodes_sh = CORE::FE::num_nodes<CORE::FE::CellType::quad4>;
  CORE::LINALG::SerialDenseMatrix xquad;
  sh->Coordinates(xquad);
  CORE::LINALG::Matrix<2, numnodes_sh> deriv(false);
  CORE::LINALG::Matrix<2, 2> metrictensor(false);
  CORE::LINALG::Matrix<3, 1> normal_side(true);
  CORE::LINALG::Matrix<3, 1> normal_bc(true);

  CORE::GEO::CUT::plain_side_set subsides;
  sh->CollectSides(subsides);
  std::vector<Teuchos::RCP<CORE::GEO::CUT::Tri3BoundaryCell>> bcs;

  CORE::LINALG::SerialDenseMatrix tcoords(3, 3);
  for (CORE::GEO::CUT::plain_side_set::iterator sit = subsides.begin(); sit != subsides.end();
       ++sit)
  {
    CORE::GEO::CUT::Side* side = *sit;
    side->Normal(CORE::LINALG::Matrix<2, 1>(true), normal_side, true);
    for (std::vector<CORE::GEO::CUT::Facet*>::const_iterator fit = side->Facets().begin();
         fit != side->Facets().end(); ++fit)
    {
      CORE::GEO::CUT::Facet* facet = *fit;
      if (facet->IsTriangulated())
      {
        for (std::size_t triangle = 0; triangle < facet->Triangulation().size(); ++triangle)
        {
          double* coord = tcoords.values();
          for (std::vector<CORE::GEO::CUT::Point*>::const_iterator tp =
                   facet->Triangulation()[triangle].begin();
               tp != facet->Triangulation()[triangle].end(); ++tp)
          {
            CORE::GEO::CUT::Point* p = *tp;
            p->Coordinates(coord);
            coord += 3;
          }
          Teuchos::RCP<CORE::GEO::CUT::Tri3BoundaryCell> tmp_bc =
              Teuchos::rcp(new CORE::GEO::CUT::Tri3BoundaryCell(
                  tcoords, facet, facet->Triangulation()[triangle]));
          tmp_bc->Normal(CORE::LINALG::Matrix<2, 1>(true), normal_bc);
          if (normal_bc.Dot(normal_side) < 0.0)
            bcs.push_back(tmp_bc);
          else
          {
            std::vector<CORE::GEO::CUT::Point*> points = facet->Triangulation()[triangle];
            std::reverse(points.begin(), points.end());
            double* coord = tcoords.values();
            for (std::vector<CORE::GEO::CUT::Point*>::const_iterator tp =
                     facet->Triangulation()[triangle].end() - 1;
                 tp != facet->Triangulation()[triangle].begin() - 1; --tp)
            {
              CORE::GEO::CUT::Point* p = *tp;
              p->Coordinates(coord);
              coord += 3;
            }
            Teuchos::RCP<CORE::GEO::CUT::Tri3BoundaryCell> tmp_bc_rev =
                Teuchos::rcp(new CORE::GEO::CUT::Tri3BoundaryCell(tcoords, facet, points));
            bcs.push_back(tmp_bc_rev);
          }
        }
      }
      else if (facet->Points().size() == 3)
      {
        facet->Coordinates(tcoords.values());
        Teuchos::RCP<CORE::GEO::CUT::Tri3BoundaryCell> tmp_bc =
            Teuchos::rcp(new CORE::GEO::CUT::Tri3BoundaryCell(tcoords, facet, facet->Points()));
        tmp_bc->Normal(CORE::LINALG::Matrix<2, 1>(true), normal_bc);
        if (normal_bc.Dot(normal_side) < 0.0)
          bcs.push_back(tmp_bc);
        else
        {
          std::vector<CORE::GEO::CUT::Point*> points = facet->Points();
          std::reverse(points.begin(), points.end());
          double* coord = tcoords.values();
          for (std::vector<CORE::GEO::CUT::Point*>::const_iterator tp = facet->Points().end() - 1;
               tp != facet->Points().begin() - 1; --tp)
          {
            CORE::GEO::CUT::Point* p = *tp;
            p->Coordinates(coord);
            coord += 3;
          }
          Teuchos::RCP<CORE::GEO::CUT::Tri3BoundaryCell> tmp_bc_rev =
              Teuchos::rcp(new CORE::GEO::CUT::Tri3BoundaryCell(tcoords, facet, points));
          bcs.push_back(tmp_bc_rev);
        }
      }
      else
      {
        std::cout << "==| Ignore facet |==" << std::endl;
        std::cout << "facet->GetSplitCells().size(): " << facet->GetSplitCells().size()
                  << std::endl;
        facet->Print(std::cout);
        if (!parallel_) dserror("Ignore Facet");
      }
    }
    if (!side->Facets().size())
    {
      if (!sh->IsunphysicalSubSide(side))
      {
        if (parallel_)  // we are on the wrong proc ...
          continue;     // return true;
        else
          dserror("There are not facets on a physical side?");
      }
      else if (side->NumNodes() == 3)
      {
        side->Coordinates(tcoords.values());
        std::vector<CORE::GEO::CUT::Point*> points;
        for (unsigned p = 0; p < side->NumNodes(); ++p) points.push_back(side->Nodes()[p]->point());
        Teuchos::RCP<CORE::GEO::CUT::Tri3BoundaryCell> tmp_bc =
            Teuchos::rcp(new CORE::GEO::CUT::Tri3BoundaryCell(tcoords, nullptr, points));
        tmp_bc->Normal(CORE::LINALG::Matrix<2, 1>(true), normal_bc);
        if (normal_bc.Dot(normal_side) < 0.0)
          bcs.push_back(tmp_bc);
        else
        {
          std::vector<CORE::GEO::CUT::Point*> tmp_points = points;
          std::reverse(tmp_points.begin(), tmp_points.end());
          double* coord = tcoords.values();
          for (std::vector<CORE::GEO::CUT::Node*>::const_iterator tp = side->Nodes().end() - 1;
               tp != side->Nodes().begin() - 1; --tp)
          {
            CORE::GEO::CUT::Point* p = (*tp)->point();
            p->Coordinates(coord);
            coord += 3;
          }
          Teuchos::RCP<CORE::GEO::CUT::Tri3BoundaryCell> tmp_bc_rev =
              Teuchos::rcp(new CORE::GEO::CUT::Tri3BoundaryCell(tcoords, nullptr, tmp_points));
          bcs.push_back(tmp_bc_rev);
        }
      }
      else
        dserror("Unphysical Subside is not a tri3?");  // Cannot happen as these sides are created
                                                       // by the SelfCut
    }
  }

  weights.clear();
  CORE::LINALG::Matrix<3, 1> x_gp_lin(true);
  CORE::LINALG::Matrix<3, 1> normal(true);
  CORE::LINALG::Matrix<2, 1> rst(true);  // local coordinates w.r.t side
  double drs = 0;
  double drs_sh = 0;
  for (std::size_t bc = 0; bc < bcs.size(); ++bc)
  {
    CORE::FE::GaussIntegration gi = bcs[bc]->gaussRule(cutwizard_->Get_BC_Cubaturedegree());
    if (gi.NumPoints())
    {
      coords.reshape(weights.size() + gi.NumPoints(), 2);
      int idx = weights.size();
      for (CORE::FE::GaussIntegration::iterator iquad = gi.begin(); iquad != gi.end(); ++iquad)
      {
        const CORE::LINALG::Matrix<2, 1> eta(iquad.Point(), false);
        XFEM::UTILS::ComputeSurfaceTransformation(drs, x_gp_lin, normal, bcs[bc].getRawPtr(), eta);

        // find element local position of gauss point
        const CORE::LINALG::Matrix<3, numnodes_sh> xquad_m(xquad.values(), true);
        Teuchos::RCP<CORE::GEO::CUT::Position> pos =
            CORE::GEO::CUT::PositionFactory::BuildPosition<3, CORE::FE::CellType::quad4>(
                xquad_m, x_gp_lin);
        pos->Compute(true);
        pos->LocalCoordinates(rst);
        coords(idx, 0) = rst(0);
        coords(idx, 1) = rst(1);

        CORE::FE::shape_function_2D_deriv1(
            deriv, coords(idx, 0), coords(idx, 1), CORE::FE::CellType::quad4);
        CORE::FE::ComputeMetricTensorForBoundaryEle<CORE::FE::CellType::quad4>(
            xquad_m, deriv, metrictensor, drs_sh, nullptr);
        weights.push_back(iquad.Weight() * drs / drs_sh);  // small tri3 to quad4 weight
        ++idx;
      }
    }
  }

  npg = weights.size() - 1;
  return;
}

void XFEM::XFluid_Contact_Comm::FillComplete_SeleMap()
{
  if (!parallel_) return;
  if (cutwizard_ == Teuchos::null)
    dserror("XFluid_Contact_Comm::FillComplete_SeleMap: CutWizard not set!");

  // We also add all unphysical sides
  for (std::size_t i = 0; i < mortarId_to_sosid_.size(); ++i)
  {
    if (mortarId_to_sosid_[i] == -1) continue;  // this entry is not set!
    CORE::GEO::CUT::SideHandle* sh = cutwizard_->GetCutSide(mortarId_to_sosid_[i]);
    if (!sh)
      dserror("Couldn't get Sidhandle for mortarId %d, soid %d!", i + min_mortar_id_,
          mortarId_to_sosid_[i]);
    if (cutwizard_->GetCutSide(mortarId_to_sosid_[i])->HasunphysicalSubSide())
    {
      my_sele_ids_.insert(GetContactEle(mortarId_to_sosid_[i])->Id());
    }
  }
  std::vector<int> my_sele_ids(my_sele_ids_.begin(), my_sele_ids_.end());
  contact_ele_rowmap_fluidownerbased_ = Teuchos::rcp(
      new Epetra_Map(-1, my_sele_ids.size(), my_sele_ids.data(), 0, fluiddis_->Comm()));
}

void XFEM::XFluid_Contact_Comm::PrepareTimeStep()
{
  higher_contact_elements_.clear();
  higher_contact_elements_comm_.clear();
}

void XFEM::XFluid_Contact_Comm::PrepareIterationStep()
{
  higher_contact_elements_comm_.clear();

  std::vector<int> src(higher_contact_elements_.begin(), higher_contact_elements_.end());

  std::vector<int> dest;
  CORE::LINALG::AllreduceVector(src, dest, fluiddis_->Comm());

  higher_contact_elements_comm_.insert(dest.begin(), dest.end());

  if (!fluiddis_->Comm().MyPID() && higher_contact_elements_comm_.size())
  {
    std::cout << "==| Interface Elements with an increased number of Contact Gausspoins:";
    for (std::set<int>::iterator sit = higher_contact_elements_comm_.begin();
         sit != higher_contact_elements_comm_.end(); ++sit)
      std::cout << " " << GetSurfSid(*sit) << "(" << *sit << ") |";
    std::cout << "==" << std::endl;
  }
}

double XFEM::XFluid_Contact_Comm::Get_fpi_pcontact_exchange_dist()
{
  return mcfpi_ps_pf_->Get_fpi_pcontact_exchange_dist();
}

double XFEM::XFluid_Contact_Comm::Get_fpi_pcontact_fullfraction()
{
  return mcfpi_ps_pf_->Get_fpi_pcontact_fullfraction();
}

void XFEM::XFluid_Contact_Comm::Create_New_Gmsh_files()
{
#ifdef WRITE_GMSH
  std::vector<std::string> sections;
  sections.push_back("Contact_Traction");           // 0
  sections.push_back("FSI_Traction");               // 1
  sections.push_back("Contact_Active");             // 2
  sections.push_back("FSI_Active");                 // 3
  sections.push_back("Contact_Traction_Solid");     // 4
  sections.push_back("Contact_Traction_Fluid");     // 5
  sections.push_back("FSI_sliplenth");              // 6
  sections.push_back("All_GPs_Contact");            // 7
  sections.push_back("Contact_PoroFlow_Active");    // 8
  sections.push_back("Contact_PoroFlow_Inactive");  // 9
  sections.push_back("FPI_PoroFlow_Ffac");          // 10

  static int counter = 0;
  if (counter)
  {
    std::stringstream str;
    str << "FSCI_" << counter << "_" << fluiddis_->Comm().MyPID() << ".pos";
    std::ofstream file(str.str().c_str());
    for (std::size_t section = 0; section < sections.size(); ++section)
    {
      CORE::GEO::CUT::OUTPUT::GmshNewSection(file, sections[section], false);
      for (std::size_t entry = 0; entry < (plot_data_[section]).size(); ++entry)
      {
        CORE::GEO::CUT::OUTPUT::GmshCoordDump(
            file, (plot_data_[section])[entry].first, (plot_data_[section])[entry].second);
      }
      CORE::GEO::CUT::OUTPUT::GmshEndSection(file, false);
    }
    file.close();
  }
  ++counter;
  if (counter > 100) counter = 1;

  plot_data_.clear();
  plot_data_.resize(sections.size());
#endif

  std::vector<int> g_sum_gps(5);
  fluiddis_->Comm().SumAll(sum_gps_.data(), g_sum_gps.data(), 5);
  if (!fluiddis_->Comm().MyPID())
  {
    std::cout << "===| Summary Contact GPs |===" << std::endl;
    // 0 ... Contact, 1 ... Contact_NoContactNoFSI, 2 ... Contact_NoContactFSI, 3 ... FSI_NoContact,
    // 4 ... FSI_Contact
    std::cout << "Contact: " << g_sum_gps[0] << ", Contact_noC_noFSI: " << g_sum_gps[1]
              << ", Contact_noC_FSI: " << g_sum_gps[2] << "("
              << g_sum_gps[0] + g_sum_gps[1] + g_sum_gps[2] << "), FSI_noC: " << g_sum_gps[3]
              << ", FSI_C: " << g_sum_gps[4] << "(" << g_sum_gps[3] + g_sum_gps[4] << ") == (Sum: "
              << g_sum_gps[0] + g_sum_gps[1] + g_sum_gps[2] + g_sum_gps[3] + g_sum_gps[4] << ")"
              << " === (Fair Sum: " << g_sum_gps[0] + g_sum_gps[1] + g_sum_gps[3] << ")"
              << std::endl;
    std::cout << "===| ------------------- |===" << std::endl;
  }
  sum_gps_.clear();
  sum_gps_.resize(5);
}

void XFEM::XFluid_Contact_Comm::Gmsh_Write(CORE::LINALG::Matrix<3, 1> x, double val, int section)
{
#ifdef WRITE_GMSH
  plot_data_[section].push_back(std::pair<CORE::LINALG::Matrix<3, 1>, double>(x, val));
#endif
}

FOUR_C_NAMESPACE_CLOSE
