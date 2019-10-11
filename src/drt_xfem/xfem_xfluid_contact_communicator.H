/*----------------------------------------------------------------------*/
/*! \file

\brief communicates between xfluid and NIT contact ... for XFSCI and XFPSCI(soon)

\level 3

\maintainer  Christoph Ager
*/
/*----------------------------------------------------------------------*/


#ifndef XFEM_XFLUID_CONTACT_COMM_H_
#define XFEM_XFLUID_CONTACT_COMM_H_


#include <vector>
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include "../drt_inpar/inpar_xfem.H"

//#define WRITE_GMSH

namespace MORTAR
{
  class MortarElement;
}
namespace CONTACT
{
  class CoElement;
  class CoNitscheStrategyFsi;
  class CoNitscheStrategyFpi;
  class CoNitscheStrategy;
}  // namespace CONTACT

namespace DRT
{
  namespace ELEMENTS
  {
    class StructuralSurface;
  }
  class Element;
}  // namespace DRT
namespace GEO
{
  namespace CUT
  {
    class SideHandle;
    class VolumeCell;
    class Facet;
    class Element;
    class ElementHandle;
    class Side;
  }  // namespace CUT
  class CutWizard;
}  // namespace GEO

namespace XFEM
{
  class ConditionManager;
  class MeshCouplingFSI;
  class MeshCoupling;
  class MeshCouplingFPI;

  class XFluid_Contact_Comm
  {
   public:
    //! constructor
    explicit XFluid_Contact_Comm(CONTACT::CoNitscheStrategy& contact_strategy)
        : fluid_init_(false),
          ele_ptrs_already_setup_(false),
          cutwizard_(Teuchos::null),
          fluiddis_(Teuchos::null),
          condition_manager_(Teuchos::null),
          mc_(std::vector<Teuchos::RCP<XFEM::MeshCoupling>>()),
          mcfpi_ps_pf_(Teuchos::null),
          mcidx_(0),
          isporo_(false),
          visc_stab_trace_estimate_(INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk),
          visc_stab_hk_(INPAR::XFEM::ViscStab_hk_ele_vol_div_by_max_ele_surf),
          nit_stab_gamma_(-1),
          is_pseudo_2D_(false),
          mass_conservation_scaling_(INPAR::XFEM::MassConservationScaling_only_visc),
          mass_conservation_combination_(INPAR::XFEM::MassConservationCombination_sum),
          dt_(-1),
          theta_(-1),
          parallel_(false),
          min_surf_id_(-1),
          min_mortar_id_(-1),
          soSurfId_to_mortar_ele_(std::vector<CONTACT::CoElement*>()),
          mortarId_to_soSurf_ele_(std::vector<DRT::ELEMENTS::StructuralSurface*>()),
          mortarId_to_somc_(std::vector<int>()),
          mortarId_to_sosid_(std::vector<int>()),
          extrapolate_to_zero_(false),
          my_sele_ids_(std::set<int>()),
          contact_ele_rowmap_fluidownerbased_(Teuchos::null),
          contact_strategy_(contact_strategy),
          contact_strategy_fsi_(NULL),
          contact_strategy_fpi_(NULL)
    {
    }

    //! destructor
    virtual ~XFluid_Contact_Comm() { return; }

    /// Initialize overall Fluid State (includes the Cut intersection information)
    void InitializeFluidState(Teuchos::RCP<GEO::CutWizard> cutwizard,
        Teuchos::RCP<DRT::Discretization> fluiddis,
        Teuchos::RCP<XFEM::ConditionManager> condition_manager,
        Teuchos::RCP<Teuchos::ParameterList> fluidparams);

    /// Reset overall Fluid State
    void ResetFluidState()
    {
      fluid_init_ = false;
      cutwizard_ = Teuchos::null;
      fluiddis_ = Teuchos::null;
    }

    /// Get the FSI traction called from contact gausspoint
    double Get_FSI_Traction(MORTAR::MortarElement* ele,  // Mortar Element
        const LINALG::Matrix<3, 1>& xsi_parent,          // local coord in the parent element
        const LINALG::Matrix<2, 1>& xsi_boundary,        // local coord on the boundary element
        const LINALG::Matrix<3, 1>& normal,              // normal for projection
        bool& FSI_integrated,
        bool& gp_on_this_proc,  // for serial run
        double* poropressure = NULL);

    /// Get the FSI traction called from contact gausspoint
    double Get_FSI_Traction(const MORTAR::MortarElement* ele,
        const LINALG::Matrix<2, 1>& xsi_parent, const LINALG::Matrix<1, 1>& xsi_boundary,
        const LINALG::Matrix<2, 1>& normal, bool& FSI_integrated,
        bool& gp_on_this_proc,  // for serial run
        double* poropressure = NULL)
    {
      dserror("no 2D xfsi with contact");
      return -1.0;
    }

    /// Get_Contact_State at gausspoint called from XFSI: return true-->evaluate FSI, return false
    /// -->evaluate NIT-Contact
    bool Get_Contact_State(int sid,                           // Solid Surface Element
        std::string mcname, const LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
        const double& full_fsi_traction,                      // stressfluid + penalty ...
        double& gap);

    /// Is this Structural surface registered in the Xfluid Contact Communicator
    bool IsRegisteredSurface(const int soSurfId)
    {
      return (soSurfId >= min_surf_id_ &&
              soSurfId < ((int)soSurfId_to_mortar_ele_.size() + min_surf_id_));
    }

    /// Get the contact element for this solid surface id
    CONTACT::CoElement* GetContactEle(const int soSurfId)
    {
      return soSurfId_to_mortar_ele_.at(soSurfId - min_surf_id_);
    }
    /// Get the solid surface element for the contact element id
    DRT::ELEMENTS::StructuralSurface* GetSurfEle(const int mortarId)
    {
      return mortarId_to_soSurf_ele_.at(mortarId - min_mortar_id_);
    }

    /// Get the mesh coupling id for the contact element id
    int GetSurfMc(const int mortarId) { return mortarId_to_somc_.at(mortarId - min_mortar_id_); }

    /// Get the solid surface element if for the contact element id
    int GetSurfSid(const int mortarId) { return mortarId_to_sosid_.at(mortarId - min_mortar_id_); }

    /// Setup Interface element connection vectors based on points
    void SetupSurfElePtrs(DRT::Discretization& contact_interface_dis);

    /// Get element size of background mesh
    double Get_h();

    /// Register Evaluation Processor rank for specific solid surface (is the fluid proc)
    void RegisterSideProc(int sid);

    /// Get the CUT integration points for this contact element (id)
    void GetCutSideIntegrationPoints(
        int sid, Epetra_SerialDenseMatrix& coords, std::vector<double>& weights, int& npg);

    /// Finalize Map of interface element owners
    void FillComplete_SeleMap();

    /// Rowmap of contact elements based on the fluid element owner
    Teuchos::RCP<Epetra_Map>& Get_ContactEleRowMap_FOwnerbased()
    {
      return contact_ele_rowmap_fluidownerbased_;
    }

    /// PrepareTimeStep
    void PrepareTimeStep();

    /// PrepareIterationStep
    void PrepareIterationStep();

    /// Register contact element for to use CUT integration points
    void RegisterContactElementforHigherIntegration(int cid)
    {
      higher_contact_elements_.insert(cid);
    }

    /// Does this contact element use CUT integration points?
    bool HigherIntegrationforContactElement(int cid)
    {
      return (higher_contact_elements_comm_.find(cid) != higher_contact_elements_comm_.end());
    }

    /// Initialize Gmsh files
    void Create_New_Gmsh_files();

    /// Write Gmsh files
    void Gmsh_Write(LINALG::Matrix<3, 1> x, double val, int section);

    /// Increment gausspoint counter
    void Inc_GP(int state) { ++sum_gps_[state]; }

    //! get distance when transition between FPSI and PSCI is started
    double Get_fpi_pcontact_exchange_dist();

    //! ration of gap/(POROCONTACTFPSI_HFRACTION*h) when full PSCI is starte
    double Get_fpi_pcontact_fullfraction();

   private:
    //! The the contact state at local coord of CoElement cele and compare to the fsi_traction,
    //! return true if contact is evaluated, reture false if FSI is evaluated
    bool CheckNitscheContactState(CONTACT::CoElement* cele,
        const LINALG::Matrix<2, 1>& xsi,  // local coord on the ele element
        const double& full_fsi_traction,  // stressfluid + penalty
        double& gap                       // gap
    );

    /// Get the fluid states at specific selexi
    void Get_States(const int fluidele_id, const std::vector<int>& fluid_nds,
        const DRT::ELEMENTS::StructuralSurface* sele, const LINALG::Matrix<2, 1>& selexsi,
        const LINALG::Matrix<3, 1>& x, DRT::Element*& fluidele, Epetra_SerialDenseMatrix& ele_xyze,
        std::vector<double>& velpres, std::vector<double>& disp, std::vector<double>& ivel,
        double& pres_m, LINALG::Matrix<3, 1>& vel_m, LINALG::Matrix<3, 1>& vel_s,
        LINALG::Matrix<3, 3>& vderxy_m, LINALG::Matrix<3, 1>& velpf_s);

    /// Get the Nitsche penalty parameter
    void Get_Penalty_Param(DRT::Element* fluidele, GEO::CUT::VolumeCell* volumecell,
        Epetra_SerialDenseMatrix& ele_xyze, const LINALG::Matrix<3, 1>& elenormal,
        double& penalty_fac, const LINALG::Matrix<3, 1>& vel_m);

    /// Get the Nitsche penalty parameter
    void Get_Penalty_Param(DRT::ELEMENTS::StructuralSurface* sele, double& penalty_fac);

    /// Get the volumecell for local coord xsi on sele
    bool GetVolumecell(DRT::ELEMENTS::StructuralSurface*& sele, LINALG::Matrix<2, 1>& xsi,
        GEO::CUT::SideHandle*& sidehandle, std::vector<int>& nds, int& eleid,
        GEO::CUT::VolumeCell*& volumecell, LINALG::Matrix<3, 1>& elenormal, LINALG::Matrix<3, 1>& x,
        bool& FSI_integrated, double& distance);

    /// Evaluate the distance of x the boundary of a side
    double DistancetoSide(
        LINALG::Matrix<3, 1>& x, GEO::CUT::Side* side, LINALG::Matrix<3, 1>& closest_x);

    /// Find the next physical interface side to x
    GEO::CUT::Side* FindnextPhysicalSide(LINALG::Matrix<3, 1>& x, GEO::CUT::Side* initSide,
        GEO::CUT::SideHandle*& sidehandle, LINALG::Matrix<2, 1>& newxsi, double& distance);

    /// Get list of potentiall next physical sides
    void Update_physical_sides(GEO::CUT::Side* side, std::set<GEO::CUT::Side*>& performed_sides,
        std::set<GEO::CUT::Side*>& physical_sides);

    /// Get neighboring sides
    std::vector<GEO::CUT::Side*> GetNewNeighboringSides(
        GEO::CUT::Side* side, std::set<GEO::CUT::Side*>& performed_sides);

    /// Get next element
    GEO::CUT::Element* GetNextElement(
        GEO::CUT::Element* ele, std::set<GEO::CUT::Element*>& performed_elements, int& lastid);

    /// access to contact/meshtying bridge
    CONTACT::CoNitscheStrategy& GetContactStrategy() { return contact_strategy_; }

    /// fluid state members initialized
    bool fluid_init_;
    /// Surface element pointers setup
    bool ele_ptrs_already_setup_;
    /// The XFluid CutWizard
    Teuchos::RCP<GEO::CutWizard> cutwizard_;
    /// The Background Fluid Discretization
    Teuchos::RCP<DRT::Discretization> fluiddis_;
    /// The XFEM Condition Manager
    Teuchos::RCP<XFEM::ConditionManager> condition_manager_;
    /// A list of all mesh coupling objects
    std::vector<Teuchos::RCP<XFEM::MeshCoupling>> mc_;
    /// In case of poro, the fluid mesh coupling object
    Teuchos::RCP<XFEM::MeshCouplingFPI> mcfpi_ps_pf_;
    /// Mesh coupling index
    int mcidx_;
    /// Is a poro problem
    bool isporo_;

    /// Viscous trace estimate for FSI-Nit-Pen
    INPAR::XFEM::ViscStab_TraceEstimate visc_stab_trace_estimate_;
    /// h-definition for FSI-Nit-Pen
    INPAR::XFEM::ViscStab_hk visc_stab_hk_;
    /// reference penalty parameter for FSI-Nit-Pen
    double nit_stab_gamma_;
    /// pseudo 2D flag for 2D simulation with one element in z-direction
    bool is_pseudo_2D_;
    /// mass conservation scaline on FSI-Nit-Pen
    INPAR::XFEM::MassConservationScaling mass_conservation_scaling_;
    /// How to combine the contribution on FSI-Nit-Pen
    INPAR::XFEM::MassConservationCombination mass_conservation_combination_;
    /// timestep
    double dt_;
    /// theta factor of OST-scheme
    double theta_;

    /// parallel computation (NumProc > 1), there are some check we can only do in serial
    bool parallel_;

    /// Min Structural Surface Id
    int min_surf_id_;
    /// Min Mortar Element Id
    int min_mortar_id_;
    /// Vector for translation of Structural Surface Id to Contact Element
    std::vector<CONTACT::CoElement*> soSurfId_to_mortar_ele_;
    /// Vector for translation of Mortar Element Id to Structural Surface
    std::vector<DRT::ELEMENTS::StructuralSurface*> mortarId_to_soSurf_ele_;
    /// Vector for translation of Mortar Element Id to Mesh Coupling Object Id
    std::vector<int> mortarId_to_somc_;
    /// Vector for translation of Mortar Element Id to Structural Surface Id
    std::vector<int> mortarId_to_sosid_;

    /// Fluid traction in extrapolation zone goes to zero
    bool extrapolate_to_zero_;

    // all sele which have a row fluid-element on this proc
    std::set<int> my_sele_ids_;
    /// contact ele romap - based on the background fluid element owners
    Teuchos::RCP<Epetra_Map> contact_ele_rowmap_fluidownerbased_;

    /// The Contact Strategy
    CONTACT::CoNitscheStrategy& contact_strategy_;

    /// The Contact Strategy casted to fsi
    CONTACT::CoNitscheStrategyFsi* contact_strategy_fsi_;

    /// The Contact Strategy casted to fpi
    CONTACT::CoNitscheStrategyFpi* contact_strategy_fpi_;

    /// Contact Elements with increased number of GPs
    std::set<int> higher_contact_elements_;
    /// Contact Elements with increased number of GPs synchronized
    std::set<int> higher_contact_elements_comm_;

    /// For Gmsh Output
    std::vector<std::vector<std::pair<LINALG::Matrix<3, 1>, double>>> plot_data_;

    /// Summarized Contact gps
    /// 0 ... Contact, 1 ... Contact_NoContactNoFSI, 2 ... Contact_NoContactFSI, 3 ...
    /// FSI_NoContact, 4 ... FSI_Contact
    std::vector<int> sum_gps_;

    /// store the last evaluted set of physical sides for solid side with id key
    std::pair<int, std::set<GEO::CUT::Side*>> last_physical_sides_;

    /// last computed element h measure with key fluidele id
    std::pair<int, double> last_ele_h_;
  };  // class XFluid_Contact_Comm
}  // namespace XFEM

#endif /*XFEM_XFLUID_CONTACT_COMM_H_*/
