/*!----------------------------------------------------------------------
\file xfluid_outputservice.cpp
\brief Service class for XFluid(Fluid)-related output, mostly GMSH.

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/

#include "xfluid_outputservice.H"

#include "../drt_lib/drt_globalproblem.H"


FLD::XFluidOutputService::XFluidOutputService(
  Teuchos::ParameterList& params_xfem,
  const Teuchos::RCP<DRT::DiscretizationXFEM>& discret,
  const Teuchos::RCP<XFEM::ConditionManager>&  cond_manager
)
  :
  gmsh_active_           (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->IOParams(),"OUTPUT_GMSH")),
  gmsh_sol_out_          (gmsh_active_ && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_SOL_OUT")),
  gmsh_debug_out_        (gmsh_active_ && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT")),
  gmsh_debug_out_screen_ (gmsh_active_ && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DEBUG_OUT_SCREEN")),
  gmsh_EOS_out_          (gmsh_active_ && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_EOS_OUT")),
  gmsh_discret_out_      (gmsh_active_ && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_DISCRET_OUT")),
  gmsh_cut_out_          (gmsh_active_ && (bool)DRT::INPUT::IntegralValue<int>(params_xfem,"GMSH_CUT_OUT")),
  gmsh_step_diff_        (500),
  discret_(discret)
{};

void FLD::XFluidOutputService::GmshOutput(
  const std::string & filename_base,  ///< name for output file
  const int step,                  ///< step number
  const int count,                 ///< counter for iterations within a global time step
  Teuchos::RCP<Epetra_Vector> vel, ///< vector holding velocity and pressure dofs
  Teuchos::RCP<Epetra_Vector> acc, ///< vector holding accelerations
  Teuchos::RCP<Epetra_Vector> disp ///< vector holding displacements
)
{
  // Todo: fill!
}

/// Gmsh output function for elements without an GEO::CUT::ElementHandle
void FLD::XFluidOutputService::GmshOutputElement(
  DRT::Discretization & discret, ///< background fluid discretization
  std::ofstream & vel_f,         ///< output file stream for velocity
  std::ofstream & press_f,       ///< output file stream for pressure
  std::ofstream & acc_f,         ///< output file stream for acceleration
  DRT::Element * actele,         ///< element
  std::vector<int> & nds,        ///< vector holding the nodal dofsets
  Teuchos::RCP<const Epetra_Vector> vel, ///< vector holding velocity and pressure dofs
  Teuchos::RCP<const Epetra_Vector> acc  ///< vector holding acceleration
)
{
  // Todo: fill!
}

/// Gmsh output function for volumecells
void FLD::XFluidOutputService::GmshOutputVolumeCell(
  DRT::Discretization & discret,    ///< background fluid discretization
  std::ofstream & vel_f,            ///< output file stream for velocity
  std::ofstream & press_f,          ///< output file stream for pressure
  std::ofstream & acc_f,            ///< output file stream for acceleration
  DRT::Element * actele,            ///< element
  GEO::CUT::ElementHandle * e,      ///<elementhandle
  GEO::CUT::VolumeCell * vc,        ///< volumecell
  const std::vector<int> & nds,     ///< vector holding the nodal dofsets
  Teuchos::RCP<const Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
  Teuchos::RCP<const Epetra_Vector> acc   ///< vector holding acceleration
)
{
  // Todo: fill!
}

/// Gmsh output function for boundarycells
void FLD::XFluidOutputService::GmshOutputBoundaryCell(
  DRT::Discretization & discret,    ///< background fluid discretization
  DRT::Discretization & cutdiscret, ///< cutter surface discretization
  std::ofstream & bound_f,          ///< output file stream for boundary mesh
  GEO::CUT::VolumeCell * vc         ///< volumecell
)
{
  // Todo: fill!
}

/// print discretization to gmsh stream
void FLD::XFluidOutputService::PrintDiscretizationToStream(
  Teuchos::RCP<DRT::Discretization> dis,
  const std::string& disname,
  const bool elements,
  const bool elecol,
  const bool nodes,
  const bool nodecol,
  const bool faces,
  const bool facecol,
  std::ostream& s,
  std::map<int, LINALG::Matrix<3,1> >* curr_pos
)
{
  // Todo: fill!
}

