/*!----------------------------------------------------------------------
\file xfluid_outputservice.H
\brief Service class for XFluid(Fluid)-related output, mostly GMSH.

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/

#ifndef XFLUID_OUTPUTSERVICE_H_
#define XFLUID_OUTPUTSERVICE_H_

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_lib/drt_discret.H"

namespace DRT
{
  class Discretization;
  class DiscretizationXFEM;
}
namespace GEO
{
  class CutWizard;

  namespace CUT
  {
    class ElementHandle;
    class VolumeCell;
  }
}

namespace XFEM
{
  class ConditionManager;
}

namespace FLD
{

/*!
 * \brief Class handles output of gmsh data
 */
class XFluidOutputService
{

public:

  XFluidOutputService(
    Teuchos::ParameterList& params_xfem,
    const Teuchos::RCP<DRT::DiscretizationXFEM>& discret,
    const Teuchos::RCP<XFEM::ConditionManager>&  cond_manager
  );

  /// Write gmsh output depending on the flags
  void GmshOutput( const std::string & filename_base,               ///< name for output file
                   const int step,                                  ///< step number
                   const int count,                                 ///< counter for iterations within a global time step
                   Teuchos::RCP<Epetra_Vector> vel,                 ///< vector holding velocity and pressure dofs
                   Teuchos::RCP<Epetra_Vector> acc  = Teuchos::null,///< vector holding accelerations
                   Teuchos::RCP<Epetra_Vector> disp = Teuchos::null ///< vector holding displacements
                   );

private:

  /// Gmsh output function for elements without an GEO::CUT::ElementHandle
  void GmshOutputElement(
      DRT::Discretization & discret, ///< background fluid discretization
      std::ofstream & vel_f,         ///< output file stream for velocity
      std::ofstream & press_f,       ///< output file stream for pressure
      std::ofstream & acc_f,         ///< output file stream for acceleration
      DRT::Element * actele,         ///< element
      std::vector<int> & nds,        ///< vector holding the nodal dofsets
      Teuchos::RCP<const Epetra_Vector> vel, ///< vector holding velocity and pressure dofs
      Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null  ///< vector holding acceleration
  );

  /// Gmsh output function for volumecells
  void GmshOutputVolumeCell(
      DRT::Discretization & discret,    ///< background fluid discretization
      std::ofstream & vel_f,            ///< output file stream for velocity
      std::ofstream & press_f,          ///< output file stream for pressure
      std::ofstream & acc_f,            ///< output file stream for acceleration
      DRT::Element * actele,            ///< element
      GEO::CUT::ElementHandle * e,      ///<elementhandle
      GEO::CUT::VolumeCell * vc,        ///< volumecell
      const std::vector<int> & nds,     ///< vector holding the nodal dofsets
      Teuchos::RCP<const Epetra_Vector> vel,  ///< vector holding velocity and pressure dofs
      Teuchos::RCP<const Epetra_Vector> acc = Teuchos::null ///< vector holding acceleration
  );

  /// Gmsh output function for boundarycells
  void GmshOutputBoundaryCell(
      DRT::Discretization & discret,    ///< background fluid discretization
      DRT::Discretization & cutdiscret, ///< cutter surface discretization
      std::ofstream & bound_f,          ///< output file stream for boundary mesh
      GEO::CUT::VolumeCell * vc         ///< volumecell
  );

  /// print discretization to gmsh stream
  void PrintDiscretizationToStream(
       Teuchos::RCP<DRT::Discretization> dis,
       const std::string& disname,
       const bool elements,
       const bool elecol,
       const bool nodes,
       const bool nodecol,
       const bool faces,
       const bool facecol,
       std::ostream& s,
       std::map<int, LINALG::Matrix<3,1> >* curr_pos = NULL
   );


  /// global flag, that decides, whether gmsh output is globally activated
  const bool gmsh_active_;

  //! @name flags for detailed gmsh output
  const bool gmsh_sol_out_;          ///< Gmsh solution output for each timestep
  const bool gmsh_debug_out_;
  const bool gmsh_debug_out_screen_; ///< print information about output to screen
  const bool gmsh_EOS_out_;          ///< output for edge-oriented stabilization and ghost-penalty stabilization
  const bool gmsh_discret_out_;
  const bool gmsh_cut_out_;
  const int  gmsh_step_diff_;
  //@}

  //! @name Managed discretizations
  //@{
  const Teuchos::RCP<DRT::DiscretizationXFEM> discret_;
  //@}


};

}

#endif
