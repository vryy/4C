/*----------------------------------------------------------------------*/
/*!
\file xfem_coupling_fpi_mesh.H

\brief manages mesh based coupling of fluid and porous media
xfluid class and the cut-library

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_COUPLING_FPI_MESH_H_
#define XFEM_COUPLING_FPI_MESH_H_

#include "xfem_coupling_mesh.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_gausspoints.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../linalg/linalg_utils.H"

namespace GEO
{
  namespace CUT
  {
    class BoundaryCell;
  }
}  // namespace GEO

namespace XFEM
{
  /*!
  specialized class for coupling of fluid with an porous media (Darcy Flow)
    \author ager
    \date 06/16
   */
  class MeshCouplingFPI : public MeshVolCoupling
  {
   public:
    enum coupled_field
    {
      ps_ps,
      ps_pf,
      pf_ps,
      pf_pf
    };

    //! constructor
    explicit MeshCouplingFPI(
        Teuchos::RCP<DRT::Discretization>& bg_dis,  ///< background discretization
        const std::string& cond_name,  ///< name of the condition, by which the derived cutter
                                       ///< discretization is identified
        Teuchos::RCP<DRT::Discretization>&
            cond_dis,           ///< discretization from which cutter discretization can be derived
        const int coupling_id,  ///< id of composite of coupling conditions
        const double time,      ///< time
        const int step,         ///< time step
        MeshCouplingFPI::coupled_field field  ///< which field is coupled to the fluid
    );

    //! destructor
    virtual ~MeshCouplingFPI() { return; }

    //! cutter dis should be loaded into the cut?
    virtual bool CutGeometry() { return (coupled_field_ == MeshCouplingFPI::ps_ps); }

    bool isBJ() const { return full_BJ_; }

    void SetFullState(
        Teuchos::RCP<const Epetra_Vector> dispnp, Teuchos::RCP<const Epetra_Vector> pres)
    {
      fulldispnp_ = dispnp;
      fullpres_ = pres;
    }

    void InitializeStrucPresMap(
        Teuchos::RCP<const Epetra_Map> pfmap, Teuchos::RCP<const Epetra_Map> psmap)
    {
      // We need to identify cutter dis dofs and pressure dofs on all processors for the whole
      // cutter_dis, as long as we don't have another ghosting strategy for the cutter_dis ...

      if (pfmap->NumMyElements() != psmap->NumMyElements())
        dserror(
            "InitializeStrucPresMap: (pfmap->NumGlobalElements() != psmap->NumGlobalElements())!");

      Teuchos::RCP<Epetra_Map> fullpfmap = LINALG::AllreduceEMap(*pfmap);
      Teuchos::RCP<Epetra_Map> fullpsmap = LINALG::AllreduceEMap(*psmap);

      if (fullpfmap->NumMyElements() != fullpsmap->NumMyElements())
        dserror(
            "InitializeStrucPresMap: (fullpfmap->NumGlobalElements() != "
            "fullpsmap->NumGlobalElements())!");

      for (int lid = 0; lid < fullpfmap->NumMyElements(); ++lid)
        lm_struct_x_lm_pres_[fullpsmap->GID(lid)] =
            fullpfmap->GID(lid) + 1;  // z-component of structure --> pressure
    }

    /*!
     Return prescribed velocities and traction vectors for a GNBC boundary condition.
     Also returns the projection matrix (to the plane of the surface) needed for the GNBC condition.
     */
    template <DRT::Element::DiscretizationType DISTYPE, class T1, class M3>
    void EvaluateCouplingConditions(T1& projection_matrix,  ///< Projection matrix
        M3& normal                                          ///< surface normal of cut element
    )
    {
      EvalProjectionMatrix<DISTYPE>(projection_matrix, normal);
      return;
    };

    // finalize the interface true residual vector
    virtual void CompleteStateVectors();

    virtual void ZeroStateVectors_FPI();

    virtual void GmshOutput(const std::string& filename_base, const int step,
        const int gmsh_step_diff, const bool gmsh_debug_out_screen);

    virtual void GmshOutputDiscretization(std::ostream& gmshfilecontent);

    virtual void LiftDrag(const int step, const double time) const;

    virtual void ReadRestart(const int step);

    // interface foces
    Teuchos::RCP<Epetra_Vector> ITrueResidual() { return itrueresidual_; }

    // for assembly of fluid interface forces
    Teuchos::RCP<Epetra_Vector> IForcecol() { return iforcecol_; }

   private:
    void Output(const int step, const double time, const bool write_restart_data);

   protected:
    //! Initializes configurationmap
    virtual void SetupConfigurationMap();

    //! set the name of the coupling object based on the field coupling
    virtual void SetCouplingName();

    //! Updates configurationmap for specific Gausspoint
    virtual void UpdateConfigurationMap_GP(double& kappa_m,  //< fluid sided weighting
        double& visc_m,                                      //< master sided dynamic viscosity
        double& visc_s,                                      //< slave sided dynamic viscosity
        double& density_m,                                   //< master sided density
        double& visc_stab_tang,         //< viscous tangential NIT Penalty scaling
        double& full_stab,              //< full NIT Penalty scaling
        const LINALG::Matrix<3, 1>& x,  //< Position x in global coordinates
        const DRT::Condition* cond,     //< Condition
        DRT::Element* ele,              //< Element
        DRT::Element* bele,             //< Boundary Element
        double* funct,  //< local shape function for Gauss Point (from fluid element)
        double* derxy,  //< local derivatives of shape function for Gauss Point (from fluid element)
        LINALG::Matrix<3, 1>& rst_slave,  //< local coord of gp on slave boundary element
        LINALG::Matrix<3, 1>& normal,     //< normal at gp
        LINALG::Matrix<3, 1>& vel_m       //< master velocity at gp
    );

    virtual void InitStateVectors();

    virtual void DoConditionSpecificSetup();

    virtual bool HasMovingInterface() { return true; }

    virtual void SetConditionSpecificParameters();

    //! Caluculate the Porosity for J,porosity pair on this FaceElement
    double CalctrPermeability(DRT::Element* ele, double& porosity, double& J);

    //! Caluculate the Porosity for this FaceElement Gausspoint
    double CalcPorosity(DRT::Element* ele, LINALG::Matrix<3, 1>& rst_slave, double& J);

    //! Compute Jacobian and extract PoroFluidPressure this FaceElement Gausspoint
    double ComputeJacobianandPressure(
        DRT::Element* ele, LINALG::Matrix<3, 1>& rst_slave, double& pres);

    //------------------------------- vectors -----------------------------
    //! @name cutter-dis state vectors
    Teuchos::RCP<Epetra_Vector>
        itrueresidual_;  //! interface forces acting on the structural surface (= -iforcenp)
    Teuchos::RCP<Epetra_Vector>
        iforcecol_;  //! interface forces acting on the fluid surface (column vector assembly)
    //@}

    //---------------------------------configuration--------------------------------
    //! @name type of poro field coupled to xfluid handled by this mesh coupling object!
    MeshCouplingFPI::coupled_field coupled_field_;

    //! Full BJ Variant or BJSaffmann?
    bool full_BJ_;
    bool Sub_tang_;

    Teuchos::RCP<const Epetra_Vector> fulldispnp_;
    Teuchos::RCP<const Epetra_Vector> fullpres_;

    //! map from structural x dof to pres dof of a node!
    std::map<int, int> lm_struct_x_lm_pres_;

    double BJ_coeff_;

    //@}
  };

}  // namespace XFEM

#endif /*XFEM_COUPLING_FPI_MESH_H_*/
