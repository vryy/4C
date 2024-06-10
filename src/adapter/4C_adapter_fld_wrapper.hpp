/*----------------------------------------------------------------------*/
/*! \file

\brief fluid wrapper

\level 1

*/
/*---------------------------------------------------------------------*/
#ifndef FOUR_C_ADAPTER_FLD_WRAPPER_HPP
#define FOUR_C_ADAPTER_FLD_WRAPPER_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Adapter
{
  /// Just wrap, do nothing new, meant to be derived from
  class FluidWrapper : public Fluid
  {
   public:
    explicit FluidWrapper(Teuchos::RCP<Fluid> fluid) : fluid_(fluid) {}

    void Init() override
    {
      fluid_->Init();
      return;
    }
    Teuchos::RCP<const Epetra_Vector> initial_guess() override { return fluid_->initial_guess(); }
    Teuchos::RCP<const Epetra_Vector> RHS() override { return fluid_->RHS(); }
    Teuchos::RCP<const Epetra_Vector> TrueResidual() override { return fluid_->TrueResidual(); }
    Teuchos::RCP<const Epetra_Vector> Velnp() override { return fluid_->Velnp(); }
    Teuchos::RCP<const Epetra_Vector> Velaf() override { return fluid_->Velaf(); }
    Teuchos::RCP<const Epetra_Vector> Veln() override { return fluid_->Veln(); }
    Teuchos::RCP<const Epetra_Vector> Velnm() override
    {
      FOUR_C_THROW("not implemented");
      return Teuchos::null;
    };
    Teuchos::RCP<const Epetra_Vector> Stepinc() override { return fluid_->Stepinc(); }
    Teuchos::RCP<const Epetra_Vector> Accnp() override { return fluid_->Accnp(); };
    Teuchos::RCP<const Epetra_Vector> Accn() override { return fluid_->Accn(); };
    Teuchos::RCP<const Epetra_Vector> Accnm() override { return fluid_->Accnm(); };
    Teuchos::RCP<const Epetra_Vector> Accam() override { return fluid_->Accam(); }
    Teuchos::RCP<const Epetra_Vector> Scaaf() override { return fluid_->Scaaf(); };
    Teuchos::RCP<const Epetra_Vector> Scaam() override
    {
      FOUR_C_THROW("not implemented");
      return Teuchos::null;
    };
    Teuchos::RCP<const Epetra_Vector> Hist() override { return fluid_->Hist(); }
    Teuchos::RCP<const Epetra_Vector> GridVel() override { return fluid_->GridVel(); }
    Teuchos::RCP<const Epetra_Vector> GridVeln() override { return fluid_->GridVeln(); }
    Teuchos::RCP<const Epetra_Vector> Dispnp() override { return fluid_->Dispnp(); }
    Teuchos::RCP<const Epetra_Vector> Dispn() override { return fluid_->Dispn(); }
    Teuchos::RCP<const Epetra_Vector> ConvectiveVel() override { return fluid_->ConvectiveVel(); }
    Teuchos::RCP<const Epetra_Vector> FsVel() override
    {
      FOUR_C_THROW("not implemented");
      return Teuchos::null;
    };
    Teuchos::RCP<Epetra_Vector> StdVeln() override
    {
      FOUR_C_THROW("not implemented");
      return Teuchos::null;
    };
    Teuchos::RCP<Epetra_Vector> StdVelnp() override
    {
      FOUR_C_THROW("not implemented");
      return Teuchos::null;
    };
    Teuchos::RCP<Epetra_Vector> StdVelaf() override
    {
      FOUR_C_THROW("not implemented");
      return Teuchos::null;
    };
    Teuchos::RCP<const Epetra_Map> dof_row_map() override { return fluid_->dof_row_map(); }
    Teuchos::RCP<const Epetra_Map> dof_row_map(unsigned nds) override
    {
      return fluid_->dof_row_map(nds);
    };
    Teuchos::RCP<Core::LinAlg::SparseMatrix> SystemMatrix() override
    {
      return fluid_->SystemMatrix();
    }
    Teuchos::RCP<Core::LinAlg::SparseMatrix> SystemSparseMatrix() override
    {
      return fluid_->SystemSparseMatrix();
    }
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> BlockSystemMatrix() override
    {
      return fluid_->BlockSystemMatrix();
    }
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> ShapeDerivatives() override
    {
      return fluid_->ShapeDerivatives();
    }
    const Teuchos::RCP<Core::FE::Discretization>& discretization() override
    {
      return fluid_->discretization();
    }
    Teuchos::RCP<const Core::DOFSets::DofSet> DofSet() override { return fluid_->DofSet(); }
    Teuchos::RCP<const Core::LinAlg::MapExtractor> GetDBCMapExtractor() override
    {
      return fluid_->GetDBCMapExtractor();
    }
    void SetInitialFlowField(
        const Inpar::FLUID::InitialField initfield, const int startfuncno) override
    {
      return fluid_->SetInitialFlowField(initfield, startfuncno);
    }
    void set_initial_porosity_field(
        const Inpar::PoroElast::InitialField initfield, const int startfuncno) override
    {
      return fluid_->set_initial_porosity_field(initfield, startfuncno);
    };
    void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext) override
    {
      return fluid_->ApplyExternalForces(fext);
    };
    void add_contribution_to_external_loads(
        const Teuchos::RCP<const Epetra_Vector> contributing_vector) override
    {
      return fluid_->add_contribution_to_external_loads(contributing_vector);
    };
    void add_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoadd) override
    {
      return fluid_->add_dirich_cond(maptoadd);
    };
    void remove_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoremove) override
    {
      return fluid_->remove_dirich_cond(maptoremove);
    };
    void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel) override
    {
      return fluid_->UpdateNewton(vel);
    };
    void set_loma_iter_scalar_fields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
        const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
        Teuchos::RCP<Core::FE::Discretization> scatradis) override
    {
      return fluid_->set_loma_iter_scalar_fields(scalaraf, scalaram, scalardtam, fsscalaraf,
          thermpressaf, thermpressam, thermpressdtaf, thermpressdtam, scatradis);
    }
    void SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<Core::FE::Discretization> scatradis, int dofset = 0) override
    {
      return fluid_->SetIterScalarFields(scalaraf, scalaram, scalardtam, scatradis, dofset);
    }
    void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp, const double thermpressnp,
        Teuchos::RCP<const Epetra_Vector> scatraresidual,
        Teuchos::RCP<Core::FE::Discretization> scatradis, const int whichscalar = -1) override
    {
      return fluid_->SetScalarFields(
          scalarnp, thermpressnp, scatraresidual, scatradis, whichscalar);
    }
    Teuchos::RCP<FLD::TurbulenceStatisticManager> turbulence_statistic_manager() override
    {
      return fluid_->turbulence_statistic_manager();
    }
    Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter() override { return fluid_->DynSmagFilter(); }
    Teuchos::RCP<FLD::Vreman> Vreman() override { return fluid_->Vreman(); }
    void set_velocity_field(Teuchos::RCP<const Epetra_Vector> velnp) override
    {
      FOUR_C_THROW("not implemented!");
      return;
    };
    //    virtual void TimeLoop()
    //    { return fluid_->TimeLoop(); }
    void Integrate() override { return fluid_->Integrate(); }
    void prepare_time_step() override { return fluid_->prepare_time_step(); }
    void increment_time_and_step() override
    {
      FOUR_C_THROW("not implemented!");
      return;
    }
    void PrepareSolve() override { fluid_->PrepareSolve(); }
    void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc) override
    {
      return fluid_->Evaluate(stepinc);
    }
    bool convergence_check(int itnum, int itmax, const double velrestol, const double velinctol,
        const double presrestol, const double presinctol) override
    {
      FOUR_C_THROW("not implemented!");
      return false;
    }
    void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment) override
    {
      FOUR_C_THROW("not implemented!");
      return;
    }
    void Update() override { return fluid_->Update(); }
    void StatisticsAndOutput() override { return fluid_->StatisticsAndOutput(); }
    void Output() override { return fluid_->Output(); }
    void StatisticsOutput() override
    {
      FOUR_C_THROW("not implemented!");
      return;
    }
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& DiscWriter() override
    {
      return fluid_->DiscWriter();
    }
    Teuchos::RCP<Core::LinAlg::MapExtractor> GetVelPressSplitter() override
    {
      return fluid_->GetVelPressSplitter();
    }
    void read_restart(int step) override { return fluid_->read_restart(step); }
    void SetRestart(const int step, const double time, Teuchos::RCP<const Epetra_Vector> readvelnp,
        Teuchos::RCP<const Epetra_Vector> readveln, Teuchos::RCP<const Epetra_Vector> readvelnm,
        Teuchos::RCP<const Epetra_Vector> readaccnp,
        Teuchos::RCP<const Epetra_Vector> readaccn) override
    {
      FOUR_C_THROW("not implemented!");
      return;
    }
    double Time() const override { return fluid_->Time(); }
    int Step() const override { return fluid_->Step(); }
    double Dt() const override { return fluid_->Dt(); }

    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    /// write access to extract velocities at \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() override { return fluid_->WriteAccessVelnp(); }

    //@}

    //! @name Time step size adaptivity in monolithic FSI
    //@{

    /*! Do one step with auxiliary time integration scheme
     *
     *  Do a single time step with the user given auxiliary time integration
     *  scheme. Result is stored in \p locerrvelnp_ and is used later to estimate
     *  the local discretization error of the marching time integration scheme.
     *
     *  \author mayr.mt \date 12/2013
     */
    void time_step_auxiliar() override{};

    /*! Indicate norms of temporal discretization error
     *
     *  \author mayr.mt \date 12/2013
     */
    void IndicateErrorNorms(
        double& err,       ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,   ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,  ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,    ///< L-inf-norm of temporal discretization error based on all DOFs
        double&
            errinfcond,  ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother  ///< L-inf-norm of temporal discretization error based on interior DOFs
        ) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    };

    //! Set fluid time step size that was computed outside
    void set_dt(const double dtnew) override { fluid_->set_dt(dtnew); }

    //! Reset last time step
    void reset_step() override { fluid_->reset_step(); }

    //! Reset time and step number of last time step, needed for time step size adaptivity an FSI
    void reset_time(const double dtold) override { fluid_->reset_time(dtold); }

    //! Set time and step
    void SetTimeStep(const double time, const int step) override
    {
      fluid_->SetTimeStep(time, step);
    }

    //@}

    double EvalTime() const override
    {
      FOUR_C_THROW("not implemented!");
      return 0.0;
    }
    void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph) override
    {
      FOUR_C_THROW("not implemented!");
      return;
    }
    void Solve() override { return fluid_->Solve(); }
    Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    Teuchos::RCP<Core::LinAlg::Solver> LinearSolver() override { return fluid_->LinearSolver(); }
    void calc_intermediate_solution() override { return fluid_->calc_intermediate_solution(); }
    Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap() override
    {
      return fluid_->InnerVelocityRowMap();
    }
    Teuchos::RCP<const Epetra_Map> VelocityRowMap() override { return fluid_->VelocityRowMap(); }
    Teuchos::RCP<const Epetra_Map> PressureRowMap() override { return fluid_->PressureRowMap(); }
    void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm, const int nds_master = 0) override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }

    /// Use residual_scaling() to convert the implemented fluid residual to an actual force with
    /// unit Newton [N]
    double residual_scaling() const override { return fluid_->residual_scaling(); }

    /// Velocity-displacement conversion at the fsi interface
    double TimeScaling() const override { return fluid_->TimeScaling(); }

    /// return time integration factor
    double TimIntParam() const override { return fluid_->TimIntParam(); }

    Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }

    Teuchos::RCP<FLD::UTILS::MapExtractor> const& FPSIInterface() const override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }
    Inpar::FLUID::TimeIntegrationScheme TimIntScheme() const override
    {
      return fluid_->TimIntScheme();
    }
    Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(
        Teuchos::RCP<const Epetra_Vector> velpres) override
    {
      return fluid_->ExtractVelocityPart(velpres);
    }
    Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(
        Teuchos::RCP<const Epetra_Vector> velpres) override
    {
      return fluid_->ExtractPressurePart(velpres);
    }
    void apply_interface_velocities(Teuchos::RCP<Epetra_Vector> ivel) override
    {
      return fluid_->apply_interface_velocities(ivel);
    }
    Teuchos::RCP<Epetra_Vector> extract_interface_velnp() override
    {
      return fluid_->extract_interface_velnp();
    }
    Teuchos::RCP<Epetra_Vector> extract_interface_veln() override
    {
      return fluid_->extract_interface_veln();
    }
    Teuchos::RCP<Epetra_Vector> extract_free_surface_veln() override
    {
      return fluid_->extract_free_surface_veln();
    }
    Teuchos::RCP<Epetra_Vector> extract_interface_forces() override
    {
      return fluid_->extract_interface_forces();
    }
    /// Apply initial mesh displacement
    void apply_initial_mesh_displacement(Teuchos::RCP<const Epetra_Vector> initfluiddisp) override
    {
      fluid_->apply_initial_mesh_displacement(initfluiddisp);
    }
    void apply_mesh_displacement(Teuchos::RCP<const Epetra_Vector> fluiddisp) override
    {
      return fluid_->apply_mesh_displacement(fluiddisp);
    }
    void apply_mesh_displacement_increment(Teuchos::RCP<const Epetra_Vector> dispstepinc) override
    {
      return fluid_->apply_mesh_displacement_increment(dispstepinc);
    }
    void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel) override
    {
      return fluid_->ApplyMeshVelocity(gridvel);
    }
    void displacement_to_velocity(Teuchos::RCP<Epetra_Vector> fcx) override
    {
      return fluid_->displacement_to_velocity(fcx);
    }
    void velocity_to_displacement(Teuchos::RCP<Epetra_Vector> fcx) override
    {
      return fluid_->velocity_to_displacement(fcx);
    }
    void free_surf_displacement_to_velocity(Teuchos::RCP<Epetra_Vector> fcx) override
    {
      return fluid_->free_surf_displacement_to_velocity(fcx);
    }
    void free_surf_velocity_to_displacement(Teuchos::RCP<Epetra_Vector> fcx) override
    {
      return fluid_->free_surf_velocity_to_displacement(fcx);
    }
    int Itemax() const override { return fluid_->Itemax(); }
    void SetItemax(int itemax) override { return fluid_->SetItemax(itemax); }
    Teuchos::RCP<Epetra_Vector> integrate_interface_shape() override
    {
      return fluid_->integrate_interface_shape();
    }
    void use_block_matrix(bool splitmatrix) override
    {
      return fluid_->use_block_matrix(splitmatrix);
    }
    Teuchos::RCP<Core::UTILS::ResultTest> CreateFieldTest() override
    {
      return fluid_->CreateFieldTest();
    }
    void Reset(bool completeReset = false, int numsteps = 1, int iter = -1) override
    {
      return fluid_->Reset(completeReset, numsteps, iter);
    };
    void SetFldGrDisp(Teuchos::RCP<Epetra_Vector> fluid_growth_disp) override
    {
      return fluid_->SetFldGrDisp(fluid_growth_disp);
    }

    /// calculate error in comparison to analytical solution
    void CalculateError() override
    {
      FOUR_C_THROW("Not implemented in the base class, may be overridden by a subclass.");
    }

    /// return physical type of fluid algorithm
    Inpar::FLUID::PhysicalType PhysicalType() const override { return fluid_->PhysicalType(); }

   protected:
    Teuchos::RCP<Fluid> fluid_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
