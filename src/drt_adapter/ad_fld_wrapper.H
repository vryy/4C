/*----------------------------------------------------------------------*/
/*! \file

\brief fluid wrapper

\level 1

\maintainer  Martin Kronbichler
*/
/*---------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_WRAPPER_H
#define ADAPTER_FLUID_WRAPPER_H

#include "ad_fld_fluid.H"

namespace ADAPTER
{
  /// Just wrap, do nothing new, meant to be derived from
  class FluidWrapper : public Fluid
  {
   public:
    explicit FluidWrapper(Teuchos::RCP<Fluid> fluid) : fluid_(fluid) {}

    virtual void Init()
    {
      fluid_->Init();
      return;
    }
    virtual Teuchos::RCP<const Epetra_Vector> InitialGuess() { return fluid_->InitialGuess(); }
    virtual Teuchos::RCP<const Epetra_Vector> RHS() { return fluid_->RHS(); }
    virtual Teuchos::RCP<const Epetra_Vector> TrueResidual() { return fluid_->TrueResidual(); }
    virtual Teuchos::RCP<const Epetra_Vector> Velnp() { return fluid_->Velnp(); }
    virtual Teuchos::RCP<const Epetra_Vector> Velaf() { return fluid_->Velaf(); }
    virtual Teuchos::RCP<const Epetra_Vector> Veln() { return fluid_->Veln(); }
    virtual Teuchos::RCP<const Epetra_Vector> Velnm()
    {
      dserror("not implemented");
      return Teuchos::null;
    };
    virtual Teuchos::RCP<const Epetra_Vector> Stepinc() { return fluid_->Stepinc(); }
    virtual Teuchos::RCP<const Epetra_Vector> Accnp() { return fluid_->Accnp(); };
    virtual Teuchos::RCP<const Epetra_Vector> Accn() { return fluid_->Accn(); };
    virtual Teuchos::RCP<const Epetra_Vector> Accnm() { return fluid_->Accnm(); };
    virtual Teuchos::RCP<const Epetra_Vector> Accam() { return fluid_->Accam(); }
    virtual Teuchos::RCP<const Epetra_Vector> Scaaf() { return fluid_->Scaaf(); };
    virtual Teuchos::RCP<const Epetra_Vector> Scaam()
    {
      dserror("not implemented");
      return Teuchos::null;
    };
    virtual Teuchos::RCP<const Epetra_Vector> Hist() { return fluid_->Hist(); }
    virtual Teuchos::RCP<const Epetra_Vector> GridVel() { return fluid_->GridVel(); }
    virtual Teuchos::RCP<const Epetra_Vector> GridVeln() { return fluid_->GridVeln(); }
    virtual Teuchos::RCP<const Epetra_Vector> Dispnp() { return fluid_->Dispnp(); }
    virtual Teuchos::RCP<const Epetra_Vector> Dispn() { return fluid_->Dispn(); }
    virtual Teuchos::RCP<const Epetra_Vector> ConvectiveVel() { return fluid_->ConvectiveVel(); }
    virtual Teuchos::RCP<const Epetra_Vector> FsVel()
    {
      dserror("not implemented");
      return Teuchos::null;
    };
    virtual Teuchos::RCP<Epetra_Vector> StdVeln()
    {
      dserror("not implemented");
      return Teuchos::null;
    };
    virtual Teuchos::RCP<Epetra_Vector> StdVelnp()
    {
      dserror("not implemented");
      return Teuchos::null;
    };
    virtual Teuchos::RCP<Epetra_Vector> StdVelaf()
    {
      dserror("not implemented");
      return Teuchos::null;
    };
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap() { return fluid_->DofRowMap(); }
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds)
    {
      return fluid_->DofRowMap(nds);
    };
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return fluid_->SystemMatrix(); }
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemSparseMatrix()
    {
      return fluid_->SystemSparseMatrix();
    }
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    {
      return fluid_->BlockSystemMatrix();
    }
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives()
    {
      return fluid_->ShapeDerivatives();
    }
    virtual const Teuchos::RCP<DRT::Discretization>& Discretization()
    {
      return fluid_->Discretization();
    }
    virtual Teuchos::RCP<const DRT::DofSet> DofSet() { return fluid_->DofSet(); }
    virtual Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return fluid_->GetDBCMapExtractor();
    }
    virtual void SetInitialFlowField(
        const INPAR::FLUID::InitialField initfield, const int startfuncno)
    {
      return fluid_->SetInitialFlowField(initfield, startfuncno);
    }
    virtual void SetInitialPorosityField(
        const INPAR::POROELAST::InitialField initfield, const int startfuncno)
    {
      return fluid_->SetInitialPorosityField(initfield, startfuncno);
    };
    virtual void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext)
    {
      return fluid_->ApplyExternalForces(fext);
    };
    virtual void AddContributionToExternalLoads(
        const Teuchos::RCP<const Epetra_Vector> contributing_vector)
    {
      return fluid_->AddContributionToExternalLoads(contributing_vector);
    };
    virtual void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
    {
      return fluid_->AddDirichCond(maptoadd);
    };
    virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
    {
      return fluid_->RemoveDirichCond(maptoremove);
    };
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> vel)
    {
      return fluid_->UpdateNewton(vel);
    };
    virtual void SetLomaIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
        const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
        Teuchos::RCP<DRT::Discretization> scatradis)
    {
      return fluid_->SetLomaIterScalarFields(scalaraf, scalaram, scalardtam, fsscalaraf,
          thermpressaf, thermpressam, thermpressdtaf, thermpressdtam, scatradis);
    }
    virtual void SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<DRT::Discretization> scatradis, int dofset = 0)
    {
      return fluid_->SetIterScalarFields(scalaraf, scalaram, scalardtam, scatradis, dofset);
    }
    virtual void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
        const double thermpressnp, Teuchos::RCP<const Epetra_Vector> scatraresidual,
        Teuchos::RCP<DRT::Discretization> scatradis, const int whichscalar = -1)
    {
      return fluid_->SetScalarFields(
          scalarnp, thermpressnp, scatraresidual, scatradis, whichscalar);
    }
    virtual Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager()
    {
      return fluid_->TurbulenceStatisticManager();
    }
    virtual Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter() { return fluid_->DynSmagFilter(); }
    virtual Teuchos::RCP<FLD::Vreman> Vreman() { return fluid_->Vreman(); }
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> velnp)
    {
      dserror("not implemented!");
      return;
    };
    //    virtual void TimeLoop()
    //    { return fluid_->TimeLoop(); }
    virtual void Integrate() { return fluid_->Integrate(); }
    virtual void PrepareTimeStep() { return fluid_->PrepareTimeStep(); }
    virtual void IncrementTimeAndStep()
    {
      dserror("not implemented!");
      return;
    }
    virtual void PrepareSolve() { fluid_->PrepareSolve(); }
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc)
    {
      return fluid_->Evaluate(stepinc);
    }
    virtual bool ConvergenceCheck(int itnum, int itmax, const double ittol)
    {
      dserror("not implemented!");
      return false;
    }
    virtual void IterUpdate(const Teuchos::RCP<const Epetra_Vector> increment)
    {
      dserror("not implemented!");
      return;
    }
    virtual void Update() { return fluid_->Update(); }
    virtual void StatisticsAndOutput() { return fluid_->StatisticsAndOutput(); }
    virtual void Output() { return fluid_->Output(); }
    virtual void StatisticsOutput()
    {
      dserror("not implemented!");
      return;
    }
    virtual const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter()
    {
      return fluid_->DiscWriter();
    }
    virtual const Teuchos::RCP<LINALG::MapExtractor> GetVelPressSplitter()
    {
      return fluid_->GetVelPressSplitter();
    }
    virtual void ReadRestart(int step) { return fluid_->ReadRestart(step); }
    virtual void SetRestart(const int step, const double time,
        Teuchos::RCP<const Epetra_Vector> readvelnp, Teuchos::RCP<const Epetra_Vector> readveln,
        Teuchos::RCP<const Epetra_Vector> readvelnm, Teuchos::RCP<const Epetra_Vector> readaccnp,
        Teuchos::RCP<const Epetra_Vector> readaccn)
    {
      dserror("not implemented!");
      return;
    }
    virtual double Time() const { return fluid_->Time(); }
    virtual int Step() const { return fluid_->Step(); }
    virtual double Dt() const { return fluid_->Dt(); }

    //! @name Write access to field solution variables at \f$t^{n+1}\f$
    //@{

    /// write access to extract velocities at \f$t^{n+1}\f$
    virtual Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() { return fluid_->WriteAccessVelnp(); }

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
    virtual void TimeStepAuxiliar(){};

    /*! Indicate norms of temporal discretization error
     *
     *  \author mayr.mt \date 12/2013
     */
    virtual void IndicateErrorNorms(
        double& err,       ///< L2-norm of temporal discretization error based on all DOFs
        double& errcond,   ///< L2-norm of temporal discretization error based on interface DOFs
        double& errother,  ///< L2-norm of temporal discretization error based on interior DOFs
        double& errinf,    ///< L-inf-norm of temporal discretization error based on all DOFs
        double&
            errinfcond,  ///< L-inf-norm of temporal discretization error based on interface DOFs
        double& errinfother  ///< L-inf-norm of temporal discretization error based on interior DOFs
    )
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    };

    //! Set fluid time step size that was computed outside
    virtual void SetDt(const double dtnew) { fluid_->SetDt(dtnew); }

    //! Reset last time step
    virtual void ResetStep() { fluid_->ResetStep(); }

    //! Reset time and step number of last time step, needed for time step size adaptivity an FSI
    virtual void ResetTime(const double dtold) { fluid_->ResetTime(dtold); }

    //! Set time and step
    virtual void SetTimeStep(const double time, const int step) { fluid_->SetTimeStep(time, step); }

    //@}

    virtual double EvalTime() const
    {
      dserror("not implemented!");
      return 0.0;
    }
    virtual void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
    {
      dserror("not implemented!");
      return;
    }
    virtual void Solve() { return fluid_->Solve(); }
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      return Teuchos::null;
    }
    virtual Teuchos::RCP<LINALG::Solver> LinearSolver() { return fluid_->LinearSolver(); }
    virtual void CalcIntermediateSolution() { return fluid_->CalcIntermediateSolution(); }
    virtual Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap()
    {
      return fluid_->InnerVelocityRowMap();
    }
    virtual Teuchos::RCP<const Epetra_Map> VelocityRowMap() { return fluid_->VelocityRowMap(); }
    virtual Teuchos::RCP<const Epetra_Map> PressureRowMap() { return fluid_->PressureRowMap(); }
    virtual void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm, const int nds_master = 0)
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }

    /// Use ResidualScaling() to convert the implemented fluid residual to an actual force with unit
    /// Newton [N]
    virtual double ResidualScaling() const { return fluid_->ResidualScaling(); }

    /// Velocity-displacement conversion at the fsi interface
    virtual double TimeScaling() const { return fluid_->TimeScaling(); }

    /// return time integration factor
    virtual double TimIntParam() const { return fluid_->TimIntParam(); }

    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }

    virtual Teuchos::RCP<FLD::UTILS::MapExtractor> const& FPSIInterface() const
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
      static Teuchos::RCP<FLD::UTILS::MapExtractor> ret = Teuchos::null;
      return ret;
    }
    virtual INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const
    {
      return fluid_->TimIntScheme();
    }
    virtual Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(
        Teuchos::RCP<const Epetra_Vector> velpres)
    {
      return fluid_->ExtractVelocityPart(velpres);
    }
    virtual Teuchos::RCP<const Epetra_Vector> ExtractPressurePart(
        Teuchos::RCP<const Epetra_Vector> velpres)
    {
      return fluid_->ExtractPressurePart(velpres);
    }
    virtual void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
    {
      return fluid_->ApplyInterfaceVelocities(ivel);
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp()
    {
      return fluid_->ExtractInterfaceVelnp();
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln()
    {
      return fluid_->ExtractInterfaceVeln();
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractFreeSurfaceVeln()
    {
      return fluid_->ExtractFreeSurfaceVeln();
    }
    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces()
    {
      return fluid_->ExtractInterfaceForces();
    }
    virtual void ApplyInitialMeshDisplacement(Teuchos::RCP<const Epetra_Vector> initfluiddisp)
    {
      fluid_->ApplyInitialMeshDisplacement(initfluiddisp);
    }
    virtual void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
    {
      return fluid_->ApplyMeshDisplacement(fluiddisp);
    }
    virtual void ApplyMeshDisplacementIncrement(Teuchos::RCP<const Epetra_Vector> dispstepinc)
    {
      return fluid_->ApplyMeshDisplacementIncrement(dispstepinc);
    }
    virtual void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
    {
      return fluid_->ApplyMeshVelocity(gridvel);
    }
    virtual void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
    {
      return fluid_->DisplacementToVelocity(fcx);
    }
    virtual void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
    {
      return fluid_->VelocityToDisplacement(fcx);
    }
    virtual void FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
    {
      return fluid_->FreeSurfDisplacementToVelocity(fcx);
    }
    virtual void FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
    {
      return fluid_->FreeSurfVelocityToDisplacement(fcx);
    }
    virtual int Itemax() const { return fluid_->Itemax(); }
    virtual void SetItemax(int itemax) { return fluid_->SetItemax(itemax); }
    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape()
    {
      return fluid_->IntegrateInterfaceShape();
    }
    virtual void UseBlockMatrix(bool splitmatrix) { return fluid_->UseBlockMatrix(splitmatrix); }
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest() { return fluid_->CreateFieldTest(); }
    virtual void Reset(bool completeReset = false, int numsteps = 1, int iter = -1)
    {
      return fluid_->Reset(completeReset, numsteps, iter);
    };
    virtual void SetFldGrDisp(Teuchos::RCP<Epetra_Vector> fluid_growth_disp)
    {
      return fluid_->SetFldGrDisp(fluid_growth_disp);
    }

    /// calculate error in comparison to analytical solution
    virtual void CalculateError()
    {
      dserror("Not implemented in the base class, may be overridden by a subclass.");
    }

    /// return physical type of fluid algorithm
    virtual INPAR::FLUID::PhysicalType PhysicalType() const { return fluid_->PhysicalType(); }

   protected:
    Teuchos::RCP<Fluid> fluid_;
  };
}  // namespace ADAPTER

#endif
