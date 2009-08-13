/*
 * adapter_fluid_projection.cpp
 *
 *  Created on: Jun 16, 2009
 *      Author: wiesner
 */

#ifdef CCADISCRET

#include "adapter_fluid_projection.H"
#include "../drt_lib/drt_condition_utils.H"

ADAPTER::FluidProjection::FluidProjection(
        Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<LINALG::Solver> solverp,
        Teuchos::RCP<ParameterList> params,
        Teuchos::RCP<IO::DiscretizationWriter> output,
        bool isale,
        bool dirichletcond)
: 	fluid_(dis,*solver,*solverp,*params,*output,isale),
dis_(dis),
solver_(solver),
solverp_(solverp),
params_(params),
output_(output)
{
            DRT::UTILS::SetupNDimExtractor(*dis,"FSICoupling",interface_);
            DRT::UTILS::SetupNDimExtractor(*dis,"FREESURFCoupling",freesurface_);

            fluid_.SetFSISurface(&interface_);
            fluid_.SetFreeSurface(&freesurface_);

            // build inner velocity map
            // dofs at the interface are excluded
            // we use only velocity dofs and only those without Dirichlet constraint
            const Teuchos::RCP<const LINALG::MapExtractor> dbcmaps = fluid_.DirichMaps();
            std::vector<Teuchos::RCP<const Epetra_Map> > maps;
            maps.push_back(interface_.OtherMap());
            maps.push_back(dbcmaps->OtherMap());
            innervelmap_ = LINALG::MultiMapExtractor::IntersectMaps(maps);

            if (dirichletcond)
            {
                // mark all interface velocities as dirichlet values
                fluid_.AddDirichCond(interface_.CondMap());
            }
}

        /*----------------------------------------------------------------------*/
        /*----------------------------------------------------------------------*/
        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::InitialGuess()
        {
            return fluid_.InitialGuess();
        }

        /*----------------------------------------------------------------------*/
        /*----------------------------------------------------------------------*/
        void ADAPTER::FluidProjection::PrepareTimeStep()
        {
            fluid_.PrepareTimeStep();
            // we add the whole fluid mesh displacement later on?
            //fluid_.Dispnp()->PutScalar(0.);
        }

        /*----------------------------------------------------------------------*/
        /*----------------------------------------------------------------------*/
        void ADAPTER::FluidProjection::Update()
        {
            fluid_.TimeUpdate();
        }

        /*----------------------------------------------------------------------*/
        /*----------------------------------------------------------------------*/
        void ADAPTER::FluidProjection::Output()
        {
            fluid_.Output();
        }

        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::RHS()
        {
            dserror("not implemented");
            return fluid_.TrueResidual();
        }

        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::TrueResidual()
        {
            return fluid_.TrueResidual();
        }

        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::Velnp()
        {
            return fluid_.Velnp();
        }

        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::Velaf()
        {
            dserror("not implemented");
            return fluid_.Velnp();
        }

        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::Veln()
        {
            return fluid_.Veln();
        }

        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::Dispnp()
        {
            return fluid_.Dispnp();
        }

        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::SgVelVisc()
        {
            dserror("not implemented");
            return Teuchos::null;
        }

        Teuchos::RCP<const Epetra_Map>    ADAPTER::FluidProjection::DofRowMap()
        {
            const Epetra_Map* dofrowmap = dis_->DofRowMap();
            return Teuchos::rcp(dofrowmap, false);
        }

        Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FluidProjection::SystemMatrix()
        {
            dserror("not implemented");
            return Teuchos::null;
        }

        Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::FluidProjection::BlockSystemMatrix()
        {
            dserror("not implemented");
            return Teuchos::null;
        }

        Teuchos::RCP<LINALG::BlockSparseMatrixBase> ADAPTER::FluidProjection::ShapeDerivatives()
        {
            dserror("not implemented");
            return Teuchos::null;
        }

        Teuchos::RCP<DRT::Discretization> ADAPTER::FluidProjection::Discretization()
        {
            return fluid_.Discretization();
        }

        void ADAPTER::FluidProjection::Evaluate(Teuchos::RCP<const Epetra_Vector> vel)
        {
            dserror("not implemented");	return;
        }

        void ADAPTER::FluidProjection::StatisticsAndOutput()
        {
            // we don't support statistics yet
            // just write output
            fluid_.Output();
        }

        void ADAPTER::FluidProjection::NonlinearSolve()
        {
            fluid_.ProjectionSolve();
        }

        void ADAPTER::FluidProjection::Predictor()
        {
            dserror("not implemented"); return;
        }

        void ADAPTER::FluidProjection::MultiCorrector()
        {
            dserror("not implemented"); return;
        }

        Teuchos::RCP<const Epetra_Map>    ADAPTER::FluidProjection::InnerVelocityRowMap()
        {
            return innervelmap_;
        }

        Teuchos::RCP<const Epetra_Map>    ADAPTER::FluidProjection::VelocityRowMap()
        {
            return fluid_.VelocityRowMap();
        }

        Teuchos::RCP<const Epetra_Map>    ADAPTER::FluidProjection::PressureRowMap()
        {
            return fluid_.PressureRowMap();
        }

        void ADAPTER::FluidProjection::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
        {
            meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
        }

        double ADAPTER::FluidProjection::ResidualScaling() const
        {
            return fluid_.ResidualScaling();
        }

        double ADAPTER::FluidProjection::TimeScaling() const
        {
            if (params_->get<bool>("interface second order"))
            {
                return 2./fluid_.Dt();
            }
            else
                return 1./fluid_.Dt();
        }

        void ADAPTER::FluidProjection::SetInitialFlowField(int whichinitialfield,int startfuncno)
        {
            dserror("not implemented! todo"); return;
        }

        void ADAPTER::FluidProjection::SetTimeLomaFields(RCP<const Epetra_Vector> densnp,
                RCP<const Epetra_Vector> densn,
                RCP<const Epetra_Vector> densnm,
                RCP<const Epetra_Vector> scatraresidual,
                const int                numscal,
                const double             eosfac)
        {
            dserror("not implemented"); return;
        }

        void ADAPTER::FluidProjection::SetIterLomaFields(RCP<const Epetra_Vector> densnp,
                RCP<const Epetra_Vector> densdtnp,
                const int                numscal,
                const double             eosfac)
        {
            dserror("not implemented"); return;
        }

        void ADAPTER::FluidProjection::ReadRestart(int step)
        {
            fluid_.ReadRestart(step);
        }

        double ADAPTER::FluidProjection::Time() const
        {
            return fluid_.Time();
        }

        int ADAPTER::FluidProjection::Step() const
        {
            return fluid_.Step();
        }

        void ADAPTER::FluidProjection::LiftDrag()
        {
            dserror("not implemented"); return;
        }

        Teuchos::RCP<Epetra_Vector> ADAPTER::FluidProjection::ExtractInterfaceForces()
        {
            return interface_.ExtractCondVector(fluid_.TrueResidual());
        }

        Teuchos::RCP<Epetra_Vector> ADAPTER::FluidProjection::ExtractInterfaceForcesRobin()
        {
            // Calculate interface force from (externally applied) Robin force and
            // velocity. This assumes the fluid solve results in
            //
            // f_int - alpha_f*u(n+1) + f_robin = 0
            //
            // where f_robin consists of structural interface force and
            // displacement. The point here is to notice non-matching interface
            // displacements in the force vector, so that a testing of interface forces
            // is sufficient as convergence check.

            /*Teuchos::RCP<Epetra_Vector> robinforce = interface_.ExtractCondVector(fluid_.RobinRHS());
	double alphaf = params_->get<double>("alpharobinf",-1.);
	Teuchos::RCP<Epetra_Vector> ivelnp = interface_.ExtractCondVector(fluid_.Velnp());

	robinforce->Update(alphaf,*ivelnp,-1.0);

	return robinforce;*/
            dserror("not implemented");
            return Teuchos::null;
        }

        Teuchos::RCP<Epetra_Vector> ADAPTER::FluidProjection::ExtractInterfaceFluidVelocity()
        {
            return interface_.ExtractCondVector(fluid_.Velnp());
        }

        Teuchos::RCP<Epetra_Vector> ADAPTER::FluidProjection::ExtractInterfaceVeln()
        {
            return interface_.ExtractCondVector(fluid_.Veln());
        }

Teuchos::RCP<Epetra_Vector> ADAPTER::FluidProjection::ExtractFreeSurfaceVeln()
{
  dserror("not implemented");
  return Teuchos::null;
}

        void ADAPTER::FluidProjection::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
        {
            interface_.InsertCondVector(ivel,fluid_.Velnp());
        }

        void ADAPTER::FluidProjection::ApplyInterfaceRobinValue(Teuchos::RCP<Epetra_Vector> ivel, Teuchos::RCP<Epetra_Vector> iforce)
        {
            // use the known parts of structure field to create the robin
            // boundary value
            // the robin boundary value consists of a linear combination of
            // interface velocity and interface forces:

            // Robin-RHS = alpha_f * structural interface velocity
            //             - interface force (form structure to fluid)

            // get linear combination parameter
            /*double alphaf = params_->get<double>("alpharobinf",-1.);
	if (alphaf<0) dserror("falscher alpharobinf-Parameter");

	// robinboundaryvalue vorerst nur interfacegeschwindigkeit
	Teuchos::RCP<Epetra_Vector> robinboundaryvalue = Teuchos::rcp(new Epetra_Vector(*ivel));

	// at the moment iforce is the force to the structure, we have to
	// multiply with -1
	robinboundaryvalue->Update(-1.,*iforce,alphaf);

	// apply robin values to fluid equations RobinRHS vector
	interface_.InsertCondVector(robinboundaryvalue,fluid_.RobinRHS());
             */

            // at this point we have to omit the setting of dirichlet values at
            // the interface
            dserror("not implemented");
            return;
        }

        void ADAPTER::FluidProjection::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp)
        {
            meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

            // new grid velocity
            fluid_.UpdateGridv();
        }

        void ADAPTER::FluidProjection::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel)
        {
            meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
        }

        void ADAPTER::FluidProjection::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
        {
            // get interface velocity at t(n)
            const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractCondVector(Veln());

            // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
            //
            // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
            //
            double timescale = TimeScaling();
            fcx->Update(-timescale*fluid_.Dt(),*veln,timescale);
        }

        void ADAPTER::FluidProjection::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
        {
            // get interface velocity at t(n)
            const Teuchos::RCP<Epetra_Vector> veln = Interface().ExtractCondVector(Veln());

            // We convert Delta u(n+1,i+1) to Delta d(n+1,i+1) here.
            //
            // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
            //
            double timescale = 1./TimeScaling();
            fcx->Update(fluid_.Dt(),*veln,timescale);
        }

void ADAPTER::FluidProjection::FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  dserror("not implemented");
}

void ADAPTER::FluidProjection::FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  dserror("not implemented");
}

        int  ADAPTER::FluidProjection::Itemax() const
        {
            return fluid_.Itemax();
        }

        void ADAPTER::FluidProjection::SetItemax(int itemax)
        {
            fluid_.SetItemax(itemax);
        }

        Teuchos::RCP<Epetra_Vector> ADAPTER::FluidProjection::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
        {
            dserror("not implemented");
            return Teuchos::null;
        }

        Teuchos::RCP<const Epetra_Vector> ADAPTER::FluidProjection::ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
        {
            dserror("not implemented");
            return Teuchos::null;
        }

        Teuchos::RCP<Epetra_Vector> ADAPTER::FluidProjection::IntegrateInterfaceShape()
        {
            dserror("not implemented");
            return Teuchos::null;
        }

        void ADAPTER::FluidProjection::UseBlockMatrix(const LINALG::MultiMapExtractor& domainmaps,
                                                      const LINALG::MultiMapExtractor& rangemaps,
                                                      std::string condname,
                                                      bool splitmatrix)
        {
            dserror("not implemented");
            return;
        }

        Teuchos::RCP<DRT::ResultTest>  ADAPTER::FluidProjection::CreateFieldTest()
        {
            return Teuchos::rcp(new FLD::FluidResultTest(fluid_));
        }

#endif
