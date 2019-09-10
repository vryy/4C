/*----------------------------------------------------------------------------*/
/*! \file

\brief Solver for fluid field on a moving ALE mesh in combination with an immersed structure

\level 3

\maintainer Jonas Eichinger

*/
/*----------------------------------------------------------------------------*/
#include "ad_fld_fluid_ale_immersed.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::FluidAleImmersed::FluidAleImmersed(
    const Teuchos::ParameterList& prbdyn, std::string condname)
    : FluidAle(prbdyn, condname)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAleImmersed::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  FluidField()->AddDirichCond(maptoadd);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAleImmersed::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  FluidField()->RemoveDirichCond(maptoremove);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidAleImmersed::Output(const int step, const double time)
{
  int upres_ = 1;
  int uprestart_ = 0;
  bool alefluid_ = true;
  bool write_eledata_every_step_ = true;

  // write standard output
  if (step == -1 and time == -1.0)
    FluidField()->Output();
  else
  {
    if (FluidField()->Discretization()->Comm().MyPID() == 0)
      std::cout << "\n   Write EXTRA Fluid Output Step=" << step << " Time=" << time << " ...   \n"
                << std::endl;

    // output of solution
    if (FluidField()->Step() % upres_ == 0)
    {
      // step number and time
      FluidField()->DiscWriter()->NewStep(step, time);

      // time step, especially necessary for adaptive dt
      FluidField()->DiscWriter()->WriteDouble("timestep", FluidField()->Dt());

      // velocity/pressure vector
      FluidField()->DiscWriter()->WriteVector("velnp", FluidField()->Velnp());
      // (hydrodynamic) pressure
      Teuchos::RCP<Epetra_Vector> pressure =
          FluidField()->GetVelPressSplitter()->ExtractCondVector(FluidField()->Velnp());
      FluidField()->DiscWriter()->WriteVector("pressure", pressure);

      if (alefluid_) FluidField()->DiscWriter()->WriteVector("dispnp", FluidField()->Dispnp());

      // write domain decomposition for visualization (only once!)
      if ((FluidField()->Step() == upres_ or FluidField()->Step() == 0) and
          !write_eledata_every_step_)
        FluidField()->DiscWriter()->WriteElementData(true);
      else
        FluidField()->DiscWriter()->WriteElementData(true);

      if (uprestart_ != 0 && FluidField()->Step() % uprestart_ == 0)  // add restart data
      {
        // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
        FluidField()->DiscWriter()->WriteVector("accnp", FluidField()->Accnp());
        FluidField()->DiscWriter()->WriteVector("accn", FluidField()->Accn());
        FluidField()->DiscWriter()->WriteVector("veln", FluidField()->Veln());
        FluidField()->DiscWriter()->WriteVector("velnm", FluidField()->Velnm());

        if (alefluid_)
        {
          FluidField()->DiscWriter()->WriteVector("dispn", FluidField()->Dispn());
          FluidField()->DiscWriter()->WriteVector("gridvn", FluidField()->GridVeln());
        }

        // write mesh in each restart step --- the elements are required since
        // they contain history variables (the time dependent subscales)
        // But never do this for step 0 (visualization of initial field) since
        // it would lead to writing the mesh twice for step 0
        // (FluidBaseAlgorithm already wrote the mesh) -> HDF5 writer will claim!
        // if (( FluidField()->Step()!=0) and ((params_->sublist("RESIDUAL-BASED
        // STABILIZATION").get<std::string>("TDS")) != "quasistatic"))
        FluidField()->DiscWriter()->WriteMesh(FluidField()->Step(), FluidField()->Time());
      }
    }
    // write restart also when uprestart_ is not a integer multiple of upres_
    else if (uprestart_ > 0 && FluidField()->Step() % uprestart_ == 0)
    {
      // step number and time
      FluidField()->DiscWriter()->NewStep(step, time);

      // time step, especially necessary for adaptive dt
      FluidField()->DiscWriter()->WriteDouble("timestep", FluidField()->Dt());

      // velocity/pressure vector
      FluidField()->DiscWriter()->WriteVector("velnp", FluidField()->Velnp());

      if (alefluid_)
      {
        FluidField()->DiscWriter()->WriteVector("dispnp", FluidField()->Dispnp());
        FluidField()->DiscWriter()->WriteVector("dispn", FluidField()->Dispn());
        FluidField()->DiscWriter()->WriteVector("gridvn", FluidField()->GridVeln());
      }

      // write mesh in each restart step --- the elements are required since
      // they contain history variables (the time dependent subscales)
      // But never do this for step 0 (visualization of initial field) since
      // it would lead to writing the mesh twice for step 0
      // (FluidBaseAlgorithm already wrote the mesh) -> HDF5 writer will claim!
      // if (( FluidField()->Step()!=0) and ((params_->sublist("RESIDUAL-BASED
      // STABILIZATION").get<std::string>("TDS")) != "quasistatic"))
      FluidField()->DiscWriter()->WriteMesh(FluidField()->Step(), FluidField()->Time());

      // acceleration vector at time n+1 and n, velocity/pressure vector at time n and n-1
      FluidField()->DiscWriter()->WriteVector("accnp", FluidField()->Accnp());
      FluidField()->DiscWriter()->WriteVector("accn", FluidField()->Accn());
      FluidField()->DiscWriter()->WriteVector("veln", FluidField()->Veln());
      FluidField()->DiscWriter()->WriteVector("velnm", FluidField()->Velnm());
    }

    return;
  }
}
