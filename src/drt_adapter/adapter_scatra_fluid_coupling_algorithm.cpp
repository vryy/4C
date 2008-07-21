/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_fluid_coupling_algorithm.cpp

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "adapter_scatra_fluid_coupling_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidCouplingAlgorithm::ScaTraFluidCouplingAlgorithm(
    Epetra_Comm& comm, 
    const Teuchos::ParameterList& prbdyn
    )
:  FluidBaseAlgorithm(prbdyn,false), // false -> no ALE in fluid algorithm
   ScaTraBaseAlgorithm(prbdyn),
   comm_(comm),
   params_(prbdyn),
   step_(0),
   time_(0.0)
{
  //print default parameters of parameter list
  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, prbdyn);

  // maximum simulation time
  maxtime_=prbdyn.get<double>("MAXTIME");
  // maximum number of timesteps
  nstep_ = prbdyn.get<int>("NUMSTEP");
  // time step size
  dt_ = prbdyn.get<double>("TIMESTEP");

  // get RCP to actual velocity field (time n+1)
  velocitynp_=FluidField().ExtractVelocityPart(FluidField().Velnp());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraFluidCouplingAlgorithm::~ScaTraFluidCouplingAlgorithm()
{
}


#endif // CCADISCRET
