/*!----------------------------------------------------------------------
\file fluid_projectionmethod.H
\brief fluid solver based on pressure correction method

<pre>
Created on: Jun 10, 2009
Maintainer: Tobias Wiesner
wiesner@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FLUID_PROJECTIONMETHOD_H_
#define FLUID_PROJECTIONMETHOD_H_

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_io/io.H"

#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_function.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>

#include "fluid_utils_mapextractor.H"

using namespace std;
using namespace Teuchos;

/// fluid algorithms
namespace FLD
{

/// Pressure correction fluid solver
/*!
This class implements pressure correction fluid solver. Choose 'FLUID_SOLVER Pressure Correction' or
'FLUID_SOLVER Pressure Correction SemiImplicit' in the 'FLUID DYNAMIC' section of your 'dat' input file.
Only Taylor Hood elements (THQ9) are supported by pressure correction fluid solver. (see e.g. f2_channel_q2q1_drt.dat)

So far FluidProjectionMethod supports only Taylor Hood (Q2Q1) elements for fluid discretization
*/
class FluidProjectionMethod
{
    friend class FluidResultTest;

public:
    /*!
    * \brief Standard Constructor
    */
    FluidProjectionMethod(RCP<DRT::Discretization> dis,
            LINALG::Solver&		solver,
            LINALG::Solver&		psolver,
            ParameterList&		params,
            IO::DiscretizationWriter& output,
            bool				alefluid=false);

    LINALG::MapExtractor& VelPresSplitter() { return velpressplitter_; }
    Teuchos::RCP<const Epetra_Map> VelocityRowMap() { return velpressplitter_.OtherMap(); }
    Teuchos::RCP<const Epetra_Map> PressureRowMap() { return velpressplitter_.CondMap(); }

  /// introduce surface split extractor object
  /*!
    This method must (and will) be called during setup with a properly
    initialized extractor object iff we are on an ale mesh.
   */
  void SetSurfaceSplitter(const UTILS::MapExtractor* surfacesplitter) { surfacesplitter_ = surfacesplitter; }


    double ResidualScaling() const { return density_/dta_/theta_; }
    double Theta() const {return theta_; }
    double Dt() const { return dta_; }

    double Time() const { return time_; }
    int Step() const { return step_; }

    /*!
    * \brief start time loop
    */
    void Integrate();

    /*!
    * \brief do time integration (time loop)
    */
    void TimeLoop();
    void PrepareTimeStep();	//!< setup the variables to do a new time step
    void TimeUpdate();		//!< update solution. current solution becomes old solution of next timestep
    void Output();			//!< output of solution

    void ProjectionSolve();	//!< projection based solving strategy

    void CalcResidual();

    Teuchos::RCP<Epetra_Vector> InitialGuess() { dserror("InitialGuess?"); /*return incvel_;*/ return velnp_; }
    Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }
    Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }
    Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }
    Teuchos::RCP<Epetra_Vector> Dispnp() { return dispnp_; }
    Teuchos::RCP<Epetra_Vector> GridVel() { return gridv_; }
    Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

    void UpdateGridv();

    // only for compatibility
    // "max nonlin iter steps" is used for SolveImpulseEqn
    int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
    void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

    const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }
    void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    void ReadRestart(int step);

protected:
    void GenAlphaIntermediateValues();

    void ProjectionSolveImplicit();			//!< pressure correction method (implicit)
    void ProjectionSolveSemiImplicit();		//!< pressure correction method (semi-implicit)
    void getMatrices();						//!< calculates gradop, massmat and inverse lumped mass matrix
    void SolveImpulseEqn();					//!< calculate impulse eqn (implicit)
    void SolveImpulseEqnSemi();				//!< calculate impulse eqn (semi-implicit)
    void ProjectionStep();					//!< projection step
    void VelUpdate();						//!< update step

    // don't want = operator and cctor
    FluidProjectionMethod operator = (const FluidProjectionMethod& old);
    FluidProjectionMethod (const FluidProjectionMethod& old);

    RCP<DRT::Discretization>	discret_;
    LINALG::Solver&				solver_;
    LINALG::Solver&				psolver_;
    ParameterList&				params_;
    IO::DiscretizationWriter&	output_;
    int							myrank_;
    bool						alefluid_;

    FLUID_TIMEINTTYPE timealgo_;

    double dta_;
    double dtp_;
    double theta_;
    double alphaM_;
    double alphaF_;
    double gamma_;

    double  time_;
    int     step_;
    int     stepmax_;
    double  maxtime_;
    int     numstasteps_;
    bool    startalgo_;

    bool extrapolationpredictor_;

    int uprestart_;
    int upres_;

    string newton_;

    // MapExtractor velocity <-> pressure
    LINALG::MapExtractor      velpressplitter_;

    double density_;

    // matrices
    RCP<LINALG::SparseMatrix> pressmat_;
    RCP<LINALG::SparseMatrix> massmat_;
    RCP<LINALG::SparseMatrix> gradop_;
    RCP<LINALG::SparseMatrix> gradopwithoutbc_;	// i know, that's not perfect...

    RCP<Epetra_Vector> 	  trueresidual_;

    // vectors
    RCP<Epetra_Vector>    accn_;
    RCP<Epetra_Vector>    accnm_;
    RCP<Epetra_Vector>    accnp_;
    RCP<Epetra_Vector>    accam_;

    RCP<Epetra_Vector>    velnp_;
    RCP<Epetra_Vector>    veln_;
    RCP<Epetra_Vector>    velnm_;
    RCP<Epetra_Vector>    velaf_;

    RCP<Epetra_Vector>    vescnp_;

    RCP<Epetra_Vector>    dispnp_;
    RCP<Epetra_Vector>    dispn_;
    RCP<Epetra_Vector>    dispnm_;

    RCP<Epetra_Vector>    gridv_;

    RCP<Epetra_Vector>    hist_;

    RCP<Epetra_Vector> veltilde_;
    RCP<Epetra_Vector> phi_;
    RCP<Epetra_Vector> lmassinvvec_;

    RCP<Epetra_Vector>    zeros_;
    RCP<Epetra_Vector>    neumann_loads_;


    RCP<LINALG::MapExtractor> dbcmaps_;

  /// row dof map extractor
  const UTILS::MapExtractor* surfacesplitter_;
};

} // end namespace FLD
#endif /* FLUID_PROJECTIONMETHOD_H_ */
#endif /* CCADISCRET */
