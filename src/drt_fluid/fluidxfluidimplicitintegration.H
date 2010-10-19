/*----------------------------------------------------------------------*/
/*!
\file fluidxfluidimplicitintegration.H


<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef FLUIDXFLUIDIMPLICITINTEGRATION_H_
#define FLUIDXFLUIDIMPLICITINTEGRATION_H_

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_fluid.H"
#include "fluid_utils_mapextractor.H" // should go away

using namespace std;
using namespace Teuchos;

//forward declarations

namespace DRT
{
  class Discretization;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MultiMapExtractor;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}

/// fluid algorithms
namespace FLD
{
  
namespace UTILS
{
  class MapExtractor;
//  class KSPMapExtractor; 
}


class FluidXFluidImplicitTimeInt
{
  friend class FluidXFluidResultTest;

public:
  
  /*!
    \brief Standard Constructor

    */
  
  FluidXFluidImplicitTimeInt(RCP<DRT::Discretization>  actdis,
                             LINALG::Solver&           solver,
                             ParameterList&            params,
                             IO::DiscretizationWriter& output,
                             bool                      alefluid=false);
 
  /*!
  \brief Destructor

  */
  
  virtual  ~FluidXFluidImplicitTimeInt();
  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate();

  /*!
  \brief Do time integration (time loop)

   */
  void TimeLoop();

  /*!
  \brief Solve stationary problem
  
  */
  void SolveStationaryProblem();

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /*!
  \brief do nonlinear iteration, e.g. full Newton, Newton like or
           Fixpoint iteration

  */
  void NonlinearSolve();

   /*!
   \brief Update the solution after convergence of the nonlinear
           iteration. Current solution becomes old solution of next
           timestep.
   */
  void TimeUpdate();

  /*!
  \brief lift'n'drag forces, statistics time sample and
           output of solution and statistics

  */
  void StatisticsAndOutput();

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief Access output object

  */
  IO::DiscretizationWriter& DiscWriter()
  {
    return output_;
  }

  /*!
  \brief set initial flow field for analytical test problems
  */
  void SetInitialFlowField(
      const INPAR::FLUID::InitialField initfield,
      const int startfuncno
      );

   /*!
   \brief calculate stresses

   */
  Teuchos::RCP<Epetra_Vector> CalcStresses();

   /*!
   \brief calculate wall shear stresses

   */
  Teuchos::RCP<Epetra_Vector> CalcWallShearStresses();

   /*!
   \brief read restart data
    */
  void ReadRestart(int step);

   //! @name access methods for composite algorithms
   /// monolithic FSI needs to access the linear fluid problem

  Teuchos::RCP<Epetra_Vector> InitialGuess() { return incvel_; }
  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }
  Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }
  Teuchos::RCP<Epetra_Vector> Velaf() { return velaf_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }
  Teuchos::RCP<Epetra_Vector> GridVel() { return gridv_; }

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// Expand the Dirichlet DOF set
  ///
  /// The method expands the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is
  /// called by the staggered FSI in which the velocities on the FSI
  /// interface are prescribed by the other fields.
  void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  Teuchos::RCP<Epetra_Vector> Dispnp() { return dispnp_; }
  Teuchos::RCP<Epetra_Vector> Dispn() { return dispn_; }
  
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }
 
  LINALG::MapExtractor& VelPresSplitter(){ return velpressplitter_;};
  Teuchos::RCP<const Epetra_Map> VelocityRowMap();
  Teuchos::RCP<const Epetra_Map> PressureRowMap(); 
   
  double ResidualScaling() const
  { if (TimIntScheme()==INPAR::FLUID::timeint_stationary) return density_;
    else return density_/dta_/theta_;}
  double Dt() const { return dta_; }

  double Time() const { return time_; }
  int Step() const { return step_; }

  int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

  /*!
  \brief return type of time integration scheme

  */
  const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return timealgo_;}

 
  /// determine grid velocity
  void UpdateGridv();

  /// calculate lift&drag forces and angular momenta
  void LiftDrag() const;
    
   /// compute flow rate
  void ComputeFlowRates() const;

   /// integrate shape functions at nodes marked by condition
   /*!
     Needed for Mortar coupling at the FSI interface
    */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  protected:

   // don't want = operator and cctor
  FluidXFluidImplicitTimeInt operator = (const FluidXFluidImplicitTimeInt& old);
  FluidXFluidImplicitTimeInt (const FluidXFluidImplicitTimeInt& old);

   /*!
   \brief compute values at intermediate time steps for gen.-alpha

   */
  void GenAlphaIntermediateValues();

    /*!
    \brief call elements to calculate system matrix/rhs and assemble

    */
  void AssembleMatAndRHS();

    /*!
    \brief update acceleration for generalized-alpha time integration

    */
  void GenAlphaUpdateAcceleration();


    //! @name general algorithm parameters
    //! fluid discretization
  RCP<DRT::Discretization>      discret_;
  LINALG::Solver&               solver_;
  //! project nullspace (useful for completely Dirichlet bounded problems)
  bool                          project_;
  //! integrated basis vector (pressure component)
  Teuchos::RCP<Epetra_Vector>   w_;
  //! basis vector for nullspace (vector of ones in the pressure component)
  Teuchos::RCP<Epetra_Vector>   c_;
  ParameterList&                params_;
  IO::DiscretizationWriter&     output_;
  //! the processor ID from the communicator
  int                           myrank_;
  //! do we move the fluid mesh and calculate the fluid on this moving mesh?
  bool                          alefluid_;
  //@}

  /// number of spatial dimensions
  int numdim_;

  //! @name time stepping variables
  double  time_;        ///< physical time
  int     step_;        ///< timestep
  int     stepmax_;     ///< maximal number of timesteps
  double  maxtime_;     ///< maximal physical computation time
  int     numstasteps_; ///< number of steps for starting algorithm
  bool    startalgo_;   ///< flag for starting algorithm
  //@}

  /// time algorithm flag
  INPAR::FLUID::TimeIntegrationScheme timealgo_;

  /// constant density extracted from element material for incompressible flow
  /// (set to 1.0 for low-Mach-number flow)
  double density_;

  /// for low-Mach-number flow solver: thermodynamic pressure at n+alpha_F/n+1
  /// and at n+alpha_M/n as well as its time derivative at n+alpha_M/n
  double thermpressaf_;
  double thermpressam_;
  double thermpressdtam_;

  //! use (or not) linearisation of reactive terms on the element
  INPAR::FLUID::LinearisationAction newton_;
   /// whether we do some extrapolation at the begin of each time step (additional switch)
  bool extrapolationpredictor_;

  //! kind of predictor used in nonlinear iteration
  string  predictor_;

  //! @name restart variables
  int uprestart_;
  int upres_;
  int writestresses_;
  int write_wall_shear_stresses_;
  //@}

  //! @name time step sizes
  double dta_;
  double dtp_;
  //@}

  //! @name time-integration-scheme factors
  double theta_;
  double omtheta_;
  double alphaM_;
  double alphaF_;
  double gamma_;
  //@}

  /// flag for physical type of fluid flow (standard: incompressible)
  enum INPAR::FLUID::PhysicalType physicaltype_;

  /// flag for potential Neumann inflow
  bool neumanninflow_;

  /// form of convective term
  string convform_;

  /// fine-scale subgrid-viscosity flag
  string fssgv_;

  /// cpu-time measures
  double  dtele_;
  double  dtfilter_;
  double  dtsolve_;

  /// (standard) system matrix
  RCP<LINALG::SparseOperator> sysmat_;

  /// maps for extracting Dirichlet and free DOF sets
  RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// the vector containing potential Neumann-type outflow terms
  RCP<Epetra_Vector>    outflow_;

  /// (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>    residual_;

  /// rhs: right hand side vector for linearised flow solution
  RCP<Epetra_Vector>    rhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  RCP<Epetra_Vector>    incvel_;

  //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
  RCP<Epetra_Vector>    accnp_;
  RCP<Epetra_Vector>    accn_;
  RCP<Epetra_Vector>    accam_;
  //@}

  //! @name velocity and pressure at time n+1, n, n-1 and n+alpha_F
  RCP<Epetra_Vector>    velnp_;
  RCP<Epetra_Vector>    veln_;
  RCP<Epetra_Vector>    velnm_;
  RCP<Epetra_Vector>    velaf_;
   //@}

  //! @name scalar at time n+alpha_F/n+1 and n+alpha_M/n
  RCP<Epetra_Vector>    scaaf_;
  RCP<Epetra_Vector>    scaam_;
  //@}

  //! @name displacemets at time n+1, n and n-1
  RCP<Epetra_Vector>    dispnp_;
  RCP<Epetra_Vector>    dispn_;
  RCP<Epetra_Vector>    dispnm_;
  //@}

  /// grid velocity (set from the adapter!)
  RCP<Epetra_Vector>    gridv_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

  /// row dof map extractor
  const UTILS::MapExtractor* surfacesplitter_;

  /// Krylov space projection map extractor
  FLD::UTILS::KSPMapExtractor kspsplitter_;

   /// connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > >     pbcmapmastertoslave_;

  /// @name special relaxation state

  bool inrelaxation_;

  Teuchos::RCP<LINALG::SparseMatrix> dirichletlines_;

  Teuchos::RCP<LINALG::SparseMatrix> meshmatrix_;

  //@}

}; // class FluidXFluidImplicitTimeInt

}

#endif /* FLUIDXFLUIDIMPLICITINTEGRATION_H_ */
#endif  // #ifdef CCADISCRET
