/*!----------------------------------------------------------------------
\file fluid_coupling_red_models.H
\brief Method to deal with coupling between 3D fluid and 0/1D vascular
problem

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
//#ifdef D_COUPLED_ARTNET
#ifdef D_ARTNET
#ifndef FLUID_COUPLING_RED_MODELS_H_
#define FLUID_COUPLING_RED_MODELS_H_


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_art_net/art_net_dyn_drt.H"
#include "../drt_io/io.H"

using namespace std;
using namespace Teuchos;

namespace FLD
{
namespace UTILS
{

//--------------------------------------------------------------------
// Wrapper class (to be called from outside) for coupling 3D/red D  bc
//--------------------------------------------------------------------

/*!
\brief coupling boundary condition wrapper
this class is meant to do some organisation stuff

\author ismail (ismail@lnm.mw.tum.de)
*/
class Fluid_couplingWrapper
{
  friend class FluidImplicitTimeInt;


public:
  /*!
  \brief Standard Constructor
  */
  Fluid_couplingWrapper( RCP<DRT::Discretization>  dis_3D,
                         RCP<DRT::Discretization>  dis_redD,
                         RCP<ART::ArtNetExplicitTimeInt>,
			 IO::DiscretizationWriter& output,
                         double dt_3d,
                         double dt_rm);

  /*!
  \brief Destructor
  */
  virtual ~Fluid_couplingWrapper();

  /*!
  \brief Wrapper for Fluid_couplingWrapper::FlowRateCalculation
  */
  void FlowRateCalculation(double time, double dta);

  /*!
  \brief Wrapper for Fluid_couplingWrapper::PressureCalculation
  */
  void PressureCalculation(double time, double dta);

  /*!
  \brief Wrapper for Fluid_couplingWrapper::OutflowBoundary
  */
  void ApplyBoundaryConditions(double time, double dta, double theta);

  /*!
  \brief Wrapper for Fluid_couplingWrapper::UpdateResidual
  */
  void UpdateResidual(RCP<Epetra_Vector>  residual );

  /*!
  \brief Wrapper for Fluid_couplingWrapper::WriteRestart
  */
  void WriteRestart( IO::DiscretizationWriter&  output );

  /*!
  \brief Wrapper for Fluid_couplingWrapper::ReadRestart
  */
  void ReadRestart( IO::DiscretizationReader& reader);


private:
  /*!
  \brief all single coupling conditions
  */
  map<const int, RCP<class Fluid_couplingBc> > coup_map3D_;

  //! map of coupling variables returned by the reduced-D model at time step n+1
  RCP<map<string, double > >                   mapRed_Dnp_;

  //! map of coupling variables returned by the reduced-D model at time step n
  RCP<map<string, double > >                   mapRed_Dn_;

  //! map of coupling variables returned by the 3-D model at time step n+1
  RCP<map<string, double > >                   map3_Dnp_;

  //! map of coupling variables returned by the 3-D model at time step n
  RCP<map<string, double > >                   map3_Dn_;

  //! 3D fluid discretization
  RCP<DRT::Discretization>                     discret3D_;

  //! Reduced-D artery network discretization
  RCP<DRT::Discretization>                     discret_redD_;

  //! Reduced-D artery network time integration
  RCP<ART::ArtNetExplicitTimeInt>              ArtExpTime_integ_;

  //! the output writer
  IO::DiscretizationWriter&                    output_;

  //! the fluid 3D time step size
  double                                       dt_f3_;

  //! the reduced model time step size
  double                                       dt_rm_;

}; // class Fluid_couplingWrapper



//--------------------------------------------------------------------
// Actual coupling bc calculation
//--------------------------------------------------------------------
/*!
\brief coupling boundary condition for vascular outflow boundaries

\author ismail (ismail@lnm.mw.tum.de)
*/
class Fluid_couplingBc
{
  friend class Fluid_couplingWrapper;

public:

  /*!
  \brief Standard Constructor
  */
  Fluid_couplingBc( RCP<DRT::Discretization>  dis3D,
                    RCP<DRT::Discretization>  dis_reD,
		    IO::DiscretizationWriter& output,
                    double dt_3d,
                    double dt_rm,
		    int condid,
		    int numcond,
                    int numcond2);

  /*!
  \brief Empty Constructor
  */
  Fluid_couplingBc();

  /*!
  \brief Destructor
  */
  virtual ~Fluid_couplingBc() {}

protected:
  /*!
  \brief write flowrates_ and flowratespos_ to result files
  */
  void WriteRestart( IO::DiscretizationWriter&  output, int condnum );

  /*!
  \brief read flowrates_ and flowratespos_
  */
  void ReadRestart( IO::DiscretizationReader& reader, int condnum );

  /*!
    \brief compute and store flow rate of all previous
    time steps belonging to one cycle
  */
  double FlowRateCalculation( double time,
                              double dta,
                              int condid );

  /*!
    \brief compute and store pressure of all previous
    time steps belonging to one cycle
  */
  double PressureCalculation( double time,
                              double dta,
                              int condid );


  /*!
    \brief compute convolution integral and apply pressure
    to elements
  */
  void OutflowBoundary(double pressure, 
                       double time,
                       double dta,
                       double theta,
                       int condid );

  void UpdateResidual( RCP<Epetra_Vector> residual );

private:
  /*!
  \brief calculate area at outflow boundary
  */
  double Area( double& density, double& viscosity, int numcond );



protected:

  RCP<Epetra_Vector> couplingbc_;

private:

  //! ID of present condition
  int condid_;

  //! 3D fluid time step size
  double dt_f3_;

  //! reduced-D time step size
  double dt_rm_;

  //! coupling error at the boundary
  double max_error_;

  //! number of maximum allowable iterations at the boundary
  double max_itr_;

  //! the processor ID from the communicator
  int                           myrank_;

  //! 3D fluid discretization
  RCP<DRT::Discretization>      discret_3D_;

  //! fluid discretization
  RCP<DRT::Discretization>      discret_redD_;

  //! the output writer
  IO::DiscretizationWriter&     output_;

  //! flow rate
  double flowrate_;

  //! pressure
  double pressure_;

};  // class Fluid_couplingBc

}
}

#endif /* D_ARTNET */
#endif /* FLUID_COUPLING_RED_MODELS_H_ */
#endif /* CCADISCRET */
