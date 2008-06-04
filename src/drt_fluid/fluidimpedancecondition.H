/*!----------------------------------------------------------------------
\file fluidimpedancecondition.H
\brief Method to deal with inpediance bc applied to vascular outflow bcs

<pre>
Maintainer: Christiane Förster
            foerster@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef FLUIDIMPEDANCECONDITION_H
#define FLUIDIMPEDANCECONDITION_H


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_elementgroup.H"

using namespace std;
using namespace Teuchos;


/*!
\brief impedance boundary condition for vascular outflow boundaries

\author chfoe (foerster@lnm.mw.tum.de)
*/
class FluidImpedanceBc
{
  friend class FluidImplicitTimeInt;

public:

  /*!
  \brief Standard Constructor
  */
  FluidImpedanceBc( RCP<DRT::Discretization>  actdis, 
		    IO::DiscretizationWriter& output,
		    double dta);

  /*!
  \brief Empty Constructor
  */
  FluidImpedanceBc();

  /*!
  \brief Destructor
  */
  virtual ~FluidImpedanceBc();

  /*!
  \brief write flowrates_ and flowratespos_ to result files
  */
  void WriteRestart( IO::DiscretizationWriter&  output );

  /*!
  \brief read flowrates_ and flowratespos_
  */
  void ReadRestart( IO::DiscretizationReader& reader );

  /*!
  \brief calculate area at outflow boundary
  */
  double Area( double& density, double& viscosity );

  /*!
    \brief compute and store flow rate of all previous 
    time steps belonging to one cycle 
  */
  void FlowRateCalculation( double time,
			    double dta);

  /*!
    \brief determine impedances in frequency domain and time domain
  */
  void Impedances( double area, double density, double viscosity );

  /*!
    \brief compute convolution integral and apply pressure
    to elements
  */
  void OutflowBoundary( double time,
			double dta,
			double theta);

  void UpdateResidual(RCP<Epetra_Vector>  residual );

private:
  /*!
  \brief compute tree and root impedance for given frequency for artery
  */
  std::complex<double> ArteryImpedance( int k, 
					int generation, 
					double radius,
					double termradius,
					double density,
					double viscosity,
					std::complex<double> zparent);

  /*!
  \brief compute tree and root impedance for zero frequency for artery
  */
  std::complex<double> DCArteryImpedance( int generation,  
					  double radius,
					  double termradius,
					  double density,
					  double viscosity,
					  std::complex<double> zparentdc);

  /*!
  \brief compute tree and root impedance for given frequency for lung
  */
  std::complex<double> LungImpedance( int k, 
					int generation,
					double radius,
					double termradius,
					double density,
					double visosity,
					std::complex<double> zparent);

  /*!
  \brief compute tree and root impedance for zero frequency for lung
  */
  std::complex<double> DCLungImpedance( int generation,
  					  double radius,
  					  double termradius,
  					  double density,
  					  double viscosity,
					  std::complex<double> zparentdc);

protected:

  RCP<Epetra_Vector> impedancetbc_;

private:

  //! number of impedance boundary condition lines in input
  int numcondlines_;
  //! number
  int numcond_;

  //! ID of present condition
  int condid_;
  //! time period of present cyclic problem
  double period_;
  //! condition type ( lung / artery / kidney / ... )
  //! implemented so far: lung, artery
  string treetype_;

  //! 'material' parameters required for artery tree
  double k1_, k2_, k3_;

  //! the processor ID from the communicator
  int                           myrank_;

  //! fluid discretization
  RCP<DRT::Discretization>      discret_;

  //! the output writer
  IO::DiscretizationWriter&     output_;

  //! Flowrate array for impedance conditions
  RCP<std::vector<double> > flowrates_;

  //! Position at which the next element should be replaced
  //! initialised to zero as the first element will be replaced first
  int flowratespos_;

  //! Number of time steps which fit (exactly) into one period
  int cyclesteps_;

  //! Real impedance values in the time domain
  vector<double> impvalues_;   // real impedance values

  //! histvector --- a linear combination of velnm, veln (BDF)
  //!                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

};  // class FluidImpedanceBc

#endif /* FLUIDIMPEDANCECONDITION_H */

#endif /* CCADISCRET */
