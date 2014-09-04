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

#ifndef FLUIDIMPEDANCECONDITION_H
#define FLUIDIMPEDANCECONDITION_H


#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RCP.hpp"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"


namespace FLD
{
namespace UTILS
{

//--------------------------------------------------------------------
// Wrapper class (to be called from outside) for impedance bc
//--------------------------------------------------------------------

/*!
\brief impedance boundary condition wrapper
this class is meant to do some organisation stuff

\author chfoe (foerster@lnm.mw.tum.de)
*/
class FluidImpedanceWrapper
{
  friend class FluidImplicitTimeInt;


public:
  /*!
  \brief Standard Constructor
  */
  FluidImpedanceWrapper( Teuchos::RCP<DRT::Discretization>  actdis,
       IO::DiscretizationWriter& output,
       double dta);

  /*!
  \brief Destructor
  */
  virtual ~FluidImpedanceWrapper();

  /*!
  \brief Wrapper for FluidImpedanceBc::FlowRateCalculation
  */
  void FlowRateCalculation(double time, double dta);


  /*!
  \brief Wrapper for FluidImpedacnceBc::OutflowBoundary
  */
  void OutflowBoundary(double time, double dta, double theta);

  /*!
  \brief Wrapper for FluidImpedacnceBc::Impedances
  */
  void Impedances();

  /*!
  \brief Wrapper for FluidImpedacnceBc::UpdateResidual
  */
  void UpdateResidual(Teuchos::RCP<Epetra_Vector>  residual );

  /*!
  \brief Wrapper for FluidImpedacnceBc::WriteRestart
  */
  void WriteRestart( IO::DiscretizationWriter&  output );

  /*!
  \brief Wrapper for FluidImpedacnceBc::ReadRestart
  */
  void ReadRestart( IO::DiscretizationReader& reader);

  /*!
  \brief return pressures of one cycle
  */
  Teuchos::RCP<std::vector<double> > getPressures(int condid);

  /*!
  \brief Wrapper for FluidImpedacnceBc::SetWindkesselParams
  */
  void SetWindkesselParams(Teuchos::ParameterList & params, int condid);

  /*!
  \brief Wrapper for FluidImpedacnceBc::GetWindkesselParams
  */
  void GetWindkesselParams(Teuchos::ParameterList & params, int condid);

 /*!
  \brief Wrapper for FluidImpedacnceBc::getResultsAtEndOfPeriod
  */
  void getResultsOfAPeriod(Teuchos::ParameterList & params);

private:
  /*!
  \brief all single impedance conditions
  */
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> > impmap_;

  //! fluid discretization
  Teuchos::RCP<DRT::Discretization>      discret_;

  //! the output writer
  IO::DiscretizationWriter&     output_;

}; // class FluidImpedanceWrapper



//--------------------------------------------------------------------
// Actual impedance bc calculation stuff
//--------------------------------------------------------------------
/*!
\brief impedance boundary condition for vascular outflow boundaries

\author chfoe (foerster@lnm.mw.tum.de)
*/
class FluidImpedanceBc
{
  friend class FluidImpedanceWrapper;

public:

  /*!
  \brief Standard Constructor
  */
  FluidImpedanceBc( Teuchos::RCP<DRT::Discretization>  actdis,
        IO::DiscretizationWriter& output,
        double dta,
        int condid,
        int numcond );

  /*!
  \brief Empty Constructor
  */
  FluidImpedanceBc();

  /*!
  \brief Destructor
  */
  virtual ~FluidImpedanceBc() {}

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
  void FlowRateCalculation( double time,
          double dta,
          int condid );

  /*!
    \brief compute convolution integral and apply pressure
    to elements
  */
  void OutflowBoundary( double time,
      double dta,
      double theta,
      int condid );

  /*!
  \brief SetWindkesselParams
  */
  void SetWindkesselParams(Teuchos::ParameterList & params);

  /*!
  \brief GetWindkesselParams
  */
  void GetWindkesselParams(Teuchos::ParameterList & params);

  /*!
  \brief getResultsOfAPeriod
  */
  void getResultsOfAPeriod(Teuchos::ParameterList & params,
                           int             condid);

  void UpdateResidual( Teuchos::RCP<Epetra_Vector> residual );

private:
  /*!
  \brief calculate area at outflow boundary
  */
  double Area( double& density, double& viscosity, int numcond );

  /*!
    \brief determine impedances in frequency domain and time domain
  */
  void Impedances( double area, double density, double viscosity );

  /*!
    \brief determine impedances in frequency domain from windkessel model
  */
  std::complex<double> WindkesselImpedance(double k);

  /*!
  \brief compute tree and root impedance for given frequency for artery
  */
  std::complex<double> ArteryImpedance( int k,
          int generation,
          double radius,
          double termradius,
          double density,
          double viscosity,
          std::map<const double,std::complex<double> > zstored);

  /*!
  \brief compute tree and root impedance for zero frequency for artery
  */
  std::complex<double> DCArteryImpedance( int generation,
            double radius,
            double termradius,
            double density,
            double viscosity,
            std::map<const double,std::complex<double> > zstored);

  /*!
  \brief compute tree and root impedance for given frequency for lung
  */
  std::complex<double> LungImpedance( int k,
          int generation,
          double radius,
          double termradius,
          double density,
          double visosity,
          std::map<const double,std::complex<double> > zstored);

  /*!
  \brief compute tree and root impedance for zero frequency for lung
  */
  std::complex<double> DCLungImpedance( int generation,
              double radius,
              double termradius,
              double density,
              double viscosity,
              std::map<const double,std::complex<double> > zstored);


  /*!
  \brief return pressures of one cycle
  */
  Teuchos::RCP<std::vector<double> > getPressures(){return pressures_;}


  /*!
  \brief interpolate Vector1 to get Vector2
  */
  void interpolate(Teuchos::RCP<std::vector<double> > V1,
                   Teuchos::RCP<std::vector<double> > V2,
                   int  index1,
                   int & index2,
                   double time);


protected:

  Teuchos::RCP<Epetra_Vector> impedancetbc_;

private:

  //! ID of present condition
  int condid_;
  //! time period of present cyclic problem
  double period_;
  //! truncation radius at leafs of the tree
  double termradius_;
  //! condition type ( lung / artery / kidney / ... )
  //! implemented so far: lung, artery, windkessel
  std::string treetype_;

  //! 'material' parameters required for artery tree
  double R1_, R2_, C_, E_, H1_, H2_, H3_;

  //! time step size
  double dta_;

  //! the processor ID from the communicator
  int                           myrank_;

  //! fluid discretization
  Teuchos::RCP<DRT::Discretization>      discret_;

  //! the output writer
  IO::DiscretizationWriter&     output_;

  //! Flowrate array for impedance conditions
  Teuchos::RCP<std::vector<double> > flowrates_;

  //! pressures array for impedance conditions
  Teuchos::RCP<std::vector<double> > pressures_;

  //! Position at which the next element should be replaced
  //! initialised to zero as the first element will be replaced first
  int flowratespos_;

  //! initialised to zero as the first element will be replaced first
  int pressurespos_;

  //! Number of time steps which fit (exactly) into one period
  int cyclesteps_;

  //! Real impedance values in the time domain
  std::vector<double> impvalues_;   // real impedance values

  //! Information associated with requency indipendent windkessel
  //! Flow in at time step n+1
  double Qin_np_;

  //! Pressure in at time step n+1
  double Pin_np_;

  //! Pressure in at time step n
  double Pin_n_;

  //! Capcitors Pressure in at time step n+1
  double Pc_np_;

  //! Capcitors Pressure in at time step n
  double Pc_n_;

  //! Flag which indicated whether precalibrating bc is defined
  bool   IsPrecalibrated_;

  //! variable discribing the value difference between pressure at (n+1)T and at (n)T
  double dP_;

  //! Flag indicating the end of a cardiac cycle
  bool   endOfCycle_;

};  // class FluidImpedanceBc

}
}

#endif /* FLUIDIMPEDANCECONDITION_H */

