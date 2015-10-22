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
  FluidImpedanceWrapper(
      const Teuchos::RCP<DRT::Discretization> actdis,
      const IO::DiscretizationWriter& output,
      const double dta);

  /*!
  \brief Destructor
  */
  virtual ~FluidImpedanceWrapper();

  /*!
  \brief Wrapper for FluidImpedanceBc::FlowRateCalculation
  */
  void FlowRateCalculation(const double time);

  /*!
  \brief Wrapper for FluidImpedacnceBc::OutflowBoundary
  */
  void OutflowBoundary(const double time, const double dta);

  /*!
    \brief Wrapper for FluidImpedacnceBc::TimeUpdateImpedance
    */
  void PrepareTimeStepImpedances(const double time);

  /*!
  \brief Wrap for time update of impedance conditions
  */
  void TimeUpdateImpedances(const double time, const double dta);

  /*!
  \brief Wrapper for FluidImpedacnceBc::Impedances
  */
  void Impedances();

  /*!
    \brief Wrapper for FluidImpedacnceBc::UseBlockMatrix
  */
  void UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool                             splitmatrix);

  /*!
  \brief Wrapper for FluidImpedacnceBc::UpdateResidual
  */
  void UpdateResidual(Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<LINALG::SparseOperator>& sysmat );

  /*!
  \brief Wrapper for FluidImpedacnceBc::WriteRestart
  */
  void WriteRestart( IO::DiscretizationWriter& output );

  /*!
  \brief Wrapper for FluidImpedacnceBc::ReadRestart
  */
  void ReadRestart( IO::DiscretizationReader& reader);

  /*!
  \brief return vector of relative pressure errors of last cycle
  */
  std::vector<double> getWKrelerrors( );

  /*!
  \brief Wrapper for FluidImpedacnceBc::SetWindkesselParams
  */
  void SetWindkesselParams(Teuchos::ParameterList & params, const int condid);

  /*!
  \brief Wrapper for FluidImpedacnceBc::GetWindkesselParams
  */
  void GetWindkesselParams(Teuchos::ParameterList & params, const int condid);

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
  const Teuchos::RCP<DRT::Discretization> discret_;

  //! the output writer
  const IO::DiscretizationWriter& output_;

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
  FluidImpedanceBc(
      const Teuchos::RCP<DRT::Discretization> actdis,
      const IO::DiscretizationWriter& output,
      const double dta,
      const int condid,
      const int numcond );

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
  void WriteRestart( IO::DiscretizationWriter&  output, const int condnum );

  /*!
  \brief read flowrates_ and flowratespos_
  */
  void ReadRestart( IO::DiscretizationReader& reader, const int condnum );

  /*!
    \brief compute and store flow rate of all previous
    time steps belonging to one cycle
  */
  void FlowRateCalculation(
      const double time,
      const int condid );

  /*!
    \brief compute convolution integral and apply pressure
    to elements
  */
  double OutflowBoundary(
      const double time,
      const double dta,
      const int condid );

  /*!
    \brief Update flowrate and pressure vector
  */
  void TimeUpdateImpedance(const double time, const double pressure);

  /*!
    \brief Time update impedance conditions
  */
  void PrepareTimeStepImpedance(const double time);

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
  void getResultsOfAPeriod(Teuchos::ParameterList & params, const int condid);

  /*!
  \brief Split linearization matrix to a BlockSparseMatrixBase
  */
  void UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool                             splitmatrix);

  /*!
  \brief Update of residual vector and its linearization
  */
  void UpdateResidual( Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<LINALG::SparseOperator>& sysmat );

private:
  /*!
  \brief calculate area at outflow boundary
  */
  double Area( double& density, double& viscosity, const int numcond );

  /*!
    \brief determine impedances in frequency domain and time domain
  */
  void Impedances( const double area, const double density, const double viscosity );

  /*!
    \brief determine impedances in frequency domain from windkessel model
  */
  std::complex<double> WindkesselImpedance(const double k);

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
  \brief return relative error of last cycle
  */
  double getWKrelerror(){return WKrelerror_;}

  /*!
  \brief interpolate Vector1 to get Vector2
  */
  void interpolate(const Teuchos::RCP<const std::vector<double> > V1,
                   Teuchos::RCP<std::vector<double> > V2,
                   const int index1,
                   int & index2,
                   const double time);


private:

  Teuchos::RCP<Epetra_Vector> impedancetbc_;

  Teuchos::RCP<LINALG::SparseOperator> impedancetbcsysmat_;

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
  int myrank_;

  //! fluid discretization
  const Teuchos::RCP<DRT::Discretization> discret_;

  //! the output writer
  const IO::DiscretizationWriter& output_;

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

  //! Information associated with frequency independent windkessel
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

  //! variable describing the relative error between pressure at (n+1)T and at (n)T
  double WKrelerror_;

  //! Flag indicating the end of a cardiac cycle
  bool   endOfCycle_;

};  // class FluidImpedanceBc

}
}

#endif /* FLUIDIMPEDANCECONDITION_H */

