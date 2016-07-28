/*!----------------------------------------------------------------------
\file fluidimpedancecondition.H

\brief Method to deal with windkessel and other vascular BCs

\level 3

<pre>
\maintainer  Moritz Thon
             thon@mhpc.mw.tum.de
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
      const double dta);

  /*!
  \brief Destructor
  */
  virtual ~FluidImpedanceWrapper(){};

  /*!
    \brief Wrapper for FluidImpedacnceBc::UseBlockMatrix
  */
  void UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool                             splitmatrix);

  /*!
    \brief Wrapper for FluidImpedacnceBc::TimeUpdateImpedance
    */
  void PrepareTimeStepImpedances(const double time);

  /*!
  \brief Calculate impedance tractions and add it to fluid residual and linearisation
  */
  void AddImpedanceBCToResidualAndSysmat(const double time, const double dta, Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<LINALG::SparseOperator>& sysmat );

  /*!
  \brief Wrap for time update of impedance conditions
  */
  void TimeUpdateImpedances(const double time, const double dta);

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

private:
  /*!
  \brief all single impedance conditions
  */
  std::map<const int, Teuchos::RCP<class FluidImpedanceBc> > impmap_;

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
      const double dta,
      const int condid,
      DRT::Condition* impedancecond );

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
      const double dta,
      const int condid );

  /*!
    \brief Update flowrate and pressure vector
  */
  void TimeUpdateImpedance(const double time);

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
  void UpdateResidualAndSysmat( Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<LINALG::SparseOperator>& sysmat );

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
  double getWKrelerror(){return dP_/P_0_;}

  /*!
  \brief interpolate Vector1 to get Vector2
  */
  void interpolate(const Teuchos::RCP<const std::vector<double> > V1,
                   Teuchos::RCP<std::vector<double> > V2,
                   const int index1,
                   int & index2,
                   const double time);


private:

  //! fluid discretization
  const Teuchos::RCP<DRT::Discretization> discret_;

  //! the processor ID from the communicator
  const int myrank_;


  Teuchos::RCP<Epetra_Vector> impedancetbc_;

  Teuchos::RCP<LINALG::SparseOperator> impedancetbcsysmat_;

  //! time period of present cyclic problem
  const double period_;
  //! truncation radius at leafs of the tree
  const double termradius_;
  //! condition type ( lung / artery / kidney / ... )
  //! implemented so far: lung, artery, windkessel
  const std::string treetype_;

  //! 'material' parameters required for artery tree
  double R1_, R2_, C_, E_, H1_, H2_, H3_;

  //! time step size
  double dta_;

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
  //! Pressure at time step n+1
  double P_np_;

  //! Pressure at time step n
  double P_n_;

  //! Capcitors Pressure at time step n+1
  double Pc_np_;

  //! Capcitors Pressure at time step n
  double Pc_n_;

  //! Flux at time step n+1
  double Q_np_;

  //! Flux at time step n
  double Q_n_;

  //! variable describing the relative error between pressure at (n+1)T and at (n)T
  double dP_;

  //! Pressure at beginning of the period
  double P_0_;

};  // class FluidImpedanceBc

}
}

#endif /* FLUIDIMPEDANCECONDITION_H */

