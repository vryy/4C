/*!----------------------------------------------------------------------
\file fluidimpedancecondition.H

\brief method to deal with three-element windkessel and other flow dependent pressure conditions

\level 3

<pre>
\maintainer  Moritz Thon
             thon@mhpc.mw.tum.de
</pre>

*----------------------------------------------------------------------*/

#ifndef FLUIDIMPEDANCECONDITION_H
#define FLUIDIMPEDANCECONDITION_H

#include "../drt_lib/drt_discret.H"

//forward declarations
namespace IO
{
  class DiscretizationReader;
}
namespace LINALG
{
  class MultiMapExtractor;
}


namespace FLD
{
namespace UTILS
{

class FluidImpedanceWrapper
{
  //friend class FLD::FluidImplicitTimeInt;

public:
  /*!
  \brief Standard Constructor
  */
  FluidImpedanceWrapper( const Teuchos::RCP<DRT::Discretization> actdis );

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
  \brief Calculate impedance tractions and add it to fluid residual and linearisation
  */
  void AddImpedanceBCToResidualAndSysmat(const double dta, Teuchos::RCP<Epetra_Vector>& residual, Teuchos::RCP<LINALG::SparseOperator>& sysmat );

  /*!
  \brief Wrap for time update of impedance conditions
  */
  void TimeUpdateImpedances(const double time);

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
  \brief Split linearization matrix to a BlockSparseMatrixBase
  */
  void UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool                             splitmatrix);

  /*!
    \brief compute and store flow rate of all previous
    time steps belonging to one cycle
  */
  void FlowRateCalculation( const int condid );

  /*!
    \brief compute convolution integral and apply pressure
    to elements
  */
  void CalculateImpedanceTractionsAndUpdateResidualAndSysmat(
      Teuchos::RCP<Epetra_Vector>& residual,
      Teuchos::RCP<LINALG::SparseOperator>& sysmat,
      const double dta,
      const int condid );

  /*!
    \brief Update flowrate and pressure vector
  */
  void TimeUpdateImpedance(const double time, const int condid);

  /*!
  \brief write flowrates_ and flowratespos_ to result files
  */
  void WriteRestart( IO::DiscretizationWriter&  output, const int condnum );

  /*!
  \brief read flowrates_ and flowratespos_
  */
  void ReadRestart( IO::DiscretizationReader& reader, const int condnum );

private:
  /*!
  \brief calculate area at outflow boundary
  */
  double Area( const int numcond );

  /*!
  \brief return relative error of last cycle
  */
  double getWKrelerror(){return WKrelerror_;}

private:

  //! fluid discretization
  const Teuchos::RCP<DRT::Discretization> discret_;

  //! the processor ID from the communicator
  const int myrank_;

  //! one-step theta time integration factor
  double theta_;

  //! condition type ( implemented so far: windkessel, ressistive )
  const std::string treetype_;

  //! time period of present cyclic problem
  const double period_;

  //! 'material' parameters required for artery tree
  const double R1_, R2_, C_;

  //! traction vector for impedance bc
  Teuchos::RCP<Epetra_Vector> impedancetbc_;

  //! linearisation of traction vector
  Teuchos::RCP<LINALG::SparseOperator> impedancetbcsysmat_;

  //! Pressure at time step n+1
  double P_np_;

  //! Pressure at time step n
  double P_n_;

  //! Flux at time step n+1
  double Q_np_;

  //! Flux at time step n
  double Q_n_;

  //! variable describing the relative error between pressure at (n+1)T and at (n)T
  double WKrelerror_;

  //! Pressure at beginning of the period
  double P_0_;
};  // class FluidImpedanceBc

}
}

#endif /* FLUIDIMPEDANCECONDITION_H */

