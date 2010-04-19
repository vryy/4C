/*!----------------------------------------------------------------------
\file statmech_time.H

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15234
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STATMECH_TIME_H
#define STATMECH_TIME_H

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_structure/strugenalpha.H"
#include "../drt_lib/drt_globalproblem.H"
/*!
\brief time integration for structural problems with statistical mechanics

\author cyron (cyron@lnm.mw.tum.de)
*/

// forward declarations
class StatMechManager;


class StatMechTime : public virtual StruGenAlpha
{
public:

  //! \brief Standard Constructor
  explicit StatMechTime(ParameterList& params,
                   DRT::Discretization& dis,
                   LINALG::Solver& solver,
                   IO::DiscretizationWriter& output);

  //! \brief Destructor
  virtual ~StatMechTime() {};

  //! \brief Do time integration
  void Integrate();

  //! \brief do consistent predictor step
  virtual void ConsistentPredictor(RCP<Epetra_MultiVector> randomnumbers);

  //!\brief do full Newton iteration
  virtual void FullNewton(RCP<Epetra_MultiVector> randomnumbers);

  //!\brief solve non-linear problem by means of PTC method
  virtual void PTC(RCP<Epetra_MultiVector> randomnumbers);

  //!\brief writing output including restart information
  virtual void Output();

  //!\brief read restart information
  virtual void ReadRestart(int step);

  /*!
  \brief adjusted DBC evaluation for periodic case

  \param params					 a parameter list
  \date 03/10
   */
  void EvaluateDirichletPeriodic(ParameterList& params);

  /*!
  \brief Determine Dirichlet condition at given time and apply its
         values to a system vector

  \param usetime
  \param time            Evaluation time
  \param nodeids				 vector containing DBC Node IDs
  \param curve					 vector containing curve values for DOFs
  \param funct					 vector containing funct values for DOFs
  \param onoff					 vector containing onoff values for DOFs
  \param val						 vector containing double multiplicator values for DOFs (amplitude)
  \date 03/10
  */
  void DoDirichletConditionPeriodic(const bool                  usetime,
																		const double                time,
																		vector<int>*								nodeids,
																		vector<int>*								curve,
																		vector<int>*								funct,
																		vector<int>*								onoff,
																		vector<double>*							val);


protected:

  // don't want = operator and cctor
  StatMechTime operator = (const StatMechTime& old);
  StatMechTime(const StatMechTime& old);

  //! \brief manager for administration of statistical mechanics stuff
  RCP<StatMechManager> statmechmanager_;

  //! \brief maximal number of random numbers to be generated in each time step per element
  int maxrandomnumbersperglobalelement_;

  //! \brief zero if in current time step the iterative solution of the nonlinear system of equations has not converged
  int   isconverged_;

  //! \brief how often a time step had to be restarted with new random numbers due to unconverged solution of nonlinear system of equations
  int   unconvergedsteps_;

private:
	// variables needed for EvaluateDirichletPeriodic()
	bool 												isinit_;
	// \brief vector used to signal when to store, change or delete a new/given reference position for DBC evaluation
	RCP<Epetra_Vector> 					dbctype_;
	// \brief incremental displacement which is to be added to the displacement in DoDirichletConditionPeriodic()
	RCP<Epetra_Vector> 					deltadbc_;
	// \brief direction of oscillation, extracted from discretization
	int 												oscdir_;
	// \brief amplitude of oscillation, extracted from discretization
	double											amp_;
	// \brief curve number, extracted from discretization
	int 												curvenumber_;

}; // class StatMechTime

#endif  // #ifndef STATMECH_H
#endif  // #ifdef CCADISCRET
