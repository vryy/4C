/*--------------------------------------------------------------------------*/
/*!
\file ehl_partitioned.H

\brief class for partitioned elastohydrodynamic lubrication (lubrication structure interaction)

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_EHL_EHL_PARTITIONED_H_
#define SRC_DRT_EHL_EHL_PARTITIONED_H_


#include "ehl_base.H"

namespace EHL {

class Partitioned : public Base
{
public:

  /// setup EHL algorithm
  explicit Partitioned(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& lubricationparams,
      const Teuchos::ParameterList& structparams,
      const std::string struct_disname,
      const std::string lubrication_disname);// Problem builder

  /// virtual destructor
  virtual ~Partitioned(){};

  /// setup of single fields (if needed)
  void SetupSystem(){};

  /// time loop of coupled problem
  void Timeloop();

protected:

  ///prepare time step of single fields
  void PrepareTimeStep();

  //! perform iteration loop between fields
  void OuterLoop();

  //! update time step and print to screen
  void UpdateAndOutput();

  //! convergence check of outer loop
  bool ConvergenceCheck(int itnum);

  ///do one inner iteration loop step of the structure
  void DoStructStep();

  ///do one inner iteration loop step of the lubrication
  void DoLubricationStep();

  /// velocity calculation given the displacements
  Teuchos::RCP<Epetra_Vector> CalcVelocity( Teuchos::RCP<const Epetra_Vector> dispnp );

  //! pressure increment of the outer loop
  Teuchos::RCP<Epetra_Vector> preincnp_;
  //! displacement increment of the outer loop
  Teuchos::RCP<Epetra_Vector> dispincnp_;

  //! maximum iteration steps
  int itmax_;
  //! convergence tolerance
  double ittol_;


};

//! class to deal with displacement relaxated partioned EHL. Relaxation parameter is constant
class Partitioned_StrToLub_Relax : public Partitioned
{
public:

  //!explicit constructor
  explicit Partitioned_StrToLub_Relax(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& lubricationparams,
      const Teuchos::ParameterList& structparams,
      const std::string struct_disname,
      const std::string lubrication_disname);

  //!virtual destructor
  virtual ~Partitioned_StrToLub_Relax(){};

protected:
  //! perform iteration loop between fields with relaxed displacements
  void OuterLoop();

  //! //! calculate relaxation parameter
  virtual void CalcOmega( double& omega , const int itnum );

  //! relaxation parameter
  double omega_;
};

//! class to deal with displacement relaxated partioned EHL. Relaxation parameter is calculated via Aitken
class Partitioned_StrToLub_Aitken : public Partitioned_StrToLub_Relax
{
public:

  //!explicit constructor
  explicit Partitioned_StrToLub_Aitken(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& lubricationparams,
      const Teuchos::ParameterList& structparams,
      const std::string struct_disname,
      const std::string lubrication_disname);

  //!virtual destructor
  virtual ~Partitioned_StrToLub_Aitken(){};

protected:

  //! Calculate relaxation parameter via Aitken
  virtual void CalcOmega( double& omega , const int itnum );

  //! displacement increment of the outer loop
  Teuchos::RCP<Epetra_Vector> del_;

  //! displacement increment of the outer loop
  Teuchos::RCP<Epetra_Vector> delhist_;

  //! maximal relaxation parameter
  double maxomega_;

  //! minimal relaxation parameter
  double minomega_;
};

//! class to deal with pressure relaxated partioned EHL. Relaxation parameter is constant
class Partitioned_LubToStr_Relax : public Partitioned
{
public:

  //!explicit constructor
  explicit Partitioned_LubToStr_Relax(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& lubricationparams,
      const Teuchos::ParameterList& structparams,
      const std::string struct_disname,
      const std::string lubrication_disname);

  //!virtual destructor
  virtual ~Partitioned_LubToStr_Relax(){};

protected:
  //! perform iteration loop between fields with relaxed pressure
  void OuterLoop();

  //! calculate relaxation parameter
  virtual void CalcOmega( double& omega , const int itnum );

  //! relaxation parameter
  double omega_;
};

//! class to deal with pressure relaxated partioned EHL. Relaxation parameter is calculated via Aitken
class Partitioned_LubToStr_Aitken : public Partitioned_LubToStr_Relax
{
public:

  //!explicit constructor
  explicit Partitioned_LubToStr_Aitken(const Epetra_Comm& comm,
      const Teuchos::ParameterList& globaltimeparams,
      const Teuchos::ParameterList& lubricationparams,
      const Teuchos::ParameterList& structparams,
      const std::string struct_disname,
      const std::string lubrication_disname);

  //!virtual destructor
  virtual ~Partitioned_LubToStr_Aitken(){};

protected:

  //! Calculate relaxation parameter via Aitken
  virtual void CalcOmega( double& omega , const int itnum );

  //! displacement increment of the outer loop
  Teuchos::RCP<Epetra_Vector> del_;

  //! displacement increment of the outer loop
  Teuchos::RCP<Epetra_Vector> delhist_;

  //! maximal relaxation parameter
  double maxomega_;

  //! minimal relaxation parameter
  double minomega_;
};

}

#endif /* SRC_DRT_EHL_EHL_PARTITIONED_H_ */
