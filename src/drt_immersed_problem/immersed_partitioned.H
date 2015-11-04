/*!----------------------------------------------------------------------
\file immersed_partitioned.H

\brief base class for all partitioned immersed algorithms

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_H
#define IMMERSED_PARTITIONED_H

#include "immersed_base.H"
#include "../drt_adapter/adapter_algorithmbase.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <NOX.H>
#include <NOX_Epetra.H>
#include <NOX_Epetra_Interface_Required.H>


namespace IMMERSED {

class ImmersedPartitioned : public ImmersedBase,
                            public ADAPTER::AlgorithmBase,
                            public NOX::Epetra::Interface::Required

{
public:

  explicit ImmersedPartitioned(const Epetra_Comm& comm);

  /// loop over time steps of partitioned scheme
  void Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface);

  /// compute residual of partitioned scheme S^{-1}(F(d)) - d
  bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// composed global nonlinear problem operator
    virtual void CouplingOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

protected:

    virtual void SetFieldDt(){};

    virtual void PrepareTimeStep(){};

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess()
  {dserror("not implemented in this class. may be overriden in derived class."); return Teuchos::null;};

  /// create linear solver framework
  Teuchos::RCP<NOX::Epetra::LinearSystem>
  CreateLinearSystem(Teuchos::ParameterList& nlParams,
                     const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface,
                     NOX::Epetra::Vector& noxSoln,
                     Teuchos::RCP<NOX::Utils> utils);

  /// create convergence tests including testing framework
  Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(Teuchos::ParameterList& nlParams,
                                                        Teuchos::RCP<NOX::Epetra::Group> grp);

  /// create convergence tests
  virtual void CreateStatusTest(Teuchos::ParameterList& nlParams,
                                Teuchos::RCP<NOX::Epetra::Group> grp,
                                Teuchos::RCP<NOX::StatusTest::Combo> converged);

  virtual void PrepareOutput();

  /// background operator
  void BackgroundOp(Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values,
                    const FillType fillFlag);

  /// immersed operator
  Teuchos::RCP<Epetra_Vector>
       ImmersedOp(Teuchos::RCP<Epetra_Vector> bdry_traction,
                  const FillType fillFlag);

  /// access to iteration counter
  std::vector<int> IterationCounter(){return counter_;};

  /// parameters handed in to NOX
  Teuchos::ParameterList noxparameterlist_;

protected:


private:

  /// setup list with default parameters
  virtual void SetDefaultParameters(const Teuchos::ParameterList& immerseddyn, Teuchos::ParameterList& list);


  /// print parameters and stuff
  /*!
    \warning This variable is only valid during when the time loop runs.
  */
  Teuchos::RCP<NOX::Utils> utils_;

  //! counters on how many times the residuum was called in a time step
  /*!
    NOX knows different types of residuum calls depending on
    circumstances (normal, finite difference, matrix free
    jacobi). It is possible to do approximations depending on the
    type.
  */
  std::vector<int> counter_;

  //! number of residuum calculations per nonlinear solve in one time step
  std::vector<int> linsolvcount_;

};// class ImmersedPartitioned
}// namespace immersed

#endif
