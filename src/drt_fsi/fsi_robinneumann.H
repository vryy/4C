/*----------------------------------------------------------------------*/
/*!
\brief Solve FSI problems using a Robin partitioning approach

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FSI_ROBINNEUMANN_H
#define FSI_ROBINNEUMANN_H

#include "fsi_partitioned.H"
#include "fsi_statustest.H"


namespace FSI {

/// Robin-Neumann coupled interface system
/*!
  This is based on force coupling. Relaxation of interface forces is
  supported.
 */
class RobinNeumann : public Partitioned
{
public:

  explicit RobinNeumann(Epetra_Comm& comm);

protected:

  /// composed FSI operator
  virtual void FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

  /// interface fluid operator
  virtual Teuchos::RCP<Epetra_Vector> FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                                              const FillType fillFlag);

  /// interface structural operator
  virtual Teuchos::RCP<Epetra_Vector> StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                               const FillType fillFlag);

  virtual Teuchos::RCP<Epetra_Vector> InitialGuess();

  /// create convergance tests
  virtual void CreateStatusTest(ParameterList& nlParams,
                                Teuchos::RCP<NOX::Epetra::Group> grp,
                                Teuchos::RCP<NOX::StatusTest::Combo> converged);

private:

#if 0
  /// calculate interface force residual
  /*!
    At the FSI interface the residual norm of the sum of fluid and
    structural residual is required to vanish.

    \author u.kue
    \date 04/08
  */
  class InterfaceForceNormF : public GenericNormF
  {
  public:

    InterfaceForceNormF(RobinNeumann& algorithm,
                        std::string name,
                        double tolerance)
      : GenericNormF(name,tolerance,Scaled),
        algorithm_(algorithm)
    {}

  protected:

    double computeNorm(const NOX::Abstract::Group& grp);

  private:
    RobinNeumann& algorithm_;
  };
#endif

  /// (structural) interface force at begin of FSI step
  Teuchos::RCP<Epetra_Vector> iforcen_;
};

}

#endif
#endif
