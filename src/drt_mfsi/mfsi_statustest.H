/*----------------------------------------------------------------------*/
/*!
\file mfsi_statustest.H

\brief NOX Thyra group enhancement

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef MFSI_STATUSTEST_H
#define MFSI_STATUSTEST_H

#include <NOX_StatusTest_Generic.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>

#include <Epetra_Map.h>
#include <Epetra_Import.h>
#include <Epetra_Vector.h>

#include "../drt_fsi/fsi_coupling.H"


namespace MFSI
{

  /// a copy of NOX::StatusTest::NormF works on a vector part only
  /*!

    This is a copy of NOX::StatusTest::NormF that allows subclasses to
    define what vector to test. This way we can specify tests on
    vectors that cover only parts of our field.

    \author u.kue
    \date 11/07
   */
  class GenericNormF : public NOX::StatusTest::Generic
  {
  public:

    //! Type that determines whether to scale the norm by the problem size.
    enum ScaleType {
      //! No norm scaling
      Unscaled,
      //! Scale the norm by the length of the vector
      Scaled
    };

    //! Type that determines whether the norm is absolute or relative to the intial guess
    enum ToleranceType {
      //! Relative to starting guess
      Relative,
      //! Absolute
      Absolute
    };

    GenericNormF(std::string name,
                 double tolerance,
                 ScaleType stype = Scaled);

    /// do the check
    virtual NOX::StatusTest::StatusType
    checkStatus(const NOX::Solver::Generic& problem,
                NOX::StatusTest::CheckType checkType);

    /// get the result
    virtual NOX::StatusTest::StatusType getStatus() const;

    /// output
    virtual std::ostream& print(std::ostream& stream, int indent = 0) const;

    /* @name Accessor Functions
       Used to query current values of variables in the status test.
    */
    //@{

    //! Returns the value of the F-norm computed in the last call to checkStatus.
    virtual double getNormF() const;

    //! Returns the true tolerance.
    virtual double getTrueTolerance() const;

    //! Returns the specified tolerance set in the constructor.
    virtual double getSpecifiedTolerance() const;

    //! Returns the initial tolerance.
    virtual double getInitialTolerance() const;

    //@}

  protected:

    /*! \brief Calculate the norm of F for the given group according to
      the scaling type, norm type, and tolerance type.

      \note Returns -1.0 if F(x) has not been calculated for the given
      grp (i.e., grp.isF() is false).
    */
    virtual double computeNorm(const NOX::Abstract::Group& grp) = 0;

    /// compute the norm of a given vector
    double computeNorm(const Epetra_Vector& v);

    /*! In the case of a relative norm calculation, initializes
      \c trueTolerance based on the F-value at the initial guess.*/
    //void relativeSetup(NOX::Abstract::Group& initialGuess);

  private:

    //! %Status
    NOX::StatusTest::StatusType status_;

    //! Type of norm to use
    NOX::Abstract::Vector::NormType normType_;

    //! Scaling to use
    ScaleType scaleType_;

    //! Tolerance type (i.e., relative or absolute)
    ToleranceType toleranceType_;

    //! Tolerance required for convergence.
    double specifiedTolerance_;

    //! Initial tolerance
    double initialTolerance_;

    //! True tolerance value, i.e., specifiedTolerance / initialTolerance
    double trueTolerance_;

    //! Norm of F to be compared to trueTolerance
    double normF_;

    std::string name_;
  };


  /// Norm of inner dofs of one of our fields
  /*!

    For FSI status tests independant norms of the inner dofs of
    structural and fluid field are required.

    \author u.kue
    \date 11/07
   */
  class PartialNormF : public GenericNormF
  {
  public:

    PartialNormF(std::string name,
                 int blocknum,
                 const Epetra_Map &blockmap,
                 const Epetra_Map &innermap,
                 double tolerance,
                 ScaleType stype = Scaled);

  protected:

    double computeNorm(const NOX::Abstract::Group& grp);

  private:
    const Epetra_Map &blockmap_;
    const Epetra_Map &innermap_;
    int blocknum_;

    Teuchos::RCP<Epetra_Import> extractor_;
  };


  /// Norm of combined fluid-structure interface dofs
  /*!

    At the FSI interface the residual norm of the sum of fluid and
    structural residual is required to vanish.

    \author u.kue
    \date 11/07
   */
  class InterfaceNormF : public GenericNormF
  {
  public:

    InterfaceNormF(double structfac,
                   const Epetra_Map &structblockmap,
                   const Epetra_Map &structinterfacemap,
                   double fluidfac,
                   const Epetra_Map &fluidblockmap,
                   const Epetra_Map &fluidinterfacemap,
                   const FSI::Coupling& coupsf,
                   double tolerance,
                   ScaleType stype = Scaled);

  protected:

    double computeNorm(const NOX::Abstract::Group& grp);

  private:

    double structfac_;
    double fluidfac_;

    const Epetra_Map &structblockmap_;
    const Epetra_Map &structinterfacemap_;
    const Epetra_Map &fluidblockmap_;
    const Epetra_Map &fluidinterfacemap_;

    const FSI::Coupling& coupsf_;

    Teuchos::RCP<Epetra_Import> structextractor_;
    Teuchos::RCP<Epetra_Import> fluidextractor_;
  };
}

#endif
#endif
