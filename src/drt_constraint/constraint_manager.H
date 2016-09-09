/*!----------------------------------------------------------------------
\file constraint_manager.H

\brief Class controlling constraints and containing the necessary data, code originally by Thomas Kloeppel

\maintainer Marc Hirschvogel

\level 2

*----------------------------------------------------------------------*/

#ifndef CONSTRAINT_MANAGER_H
#define CONSTRAINT_MANAGER_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Teuchos_ParameterList.hpp>

namespace IO
{
  class DiscretizationReader;
}

namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class SparseOperator;
  class MultiMapExtractor;
}

namespace UTILS
{

  //forward declarations
  class Constraint;
  class ConstraintPenalty;
  class MPConstraint3;
  class MPConstraint3Penalty;
  class MPConstraint2;
  class Monitor;
  class ConstraintDofSet;

  /*!
  \brief Class controlling constraints and containing the necessary data

  \author Kloeppel (kloeppel@lnm.mw.tum.de)
  */
  class ConstrManager
  {
  public:

    //! Constructor of constraint manager, allocating the constraints
    ConstrManager();

     //! Destructor
    ~ConstrManager()
    {return;};

    //! initialize this class
    void Init(
        Teuchos::RCP<DRT::Discretization> discr,
        Teuchos::ParameterList params
        );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning Here Setup() needs to get an input argument because of the faulty implementation
             of this manager class. This needs to be fixed.

    \warning none
    \return void
    \date 09/16
    \author rauch  */
    void Setup(
        Teuchos::RCP<const Epetra_Vector> disp,
        Teuchos::ParameterList params);

    /*!
      \brief Change stiffness matrix and force vector according to the constraints.
      Values of lagrange multiplier are taken from intern variable.
      Difference between current and prescribed values is calculated and stored as well.
    */
    void EvaluateForceStiff
    (
      const double time,  ///< time at end of time step
      Teuchos::RCP<const Epetra_Vector> displast,  ///< displacement at beginning of time step
      Teuchos::RCP<const Epetra_Vector> disp,  ///< displacement at end of time step
      Teuchos::RCP<Epetra_Vector> fint,  ///< vector of internal forces
      Teuchos::RCP<LINALG::SparseOperator> stiff,  ///< stiffness matrix
      Teuchos::ParameterList scalelist
    );

    /*!
     \brief Return norm of difference between actual and constraint values
    */
    double GetErrorNorm() const
    {
      double foo;
      constrainterr_->Norm2(&foo);
      return foo;
    };

    /*!
         \brief Return number of constraints
    */
    int GetNumberOfConstraints() const
    {
        return numConstrID_;
    };

    /*!
     \brief Scale all lagrange multipliers by a double d
    */
    void ScaleLagrMult
    (
      double d   ///< scale factor
    )
    {
        lagrMultVec_->Scale(d);
        return;
    };

    /*!
         \brief Update constraint variables
    */
    void Update();

    /*!
         \brief Update lagrange multiplier \f$\lambda_{n+1}=\lambda_{n}+factor*\f$(volerr)
    */
    void UpdateLagrMult
    (
      double factor
    );

    /// Add a vector as residual increment to the vector of Lagrange multipliers
    void UpdateLagrMult
    (
      Teuchos::RCP<Epetra_Vector> vect  ///< vector to add
    );

    /// Add a vector as total increment to the vector of Lagrange multipliers
    void UpdateTotLagrMult
    (
      Teuchos::RCP<Epetra_Vector> vect  ///< vector to add
    );

    /*!
         \brief Compute difference between current and prescribed values at a given time and a given displacement
    */
    void ComputeError
    (
      double time,   ///< time, at which the error is to compute at
      Teuchos::RCP<Epetra_Vector> disp  ///< displacement vector at the given time
    );

    /*!
         \brief Return differences between prescribed and actual value of constraint number i
    */
    double GetError
    (
      int i ///< ID of constraint of interest
    ) const
    {
        return (*constrainterr_)[i];
    }

    /// return vector of differences between prescribed and actual values
    Teuchos::RCP<Epetra_Vector> GetError() const
    {
      return constrainterr_;
    }

    /*!
     \brief Return EpetraMap that determined distribution of constraints and lagrange
     multiplier over processors
    */
    Teuchos::RCP<Epetra_Map> GetConstraintMap() const
    {
      return constrmap_;
    };

    //! Return the additional rectangular matrix, constructed for lagrange multiplier evaluation
    Teuchos::RCP<LINALG::SparseOperator> GetConstrMatrix() //const
    {
      return constrMatrix_;
    };

    /*!
      \brief Return lagrange multiplier for constraint i
    */
    double GetLagrMult
    (
      int i ///< ID of constraint of interest
    ) const
    {
      return (*lagrMultVec_)[i];
    };

    /*!
      \brief Return lagrange multiplier vector
    */
    Teuchos::RCP<Epetra_Vector> GetLagrMultVector() const
    {
      return lagrMultVec_;
    };

    /*!
      \brief Return lagrange multiplier of last converged step
    */
    Teuchos::RCP<Epetra_Vector> GetLagrMultVectorOld() const
    {
      return lagrMultVecOld_;
    };

    /*!
     \brief Return if there are constraints
    */
    bool HaveConstraint() const
    {
        return haveconstraint_;
    };

    /*!
     \brief Return if there are constraints
    */
    bool HaveConstraintLagr() const
    {
        return havelagrconstr_;
    };

    /*!
     \brief Return if there are constraints
    */
    bool HaveConstraintPen() const
    {
        return havepenaconstr_;
    };

    /*!
       \brief Return if there are monitors
    */
    bool HaveMonitor() const
    {
        return havemonitor_;
    };

    /*!
     \brief Read restart information
    */
    void ReadRestart(
      IO::DiscretizationReader& reader,
      const double& time);

    /*!
     \brief Return current value
    */
    double GetCurrValue
    (
      int i ///< ID of constraint of interest
    ) const
    {
        return (*actvalues_)[i];
    };

    /*!
         \brief Print out the values of current monitor values
     */
    void PrintMonitorValues() const;

    /*!
       \brief Compute values described by a monitor boundary condition
    */
    void ComputeMonitorValues
    (
      Teuchos::RCP<Epetra_Vector> disp ///< current displacement
    );

    /*!
       \brief Compute values described by a monitor boundary condition
    */
    void ComputeMonitorValues
    (
      Teuchos::RCP<const Epetra_Vector> disp ///< current displacement
    );

    /// Reset reference base values for restart computations
    void SetRefBaseValues
    (
      Teuchos::RCP<Epetra_Vector> newrefvals,  ///< new reference base values
      const double& time   ///< current time
    );

    /// Reset lagrange multipliers
    void SetLagrMultVector
    (
      Teuchos::RCP<Epetra_Vector> newlagrmult  ///< new lagrange multipliers
    )
    {
      lagrMultVec_->Update(1.0,*newlagrmult,0.0);
      lagrMultVecOld_->Update(1.0,*newlagrmult,0.0);
      return;
    }

    /// Return Reference base values to write restart
    Teuchos::RCP<Epetra_Vector> GetRefBaseValues() const
    {
      return refbasevalues_;
    }

    //! switch constraint matrix to block matrix
    void UseBlockMatrix(Teuchos::RCP<const LINALG::MultiMapExtractor> domainmaps,
                        Teuchos::RCP<const LINALG::MultiMapExtractor> rangemaps);


  private:

    // don't want = operator, cctor and destructor
    ConstrManager operator = (const ConstrManager& old);
    ConstrManager(const ConstrManager& old);

    ///Build Monitor type Vector
    void BuildMoniType();

    Teuchos::RCP<DRT::Discretization> actdisc_; ///< discretization, elements to constraint live in
    Teuchos::RCP<ConstraintDofSet> constrdofset_; ///< degrees of freedom of lagrange multipliers
    Teuchos::RCP<Epetra_Map> constrmap_;  ///< unique map of constraint values
    Teuchos::RCP<Epetra_Map> redconstrmap_;   ///< fully redundant map of constraint values
    Teuchos::RCP<Epetra_Export> conimpo_;  ///< importer for fully redundant constraint vector into distributed one
    Teuchos::RCP<Epetra_Map> monitormap_;  ///< unique map of monitor values
    Teuchos::RCP<Epetra_Map> redmonmap_;  ///< fully redundant map of monitor values
    Teuchos::RCP<Epetra_Export> monimpo_;  ///< importer for fully redundant monitor vector into distributed one
    Teuchos::RCP<Epetra_Vector> referencevalues_;  ///< reference at current time step to constrain values to
    Teuchos::RCP<Epetra_Vector> refbasevalues_;   ///< reference base values at activation time of constrained structures
    Teuchos::RCP<Epetra_Vector> actvalues_;   ///< current values of constrained structures
    Teuchos::RCP<Epetra_Vector> constrainterr_;   ///< vector with deflection between reference and current values
    Teuchos::RCP<Epetra_Vector> monitorvalues_;  ///< current values of monitored structures
    Teuchos::RCP<Epetra_Vector> initialmonvalues_;   ///< initial values of monitored structures
    Teuchos::RCP<Epetra_Vector> monitortypes_;  ///< vector containing type of monitors
    int offsetID_;  ///< smallest constraint boundary condition ID
    int numConstrID_;  ///< number of constraint boundary conditions
    int numMonitorID_;  ///< smallest monitor boundary condition ID
    int minMonitorID_;  ///< number monitor boundary condition ID
    Teuchos::RCP<Epetra_Vector> fact_;  ///< vector with current time curve values
    Teuchos::RCP<Epetra_Vector> lagrMultVec_;  ///< lagrange multipliers
    Teuchos::RCP<Epetra_Vector> lagrMultVecOld_;  ///< lagrange multipliers
    Teuchos::RCP<LINALG::SparseOperator> constrMatrix_;  ///< additional rectangular matrix
    bool haveconstraint_;  ///< are there constraints at all?
    bool havelagrconstr_;  ///< are there constraints controlled by Lagrange multiplier?
    bool havepenaconstr_;  ///< are there constraints controlled by Penalty approach?
    bool havemonitor_;  ///< are there monitor conditions?
    double uzawaparam_;  ///< parameter of Uzawa algorithm (only for the case the linear uzawa is not used)

    Teuchos::RCP<Constraint> volconstr3d_;  ///< 3d volume constraints defined on surfaces
    Teuchos::RCP<Constraint> areaconstr3d_;  ///< 3d area constraints defined on surfaces
    Teuchos::RCP<Constraint> areaconstr2d_;  ///< 2d area constraints defined on lines
    Teuchos::RCP<ConstraintPenalty> volconstr3dpen_;
    Teuchos::RCP<ConstraintPenalty> areaconstr3dpen_;
    Teuchos::RCP<MPConstraint3> mpconplane3d_;  ///< 3d multipoint constraint prescribing the motion of a node relatively to a plane
    Teuchos::RCP<MPConstraint3> mpcnormcomp3d_;  ///< 3d multipoint constraint prescribing the motion of a node to a plane masternode
    Teuchos::RCP<MPConstraint2> mpconline2d_; ///< 2d multipoint constraint prescribing the motion of a node relatively to a straight line
    Teuchos::RCP<MPConstraint3Penalty> mpcnormcomp3dpen_;


    Teuchos::RCP<Monitor> volmonitor3d_; ///< 3d volume monitors defined on surfaces
    Teuchos::RCP<Monitor> areamonitor3d_;  ///< 3d area monitors defined on surfaces
    Teuchos::RCP<Monitor> areamonitor2d_;  ///< 2d area monitors defined on lines


  private:

    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

  protected:

    //! returns true if Setup() was called and is still valid
    bool IsSetup()
    {return issetup_;};

    //! returns true if Init(..) was called and is still valid
    bool IsInit()
    {return isinit_;};

    //! check if \ref Setup() was called
    void CheckIsSetup()
    {if(not IsSetup()) dserror("Setup() was not called.");};

    //! check if \ref Init() was called
    void CheckIsInit()
    {if(not IsInit()) dserror("Init(...) was not called.");};

  public:

    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse){issetup_ = trueorfalse;};

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse){isinit_ = trueorfalse;};

  }; //class
}
#endif /*CONSTRAINTMANAGER_H*/
