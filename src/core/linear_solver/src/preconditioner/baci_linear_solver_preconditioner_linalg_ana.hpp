/*----------------------------------------------------------------------*/
/*! \file

\brief A family of abstract nice algebra operations (ANA)

\level 1
*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_LINALG_ANA_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_LINALG_ANA_HPP

// Trilinos includes
#include "baci_config.hpp"

#include "baci_linalg_sparsematrix.hpp"
#include "baci_linear_solver_method.hpp"
#include "baci_linear_solver_method_linalg.hpp"
#include "baci_utils_exceptions.hpp"
#include "baci_utils_parameter_list.hpp"

#include <functional>

#define DEBUGGING_ANA 0  // turn on to get debugging printouts

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  namespace ANA
  {
    // forward declarations
    class LCSTimesVec;
    class Vector;

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A virtual class implementing an abstract implicit linear combination
           on a vector-valued quantity.

    */
    class LCBase
    {
     public:
      /// ctor
      LCBase() {}

      /// dtor
      virtual ~LCBase() = default;
      /*!
      \brief Return the range space of the result of the linear combination
      */
      virtual const Epetra_BlockMap& RangeMap() const = 0;

      /*!
      \brief Perform " v += scale * " operation

      \param v    (out): A vector with range map this->RangeMap()
                         with values to be updated on output
      \param scale (in): a scaling factor for the linear combination
                         (usually -1.0 or 1.0, used for sign changes)
      */
      virtual void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const = 0;

      /*!
      \brief Perform " v = scale * " operation

      \param v    (out): A vector with range map this->RangeMap() with values to be set on output
      \param scale (in): a scaling factor for the linear combination (usually -1.0 or 1.0, used for
      sign changes)
      */
      virtual void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const = 0;

    };  // class LCBase


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A pure virtual light weight wrapper for a (heavy-weight) Epetra_Operator

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

    */
    class LightWeightOperatorBase
    {
     public:
      /// ctor
      LightWeightOperatorBase() {}

      /// cctor
      LightWeightOperatorBase(const LightWeightOperatorBase& old) {}

      /// dtor
      virtual ~LightWeightOperatorBase() = default;
      /*!
      \brief The derived class shall return a clone of itself by calling its own copy-ctor
      */
      virtual Teuchos::RCP<LightWeightOperatorBase> Clone() const = 0;

      /*!
      \brief Use transpose of operator
      */
      virtual int SetUseTranspose(bool UseTranspose) = 0;

      /*!
      \brief Apply operator to X and return result in Y

      \note X and Y might be in-place pointing to the same physical Epetra_MultiVector
      */
      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

      /*!
      \brief Apply the inverse of the operator to X and return result in Y

      \note X and Y might be in-place pointing to the same physical Epetra_MultiVector
      */
      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

      /*!
      \brief return inf-norm of operator
      */
      virtual double NormInf() const = 0;

      /*!
      \brief return label of operator
      */
      virtual const char* Label() const = 0;

      /*!
      \brief return flag indicating whether transposed of operator is used in Apply and ApplyInverse
      */
      virtual bool UseTranspose() const = 0;

      /*!
      \brief return flag indicating whether operator supports inf-norm
      */
      virtual bool HasNormInf() const = 0;

      /*!
      \brief return communicator
      */
      virtual const Epetra_Comm& Comm() const = 0;

      /*!
      \brief return domain map of operator
      */
      virtual const Epetra_Map& OperatorDomainMap() const = 0;

      /*!
      \brief return range map of operator
      */
      virtual const Epetra_Map& OperatorRangeMap() const = 0;

     private:
    };  // class LightWeightOperatorBase



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A lightweight wrapper for a true heavy Epetra_Operator

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

     \sa LightWeightOperatorBase

    */
    class LightWeightOperator : public LightWeightOperatorBase
    {
     public:
      LightWeightOperator(const Epetra_Operator& op) : op_(op) {}

      LightWeightOperator(const LightWeightOperator& old)
          : LightWeightOperatorBase(old), op_(old.op_)
      {
      }

      Teuchos::RCP<LightWeightOperatorBase> Clone() const override
      {
        return Teuchos::rcp(new LightWeightOperator(*this));
      }


      int SetUseTranspose(bool UseTranspose) override
      {
        return const_cast<Epetra_Operator&>(op_).SetUseTranspose(UseTranspose);
      }

      int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
      {
        return op_.Apply(X, Y);
      }

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
      {
        return op_.ApplyInverse(X, Y);
      }

      double NormInf() const override { return op_.NormInf(); }

      const char* Label() const override { return "CORE::LINALG::ANA::LightWeightOperator"; }

      bool UseTranspose() const override { return op_.UseTranspose(); }

      bool HasNormInf() const override { return op_.HasNormInf(); }

      const Epetra_Comm& Comm() const override { return op_.Comm(); }

      const Epetra_Map& OperatorDomainMap() const override { return op_.OperatorDomainMap(); }

      const Epetra_Map& OperatorRangeMap() const override { return op_.OperatorRangeMap(); }

     private:
      const Epetra_Operator& op_;

    };  // class LightWeightOperator


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A lightweight wrapper implementing the transposedof a LightWeightOperatorBase

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

     \sa LightWeightOperatorBase

    */
    class OperatorTransposed : public LightWeightOperatorBase
    {
     public:
      OperatorTransposed(const LightWeightOperatorBase& op) : op_(op.Clone()) {}

      OperatorTransposed(const OperatorTransposed& old) : LightWeightOperatorBase(old), op_(old.op_)
      {
      }


      Teuchos::RCP<LightWeightOperatorBase> Clone() const override
      {
        return Teuchos::rcp(new OperatorTransposed(*this));
      }

      int SetUseTranspose(bool UseTranspose) override
      {
        // we are transposing the transposed operator
        return const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!UseTranspose);
      }

      int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
      {  // apply the transposed
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        int err = op_->Apply(X, Y);
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        return err;
      }

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
      {
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        int err = op_->ApplyInverse(X, Y);
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        return err;
      }

      double NormInf() const override
      {
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        double out = op_->NormInf();
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        return out;
      }

      const char* Label() const override { return "CORE::LINALG::ANA::OperatorTransposed"; }

      bool UseTranspose() const override { return (!op_->UseTranspose()); }

      bool HasNormInf() const override { return op_->HasNormInf(); }

      const Epetra_Comm& Comm() const override { return op_->Comm(); }

      const Epetra_Map& OperatorDomainMap() const override { return op_->OperatorRangeMap(); }

      const Epetra_Map& OperatorRangeMap() const override { return op_->OperatorDomainMap(); }

     private:
      const Teuchos::RCP<LightWeightOperatorBase> op_;

    };  // class OperatorTransposed

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A lightweight wrapper implementing a scalar-scaled LightWeightOperatorBase

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

     \sa LightWeightOperatorBase

    */
    class OperatorScaled : public LightWeightOperatorBase
    {
     public:
      OperatorScaled(const LightWeightOperatorBase& op, const double& scalar)
          : op_(op.Clone()), scalar_(scalar)
      {
      }

      OperatorScaled(const OperatorScaled& old)
          : LightWeightOperatorBase(old), op_(old.op_), scalar_(old.scalar_)
      {
      }


      Teuchos::RCP<LightWeightOperatorBase> Clone() const override
      {
        return Teuchos::rcp(new OperatorScaled(*this));
      }

      int SetUseTranspose(bool UseTranspose) override
      {
        return const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(UseTranspose);
      }

      int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
      {  // apply the transposed
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        int err = op_->Apply(X, Y);
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        Y.Scale(scalar_);
        return err;
      }

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
      {
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        int err = op_->ApplyInverse(X, Y);
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        Y.Scale(1. / scalar_);
        return err;
      }

      double NormInf() const override { return scalar_ * op_->NormInf(); }

      const char* Label() const override { return "CORE::LINALG::ANA::OperatorScaled"; }

      bool UseTranspose() const override { return (op_->UseTranspose()); }

      bool HasNormInf() const override { return op_->HasNormInf(); }

      const Epetra_Comm& Comm() const override { return op_->Comm(); }

      const Epetra_Map& OperatorDomainMap() const override { return op_->OperatorDomainMap(); }

      const Epetra_Map& OperatorRangeMap() const override { return op_->OperatorRangeMap(); }

     private:
      const Teuchos::RCP<LightWeightOperatorBase> op_;
      const double scalar_;

    };  // class OperatorScaled

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A lightweight wrapper implementing an implicit product of 2 LightWeightOperatorBase
    classes

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

     \sa LightWeightOperatorBase

    */
    class OperatorProduct : public LightWeightOperatorBase
    {
     public:
      OperatorProduct(const LightWeightOperatorBase& left, const LightWeightOperatorBase& right)
          : usetransposed_(false), left_(left.Clone()), right_(right.Clone())
      {
      }

      OperatorProduct(const OperatorProduct& old)
          : LightWeightOperatorBase(old),
            usetransposed_(old.usetransposed_),
            left_(old.left_),
            right_(old.right_)
      {
      }


      Teuchos::RCP<LightWeightOperatorBase> Clone() const override
      {
        return Teuchos::rcp(new OperatorProduct(*this));
      }

      int SetUseTranspose(bool UseTranspose) override
      {
        usetransposed_ = UseTranspose;
        return 0;
      }

      int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

      double NormInf() const override
      {
        FOUR_C_THROW(
            "CORE::LINALG::ANA::OperatorProduct does not implement "
            "LightWeightOperatorBase::NormInf()");
        return -1.0;
      }

      const char* Label() const override { return "CORE::LINALG::ANA::OperatorProduct"; }

      bool UseTranspose() const override { return usetransposed_; }

      bool HasNormInf() const override { return false; }

      const Epetra_Comm& Comm() const override { return left_->Comm(); }

      const Epetra_Map& OperatorDomainMap() const override
      {
        if (usetransposed_)
          return left_->OperatorRangeMap();
        else
          return right_->OperatorDomainMap();
      }

      const Epetra_Map& OperatorRangeMap() const override
      {
        if (usetransposed_)
          return right_->OperatorDomainMap();
        else
          return left_->OperatorRangeMap();
      }

     private:
      bool usetransposed_;
      const Teuchos::RCP<LightWeightOperatorBase> left_;
      const Teuchos::RCP<LightWeightOperatorBase> right_;

    };  // class OperatorProduct


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A lightweight wrapper implementing an implicit sum of 2 LightWeightOperatorBase classes

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

     \sa LightWeightOperatorBase

    */
    class OperatorSum : public LightWeightOperatorBase
    {
     public:
      OperatorSum(
          const LightWeightOperatorBase& left, const LightWeightOperatorBase& right, const int sign)
          : sign_(sign), usetransposed_(false), left_(left.Clone()), right_(right.Clone())
      {
        if (sign != 1 && sign != -1) FOUR_C_THROW("sign parameter has to be 1 or -1");
      }

      OperatorSum(const OperatorSum& old)
          : LightWeightOperatorBase(old),
            sign_(old.sign_),
            usetransposed_(old.usetransposed_),
            left_(old.left_),
            right_(old.right_)
      {
      }


      Teuchos::RCP<LightWeightOperatorBase> Clone() const override
      {
        return Teuchos::rcp(new OperatorSum(*this));
      }

      int SetUseTranspose(bool UseTranspose) override
      {
        usetransposed_ = UseTranspose;
        return 0;
      }

      int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
      {
        FOUR_C_THROW(
            "CORE::LINALG::ANA::OperatorSum does not implement "
            "LightWeightOperatorBase::ApplyInverse");
        return -1;
      }

      double NormInf() const override
      {
        FOUR_C_THROW(
            "CORE::LINALG::ANA::OperatorSum does not implement LightWeightOperatorBase::NormInf()");
        return -1.0;
      }

      const char* Label() const override { return "CORE::LINALG::ANA::OperatorSum"; }

      bool UseTranspose() const override { return usetransposed_; }

      bool HasNormInf() const override { return false; }

      const Epetra_Comm& Comm() const override { return left_->Comm(); }

      const Epetra_Map& OperatorDomainMap() const override
      {
        if (usetransposed_)
          return left_->OperatorRangeMap();
        else
          return left_->OperatorDomainMap();
      }

      const Epetra_Map& OperatorRangeMap() const override
      {
        if (usetransposed_)
          return left_->OperatorDomainMap();
        else
          return left_->OperatorRangeMap();
      }

     private:
      int sign_;
      bool usetransposed_;
      const Teuchos::RCP<LightWeightOperatorBase> left_;
      const Teuchos::RCP<LightWeightOperatorBase> right_;

    };  // class OperatorSum



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A lightweight wrapper implementing an implicit inverse of a LightWeightOperatorBase

    A supplied solver can be used to implement the inverse. If no solver is supplied, this operator
    uses a CORE::LINALG::Solver with no parameter list that defaults to Amesos_KLU.

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

     \sa LightWeightOperatorBase, CORE::LINALG::ANA::inverse

    */
    class OperatorInverse : public LightWeightOperatorBase
    {
     public:
      OperatorInverse(const Epetra_Operator& op, CORE::LINALG::Solver& solver, bool reset = true)
          : reset_(reset), solver_(solver), op_(op)
      {
      }

      OperatorInverse(
          const CORE::LINALG::SparseOperator& op, CORE::LINALG::Solver& solver, bool reset = true)
          : reset_(reset),
            solver_(solver),
            op_(*(const_cast<CORE::LINALG::SparseOperator&>(op).EpetraOperator()))
      {
      }

      OperatorInverse(const Epetra_Operator& op)
          : reset_(true),
            defaultsolver_(std::invoke(
                [&]()
                {
                  Teuchos::ParameterList solvparams;
                  CORE::UTILS::AddEnumClassToParameterList<CORE::LINEAR_SOLVER::SolverType>(
                      "SOLVER", CORE::LINEAR_SOLVER::SolverType::umfpack, solvparams);
                  return Teuchos::rcp(new CORE::LINALG::Solver(solvparams, op.Comm()));
                })),
            solver_(*defaultsolver_),
            op_(op)
      {
      }

      OperatorInverse(const CORE::LINALG::SparseOperator& op)
          : reset_(true),
            defaultsolver_(std::invoke(
                [&]()
                {
                  Teuchos::ParameterList solvparams;
                  CORE::UTILS::AddEnumClassToParameterList<CORE::LINEAR_SOLVER::SolverType>(
                      "SOLVER", CORE::LINEAR_SOLVER::SolverType::umfpack, solvparams);
                  return Teuchos::rcp(new CORE::LINALG::Solver(solvparams, op.Comm()));
                })),
            solver_(*defaultsolver_),
            op_(*(const_cast<CORE::LINALG::SparseOperator&>(op).EpetraOperator()))
      {
      }

      OperatorInverse(const OperatorInverse& old)
          : LightWeightOperatorBase(old),
            reset_(old.reset_),
            defaultsolver_(old.defaultsolver_),
            solver_(old.solver_),
            op_(old.op_)
      {
      }


      Teuchos::RCP<LightWeightOperatorBase> Clone() const override
      {
        return Teuchos::rcp(new OperatorInverse(*this));
      }

      int SetUseTranspose(bool UseTranspose) override
      {
        FOUR_C_THROW("CORE::LINALG::ANA::OperatorInverse does not support transpose");
        return -1;
      }

      int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override
      {
        FOUR_C_THROW(
            "CORE::LINALG::ANA::OperatorInverse does not support inverse of inverse of underlying "
            "operator, use Apply instead");
        return -1;
      }

      double NormInf() const override
      {
        FOUR_C_THROW(
            "CORE::LINALG::ANA::OperatorInverse does not support NormInf of inverse of operator");
        return -1.0;
      }

      const char* Label() const override { return "CORE::LINALG::ANA::OperatorInverse"; }

      bool UseTranspose() const override { return false; }

      bool HasNormInf() const override { return false; }

      const Epetra_Comm& Comm() const override { return op_.Comm(); }

      const Epetra_Map& OperatorRangeMap() const override  // no, this is NOT a bug
      {
        return op_.OperatorDomainMap();
      }

      const Epetra_Map& OperatorDomainMap() const override  // no, this is NOT a bug
      {
        return op_.OperatorRangeMap();
      }

     private:
      const bool reset_;
      Teuchos::RCP<CORE::LINALG::Solver> defaultsolver_;
      CORE::LINALG::Solver& solver_;
      const Epetra_Operator& op_;

    };  // class OperatorInverse


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A distributed vector class that implements Epetra_Vector

    All Epetra_Vector functionality can be used. Additionally, this class overloads a series
    of operators used in ANA linear algebra expressions.

    \sa LightWeightOperatorBase, Epetra_Vector

    */
    class Vector : public Epetra_Vector
    {
     public:
      /// Implements Epetra_Vector ctor
      Vector(const Epetra_BlockMap& m, bool init = true) : Epetra_Vector(m, init) {}

      /// Implements Epetra_Vector ctor
      Vector(const Vector& Source) : Epetra_SrcDistObject(Source), Epetra_Vector(Source) {}

      /// Implements Epetra_Vector ctor
      Vector(Epetra_DataAccess CV, const Epetra_BlockMap& m, double* V) : Epetra_Vector(CV, m, V) {}

      /// Implements Epetra_Vector ctor
      Vector(Epetra_DataAccess CV, const Epetra_MultiVector& mv, int i) : Epetra_Vector(CV, mv, i)
      {
      }

      /*!
      \brief Initialize this Vector from a scalar

      \param rhs (in): Scalar value to init this vector with
      */
      inline void operator=(const double& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator = (const double& rhs)" << endl;
        fflush(stdout);
#endif
        PutScalar(rhs);
      }

      /*!
      \brief Initialize this Vector from another Vector (deep copy)

      \param rhs (in): Vector to init this vector with
      */
      inline void operator=(const CORE::LINALG::ANA::Vector& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator = (const CORE::LINALG::ANA::Vector& rhs)" << endl;
        fflush(stdout);
#endif
        Update(1.0, rhs, 0.0);
      }

      /// Teuchos::RCP version of the above method
      inline void operator=(const Teuchos::RCP<CORE::LINALG::ANA::Vector>& rhs) { *this = *rhs; }

      /*!
      \brief Update this Vector with another Vector

      \param rhs (in): Vector to update this vector with
      */
      inline void operator+=(const CORE::LINALG::ANA::Vector& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator += (const CORE::LINALG::ANA::Vector& rhs)" << endl;
        fflush(stdout);
#endif
        Update(1.0, rhs, 1.0);
      }

      /// Teuchos::RCP version of the above method
      inline void operator+=(const Teuchos::RCP<CORE::LINALG::ANA::Vector>& rhs) { *this += *rhs; }

      /*!
      \brief Update this Vector with negative of another Vector

      \param rhs (in): Vector to update this vector with
      */
      inline void operator-=(const CORE::LINALG::ANA::Vector& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator -= (const CORE::LINALG::ANA::Vector& rhs)" << endl;
        fflush(stdout);
#endif
        Update(-1.0, rhs, 1.0);
      }

      /// Teuchos::RCP version of the above method
      inline void operator-=(const Teuchos::RCP<CORE::LINALG::ANA::Vector>& rhs) { *this -= *rhs; }

      /*!
      \brief Scale this Vector with a scalar

      \param scalar (in): Scalar the vector is scaled with
      */
      inline void operator*=(const double& scalar)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator *= (const double& scalar)" << endl;
        fflush(stdout);
#endif
        Scale(scalar);
      }

      /*!
      \brief Scale this Vector with the inverse of a scalar

      \param scalar (in): Scalar the vector is inverse-scaled with
      */
      inline void operator/=(const double& scalar)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator /= (const double& scalar)" << endl;
        fflush(stdout);
#endif
        Scale(1.0 / scalar);
      }

      /*!
      \brief Set this Vector to the result of a linear combination

      \param rhs (in): Linear combination of which the result is put into this Vector
      */
      inline void operator=(const CORE::LINALG::ANA::LCBase& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator = (const CORE::LINALG::ANA::LCBase& rhs)" << endl;
        fflush(stdout);
#endif
        rhs.Set(*this, 1.0);
      }

      /*!
      \brief Update this Vector with the result of a linear combination

      \param rhs (in): Linear combination of which the result is used to update this Vector
      */
      inline void operator+=(const CORE::LINALG::ANA::LCBase& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator += (CORE::LINALG::ANA::LCBase& rhs)" << endl;
        fflush(stdout);
#endif
        rhs.Update(*this, 1.0);
      }

      /*!
      \brief Update this Vector with the negative of the result of a linear combination

      \param rhs (in): Linear combination of which the result is used to negatively update this
      Vector
      */
      inline void operator-=(const CORE::LINALG::ANA::LCBase& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator -= (CORE::LINALG::ANA::LCBase& rhs)" << endl;
        fflush(stdout);
#endif
        rhs.Update(*this, -1.0);
      }

     private:
    };  // class  Vector



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Sum of 2 (generic) linear combinations

    \sa LCBase

    */
    class LcLcPlusLc : public LCBase
    {
     public:
      LcLcPlusLc(const LCBase& left, const LCBase& right) : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lc_plus_lc(const LCBase& left, const LCBase& right)" << endl;
        fflush(stdout);
#endif
      }

      ~LcLcPlusLc() override
      {
#if DEBUGGING_ANA
        cout << "~LC_lc_plus_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        left_.Update(v, scale);
        right_.Update(v, scale);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        left_.Set(v, scale);
        right_.Update(v, scale);
      }

     private:
      const LCBase& left_;
      const LCBase& right_;

    };  // class LC_lc_plus_lc


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of 2 (generic) linear combinations

    \sa LCBase

    */
    class LcLcMinusLc : public LCBase
    {
     public:
      LcLcMinusLc(const LCBase& left, const LCBase& right) : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lc_minus_lc(const LCBase& left, const LCBase& right)" << endl;
        fflush(stdout);
#endif
      }

      ~LcLcMinusLc() override
      {
#if DEBUGGING_ANA
        cout << "~LC_lc_minus_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        left_.Update(v, scale);
        right_.Update(v, -scale);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        left_.Set(v, scale);
        right_.Update(v, -scale);
      }

     private:
      const LCBase& left_;
      const LCBase& right_;

    };  // class LC_lc_minus_lc

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Sum of a Vector and a generic linear combination

    \sa LCBase

    */
    class LcVecPlusLc : public LCBase
    {
     public:
      LcVecPlusLc(const CORE::LINALG::ANA::Vector& vec, const LCBase& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_plus_lc(const CORE::LINALG::ANA::Vector& vec, const LCBase& right)" << endl;
        fflush(stdout);
#endif
      }

      ~LcVecPlusLc() override
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_plus_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale, vec_, 1.0);
        right_.Update(v, scale);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale, vec_, 0.0);
        right_.Update(v, scale);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec_;
      const LCBase& right_;

    };  // class LC_vec_plus_lc

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of a Vector and a generic linear combination

    \sa LCBase

    */
    class LcVecMinusLc : public LCBase
    {
     public:
      LcVecMinusLc(const CORE::LINALG::ANA::Vector& vec, const LCBase& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_minus_lc(const CORE::LINALG::ANA::Vector& vec, const LCBase& right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcVecMinusLc() override
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_minus_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale, vec_, 1.0);
        right_.Update(v, -scale);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale, vec_, 0.0);
        right_.Update(v, -scale);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec_;
      const LCBase& right_;

    };  // class LC_vec_minus_lc

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of a generic linear combination and a Vector

    \sa LCBase

    */
    class LcLcMinusVec : public LCBase
    {
     public:
      LcLcMinusVec(const LCBase& left, const CORE::LINALG::ANA::Vector& vec)
          : LCBase(), vec_(vec), left_(left)
      {
#if DEBUGGING_ANA
        cout << "LC_lc_minus_vec(const LCBase& left, const CORE::LINALG::ANA::Vector& vec)" << endl;
        fflush(stdout);
#endif
      }

      ~LcLcMinusVec() override
      {
#if DEBUGGING_ANA
        cout << "~LC_lc_minus_vec()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        left_.Update(v, scale);
        v.Update(-scale, vec_, 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        left_.Set(v, scale);
        v.Update(-scale, vec_, 1.0);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec_;
      const LCBase& left_;

    };  // class LC_lc_minus_vec

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Linear combination of a scalar with a Vector

    \sa LCBase

    */
    class LCSTimesVec : public LCBase
    {
     public:
      LCSTimesVec(const double scalar, const CORE::LINALG::ANA::Vector& vec)
          : LCBase(), scalar_(scalar), vec_(vec)
      {
#if DEBUGGING_ANA
        cout << "LCSTimesVec(const double scalar, const CORE::LINALG::ANA::Vector& vec)" << endl;
        fflush(stdout);
#endif
      }

      ~LCSTimesVec() override
      {
#if DEBUGGING_ANA
        cout << "~LCSTimesVec() " << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale * scalar_, vec_, 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale * scalar_, vec_, 0.0);
      }

      // give access to the scalar (for specialization LCs)
      inline const double& Scalar() const { return scalar_; }
      // give access to the vector (for specialization LCs)
      inline const CORE::LINALG::ANA::Vector& Vector() const { return vec_; }

     private:
      const double scalar_;
      const CORE::LINALG::ANA::Vector& vec_;

    };  // class LCSTimesVec

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Linear combination of a scalar with a generic linear combination

    \sa LCBase

    */
    class LcSTimesLc : public LCBase
    {
     public:
      LcSTimesLc(const double scalar, const CORE::LINALG::ANA::LCBase& right)
          : LCBase(), scalar_(scalar), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_s_times_lc(const double scalar, const CORE::LINALG::ANA::LCBase& right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcSTimesLc() override
      {
#if DEBUGGING_ANA
        cout << "~LC_s_times_lc() " << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return right_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        right_.Update(v, scale * scalar_);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        right_.Set(v, scale * scalar_);
      }

     private:
      const double scalar_;
      const CORE::LINALG::ANA::LCBase& right_;

    };  // class LC_s_times_lc


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Sum of 2 linear combinations, specialization of LC_lc_plus_lc
           for performance reasons

    \sa LCBase, LC_lc_plus_lc

    */
    class LcLcsvPlusLcsv : public LCBase
    {
     public:
      LcLcsvPlusLcsv(const LCSTimesVec& left, const LCSTimesVec& right)
          : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lcsv_plus_lcsv(const LCSTimesVec& left, const LCSTimesVec& right)" << endl;
        fflush(stdout);
#endif
      }

      ~LcLcsvPlusLcsv() override
      {
#if DEBUGGING_ANA
        cout << "~LC_lcsv_plus_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(
            scale * left_.Scalar(), left_.Vector(), scale * right_.Scalar(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(
            scale * left_.Scalar(), left_.Vector(), scale * right_.Scalar(), right_.Vector(), 0.0);
      }

     private:
      const LCSTimesVec left_;
      const LCSTimesVec right_;

    };  // class LC_lcsv_plus_lcsv


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Sum of Vector and a linear combinations, specialization of LC_vec_plus_lc
           for performance reasons

    \sa LCBase, LC_vec_plus_lc

    */
    class LcVecPlusLcsv : public LCBase
    {
     public:
      LcVecPlusLcsv(const CORE::LINALG::ANA::Vector& vec, const LCSTimesVec& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_plus_lcsv(const CORE::LINALG::ANA::Vector& vec, const LCSTimesVec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcVecPlusLcsv() override
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_plus_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale, vec_, scale * right_.Scalar(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale, vec_, scale * right_.Scalar(), right_.Vector(), 0.0);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec_;
      const LCSTimesVec right_;

    };  // class LC_vec_plus_lcsv

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of 2 linear combinations, specialization of LC_lc_minus_lc
           for performance reasons

    \sa LCBase, LC_lc_minus_lc

    */
    class LcLcsvMinusLcsv : public LCBase
    {
     public:
      LcLcsvMinusLcsv(const LCSTimesVec& left, const LCSTimesVec& right)
          : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lcsv_minus_lcsv(const LCSTimesVec& left, const LCSTimesVec& right)" << endl;
        fflush(stdout);
#endif
      }

      ~LcLcsvMinusLcsv() override
      {
#if DEBUGGING_ANA
        cout << "~LC_lcsv_minus_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(
            scale * left_.Scalar(), left_.Vector(), -scale * right_.Scalar(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(
            scale * left_.Scalar(), left_.Vector(), -scale * right_.Scalar(), right_.Vector(), 0.0);
      }

     private:
      const LCSTimesVec left_;
      const LCSTimesVec right_;

    };  // class LC_lcsv_minus_lcsv

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of a Vector and a linear combination, specialization of LC_vec_minus_lc
           for performance reasons

    \sa LCBase, LC_vec_minus_lc

    */
    class LcVecMinusLcsv : public LCBase
    {
     public:
      LcVecMinusLcsv(const CORE::LINALG::ANA::Vector& vec, const LCSTimesVec& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_minus_lcsv(const CORE::LINALG::ANA::Vector& vec, const LCSTimesVec& "
                "right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcVecMinusLcsv() override
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_minus_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale, vec_, -scale * right_.Scalar(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(scale, vec_, -scale * right_.Scalar(), right_.Vector(), 0.0);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec_;
      const LCSTimesVec right_;

    };  // class LC_vec_minus_lcsv

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of a linear combination and a Vector, specialization of LC_lc_minus_vec
           for performance reasons

    \sa LCBase, LC_lc_minus_vec

    */
    class LcLcsvMinusVec : public LCBase
    {
     public:
      LcLcsvMinusVec(const LCSTimesVec& left, const CORE::LINALG::ANA::Vector& vec)
          : LCBase(), vec_(vec), left_(left)
      {
#if DEBUGGING_ANA
        cout << "LC_lcsv_minus_vec(const LCSTimesVec& left, const CORE::LINALG::ANA::Vector& vec)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcLcsvMinusVec() override
      {
#if DEBUGGING_ANA
        cout << "~LC_lcsv_minus_vec()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(-scale, vec_, scale * left_.Scalar(), left_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Update(-scale, vec_, scale * left_.Scalar(), left_.Vector(), 0.0);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec_;
      const LCSTimesVec left_;

    };  // class LC_lcsv_minus_vec



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of 2 Vector s, result is a Vector, where
           Result[i] = vec1[i]*vec2[i]

    \sa LCBase

    */
    class LcVecPointwiseVec : public LCBase
    {
     public:
      LcVecPointwiseVec(
          const CORE::LINALG::ANA::Vector& vec1, const CORE::LINALG::ANA::Vector& vec2)
          : LCBase(), vec1_(vec1), vec2_(vec2)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_pointwise_vec(const CORE::LINALG::ANA::Vector& vec1, const "
                "CORE::LINALG::ANA::Vector& "
                "vec2)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcVecPointwiseVec() override
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_pointwise_vec()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return vec1_.Map(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Multiply(scale, vec1_, vec2_, 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Multiply(scale, vec1_, vec2_, 0.0);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec1_;
      const CORE::LINALG::ANA::Vector& vec2_;

    };  // class LC_vec_pointwise_vec


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of Vector and generic linear combination,
           result is a Vector, where Result[i] = vec[i]*right[i]

    \sa LCBase, LC_vec_pointwise_vec

    */
    class LcVecPointwiseLc : public LCBase
    {
     public:
      LcVecPointwiseLc(const CORE::LINALG::ANA::Vector& vec, const CORE::LINALG::ANA::LCBase& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_pointwise_lc(const CORE::LINALG::ANA::Vector& vec1, const "
                "CORE::LINALG::ANA::LCBase& "
                "right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcVecPointwiseLc() override
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_pointwise_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return vec_.Map(); }
      // perform 'v +=' operations
      void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override;
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        right_.Set(v, 1.0);
        v.Multiply(scale, vec_, v, 0.0);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec_;
      const CORE::LINALG::ANA::LCBase& right_;

    };  // class LC_vec_pointwise_lc



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of Vector and linear combination,
           result is a Vector, where Result[i] = vec[i]*right[i].
           Specialization of LC_vec_pointwise_lc for performance reasons

    \sa LCBase, LC_vec_pointwise_vec, LC_vec_pointwise_lc

    */
    class LcVecPointwiseLcsv : public LCBase
    {
     public:
      LcVecPointwiseLcsv(
          const CORE::LINALG::ANA::Vector& vec, const CORE::LINALG::ANA::LCSTimesVec& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_pointwise_lcsv(const CORE::LINALG::ANA::Vector& vec, const "
                "CORE::LINALG::ANA::LCSTimesVec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcVecPointwiseLcsv() override
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_pointwise_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Multiply(scale * right_.Scalar(), vec_, right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Multiply(scale * right_.Scalar(), vec_, right_.Vector(), 0.0);
      }

     private:
      const CORE::LINALG::ANA::Vector& vec_;
      const CORE::LINALG::ANA::LCSTimesVec right_;

    };  // class LC_vec_pointwise_lc



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of 2 linear combinations,
           result is a Vector, where Result[i] = left[i]*right[i].

    \sa LCBase, LC_vec_pointwise_vec

    */
    class LcLcPointwiseLc : public LCBase
    {
     public:
      LcLcPointwiseLc(const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::LCBase& right)
          : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lc_pointwise_lc(const CORE::LINALG::ANA::LCBase& left, const "
                "CORE::LINALG::ANA::LCBase& "
                "right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcLcPointwiseLc() override
      {
#if DEBUGGING_ANA
        cout << "~LC_lc_pointwise_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return left_.RangeMap(); }
      // perform 'v +=' operations
      void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override;
      // perform 'v =' operations
      void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override;

     private:
      const CORE::LINALG::ANA::LCBase& left_;
      const CORE::LINALG::ANA::LCBase& right_;

    };  // class LC_vec_pointwise_lc


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of 2 linear combinations,
           result is a Vector, where Result[i] = left[i]*right[i].
           Specialization of LC_lc_pointwise_lc for performance reasons

    \sa LCBase, LC_lc_pointwise_lc

    */
    class LcLcsvPointwiseLcsv : public LCBase
    {
     public:
      LcLcsvPointwiseLcsv(
          const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::LCSTimesVec& right)
          : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lcsv_pointwise_lcsv(const CORE::LINALG::ANA::LCSTimesVec& left, const "
                "CORE::LINALG::ANA::LCSTimesVec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcLcsvPointwiseLcsv() override
      {
#if DEBUGGING_ANA
        cout << "~LC_lcsv_pointwise_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Multiply(scale * left_.Scalar() * right_.Scalar(), left_.Vector(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override
      {
        v.Multiply(scale * left_.Scalar() * right_.Scalar(), left_.Vector(), right_.Vector(), 0.0);
      }

     private:
      const CORE::LINALG::ANA::LCSTimesVec left_;
      const CORE::LINALG::ANA::LCSTimesVec right_;

    };  // class LC_lcsv_pointwise_lcsv

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Apply linear operator to linear combination.
           Specialization of LC_Operator_times_lc for performance reasons.

    \sa LCBase, LC_Operator_times_lc

    */
    class LcOperatorTimesLcsv : public LCBase
    {
     public:
      LcOperatorTimesLcsv(
          const LightWeightOperatorBase& op, const CORE::LINALG::ANA::LCSTimesVec& right)
          : LCBase(), op_(op.Clone()), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_Operator_times_lcsv(const LightWeightOperatorBase& op, const "
                "CORE::LINALG::ANA::LCSTimesVec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcOperatorTimesLcsv() override
      {
#if DEBUGGING_ANA
        cout << "~LC_Operator_times_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return op_->OperatorRangeMap(); }
      // perform 'v +=' operations
      void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override;
      // perform 'v =' operations
      void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override;

     private:
      const Teuchos::RCP<LightWeightOperatorBase> op_;
      const CORE::LINALG::ANA::LCSTimesVec right_;

    };  // class LC_Operator_times_lcsv


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Apply linear operator to linear combination.

    \sa LCBase, LC_Operator_times_lcsv

    */
    class LcOperatorTimesLc : public LCBase
    {
     public:
      LcOperatorTimesLc(const LightWeightOperatorBase& op, const CORE::LINALG::ANA::LCBase& right)
          : LCBase(), op_(op.Clone()), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_Operator_times_lc(const LightWeightOperatorBase& op, const "
                "CORE::LINALG::ANA::LCBase& right)"
             << endl;
        fflush(stdout);
#endif
      }

      ~LcOperatorTimesLc() override
      {
#if DEBUGGING_ANA
        cout << "~LC_Operator_times_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      inline const Epetra_BlockMap& RangeMap() const override { return op_->OperatorRangeMap(); }
      // perform 'v +=' operations
      void Update(CORE::LINALG::ANA::Vector& v, const double& scale) const override;
      // perform 'v =' operations
      void Set(CORE::LINALG::ANA::Vector& v, const double& scale) const override;

     private:
      const Teuchos::RCP<LightWeightOperatorBase> op_;
      const CORE::LINALG::ANA::LCBase& right_;

    };  // class LC_Operator_times_lc



    /*----------------------------------------------------------------------*
       static (local) little helper method wrapping an Epetra_Operator
     *----------------------------------------------------------------------*/
    static inline CORE::LINALG::ANA::LightWeightOperator lw(const Epetra_Operator& op)
    {
      return CORE::LINALG::ANA::LightWeightOperator(op);
    }

    /*----------------------------------------------------------------------*
       scalar, vector and LC operations
     *----------------------------------------------------------------------*/
    inline CORE::LINALG::ANA::LCSTimesVec operator*(
        const double& scalar, const CORE::LINALG::ANA::Vector& vec)
    {
      return CORE::LINALG::ANA::LCSTimesVec(scalar, vec);
    }
    inline CORE::LINALG::ANA::LCSTimesVec operator*(
        const CORE::LINALG::ANA::Vector& vec, const double& scalar)
    {
      return CORE::LINALG::ANA::LCSTimesVec(scalar, vec);
    }
    inline CORE::LINALG::ANA::LCSTimesVec operator/(
        const CORE::LINALG::ANA::Vector& vec, const double& scalar)
    {
      return CORE::LINALG::ANA::LCSTimesVec(1. / scalar, vec);
    }
    inline CORE::LINALG::ANA::LcLcPlusLc operator+(
        const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcLcPlusLc(left, right);
    }
    inline CORE::LINALG::ANA::LcLcMinusLc operator-(
        const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcLcMinusLc(left, right);
    }
    inline CORE::LINALG::ANA::LcVecPlusLc operator+(
        const CORE::LINALG::ANA::Vector& vec, const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcVecPlusLc(vec, right);
    }
    inline CORE::LINALG::ANA::LcVecPlusLc operator+(
        const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::Vector& vec)
    {
      return CORE::LINALG::ANA::LcVecPlusLc(vec, left);
    }
    inline CORE::LINALG::ANA::LcVecMinusLc operator-(
        const CORE::LINALG::ANA::Vector& vec, const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcVecMinusLc(vec, right);
    }
    inline CORE::LINALG::ANA::LcLcMinusVec operator-(
        const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::Vector& vec)
    {
      return CORE::LINALG::ANA::LcLcMinusVec(left, vec);
    }
    inline CORE::LINALG::ANA::LcLcsvPlusLcsv operator+(
        const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::LcLcsvPlusLcsv(left, right);
    }
    inline CORE::LINALG::ANA::LcLcsvMinusLcsv operator-(
        const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::LcLcsvMinusLcsv(left, right);
    }
    inline CORE::LINALG::ANA::LcVecPlusLcsv operator+(
        const CORE::LINALG::ANA::Vector& vec, const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::LcVecPlusLcsv(vec, right);
    }
    inline CORE::LINALG::ANA::LcVecPlusLcsv operator+(
        const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::Vector& vec)
    {
      return CORE::LINALG::ANA::LcVecPlusLcsv(vec, left);
    }
    inline CORE::LINALG::ANA::LcVecMinusLcsv operator-(
        const CORE::LINALG::ANA::Vector& vec, const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::LcVecMinusLcsv(vec, right);
    }
    inline CORE::LINALG::ANA::LcLcsvMinusVec operator-(
        const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::Vector& vec)
    {
      return CORE::LINALG::ANA::LcLcsvMinusVec(left, vec);
    }
    inline CORE::LINALG::ANA::LcLcsvPlusLcsv operator+(
        const CORE::LINALG::ANA::Vector& left, const CORE::LINALG::ANA::Vector& right)
    {
      return (1.0 * left + 1.0 * right);
    }
    inline CORE::LINALG::ANA::LcLcsvMinusLcsv operator-(
        const CORE::LINALG::ANA::Vector& left, const CORE::LINALG::ANA::Vector& right)
    {
      return (1.0 * left - 1.0 * right);
    }
    inline CORE::LINALG::ANA::LcSTimesLc operator*(
        const double& scalar, const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcSTimesLc(scalar, right);
    }
    inline CORE::LINALG::ANA::LcSTimesLc operator*(
        const CORE::LINALG::ANA::LCBase& left, const double& scalar)
    {
      return (scalar * left);
    }
    // Teuchos::RCP versions of the above operations
    inline CORE::LINALG::ANA::LCSTimesVec operator*(
        const double& scalar, const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return (scalar * (*vec));
    }
    inline CORE::LINALG::ANA::LCSTimesVec operator*(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec, const double& scalar)
    {
      return (scalar * (*vec));
    }
    inline CORE::LINALG::ANA::LCSTimesVec operator/(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec, const double& scalar)
    {
      return ((*vec) / scalar);
    }
    inline CORE::LINALG::ANA::LcVecPlusLc operator+(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec,
        const CORE::LINALG::ANA::LCBase& right)
    {
      return ((*vec) + right);
    }
    inline CORE::LINALG::ANA::LcVecPlusLc operator+(const CORE::LINALG::ANA::LCBase& left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return (left + (*vec));
    }
    inline CORE::LINALG::ANA::LcVecMinusLc operator-(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec,
        const CORE::LINALG::ANA::LCBase& right)
    {
      return (*vec - right);
    }
    inline CORE::LINALG::ANA::LcLcMinusVec operator-(const CORE::LINALG::ANA::LCBase& left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return (left - (*vec));
    }
    inline CORE::LINALG::ANA::LcVecPlusLcsv operator+(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec,
        const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return (*vec + right);
    }
    inline CORE::LINALG::ANA::LcVecPlusLcsv operator+(const CORE::LINALG::ANA::LCSTimesVec& left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return (left + (*vec));
    }
    inline CORE::LINALG::ANA::LcVecMinusLcsv operator-(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec,
        const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return (*vec - right);
    }
    inline CORE::LINALG::ANA::LcLcsvMinusVec operator-(const CORE::LINALG::ANA::LCSTimesVec& left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return (left - (*vec));
    }
    inline CORE::LINALG::ANA::LcLcsvPlusLcsv operator+(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> right)
    {
      return ((*left) + (*right));
    }
    inline CORE::LINALG::ANA::LcLcsvMinusLcsv operator-(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> right)
    {
      return (*left - (*right));
    }


    /*----------------------------------------------------------------------*
       LightWeightOperatorBase and  Epetra_Operator operations
     *----------------------------------------------------------------------*/
    inline CORE::LINALG::ANA::LcOperatorTimesLc operator*(
        const CORE::LINALG::ANA::LightWeightOperatorBase& op,
        const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcOperatorTimesLc(op, right);
    }
    inline CORE::LINALG::ANA::LcOperatorTimesLc operator*(
        const Epetra_Operator& op, const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcOperatorTimesLc(CORE::LINALG::ANA::lw(op), right);
    }
    inline CORE::LINALG::ANA::LcOperatorTimesLcsv operator*(
        const CORE::LINALG::ANA::LightWeightOperatorBase& op,
        const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::LcOperatorTimesLcsv(op, right);
    }
    inline CORE::LINALG::ANA::LcOperatorTimesLcsv operator*(
        const CORE::LINALG::ANA::LightWeightOperatorBase& op, const CORE::LINALG::ANA::Vector& vec)
    {
      return op * (1.0 * vec);
    }
    inline CORE::LINALG::ANA::LcOperatorTimesLcsv operator*(
        const Epetra_Operator& op, const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::LcOperatorTimesLcsv(CORE::LINALG::ANA::lw(op), right);
    }
    inline CORE::LINALG::ANA::LcOperatorTimesLcsv operator*(
        const Epetra_Operator& op, const CORE::LINALG::ANA::Vector& vec)
    {
      return CORE::LINALG::ANA::lw(op) * (1.0 * vec);
    }
    inline CORE::LINALG::ANA::OperatorScaled operator*(
        const CORE::LINALG::ANA::LightWeightOperatorBase& op, const double& scalar)
    {
      return CORE::LINALG::ANA::OperatorScaled(op, scalar);
    }
    inline CORE::LINALG::ANA::OperatorScaled operator*(
        const double& scalar, const CORE::LINALG::ANA::LightWeightOperatorBase& op)
    {
      return CORE::LINALG::ANA::OperatorScaled(op, scalar);
    }
    inline CORE::LINALG::ANA::OperatorScaled operator*(
        const Epetra_Operator& op, const double& scalar)
    {
      return CORE::LINALG::ANA::OperatorScaled(CORE::LINALG::ANA::lw(op), scalar);
    }
    inline CORE::LINALG::ANA::OperatorScaled operator*(
        const double& scalar, const Epetra_Operator& op)
    {
      return CORE::LINALG::ANA::OperatorScaled(CORE::LINALG::ANA::lw(op), scalar);
    }
    inline CORE::LINALG::ANA::OperatorProduct operator*(
        const CORE::LINALG::ANA::LightWeightOperatorBase& left,
        const CORE::LINALG::ANA::LightWeightOperatorBase& right)
    {
      return CORE::LINALG::ANA::OperatorProduct(left, right);
    }
    inline CORE::LINALG::ANA::OperatorProduct operator*(
        const Epetra_Operator& left, const Epetra_Operator& right)
    {
      return CORE::LINALG::ANA::OperatorProduct(
          CORE::LINALG::ANA::lw(left), CORE::LINALG::ANA::lw(right));
    }
    inline CORE::LINALG::ANA::OperatorProduct operator*(
        const CORE::LINALG::ANA::LightWeightOperatorBase& left, const Epetra_Operator& right)
    {
      return CORE::LINALG::ANA::OperatorProduct(left, CORE::LINALG::ANA::lw(right));
    }
    inline CORE::LINALG::ANA::OperatorProduct operator*(
        const Epetra_Operator& left, const CORE::LINALG::ANA::LightWeightOperatorBase& right)
    {
      return CORE::LINALG::ANA::OperatorProduct(CORE::LINALG::ANA::lw(left), right);
    }
    inline CORE::LINALG::ANA::OperatorSum operator+(
        const CORE::LINALG::ANA::LightWeightOperatorBase& left,
        const CORE::LINALG::ANA::LightWeightOperatorBase& right)
    {
      return CORE::LINALG::ANA::OperatorSum(left, right, 1);
    }
    inline CORE::LINALG::ANA::OperatorSum operator+(
        const Epetra_Operator& left, const Epetra_Operator& right)
    {
      return CORE::LINALG::ANA::OperatorSum(
          CORE::LINALG::ANA::lw(left), CORE::LINALG::ANA::lw(right), 1);
    }
    inline CORE::LINALG::ANA::OperatorSum operator+(
        const CORE::LINALG::ANA::LightWeightOperatorBase& left, const Epetra_Operator& right)
    {
      return CORE::LINALG::ANA::OperatorSum(left, CORE::LINALG::ANA::lw(right), 1);
    }
    inline CORE::LINALG::ANA::OperatorSum operator+(
        const Epetra_Operator& left, const CORE::LINALG::ANA::LightWeightOperatorBase& right)
    {
      return CORE::LINALG::ANA::OperatorSum(CORE::LINALG::ANA::lw(left), right, 1);
    }
    inline CORE::LINALG::ANA::OperatorSum operator-(
        const CORE::LINALG::ANA::LightWeightOperatorBase& left,
        const CORE::LINALG::ANA::LightWeightOperatorBase& right)
    {
      return CORE::LINALG::ANA::OperatorSum(left, right, -1);
    }
    inline CORE::LINALG::ANA::OperatorSum operator-(
        const Epetra_Operator& left, const Epetra_Operator& right)
    {
      return CORE::LINALG::ANA::OperatorSum(
          CORE::LINALG::ANA::lw(left), CORE::LINALG::ANA::lw(right), -1);
    }
    inline CORE::LINALG::ANA::OperatorSum operator-(
        const CORE::LINALG::ANA::LightWeightOperatorBase& left, const Epetra_Operator& right)
    {
      return CORE::LINALG::ANA::OperatorSum(left, CORE::LINALG::ANA::lw(right), -1);
    }
    inline CORE::LINALG::ANA::OperatorSum operator-(
        const Epetra_Operator& left, const CORE::LINALG::ANA::LightWeightOperatorBase& right)
    {
      return CORE::LINALG::ANA::OperatorSum(CORE::LINALG::ANA::lw(left), right, -1);
    }
    // Teuchos::RCP versions of the above operations
    inline CORE::LINALG::ANA::LcOperatorTimesLcsv operator*(
        const CORE::LINALG::ANA::LightWeightOperatorBase& op,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return op * (*vec);
    }
    inline CORE::LINALG::ANA::LcOperatorTimesLcsv operator*(
        const Epetra_Operator& op, const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return op * (*vec);
    }


    /*----------------------------------------------------------------------*
      dot products (result is scalar) (general and specialization versions)
    *----------------------------------------------------------------------*/
    double operator*(const CORE::LINALG::ANA::Vector& vec1, const CORE::LINALG::ANA::Vector& vec2);
    double operator*(const CORE::LINALG::ANA::Vector& vec1, const CORE::LINALG::ANA::LCBase& right);
    double operator*(
        const CORE::LINALG::ANA::Vector& vec1, const CORE::LINALG::ANA::LCSTimesVec& right);
    double operator*(const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::LCBase& right);
    double operator*(
        const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::LCSTimesVec& right);
    inline double operator*(
        const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::Vector& vec1)
    {
      return vec1 * left;
    }
    inline double operator*(
        const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::Vector& vec1)
    {
      return vec1 * left;
    }
    // Teuchos::RCP versions of the above operations
    inline double operator*(const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec1,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec2)
    {
      return (*vec1) * (*vec2);
    }
    inline double operator*(const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec1,
        const CORE::LINALG::ANA::LCBase& right)
    {
      return (*vec1) * right;
    }
    inline double operator*(const CORE::LINALG::ANA::LCBase& left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec1)
    {
      return (*vec1) * left;
    }
    inline double operator*(const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec1,
        const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return (*vec1) * right;
    }
    inline double operator*(const CORE::LINALG::ANA::LCSTimesVec& left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec1)
    {
      return (*vec1) * left;
    }

    /*----------------------------------------------------------------------*
       pointwise multiplications of vectors (result is a vector)
     *----------------------------------------------------------------------*/
    inline CORE::LINALG::ANA::LcVecPointwiseVec pw(
        const CORE::LINALG::ANA::Vector& vec1, const CORE::LINALG::ANA::Vector& vec2)
    {
      return CORE::LINALG::ANA::LcVecPointwiseVec(vec1, vec2);
    }
    inline CORE::LINALG::ANA::LcVecPointwiseLc pw(
        const CORE::LINALG::ANA::Vector& vec, const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcVecPointwiseLc(vec, right);
    }
    inline CORE::LINALG::ANA::LcVecPointwiseLc pw(
        const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::Vector& vec)
    {
      return CORE::LINALG::ANA::pw(vec, left);
    }
    inline CORE::LINALG::ANA::LcVecPointwiseLcsv pw(
        const CORE::LINALG::ANA::Vector& vec, const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::LcVecPointwiseLcsv(vec, right);
    }
    inline CORE::LINALG::ANA::LcVecPointwiseLcsv pw(
        const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::Vector& vec)
    {
      return CORE::LINALG::ANA::pw(vec, left);
    }
    inline CORE::LINALG::ANA::LcLcPointwiseLc pw(
        const CORE::LINALG::ANA::LCBase& left, const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::LcLcPointwiseLc(left, right);
    }
    inline CORE::LINALG::ANA::LcLcsvPointwiseLcsv pw(
        const CORE::LINALG::ANA::LCSTimesVec& left, const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::LcLcsvPointwiseLcsv(left, right);
    }
    // Teuchos::RCP versions of the above operations
    inline CORE::LINALG::ANA::LcVecPointwiseVec pw(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec1,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec2)
    {
      return CORE::LINALG::ANA::pw(*vec1, *vec2);
    }
    inline CORE::LINALG::ANA::LcVecPointwiseLc pw(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec,
        const CORE::LINALG::ANA::LCBase& right)
    {
      return CORE::LINALG::ANA::pw(*vec, right);
    }
    inline CORE::LINALG::ANA::LcVecPointwiseLc pw(const CORE::LINALG::ANA::LCBase& left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return CORE::LINALG::ANA::pw(left, *vec);
    }
    inline CORE::LINALG::ANA::LcVecPointwiseLcsv pw(
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec,
        const CORE::LINALG::ANA::LCSTimesVec& right)
    {
      return CORE::LINALG::ANA::pw(*vec, right);
    }
    inline CORE::LINALG::ANA::LcVecPointwiseLcsv pw(const CORE::LINALG::ANA::LCSTimesVec& left,
        const Teuchos::RCP<const CORE::LINALG::ANA::Vector> vec)
    {
      return CORE::LINALG::ANA::pw(left, *vec);
    }
    /*----------------------------------------------------------------------*
       implicit transpose of an CORE::LINALG::ANA::LightWeightOperatorBase / Epetra_Operator
     *----------------------------------------------------------------------*/
    inline CORE::LINALG::ANA::OperatorTransposed trans(
        const CORE::LINALG::ANA::LightWeightOperatorBase& op)
    {
      return CORE::LINALG::ANA::OperatorTransposed(op);
    }
    inline CORE::LINALG::ANA::OperatorTransposed trans(const Epetra_Operator& op)
    {
      return CORE::LINALG::ANA::OperatorTransposed(CORE::LINALG::ANA::lw(op));
    }
    /*----------------------------------------------------------------------*
       implicit inverse of an operator
     *----------------------------------------------------------------------*/
    inline CORE::LINALG::ANA::OperatorInverse inverse(
        const Epetra_Operator& op, CORE::LINALG::Solver& solver, bool reset = true)
    {
      return CORE::LINALG::ANA::OperatorInverse(op, solver, reset);
    }
    inline CORE::LINALG::ANA::OperatorInverse inverse(
        const CORE::LINALG::SparseOperator& op, CORE::LINALG::Solver& solver, bool reset = true)
    {
      return CORE::LINALG::ANA::OperatorInverse(op, solver, reset);
    }
    inline CORE::LINALG::ANA::OperatorInverse inverse(const Epetra_Operator& op)
    {
      return CORE::LINALG::ANA::OperatorInverse(op);
    }
    inline CORE::LINALG::ANA::OperatorInverse inverse(const CORE::LINALG::SparseOperator& op)
    {
      return CORE::LINALG::ANA::OperatorInverse(op);
    }
    /*----------------------------------------------------------------------*
       norms
     *----------------------------------------------------------------------*/
    inline double norm2(const CORE::LINALG::ANA::Vector& vec)
    {
      double norm;
      vec.Norm2(&norm);
      return norm;
    }
    inline double norm1(const CORE::LINALG::ANA::Vector& vec)
    {
      double norm;
      vec.Norm1(&norm);
      return norm;
    }
    inline double norminf(const CORE::LINALG::ANA::Vector& vec)
    {
      double norm;
      vec.NormInf(&norm);
      return norm;
    }
    inline double norm2(const Teuchos::RCP<const CORE::LINALG::ANA::Vector>& vec)
    {
      return CORE::LINALG::ANA::norm2(*vec);
    }
    inline double norm1(const Teuchos::RCP<const CORE::LINALG::ANA::Vector>& vec)
    {
      return CORE::LINALG::ANA::norm1(*vec);
    }
    inline double norminf(const Teuchos::RCP<const CORE::LINALG::ANA::Vector>& vec)
    {
      return CORE::LINALG::ANA::norminf(*vec);
    }
    double norm2(const CORE::LINALG::ANA::LCBase& lc);
    double norm1(const CORE::LINALG::ANA::LCBase& lc);
    double norminf(const CORE::LINALG::ANA::LCBase& lc);
    inline double norm2(const CORE::LINALG::ANA::LCSTimesVec& lc)
    {
      return lc.Scalar() * norm2(lc.Vector());
    }
    inline double norm1(const CORE::LINALG::ANA::LCSTimesVec& lc)
    {
      return lc.Scalar() * norm1(lc.Vector());
    }
    inline double norminf(const CORE::LINALG::ANA::LCSTimesVec& lc)
    {
      return lc.Scalar() * norminf(lc.Vector());
    }

  }  // namespace ANA
}  // namespace CORE::LINALG



FOUR_C_NAMESPACE_CLOSE

#endif
