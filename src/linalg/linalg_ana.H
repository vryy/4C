/*!----------------------------------------------------------------------
\file linalg_ana.H

\brief A family of abstract nice algebra operations (ANA)

<pre>
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifndef LINALG_ANA_H
#define LINALG_ANA_H

// Trilinos includes
#include "../drt_lib/drt_dserror.H"
#include "linalg_sparsematrix.H"
#include "linalg_solver.H"

#define DEBUGGING_ANA 0  // turn on to get debugging printouts


namespace LINALG
{
  /*!
  \brief ANA: namespace of the abstract nice algorithm methods and classes

  ANA implements a series of classes, functions and overloaded operators (mainly the latter)
  that allow writing nicer linear algebra code with sparse matrices and distributed vectors.
  Instead of Epetra_Vector and Epetra_CrsMatrix one should use Baci's SparseOperator and
  LINALG::ANA::Vector (which fully implements Epetra_Vector) to write control routines.
  ANA comes with no performance loss and can therefore be used wherever sparse parallel linear
  algebra is needed. It is also good for rapid prototyping of some methods, as it has a very clear
  and easy syntax. Especially for younger collegues, it eliminates the need to tediously learn or
  lookup Epetra syntax. Using the statement \code using namespace LINALG::ANA; \endcode inside a
  function allows for an extremely simplified linear algebra syntax:

  <b>Example 1:</b>

  The user wants to implement \f$ x = x - M^{-1} ( b - A^T * y ) \f$,<br>
  where \f$ x , y , b \f$ are LINALG::ANA:Vector and \f$ M , A \f$ are some<br>
  implementation of Epetra_Operator.

  In classical LINALG / Epetra style this would look like this:
  \code
  Epetra_Vector tmp(A.OperatorDomainMap(),false);
  A.SetUseTranspose(true);
  A.Apply(y,tmp);
  A.SetUseTranspose(false);
  b.Update(-1.0,tmp,1.0);               // note: this ruins b!
  LINALg::Solver solver(M.Comm());
  Epetra_Vector tmp2(M.OperatorDomainMap(),true);
  solver.Solve(M,tmp2,b,true,true);
  x.Update(-1.0,tmp2,1.0);
  \endcode

  In LINALG::ANA this would look like this:
  \code
  using namespace LINALG::ANA;
  x -= inverse(M) * ( b - trans(A) * y ) // note: this does NOT change b!
  \endcode

  If you have, you can incorporate a precomputed LINALG::Solver in here:
  \code
  using namespace LINALG::ANA;
  x -= inverse(M,mysolver,true) * ( b - trans(A) * y ) // note: this does NOT change b!
  \endcode


  <b>Example 2:</b>

  The user wants to implement a vector update as e.g. needed in
  some time integration:<br>
  \f$ d_n = d_n + \Delta d \f$ <br>
  \f$ d_m = ( 1 - \alpha_f ) d_n + \alpha_f * d_o \f$ <br>

  In classical LINALG / Epetra style this would look like this:
  \code
  disn.Update(1.0,disi,1.0);
  dism.Update(1.-alphaf,disn,alphaf,dis,0.0);
  \endcode

  In LINALG::ANA this would look like this:
  \code
  using namespace LINALG::ANA;
  disn += disi;
  dism = ( 1. - alphaf ) * disn + alphaf * dis;
  \endcode

  Here is a short, incomplete list of what LINALG::ANA can do
  (vectors \f$ a,b,c \f$, Operators \f$ A,B,C \f$, scalar \f$s\f$):
  \code
  using namespace LINALG::ANA;
  c += ...                            // summation
  c -= ...                            // subtraction
  c = 0.0;                            // initialization from a scalar
  c *= 2.0;                           // vector scaling
  c /= 2.0;                           // vector scaling
  c = 3.5 * a;                        // vector update
  c = 3.5 * a + 2.4 * b;              // twice vector update
  c = 3.5 * a + 2.4 * b + 2.5 d + ...;// multiple vector update
  c += A*a + B*d - 3.5 * b;           // any combination of linear combinations
  c -= (A + B) * a;                   // implicit Operator sum
  c =  ( A * B ) * a;                 // implicit Operator product
  c = trans(A) * b;                   // implicit transposed Operator
  c = inverse(A) * b;                 // Operator inverse (uses serial Amesos_KLU)
  c = inverse(A,solver,true) * b;     // Operator inverse (with user supplied LINALG::Solver
  instance) c = 2.0 * B * b;                    // implicit Operator scaling s = c * c; // inner
  product (results is scalar) s = 0.5 * c * A * c;                // inner product with linear
  combinations s += 0.5 * c * trans(A+B) * c;      // combine stuff as you wish.... s = norm2(c); //
  norms of vectors (norm2, norm1, norminf) s = norm2(2.0*A*c);                 // norms of any
  linear combinations (norm2, norm1, norminf) c = pw(a,b);                        // pointwise
  multiplication (results is vector) c = pw( a , A * b);                 // pointwise multiplication
  (results is vector) c = pw( 2.0*a + b , A * b);         // pointwise multiplication of generic
  linear combinations \endcode

  <b>Important issues:</b>

  The LINALG::ANA layer on top of algorithms is of cost \f$ o \left( 1 \right) \f$,
  where \f$ n \f$ is the size of the problem. This means, the extra cost is size independent
  and neglectable for \f$ n >> 1 \f$.<br>
  However, one should keep a couple of things in mind (also for performance reasons)
  when writing algorithms with ANA:

  - ANA is an abstract, implicit layer on top of Epetra. All ANA operations are reduced
    to Epetra commands.

  - Though carefully tested, ANA still is a \f$ \beta \f$ - feature and might not be bug-free.
    If you experience strange behavior, tell me.

  - Implicit means, a statement "A * x" or "2.0 * b" or "b * trans(A) * b" does NOT
    compute anything. The only statement that actually computes something is the
    "=" operator. The "=" operation is explicit.

  - There is an explicit "double =" and an explicit "LINALG::ANA::Vector = " operation.
    There intentionally is no explicit "Operator =" operation.
    This means you can not become explicit on matrices.<br>
    The following won't work:
    \code
    A = B;
    A = B+C;
    A = trans(C);
    A = inverse(B);
    \endcode
    but the compiler will tell you so as well....
    Use classical LINALG and Epetra methods (on the same objects) instead.

  - ANA handles rectangular (non-square) Operators correctly.

  - ANA works with scalar double, vectorial LINALG::ANA::Vector and Operator Epetra_Operator.
    This means, you can use it with Epetra_CrsMatrix, Epetra_Operator, SparseMatrix, SparseOperator,
    BlockSparseMatrix etc. You can not use it with Epetra_Vector or Epetra_MultiVector.

  - You can use LINALG::ANA::Vector everywhere instead of Epetra_Vector as it completely implements
    Epetra_Vector including all of its constructors.

  - ANA respects + - before * / and brackets ().

  - Of course, ANA works in parallel.

  - ANA never changes objects on the right hand side of the "=" operator (Yes, they are all 'const'
  ;-))

  - The statements
    \code
    c = 2.0*a + A*b - 3.0*d;
    \endcode
    and
    \code
    c = 2.0*a - 3.0*d + A*b;
    \endcode
    do have the same result, but the latter is slightly faster. The reason is, that ANA tries to
  squeeze as large of a junk of the operation into ONE Epetra call. The first is three Epetra
  methods, the second is two because ANA recognizes \code 2.0*a - 3.0*d \endcode as one Epetra
  operation. ANA can only recognize operations that are directly next to each other to be squeezable
  into one Epetra call. If you are an Epetra expert, you can help ANA by being more explicit: \code
    c = (2.0*a - 3.0*d) + A*b; // ANA would have done this anyway even without brackets
    c = A*b + (2.0*a - 3.0*d); // Here, ANA actually benefits from your help (hard to understand
  why) \endcode A statement like \code c = 2.0*a - (3.0*d - A*b); \endcode would still have the same
  result, but would be slightly slower.

  - ANA never allocates extra memory unless needed by Epetra to perform the operation.<br>
    This does not need extra memory:
    \code
    c = A*b;
    \endcode
    This DOES need extra memory:
    \code
    c += A*b;
    \endcode
    A good rule of thumb is, that if Epetra can do it without extra temporary Vector objects,
    ANA can as well.

  - ANA works with Teuchos::RCPed Vectors on the right hand side of the "=" operator as well.
    It does not work with Teuchos::RCPed Operators and does not accept Teuchos::RCPed quantities on
  the left hand side of the "=" operator: \code Teuchos::RCP<LINALG::ANA::Vector>  a =
  Teuchos::rcp(new LINALG::ANA::Vector(A.RowMap(),false)); Teuchos::RCP<LINALG::ANA::Vector>  b =
  Teuchos::rcp(new LINALG::ANA::Vector(A.RowMap(),false)); Teuchos::RCP<LINALG::ANA::Vector>  c =
  Teuchos::rcp(new LINALG::ANA::Vector(A.RowMap(),false)); Teuchos::RCP<LINALG::SparseMatrix> A =
  Teuchos::rcp(Amatrix)); *a = 3.141;              // Always need to dereference on the left of "="
    *b = 1.0;
    *c = 2.0*a + 3.0*b;      // no need to dereference on the right side of "="
    *c += 2.0 * (*A) * b;    // ANA will not accept Teuchos::RCPed Operator objects -> always need
  to dereference \endcode The compiler will tell you when its ok to use Teuchos::RCPed objects.

    - Of course you are free to mix with classical Epetra calls (not on the same line, though).
      You can also continue using classical Epetra style programming, you can use ANA but you don't
  have to.

    - Obey the important warning below.

  \warning Unless you are an Epetra and ANA expert user, your objects on the left and right side of
  the
           "=" operator should not coincide.
           A construct
           \code
           c += 2.0*c + A*c; // this is NOT deterministic
           \endcode
           might or might NOT lead to the desired result...<br>
           It is very hard to fully understand the order in which operations are done and whether
           they are ok to be performed 'in place' without going through the mines of Moria:
           \code
           c = A*c;            // this is ok, but hard to understand why
           c = 2.0*c;          // this is ok as well.
           c += 2.0*c;         // this is unclear
           c += B*c - 2.0*A*c; // this is definitely going to lead to wrong results without message!
           \endcode
           ANA does not check or safeguard against illegal in-place operations! You will NOT
           receive notification of failure but just wrong results.
           Safeguarding against illegal in-place operations is yet to be implemented.

  <pre>
  \level 1
  \maintainer Martin Kronbichler
              http://www.lnm.mw.tum.de
              089 - 289-15235
  </pre>

  */
  namespace ANA
  {
    // forward declarations
    class LC_s_times_vec;
    class Vector;

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A virtual class implementing an abstract implicit linear combination
           on a vector-valued quantity.

    \author gee (gee@lnm.mw.tum.de)
    */
    class LCBase
    {
     public:
      /// ctor
      LCBase() {}

      /// dtor
      virtual ~LCBase() {}

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
      virtual void Update(LINALG::ANA::Vector& v, const double& scale) const = 0;

      /*!
      \brief Perform " v = scale * " operation

      \param v    (out): A vector with range map this->RangeMap() with values to be set on output
      \param scale (in): a scaling factor for the linear combination (usually -1.0 or 1.0, used for
      sign changes)
      */
      virtual void Set(LINALG::ANA::Vector& v, const double& scale) const = 0;

    };  // class LCBase


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A pure virtual light weight wrapper for a (heavy-weight) Epetra_Operator

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

    \author gee (gee@lnm.mw.tum.de)
    */
    class LightWeightOperatorBase
    {
     public:
      /// ctor
      LightWeightOperatorBase() {}

      /// cctor
      LightWeightOperatorBase(const LightWeightOperatorBase& old) {}

      /// dtor
      virtual ~LightWeightOperatorBase() {}

      /*!
      \brief The derived class shall return a clone of itself by calling its own copy-ctor
      */
      virtual const Teuchos::RCP<LightWeightOperatorBase> Clone() const = 0;

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

    \author gee (gee@lnm.mw.tum.de)
    */
    class LightWeightOperator : public LightWeightOperatorBase
    {
     public:
      LightWeightOperator(const Epetra_Operator& op) : op_(op) {}

      LightWeightOperator(const LightWeightOperator& old)
          : LightWeightOperatorBase(old), op_(old.op_)
      {
      }

      virtual const Teuchos::RCP<LightWeightOperatorBase> Clone() const
      {
        return Teuchos::rcp(new LightWeightOperator(*this));
      }

      virtual ~LightWeightOperator() {}

      virtual int SetUseTranspose(bool UseTranspose)
      {
        return const_cast<Epetra_Operator&>(op_).SetUseTranspose(UseTranspose);
      }

      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {
        return op_.Apply(X, Y);
      }

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {
        return op_.ApplyInverse(X, Y);
      }

      virtual double NormInf() const { return op_.NormInf(); }

      virtual const char* Label() const { return "LINALG::ANA::LightWeightOperator"; }

      virtual bool UseTranspose() const { return op_.UseTranspose(); }

      virtual bool HasNormInf() const { return op_.HasNormInf(); }

      virtual const Epetra_Comm& Comm() const { return op_.Comm(); }

      virtual const Epetra_Map& OperatorDomainMap() const { return op_.OperatorDomainMap(); }

      virtual const Epetra_Map& OperatorRangeMap() const { return op_.OperatorRangeMap(); }

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

    \author gee (gee@lnm.mw.tum.de)
    */
    class OperatorTransposed : public LightWeightOperatorBase
    {
     public:
      OperatorTransposed(const LightWeightOperatorBase& op) : op_(op.Clone()) {}

      OperatorTransposed(const OperatorTransposed& old) : LightWeightOperatorBase(old), op_(old.op_)
      {
      }

      virtual ~OperatorTransposed() {}

      virtual const Teuchos::RCP<LightWeightOperatorBase> Clone() const
      {
        return Teuchos::rcp(new OperatorTransposed(*this));
      }

      virtual int SetUseTranspose(bool UseTranspose)
      {
        // we are transposing the transposed operator
        return const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!UseTranspose);
      }

      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {  // apply the transposed
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        int err = op_->Apply(X, Y);
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        return err;
      }

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        int err = op_->ApplyInverse(X, Y);
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        return err;
      }

      virtual double NormInf() const
      {
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        double out = op_->NormInf();
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        return out;
      }

      virtual const char* Label() const { return "LINALG::ANA::OperatorTransposed"; }

      virtual bool UseTranspose() const { return (!op_->UseTranspose()); }

      virtual bool HasNormInf() const { return op_->HasNormInf(); }

      virtual const Epetra_Comm& Comm() const { return op_->Comm(); }

      virtual const Epetra_Map& OperatorDomainMap() const { return op_->OperatorRangeMap(); }

      virtual const Epetra_Map& OperatorRangeMap() const { return op_->OperatorDomainMap(); }

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

    \author gee (gee@lnm.mw.tum.de)
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

      virtual ~OperatorScaled() {}

      virtual const Teuchos::RCP<LightWeightOperatorBase> Clone() const
      {
        return Teuchos::rcp(new OperatorScaled(*this));
      }

      virtual int SetUseTranspose(bool UseTranspose)
      {
        return const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(UseTranspose);
      }

      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {  // apply the transposed
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        int err = op_->Apply(X, Y);
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        Y.Scale(scalar_);
        return err;
      }

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        int err = op_->ApplyInverse(X, Y);
        const_cast<LightWeightOperatorBase&>(*op_).SetUseTranspose(!op_->UseTranspose());
        Y.Scale(1. / scalar_);
        return err;
      }

      virtual double NormInf() const { return scalar_ * op_->NormInf(); }

      virtual const char* Label() const { return "LINALG::ANA::OperatorScaled"; }

      virtual bool UseTranspose() const { return (op_->UseTranspose()); }

      virtual bool HasNormInf() const { return op_->HasNormInf(); }

      virtual const Epetra_Comm& Comm() const { return op_->Comm(); }

      virtual const Epetra_Map& OperatorDomainMap() const { return op_->OperatorDomainMap(); }

      virtual const Epetra_Map& OperatorRangeMap() const { return op_->OperatorRangeMap(); }

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

    \author gee (gee@lnm.mw.tum.de)
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

      virtual ~OperatorProduct() {}

      virtual const Teuchos::RCP<LightWeightOperatorBase> Clone() const
      {
        return Teuchos::rcp(new OperatorProduct(*this));
      }

      virtual int SetUseTranspose(bool UseTranspose)
      {
        usetransposed_ = UseTranspose;
        return 0;
      }

      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      virtual double NormInf() const
      {
        dserror(
            "LINALG::ANA::OperatorProduct does not implement LightWeightOperatorBase::NormInf()");
        return -1.0;
      }

      virtual const char* Label() const { return "LINALG::ANA::OperatorProduct"; }

      virtual bool UseTranspose() const { return usetransposed_; }

      virtual bool HasNormInf() const { return false; }

      virtual const Epetra_Comm& Comm() const { return left_->Comm(); }

      virtual const Epetra_Map& OperatorDomainMap() const
      {
        if (usetransposed_)
          return left_->OperatorRangeMap();
        else
          return right_->OperatorDomainMap();
      }

      virtual const Epetra_Map& OperatorRangeMap() const
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

    \author gee (gee@lnm.mw.tum.de)
    */
    class OperatorSum : public LightWeightOperatorBase
    {
     public:
      OperatorSum(
          const LightWeightOperatorBase& left, const LightWeightOperatorBase& right, const int sign)
          : sign_(sign), usetransposed_(false), left_(left.Clone()), right_(right.Clone())
      {
        if (sign != 1 && sign != -1) dserror("sign parameter has to be 1 or -1");
      }

      OperatorSum(const OperatorSum& old)
          : LightWeightOperatorBase(old),
            sign_(old.sign_),
            usetransposed_(old.usetransposed_),
            left_(old.left_),
            right_(old.right_)
      {
      }

      virtual ~OperatorSum() {}

      virtual const Teuchos::RCP<LightWeightOperatorBase> Clone() const
      {
        return Teuchos::rcp(new OperatorSum(*this));
      }

      virtual int SetUseTranspose(bool UseTranspose)
      {
        usetransposed_ = UseTranspose;
        return 0;
      }

      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {
        dserror(
            "LINALG::ANA::OperatorSum does not implement LightWeightOperatorBase::ApplyInverse");
        return -1;
      }

      virtual double NormInf() const
      {
        dserror("LINALG::ANA::OperatorSum does not implement LightWeightOperatorBase::NormInf()");
        return -1.0;
      }

      virtual const char* Label() const { return "LINALG::ANA::OperatorSum"; }

      virtual bool UseTranspose() const { return usetransposed_; }

      virtual bool HasNormInf() const { return false; }

      virtual const Epetra_Comm& Comm() const { return left_->Comm(); }

      virtual const Epetra_Map& OperatorDomainMap() const
      {
        if (usetransposed_)
          return left_->OperatorRangeMap();
        else
          return left_->OperatorDomainMap();
      }

      virtual const Epetra_Map& OperatorRangeMap() const
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
    uses a LINALG::Solver with no parameter list that defaults to Amesos_KLU.

    \note Intentionally this class does NOT implement Epetra_Operator though it has all methods
          an Epetra_Operator has as well.

     \sa LightWeightOperatorBase, LINALG::ANA::inverse

    \author gee (gee@lnm.mw.tum.de)
    */
    class OperatorInverse : public LightWeightOperatorBase
    {
     public:
      OperatorInverse(const Epetra_Operator& op, LINALG::Solver& solver, bool reset = true)
          : reset_(reset), solver_(solver), op_(op)
      {
      }

      OperatorInverse(const LINALG::SparseOperator& op, LINALG::Solver& solver, bool reset = true)
          : reset_(reset),
            solver_(solver),
            op_(*(const_cast<LINALG::SparseOperator&>(op).EpetraOperator()))
      {
      }

      OperatorInverse(const Epetra_Operator& op)
          : reset_(true),
            defaultsolver_(Teuchos::rcp(new LINALG::Solver(op.Comm()))),
            solver_(*defaultsolver_),
            op_(op)
      {
      }

      OperatorInverse(const LINALG::SparseOperator& op)
          : reset_(true),
            defaultsolver_(Teuchos::rcp(new LINALG::Solver(op.Comm()))),
            solver_(*defaultsolver_),
            op_(*(const_cast<LINALG::SparseOperator&>(op).EpetraOperator()))
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

      virtual ~OperatorInverse() {}

      virtual const Teuchos::RCP<LightWeightOperatorBase> Clone() const
      {
        return Teuchos::rcp(new OperatorInverse(*this));
      }

      virtual int SetUseTranspose(bool UseTranspose)
      {
        dserror("LINALG::ANA::OperatorInverse does not support transpose");
        return -1;
      }

      virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

      virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {
        dserror(
            "LINALG::ANA::OperatorInverse does not support inverse of inverse of underlying "
            "operator, use Apply instead");
        return -1;
      }

      virtual double NormInf() const
      {
        dserror("LINALG::ANA::OperatorInverse does not support NormInf of inverse of operator");
        return -1.0;
      }

      virtual const char* Label() const { return "LINALG::ANA::OperatorInverse"; }

      virtual bool UseTranspose() const { return false; }

      virtual bool HasNormInf() const { return false; }

      virtual const Epetra_Comm& Comm() const { return op_.Comm(); }

      virtual const Epetra_Map& OperatorRangeMap() const  // no, this is NOT a bug
      {
        return op_.OperatorDomainMap();
      }

      virtual const Epetra_Map& OperatorDomainMap() const  // no, this is NOT a bug
      {
        return op_.OperatorRangeMap();
      }

     private:
      const bool reset_;
      Teuchos::RCP<LINALG::Solver> defaultsolver_;
      LINALG::Solver& solver_;
      const Epetra_Operator& op_;

    };  // class OperatorInverse


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief A distributed vector class that implements Epetra_Vector

    All Epetra_Vector functionality can be used. Additionally, this class overloads a series
    of operators used in ANA linear algebra expressions.

    \sa LightWeightOperatorBase, Epetra_Vector

    \author gee (gee@lnm.mw.tum.de)
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

      /// dtor
      virtual ~Vector() {}

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
      inline void operator=(const LINALG::ANA::Vector& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator = (const LINALG::ANA::Vector& rhs)" << endl;
        fflush(stdout);
#endif
        Update(1.0, rhs, 0.0);
      }

      /// Teuchos::RCP version of the above method
      inline void operator=(const Teuchos::RCP<LINALG::ANA::Vector>& rhs) { *this = *rhs; }

      /*!
      \brief Update this Vector with another Vector

      \param rhs (in): Vector to update this vector with
      */
      inline void operator+=(const LINALG::ANA::Vector& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator += (const LINALG::ANA::Vector& rhs)" << endl;
        fflush(stdout);
#endif
        Update(1.0, rhs, 1.0);
      }

      /// Teuchos::RCP version of the above method
      inline void operator+=(const Teuchos::RCP<LINALG::ANA::Vector>& rhs) { *this += *rhs; }

      /*!
      \brief Update this Vector with negative of another Vector

      \param rhs (in): Vector to update this vector with
      */
      inline void operator-=(const LINALG::ANA::Vector& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator -= (const LINALG::ANA::Vector& rhs)" << endl;
        fflush(stdout);
#endif
        Update(-1.0, rhs, 1.0);
      }

      /// Teuchos::RCP version of the above method
      inline void operator-=(const Teuchos::RCP<LINALG::ANA::Vector>& rhs) { *this -= *rhs; }

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
      inline void operator=(const LINALG::ANA::LCBase& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator = (const LINALG::ANA::LCBase& rhs)" << endl;
        fflush(stdout);
#endif
        rhs.Set(*this, 1.0);
      }

      /*!
      \brief Update this Vector with the result of a linear combination

      \param rhs (in): Linear combination of which the result is used to update this Vector
      */
      inline void operator+=(const LINALG::ANA::LCBase& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator += (LINALG::ANA::LCBase& rhs)" << endl;
        fflush(stdout);
#endif
        rhs.Update(*this, 1.0);
      }

      /*!
      \brief Update this Vector with the negative of the result of a linear combination

      \param rhs (in): Linear combination of which the result is used to negatively update this
      Vector
      */
      inline void operator-=(const LINALG::ANA::LCBase& rhs)
      {
#if DEBUGGING_ANA
        cout << "Vector::operator -= (LINALG::ANA::LCBase& rhs)" << endl;
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

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_lc_plus_lc : public LCBase
    {
     public:
      LC_lc_plus_lc(const LCBase& left, const LCBase& right) : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lc_plus_lc(const LCBase& left, const LCBase& right)" << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_lc_plus_lc()
      {
#if DEBUGGING_ANA
        cout << "~LC_lc_plus_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        left_.Update(v, scale);
        right_.Update(v, scale);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
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

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_lc_minus_lc : public LCBase
    {
     public:
      LC_lc_minus_lc(const LCBase& left, const LCBase& right) : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lc_minus_lc(const LCBase& left, const LCBase& right)" << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_lc_minus_lc()
      {
#if DEBUGGING_ANA
        cout << "~LC_lc_minus_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        left_.Update(v, scale);
        right_.Update(v, -scale);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
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

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_vec_plus_lc : public LCBase
    {
     public:
      LC_vec_plus_lc(const LINALG::ANA::Vector& vec, const LCBase& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_plus_lc(const LINALG::ANA::Vector& vec, const LCBase& right)" << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_vec_plus_lc()
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_plus_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale, vec_, 1.0);
        right_.Update(v, scale);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale, vec_, 0.0);
        right_.Update(v, scale);
      }

     private:
      const LINALG::ANA::Vector& vec_;
      const LCBase& right_;

    };  // class LC_vec_plus_lc

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of a Vector and a generic linear combination

    \sa LCBase

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_vec_minus_lc : public LCBase
    {
     public:
      LC_vec_minus_lc(const LINALG::ANA::Vector& vec, const LCBase& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_minus_lc(const LINALG::ANA::Vector& vec, const LCBase& right)" << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_vec_minus_lc()
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_minus_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale, vec_, 1.0);
        right_.Update(v, -scale);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale, vec_, 0.0);
        right_.Update(v, -scale);
      }

     private:
      const LINALG::ANA::Vector& vec_;
      const LCBase& right_;

    };  // class LC_vec_minus_lc

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of a generic linear combination and a Vector

    \sa LCBase

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_lc_minus_vec : public LCBase
    {
     public:
      LC_lc_minus_vec(const LCBase& left, const LINALG::ANA::Vector& vec)
          : LCBase(), vec_(vec), left_(left)
      {
#if DEBUGGING_ANA
        cout << "LC_lc_minus_vec(const LCBase& left, const LINALG::ANA::Vector& vec)" << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_lc_minus_vec()
      {
#if DEBUGGING_ANA
        cout << "~LC_lc_minus_vec()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        left_.Update(v, scale);
        v.Update(-scale, vec_, 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        left_.Set(v, scale);
        v.Update(-scale, vec_, 1.0);
      }

     private:
      const LINALG::ANA::Vector& vec_;
      const LCBase& left_;

    };  // class LC_lc_minus_vec

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Linear combination of a scalar with a Vector

    \sa LCBase

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_s_times_vec : public LCBase
    {
     public:
      LC_s_times_vec(const double scalar, const LINALG::ANA::Vector& vec)
          : LCBase(), scalar_(scalar), vec_(vec)
      {
#if DEBUGGING_ANA
        cout << "LC_s_times_vec(const double scalar, const LINALG::ANA::Vector& vec)" << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_s_times_vec()
      {
#if DEBUGGING_ANA
        cout << "~LC_s_times_vec() " << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale * scalar_, vec_, 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale * scalar_, vec_, 0.0);
      }

      // give access to the scalar (for specialization LCs)
      inline const double& Scalar() const { return scalar_; }
      // give access to the vector (for specialization LCs)
      inline const LINALG::ANA::Vector& Vector() const { return vec_; }

     private:
      const double scalar_;
      const LINALG::ANA::Vector& vec_;

    };  // class LC_s_times_vec

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Linear combination of a scalar with a generic linear combination

    \sa LCBase

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_s_times_lc : public LCBase
    {
     public:
      LC_s_times_lc(const double scalar, const LINALG::ANA::LCBase& right)
          : LCBase(), scalar_(scalar), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_s_times_lc(const double scalar, const LINALG::ANA::LCBase& right)" << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_s_times_lc()
      {
#if DEBUGGING_ANA
        cout << "~LC_s_times_lc() " << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return right_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        right_.Update(v, scale * scalar_);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        right_.Set(v, scale * scalar_);
      }

     private:
      const double scalar_;
      const LINALG::ANA::LCBase& right_;

    };  // class LC_s_times_lc


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Sum of 2 linear combinations, specialization of LC_lc_plus_lc
           for performance reasons

    \sa LCBase, LC_lc_plus_lc

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_lcsv_plus_lcsv : public LCBase
    {
     public:
      LC_lcsv_plus_lcsv(const LC_s_times_vec& left, const LC_s_times_vec& right)
          : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lcsv_plus_lcsv(const LC_s_times_vec& left, const LC_s_times_vec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_lcsv_plus_lcsv()
      {
#if DEBUGGING_ANA
        cout << "~LC_lcsv_plus_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(
            scale * left_.Scalar(), left_.Vector(), scale * right_.Scalar(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(
            scale * left_.Scalar(), left_.Vector(), scale * right_.Scalar(), right_.Vector(), 0.0);
      }

     private:
      const LC_s_times_vec left_;
      const LC_s_times_vec right_;

    };  // class LC_lcsv_plus_lcsv


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Sum of Vector and a linear combinations, specialization of LC_vec_plus_lc
           for performance reasons

    \sa LCBase, LC_vec_plus_lc

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_vec_plus_lcsv : public LCBase
    {
     public:
      LC_vec_plus_lcsv(const LINALG::ANA::Vector& vec, const LC_s_times_vec& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_plus_lcsv(const LINALG::ANA::Vector& vec, const LC_s_times_vec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_vec_plus_lcsv()
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_plus_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale, vec_, scale * right_.Scalar(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale, vec_, scale * right_.Scalar(), right_.Vector(), 0.0);
      }

     private:
      const LINALG::ANA::Vector& vec_;
      const LC_s_times_vec right_;

    };  // class LC_vec_plus_lcsv

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of 2 linear combinations, specialization of LC_lc_minus_lc
           for performance reasons

    \sa LCBase, LC_lc_minus_lc

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_lcsv_minus_lcsv : public LCBase
    {
     public:
      LC_lcsv_minus_lcsv(const LC_s_times_vec& left, const LC_s_times_vec& right)
          : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lcsv_minus_lcsv(const LC_s_times_vec& left, const LC_s_times_vec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_lcsv_minus_lcsv()
      {
#if DEBUGGING_ANA
        cout << "~LC_lcsv_minus_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(
            scale * left_.Scalar(), left_.Vector(), -scale * right_.Scalar(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(
            scale * left_.Scalar(), left_.Vector(), -scale * right_.Scalar(), right_.Vector(), 0.0);
      }

     private:
      const LC_s_times_vec left_;
      const LC_s_times_vec right_;

    };  // class LC_lcsv_minus_lcsv

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of a Vector and a linear combination, specialization of LC_vec_minus_lc
           for performance reasons

    \sa LCBase, LC_vec_minus_lc

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_vec_minus_lcsv : public LCBase
    {
     public:
      LC_vec_minus_lcsv(const LINALG::ANA::Vector& vec, const LC_s_times_vec& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_minus_lcsv(const LINALG::ANA::Vector& vec, const LC_s_times_vec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_vec_minus_lcsv()
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_minus_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale, vec_, -scale * right_.Scalar(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(scale, vec_, -scale * right_.Scalar(), right_.Vector(), 0.0);
      }

     private:
      const LINALG::ANA::Vector& vec_;
      const LC_s_times_vec right_;

    };  // class LC_vec_minus_lcsv

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Difference of a linear combination and a Vector, specialization of LC_lc_minus_vec
           for performance reasons

    \sa LCBase, LC_lc_minus_vec

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_lcsv_minus_vec : public LCBase
    {
     public:
      LC_lcsv_minus_vec(const LC_s_times_vec& left, const LINALG::ANA::Vector& vec)
          : LCBase(), vec_(vec), left_(left)
      {
#if DEBUGGING_ANA
        cout << "LC_lcsv_minus_vec(const LC_s_times_vec& left, const LINALG::ANA::Vector& vec)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_lcsv_minus_vec()
      {
#if DEBUGGING_ANA
        cout << "~LC_lcsv_minus_vec()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(-scale, vec_, scale * left_.Scalar(), left_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Update(-scale, vec_, scale * left_.Scalar(), left_.Vector(), 0.0);
      }

     private:
      const LINALG::ANA::Vector& vec_;
      const LC_s_times_vec left_;

    };  // class LC_lcsv_minus_vec



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of 2 Vector s, result is a Vector, where
           Result[i] = vec1[i]*vec2[i]

    \sa LCBase

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_vec_pointwise_vec : public LCBase
    {
     public:
      LC_vec_pointwise_vec(const LINALG::ANA::Vector& vec1, const LINALG::ANA::Vector& vec2)
          : LCBase(), vec1_(vec1), vec2_(vec2)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_pointwise_vec(const LINALG::ANA::Vector& vec1, const LINALG::ANA::Vector& "
                "vec2)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_vec_pointwise_vec()
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_pointwise_vec()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return vec1_.Map(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Multiply(scale, vec1_, vec2_, 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Multiply(scale, vec1_, vec2_, 0.0);
      }

     private:
      const LINALG::ANA::Vector& vec1_;
      const LINALG::ANA::Vector& vec2_;

    };  // class LC_vec_pointwise_vec


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of Vector and generic linear combination,
           result is a Vector, where Result[i] = vec[i]*right[i]

    \sa LCBase, LC_vec_pointwise_vec

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_vec_pointwise_lc : public LCBase
    {
     public:
      LC_vec_pointwise_lc(const LINALG::ANA::Vector& vec, const LINALG::ANA::LCBase& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_pointwise_lc(const LINALG::ANA::Vector& vec1, const LINALG::ANA::LCBase& "
                "right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_vec_pointwise_lc()
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_pointwise_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return vec_.Map(); }
      // perform 'v +=' operations
      void Update(LINALG::ANA::Vector& v, const double& scale) const;
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        right_.Set(v, 1.0);
        v.Multiply(scale, vec_, v, 0.0);
      }

     private:
      const LINALG::ANA::Vector& vec_;
      const LINALG::ANA::LCBase& right_;

    };  // class LC_vec_pointwise_lc



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of Vector and linear combination,
           result is a Vector, where Result[i] = vec[i]*right[i].
           Specialization of LC_vec_pointwise_lc for performance reasons

    \sa LCBase, LC_vec_pointwise_vec, LC_vec_pointwise_lc

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_vec_pointwise_lcsv : public LCBase
    {
     public:
      LC_vec_pointwise_lcsv(
          const LINALG::ANA::Vector& vec, const LINALG::ANA::LC_s_times_vec& right)
          : LCBase(), vec_(vec), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_vec_pointwise_lcsv(const LINALG::ANA::Vector& vec, const "
                "LINALG::ANA::LC_s_times_vec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_vec_pointwise_lcsv()
      {
#if DEBUGGING_ANA
        cout << "~LC_vec_pointwise_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return vec_.Map(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Multiply(scale * right_.Scalar(), vec_, right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Multiply(scale * right_.Scalar(), vec_, right_.Vector(), 0.0);
      }

     private:
      const LINALG::ANA::Vector& vec_;
      const LINALG::ANA::LC_s_times_vec right_;

    };  // class LC_vec_pointwise_lc



    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of 2 linear combinations,
           result is a Vector, where Result[i] = left[i]*right[i].

    \sa LCBase, LC_vec_pointwise_vec

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_lc_pointwise_lc : public LCBase
    {
     public:
      LC_lc_pointwise_lc(const LINALG::ANA::LCBase& left, const LINALG::ANA::LCBase& right)
          : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lc_pointwise_lc(const LINALG::ANA::LCBase& left, const LINALG::ANA::LCBase& "
                "right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_lc_pointwise_lc()
      {
#if DEBUGGING_ANA
        cout << "~LC_lc_pointwise_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return left_.RangeMap(); }
      // perform 'v +=' operations
      void Update(LINALG::ANA::Vector& v, const double& scale) const;
      // perform 'v =' operations
      void Set(LINALG::ANA::Vector& v, const double& scale) const;

     private:
      const LINALG::ANA::LCBase& left_;
      const LINALG::ANA::LCBase& right_;

    };  // class LC_vec_pointwise_lc


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Element by element product of 2 linear combinations,
           result is a Vector, where Result[i] = left[i]*right[i].
           Specialization of LC_lc_pointwise_lc for performance reasons

    \sa LCBase, LC_lc_pointwise_lc

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_lcsv_pointwise_lcsv : public LCBase
    {
     public:
      LC_lcsv_pointwise_lcsv(
          const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::LC_s_times_vec& right)
          : LCBase(), left_(left), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_lcsv_pointwise_lcsv(const LINALG::ANA::LC_s_times_vec& left, const "
                "LINALG::ANA::LC_s_times_vec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_lcsv_pointwise_lcsv()
      {
#if DEBUGGING_ANA
        cout << "~LC_lcsv_pointwise_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return left_.RangeMap(); }
      // perform 'v +=' operations
      inline void Update(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Multiply(scale * left_.Scalar() * right_.Scalar(), left_.Vector(), right_.Vector(), 1.0);
      }
      // perform 'v =' operations
      inline void Set(LINALG::ANA::Vector& v, const double& scale) const
      {
        v.Multiply(scale * left_.Scalar() * right_.Scalar(), left_.Vector(), right_.Vector(), 0.0);
      }

     private:
      const LINALG::ANA::LC_s_times_vec left_;
      const LINALG::ANA::LC_s_times_vec right_;

    };  // class LC_lcsv_pointwise_lcsv

    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Apply linear operator to linear combination.
           Specialization of LC_Operator_times_lc for performance reasons.

    \sa LCBase, LC_Operator_times_lc

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_Operator_times_lcsv : public LCBase
    {
     public:
      LC_Operator_times_lcsv(
          const LightWeightOperatorBase& op, const LINALG::ANA::LC_s_times_vec& right)
          : LCBase(), op_(op.Clone()), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_Operator_times_lcsv(const LightWeightOperatorBase& op, const "
                "LINALG::ANA::LC_s_times_vec& right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_Operator_times_lcsv()
      {
#if DEBUGGING_ANA
        cout << "~LC_Operator_times_lcsv()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return op_->OperatorRangeMap(); }
      // perform 'v +=' operations
      void Update(LINALG::ANA::Vector& v, const double& scale) const;
      // perform 'v =' operations
      void Set(LINALG::ANA::Vector& v, const double& scale) const;

     private:
      const Teuchos::RCP<LightWeightOperatorBase> op_;
      const LINALG::ANA::LC_s_times_vec right_;

    };  // class LC_Operator_times_lcsv


    /*----------------------------------------------------------------------*
     *----------------------------------------------------------------------*/
    /*!
    \brief Apply linear operator to linear combination.

    \sa LCBase, LC_Operator_times_lcsv

    \author gee (gee@lnm.mw.tum.de)
    */
    class LC_Operator_times_lc : public LCBase
    {
     public:
      LC_Operator_times_lc(const LightWeightOperatorBase& op, const LINALG::ANA::LCBase& right)
          : LCBase(), op_(op.Clone()), right_(right)
      {
#if DEBUGGING_ANA
        cout << "LC_Operator_times_lc(const LightWeightOperatorBase& op, const "
                "LINALG::ANA::LCBase& right)"
             << endl;
        fflush(stdout);
#endif
      }

      virtual ~LC_Operator_times_lc()
      {
#if DEBUGGING_ANA
        cout << "~LC_Operator_times_lc()" << endl;
        fflush(stdout);
#endif
      }

      // return the range space of the result of the linear combination
      virtual inline const Epetra_BlockMap& RangeMap() const { return op_->OperatorRangeMap(); }
      // perform 'v +=' operations
      void Update(LINALG::ANA::Vector& v, const double& scale) const;
      // perform 'v =' operations
      void Set(LINALG::ANA::Vector& v, const double& scale) const;

     private:
      const Teuchos::RCP<LightWeightOperatorBase> op_;
      const LINALG::ANA::LCBase& right_;

    };  // class LC_Operator_times_lc



    /*----------------------------------------------------------------------*
       static (local) little helper method wrapping an Epetra_Operator
     *----------------------------------------------------------------------*/
    static inline LINALG::ANA::LightWeightOperator lw(const Epetra_Operator& op)
    {
      return LINALG::ANA::LightWeightOperator(op);
    }

    /*----------------------------------------------------------------------*
       scalar, vector and LC operations
     *----------------------------------------------------------------------*/
    inline LINALG::ANA::LC_s_times_vec operator*(
        const double& scalar, const LINALG::ANA::Vector& vec)
    {
      return LINALG::ANA::LC_s_times_vec(scalar, vec);
    }
    inline LINALG::ANA::LC_s_times_vec operator*(
        const LINALG::ANA::Vector& vec, const double& scalar)
    {
      return LINALG::ANA::LC_s_times_vec(scalar, vec);
    }
    inline LINALG::ANA::LC_s_times_vec operator/(
        const LINALG::ANA::Vector& vec, const double& scalar)
    {
      return LINALG::ANA::LC_s_times_vec(1. / scalar, vec);
    }
    inline LINALG::ANA::LC_lc_plus_lc operator+(
        const LINALG::ANA::LCBase& left, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_lc_plus_lc(left, right);
    }
    inline LINALG::ANA::LC_lc_minus_lc operator-(
        const LINALG::ANA::LCBase& left, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_lc_minus_lc(left, right);
    }
    inline LINALG::ANA::LC_vec_plus_lc operator+(
        const LINALG::ANA::Vector& vec, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_vec_plus_lc(vec, right);
    }
    inline LINALG::ANA::LC_vec_plus_lc operator+(
        const LINALG::ANA::LCBase& left, const LINALG::ANA::Vector& vec)
    {
      return LINALG::ANA::LC_vec_plus_lc(vec, left);
    }
    inline LINALG::ANA::LC_vec_minus_lc operator-(
        const LINALG::ANA::Vector& vec, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_vec_minus_lc(vec, right);
    }
    inline LINALG::ANA::LC_lc_minus_vec operator-(
        const LINALG::ANA::LCBase& left, const LINALG::ANA::Vector& vec)
    {
      return LINALG::ANA::LC_lc_minus_vec(left, vec);
    }
    inline LINALG::ANA::LC_lcsv_plus_lcsv operator+(
        const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::LC_lcsv_plus_lcsv(left, right);
    }
    inline LINALG::ANA::LC_lcsv_minus_lcsv operator-(
        const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::LC_lcsv_minus_lcsv(left, right);
    }
    inline LINALG::ANA::LC_vec_plus_lcsv operator+(
        const LINALG::ANA::Vector& vec, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::LC_vec_plus_lcsv(vec, right);
    }
    inline LINALG::ANA::LC_vec_plus_lcsv operator+(
        const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::Vector& vec)
    {
      return LINALG::ANA::LC_vec_plus_lcsv(vec, left);
    }
    inline LINALG::ANA::LC_vec_minus_lcsv operator-(
        const LINALG::ANA::Vector& vec, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::LC_vec_minus_lcsv(vec, right);
    }
    inline LINALG::ANA::LC_lcsv_minus_vec operator-(
        const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::Vector& vec)
    {
      return LINALG::ANA::LC_lcsv_minus_vec(left, vec);
    }
    inline LINALG::ANA::LC_lcsv_plus_lcsv operator+(
        const LINALG::ANA::Vector& left, const LINALG::ANA::Vector& right)
    {
      return (1.0 * left + 1.0 * right);
    }
    inline LINALG::ANA::LC_lcsv_minus_lcsv operator-(
        const LINALG::ANA::Vector& left, const LINALG::ANA::Vector& right)
    {
      return (1.0 * left - 1.0 * right);
    }
    inline LINALG::ANA::LC_s_times_lc operator*(
        const double& scalar, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_s_times_lc(scalar, right);
    }
    inline LINALG::ANA::LC_s_times_lc operator*(
        const LINALG::ANA::LCBase& left, const double& scalar)
    {
      return (scalar * left);
    }
    // Teuchos::RCP versions of the above operations
    inline LINALG::ANA::LC_s_times_vec operator*(
        const double& scalar, const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return (scalar * (*vec));
    }
    inline LINALG::ANA::LC_s_times_vec operator*(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec, const double& scalar)
    {
      return (scalar * (*vec));
    }
    inline LINALG::ANA::LC_s_times_vec operator/(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec, const double& scalar)
    {
      return ((*vec) / scalar);
    }
    inline LINALG::ANA::LC_vec_plus_lc operator+(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec, const LINALG::ANA::LCBase& right)
    {
      return ((*vec) + right);
    }
    inline LINALG::ANA::LC_vec_plus_lc operator+(
        const LINALG::ANA::LCBase& left, const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return (left + (*vec));
    }
    inline LINALG::ANA::LC_vec_minus_lc operator-(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec, const LINALG::ANA::LCBase& right)
    {
      return (*vec - right);
    }
    inline LINALG::ANA::LC_lc_minus_vec operator-(
        const LINALG::ANA::LCBase& left, const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return (left - (*vec));
    }
    inline LINALG::ANA::LC_vec_plus_lcsv operator+(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec, const LINALG::ANA::LC_s_times_vec& right)
    {
      return (*vec + right);
    }
    inline LINALG::ANA::LC_vec_plus_lcsv operator+(
        const LINALG::ANA::LC_s_times_vec& left, const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return (left + (*vec));
    }
    inline LINALG::ANA::LC_vec_minus_lcsv operator-(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec, const LINALG::ANA::LC_s_times_vec& right)
    {
      return (*vec - right);
    }
    inline LINALG::ANA::LC_lcsv_minus_vec operator-(
        const LINALG::ANA::LC_s_times_vec& left, const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return (left - (*vec));
    }
    inline LINALG::ANA::LC_lcsv_plus_lcsv operator+(
        const Teuchos::RCP<const LINALG::ANA::Vector> left,
        const Teuchos::RCP<const LINALG::ANA::Vector> right)
    {
      return ((*left) + (*right));
    }
    inline LINALG::ANA::LC_lcsv_minus_lcsv operator-(
        const Teuchos::RCP<const LINALG::ANA::Vector> left,
        const Teuchos::RCP<const LINALG::ANA::Vector> right)
    {
      return (*left - (*right));
    }


    /*----------------------------------------------------------------------*
       LightWeightOperatorBase and  Epetra_Operator operations
     *----------------------------------------------------------------------*/
    inline LINALG::ANA::LC_Operator_times_lc operator*(
        const LINALG::ANA::LightWeightOperatorBase& op, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_Operator_times_lc(op, right);
    }
    inline LINALG::ANA::LC_Operator_times_lc operator*(
        const Epetra_Operator& op, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_Operator_times_lc(LINALG::ANA::lw(op), right);
    }
    inline LINALG::ANA::LC_Operator_times_lcsv operator*(
        const LINALG::ANA::LightWeightOperatorBase& op, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::LC_Operator_times_lcsv(op, right);
    }
    inline LINALG::ANA::LC_Operator_times_lcsv operator*(
        const LINALG::ANA::LightWeightOperatorBase& op, const LINALG::ANA::Vector& vec)
    {
      return op * (1.0 * vec);
    }
    inline LINALG::ANA::LC_Operator_times_lcsv operator*(
        const Epetra_Operator& op, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::LC_Operator_times_lcsv(LINALG::ANA::lw(op), right);
    }
    inline LINALG::ANA::LC_Operator_times_lcsv operator*(
        const Epetra_Operator& op, const LINALG::ANA::Vector& vec)
    {
      return LINALG::ANA::lw(op) * (1.0 * vec);
    }
    inline LINALG::ANA::OperatorScaled operator*(
        const LINALG::ANA::LightWeightOperatorBase& op, const double& scalar)
    {
      return LINALG::ANA::OperatorScaled(op, scalar);
    }
    inline LINALG::ANA::OperatorScaled operator*(
        const double& scalar, const LINALG::ANA::LightWeightOperatorBase& op)
    {
      return LINALG::ANA::OperatorScaled(op, scalar);
    }
    inline LINALG::ANA::OperatorScaled operator*(const Epetra_Operator& op, const double& scalar)
    {
      return LINALG::ANA::OperatorScaled(LINALG::ANA::lw(op), scalar);
    }
    inline LINALG::ANA::OperatorScaled operator*(const double& scalar, const Epetra_Operator& op)
    {
      return LINALG::ANA::OperatorScaled(LINALG::ANA::lw(op), scalar);
    }
    inline LINALG::ANA::OperatorProduct operator*(const LINALG::ANA::LightWeightOperatorBase& left,
        const LINALG::ANA::LightWeightOperatorBase& right)
    {
      return LINALG::ANA::OperatorProduct(left, right);
    }
    inline LINALG::ANA::OperatorProduct operator*(
        const Epetra_Operator& left, const Epetra_Operator& right)
    {
      return LINALG::ANA::OperatorProduct(LINALG::ANA::lw(left), LINALG::ANA::lw(right));
    }
    inline LINALG::ANA::OperatorProduct operator*(
        const LINALG::ANA::LightWeightOperatorBase& left, const Epetra_Operator& right)
    {
      return LINALG::ANA::OperatorProduct(left, LINALG::ANA::lw(right));
    }
    inline LINALG::ANA::OperatorProduct operator*(
        const Epetra_Operator& left, const LINALG::ANA::LightWeightOperatorBase& right)
    {
      return LINALG::ANA::OperatorProduct(LINALG::ANA::lw(left), right);
    }
    inline LINALG::ANA::OperatorSum operator+(const LINALG::ANA::LightWeightOperatorBase& left,
        const LINALG::ANA::LightWeightOperatorBase& right)
    {
      return LINALG::ANA::OperatorSum(left, right, 1);
    }
    inline LINALG::ANA::OperatorSum operator+(
        const Epetra_Operator& left, const Epetra_Operator& right)
    {
      return LINALG::ANA::OperatorSum(LINALG::ANA::lw(left), LINALG::ANA::lw(right), 1);
    }
    inline LINALG::ANA::OperatorSum operator+(
        const LINALG::ANA::LightWeightOperatorBase& left, const Epetra_Operator& right)
    {
      return LINALG::ANA::OperatorSum(left, LINALG::ANA::lw(right), 1);
    }
    inline LINALG::ANA::OperatorSum operator+(
        const Epetra_Operator& left, const LINALG::ANA::LightWeightOperatorBase& right)
    {
      return LINALG::ANA::OperatorSum(LINALG::ANA::lw(left), right, 1);
    }
    inline LINALG::ANA::OperatorSum operator-(const LINALG::ANA::LightWeightOperatorBase& left,
        const LINALG::ANA::LightWeightOperatorBase& right)
    {
      return LINALG::ANA::OperatorSum(left, right, -1);
    }
    inline LINALG::ANA::OperatorSum operator-(
        const Epetra_Operator& left, const Epetra_Operator& right)
    {
      return LINALG::ANA::OperatorSum(LINALG::ANA::lw(left), LINALG::ANA::lw(right), -1);
    }
    inline LINALG::ANA::OperatorSum operator-(
        const LINALG::ANA::LightWeightOperatorBase& left, const Epetra_Operator& right)
    {
      return LINALG::ANA::OperatorSum(left, LINALG::ANA::lw(right), -1);
    }
    inline LINALG::ANA::OperatorSum operator-(
        const Epetra_Operator& left, const LINALG::ANA::LightWeightOperatorBase& right)
    {
      return LINALG::ANA::OperatorSum(LINALG::ANA::lw(left), right, -1);
    }
    // Teuchos::RCP versions of the above operations
    inline LINALG::ANA::LC_Operator_times_lcsv operator*(
        const LINALG::ANA::LightWeightOperatorBase& op,
        const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return op * (*vec);
    }
    inline LINALG::ANA::LC_Operator_times_lcsv operator*(
        const Epetra_Operator& op, const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return op * (*vec);
    }


    /*----------------------------------------------------------------------*
      dot products (result is scalar) (general and specialization versions)
    *----------------------------------------------------------------------*/
    double operator*(const LINALG::ANA::Vector& vec1, const LINALG::ANA::Vector& vec2);
    double operator*(const LINALG::ANA::Vector& vec1, const LINALG::ANA::LCBase& right);
    double operator*(const LINALG::ANA::Vector& vec1, const LINALG::ANA::LC_s_times_vec& right);
    double operator*(const LINALG::ANA::LCBase& left, const LINALG::ANA::LCBase& right);
    double operator*(
        const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::LC_s_times_vec& right);
    inline double operator*(const LINALG::ANA::LCBase& left, const LINALG::ANA::Vector& vec1)
    {
      return vec1 * left;
    }
    inline double operator*(
        const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::Vector& vec1)
    {
      return vec1 * left;
    }
    // Teuchos::RCP versions of the above operations
    inline double operator*(const Teuchos::RCP<const LINALG::ANA::Vector> vec1,
        const Teuchos::RCP<const LINALG::ANA::Vector> vec2)
    {
      return (*vec1) * (*vec2);
    }
    inline double operator*(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec1, const LINALG::ANA::LCBase& right)
    {
      return (*vec1) * right;
    }
    inline double operator*(
        const LINALG::ANA::LCBase& left, const Teuchos::RCP<const LINALG::ANA::Vector> vec1)
    {
      return (*vec1) * left;
    }
    inline double operator*(const Teuchos::RCP<const LINALG::ANA::Vector> vec1,
        const LINALG::ANA::LC_s_times_vec& right)
    {
      return (*vec1) * right;
    }
    inline double operator*(
        const LINALG::ANA::LC_s_times_vec& left, const Teuchos::RCP<const LINALG::ANA::Vector> vec1)
    {
      return (*vec1) * left;
    }

    /*----------------------------------------------------------------------*
       pointwise multiplications of vectors (result is a vector)
     *----------------------------------------------------------------------*/
    inline LINALG::ANA::LC_vec_pointwise_vec pw(
        const LINALG::ANA::Vector& vec1, const LINALG::ANA::Vector& vec2)
    {
      return LINALG::ANA::LC_vec_pointwise_vec(vec1, vec2);
    }
    inline LINALG::ANA::LC_vec_pointwise_lc pw(
        const LINALG::ANA::Vector& vec, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_vec_pointwise_lc(vec, right);
    }
    inline LINALG::ANA::LC_vec_pointwise_lc pw(
        const LINALG::ANA::LCBase& left, const LINALG::ANA::Vector& vec)
    {
      return LINALG::ANA::pw(vec, left);
    }
    inline LINALG::ANA::LC_vec_pointwise_lcsv pw(
        const LINALG::ANA::Vector& vec, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::LC_vec_pointwise_lcsv(vec, right);
    }
    inline LINALG::ANA::LC_vec_pointwise_lcsv pw(
        const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::Vector& vec)
    {
      return LINALG::ANA::pw(vec, left);
    }
    inline LINALG::ANA::LC_lc_pointwise_lc pw(
        const LINALG::ANA::LCBase& left, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::LC_lc_pointwise_lc(left, right);
    }
    inline LINALG::ANA::LC_lcsv_pointwise_lcsv pw(
        const LINALG::ANA::LC_s_times_vec& left, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::LC_lcsv_pointwise_lcsv(left, right);
    }
    // Teuchos::RCP versions of the above operations
    inline LINALG::ANA::LC_vec_pointwise_vec pw(const Teuchos::RCP<const LINALG::ANA::Vector> vec1,
        const Teuchos::RCP<const LINALG::ANA::Vector> vec2)
    {
      return LINALG::ANA::pw(*vec1, *vec2);
    }
    inline LINALG::ANA::LC_vec_pointwise_lc pw(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec, const LINALG::ANA::LCBase& right)
    {
      return LINALG::ANA::pw(*vec, right);
    }
    inline LINALG::ANA::LC_vec_pointwise_lc pw(
        const LINALG::ANA::LCBase& left, const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return LINALG::ANA::pw(left, *vec);
    }
    inline LINALG::ANA::LC_vec_pointwise_lcsv pw(
        const Teuchos::RCP<const LINALG::ANA::Vector> vec, const LINALG::ANA::LC_s_times_vec& right)
    {
      return LINALG::ANA::pw(*vec, right);
    }
    inline LINALG::ANA::LC_vec_pointwise_lcsv pw(
        const LINALG::ANA::LC_s_times_vec& left, const Teuchos::RCP<const LINALG::ANA::Vector> vec)
    {
      return LINALG::ANA::pw(left, *vec);
    }
    /*----------------------------------------------------------------------*
       implicit transpose of an LINALG::ANA::LightWeightOperatorBase / Epetra_Operator
     *----------------------------------------------------------------------*/
    inline LINALG::ANA::OperatorTransposed trans(const LINALG::ANA::LightWeightOperatorBase& op)
    {
      return LINALG::ANA::OperatorTransposed(op);
    }
    inline LINALG::ANA::OperatorTransposed trans(const Epetra_Operator& op)
    {
      return LINALG::ANA::OperatorTransposed(LINALG::ANA::lw(op));
    }
    /*----------------------------------------------------------------------*
       implicit inverse of an operator
     *----------------------------------------------------------------------*/
    inline LINALG::ANA::OperatorInverse inverse(
        const Epetra_Operator& op, LINALG::Solver& solver, bool reset = true)
    {
      return LINALG::ANA::OperatorInverse(op, solver, reset);
    }
    inline LINALG::ANA::OperatorInverse inverse(
        const LINALG::SparseOperator& op, LINALG::Solver& solver, bool reset = true)
    {
      return LINALG::ANA::OperatorInverse(op, solver, reset);
    }
    inline LINALG::ANA::OperatorInverse inverse(const Epetra_Operator& op)
    {
      return LINALG::ANA::OperatorInverse(op);
    }
    inline LINALG::ANA::OperatorInverse inverse(const LINALG::SparseOperator& op)
    {
      return LINALG::ANA::OperatorInverse(op);
    }
    /*----------------------------------------------------------------------*
       norms
     *----------------------------------------------------------------------*/
    inline double norm2(const LINALG::ANA::Vector& vec)
    {
      double norm;
      vec.Norm2(&norm);
      return norm;
    }
    inline double norm1(const LINALG::ANA::Vector& vec)
    {
      double norm;
      vec.Norm1(&norm);
      return norm;
    }
    inline double norminf(const LINALG::ANA::Vector& vec)
    {
      double norm;
      vec.NormInf(&norm);
      return norm;
    }
    inline double norm2(const Teuchos::RCP<const LINALG::ANA::Vector>& vec)
    {
      return LINALG::ANA::norm2(*vec);
    }
    inline double norm1(const Teuchos::RCP<const LINALG::ANA::Vector>& vec)
    {
      return LINALG::ANA::norm1(*vec);
    }
    inline double norminf(const Teuchos::RCP<const LINALG::ANA::Vector>& vec)
    {
      return LINALG::ANA::norminf(*vec);
    }
    double norm2(const LINALG::ANA::LCBase& lc);
    double norm1(const LINALG::ANA::LCBase& lc);
    double norminf(const LINALG::ANA::LCBase& lc);
    inline double norm2(const LINALG::ANA::LC_s_times_vec& lc)
    {
      return lc.Scalar() * norm2(lc.Vector());
    }
    inline double norm1(const LINALG::ANA::LC_s_times_vec& lc)
    {
      return lc.Scalar() * norm1(lc.Vector());
    }
    inline double norminf(const LINALG::ANA::LC_s_times_vec& lc)
    {
      return lc.Scalar() * norminf(lc.Vector());
    }

  }  // namespace ANA
}  // namespace LINALG



#endif  // LINALG_ANA_H
