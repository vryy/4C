/*!----------------------------------------------------------------------
\file linalg_mlapi_operator.H

\brief Declaration of MLAPI interface to ML

\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
*/
/*----------------------------------------------------------------------*/
#ifndef LINALG_MLAPI_OPERATOR_H
#define LINALG_MLAPI_OPERATOR_H

// Trilinos includes
#include <Epetra_CrsMatrix.h>
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_genP.h"
#include "MLAPI_Error.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"

#include "../drt_lib/drt_dserror.H"


namespace LINALG
{
  /*!
  \brief A multipurpose experimental multigrid operator

  This operator based on the ml advanced programming interface is a
  multipurpose development object for amg ideas that shall be tested in
  the baci framework

  \author gee (gee@lnm.mw.tum.de)
  */
  class AMG_Operator : public virtual Epetra_Operator
  {
   public:
    /*!
    \brief Standard Constructor
    */
    explicit AMG_Operator(
        Teuchos::RCP<Epetra_RowMatrix> A, Teuchos::ParameterList& params, const bool compute);

    /*!
    \brief Destructor
    */
    virtual ~AMG_Operator() {}

    /*!
    \brief Returns the label of this class.
    */
    const char* Label() const { return &label_[0]; }

    /*!
    \brief get Comm of this class

    Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class

    */
    const Epetra_Comm& Comm() const { return (Ainput_->Comm()); }

    /*!
    \brief get reference to parameter list

    */
    Teuchos::ParameterList& Params() const { return params_; }

    /*!
    \brief get input matrix

    */
    Teuchos::RCP<Epetra_CrsMatrix> A() { return Ainput_; }

    /*!
    \brief Get fine level OperatorDomainMap

    Derived from Epetra_Operator, get fine level OperatorDomainMap

    */
    const Epetra_Map& OperatorDomainMap() const { return Ainput_->OperatorDomainMap(); }

    /*!
    \brief Get fine level OperatorRangeMap

    Derived from Epetra_Operator, get fine level OperatorRangeMap

    */
    const Epetra_Map& OperatorRangeMap() const { return Ainput_->OperatorRangeMap(); }

    /*!
    \brief ApplyInverse the preconditioner

     ApplyInverse the preconditioner. Method is derived from Epetra_Operator.


     \param X   (In) : Epetra_MultiVector matching the fine level map of this
                       preconditioner
     \param Y (Out)  : Epetra_MultiVector containing the result on output
    */
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    /*!
    \brief not implemented
    */
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
      dserror("Apply does not make sense for LINALG::AMG_Operator");
      return (-1);
    }

    /*!
    \brief not implemented
    */
    int SetUseTranspose(bool UseTranspose)
    {
      dserror("SetUseTranspose not impl.");
      return -1;
    }

    /*!
    \brief not implemented
    */
    double NormInf() const
    {
      dserror("NormInf not impl.");
      return (-1.0);
    }


    /*!
    \brief not implemented
    */
    bool UseTranspose() const
    {
      dserror("UseTranspose not impl.");
      return false;
    }

    /*!
    \brief not implemented
    */
    bool HasNormInf() const
    {
      dserror("HasNormInf not impl.");
      return false;
    }

   private:
    // don't want copy-ctor and = operator
    AMG_Operator(AMG_Operator& old);
    AMG_Operator operator=(const AMG_Operator& old);

    /*!
    \brief Compute the preconditioner
    */
    void SetupNonSymStab();

    /*!
    \brief Compute row wise 1-norm of a matrix
    */
    MLAPI::MultiVector Row1Norm(MLAPI::Operator& A);

    /*!
    \brief Get number of levels in hierarchy
    */
    inline int NLevel() const { return nlevel_; }

    /*!
    \brief Multigrid vcycle
    */
    void Vcycle(const MLAPI::MultiVector& b_f, MLAPI::MultiVector& x_f, const int level) const;

    /*!
    \brief Get smoother on a certain level
    */
    inline const MLAPI::InverseOperator& S(const int level) const { return S_[level]; }

    /*!
    \brief Get prolongation on a certain level
    */
    inline const MLAPI::Operator& P(const int level) const { return P_[level]; }

    /*!
    \brief Get restruction on a certain level
    */
    inline const MLAPI::Operator& R(const int level) const { return R_[level]; }

    /*!
    \brief Get matrix on a certain level
    */
    inline const MLAPI::Operator& A(const int level) const { return A_[level]; }

    /*!
    \brief Get stabilization matrix on a certain level
    */
    inline const MLAPI::Operator& Astab(const int level) const { return Astab_[level]; }

    /*!
    \brief Get original nonstab matrix on a certain level
    */
    inline const MLAPI::Operator& Anonstab(const int level) const { return Anonstab_[level]; }

    std::string label_;                      // name of this operator
    Teuchos::ParameterList& params_;         // ml parameter list (intput)
    Teuchos::RCP<Epetra_CrsMatrix> Ainput_;  // input matrix
    int nlevel_;                             // no. of grids in hierarchy
    std::vector<MLAPI::Operator> P_;         // the matrix hierarchy
    std::vector<MLAPI::Operator> R_;         // the matrix hierarchy
    std::vector<MLAPI::Operator> A_;         // the matrix hierarchy
    std::vector<MLAPI::Operator> Astab_;     // the matrix hierarchy
    std::vector<MLAPI::Operator> Anonstab_;  // the matrix hierarchy
    std::vector<MLAPI::InverseOperator> S_;  // the matrix hierarchy


  };  // class  AMG_Operator : public virtual Epetra_Operator
}  // namespace LINALG


#endif  // LINALG_MLAPI_OPERATOR_H
