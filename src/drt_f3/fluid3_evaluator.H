#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef FLUID3_EVALUATOR_H
#define FLUID3_EVALUATOR_H

#include <string>

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_elementgroup.H"

#include "fluid3.H"

#include "../drt_mat/newtonianfluid.H"


namespace DRT
{
  namespace ELEMENTS
  {

    /// Fluid3 specific evaluator
    /*!

      The system matrix and RHS construction of Fluid3 elements.

      \author u.kue
      \date 03/08
     */
    class Fluid3SystemEvaluator : public DRT::EGROUP::SystemEvaluatorBase
    {
    public:
      Fluid3SystemEvaluator(Teuchos::RCP<DRT::Discretization> dis,
                            const Teuchos::ParameterList& params,
                            Teuchos::RCP<LINALG::SparseOperator> systemmatrix,
                            Teuchos::RCP<Epetra_Vector> systemvector);

      virtual void Evaluate(DRT::EGROUP::OtherElements& elements);

      virtual void Evaluate(DRT::EGROUP::AlignedHex8& elements);

    private:

      /// Evaluate one element
      void ElementEvaluation(DRT::Element* ele,
                             Epetra_SerialDenseMatrix& elematrix1,
                             Epetra_SerialDenseVector& elevector1,
                             std::vector<int>& lm);

      /// Simple element by element evaluation
      void PlainEvaluate(DRT::EGROUP::Group& elements);

      /// Evaluate elements in block to keep element data in cache
      void BlockEvaluate(DRT::EGROUP::Group& elements, int blocksize=20);

      Fluid3::StabilisationAction ConvertStringToStabAction(const string& action) const;

      /** \name Smagorinsky variables */
      //@{

      Teuchos::RCP<const Epetra_Vector> velnp_;
      Teuchos::RCP<const Epetra_Vector> hist_;
      Teuchos::RCP<const Epetra_Vector> dispnp_;
      Teuchos::RCP<const Epetra_Vector> gridv_;
      Teuchos::RCP<const Epetra_Vector> fsvelnp_;

      //@}

      /// fine scale something
      Fluid3::StabilisationAction fssgv_;

      /** \name Smagorinsky variables */
      //@{

      Fluid3::TurbModelAction turb_mod_action_;

      double Cs_;
      double Cs_delta_sq_;
      double l_tau_;
      double visceff_;

      //@}

      /// map that converts strings to actions (stabilisation)
      static std::map<std::string,Fluid3::StabilisationAction> stabstrtoact_;

      MATERIAL* actmat_;

      double time_;
      double timefac_;

      bool newton_;

      Fluid3::StabilisationAction pspg_;
      Fluid3::StabilisationAction supg_;
      Fluid3::StabilisationAction vstab_;
      Fluid3::StabilisationAction cstab_;
      Fluid3::StabilisationAction cross_;
      Fluid3::StabilisationAction reynolds_;

      static Teuchos::RCP<Teuchos::Time> otherelementtime_;
      static Teuchos::RCP<Teuchos::Time> alignedhex8time_;
    };
  }
}

#endif
#endif
#endif
