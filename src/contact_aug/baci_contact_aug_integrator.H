/*---------------------------------------------------------------------*/
/*! \file
\brief A class to perform integrations of Mortar matrices on the overlap
       of two MORTAR::Elements in 1D and 2D (derived version for
       augmented contact)

\level 2

*/
/*---------------------------------------------------------------------*/

#ifndef BACI_CONTACT_AUG_INTEGRATOR_H
#define BACI_CONTACT_AUG_INTEGRATOR_H

#include "baci_config.H"

#include "baci_contact_aug_integrator_policy.H"
#include "baci_contact_integrator.H"
#include "baci_discretization_fem_general_utils_integration.H"

BACI_NAMESPACE_OPEN

namespace MORTAR
{
  enum ActionType : int;
}

namespace CONTACT
{
  namespace AUG
  {
    // forward declaration
    class IntegratorGeneric;

    /*--------------------------------------------------------------------------*/
    /** \brief Intermediate class necessary to handle a vector of master elements
     *  with arbitrary shapes
     *
     *  \author hiermeier \date 03/17 */
    class IntegrationWrapper : public Integrator
    {
      friend class IntegratorGeneric;

     public:
      /// constructor
      IntegrationWrapper(
          Teuchos::ParameterList& params, CORE::FE::CellType eletype, const Epetra_Comm& comm);

      //! @name 2D augmented Lagrange integration methods
      //! @{
      //! [derived]
      void IntegrateDerivSegment2D(MORTAR::Element& sele, double& sxia, double& sxib,
          MORTAR::Element& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
          const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr) override;

      //! [derived]
      void IntegrateDerivEle2D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
          bool* boundary_ele, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr) override;

      //! @}

      //! @name 3D augmented Lagrange integration methods
      //! @{
      //! [derived]
      void IntegrateDerivCell3DAuxPlane(MORTAR::Element& sele, MORTAR::Element& mele,
          Teuchos::RCP<MORTAR::IntCell> cell, double* auxn, const Epetra_Comm& comm,
          const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr) override;

      //! [derived]
      void IntegrateDerivEle3D(MORTAR::Element& sele, std::vector<MORTAR::Element*> meles,
          bool* boundary_ele, bool* proj_, const Epetra_Comm& comm,
          const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr) override;

      //! @}

      //! Build remaining integrals and linearizations on the 1-D/2-D slave interface
      void IntegrateDerivSlaveElement(MORTAR::Element& sele, const Epetra_Comm& comm,
          const Teuchos::RCP<MORTAR::ParamsInterface>& mparams_ptr);
      void IntegrateDerivSlaveElement(MORTAR::Element& sele, const Epetra_Comm& comm,
          const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

     private:
      IntegratorGeneric* integrator_;

      static INTEGRATOR::UniqueProjInfoPair projInfo_;
    };  // class IntegratorWrapper


    /*--------------------------------------------------------------------------*/
    class IntegratorGeneric
    {
     public:
      static CONTACT::AUG::IntegratorGeneric* Create(int probdim, CORE::FE::CellType slavetype,
          CORE::FE::CellType mastertype, CONTACT::ParamsInterface& cparams,
          CONTACT::Integrator* wrapper);

     private:
      static CONTACT::AUG::IntegratorGeneric* Create2D(CORE::FE::CellType slavetype,
          CORE::FE::CellType mastertype, CONTACT::ParamsInterface& cparams,
          CONTACT::Integrator* wrapper);

      template <CORE::FE::CellType slavetype>
      static CONTACT::AUG::IntegratorGeneric* Create2D(CORE::FE::CellType mastertype,
          CONTACT::ParamsInterface& cparams, CONTACT::Integrator* wrapper);

      template <CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
      static CONTACT::AUG::IntegratorGeneric* Create2D(
          CONTACT::ParamsInterface& cparams, CONTACT::Integrator* wrapper);

      static CONTACT::AUG::IntegratorGeneric* Create3D(CORE::FE::CellType slavetype,
          CORE::FE::CellType mastertype, CONTACT::ParamsInterface& cparams,
          CONTACT::Integrator* wrapper);

      template <CORE::FE::CellType slavetype>
      static CONTACT::AUG::IntegratorGeneric* Create3D(CORE::FE::CellType mastertype,
          CONTACT::ParamsInterface& cparams, CONTACT::Integrator* wrapper);

      template <CORE::FE::CellType slavetype, CORE::FE::CellType mastertype>
      static CONTACT::AUG::IntegratorGeneric* Create3D(
          CONTACT::ParamsInterface& cparams, CONTACT::Integrator* wrapper);

     protected:
      /// empty constructor for singleton construction
      IntegratorGeneric() : cparams_(nullptr), wrapper_(nullptr)
      { /* empty */
      }

      /// initialize the singleton instance
      void Init(CONTACT::ParamsInterface* cparams, CONTACT::Integrator* wrapper)
      {
        if (not cparams or not wrapper)
          dserror(
              "The initialization list is not properly filled. No new"
              "instance can be created.");

        cparams_ = cparams;
        wrapper_ = dynamic_cast<IntegrationWrapper*>(wrapper);
        if (not wrapper_) dserror("dynamic cast failed!");
      }

     public:
      IntegratorGeneric(CONTACT::ParamsInterface& cparams, CONTACT::Integrator& wrapper)
          : cparams_(&cparams), wrapper_(dynamic_cast<IntegrationWrapper*>(&wrapper)){};

      virtual ~IntegratorGeneric() = default;

      //! @name segment based integration methods
      //! @{

      virtual void IntegrateDerivSegment2D(MORTAR::Element& sele, double& sxia, double& sxib,
          MORTAR::Element& mele, double& mxia, double& mxib) = 0;

      //! [derived]
      virtual void IntegrateDerivCell3DAuxPlane(
          MORTAR::Element& sele, MORTAR::Element& mele, MORTAR::IntCell& cell, double* auxn) = 0;
      //! @}

      virtual void Evaluate(MORTAR::Element& sele, MORTAR::Element& mele, bool boundary_ele,
          const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo) = 0;


      //! Build the remaining integrals and linearizations on the 1-D/2-D slave interface
      virtual void IntegrateDerivSlaveElement(MORTAR::Element& sele) = 0;

     protected:
      INPAR::MORTAR::ShapeFcn ShapeFcn() { return Wrapper().shapefcn_; }

      CONTACT::ParamsInterface& CParams()
      {
        if (not cparams_) dserror("cparams_ seems no longer valid!");

        return *cparams_;
      }

      const CONTACT::ParamsInterface& CParams() const
      {
        if (not cparams_) dserror("cparams_ seems no longer valid!");

        return *cparams_;
      }

      IntegrationWrapper& Wrapper()
      {
        if (not wrapper_) dserror("wrapper_ seems no longer valid!");

        return *wrapper_;
      }

     private:
      CONTACT::ParamsInterface* cparams_;

      CONTACT::AUG::IntegrationWrapper* wrapper_;
    };  //

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype,
        class IntPolicy>
    class Integrator : public IntegratorGeneric, public IntPolicy
    {
      // used for static const declarations
      typedef IntPolicy my;

     private:
      /// empty constructor for singleton construction
      Integrator();

     public:
      /** \brief Create a new instance or return a already existing one.
       *
       *  \author hiermeier \date 06/17 */
      static Integrator<probdim, slavetype, mastertype, IntPolicy>* Instance(
          CONTACT::ParamsInterface* cparams, CONTACT::Integrator* wrapper);

     public:
      //! Default constructor without active set information
      Integrator(CONTACT::ParamsInterface& cparams, CONTACT::Integrator& wrapper);

      // don't want = operator and cctor
      Integrator operator=(const Integrator& old) = delete;
      Integrator(const Integrator& old) = delete;

     protected:
      //! @name segment based integration methods
      //! @{

      void IntegrateDerivSegment2D(MORTAR::Element& sele, double& sxia, double& sxib,
          MORTAR::Element& mele, double& mxia, double& mxib) override;

      //! [derived]
      void IntegrateDerivCell3DAuxPlane(MORTAR::Element& sele, MORTAR::Element& mele,
          MORTAR::IntCell& cell, double* auxn) override;

      //! element based evaluation
      void Evaluate(MORTAR::Element& sele, MORTAR::Element& mele, bool boundary_ele,
          const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo) override;

      //! @}

      /**  \brief Build the remaining integrals and linearizations on the 1-D/2-D
       *          slave interface
       *
       *  \author hiermeier */
      void IntegrateDerivSlaveElement(MORTAR::Element& sele) override;

     private:
      class Evaluator;
      class EvaluatorDeriv1stOnly;
      class EvaluatorFull;

      void SetEvaluator(const enum MORTAR::ActionType action);

      //! element based integration
      void IntegrateDerivEle(MORTAR::Element& sele, MORTAR::Element& mele, bool boundary_ele,
          const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo);

      void IntegrateWeightedGap(MORTAR::Element& sele, MORTAR::Element& mele, bool boundary_ele,
          const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo);

      //! integrate the weighted gap gradient error
      void IntegrateWeightedGapGradientError(MORTAR::Element& sele, MORTAR::Element& mele,
          bool boundary_ele, const CONTACT::INTEGRATOR::UniqueProjInfo& projInfo);

      void ExtractActiveSlaveNodeLIDs(
          std::vector<unsigned>& active_nlids, const MORTAR::Element& sele) const;

      int GetLinSize(MORTAR::Element& sele) const;

      /// reset member variables at the beginning of each gauss point loop
      void HardReset(const unsigned linsize);

      void WeakReset(const unsigned linsize);

      //! evaluate the scaling factor kappa at gp
      void GP_kappa(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          double wgt, double jac) const;

      /// \brief evaluate the lin. of the scaling factor kappa at gp (3-D)
      void Get_Deriv1st_Kappa(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, double wgt,
          const Deriv1stMap& d_jac);


      //! evaluate the weighted element Area at gp
      void GP_AugA(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          double wgt, double jac) const;

      //! evaluate the linearization of weighted element area at gp
      void Get_Deriv1st_AugA(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, double wgt, double jac,
          const Deriv1stMap& derivjac) const;

      /// fill unified gauss point normal
      void GP_Normal(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          double* gpn) const;

      /** \brief pure element based variant
       *
       *  Note that it is crucial that we unify the gauss point normal once more
       *  since we neglect the variational term of the smooth normal during the
       *  variation of the normal gauss point gap. Without this normalization
       *  the neglected term wouldn't be equal to zero.
       *
       *  \author hiermeier */
      void GP_Normal_DerivNormal(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval, double* gpn,
          Deriv1stVecMap& dn_non_unit, Deriv2ndVecMap& ddn_non_unit, Deriv1stVecMap& dn_unit,
          Deriv2ndVecMap& ddn_unit);

      //! evaluate the linearization of weighted element area at gp
      void Get_Deriv2nd_Kappa(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv2ndMap& dd_jac) const;

      /** \brief Evaluate the discrete gap at the current GP
       *
       *  \author hiermeier \date 06/17 */
      void GapN(MORTAR::Element& sele, MORTAR::Element& mele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval, const double* gpn,
          double& gapn_sl, double& gapn_ma) const;

      void GP_WGap(MORTAR::Element& sele, const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval,
          const double gapn_sl, const double gapn_ma, const double wg, const double jac) const;

     private:
      // vectors for shape fct. evaluation
      CORE::LINALG::Matrix<my::SLAVENUMNODE, 1> sval_;
      CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE> sderiv_;
      CORE::LINALG::Matrix<my::MASTERNUMNODE, 1> mval_;
      CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE> mderiv_;
      CORE::LINALG::Matrix<3, my::MASTERNUMNODE> mderiv2nd_;
      CORE::LINALG::Matrix<my::SLAVENUMNODE, 1> lmval_;
      CORE::LINALG::Matrix<my::SLAVENUMNODE, my::SLAVEDIM> lmderiv_;

      // create empty vectors for shape fct. evaluation
      CORE::LINALG::Matrix<my::SLAVENUMNODE, 3> ssecderiv_;

      // slave and master nodal coords for Jacobian / GP evaluation
      CORE::LINALG::Matrix<3, my::SLAVENUMNODE> scoord_;
      CORE::LINALG::Matrix<3, my::MASTERNUMNODE> mcoord_;

      // directional derivatives of sxia, sxib, mxia, mxib
      Deriv1stVecMap ximaps_;

      // smooth gauss point normal
      double gpn_[3];

      // 1-st order derivative of the smooth unit gauss point normal
      Deriv1stVecMap dn_unit_;

      // 2-nd order derivative of the smooth unit gauss point normal
      Deriv2ndVecMap ddn_unit_;

      // 1-st order derivative of the smooth non-unit gauss point normal
      Deriv1stVecMap dn_non_unit_;

      // 2-nd order derivative of the smooth non-unit gauss point normal
      Deriv2ndVecMap ddn_non_unit_;

      // GP slave coordinate derivatives
      Deriv1stVecMap dsxigp_;

      // GP master coordinate derivatives (1-st order)
      Deriv1stVecMap dmxigp_;

      // GP auxiliary distance derivatives (1-st order)
      Deriv1stMap dalpha_;

      // GP master coordinate derivatives (2-nd order)
      Deriv2ndVecMap ddmxigp_;

      // Jacobian derivative (1-st order)
      Deriv1stMap derivjac_;

      CORE::GEN::pairedvector<int, CORE::LINALG::Matrix<3, 1>> lingp_;

      // directional derivative of slave GP normal (non-unit)
      Deriv1stMap dmap_nxsl_gp_;
      Deriv1stMap dmap_nysl_gp_;
      Deriv1stMap dmap_nzsl_gp_;

      // 2-nd order derivative of the jacobian matrix
      Deriv2ndMap deriv2ndjac_;

      Deriv1stMap deriv_gapn_sl_;
      Deriv1stMap deriv_gapn_ma_;

      Teuchos::RCP<Evaluator> evaluator_ = Teuchos::null;
    };  // class AugmentedIntegrator

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype,
        class IntPolicy>
    class Integrator<probdim, slavetype, mastertype, IntPolicy>::Evaluator
    {
     protected:
      typedef Integrator<probdim, slavetype, mastertype, IntPolicy> parent_type;
      typedef IntPolicy my;  // parent_type::my my;
     private:
      typedef parent_type::Evaluator class_type;

     public:
      enum class Type
      {
        deriv1st_only,
        full
      };

      Evaluator(parent_type& parent) : parent_(parent){};

      Evaluator(const class_type& source) = delete;

      virtual ~Evaluator() = default;

      virtual Type GetType() const = 0;

      virtual void Deriv_Jacobian(MORTAR::Element& ele, const double* xi,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau) = 0;

      virtual void Deriv_MXiGP(MORTAR::Element& sele, MORTAR::Element& mele, const double* sxi,
          const double* mxi, const double alpha,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, 2>& mtau) = 0;

      virtual void Get_Deriv2nd_AugA(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv2ndMap& dd_jac) const = 0;

     protected:
      parent_type& parent_;
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype,
        class IntPolicy>
    class Integrator<probdim, slavetype, mastertype, IntPolicy>::EvaluatorDeriv1stOnly
        : public Integrator<probdim, slavetype, mastertype, IntPolicy>::Evaluator
    {
     protected:
      typedef Integrator<probdim, slavetype, mastertype, IntPolicy>::Evaluator base_type;
      typedef typename base_type::parent_type parent_type;
      typedef IntPolicy my;  // typename base_type::my my;
     private:
      typedef typename parent_type::EvaluatorDeriv1stOnly class_type;

     public:
      EvaluatorDeriv1stOnly(parent_type& parent) : base_type(parent){};

      EvaluatorDeriv1stOnly(const class_type& source) = delete;

     protected:
      typename base_type::Type GetType() const override { return base_type::Type::deriv1st_only; };

      void Deriv_Jacobian(MORTAR::Element& ele, const double* xi,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau) override;

      void Deriv1st_Jacobian(MORTAR::Element& ele, const double* xi,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau,
          CORE::LINALG::Matrix<probdim, my::SLAVENUMNODE, int>& nodal_dofs,
          CORE::LINALG::Matrix<probdim, 1>& unit_normal, double& length_n_inv,
          Deriv1stVecMap& d_non_unit_normal);

      void Deriv_MXiGP(MORTAR::Element& sele, MORTAR::Element& mele, const double* sxi,
          const double* mxi, const double alpha,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, 2>& mtau) override;

      void Deriv1st_MXiGP(MORTAR::Element& sele, MORTAR::Element& mele, const double* sxi,
          const double* mxi, const double alpha,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, 2>& mtau, CORE::LINALG::Matrix<probdim, probdim>& lmat_inv);

      void Get_Deriv2nd_AugA(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv2ndMap& dd_jac) const override{/* do nothing */};
    };

    /*--------------------------------------------------------------------------*/
    template <unsigned probdim, CORE::FE::CellType slavetype, CORE::FE::CellType mastertype,
        class IntPolicy>
    class Integrator<probdim, slavetype, mastertype, IntPolicy>::EvaluatorFull
        : public Integrator<probdim, slavetype, mastertype, IntPolicy>::EvaluatorDeriv1stOnly
    {
     private:
      typedef Integrator<probdim, slavetype, mastertype, IntPolicy>::EvaluatorDeriv1stOnly
          base_type;
      typedef typename base_type::parent_type parent_type;
      typedef IntPolicy my;  // typename base_type::my my;
     private:
      typedef typename parent_type::EvaluatorFull class_type;

     public:
      EvaluatorFull(parent_type& parent) : base_type(parent){};

      EvaluatorFull(const class_type& source) = delete;

     private:
      typename base_type::Type GetType() const override { return base_type::Type::full; };

      void Deriv_Jacobian(MORTAR::Element& ele, const double* xi,
          const CORE::LINALG::Matrix<my::SLAVEDIM, my::SLAVENUMNODE>& sderiv,
          const CORE::LINALG::Matrix<3, 2>& stau) override;

      void Deriv_MXiGP(MORTAR::Element& sele, MORTAR::Element& mele, const double* sxi,
          const double* mxi, const double alpha,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& sval,
          const CORE::LINALG::Matrix<my::MASTERNUMNODE, 1>& mval,
          const CORE::LINALG::Matrix<my::MASTERDIM, my::MASTERNUMNODE>& mderiv,
          const CORE::LINALG::Matrix<3, 2>& mtau) override;

      void Get_Deriv2nd_AugA(MORTAR::Element& sele,
          const CORE::LINALG::Matrix<my::SLAVENUMNODE, 1>& lmval, const double wgt,
          const Deriv2ndMap& dd_jac) const override;
    };
  }  // namespace AUG
}  // namespace CONTACT


BACI_NAMESPACE_CLOSE

#endif  // CONTACT_AUG_INTEGRATOR_H
