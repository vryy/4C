/*---------------------------------------------------------------------------*/
/*!
\brief Integrator class for the inequality constraint level-set approach
       (originally created by student Michael Hofer)

\level 3

\maintainer Matthias Mayr
*/
/*---------------------------------------------------------------------------*/

#ifndef XCONTACT_INTEGRATOR_H
#define XCONTACT_INTEGRATOR_H

#include "../drt_contact/contact_integrator.H"

namespace XCONTACT
{
  // forward declarations
  template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
      unsigned numNodeEle>
  class GenericEvaluator;

  template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
      unsigned numNodeEle>
  class GapEvaluator;

  template <unsigned probdim, DRT::Element::DiscretizationType eletype, bool wgap_only,
      unsigned dim, unsigned numNodeEle>
  class ForceEvaluator;

  template <unsigned probdim, DRT::Element::DiscretizationType eletype, unsigned dim,
      unsigned numNodeEle>
  class ForceTangentEvaluator;

  /** \brief Extended contact integration class
   *
   *  Here you can find the actual calculation of the different contact terms.
   *
   *  Information to the template definition:
   *
   *  probdim:    Dimension of domain for current problem (2D or 3D)
   *  eletype:    Discretization type of slave element
   *  dim:        Dimension of boundary for current problem (1D or 2D)
   *  numNodeEle: Number of slave element nodes
   *
   *  \remark Currently the template considers only the slave element, assuming that
   *          all master elements have the same element type */
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  class Integrator : public CONTACT::CoIntegrator
  {
    // friend classes
    friend class ForceEvaluator<probdim, eletype, false, dim, numNodeEle>;
    friend class ForceEvaluator<probdim, eletype, true, dim, numNodeEle>;
    friend class ForceTangentEvaluator<probdim, eletype, dim, numNodeEle>;
    friend class GapEvaluator<probdim, eletype, dim, numNodeEle>;

   public:
    //! @name Constructors and destructors and related methods
    //! @{

    //! Constructor
    Integrator(Teuchos::ParameterList& params, const Epetra_Comm& comm);

    //! Destructor
    virtual ~Integrator(){};

    //! @}

    //! @name Derived functions
    //! @{

    //! @name currently unsupported derived methods
    //! @{

    virtual void IntegrateDerivSegment2D(MORTAR::MortarElement& sele, double& sxia, double& sxib,
        MORTAR::MortarElement& mele, double& mxia, double& mxib, const Epetra_Comm& comm,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
    {
      dserror("Segment based integration is currently unsupported!");
    }

    virtual void IntegrateDerivEle3D(MORTAR::MortarElement& sele,
        std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele, bool* proj_,
        const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
    {
      dserror("Element based integration in 3D is currently unsupported!");
    }

    virtual void IntegrateDerivCell3DAuxPlane(MORTAR::MortarElement& sele,
        MORTAR::MortarElement& mele, Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
        const Epetra_Comm& comm, const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr)
    {
      dserror("The auxiliary plane 3-D coupling integration case is currently unsupported!");
    }

    //! @}

    /*! Compute contact residual and tangent matrix with H1 duality pairing
     *  for 2D case using element based integration on slave element and
     *  closest point projection
     *
     *  \author hofer
     *  \date 07/16 */
    virtual void IntegrateDerivEle2D(MORTAR::MortarElement& sele,
        std::vector<MORTAR::MortarElement*> meles, bool* boundary_ele,
        const Teuchos::RCP<CONTACT::ParamsInterface>& cparams_ptr);

    //! @}

   private:
    //! @name Internal functions
    //! @{

    /// Build the desired evaluator object
    Teuchos::RCP<const XCONTACT::GenericEvaluator<probdim, eletype, dim, numNodeEle>>
    BuildEvaluator(const CONTACT::ParamsInterface& cparams) const;

    //! @}
   private:
#if __cplusplus < 201103L
    /// get the number of nodes per slave element
    unsigned SlNumNodeEle() const { return numNodeEle; }

    /// get the number of nodes per master element
    unsigned MaNumNodeEle() const { return numNodeEle; }
#else
    /// get the number of nodes per slave element
    constexpr unsigned SlNumNodeEle() const { return numNodeEle; }

    /// get the number of nodes per master element
    constexpr unsigned MaNumNodeEle() const { return numNodeEle; }
#endif

    const bool& IsConstantNormal() const { return is_const_normal_; };

    const bool& IsL2VarJacobi() const { return is_l2_var_jacobi_; }

    const bool& IsH1() const { return is_h1_; }

   private:
    //! @name Member variables
    //! @{

    bool is_const_normal_;

    //
    bool is_l2_var_jacobi_;

    //
    bool is_h1_;

    //! @}

  };  // class Integrator

  // abstract helper class
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  class GenericEvaluator
  {
   public:
    // constructor
    GenericEvaluator(){/* intentionally left blank */};

    // destructor
    virtual ~GenericEvaluator(){};

    /** \brief Abstract evaluate routine
     *
     *  \author hiermeier \date 11/16 */
    virtual void Evaluate(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const LINALG::Matrix<1, numNodeEle>& sval, const LINALG::Matrix<1, numNodeEle>& mval,
        const LINALG::Matrix<1, numNodeEle>& lmval, const LINALG::Matrix<dim, numNodeEle>& sderiv,
        const LINALG::Matrix<dim, numNodeEle>& mderiv,
        const LINALG::Matrix<dim, numNodeEle>& lmderiv, const LINALG::Matrix<dim, 1>& sxi,
        const LINALG::Matrix<dim, 1>& mxi, const LINALG::Matrix<probdim, 1>& sx,
        const LINALG::Matrix<probdim, 1>& mx, const LINALG::Matrix<probdim, probdim>& sx_sx,
        const LINALG::Matrix<dim, probdim>& mg, const LINALG::Matrix<probdim, 1>& mn,
        const double& gN, const LINALG::Matrix<1, dim>& gN_sxi, const double& lm,
        const LINALG::Matrix<1, dim>& lm_sxi, const double& sj, const double& sjinv,
        const double& sjc, const double& binv, const double& wgt,
        const XCONTACT::Integrator<probdim, eletype, dim, numNodeEle>& integrator) const = 0;

    /*! \brief Project spatial point from slave element onto master element
     *  via closest point projection
     *
     *  \author Hofer
     *  \date 08/16 */
    void ClosestPointProjection(const LINALG::Matrix<probdim, 1>& sx,
        const LINALG::Matrix<dim, probdim>& sg, MORTAR::MortarElement& mele,
        LINALG::Matrix<dim, 1>& mxi, LINALG::Matrix<1, numNodeEle>& mval,
        LINALG::Matrix<dim, numNodeEle>& mderiv, LINALG::Matrix<probdim, 1>& mx,
        LINALG::Matrix<dim, probdim>& mg, double& binv) const;

    /*! \brief Evaluate gap in normal direction and parametric derivatives based on closest point
     *  projection
     *
     *  \author hofer
     *  \date 08/16 */
    void NormalGapCPP(const LINALG::Matrix<probdim, 1>& sx, const LINALG::Matrix<probdim, 1>& mx,
        const LINALG::Matrix<probdim, 1>& mn, const LINALG::Matrix<dim, probdim>& sg, double& gN,
        LINALG::Matrix<1, dim>& gN_sxi) const;
  };  // class GenericEvaluator

  // force evaluator class object
  template <unsigned probdim, DRT::Element::DiscretizationType eletype, bool wgap_only = false,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  class ForceEvaluator : public GenericEvaluator<probdim, eletype, dim, numNodeEle>
  {
   public:
    // constructor
    ForceEvaluator(){/* intentionally left blank */};

    // destructor
    virtual ~ForceEvaluator(){};

    /** \brief Add contributions of Gauss point to contact residual
     *
     *  \author hofer \date 11/16 */
    virtual void Evaluate(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const LINALG::Matrix<1, numNodeEle>& sval, const LINALG::Matrix<1, numNodeEle>& mval,
        const LINALG::Matrix<1, numNodeEle>& lmval, const LINALG::Matrix<dim, numNodeEle>& sderiv,
        const LINALG::Matrix<dim, numNodeEle>& mderiv,
        const LINALG::Matrix<dim, numNodeEle>& lmderiv, const LINALG::Matrix<dim, 1>& sxi,
        const LINALG::Matrix<dim, 1>& mxi, const LINALG::Matrix<probdim, 1>& sx,
        const LINALG::Matrix<probdim, 1>& mx, const LINALG::Matrix<probdim, probdim>& sx_sx,
        const LINALG::Matrix<dim, probdim>& mg, const LINALG::Matrix<probdim, 1>& mn,
        const double& gN, const LINALG::Matrix<1, dim>& gN_sxi, const double& lm,
        const LINALG::Matrix<1, dim>& lm_sxi, const double& sj, const double& sjinv,
        const double& sjc, const double& binv, const double& wgt,
        const XCONTACT::Integrator<probdim, eletype, dim, numNodeEle>& integrator) const;
  };  // class ForceEvaluator

  // gap evaluator class object
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  class GapEvaluator : public ForceEvaluator<probdim, eletype, true, dim, numNodeEle>
  {
   public:
    // constructor
    GapEvaluator(){/* intentionally left blank */};

    // destructor
    virtual ~GapEvaluator(){};

  };  // class GapEvaluator

  // force evaluator class object
  template <unsigned probdim, DRT::Element::DiscretizationType eletype,
      unsigned dim = DRT::UTILS::DisTypeToDim<eletype>::dim,
      unsigned numNodeEle = DRT::UTILS::DisTypeToNumNodePerEle<eletype>::numNodePerElement>
  class ForceTangentEvaluator : public ForceEvaluator<probdim, eletype, false, dim, numNodeEle>
  {
   public:
    // constructor
    ForceTangentEvaluator(){/* intentionally left blank */};

    // destructor
    virtual ~ForceTangentEvaluator(){};

    /** \brief Add contributions of Gauss point to contact residual and
     *  tangential matrix
     *
     *  \author hofer \date 11/16 */
    virtual void Evaluate(MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const LINALG::Matrix<1, numNodeEle>& sval, const LINALG::Matrix<1, numNodeEle>& mval,
        const LINALG::Matrix<1, numNodeEle>& lmval, const LINALG::Matrix<dim, numNodeEle>& sderiv,
        const LINALG::Matrix<dim, numNodeEle>& mderiv,
        const LINALG::Matrix<dim, numNodeEle>& lmderiv, const LINALG::Matrix<dim, 1>& sxi,
        const LINALG::Matrix<dim, 1>& mxi, const LINALG::Matrix<probdim, 1>& sx,
        const LINALG::Matrix<probdim, 1>& mx, const LINALG::Matrix<probdim, probdim>& sx_sx,
        const LINALG::Matrix<dim, probdim>& mg, const LINALG::Matrix<probdim, 1>& mn,
        const double& gN, const LINALG::Matrix<1, dim>& gN_sxi, const double& lm,
        const LINALG::Matrix<1, dim>& lm_sxi, const double& sj, const double& sjinv,
        const double& sjc, const double& binv, const double& wgt,
        const XCONTACT::Integrator<probdim, eletype, dim, numNodeEle>& integrator) const;
  };  // class ForceTangentEvaluator

}  // namespace XCONTACT

#endif  // #ifndef XCONTACT_INTEGRATOR_H
