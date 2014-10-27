/*!----------------------------------------------------------------------
\file contact_augmented_integrator.H

<pre>
Created on: Apr 28, 2014

Maintainer: Michael Hiermeier
            hiermeier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15268
</pre>

*----------------------------------------------------------------------*/

#ifndef CONTACT_AUGMENTED_INTEGRATOR_H_
#define CONTACT_AUGMENTED_INTEGRATOR_H_

# include "../drt_contact/contact_integrator.H"
# include "../drt_fem_general/drt_utils_integration.H"

namespace CONTACT
{
class AugmentedIntegrator : public CoIntegrator
{
public:
  /*!
  \brief Constructor with shape function specification

  */
  AugmentedIntegrator(Teuchos::ParameterList& params,
                      DRT::Element::DiscretizationType eletype,
                      const Epetra_Comm& comm,
                      const Teuchos::RCP<Epetra_Map> augActiveSlaveNodes);
  /*!
  \brief Destructor

  */
  virtual ~AugmentedIntegrator() {}

  //! @name 2D augmented Lagrange integration methods

  /*!
  \brief Build integrals on a 1D slave/master overlap
  (segment based integration)

  */
  void IntegrateDerivSegment2D(
       MORTAR::MortarElement& sele, double& sxia, double& sxib,
       MORTAR::MortarElement& mele, double& mxia, double& mxib,
       const Epetra_Comm& comm);

  /*!
  \brief Build integrals on a 1D slave/master overlap
  (element based integration)

  */
  void IntegrateDerivEle2D(
      MORTAR::MortarElement& sele,
      std::vector<MORTAR::MortarElement*> meles,
      bool *boundary_ele);

  /*!
  \brief Build remaining integrals and linearizations on the 1D slave interface

  */
  void IntegrateDerivSlEle2D(
       MORTAR::MortarElement& sele,const Epetra_Comm& comm);

  /*!
  \brief Build integrals on a 2D slave/master overlap
  (segment based integration)

  */
  void IntegrateDerivCell3DAuxPlane(
       MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
       Teuchos::RCP<MORTAR::IntCell> cell, double* auxn,
       const Epetra_Comm& comm);

  /*!
  \brief Build integrals on a 2D slave/master overlap
  (element based integration)

  */
  void IntegrateDerivEle3D(
       MORTAR::MortarElement& sele,
       std::vector<MORTAR::MortarElement*> meles,
       bool *boundary_ele, bool *proj_,
       const Epetra_Comm& comm);

  /*!
  \brief Build the remaining integrals and linearizations on the 2D slave interface

  */
  void IntegrateDerivSlEle3D(
       MORTAR::MortarElement& sele,const Epetra_Comm& comm);

  /*!
  \brief evaluate averaged weighted gap at gp

  */
  void inline GP_2D_wG(
       MORTAR::MortarElement& sele,
       MORTAR::MortarElement& mele,
       LINALG::SerialDenseVector& sval,
       LINALG::SerialDenseVector& mval,
       LINALG::SerialDenseVector& lmval,
       LINALG::SerialDenseMatrix& scoord,
       LINALG::SerialDenseMatrix& mcoord,
       double* awgap,double* gpn,
       double* lengthn,
       double& dxdsxi,double& dsxideta,
       double& wgt);

  /*!
  \brief evaluate the scaling factor kappa at gp

  */
  void inline GP_kappa(
       MORTAR::MortarElement& sele,
       LINALG::SerialDenseVector& lmval,
       double& wgt, double& jac);

  /*!
  \brief evaluate the weighted element Area at gp

  */
  void inline GP_AugA(
      int& it,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      double& wgt,
      double& jac);

  /*!
  \brief evaluate the normal and the derivative of the normal at gp

  */
  void inline GP_Normal_DerivNormal(
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector sval,
      LINALG::SerialDenseMatrix sderiv,
      const std::vector<GEN::pairedvector<int,double> >& dsxigp,
      double* gpn,
      std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
      int& linsize);

  /*!
  \brief evaluate the variation of the weighted gap (inner integral)
  ---> SLAVE and MASTER side

  */
  void inline GP_VarWGap(
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      double* gpn,
      double& wgt, double& jac);

  /*!
  \brief evaluate the lin. of the scaling factor kappa at gp (2-D)
  (segment based integration)

  */
  void inline GP_2D_kappa_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dsxideta, double& dxdsxi,
      double& dxdsxidsxi,
      double& wgt,
      const GEN::pairedvector<int,double>& dsxigp,
      const GEN::pairedvector<int,double>& derivjac,
      const std::vector<GEN::pairedvector<int,double> >& ximaps);

  /*!
  \brief evaluate the lin. of the scaling factor kappa at gp (2-D)
  (Element based integration)

  */
  void inline GP_2D_kappa_Ele_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dxdsxi,double& wgt,
      const GEN::pairedvector<int,double>& derivjac);

  /*!
  \brief evaluate the lin. of the scaling factor kappa at gp (3-D)

  */
  void inline GP_3D_kappa_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseMatrix& lmderiv,
      double& wgt, double& jac,
      const std::vector<GEN::pairedvector<int,double> >& dsxigp,
      const GEN::pairedvector<int,double>& jacintcellmap);

  /*!
  \brief evaluate the linearization of weighted element area at gp

  */
  void inline GP_AugA_Lin(
      int& it,
      MORTAR::MortarElement& sele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      double& wgt,
      double& jac,
      const GEN::pairedvector<int,double>& derivjac);

  /*!
  \brief evaluate linearization of the variation of the weighted gap (2-D)
  ---> SLAVE and MASTER side

  */
  void inline GP_2D_VarWGap_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      double* gpn,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dsxideta, double& dxdsxi,
      double& dxdsxidsxi,
      double& wgt,
      const GEN::pairedvector<int,double>& dsxigp,
      const GEN::pairedvector<int,double>& dmxigp,
      const GEN::pairedvector<int,double>& derivjac,
      const std::vector<GEN::pairedvector<int,double> >& dnmap_unit,
      const std::vector<GEN::pairedvector<int,double> >& ximaps);

  /*!
  \brief evaluate linearization of the variation of the weighted gap (2-D)
  ---> SLAVE and MASTER side
  (element based integration)

  */
  void inline GP_2D_VarWGap_Ele_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      double* gpn,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& dxdsxi, double& wgt,
      const GEN::pairedvector<int,double>& dmxigp,
      const GEN::pairedvector<int,double>&derivjac,
      const std::vector<GEN::pairedvector<int,double> >& dnmap_unit);

  /*!
  \brief evaluate linearization of the variation of the weighted gap (3-D)
  ---> SLAVE and MASTER side

  */
  void inline GP_3D_VarWGap_Lin(
      int& iter,
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseVector& lmval,
      double* gpn,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      double& wgt, double& jac,
      const std::vector<GEN::pairedvector<int,double> >& dsxigp,
      const std::vector<GEN::pairedvector<int,double> >& dmxigp,
      const GEN::pairedvector<int,double>& jacintcellmap,
      const std::vector<GEN::pairedvector<int,double> >& dnmap_unit);
  //@}

  //! @name Return methods

  /*!
  \brief Return the solution type

  */
  INPAR::CONTACT::SolvingStrategy SolType()
     { return DRT::INPUT::IntegralValue<INPAR::CONTACT::SolvingStrategy>(imortar_,"STRATEGY"); }

  //@}

protected:
  // don't want = operator and cctor
  AugmentedIntegrator operator = (const AugmentedIntegrator& old);
  AugmentedIntegrator(const AugmentedIntegrator& old);

  const Teuchos::RCP<Epetra_Map> augActiveSlaveNodes_;            // Map of all active augmented slave nodes of the current interface
};
}

#endif /* CONTACT_AUGMENTED_INTEGRATOR_H_ */
