/*---------------------------------------------------------------------*/
/*!
\file contact_nitsche_integrator.H

\brief A class to perform integrations of nitsche related terms

\level 3

\maintainer Alexander Seitz

*/
/*---------------------------------------------------------------------*/
#ifndef CONTACT_INTEGRATOR_NITSCHE_H
#define CONTACT_INTEGRATOR_NITSCHE_H

#include "contact_integrator.H"
#include "../headers/pairedvector.H"
#include <Epetra_FEVector.h>
#include <Epetra_CrsMatrix.h>

// forward declarations
namespace LINALG{class SerialDenseVector;}
/*!
 \brief CONTACT: namespace of the BACI contact and meshtying algorithms

 */
namespace CONTACT
{
template<int dim>
 struct DimToNumStr {};
 template<> struct DimToNumStr<1>     {static const int NumStr =  1;};
 template<> struct DimToNumStr<2>     {static const int NumStr =  3;};
 template<> struct DimToNumStr<3>     {static const int NumStr =  6;};

class CoIntegratorNitsche : public CONTACT::CoIntegrator
{
public:

  /*!
   \brief Constructor  with shape function specification

   Constructs an instance of this class using a specific type of shape functions.<br>
   Note that this is \b not a collective call as overlaps are
   integrated in parallel by individual processes.<br>
   Note also that this constructor relies heavily on the
   DRT::UTILS::IntegrationPoints structs to get Gauss points
   and corresponding weights.

   */
   CoIntegratorNitsche(Teuchos::ParameterList& params,
      DRT::Element::DiscretizationType eletype, const Epetra_Comm& comm):
        CoIntegrator(params,eletype,comm),
        fc_(params.get<Teuchos::RCP<Epetra_FEVector> >("GLOBAL_CONTACT_FORCE",Teuchos::null)),
        kc_(params.get<Teuchos::RCP<LINALG::SparseMatrix> >("GLOBAL_CONTACT_STIFF",Teuchos::null)),
        theta_(params.get<double>("NITSCHE_THETA")),
        nit_wgt_(DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheWeighting>(params,"NITSCHE_WEIGHTING")){}

  /*!
   \brief Destructor

   */
  virtual ~CoIntegratorNitsche()
  {
  }

  /*!
   \brief Perform integration at GP
          This is where the distinction between methods should be,
          i.e. mortar, augmented, gpts,...
   */
  virtual void IntegrateGP_2D(
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      GEN::pairedvector<int,Epetra_SerialDenseMatrix>& dualmap,
      double& wgt,
      double& jac,
      GEN::pairedvector<int, double>& derivjac,
      double* normal,
      std::vector<GEN::pairedvector<int, double> >& dnmap_unit,
      double& gap,
      GEN::pairedvector<int, double>& deriv_gap,
      double* sxi,
      double* mxi,
      std::vector<GEN::pairedvector<int,double> >& derivsxi,
      std::vector<GEN::pairedvector<int,double> >& derivmxi
      );

  /*!
   \brief Perform integration at GP
          This is where the distinction between methods should be,
          i.e. mortar, augmented, gpts,...
   */
  void IntegrateGP_3D(
      MORTAR::MortarElement& sele,
      MORTAR::MortarElement& mele,
      LINALG::SerialDenseVector& sval,
      LINALG::SerialDenseVector& lmval,
      LINALG::SerialDenseVector& mval,
      LINALG::SerialDenseMatrix& sderiv,
      LINALG::SerialDenseMatrix& mderiv,
      LINALG::SerialDenseMatrix& lmderiv,
      GEN::pairedvector<int,Epetra_SerialDenseMatrix>& dualmap,
      double& wgt,
      double& jac,
      GEN::pairedvector<int, double>& derivjac,
      double* normal,
      std::vector<GEN::pairedvector<int, double> >& dnmap_unit,
      double& gap,
      GEN::pairedvector<int, double>& deriv_gap,
      double* sxi,
      double* mxi,
      std::vector<GEN::pairedvector<int,double> >& derivsxi,
      std::vector<GEN::pairedvector<int,double> >& derivmxi
      );


  /*!
  \brief evaluate GPTS forces and linearization at this gp
  */
    template <int dim>
    void GPTS_forces(
        MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const LINALG::SerialDenseVector& sval, const LINALG::SerialDenseMatrix& sderiv,const std::vector<GEN::pairedvector<int,double> >& dsxi,
        const LINALG::SerialDenseVector& mval, const LINALG::SerialDenseMatrix& mderiv,const std::vector<GEN::pairedvector<int,double> >& dmxi,
        const double jac,const GEN::pairedvector<int,double>& jacintcellmap, const double wgt,
        const double gap, const GEN::pairedvector<int,double>& dgapgp,
        double* gpn, std::vector<GEN::pairedvector<int,double> >& dnmap_unit,double* sxi, double* mxi);

    template <int dim>
    void IntegrateNormalTest(
        const double fac,
        MORTAR::MortarElement& sele, MORTAR::MortarElement& mele,
        const LINALG::SerialDenseVector& sval, const LINALG::SerialDenseMatrix& sderiv,const std::vector<GEN::pairedvector<int,double> >& dsxi,
        const LINALG::SerialDenseVector& mval, const LINALG::SerialDenseMatrix& mderiv,const std::vector<GEN::pairedvector<int,double> >& dmxi,
        const double jac,const GEN::pairedvector<int,double>& jacintcellmap, const double wgt,
        const double gap, const GEN::pairedvector<int,double>& dgapgp,
        double* gpn, std::vector<GEN::pairedvector<int,double> >& dnmap_unit);

    template<int dim>
        void BuildNormalAdjointTest(
        MORTAR::MortarElement& moEle,
        const double fac,
        const Epetra_SerialDenseMatrix& dsnndd,
        const Epetra_SerialDenseMatrix& d2snndd2,
        const Epetra_SerialDenseMatrix& d2snnDdDn,
        const Epetra_SerialDenseMatrix& d2snnDdDpxi,
        const std::vector<GEN::pairedvector<int,double> > boundary_gpcoord_lin,
        LINALG::Matrix<dim,dim> derivtravo_slave,
        const std::vector<GEN::pairedvector<int,double> >& contactN_deriv,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test);

    template<int dim>
    void IntegrateNormalAdjointTest(
        const double fac,
        const double jac, const GEN::pairedvector<int,double>& jacintcellmap, const double wgt,
        const double test, const GEN::pairedvector<int,double>& deriv_test,
        MORTAR::MortarElement& moEle,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test
        );



    template <DRT::Element::DiscretizationType parentdistype, int dim>
    void inline SoEleGP(
        MORTAR::MortarElement& sele,
        const double wgt,
        const double* gpcoord,
        LINALG::Matrix<dim,1>& pxsi,
        LINALG::Matrix<dim,dim>& derivtrafo
  );

    template <int dim>
    void SoEleCauchy(
        MORTAR::MortarElement& moEle,
        double* boundary_gpcoord,
        std::vector<GEN::pairedvector<int,double> > boundary_gpcoord_lin,
        const double gp_wgt,
        const double* contactN,
        std::vector<GEN::pairedvector<int,double> >& contactN_deriv,
        const bool useEleN,
        const double w,
        double& cauchy_nn,
        GEN::pairedvector<int,double>& deriv_sigma_nn,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test);

protected:

    Teuchos::RCP<Epetra_FEVector> fc_;
    Teuchos::RCP<LINALG::SparseMatrix> kc_;
    double theta_;
    INPAR::CONTACT::NitscheWeighting nit_wgt_;
};
}
#endif // CONTACT_INTEGRATOR_NITSCHE_H
