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
        theta_(params.get<double>("NITSCHE_THETA")),
        nit_wgt_(DRT::INPUT::IntegralValue<INPAR::CONTACT::NitscheWeighting>(params,"NITSCHE_WEIGHTING")),
        ppn_(imortar_.get<double>("PENALTYPARAM")),
        ppt_(imortar_.get<double>("PENALTYPARAMTAN")),
        frcoeff_(imortar_.get<double>("FRCOEFF",-1.)),
        frbound_(imortar_.get<double>("FRBOUND",-1.)),
        frtype_(DRT::INPUT::IntegralValue<INPAR::CONTACT::FrictionType>(imortar_,"FRICTION")),
        dt_(imortar_.get<double>("TIMESTEP"))
   {}

  /*!
   \brief Destructor

   */
  virtual ~CoIntegratorNitsche()
  {
  }

protected:

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
  virtual void IntegrateGP_3D(
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

private:

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
        const LINALG::Matrix<dim,1>& normal,
        std::vector<GEN::pairedvector<int,double> >& normal_deriv,
        const LINALG::Matrix<dim,1>& direction,
        std::vector<GEN::pairedvector<int,double> >& direction_deriv,
        const double w,
        double& cauchy_nt,
        GEN::pairedvector<int,double>& deriv_sigma_nt,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test);

protected:

    template<int dim>
        void BuildAdjointTest(
        MORTAR::MortarElement& moEle,
        const double fac,
        const Epetra_SerialDenseMatrix& dsntdd,
        const Epetra_SerialDenseMatrix& d2sntdd2,
        const Epetra_SerialDenseMatrix& d2sntDdDn,
        const Epetra_SerialDenseMatrix& d2sntDdDt,
        const Epetra_SerialDenseMatrix& d2sntDdDpxi,
        const std::vector<GEN::pairedvector<int,double> > boundary_gpcoord_lin,
        LINALG::Matrix<dim,dim> derivtravo_slave,
        const std::vector<GEN::pairedvector<int,double> >& normal_deriv,
        const std::vector<GEN::pairedvector<int,double> >& direction_deriv,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test);

    template <int dim>
    void IntegrateTest(
        const double fac,
        MORTAR::MortarElement& ele,
        const LINALG::SerialDenseVector& shape,
        const LINALG::SerialDenseMatrix& deriv,
        const std::vector<GEN::pairedvector<int,double> >& dxi,
        const double jac,const GEN::pairedvector<int,double>& jacintcellmap, const double wgt,
        const double test_val, const GEN::pairedvector<int,double>& test_deriv,
        const LINALG::Matrix<dim,1>& test_dir, const std::vector<GEN::pairedvector<int,double> >& test_dir_deriv
        );

    template<int dim>
    void IntegrateAdjointTest(
        const double fac,
        const double jac, const GEN::pairedvector<int,double>& jacintcellmap, const double wgt,
        const double test, const GEN::pairedvector<int,double>& deriv_test,
        MORTAR::MortarElement& moEle,
        LINALG::SerialDenseVector& adjoint_test,
        GEN::pairedvector<int,LINALG::SerialDenseVector>& deriv_adjoint_test
        );


    double theta_;
    INPAR::CONTACT::NitscheWeighting nit_wgt_;
    double ppn_;
    double ppt_;
    double frcoeff_; // coulomb
    double frbound_; // tresca
    INPAR::CONTACT::FrictionType frtype_;
    double dt_;

};

namespace UTILS
{
template <int dim>

// actually not the velocity but the displacement increment
void RelVelInvariant(
    MORTAR::MortarElement& sele,
    const double* sxi,
    const std::vector<GEN::pairedvector<int,double> >& derivsxi,
    const LINALG::SerialDenseVector& sval,
    const LINALG::SerialDenseMatrix& sderiv,
    MORTAR::MortarElement& mele,
    const double* mxi,
    const std::vector<GEN::pairedvector<int,double> >& derivmxi,
    const LINALG::SerialDenseVector& mval,
    const LINALG::SerialDenseMatrix& mderiv,
    const double& gap,
    const GEN::pairedvector<int, double>& deriv_gap,
    LINALG::Matrix<dim,1>& relVel,
    std::vector<GEN::pairedvector<int,double> >& relVel_deriv,
    const double fac=1.
);

template <int dim>
void RelVel(
    MORTAR::MortarElement& ele,
    const LINALG::SerialDenseVector& shape,
    const LINALG::SerialDenseMatrix& deriv,
    const std::vector<GEN::pairedvector<int,double> >& dxi,
    const double fac,
    LINALG::Matrix<dim,1>& relVel,
    std::vector<GEN::pairedvector<int,double> >& relVel_deriv
    );

template <int dim>
void VectorScalarProduct(
    const LINALG::Matrix<dim,1>& v1,
    const std::vector<GEN::pairedvector<int,double> >& v1d,
    const LINALG::Matrix<dim,1>& v2,
    const std::vector<GEN::pairedvector<int,double> >& v2d,
    double& val,
    GEN::pairedvector<int,double>& val_deriv
);

void BuildTangentVectors_3D(
    const double* np,
    const std::vector<GEN::pairedvector<int,double> >& dn,
    double* t1p,
    std::vector<GEN::pairedvector<int,double> >& dt1,
    double* t2p,
    std::vector<GEN::pairedvector<int,double> >& dt2
    );

template <int dim>
void BuildTangentVectors(
    const double* np,
    const std::vector<GEN::pairedvector<int,double> >& dn,
    double* t1p,
    std::vector<GEN::pairedvector<int,double> >& dt1,
    double* t2p,
    std::vector<GEN::pairedvector<int,double> >& dt2
    );
}
}
#endif // CONTACT_INTEGRATOR_NITSCHE_H
