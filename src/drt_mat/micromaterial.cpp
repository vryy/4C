#ifdef D_SOH8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "micromaterial.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"

using namespace std;
using namespace Teuchos;

/// construct an instance of MicroMaterial for a given Gauss point and
/// microscale discretization

MAT::MicroMaterial::MicroMaterial(int gp, int microdis_num)
  : gp_(gp)
{
  RefCountPtr<DRT::Problem> microproblem = DRT::Problem::Instance(microdis_num);
  microdis_ = microproblem->Dis(0, 0);
  disp_ = LINALG::CreateVector(*microdis_->DofRowMap(),true);
  return;
}

/// destructor

MAT::MicroMaterial::~MicroMaterial()
{ }


/// test routine for calculating stresses, constitutive matrix and density in
/// case of St Venant Kirchhoff material

void MAT::MicroMaterial::CalcStressStiffDens (Epetra_SerialDenseVector* stress,
                          Epetra_SerialDenseMatrix* cmat,
                          double* density,
                          const Epetra_SerialDenseVector* glstrain)
{
  double Emod = 1.0;    // Young's modulus (modulus of elasticity)
  double nu = 0.0;      // Poisson's ratio (Querdehnzahl)
  (*density) = 0.0;     //
                                            // density, returned to evaluate mass matrix
  double mfac = Emod/((1.0+nu)*(1.0-2.0*nu));  /* factor */
  /* write non-zero components */
  (*cmat)(0,0) = mfac*(1.0-nu);
  (*cmat)(0,1) = mfac*nu;
  (*cmat)(0,2) = mfac*nu;
  (*cmat)(1,0) = mfac*nu;
  (*cmat)(1,1) = mfac*(1.0-nu);
  (*cmat)(1,2) = mfac*nu;
  (*cmat)(2,0) = mfac*nu;
  (*cmat)(2,1) = mfac*nu;
  (*cmat)(2,2) = mfac*(1.0-nu);
  /* ~~~ */
  (*cmat)(3,3) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(4,4) = mfac*0.5*(1.0-2.0*nu);
  (*cmat)(5,5) = mfac*0.5*(1.0-2.0*nu);

  // evaluate stresses
  (*cmat).Multiply('N',(*glstrain),(*stress));   // sigma = C . epsilon
}

#endif
#endif
#endif
