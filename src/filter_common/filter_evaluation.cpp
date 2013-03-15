/*----------------------------------------------------------------------*/
/*!
\file post_drt_evaluation.cpp

\brief compatibility definitions

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Some discretization functions cannot be included in the filter build
because they use ccarat facilities that are not available inside the
filter. But to link the filter stubs of these functions are needed.

*/
/*----------------------------------------------------------------------*/



#include "../drt_mat/micromaterial.H"


namespace MAT
{

class MicroMaterialGP
{
};

}

void MAT::MicroMaterial::Evaluate(const LINALG::Matrix<3,3>* defgrd,
                                          const LINALG::Matrix<6,1>* glstrain,
                                          Teuchos::ParameterList& params,
                                          LINALG::Matrix<6,1>* stress,
                                          LINALG::Matrix<6,6>* cmat)
{
  dserror("MAT::MicroMaterial::Evaluate not available");
}

double MAT::MicroMaterial::Density() const
{
  dserror("MAT::MicroMaterial::Density not available");
  return 0.0;
}

void MAT::MicroMaterial::PrepareOutput()
{
  dserror("MAT::MicroMaterial::PrepareOutput not available");
}

void MAT::MicroMaterial::Output()
{
  dserror("MAT::MicroMaterial::Output not available");
}

void MAT::MicroMaterial::Update()
{
  dserror("MAT::MicroMaterial::Update not available");
}

void MAT::MicroMaterial::ReadRestart(const int gp, const int eleID, const bool eleowner)
{
  dserror("MAT::MicroMaterial::ReadRestart not available");
}

void MAT::MicroMaterial::InvAnaInit(bool eleowner, int eleID)
{
  dserror("Mat::MicroMaterial::InvAna_Init not available");
}

