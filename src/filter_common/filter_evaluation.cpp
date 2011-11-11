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

#ifdef CCADISCRET

#include "../drt_s8/shell8.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/micromaterial.H"

struct _PAR     par;
struct _GENPROB genprob;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

#if 0
#ifdef D_SHELL8
int ELEMENTS::Shell8Type::Initialize(Discretization&)
{
  return 0;
}

bool DRT::ELEMENTS::Shell8::ReadElement(const std::string& eletype,
                                        const std::string& distype,
                                        DRT::INPUT::LineDefinition* linedef)
{
  dserror("ELEMENTS::Shell8::ReadElement undefined");
  return false;
}

int ELEMENTS::Shell8::Evaluate(ParameterList&,
                               Discretization&,
                               vector<int>&,
                               Epetra_SerialDenseMatrix&,
                               Epetra_SerialDenseMatrix&,
                               Epetra_SerialDenseVector&,
                               Epetra_SerialDenseVector&,
                               Epetra_SerialDenseVector&)
{
  dserror("ELEMENTS::Shell8::Evaluate undefined");
  return 0;
}

int ELEMENTS::Shell8::EvaluateNeumann(ParameterList&, Discretization&, Condition&, vector<int>&, Epetra_SerialDenseVector&, Epetra_SerialDenseMatrix*)
{
  dserror("ELEMENTS::Shell8::EvaluateNeumann undefined");
  return 0;
}

int ELEMENTS::Shell8Line::EvaluateNeumann(ParameterList& params,
                                          Discretization&      discretization,
                                          Condition&           condition,
                                          vector<int>&              lm,
                                          Epetra_SerialDenseVector& elevec1,
                                          Epetra_SerialDenseMatrix* elemat1)
{
  dserror("ELEMENTS::Shell8Line::EvaluateNeumann undefined");
  return 0;
}

void DRT::ELEMENTS::Shell8::VisNames(map<string,int>& names)
{
}

bool DRT::ELEMENTS::Shell8::VisData(const string& name, vector<double>& data)
{
  return false;
}

#endif
#endif

namespace MAT
{

class MicroMaterialGP
{
};

}

void MAT::MicroMaterial::Evaluate(LINALG::Matrix<3,3>* defgrd,
                                  LINALG::Matrix<6,6>* cmat,
                                  LINALG::Matrix<6,1>* stress,
                                  double* density,
                                  const int gp,
                                  const int ele_ID)
{
  dserror("MAT::MicroMaterial::Evaluate not available");
}

void PrepareOutput()
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

void MAT::MicroMaterial::InvAnaInit(const bool eleowner)
{
  dserror("Mat::MicroMaterial::InvAna_Init not available");
}

#endif
