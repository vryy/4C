/*----------------------------------------------------------------------*/
/*!
\file micromaterial_evaluate.cpp

\brief class for handling of micro-macro transitions

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "micromaterial.H"
//#include "micromaterialgp.H"
#include "micromaterialgp_static.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"

#include "../drt_stru_multi/microstatic.H"

using namespace std;
using namespace Teuchos;
using namespace IO;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


// This function has to be separated from the remainder of the
// MicroMaterial class. MicroMaterialGP is NOT a member of
// FILTER_OBJECTS hence the MicroMaterial::Evaluate function that
// builds the connection to MicroMaterialGP is not either. In
// post_drt_evaluation.cpp this function is defined to content the
// compiler. If during postprocessing the MicroMaterial::Evaluate
// function should be called, an error is invoked.
//
// -> see also Makefile.objects and setup-objects.sh
//
// In case of any changes of the function prototype make sure that the
// corresponding prototype in src/filter_common/filter_evaluation.cpp is adapted, too!!

void MAT::MicroMaterial::Evaluate(LINALG::Matrix<3,3>* defgrd,
                                  LINALG::Matrix<6,6>* cmat,
                                  LINALG::Matrix<6,1>* stress,
                                  double* density,
                                  const int gp,
                                  const int ele_ID)
{
  // activate microscale material

  int microdisnum = MicroDisNum();
  double V0 = InitVol();
  RefCountPtr<DRT::Problem> micro_problem = DRT::Problem::Instance(microdisnum);
  DRT::Problem::Instance()->Materials()->SetReadFromProblem(microdisnum);

  // avoid writing output also for ghosted elements
  const bool eleowner = DRT::Problem::Instance(0)->Dis(genprob.numsf,0)->ElementRowMap()->MyGID(ele_ID);

  if (matgp_.find(gp) == matgp_.end())
  {
    matgp_[gp] = rcp(new MicroMaterialGP(gp, ele_ID, eleowner, microdisnum, V0));
  }

  RefCountPtr<MicroMaterialGP> actmicromatgp = matgp_[gp];

  // perform microscale simulation and homogenization (if fint and stiff/mass or stress calculation is required)
  actmicromatgp->PerformMicroSimulation(defgrd, stress, cmat, density);

  // reactivate macroscale material
  DRT::Problem::Instance()->Materials()->ResetReadFromProblem();

  return;
}


void MAT::MicroMaterial::Update()
{
  std::map<int, RefCountPtr<MicroMaterialGP> >::iterator it;
  for (it=matgp_.begin(); it!=matgp_.end(); ++it)
  {
    RefCountPtr<MicroMaterialGP> actmicromatgp = (*it).second;
    actmicromatgp->Update();
  }
}


void MAT::MicroMaterial::PrepareOutput()
{
  std::map<int, RefCountPtr<MicroMaterialGP> >::iterator it;
  for (it=matgp_.begin(); it!=matgp_.end(); ++it)
  {
    RefCountPtr<MicroMaterialGP> actmicromatgp = (*it).second;
    actmicromatgp->PrepareOutput();
  }
}


void MAT::MicroMaterial::Output()
{
  std::map<int, RefCountPtr<MicroMaterialGP> >::iterator it;
  for (it=matgp_.begin(); it!=matgp_.end(); ++it)
  {
    RefCountPtr<MicroMaterialGP> actmicromatgp = (*it).second;
    actmicromatgp->Output();
  }
}


void MAT::MicroMaterial::ReadRestart(const int gp, const int eleID, const bool eleowner)
{
  if (matgp_.find(gp) == matgp_.end())
  {
    int microdisnum = MicroDisNum();
    double V0 = InitVol();
    matgp_[gp] = rcp(new MicroMaterialGP(gp, eleID, eleowner, microdisnum, V0));
  }

  RefCountPtr<MicroMaterialGP> actmicromatgp = matgp_[gp];
  actmicromatgp->ReadRestart();
}


void MAT::MicroMaterial::InvAnaInit(const bool eleowner)
{
  std::map<int, RefCountPtr<MicroMaterialGP> >::iterator it;
  for (it=matgp_.begin(); it!=matgp_.end(); ++it)
  {
    RefCountPtr<MicroMaterialGP> actmicromatgp = (*it).second;
    actmicromatgp->ResetTimeAndStep();
    std::string newfilename;
    actmicromatgp->NewResultFile(eleowner, newfilename);
  }
}
#endif
