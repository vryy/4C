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
#include "../drt_lib/linalg_utils.H"

//#include "../drt_stru_multi/microstrugenalpha.H"
#include "../drt_stru_multi/microstatic.H"
#include "../drt_io/io_micro.H"

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
// corresponding prototype in post_drt_evaluation.cpp is adapted, too!!

void MAT::MicroMaterial::Evaluate(const Epetra_SerialDenseMatrix* defgrd,
                                  Epetra_SerialDenseMatrix* cmat,
                                  Epetra_SerialDenseVector* stress,
                                  double* density,
                                  const int gp,
                                  const int ele_ID,
                                  const double time,
                                  const string action)
{
  // activate microscale material

  int microdisnum = matdata_->m.struct_multiscale->microdis;
  RefCountPtr<DRT::Problem> micro_problem = DRT::Problem::Instance(microdisnum);
  micro_problem->ActivateMaterial();

  // avoid writing output also for ghosted elements
  const bool eleowner = DRT::Problem::Instance(0)->Dis(genprob.numsf,0)->ElementRowMap()->MyGID(ele_ID);

  if (gp > static_cast<int>(matgp_.size())-1)
  {
    matgp_.resize(gp+1);
    matgp_[gp] = rcp(new MicroMaterialGP(gp, ele_ID, eleowner, time, microdisnum));
  }

  RefCountPtr<MicroMaterialGP> actmicromatgp = matgp_[gp];

  // read restart if necessary
  if (action == "multi_readrestart")
  {
    actmicromatgp->ReadRestart();
  }
  // perform microscale simulation and homogenization (if fint and stiff/mass or stress calculation is required)
  else
  {
    actmicromatgp->PerformMicroSimulation(defgrd, stress, cmat, density, time, eleowner);
  }

  // reactivate macroscale material
  RefCountPtr<DRT::Problem> macro_problem = DRT::Problem::Instance(0);
  macro_problem->ActivateMaterial();

  return;
}

#endif
