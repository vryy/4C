/*!----------------------------------------------------------------------
\file microstrugenalpha.H

\class MicroStruGenAlpha

\brief Generalized Alpha time integration for microstructural problems

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#ifndef STRU_MICRO_H
#define STRU_MICRO_H

#ifdef CCADISCRET

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "ml_common.h"
#include "ml_include.h"
#include "nlnml_preconditioner.H"
#include "nlnml_linearsystem.H"
#include "nlnml_finelevelnoxinterface.H"
#include <Epetra_Map.h>
#include <Teuchos_RefCountPtr.hpp>

#include "../io/io_drt_micro.H"

using namespace Teuchos;
using namespace std;
using namespace IO;

/*!
\class MicroStruGenAlpha

\brief Generalized Alpha time integration for microstructural problems
in case of microscale analyses

\author gee (gee@lnm.mw.tum.de)
*/

class MicroStruGenAlpha
{
public:

  /*!
  \brief Standard Constructor -> currently no output for microscale possible

  */
  MicroStruGenAlpha(RefCountPtr<ParameterList> params,
                    RefCountPtr<DRT::Discretization> dis,
                    RefCountPtr<LINALG::Solver> solver);

  /*!
  \brief Destructor

  */
  virtual ~MicroStruGenAlpha();

  /*!
  \brief Create default parameter list

  */
  static void SetDefaults(ParameterList& params);

  /*!
  \brief Do time integration

  */
  void Integrate();

  /*!
  \brief do constant predictor step

  */
  void ConstantPredictor(const Epetra_SerialDenseMatrix* defgrd);

  /*!
  \brief do consistent constant predictor step

  */
  // void ConsistentPredictor();

  /*!
  \brief do full Newton iteration

  */
  void FullNewton();

  /*!
  \brief update configuration

  */
  void Update();


  /*!
  \brief write update

  */
  void Output(RefCountPtr<MicroDiscretizationWriter> output,
              const double time,
              const int istep);


  /// new stuff for multiscale analyses

  void DetermineToggle();

  void EvaluateMicroBC(const Epetra_SerialDenseMatrix* defgrd);

  void SetOldState(RefCountPtr<Epetra_Vector> disp,
                   RefCountPtr<Epetra_Vector> vel,
                   RefCountPtr<Epetra_Vector> acc,
                   RefCountPtr<Epetra_Vector> disi);

  void SetTime(double timen, int istep);

  RefCountPtr<Epetra_Vector> ReturnNewDisp();

  RefCountPtr<Epetra_Vector> ReturnNewVel();

  RefCountPtr<Epetra_Vector> ReturnNewAcc();

  RefCountPtr<Epetra_Vector> ReturnNewResDisp();

  void ClearState();

  void SetUpHomogenization();

  void Homogenization(Epetra_SerialDenseVector* stress,
                      Epetra_SerialDenseMatrix* cmat,
                      double *density,
                      const Epetra_SerialDenseMatrix* defgrd,
                      const string action);

  void StaticHomogenization(Epetra_SerialDenseVector* stress,
                            Epetra_SerialDenseMatrix* cmat,
                            double *density,
                            const Epetra_SerialDenseMatrix* defgrd);

  void calc_cmat(const Epetra_MultiVector& K_M,
                 const Epetra_SerialDenseMatrix& inv_defgrd,
                 const Epetra_SerialDenseVector& S,
                 Epetra_SerialDenseMatrix* cmat,
                 const Epetra_SerialDenseMatrix* F);


protected:

  // don't want = operator and cctor
  MicroStruGenAlpha operator = (const MicroStruGenAlpha& old);
  MicroStruGenAlpha(const MicroStruGenAlpha& old);

  RefCountPtr<ParameterList>                 params_;
  RefCountPtr<DRT::Discretization>           discret_;
  RefCountPtr<LINALG::Solver>                solver_;
  int                           myrank_;
  int                           maxentriesperrow_;
  double                        norm_;

  RefCountPtr<Epetra_CrsMatrix> stiff_;
  RefCountPtr<Epetra_CrsMatrix> mass_;
  RefCountPtr<Epetra_CrsMatrix> damp_;

  RefCountPtr<Epetra_Vector>    dirichtoggle_;
  RefCountPtr<Epetra_Vector>    invtoggle_;
  RefCountPtr<Epetra_Vector>    zeros_;
  RefCountPtr<Epetra_Vector>    dis_;
  RefCountPtr<Epetra_Vector>    vel_;
  RefCountPtr<Epetra_Vector>    acc_;
  RefCountPtr<Epetra_Vector>    disn_;
  RefCountPtr<Epetra_Vector>    veln_;
  RefCountPtr<Epetra_Vector>    accn_;
  RefCountPtr<Epetra_Vector>    dism_;
  RefCountPtr<Epetra_Vector>    velm_;
  RefCountPtr<Epetra_Vector>    accm_;
  RefCountPtr<Epetra_Vector>    disi_;
  RefCountPtr<Epetra_Vector>    fint_;
  RefCountPtr<Epetra_Vector>    fext_;
  RefCountPtr<Epetra_Vector>    fextm_;
  RefCountPtr<Epetra_Vector>    fextn_;
  RefCountPtr<Epetra_Vector>    fresm_;
  RefCountPtr<Epetra_Vector>    fresm_dirich_;

  double                         V0_;  // initial volume of RVE

  int                          ndof_;  // number of dofs overall
  int                            np_;  // number of boundary dofs
  RefCountPtr<Epetra_Vector>     Xp_;  // vector containing material
                                       // coordinates of boundary nodes
  RefCountPtr<Epetra_Map>      pdof_;  // prescribed dofs
  RefCountPtr<Epetra_Map>      fdof_;  // free dofs
  RefCountPtr<Epetra_Export> export_;
};

#endif
#endif
