/*!----------------------------------------------------------------------
\file microstatic.H

\class MicroStatic

\brief Static control for microstructural problems

<pre>
Maintainer: Lena Wiechert
            wiechert@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>

*----------------------------------------------------------------------*/

#ifndef STRU_MICRO_STAT_H
#define STRU_MICRO_STAT_H

#ifdef CCADISCRET

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_surfstress/drt_surfstress_manager.H"

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
\class MicroStatic

\brief Static control for microstructural problems
in case of microscale analyses

\author wiechert (wiechert@lnm.mw.tum.de)
*/

class MicroStatic
{
public:

  /*!
  \brief Standard Constructor -> currently no output for microscale possible

  */
  MicroStatic(RCP<ParameterList> params,
              RCP<DRT::Discretization> dis,
              RCP<LINALG::Solver> solver);

  /*!
  \brief Destructor

  */
  virtual ~MicroStatic();

  /*!
  \brief Create default parameter list

  */
  static void SetDefaults(ParameterList& params);

  /*!
  \brief read restart

  */
  void ReadRestart(int step);

  /*!
  \brief do constant predictor step

  */
  void Predictor(const Epetra_SerialDenseMatrix* defgrd);

  /*!
  \brief do full Newton iteration

  */
  void FullNewton();

  /*!
  \brief update for a new time step

  */
  void UpdateNewTimeStep(RefCountPtr<Epetra_Vector> dis,
                         RefCountPtr<Epetra_Vector> dism);

  /*!
  \brief write update

  */
  void Output(RefCountPtr<MicroDiscretizationWriter> output,
              const double time,
              const int istep,
              const double dt);


  /// new stuff for multiscale analyses

  void DetermineToggle();

  void EvaluateMicroBC(const Epetra_SerialDenseMatrix* defgrd);

  void SetOldState(RefCountPtr<Epetra_Vector> dis,
                   RefCountPtr<Epetra_Vector> dism);

  void SetTime(double timen, int istep);

  RefCountPtr<Epetra_Vector> ReturnNewDism();

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

  void ConvertMat(const Epetra_MultiVector& cmatpf,
                  const Epetra_SerialDenseMatrix& inv_defgrd,
                  const Epetra_SerialDenseVector& S,
                  Epetra_SerialDenseMatrix& cmat);

  void calc_cmat(const Epetra_MultiVector& K_M,
                 const Epetra_SerialDenseMatrix& inv_defgrd,
                 const Epetra_SerialDenseVector& S,
                 Epetra_SerialDenseMatrix* cmat,
                 const Epetra_SerialDenseMatrix* F);

  bool Converged(const string type, const double disinorm,
                 const double resnorm, const double toldisp,
                 const double tolres);

  void CalcRefNorms();

  void PrintNewton(bool printscreen, bool print_unconv,
                   Epetra_Time timer, int numiter, int maxiter,
                   double fresmnorm, double disinorm, string convcheck);

   void PrintPredictor(string convcheck, double fresmnorm);



protected:

  // don't want = operator and cctor
  MicroStatic operator = (const MicroStatic& old);
  MicroStatic(const MicroStatic& old);

  RefCountPtr<DRT::SurfStressManager>   surf_stress_man_;

  RCP<ParameterList>            params_;
  RCP<DRT::Discretization>      discret_;
  RCP<LINALG::Solver>           solver_;
  int                           myrank_;
  int                           maxentriesperrow_;

  double                        ref_fnorm_;
  double                        ref_disnorm_;

  RefCountPtr<LINALG::SparseMatrix> stiff_;
  RefCountPtr<LINALG::SparseMatrix> stiff_dirich_;

  RefCountPtr<Epetra_Vector>    dirichtoggle_;
  RefCountPtr<Epetra_Vector>    invtoggle_;
  RefCountPtr<Epetra_Vector>    zeros_;
  RefCountPtr<Epetra_Vector>    dis_;
  RefCountPtr<Epetra_Vector>    disn_;
  RefCountPtr<Epetra_Vector>    dism_;
  RefCountPtr<Epetra_Vector>    disi_;
  RefCountPtr<Epetra_Vector>    fint_;
  RefCountPtr<Epetra_Vector>    fext_;
  RefCountPtr<Epetra_Vector>    fextm_;
  RefCountPtr<Epetra_Vector>    fextn_;
  RefCountPtr<Epetra_Vector>    fresm_;
  RefCountPtr<Epetra_Vector>    fresm_dirich_;

  RefCountPtr<Epetra_MultiVector> D_;  // D Matrix following Miehe et al., 2002
  RefCountPtr<Epetra_MultiVector> rhs_; // exported transpose of D (pdof -> dofrowmap)

  double                         V0_;  // initial volume of RVE
  double                    density_;  // initial density of RVE

  int                          ndof_;  // number of dofs overall
  int                            np_;  // number of boundary dofs
  RefCountPtr<Epetra_Vector>     Xp_;  // vector containing material
                                       // coordinates of boundary nodes
  RefCountPtr<Epetra_Map>      pdof_;  // prescribed dofs
  RefCountPtr<Epetra_Map>      fdof_;  // free dofs
  RefCountPtr<Epetra_Import> importp_;
  RefCountPtr<Epetra_Import> importf_;
};

#endif
#endif
