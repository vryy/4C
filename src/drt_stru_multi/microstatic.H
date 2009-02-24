/*!----------------------------------------------------------------------
\file microstatic.H

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
#include "../drt_lib/linalg_sparsematrix.H"
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

#include "../drt_io/io.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_structure/str_aux.H"


using namespace Teuchos;
using namespace std;
using namespace IO;

/*!
\brief Namespace for structural multi-scale related algorithms

\author wiechert (wiechert@lnm.mw.tum.de)
*/
namespace STRUMULTI
{

/*!
\brief Quasi-static control for microstructural analysis
in case of multi-scale problems

Note that implementation currently only holds for imr-like generalized
alpha time integration. Corresponding functions (e.g. UpdateNewTimeStep,
but also calls to SurfaceStressManager!) need to be adapted accordingly
if usage of other time integration schemes should be enabled.

\author wiechert (wiechert@lnm.mw.tum.de)
*/

class MicroStatic
{
public:

  /*!
  \brief Standard Constructor

  */
  MicroStatic(const int microdisnum, const double V0);

  /*!
  \brief Destructor

  */
  virtual ~MicroStatic();

  /*!
  \brief Read restart

  */
  void ReadRestart(int step,
                   RCP<Epetra_Vector> dis,
                   RCP<std::map<int, RCP<Epetra_SerialDenseMatrix> > > lastalpha,
                   RCP<UTILS::SurfStressManager> surf_stress_man,
                   string name);

  /*!
  \brief Return time from parameter list

  */
  double GetTime()
  {
    return time_;
  }

  /*!
  \brief Predictor step

  */
  void Predictor(LINALG::Matrix<3,3>* defgrd);

  /*!
  \brief Predictor step

  */
  void PredictConstDis(LINALG::Matrix<3,3>* defgrd);

   /*!
  \brief Predictor step

  */
  void PredictTangDis(LINALG::Matrix<3,3>* defgrd);

  /*!
  \brief Full Newton iteration

  */
  void FullNewton();

  /*!
  \brief Update for a new time step in case of imr-like generalized
  alpha time integration

  */
  void UpdateNewTimeStep(RefCountPtr<Epetra_Vector> dis,
                         RefCountPtr<Epetra_Vector> dism,
                         RefCountPtr<Epetra_Vector> disn,
                         RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldalpha,
                         RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > alpha,
                         RefCountPtr<UTILS::SurfStressManager> surf_stress_man);

  /*!
  \brief Write output and (possibly) restart

  */
  void Output(RefCountPtr<DiscretizationWriter> output,
              const double time,
              const int istep,
              const double dt);

  /*!
  \brief Determine toggle vector identifying prescribed boundary dofs

  */
  void DetermineToggle();

  /*!
  \brief Evaluate microscale boundary displacement according to
  associated macroscale deformation gradient

  */
  void EvaluateMicroBC(LINALG::Matrix<3,3>* defgrd, RefCountPtr<Epetra_Vector> disp);

  /*!
  \brief Set old state given from micromaterialgp

  */
  void SetOldState(RefCountPtr<Epetra_Vector> dis,
                   RefCountPtr<Epetra_Vector> dism,
                   RefCountPtr<Epetra_Vector> disn,
                   RefCountPtr<UTILS::SurfStressManager> surfman,
                   RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > lastalpha,
                   RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldalpha,
                   RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldfeas,
                   RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKaainv,
                   RefCountPtr<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKda);

  /*!
  \brief Set time and step

  */
  void SetTime(const double timen, const double dt, const int istep);

  /*!
  \brief Clear all displacement states

  */
  void ClearState();

  /*!
  \brief Set up everything for homogenization
  (e.g. calculation of matrix D containing reference boundary coordinates)

  */
  void SetUpHomogenization();

  /*!
  \brief Perform homogenization, i.e. calculate second Piola-Kirchhoff
  stresses and constitutive tensor by averaging over RVE

  */
  void StaticHomogenization(LINALG::Matrix<6,1>* stress,
                            LINALG::Matrix<6,6>* cmat,
                            double *density,
                            LINALG::Matrix<3,3>* defgrd,
                            const bool mod_newton,
                            bool& build_stiff);

  /*!
  \brief Convert constitutive tensor relating first Piola-Kirchhoff
  stresses and deformation gradient to tensor relating second
  Piola-Kirchhoff stresses and Green-Lagrange strains

  For details cf.

  Marsden and Hughes, Mathematical Foundations of Elasticity,
  Dover, pg. 215
  */
  void ConvertMat(const Epetra_MultiVector& cmatpf,
                  const LINALG::Matrix<3,3>& inv_defgrd,
                  const LINALG::Matrix<6,1>& S,
                  LINALG::Matrix<6,6>& cmat);


  /*!
  \brief Check for Newton convergence
  */

  bool Converged();

  /*!
  \brief Calculate reference norms for relative convergence checks

  */
  void CalcRefNorms();

  /*!
  \brief Output of Newton details

  Note that this is currently disabled for the sake of clearness
  */
  void PrintNewton(bool print_unconv, Epetra_Time timer);

  /*!
  \brief Output of predictor details

  Note that this is currently disabled for the sake of clearness
  */
  void PrintPredictor();



protected:

  // don't want = operator and cctor
  MicroStatic operator = (const MicroStatic& old);
  MicroStatic(const MicroStatic& old);

  RefCountPtr<UTILS::SurfStressManager>   surf_stress_man_;

  RCP<DRT::Discretization>      discret_;
  RCP<LINALG::Solver>           solver_;
  int                           myrank_;
  int                           maxentriesperrow_;

  double                        beta_;
  double                        gamma_;
  double                        alphaf_;
  double                        alpham_;

  double                        dt_;
  double                        time_;

  INPAR::STR::PredEnum          pred_;  //!< predictor

  bool                          isadapttol_;
  double                        adaptolbetter_;

  int                           maxiter_;
  int                           numiter_;
  int                           numstep_;
  int                           step_;

  bool                          iodisp_;
  int                           resevrydisp_;
  INPAR::STR::StressType        iostress_;
  int                           resevrystrs_;
  INPAR::STR::StrainType        iostrain_;
  bool                          iosurfactant_;
  int                           restart_;
  int                           restartevry_;
  bool                          printscreen_;

  INPAR::STR::ConvCheck         convcheck_;
  INPAR::STR::VectorNorm        iternorm_;
  double                        tolres_;
  double                        toldis_;

  double                        ref_resnorm_;
  double                        resnorm_;
  double                        ref_disinorm_;
  double                        disinorm_;

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
  RefCountPtr<Epetra_Vector>    fresm_;
  RefCountPtr<Epetra_Vector>    fresm_dirich_;

  //EAS history data
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > lastalpha_;
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldalpha_;
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldfeas_;
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKaainv_;
  RCP<std::map<int, RefCountPtr<Epetra_SerialDenseMatrix> > > oldKda_;

  RefCountPtr<Epetra_MultiVector> D_;  //!< D Matrix following Miehe et al., 2002
  RefCountPtr<Epetra_MultiVector> rhs_; //!< exported transpose of D (pdof -> dofrowmap)

  int                   microdisnum_;  //!< number of RVE

  double                         V0_;  //!< initial volume of RVE
  double                    density_;  //!< initial density of RVE

  int                          ndof_;  //!< number of dofs overall
  int                            np_;  //!< number of boundary dofs
  RefCountPtr<Epetra_Vector>     Xp_;  //!< vector containing material
                                       //!< coordinates of boundary nodes
  RefCountPtr<Epetra_Map>      pdof_;  //!< prescribed dofs
  RefCountPtr<Epetra_Map>      fdof_;  //!< free dofs
  RefCountPtr<Epetra_Import> importp_;
  RefCountPtr<Epetra_Import> importf_;
};



} // namespace STRUMULTI
#endif
#endif
