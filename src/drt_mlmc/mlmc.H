/*!----------------------------------------------------------------------
\file  mlmc.H
\brief Class for performing Multi Level Monte Carlo (MLMC)analysis of structure


 <pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
 *!----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MLMC_H
#define MLMC_H

#include "Teuchos_RefCountPtr.hpp"
#include "../drt_adapter/adapter_structure.H"

namespace STR
{

/// Multi Level Monte Carlo Simulation of Structures with
/// some sort of stochastic property


class MLMC {

private:

  /// my structure discretisation to work at
  Teuchos::RCP<DRT::Discretization> discret_;

  /// linear algebraic solver
  Teuchos::RCP<LINALG::Solver> solver_;

  /// context for output and restart
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  ///  time integrator for structural dynamics
  Teuchos::RCP<STR::TimInt> sti_;

  // prolongator which maps all nodal results to finest grid

  Teuchos::RCP<Epetra_MultiVector> prolongator_disp_;
  // for stresses
  Teuchos::RCP<Epetra_MultiVector> prolongator_stress_;
  // Parameters from input file
  // Get number of Newton iterations
   int num_newton_it_;
   // Get convergence tolerance
   double convtol_;

   double InEleRange_ ;

   // number of timesteps
   int tsteps_;

  ///  number of measurment points the
  double tstep_;                    // timestep

   // number of run in MLMC
  int numb_run_;

  int reset_out_count_;


public:

  /// standard constructor
  MLMC(
    Teuchos::RCP<DRT::Discretization> dis,  ///< the discretisation
    Teuchos::RCP<LINALG::Solver> solver,  ///< algebraic solver
    Teuchos::RCP<IO::DiscretizationWriter> output  ///< outputer
    );

  // Setup prolongator to transfer results to finest grid
  void SetupProlongator();
  // Map results from one Discretization to another using the prolongator
  void ProlongateResults();
  ///
  void Integrate();

  // Find backgroundelement of a node
  void FindBackgroundElement(DRT::Node node, Teuchos::RCP<DRT::Discretization> background_dis, int* bg_ele_id, double* xsi);

  // check if a node lies within an element
  // returns double based on minimum value of local coordinates abscissa
  double CheckIfNodeInElement(DRT::Node& node ,DRT::Element& ele, double* xsi);

  bool EvaluateF(double* f, DRT::Node& node, DRT::Element& ele,const double* xsi);

  bool EvaluateGradF(LINALG::Matrix<3,3>& fgrad,DRT::Node& node, DRT::Element& ele,const double* xsi);

  // setup stoch mat law
  void SetupStochMat(unsigned int random_seed);


};  // class MLMC

}  // namespace STR

#endif /*MLMC_H_*/
#endif  // #ifdef CCADISCRET
