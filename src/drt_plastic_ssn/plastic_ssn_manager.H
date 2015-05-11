/*!----------------------------------------------------------------------
\file plastic_ssn_manager.H

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#ifndef PLASTIC_SSN_MANAGER_H_
#define PLASTIC_SSN_MANAGER_H_

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include "../drt_lib/drt_globalproblem_enums.H"
#include "../drt_inpar/inpar_structure.H"

// forward declarations
namespace DRT
{
class Discretization;

namespace UTILS
{

/*!
  \brief struct containing plasticity algorithm info needed by every element
 */
struct PlastSsnData {
  bool unconverged_active_set_;
  double pl_inc_;
  double pl_res_;
  bool eas_;
  double eas_res_;
  double eas_inc_;
  int num_active_;
  INPAR::STR::PredEnum pred_type_;
  bool no_pl_condensation_;
  bool no_recovery_;
  double scale_timint_;
  double dt_;
  bool split_res_;
  double cond_rhs_norm_;
  bool ls_;
  double alpha_ls_;
  int myPID;
};

/*!
\brief Class controlling semi-smooth Newton strategy for plasticity
*/
class PlastSsnManager
{
public:

  /*!
    \brief Standard Constructor
  */
  PlastSsnManager(Teuchos::RCP<DRT::Discretization> discret);

  /*!
    \brief Destructor
  */
  ~PlastSsnManager() {};

  //@}

  //! @name Evaluation methods

  /*!
  \brief Read and check plasticity input parameters

  All specified plasticity-related input parameters are read from the
  DRT::Problem::Instance() and stored into a local variable of
  type Teuchos::ParameterList. Invalid parameter combinations are
  sorted out and throw a dserror.

  */
  void ReadAndCheckInput();

  /*
  \brief set plastic parameters before element evaluation
  */
  void SetPlasticParams(Teuchos::ParameterList& params);

  /*
   \brief get plastic parameters from elements after evaluation
   */
  void GetPlasticParams(Teuchos::ParameterList& params);

  /*
  \brief check active set for convergence
  */
  bool ActiveSetConverged() {return !unconvergedactiveset_;};

  /*
  \brief check constraint residual for convergence
  */
  bool ConstraintConverged() {return lp_residual_norm_global_<pl_res_tol_;} ;

  /*
  \brief check plastic flow increment for convergence
  */
  bool IncrementConverged() {return lp_increment_norm_global_<pl_incr_tol_;};

  /*
  \brief check constraint residual for convergence
  */
  bool EasResConverged() {return eas_residual_norm_global_<eas_res_tol_;} ;

  /*
  \brief check plastic flow increment for convergence
  */
  bool EasIncrConverged() {return eas_increment_norm_global_<eas_incr_tol_;};

  /*
  \brief return total number of active Gauss points
  */
  int NumActivePlasticGP(){return numactive_global_;};

  /*
  \brief return L2 norm of plastic flow increment
  */
  double DeltaLp_increment_norm() {return lp_increment_norm_global_;};

  /*
  \brief return L2 norm of constraint equations
  */
  double DeltaLp_residual_norm() {return lp_residual_norm_global_;};

  /*
  \brief return L2 norm of plastic flow increment
  */
  double EAS_increment_norm() {return eas_increment_norm_global_;};

  /*
  \brief return L2 norm of constraint equations
  */
  double EAS_residual_norm() {return eas_residual_norm_global_;};

  /*
  \brief return parameter list of plastic algorithmic parameters
  */
  Teuchos::RCP<Teuchos::ParameterList> Params() {return plparams_; };

  /*
  \brief Do some plastic elements use EAS?
  */
  bool EAS() {return have_eas; };

  /*
  \brief Is this a TSI problem
  */
  bool TSI() {return probtype_==prb_tsi; };

  /*
  \brief Return the data container for the Semi-smooth Newton data
  */
  DRT::UTILS::PlastSsnData& SetData() { return data_; };

  /*
  \brief Print the data container for the Semi-smooth Newton data
  */
  void PrintData()
  {
    std::cout << "unconverged_active_set_: " <<  data_.unconverged_active_set_  << std::endl;
    std::cout << "pl_inc_: " <<  data_.pl_inc_  << std::endl;
    std::cout << "pl_res_: " <<  data_.pl_res_  << std::endl;
    std::cout << "eas_: " << data_.eas_   << std::endl;
    std::cout << "eas_res_: " <<  data_.eas_res_  << std::endl;
    std::cout << "eas_inc_: " <<  data_.eas_inc_  << std::endl;
    std::cout << "num_active_: " <<  data_.num_active_  << std::endl;
    std::cout << "pred_type_: " <<  data_.pred_type_  << std::endl;
    std::cout << "no_pl_condensation_: " << data_.no_pl_condensation_   << std::endl;
    std::cout << "no_recovery_: " <<   data_.no_recovery_ << std::endl;
    std::cout << "scale_timint_: " <<   data_.scale_timint_ << std::endl;
    std::cout << "dt_: " <<    data_.dt_<< std::endl;
    std::cout << "split_res_: " << data_.split_res_   << std::endl;
    std::cout << "cond_rhs_norm_: " <<   data_.cond_rhs_norm_ << std::endl;
    std::cout << "ls_: " << data_.ls_   << std::endl;
    std::cout << "alpha_ls_: " <<  data_.alpha_ls_  << std::endl;
    std::cout << "myPID: " <<  data_.myPID  << std::endl;
  }

private:

  Teuchos::RCP<DRT::Discretization> discret_; // current discretization
  Teuchos::RCP<Teuchos::ParameterList> plparams_; // plasticity algorithm input parameters
  PROBLEM_TYP probtype_;                          // problem type TSI or Structure

  int numactive_global_;// number of active plastic Gauss points on all procs
  bool unconvergedactiveset_; // true, if the active set has not changed in the last step

  double lp_increment_norm_global_; // L2-norm of global vector of plastic flow
  double lp_residual_norm_global_; // L2-norm of global vector of plastic constraint residuals

  bool have_eas;                    // is EAS technology used
  double eas_increment_norm_global_; // L2-norm of global vector EAS variables
  double eas_residual_norm_global_; // L2-norm of global vector of EAS residuals

  double pl_incr_tol_;      // convergence tolerance for plastic flow increment
  double pl_res_tol_;       // convergence tolerance for plastic constraints
  double eas_incr_tol_;     // convergence tolerance for EAS increment
  double eas_res_tol_;       // convergence tolerance for EAS residual

  struct DRT::UTILS::PlastSsnData data_;
};// class SsnPlastManager

}//namespace UTILS
} // namespace DRT

#endif // PLASTIC_SSN_MANAGER_H_
