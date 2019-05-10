/*-----------------------------------------------------------*/
/*!
\file contact_aug_parallel_distribution_controller.H

\maintainer Michael Hiermeier

\date May 25, 2018

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_PARALLEL_DISTRIBUTION_CONTROLLER_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_PARALLEL_DISTRIBUTION_CONTROLLER_H_

#include "contact_aug_timemonitor.H"
#include "../drt_mortar/mortar_paramsinterface.H"

#include <Teuchos_RCP.hpp>

class Epetra_Vector;
namespace CONTACT
{
  class ParamsInterface;
  namespace AUG
  {
    class DataContainer;
    class Strategy;

    /*-------------------------------------------------------------------------*/
    /** \brief Parallel distribution controller
     *
     *  This class controls the dynamic parallel redistribution in case of the
     *  contact augmented framework and tries to achieve a load balance over all
     *  processors based on the time measurements by the used time monitor object.
     *
     *  \athour hiermeier \date 05/18 */
    class ParallelDistributionController
    {
     public:
      /// constructor
      ParallelDistributionController(CONTACT::AUG::Strategy& strat,
          CONTACT::AUG::DataContainer& data, const int par_redist_interval)
          : strat_(strat),
            data_(data),
            interval_((par_redist_interval > 0 ? par_redist_interval
                                               : std::numeric_limits<int>::max())){/* empty */};

      /// destructor
      ~ParallelDistributionController() = default;

      /// setup the controller
      void setup(CONTACT::ParamsInterface& cparams);

      /** \brief check the balance of the current parallel distribution
       *
       *  This check is based on the element evaluation times collected by the
       *  used time monitor object. The check is performed during here specified
       *  evaluation actions while the respective time measurements are averaged.
       *
       *  \note that these evaluations might correspond to quite different displacement
       *  fields: While force evaluation takes place at the end of a successful
       *  Newton iterate and, therefore, the already updated displacement field is
       *  considered, the stiffness matrix is evaluated at the beginning and, hence,
       *  corresponds to the previously accepted iterate.
       *
       *  \author hiermeier */
      void check(CONTACT::ParamsInterface& cparams);

      /** \brief Initiate a redistribution if the set \c interval_ suggests it
       *
       *  \return TRUE if a redistribution took place. Return FALSE otherwise. */
      bool redistribute(const Teuchos::RCP<const Epetra_Vector>& dis, const int nlniter);

     private:
      /// reference to the surrounding contact strategy
      Strategy& strat_;

      /// reference to the data container of this strategy
      DataContainer& data_;

      /// global timer object used for the time measurements
      GlobalTimeMonitor global_timer_;

      /** iteration interval which specifies the check interval for the
       *  dynamic parallel redistribution during the Newton iterations. */
      const int interval_;

      /// current action type
      MORTAR::ActionType acttype_ = MORTAR::ActionType::eval_none;

      /// pointer containing the evaluation times of each slave element
      Teuchos::RCP<Epetra_Vector> sele_eval_times_ = Teuchos::null;
    };
  }  // namespace AUG
}  // namespace CONTACT



#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_PARALLEL_DISTRIBUTION_CONTROLLER_H_ */
