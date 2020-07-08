/*----------------------------------------------------------------------------*/
/** \file
\brief base class of the inequality level-set approach algorithm for contact
       problems a.k.a. xcontact or extended contact


\level 3
*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_ALGORITHM_BASE_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_ALGORITHM_BASE_H_


#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_lib/drt_dserror.H"

namespace DRT
{
  class DiscretizationInterface;
  class Discretization;
  class Condition;
}  // namespace DRT
namespace ADAPTER
{
  class StructureXContact;
  class Coupling;
}  // namespace ADAPTER
namespace XCONTACT
{
  class MultiDiscretizationWrapper;
  namespace LEVELSET
  {
    class Algorithm;
  }  // namespace LEVELSET
  namespace ALGORITHM
  {
    /*! \brief Base algorithm for the extended contact formulation
     *
     *  \author hiermeier
     *  \date 07/16 */
    class Base : public ADAPTER::AlgorithmBase
    {
     public:
      //! constructor
      Base();

      //! destructor
      virtual ~Base(){};

      /*! \brief Initialize class variables
       *
       *  \param p_xcontact_dyn  (in): X-Contact Dynamics parameter-list from dat-file (read-only)
       *  \param p_structure_dyn (in): Structural Dynamics parameter-list from dat-file (read-only)
       *  \param p_scatra_dyn    (in): ScaTra Dynamics parameter-list from dat-file (read-only)
       *  \param structdis       (in): full standard structure discretization
       *
       *  \date 06/16
       *  \author hiermeier
       */
      void Init(const Teuchos::ParameterList& p_xcontact_dyn,
          const Teuchos::ParameterList& p_structure_dyn, const Teuchos::ParameterList& p_scatra_dyn,
          const Teuchos::ParameterList& p_xfem_general,
          const Teuchos::RCP<DRT::Discretization>& structdis);

      //! setup class variables
      virtual void Setup();

      //! timeloop of coupled problem
      virtual void Timeloop();

      virtual void ReadRestart(int step) { dserror("ReadRestart() is not yet implemented!"); };

     protected:
      void SetupStructure();

      void SetupScaTra();

      DRT::Condition* ExtractSlaveCondition(const DRT::DiscretizationInterface& discret) const;

      void SetInitialScaTraParams(Teuchos::ParameterList& p_scatra_dyn) const;

      //! @name General time loop routines
      //! @{
      //! output the initial field at the beginning of the simulation
      void OutputInitialField();

      //! Prepare the different fields for the next time step
      virtual void PrepareTimeStep();

      /*! \brief Call the update routines in the different fields.
       *
       *  Basically we store the converged values of the current time
       *  step and update the step and time counters. \f$(n \leftarrow n+1)\f$
       */
      virtual void Update();

      /*! \brief Do outer iteration loop for particular type of algorithm
       *
       *  The corresponding implementation can be found in the derived algorithms. */
      virtual void OuterLoop() = 0;

      //! important for the structural output routine. This has to be done, before the
      //  Update() routines are called.
      void PrepareOutput();

      //! write output of the different fields
      virtual void Output();

      //! @}

      //! @ Set routines
      //! @{
      //! set the scaTra values in the structure field
      void SetScaTraValuesInStructure();

      //! set the structure values in the scaTra field
      bool SetStructureValuesInScaTra();
      //! @}

     protected:
      const XCONTACT::MultiDiscretizationWrapper& MultiDiscret() const
      {
        CheckInit();
        return *multi_discret_ptr_;
      }

      XCONTACT::MultiDiscretizationWrapper& MultiDiscret()
      {
        CheckInit();
        return *multi_discret_ptr_;
      }

      const Teuchos::RCP<DRT::DiscretizationInterface>& StructDiscretPtr();

      DRT::DiscretizationInterface& StructDiscret();

      const Teuchos::RCP<DRT::DiscretizationInterface>& XStructDiscretPtr();

      DRT::DiscretizationInterface& XStructDiscret();

      const Teuchos::RCP<DRT::Discretization>& ScaTraDiscretPtr();

      DRT::Discretization& ScaTraDiscret();

      ADAPTER::StructureXContact& StructureField()
      {
        CheckInitSetup();
        return *structure_ptr_;
      }

      XCONTACT::LEVELSET::Algorithm& ScaTraField()
      {
        CheckInitSetup();
        return *scatra_ptr_;
      }

      inline Teuchos::ParameterList& XContactDynParams()
      {
        if (p_xcontact_dyn_ptr_.is_null()) dserror("The parameter list ptr is NULL!");

        return *p_xcontact_dyn_ptr_;
      }

      inline Teuchos::ParameterList& StructDynParams()
      {
        if (p_struct_dyn_ptr_.is_null()) dserror("The parameter list ptr is NULL!");

        return *p_struct_dyn_ptr_;
      }

      inline Teuchos::ParameterList& ScaTraDynParams()
      {
        if (p_scatra_dyn_ptr_.is_null()) dserror("The parameter list ptr is NULL!");

        return *p_scatra_dyn_ptr_;
      }

      inline Teuchos::ParameterList& XFEMGeneralParams()
      {
        if (p_xfem_general_ptr_.is_null()) dserror("The parameter list ptr is NULL!");

        return *p_xfem_general_ptr_;
      }

      //! get the maximum number of steps
      inline const int& NumStep() const
      {
        CheckInit();
        return numstep_;
      }

      //! get the number of dofs per node
      inline const int& NumDofPerNode() const
      {
        if (num_dof_per_node_ < 1)
          dserror("The number of dofs per node is not meaningful! (%d<1)", num_dof_per_node_);
        return num_dof_per_node_;
      }

      inline const int& MaxNumDofSets() const
      {
        if (max_num_dof_sets_ < 1)
          dserror("The maximal number of dofsets is not meaningful! (%d<1)", max_num_dof_sets_);
        return max_num_dof_sets_;
      }

      //! check the initialization indicator
      inline void CheckInit() const
      {
        if (not isinit_) dserror("Call Init() first!");
      }

      //! check the initialization and setup indicators
      inline void CheckInitSetup() const
      {
        if (not issetup_ or not isinit_) dserror("Call Init() and Setup() first!");
      }

     protected:
      //! initialization flag
      bool isinit_;

      //! setup flag
      bool issetup_;

     private:
      //! copy of the x-contact dynamics parameter-list
      Teuchos::RCP<Teuchos::ParameterList> p_xcontact_dyn_ptr_;

      //! copy of the structural dynamics parameter-list
      Teuchos::RCP<Teuchos::ParameterList> p_struct_dyn_ptr_;

      //! copy of the scalar transport dynamics parameter-list
      Teuchos::RCP<Teuchos::ParameterList> p_scatra_dyn_ptr_;

      //! copy of the general xfem parameter-list
      Teuchos::RCP<Teuchos::ParameterList> p_xfem_general_ptr_;

      //! std::map containing pointers to the different slave interface scatra discretizations
      std::map<int, Teuchos::RCP<DRT::Discretization>> scatra_dis_map_;

      //! pointer to the structural x-contact wrapper
      Teuchos::RCP<ADAPTER::StructureXContact> structure_ptr_;

      //! pointer to the scatra/level-set algorithm
      Teuchos::RCP<XCONTACT::LEVELSET::Algorithm> scatra_ptr_;

      Teuchos::RCP<XCONTACT::MultiDiscretizationWrapper> multi_discret_ptr_;

      Teuchos::RCP<ADAPTER::Coupling> i_scatra_contact_coupling_;

      int numstep_;

      int num_dof_per_node_;

      int max_num_dof_sets_;
    };  // class Base
  }     // namespace ALGORITHM
}  // namespace XCONTACT


#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_ALGORITHM_BASE_H_ */
