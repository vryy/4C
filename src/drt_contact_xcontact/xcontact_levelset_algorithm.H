/*----------------------------------------------------------------------------*/
/**
\file xcontact_levelset_algorithm.H

\brief level-set algorithm for the extended contact formulation

\maintainer Michael Hiermeier

\date Nov 22, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_ALGORITHM_H_
#define SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_ALGORITHM_H_

#include "../drt_levelset/levelset_algorithm.H"

// forward declarations
namespace LINALG
{
  class MultiMapExtractor;
}  // namespace LINALG
namespace XCONTACT
{
  namespace LEVELSET
  {
    namespace REINIT
    {
      class Generic;
      class Elliptic;
    }  // namespace REINIT

    class Intersection;

    enum MapPart
    {
      active = 0,
      inactive = 1
    };

    class Algorithm : public virtual SCATRA::LevelSetAlgorithm
    {
      friend class REINIT::Elliptic;

     public:
      /** \brief constructor
       *
       *  \param dis             (in) : ScaTra discretization
       *  \param solver          (in) : linear solver object
       *  \param ls_params       (in) : level-set control parameters
       *  \param scatradynparams (in) : scalar transport dynamic parameters
       *  \param extraparams     (in) : collection of parameters of other parameter lists
       *                                ( see adapter_scatra_base_algorithm for more info )
       *  \param output          (in) : Discretization writer object of the ScaTra discretization
       *
       *  \author hiermeier \date 11/16 */
      Algorithm(const Teuchos::RCP<DRT::Discretization>& dis,
          const Teuchos::RCP<LINALG::Solver>& solver,
          const Teuchos::RCP<Teuchos::ParameterList>& level_set_params,
          const Teuchos::RCP<Teuchos::ParameterList>& scatradynparams,
          const Teuchos::RCP<Teuchos::ParameterList>& extraparams,
          const Teuchos::RCP<IO::DiscretizationWriter>& output);

      /// destructor
      virtual ~Algorithm(){};

      /// setup member variables
      virtual void Setup();

      /// reinitialization of the level-set properties
      virtual void Reinitialization();

     protected:
      /** \brief Calculate the l2-projection of the gradient field if desired
       *
       *  \author hiermeier \date 12/16 */
      Teuchos::RCP<Epetra_MultiVector> ComputeNodalGradientL2Projection(const Epetra_Vector& phinp,
          const enum INPAR::SCATRA::L2ProjectionSystemType& l2_proj_system);

      /** \brief derived
       *
       *  Currently we support no initial contact states in the XCONTACT framework! */
      virtual void GetInitialVolumeOfMinusDomain(const Teuchos::RCP<const Epetra_Vector>& phinp,
          const Teuchos::RCP<const DRT::Discretization>& scatradis,
          double& volumedomainminus) const;

      /// does currently nothing for the xcontact case
      virtual void OutputOfLevelSetSpecificValues();

      void CreateActiveMaps(const Epetra_Vector& phinp);

      void CreateActiveMaps(const DRT::Discretization& dis, const Epetra_Vector& phinp,
          Teuchos::RCP<Epetra_Map>& active_row_elements, Teuchos::RCP<Epetra_Map>& active_row_nodes,
          Teuchos::RCP<Epetra_Map>& active_col_elements,
          Teuchos::RCP<LINALG::MultiMapExtractor>& active_dof_map_extractor) const;

      void CreateNewMap(const std::vector<int>& my_entries, const Epetra_Comm& comm,
          Teuchos::RCP<Epetra_Map>& new_map) const;

      /// set the solution state in the current active time integration
      virtual void SetState(const Epetra_Vector& sol_state) = 0;

      /// @name accessors
      /// @{

      /// \brief access the zero iso line as rc-pointer
      Teuchos::RCP<const std::map<int, GEO::BoundaryIntCellPtrs>> GetZeroIsoLinePtr() const
      {
        return Teuchos::rcp(&zero_iso_line_, false);
      }

      const Epetra_Map& ActiveRowDofMap() const;

      inline Teuchos::RCP<Epetra_Map>& ActiveRowEleMapPtr() { return active_row_elements_; }

      inline const Epetra_Map& ActiveRowEleMap() const
      {
        dsassert(active_row_elements_ != Teuchos::null, "Active element map was not initialized!");
        return *active_row_elements_;
      }

      inline Teuchos::RCP<Epetra_Map>& ActiveRowNodeMapPtr() { return active_row_nodes_; }

      inline const Epetra_Map& ActiveRowNodeMap() const
      {
        dsassert(active_row_nodes_ != Teuchos::null, "Active node map was not initialized!");
        return *active_row_nodes_;
      }

      inline Teuchos::RCP<Epetra_Map>& ActiveColEleMapPtr() { return active_col_elements_; }

      inline const Epetra_Map& ActiveColEleMap() const
      {
        dsassert(active_col_elements_ != Teuchos::null, "Active element map was not initialized!");
        return *active_col_elements_;
      }

      inline const LINALG::MultiMapExtractor& ActiveDofMapExtractor() const
      {
        dsassert(active_dof_map_extractor_ != Teuchos::null,
            "Active dof map extractor was not initialized!");
        return *active_dof_map_extractor_;
      }

      inline Teuchos::RCP<const Epetra_Vector> L2ProjSysMatDiagonalPtr() const
      {
        return l2_sys_mat_diagonal_.getConst();
      }


      /// @}

     private:
      /** \brief capture the zero iso line
       *
       *  This method calls the level-set intersection routine and performs a cut
       *  on the embedded contact interface. This is necessary to generate the boundary
       *  cells for the integration procedure.
       *
       *  \author hiermeier \date 11/16    */
      void CaptureInterface(std::map<int, GEO::BoundaryIntCellPtrs>& zero_iso_line);

      inline XCONTACT::LEVELSET::REINIT::Generic& ReinitStrategy()
      {
        if (reinit_strategy_.is_null())
          dserror("The reinitialization strategy was not initialized properly!");

        return *reinit_strategy_;
      }

      inline const XCONTACT::LEVELSET::REINIT::Generic& ReinitStrategy() const
      {
        if (reinit_strategy_.is_null())
          dserror("The reinitialization strategy was not initialized properly!");

        return *reinit_strategy_;
      }

     private:
      Teuchos::RCP<XCONTACT::LEVELSET::Intersection> ls_intersect_;

      std::map<int, GEO::BoundaryIntCellPtrs> zero_iso_line_;

      Teuchos::RCP<XCONTACT::LEVELSET::REINIT::Generic> reinit_strategy_;

      Teuchos::RCP<Epetra_Map> active_row_elements_;

      Teuchos::RCP<Epetra_Map> active_row_nodes_;

      Teuchos::RCP<Epetra_Map> active_col_elements_;

      Teuchos::RCP<LINALG::MultiMapExtractor> active_dof_map_extractor_;

      Teuchos::RCP<Epetra_Vector> l2_sys_mat_diagonal_;

    };  // class Algorithm
  }     // namespace LEVELSET
}  // namespace XCONTACT


#endif /* SRC_DRT_CONTACT_XCONTACT_XCONTACT_LEVELSET_ALGORITHM_H_ */
