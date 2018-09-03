/*!----------------------------------------------------------------------
\file immersed_partitioned_multiphysics.H

\brief base class for all multifield partitioned immersed algorithms

\level 2

<pre>
\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
 *----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_MULTIPHYSICS_H
#define IMMERSED_PARTITIONED_MULTIPHYSICS_H

#include "immersed_partitioned.H"


namespace ADAPTER
{
  class FSIStructureWrapperImmersed;
}

namespace POROELAST
{
  class PoroScatraBase;
}

namespace SSI
{
  class SSI_Part2WC;
}


namespace IMMERSED
{
  class ImmersedPartitionedAdhesionTraction;
  class ImmersedPartitionedFlowCellInteraction;

  class ImmersedPartitionedMultiphysics : public ImmersedPartitioned
  {
   public:
    explicit ImmersedPartitionedMultiphysics(const Epetra_Comm& comm);

    /*! \brief Initialize this object

    \warning none
    \return int
    \date 08/16
    \author rauch  */
    virtual int Init(const Teuchos::ParameterList& params);

    /*! \brief Setup all class internal objects and members

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

    /// read restart
    virtual void ReadRestart(int step);

    /// solve the multiphysics problem
    virtual void Solve();

   protected:
    /// perform one time step solver cascade
    virtual void TimeStep();

    /// prepare time step
    virtual void PrepareTimeStep();

    /// prepare output
    virtual void PrepareOutput();

    /// update counters and variables
    virtual void Update();

    /// write output
    virtual void Output();


   private:
    /// solve the cell-flow interaction module
    virtual void CFIOperator();

    /// solve the cell-ecm interaction module
    virtual void CEIOperator();

    /// my processor id
    int myrank_;

    /// pointer to global problem
    DRT::Problem* globalproblem_;

    /// immersed discretization
    Teuchos::RCP<DRT::Discretization> immerseddis_;

    /// discretization of structural phase in ECM
    Teuchos::RCP<DRT::Discretization> backgroundstructuredis_;

    /// discretization of liquid phase in ECM
    Teuchos::RCP<DRT::Discretization> backgroundfluiddis_;

    /// cell-flow interaction module
    Teuchos::RCP<IMMERSED::ImmersedPartitionedFlowCellInteraction> cfi_module_;

    /// cell-ECM adhesion interaction module
    Teuchos::RCP<IMMERSED::ImmersedPartitionedAdhesionTraction> adh_module_;

    /// cell structure of the immersed problem
    Teuchos::RCP<::ADAPTER::FSIStructureWrapperImmersed> cellstructure_;

    /// pointer to cellscatra subproblem (structure-scatra interaction)
    Teuchos::RCP<SSI::SSI_Part2WC> cellscatra_subproblem_;

    /// pointer to poroelast subproblem (poro-scatra interaction)
    Teuchos::RCP<POROELAST::PoroScatraBase> poroscatra_subproblem_;

  };  // class ImmersedPartitionedMultiphysics
}  // namespace IMMERSED

#endif
