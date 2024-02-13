/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with algorithmic routines for partitioned
     solution approaches to fluid-porous-structure-scalar-scalar interaction
     (FPS3I) specifically related to one-way-coupled problem
     configurations

\level 3



*----------------------------------------------------------------------*/


#ifndef BACI_FS3I_FPS3I_PARTITIONED_1WC_HPP
#define BACI_FS3I_FPS3I_PARTITIONED_1WC_HPP


#include "baci_config.hpp"

#include "baci_fs3i_fps3i_partitioned.hpp"

BACI_NAMESPACE_OPEN


namespace FS3I
{
  class PartFPS3I_1WC : public PartFPS3I
  {
   public:
    //! constructor of one-way coupled FPS3I
    PartFPS3I_1WC(const Epetra_Comm& comm);

    //! initialize this class
    void Init() override;

    //! setup this class
    void Setup() override;

    /// timeloop of coupled problem
    void Timeloop() override;

    /// FPSI step
    void DoFPSIStep();

    /// Scatra step
    void DoScatraStep();

    //! routine for preparing time step
    virtual void PrepareTimeStep();

    //! check convergence of monolithic ScaTra problem
    virtual bool ScatraConvergenceCheck(int itnum);
  };
}  // namespace FS3I

BACI_NAMESPACE_CLOSE

#endif
