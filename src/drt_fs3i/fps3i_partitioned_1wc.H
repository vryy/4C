/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with algorithmic routines for partitioned
     solution approaches to fluid-porous-structure-scalar-scalar interaction
     (FPS3I) specifically related to one-way-coupled problem
     configurations

\level 3

\maintainer Christoph Schmidt


*----------------------------------------------------------------------*/


#ifndef FPS3I_PART_1WC_H
#define FPS3I_PART_1WC_H


#include "fps3i_partitioned.H"


namespace FS3I
{
  class PartFPS3I_1WC : public PartFPS3I
  {
   public:
    //! constructor of one-way coupled FPS3I
    PartFPS3I_1WC(const Epetra_Comm& comm);

    //! Destructor of one-way coupled FPS3I
    virtual ~PartFPS3I_1WC(){};

    //! initialize this class
    virtual void Init();

    //! setup this class
    virtual void Setup();

    /// timeloop of coupled problem
    virtual void Timeloop();

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

#endif
