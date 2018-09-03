/*!----------------------------------------------------------------------
\file immersed_partitioned_fsi_dirichletneumann_membrane.H

\brief partitioned immersed fsi algorithm for dirichlet-neumann coupling with membrane finite
elements

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/
#ifndef IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_MEMBRANE_H
#define IMMERSED_PARTITIONED_FSI_DIRICHLETNEUMANN_MEMBRANE_H

#include "immersed_partitioned_fsi_dirichletneumann.H"


namespace IMMERSED
{
  class ImmersedPartitionedFSIDirichletNeumannMembrane
      : public ImmersedPartitionedFSIDirichletNeumann
  {
   public:
    /// setup partitioned ImmersedFSI algorithm
    explicit ImmersedPartitionedFSIDirichletNeumannMembrane(const Epetra_Comm& comm);

   protected:
    /*!
    \brief calc the current artificial velocity by method of least squares from structure velocity
    onto fluid considering div*u=0 on all fluid nodes

    Calc new velocity in artificial fluid domain. After leaving this method the vector
    fluid_artificial_velocity_ contains the current velocity approximated from the structure onto
    the fluid.

    \note Artificial velocity is only approximated if the member artificial_velocity_isvalid_
    is=false on entry. if artificial_velocity_isvalid_=true the current artificial velocity stored
    in fluid_artificial_velocity_ is returned. This is important when restart is requested since
    then the we need to perform the projection in the first restart step. This is done in \ref
    InitialGuess(). In case of no restart, we can just return the velocity since it had recently
    been projected after the last structural solve. This work would be done twice, else.
    */
    virtual Teuchos::RCP<Epetra_Vector> CalcArtificialVelocity();

    /// improve quality of solution near the interface
    virtual void CorrectInterfaceVelocity();

    /// calc the fluid tractions on the structure
    virtual void CalcFluidTractionsOnStructure();

  };  // class ImmersedPartitionedFSIDirichletNeumannMembrane
}  // namespace IMMERSED

#endif
