#ifdef CCADISCRET

#ifndef ADAPTER_FLUID_XFEM_H
#define ADAPTER_FLUID_XFEM_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "adapter_general_fluid.H"
#include "adapter_fluid_base_algorithm.H"
#include "../drt_fsi/fsi_coupling.H"

namespace ADAPTER
{

  /// fluid with moving interfaces implemented by the XFEM
  class FluidXFEMAdapter : public GeneralFluid
  {
  public:

    FluidXFEMAdapter(const Teuchos::ParameterList& prbdyn,
                    std::string condname);

    /// fluid field
    Fluid& FluidField() { return fluid_.FluidField(); }
    const Fluid& FluidField() const { return fluid_.FluidField(); }

    virtual Teuchos::RCP<DRT::Discretization> Discretization();
    virtual const LINALG::MapExtractor& Interface() const;

    virtual void PrepareTimeStep();
    virtual void Update();
    virtual void Output();
    virtual double ReadRestart(int step);

    virtual void NonlinearSolve(Teuchos::RCP<Epetra_Vector> idisp,Teuchos::RCP<Epetra_Vector> ivel);
    virtual Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> idisp,
                                                        double dt);

    virtual Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces();

    virtual int Itemax() const { return FluidField().Itemax(); }
    virtual void SetItemax(int itemax) { FluidField().SetItemax(itemax); }

    virtual Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape();

    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

  protected:

//    //! @name Transfer helpers
//
//    /// field transform
//    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const;
//
//    /// field transform
//    virtual Teuchos::RCP<Epetra_Vector> AleToFluidField(Teuchos::RCP<const Epetra_Vector> iv) const;
//
//    /// interface transform
//    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const;
//
//    /// interface transform
//    virtual Teuchos::RCP<Epetra_Vector> FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const;
//
//    //@}
    
  private:

    FluidBaseAlgorithm fluid_;

    /// coupling of fluid and ale (interface only)
    FSI::Coupling               icoupfa_;

  };

}

#endif
#endif
