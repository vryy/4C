/*----------------------------------------------------------------------*/
/*!
\file fsi_algorithm.H

\brief Basis of all FSI algorithms

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifndef FSI_ALGORITHM_H
#define FSI_ALGORITHM_H

#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_fluid_ale.H"
#include "fsi_coupling.H"


namespace FSI
{

  /// FSI algorithm base
  /*!

    Base class of FSI algorithms with generalized fluid field. There can (and
    will) be different subclasses that implement different coupling schemes.

    \note The generalized fluid field hides any ale or xfem handling of the
    variable fluid domain. This is the appropriate base class if direct access
    to ale or xfem is not required. If a coupling algorithm needs an ale
    field, MonolithicBase is the better choice for a base class.

    \warning The order of calling the three BaseAlgorithm-constructors (that
    is the order in which we list the base classes) is important here! In the
    constructors control file entries are written. And these entries define
    the order in which the filters handle the Discretizations, which in turn
    defines the dof number ordering of the Discretizations... Don't get
    confused. Just always list structure, fluid, ale. In that order.

    \author u.kue
    \date 02/08
   */
  class Algorithm : public ADAPTER::StructureBaseAlgorithm,
                    public ADAPTER::GeneralFluidBaseAlgorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit Algorithm(Epetra_Comm& comm);

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// read restart data
    void ReadRestart(int step);

    /// communicator
    const Epetra_Comm& Comm() const { return comm_; }

  protected:

    /// time step size
    double Dt() const { return dt_; }

    /// step number
    int Step() const { return step_; }

    //! @name Time loop building blocks

    /// tests if there are more time steps to do
    bool NotFinished() { return step_ < nstep_ and time_ <= maxtime_; }

    /// start a new time step
    void PrepareTimeStep();

    /// take current results for converged and save for next time step
    void Update();

    /// write output
    void Output();

    /// set method name for screen output
    void SetMethod(std::string method) { method_ = method; }

    //@}

    //! @name Transfer helpers

    //virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<Epetra_Vector> iv) const;
    //virtual Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const;

    //virtual Teuchos::RCP<Epetra_Vector> StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const;
    //virtual Teuchos::RCP<Epetra_Vector> AleToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const;
    virtual Teuchos::RCP<Epetra_Vector> FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

    Coupling& StructureFluidCoupling() { return coupsf_; }
    //Coupling& StructureAleCoupling()   { return coupsa_; }

    const Coupling& StructureFluidCoupling() const { return coupsf_; }
    //const Coupling& StructureAleCoupling()   const { return coupsa_; }

  private:

    /// coupling of structure and fluid at the interface
    Coupling               coupsf_;

    /// coupling of structure and ale at the interface
    //Coupling               coupsa_;

    /// comunication (mainly for screen output)
    Epetra_Comm& comm_;

    /// method name prepared for output
    std::string method_;

    //! @name Time stepping variables
    int step_;
    int nstep_;
    double time_;
    double maxtime_;
    double dt_;
    //@}
  };
}

#endif
#endif
