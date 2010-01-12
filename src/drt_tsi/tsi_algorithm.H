/*----------------------------------------------------------------------*/
/*!
\file tsi_algorithm.H

\brief Basis of all TSI algorithms that perform a coupling between the structure field
       equation and temperature field equations

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifndef TSI_ALGORITHM_H
#define TSI_ALGORITHM_H

/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
#include <Epetra_Comm.h>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "../drt_adapter/adapter_structure.H"
#include "../drt_adapter/adapter_thermo.H"
#include "../drt_adapter/adapter_algorithmbase.H"
//#include "../drt_adapter/adapter_coupling.H"


/*----------------------------------------------------------------------*
 | headers                                                   dano 12/09 |
 *----------------------------------------------------------------------*/
namespace TSI
{
  /// TSI algorithm base
  //!
  //!
  //!  Base class of TSI algorithms. Derives from StructureBaseAlgorithm and
  //!  ThermoBaseAlgorithm with temperature field.
  //!  There can (and will) be different subclasses that implement different
  //!  coupling schemes.
  //!
  //!  \warning The order of calling the two BaseAlgorithm-constructors (that
  //!  is the order in which we list the base classes) is important here! In the
  //!  constructors control file entries are written. And these entries define
  //!  the order in which the filters handle the Discretizations, which in turn
  //!  defines the dof number ordering of the Discretizations... Don't get
  //!  confused. Just always list structure, thermo. In that order.
  //!
  //!  \author u.kue
  //!  \date 02/08
  class Algorithm : public ADAPTER::AlgorithmBase,
                    public ADAPTER::StructureBaseAlgorithm,
                    public ADAPTER::ThermoBaseAlgorithm
  {
  public:

    /// create using a Epetra_Comm
    explicit Algorithm(Epetra_Comm& comm);

//    /// constructor using a Epetra_Comm
//    explicit Algorithm(
//      Epetra_Comm& comm,  ///< communicator
//      const Teuchos::ParameterList& prbdyn  ///< problem-specific parameters
//      );

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// outer level time loop (to be implemented by deriving classes)
    void TimeLoop();

    /// read restart data
    virtual void ReadRestart(
      int step  ///< step number where the calculation is continued
      );

  protected:

    //! @name Time loop building blocks

    /// start a new time step
    virtual void PrepareTimeStep();

    /// solve structure equations for current time step
    void DoStructureStep();

    /// solve tempterature equations for current time step
    void DoThermoStep();

    /// take current results for converged and save for next time step
    virtual void Update();

    /// write output
    virtual void Output();

    //@}

    //! @name Transfer helpers

//    virtual Teuchos::RCP<Epetra_Vector> StructToThermo(Teuchos::RCP<Epetra_Vector> iv) const;
//    virtual Teuchos::RCP<Epetra_Vector> ThermoToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
//
//    virtual Teuchos::RCP<Epetra_Vector> StructToThermo(Teuchos::RCP<const Epetra_Vector> iv) const;
//    virtual Teuchos::RCP<Epetra_Vector> ThermoToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;

    //@}

//    ADAPTER::TSCoupling& StructureThermoCoupling() { return coupst_; }
//
//    const ADAPTER::TSCoupling& StructureThermoCoupling() const { return coupst_; }

  private:

    // 14.12.09 not yet used!
//    /// coupling of structure and temperature
//    TSI::Algorithm  coupst_;
  };
} // namespace TSI

/*----------------------------------------------------------------------*/
#endif  // TSI_ALGORITHM_H
#endif  // CCADISCRET
