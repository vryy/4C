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
 |                                                           dano 12/09 |
 *----------------------------------------------------------------------*/
/// TSI: Thermo-Structure Interaction
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

    /// virtual destructor to support polymorph destruction
    virtual ~Algorithm();

    /// outer level time loop (to be implemented by deriving classes)
    void TimeLoop();

    /// time loop for TSI algorithm with one-way coupling
    void TimeLoopOneWay();

    /// time loop for TSI algorithm with iteration between fields (full coupling)
    void TimeLoopFull();

    /// read restart data
    virtual void ReadRestart(
      int step  ///< step number where the calculation is continued
      );

  protected:

    //! @name Time loop building blocks

    // one-way coupling

    /// solve tempterature equations for current time step
    /// and return temperatures (coupling towards structure field)
    Teuchos::RCP<Epetra_Vector> DoThermoStep();

    /// solve structure equations for current time step
    void DoStructureStep(Teuchos::RCP<Epetra_Vector> temp);

    /// solve structure equation and return displacements
    /// (coupling towards thermo field)
    Teuchos::RCP<Epetra_Vector> DoStructureStep();

    /// solve tempterature equations for current time step with respect to
    /// current velocities
    void DoThermoStep(
      Teuchos::RCP<Epetra_Vector> disp,
      Teuchos::RCP<Epetra_Vector> vel
      );

    // iterative (two-way) coupling

    /// Outer Iteration loop
    void OuterIterationLoop();

    /// solve structure equations for current time step
    Teuchos::RCP<Epetra_Vector> DoFullStructureStep(Teuchos::RCP<Epetra_Vector> temp);

    /// solve tempterature equations for current time step with respect to
    /// current velocities
    Teuchos::RCP<Epetra_Vector> DoFullThermoStep(
      Teuchos::RCP<Epetra_Vector> disp,
      Teuchos::RCP<Epetra_Vector> vel
      );

//    /// solve tempterature equations for current time step
//    void DoThermoStep(
//      Teuchos::RCP<Epetra_Vector> idisp
//      );


    //@}

    //! @name Access methods

    /// velocity calculation given the displacements (like in FSI)
   Teuchos::RCP<Epetra_Vector> CalcVelocity(
    Teuchos::RCP<const Epetra_Vector> idispnp
    ) const;

    //@}

    /// take current results for converged and save for next time step
    virtual void Update();

    /// write output
    virtual void Output();

    /// convergence check for natural convection solver
    bool ConvergenceCheck(int itnum, int itmax, double ittol);

    // apply modifications towards thermal contact
    void ApplyThermoContact(RCP<LINALG::SparseMatrix>& stiff,
                            RCP<Epetra_Vector>& rhs,
                            RCP<MORTAR::ManagerBase> cmtman);

    // convert maps form structure discretization to thermo discretization
    void ConvertMaps(RCP<Epetra_Map>& slavedofs,
                     RCP<Epetra_Map>& activedofs,
                     RCP<Epetra_Map>& masterdofs,
                     RCP<MORTAR::ManagerBase> cmtman);

    // apply mortar matrices for nodes in contact
    void AssembleDM(LINALG::SparseMatrix& dmatrix,
                    LINALG::SparseMatrix& mmatrix,
                    RCP<MORTAR::ManagerBase> cmtman);

    // apply mortar matrices for nodes in contact
    void AssembleMechDiss(Epetra_Vector& mechdiss,
                          RCP<MORTAR::ManagerBase> cmtman);

    // assemble thermal contact condition
    void AssembleThermContCondition(LINALG::SparseMatrix& thermcontLM,
                                    LINALG::SparseMatrix& thermcontTEMP,
                                    Epetra_Vector& thermcontRHS,
                                    LINALG::SparseMatrix& dmatrix,
                                    LINALG::SparseMatrix& mmatrix,
                                    RCP<Epetra_Map> activedofs,
                                    RCP<MORTAR::ManagerBase> cmtman);
    //@}

    /// maximum iteration steps
    int itmax_;
    /// temperature tolerance
    double ittol_;

    /// temperature increment of the outer loop
    Teuchos::RCP<Epetra_Vector> tempincnp_;

    //! @name Transfer helpers

////    virtual Teuchos::RCP<Epetra_Vector> StructToThermo(Teuchos::RCP<Epetra_Vector> iv) const;
//    virtual Teuchos::RCP<Epetra_Vector> ThermoToStruct(Teuchos::RCP<Epetra_Vector> iv) const;
////
////    virtual Teuchos::RCP<Epetra_Vector> StructToThermo(Teuchos::RCP<const Epetra_Vector> iv) const;
//    virtual Teuchos::RCP<Epetra_Vector> ThermoToStruct(Teuchos::RCP<const Epetra_Vector> iv) const;
//
//    //@}
//
//    TSI::Algorithm& StructureThermoCoupling() { return coupst_; }
//
//    const TSI::Algorithm& StructureThermoCoupling() const { return coupst_; }
//

    //const Epetra_Comm&          comm_;          // communicator

  private:

    /// displacements from time step begin (t_n)
    Teuchos::RCP<Epetra_Vector> dispn_;
    /// velocities from time step begin (t_n)
    Teuchos::RCP<Epetra_Vector> veln_;
    /// velocities from time step begin (t_n+1)
    Teuchos::RCP<Epetra_Vector> velnp_;

    /// temperature from time step begin (t_n)
    Teuchos::RCP<Epetra_Vector> tempn_;
//    /// temperature from time step begin (t_n+1)
//    Teuchos::RCP<Epetra_Vector> tempnp_;

    //! we couple based on displacements
    bool displacementcoupling_;

//    /// coupling of structure and temperature
//    TSI::Algorithm  coupst_;
  };
} // namespace TSI

/*----------------------------------------------------------------------*/
#endif  // TSI_ALGORITHM_H
#endif  // CCADISCRET
