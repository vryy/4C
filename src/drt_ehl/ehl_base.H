/*--------------------------------------------------------------------------*/
/*!
\file ehl_base.H

\brief base class for all elastohydrodynamic lubrication (lubrication structure interaction) algorithms

\level 3

<pre>
\maintainer Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_EHL_EHL_BASE_H_
#define SRC_DRT_EHL_EHL_BASE_H_


#include "../drt_adapter/adapter_algorithmbase.H"

#include <Epetra_Vector.h>

#include "../drt_inpar/inpar_ehl.H"

namespace ADAPTER
{
  class LubricationBaseAlgorithm;
  class Structure;
  class CouplingMortar;
  class MortarVolCoupl;
}

namespace LINALG
{
  class MapExtractor;
}

namespace EHL
{

  class Base : public ADAPTER::AlgorithmBase
  {
  public:

    /// create using a Epetra_Comm
    explicit Base(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& lubricationparams,
        const Teuchos::ParameterList& structparams,
        const std::string struct_disname,
        const std::string lubrication_disname);// Problem builder

    /// virtual destructor
    virtual ~Base(){};

    /// setup
    virtual void SetupSystem() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// test results (if necessary)
    void TestResults(const Epetra_Comm& comm);

    /// read restart
    void ReadRestart( int restart );

    /// read restart from given time
    void ReadRestartfromTime( double restarttime );

    //! access to structural field
    const Teuchos::RCP<ADAPTER::Structure>& StructureField(){return structure_;}

    //! access to lubrication field
//    const Teuchos::RCP<ADAPTER::LubricationBaseAlgorithm>&  LubricationField(){return lubrication_;}

    /// set structure solution on lubrication field
    void SetStructSolution( Teuchos::RCP<const Epetra_Vector> disp,
                            Teuchos::RCP<const Epetra_Vector> vel );

    /// set lubrication solution on structure field
    void SetLubricationSolution( Teuchos::RCP<const Epetra_Vector> quantity );

  protected:

    /// underlying structure of the EHL problem
    Teuchos::RCP< ::ADAPTER::Structure>    structure_;

    /// underlying lubrication problem of the EHL problem
    Teuchos::RCP<ADAPTER::LubricationBaseAlgorithm> lubrication_;

    //! Type of coupling strategy between the two fields of the EHL problems
    const INPAR::EHL::FieldCoupling fieldcoupling_;

    /// setup adapters for EHL on boundary
    void SetupFieldCoupling(const std::string struct_disname, const std::string lubrication_disname);

    //! take current results for converged and save for next time step
    virtual void Update();

    //! write output
    virtual void Output(bool forced_writerestart = false);

    /// velocity calculation given the displacements
    Teuchos::RCP<Epetra_Vector> CalcVelocity( Teuchos::RCP<const Epetra_Vector> dispnp );

  private:
    /// setup discretizations and dofsets
    void SetupDiscretizations(const Epetra_Comm& comm, const std::string struct_disname, const std::string lubrication_disname);

    /// set structure mesh displacement on lubrication field
    void SetMeshDisp( Teuchos::RCP<const Epetra_Vector> disp ) ;

    /// set structure velocity fields on lubrication field
    void SetVelocityFields( Teuchos::RCP<const Epetra_Vector> quantity ) ;

    /// Function for checking that the different time steps are a
    /// multiplicative of each other
    //void CheckTimeStepping(double dt1, double dt2);
  };
}


#endif /* SRC_DRT_EHL_EHL_BASE_H_ */
