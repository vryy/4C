/*----------------------------------------------------------------------*/
/*!
\file fsi_structureale.cpp

\brief Solve structure only using FSI framework

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "fsi_structureale.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

#include "../drt_lib/drt_colors.H"

#include <string>
#include <Epetra_Time.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::StructureALE::StructureALE(Epetra_Comm& comm)
  : Algorithm(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  ADAPTER::Coupling& coupsf = StructureFluidCoupling();

  if (Teuchos::getIntegralValue<int>(fsidyn,"COUPMETHOD"))
  {
    matchingnodes_ = true;
    coupsf.SetupConditionCoupling(*StructureField().Discretization(),
                                   StructureField().Interface(),
                                  *MBFluidField().Discretization(),
                                   MBFluidField().Interface(),
                                  "FSICoupling");

    // In the following we assume that both couplings find the same dof
    // map at the structural side. This enables us to use just one
    // interface dof map for all fields and have just one transfer
    // operator from the interface map to the full field map.
//     if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
//       dserror("structure interface dof maps do not match");

    if (coupsf.MasterDofMap()->NumGlobalElements()==0)
    {
      dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");
    }
  }
  else
  {
    matchingnodes_ = false;
    coupsfm_.Setup( *StructureField().Discretization(),
                    *MBFluidField().Discretization(),
                    comm );
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureALE::Timeloop()
{
  while (NotFinished())
  {
    PrepareTimeStep();
    Solve();
    Update();
    Output();
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::StructureALE::Solve()
{
  StructureField().Solve();
  MBFluidField().NonlinearSolve(StructToFluid(StructureField().ExtractInterfaceDispnp()));
}


Teuchos::RCP<Epetra_Vector> FSI::StructureALE::StructToFluid(Teuchos::RCP<Epetra_Vector> iv)
{
  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  if (matchingnodes_)
  {
    return coupsf.MasterToSlave(iv);
  }
  else
  {
    return coupsfm_.MasterToSlave(iv);
  }
}

#endif
