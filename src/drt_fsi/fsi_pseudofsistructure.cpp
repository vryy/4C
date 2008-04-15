
#ifdef CCADISCRET

#include "fsi_pseudofsistructure.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

#include "../drt_lib/drt_colors.H"

#include <string>
#include <Epetra_Time.h>
#include <Teuchos_StandardParameterEntryValidators.hpp>


/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::PseudoFSIStructure::PseudoFSIStructure(Epetra_Comm& comm)
  : Algorithm(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  FSI::Coupling& coupsf = StructureFluidCoupling();

  if (Teuchos::getIntegralValue<int>(fsidyn,"COUPMETHOD"))
  {
    matchingnodes_ = true;
    coupsf.SetupConditionCoupling(*StructureField().Discretization(),
                                   StructureField().Interface(),
                                  *FluidField().Discretization(),
                                   FluidField().Interface(),
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
                    *FluidField().Discretization(),
                    comm );
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PseudoFSIStructure::Timeloop()
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
void FSI::PseudoFSIStructure::Solve()
{
  StructureField().Solve();  
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::PseudoFSIStructure::Update()
{
  StructureField().Update();
}

#endif
