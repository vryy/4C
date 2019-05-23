/*----------------------------------------------------------------------*/
/*!

\brief testing of FSI calculation results

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#include <string>
#include "fsi_resulttest.H"
#include "fsi_monolithicstructuresplit.H"
#include "fsi_monolithicfluidsplit.H"
#include "fsi_mortarmonolithic_structuresplit.H"
#include "fsi_mortarmonolithic_fluidsplit.H"
#include "fsi_fluidfluidmonolithic_structuresplit_nonox.H"
#include "fsi_fluidfluidmonolithic_fluidsplit_nonox.H"
#include "fsi_slidingmonolithic_fluidsplit.H"
#include "fsi_slidingmonolithic_structuresplit.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FSIResultTest::FSIResultTest(
    Teuchos::RCP<FSI::Monolithic>& fsi, const Teuchos::ParameterList& fsidyn)
    : DRT::ResultTest("FSI"), fsi_(fsi)
{
  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_fluidfluid_monolithicfluidsplit:
    {
      const Teuchos::RCP<FSI::MonolithicFluidSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MonolithicFluidSplit>(fsi);

      if (fsiobject == Teuchos::null) dserror("Cast to FSI::MonolithicFluidSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->FluidField()->Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_monolithicstructuresplit:
    case fsi_iter_fluidfluid_monolithicstructuresplit:
    {
      const Teuchos::RCP<FSI::MonolithicStructureSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MonolithicStructureSplit>(fsi);

      if (fsiobject == Teuchos::null) dserror("Cast to FSI::MonolithicStructureSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->StructureField()->Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_mortar_monolithicfluidsplit:
    {
      const Teuchos::RCP<FSI::MortarMonolithicFluidSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MortarMonolithicFluidSplit>(fsi);

      if (fsiobject == Teuchos::null) dserror("Cast to FSI::MortarMonolithicFluidSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->FluidField()->Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_mortar_monolithicstructuresplit:
    {
      const Teuchos::RCP<FSI::MortarMonolithicStructureSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::MortarMonolithicStructureSplit>(fsi);

      if (fsiobject == Teuchos::null)
        dserror("Cast to FSI::MortarMonolithicStructureSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->StructureField()->Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_sliding_monolithicfluidsplit:
    {
      const Teuchos::RCP<FSI::SlidingMonolithicFluidSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::SlidingMonolithicFluidSplit>(fsi);

      if (fsiobject == Teuchos::null) dserror("Cast to FSI::SlidingMonolithicFluidSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->FluidField()->Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_sliding_monolithicstructuresplit:
    {
      const Teuchos::RCP<FSI::SlidingMonolithicStructureSplit>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::SlidingMonolithicStructureSplit>(fsi);

      if (fsiobject == Teuchos::null)
        dserror("Cast to FSI::SlidingMonolithicStructureSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->StructureField()->Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    default:
    {
      slavedisc_ = Teuchos::null;
      fsilambda_ = Teuchos::null;

      std::cout << "\nNo FSI test routines implemented for this coupling algorithm." << std::endl;

      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FSIResultTest::FSIResultTest(
    Teuchos::RCP<FSI::MonolithicNoNOX> fsi, const Teuchos::ParameterList& fsidyn)
    : DRT::ResultTest("FSI")
{
  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_fluidfluid_monolithicstructuresplit_nonox:
    {
      // Lagrange multipliers live on the slave field
      slavedisc_ = fsi->StructureField()->Discretization();

      const Teuchos::RCP<FSI::FluidFluidMonolithicStructureSplitNoNOX>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::FluidFluidMonolithicStructureSplitNoNOX>(fsi);

      if (fsiobject == Teuchos::null)
        dserror("Cast to FSI::FluidFluidMonolithicStructureSplitNoNOX failed.");
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_fluidfluid_monolithicfluidsplit_nonox:
    {
      // Lagrange multiplier lives on the slave field (fluid in this case!)
      slavedisc_ = fsi->FluidField()->Discretization();

      const Teuchos::RCP<FSI::FluidFluidMonolithicFluidSplitNoNOX>& fsiobject =
          Teuchos::rcp_dynamic_cast<FSI::FluidFluidMonolithicFluidSplitNoNOX>(fsi);

      if (fsiobject == Teuchos::null)
        dserror("Cast to FSI::FluidFluidMonolithicFluidSplitNoNOX failed.");

      fsilambda_ = fsiobject->lambda_;

      break;
    }
    default:
    {
      slavedisc_ = Teuchos::null;
      fsilambda_ = Teuchos::null;

      std::cout << "\nNo FSI test routines implemented for this coupling algorithm." << std::endl;

      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FSIResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(slavedisc_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  slavedisc_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", node + 1, slavedisc_->Name().c_str());
  }
  else
  {
    if (slavedisc_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = slavedisc_->gNode(node);

      // Strange! It seems we might actually have a global node around
      // even if it does not belong to us. But here we are just
      // interested in our nodes!
      if (actnode->Owner() != slavedisc_->Comm().MyPID()) return;

      std::string quantity;
      res.ExtractString("QUANTITY", quantity);
      bool unknownquantity = true;  // make sure the result value std::string can be handled
      double result = 0.0;          // will hold the actual result of run

      // test Lagrange multipliers
      if (fsilambda_ != Teuchos::null)
      {
        const Epetra_BlockMap& fsilambdamap = fsilambda_->Map();
        if (quantity == "lambdax")
        {
          unknownquantity = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->Dof(0, actnode, 0))];
        }
        else if (quantity == "lambday")
        {
          unknownquantity = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->Dof(0, actnode, 1))];
        }
        else if (quantity == "lambdaz")
        {
          unknownquantity = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->Dof(0, actnode, 2))];
        }
      }

      // catch quantity strings, which are not handled by fsi result test
      if (unknownquantity) dserror("Quantity '%s' not supported in fsi testing", quantity.c_str());

      // compare values
      const int err = CompareValues(result, "NODE", res);
      nerr += err;
      test_count++;
    }
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FSIResultTest::TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  dserror("FSI ELEMENT test not implemented, yet.");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FSIResultTest::TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  std::string quantity;
  res.ExtractString("QUANTITY", quantity);
  bool unknownquantity = true;  // make sure the result value std::string can be handled
  double result = 0.0;          // will hold the actual result of run

  // test for time step size
  if (quantity == "dt")
  {
    unknownquantity = false;
    result = fsi_->Dt();
  }

  // test for number of repetitions of time step in case of time step size adaptivity
  if (quantity == "adasteps")
  {
    unknownquantity = false;
    result = fsi_->GetNumAdaptSteps();
  }

  // test for simulation time in case of time step size adaptivity
  if (quantity == "time")
  {
    unknownquantity = false;
    result = fsi_->Time();
  }

  // catch quantity strings, which are not handled by fsi result test
  if (unknownquantity) dserror("Quantity '%s' not supported in fsi testing", quantity.c_str());

  // compare values
  const int err = CompareValues(result, "SPECIAL", res);
  nerr += err;
  test_count++;

  return;
}
