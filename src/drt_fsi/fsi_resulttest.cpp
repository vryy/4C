/*----------------------------------------------------------------------*/
/*!
\file fsi_resulttest.cpp

\brief testing of FSI calculation results

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/
/*----------------------------------------------------------------------*/

#include <string>
#include "fsi_resulttest.H"
#include "fsi_monolithicstructuresplit.H"
#include "fsi_monolithicfluidsplit.H"
#include "fsi_mortarmonolithic_structuresplit.H"
#include "fsi_mortarmonolithic_fluidsplit.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FSIResultTest::FSIResultTest(Teuchos::RCP<FSI::Monolithic> fsi,
                                  const Teuchos::ParameterList& fsidyn)
  : DRT::ResultTest("FSI")
{
  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    {
      const Teuchos::RCP<FSI::MonolithicFluidSplit>& fsiobject
        = Teuchos::rcp_dynamic_cast<FSI::MonolithicFluidSplit>(fsi);

      if (fsiobject == Teuchos::null)
        dserror("Cast to FSI::MonolithicFluidSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->FluidField().Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_monolithicstructuresplit:
    {
      const Teuchos::RCP<FSI::MonolithicStructureSplit>& fsiobject
        = Teuchos::rcp_dynamic_cast<FSI::MonolithicStructureSplit>(fsi);

      if (fsiobject == Teuchos::null)
        dserror("Cast to FSI::MonolithicStructureSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->StructureField()->Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_mortar_monolithicfluidsplit:
    {
      const Teuchos::RCP<FSI::MortarMonolithicFluidSplit>& fsiobject
        = Teuchos::rcp_dynamic_cast<FSI::MortarMonolithicFluidSplit>(fsi);

      if (fsiobject == Teuchos::null)
        dserror("Cast to FSI::MortarMonolithicFluidSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->FluidField().Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    case fsi_iter_mortar_monolithicstructuresplit:
    {
      const Teuchos::RCP<FSI::MortarMonolithicStructureSplit>& fsiobject
        = Teuchos::rcp_dynamic_cast<FSI::MortarMonolithicStructureSplit>(fsi);

      if (fsiobject == Teuchos::null)
        dserror("Cast to FSI::MortarMonolithicStructureSplit failed.");

      // Lagrange multipliers live on the slave field
      slavedisc_ = fsiobject->StructureField()->Discretization();
      fsilambda_ = fsiobject->lambda_;

      break;
    }
    default:
    {
      slavedisc_ = Teuchos::null;
      fsilambda_ = Teuchos::null;

      cout << "\nNo FSI test routines implemented for this coupling algorithm." << endl;

      break;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FSIResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(slavedisc_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  slavedisc_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node+1,slavedisc_->Name().c_str());
  }
  else
  {
    if (slavedisc_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = slavedisc_->gNode(node);

      // Strange! It seems we might actually have a global node around
      // even if it does not belong to us. But here we are just
      // interested in our nodes!
      if (actnode->Owner() != slavedisc_->Comm().MyPID())
        return;

      std::string position;
      res.ExtractString("QUANTITY",position);
      bool unknownpos = true; // make sure the result value string can be handled
      double result = 0.0;    // will hold the actual result of run

      // test Lagrange multipliers
      if (fsilambda_ != Teuchos::null)
      {
        const Epetra_BlockMap& fsilambdamap = fsilambda_->Map();
        if (position=="lambdax")
        {
          unknownpos = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->Dof(0,actnode,0))];
        }
        else if (position=="lambday")
        {
          unknownpos = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->Dof(0,actnode,1))];
        }
        else if (position=="lambdaz")
        {
          unknownpos = false;
          result = (*fsilambda_)[fsilambdamap.LID(slavedisc_->Dof(0,actnode,2))];
        }
      }
      // catch position strings, which are not handled by fsi result test
      if (unknownpos)
        dserror("Quantity '%s' not supported in fsi testing", position.c_str());

      // compare values
      const int err = CompareValues(result, res);
      nerr += err;
      test_count++;
    }
  }
}
