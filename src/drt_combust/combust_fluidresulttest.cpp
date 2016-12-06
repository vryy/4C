/*----------------------------------------------------------------------*/
/*!
\file combust_fluidresulttest.cpp

\brief tesing of fluid calculation results

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/

#include "combust_fluidimplicitintegration.H"
#include "combust_fluidresulttest.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dofset_independent_pbc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_xfem/dof_management.H"
#include "../linalg/linalg_blocksparsematrix.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::CombustFluidResultTest::CombustFluidResultTest(CombustFluidImplicitTimeInt& fluid)
  : DRT::ResultTest("FLUID")
{
  fluiddis_ = fluid.discret_;
  fluidstddofset_ = fluid.standarddofset_;

  std::set<XFEM::PHYSICS::Field> outputfields;
  outputfields.insert(XFEM::PHYSICS::Velx);
  outputfields.insert(XFEM::PHYSICS::Vely);
  outputfields.insert(XFEM::PHYSICS::Velz);
  outputfields.insert(XFEM::PHYSICS::Pres);

  // transform XFEM velnp vector to (standard FEM) velnp vector
  Teuchos::RCP<Epetra_Vector> velnp_std = fluid.dofmanagerForOutput_->transformXFEMtoStandardVector(
      *fluid.state_.velnp_, *fluid.standarddofset_, fluid.state_.nodalDofDistributionMap_, outputfields);
  mysol_ = velnp_std;

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::CombustFluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != fluiddis_->Name())
    return;

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(fluiddis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  fluiddis_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretization %s",node+1,fluiddis_->Name().c_str());
  }
  else
  {
    if (fluiddis_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = fluiddis_->gNode(node);

      // Test only, if actnode is a row node
      if (actnode->Owner() != fluiddis_->Comm().MyPID())
        return;

      double result = 0.;

      // get map of the standard fluid dofset (no XFEM dofs)
      const Epetra_Map& velnpmap = *fluidstddofset_->DofRowMap();

      const int numdim = DRT::Problem::Instance()->NDim();

      std::string position;
      res.ExtractString("QUANTITY",position);
      if (position=="velx")
      {
        result = (*mysol_)[velnpmap.LID(fluidstddofset_->Dof(actnode,0))];
      }
      else if (position=="vely")
      {
        result = (*mysol_)[velnpmap.LID(fluidstddofset_->Dof(actnode,1))];
      }
      else if (position=="velz")
      {
        if (numdim==2)
          dserror("Cannot test result for velz in 2D case.");
        result = (*mysol_)[velnpmap.LID(fluidstddofset_->Dof(actnode,2))];
      }
      else if (position=="pressure")
      {
        if (numdim==2)
        {
          if (fluiddis_->NumDof(actnode)<3)
            dserror("too few dofs at node %d for pressure testing",actnode->Id());
          result = (*mysol_)[velnpmap.LID(fluidstddofset_->Dof(actnode,2))];
        }
        else
        {
          if (fluiddis_->NumDof(actnode)<4)
            dserror("too few dofs at node %d for pressure testing",actnode->Id());
          result = (*mysol_)[velnpmap.LID(fluidstddofset_->Dof(actnode,3))];
        }
      }
      else
      {
        dserror("Quantity '%s' not supported in fluid testing", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}
