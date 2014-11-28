
/*----------------------------------------------------------------------*/
/*!
\file  xfluidfluidresulttest.cpp

\brief tesing of fluid calculation results

<pre>
Maintainer: Shadan Shahmiri
            shahmiri@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*/
/*----------------------------------------------------------------------*/


#include <string>

#include "xfluidfluidresulttest.H"
#include "xfluidfluid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::XFluidFluidResultTest::XFluidFluidResultTest(XFluidFluid& fluid)
  : DRT::ResultTest("FLUID")
{
  embfluiddis_= fluid.embdis_;
  bgfluiddis_= fluid.bgdis_;
  embfluidsol_   = fluid.velnp_ ;
  bgfluidsol_   = fluid.state_->velnp_ ;
//  mytraction_ = fluid.calcStresses();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::XFluidFluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if ((dis != bgfluiddis_->Name()) and (dis != embfluiddis_->Name()))
    return;

   int node;
   res.ExtractInt("NODE",node);
   //node -= 1;
   //ToDo  activate the line above

   if (dis == bgfluiddis_->Name())
   {
     int havenode(bgfluiddis_->HaveGlobalNode(node));
     int isnodeofanybody(0);
     bgfluiddis_->Comm().SumAll(&havenode,&isnodeofanybody,1);

     if (isnodeofanybody==0)
     {
       dserror("Node %d does not belong to discretization %s",node+1,bgfluiddis_->Name().c_str());
     }

     if (bgfluiddis_->HaveGlobalNode(node))
     {
         const DRT::Node* actnode = bgfluiddis_->gNode(node);

         // Test only, if actnode is a row node
         if (actnode->Owner() != bgfluiddis_->Comm().MyPID())
             return;

         double result = 0.;

         const Epetra_BlockMap& velnpmap = bgfluidsol_->Map();

         const int numdim = DRT::Problem::Instance()->NDim();

     std::string position;
     res.ExtractString("QUANTITY",position);
         if (position=="velx")
         {
             result = (*bgfluidsol_)[velnpmap.LID(bgfluiddis_->Dof(actnode,0))];
         }
         else if (position=="vely")
         {
             result = (*bgfluidsol_)[velnpmap.LID(bgfluiddis_->Dof(actnode,1))];
         }
         else if (position=="velz")
         {
             if (numdim==2)
                 dserror("Cannot test result for velz in 2D case.");
             result = (*bgfluidsol_)[velnpmap.LID(bgfluiddis_->Dof(actnode,2))];
         }
         else if (position=="pressure")
         {
             if (numdim==2)
             {
                 if (bgfluiddis_->NumDof(actnode)<3)
                     dserror("too few dofs at node %d for pressure testing",actnode->Id());
                 result = (*bgfluidsol_)[velnpmap.LID(bgfluiddis_->Dof(actnode,2))];
             }
             else
             {
                 if (bgfluiddis_->NumDof(actnode)<4)
                     dserror("too few dofs at node %d for pressure testing",actnode->Id());
                 result = (*bgfluidsol_)[velnpmap.LID(bgfluiddis_->Dof(actnode,3))];
             }
         }
//          else if (position=="tractionx")
//              result = (*mytraction_)[(mytraction_->Map()).LID(bgfluiddis_->Dof(actnode,0))];
//          else if (position=="tractiony")
//              result = (*mytraction_)[(mytraction_->Map()).LID(bgfluiddis_->Dof(actnode,1))];
//          else if (position=="tractionz")
//          {
//              if (numdim==2)
//                  dserror("Cannot test result for tractionz in 2D case.");
//              result = (*mytraction_)[(mytraction_->Map()).LID(bgfluiddis_->Dof(actnode,2))];
//          }
         else
         {
             dserror("Quantity '%s' not supported in fluid testing", position.c_str());
         }

         nerr += CompareValues(result, "NODE", res);
         test_count++;
     }
   }
   else if(dis == embfluiddis_->Name())
   {
     int havenode(embfluiddis_->HaveGlobalNode(node));
     int isnodeofanybody(0);
     embfluiddis_->Comm().SumAll(&havenode,&isnodeofanybody,1);

     if (isnodeofanybody==0)
     {
       dserror("Node %d does not belong to discretization %s",node+1,embfluiddis_->Name().c_str());
     }

     if (embfluiddis_->HaveGlobalNode(node))
     {
       const DRT::Node* actnode = embfluiddis_->gNode(node);

       // Test only, if actnode is a row node
       if (actnode->Owner() != embfluiddis_->Comm().MyPID())
         return;

       double result = 0.;

       const Epetra_BlockMap& velnpmap = embfluidsol_->Map();

       const int numdim = DRT::Problem::Instance()->NDim();

       std::string position;
       res.ExtractString("QUANTITY",position);
       if (position=="velx")
       {
         result = (*embfluidsol_)[velnpmap.LID(embfluiddis_->Dof(actnode,0))];
       }
       else if (position=="vely")
       {
         result = (*embfluidsol_)[velnpmap.LID(embfluiddis_->Dof(actnode,1))];
       }
       else if (position=="velz")
       {
         if (numdim==2)
           dserror("Cannot test result for velz in 2D case.");
         result = (*embfluidsol_)[velnpmap.LID(embfluiddis_->Dof(actnode,2))];
       }
       else if (position=="pressure")
       {
         if (numdim==2)
         {
           if (embfluiddis_->NumDof(actnode)<3)
             dserror("too few dofs at node %d for pressure testing",actnode->Id());
           result = (*embfluidsol_)[velnpmap.LID(embfluiddis_->Dof(actnode,2))];
         }
         else
         {
           if (embfluiddis_->NumDof(actnode)<4)
             dserror("too few dofs at node %d for pressure testing",actnode->Id());
           result = (*embfluidsol_)[velnpmap.LID(embfluiddis_->Dof(actnode,3))];
         }
       }
//          else if (position=="tractionx")
//             result = (*mytraction_)[(mytraction_->Map()).LID(xfluiddis_->Dof(actnode,0))];
//          else if (position=="tractiony")
//              result = (*mytraction_)[(mytraction_->Map()).LID(xfluiddis_->Dof(actnode,1))];
//          else if (position=="tractionz")
//          {
//              if (numdim==2)
//                  dserror("Cannot test result for tractionz in 2D case.");
//              result = (*mytraction_)[(mytraction_->Map()).LID(xfluiddis_->Dof(actnode,2))];
//          }
         else
         {
             dserror("Quantity '%s' not supported in fluid testing", position.c_str());
         }

         nerr += CompareValues(result, "NODE", res);
         test_count++;
     }
   }
}
