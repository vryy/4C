/*!----------------------------------------------------------------------
\file immersed_discretization.cpp

\brief immersed specific background discretization

</pre>

<pre>
\level 1

\maintainer Andreas Rauch
</pre>

*----------------------------------------------------------------------*/

#include "immersed_discretization.H"
#include "immersed_node.H"
#include "immersed_field_exchange_manager.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_fluid_ele/fluid_ele_poro_immersed.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                    rauch 05/17 |
 *----------------------------------------------------------------------*/
DRT::DiscretizationImmersed::DiscretizationImmersed(
    const std::string name,
    Teuchos::RCP<Epetra_Comm> comm)
  :
  Discretization(name, comm) // use base class constructor
{};


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DiscretizationImmersed::Evaluate(
                        Teuchos::ParameterList& params,
                        DRT::AssembleStrategy&  strategy )
{
  bool selectiveassembly = false;
  DRT::ImmersedFieldExchangeManager* exchangemanager =
      DRT::ImmersedFieldExchangeManager::Instance();

  if(not exchangemanager->GetPseudoBoundarySwitch())
  {
    DRT::Discretization::Evaluate(params,strategy);
  }
  else
  {
    std::vector<int>* immersedeles = exchangemanager -> GetEles();

    int mysize = immersedeles->size();
    int maxsize = -1;
    Comm().MaxAll(&mysize,&maxsize,1);

    if (maxsize == 0) dserror("selective assembly requested but no elements provided");
    if (!Filled()) dserror("FillComplete() was not called");
    if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    int row = strategy.FirstDofSet();
    int col = strategy.SecondDofSet();

    // call the element's register class preevaluation method
    // for each type of element
    // for most element types, just the base class dummy is called
    // that does nothing
    ParObjectFactory::Instance().PreEvaluate(
        *this,
        params,
        strategy.Systemmatrix1(),
        strategy.Systemmatrix2(),
        strategy.Systemvector1(),
        strategy.Systemvector2(),
        strategy.Systemvector3());

    Element::LocationArray la(dofsets_.size());
    Element::LocationArray la_selective(dofsets_.size());

    // loop over column elements
    const int numcolele = NumMyColElements();
    for (int i=0; i<numcolele; ++i)
    {
      selectiveassembly = false;
      DRT::Element* actele = lColElement(i);

      // get element location vector, dirichlet flags and ownerships
      actele->LocationVector(*this,la,false);

      // get dimension of element matrices and vectors
      // Reshape element matrices and vectors and init to zero
      strategy.ClearElementStorage( la[row].Size(), la[col].Size() );


      // call the element evaluate method
      int err = actele->Evaluate(
          params,
          *this,
          la,
          strategy.Elematrix1(),
          strategy.Elematrix2(),
          strategy.Elevector1(),
          strategy.Elevector2(),
          strategy.Elevector3());
      if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);

      if( exchangemanager->GetPseudoBoundarySwitch() )
      {
        if( immersedeles->size()>0 )
        {
          for(int iele=0;iele<(int)immersedeles->size();iele++)
          {
            if(actele->Id()==immersedeles->at(iele))
            {
              selectiveassembly = true;
              break;
            }
          }

          if(selectiveassembly)
          {
            LocationVectorSelectiveAssembly(actele,*this,la_selective);
          }
        }
      } // if selective assembly is switched on externally

      if(selectiveassembly == false)
      {
        int eid = actele->Id();
        strategy.AssembleMatrix1( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
        strategy.AssembleMatrix2( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
        strategy.AssembleVector1( la[row].lm_, la[row].lmowner_ );
        strategy.AssembleVector2( la[row].lm_, la[row].lmowner_ );
        strategy.AssembleVector3( la[row].lm_, la[row].lmowner_ );
      }
      else
      {
        int eid = actele->Id();
        strategy.AssembleMatrix1( eid, la_selective[row].lm_, la[col].lm_, la_selective[row].lmowner_, la[col].stride_ );
        strategy.AssembleMatrix2( eid, la_selective[row].lm_, la[col].lm_, la_selective[row].lmowner_, la[col].stride_ );
        strategy.AssembleVector1( la_selective[row].lm_, la_selective[row].lmowner_ );
        strategy.AssembleVector2( la_selective[row].lm_, la_selective[row].lmowner_ );
        strategy.AssembleVector3( la_selective[row].lm_, la_selective[row].lmowner_ );
      }

    } // for (int i=0; i<numcolele; ++i)
  } // if selective assembly requested

  return;
}// DRT::DiscretizationImmersed::Evaluate


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DiscretizationImmersed::EvaluateCondition
(
  Teuchos::ParameterList& params,
  DRT::AssembleStrategy & strategy,
  const std::string& condstring,
  const int condid
)
{
  bool selectiveassembly = false;
  DRT::ImmersedFieldExchangeManager* exchangemanager =
      DRT::ImmersedFieldExchangeManager::Instance();

  if(not exchangemanager->GetPseudoBoundarySwitch())
  {
    DRT::Discretization::EvaluateCondition(params,strategy,condstring,condid);
  }
  else
  {
    std::vector<int>* immersedeles = exchangemanager -> GetEles();

    int mysize = immersedeles->size();
    int maxsize = -1;
    Comm().MaxAll(&mysize,&maxsize,1);

    if (maxsize == 0) dserror("selective assembly requested but no elements provided");
    if (!Filled()) dserror("FillComplete() was not called");
    if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

    int row = strategy.FirstDofSet();
    int col = strategy.SecondDofSet();

    // get the current time
    bool usetime = true;
    const double time = params.get("total time",-1.0);
    if (time<0.0) usetime = false;

    Element::LocationArray la(dofsets_.size());
    Element::LocationArray la_selective(dofsets_.size());


    std::multimap<std::string,Teuchos::RCP<Condition> >::iterator fool;

    //----------------------------------------------------------------------
    // loop through conditions and evaluate them if they match the criterion
    //----------------------------------------------------------------------
    for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
    {
      if (fool->first == condstring)
      {
        DRT::Condition& cond = *(fool->second);
        if (condid == -1 || condid ==cond.GetInt("ConditionID"))
        {
          std::map<int,Teuchos::RCP<DRT::Element> >& geom = cond.Geometry();
          // if (geom.empty()) dserror("evaluation of condition with empty geometry");
          // no check for empty geometry here since in parallel computations
          // can exist processors which do not own a portion of the elements belonging
          // to the condition geometry
          std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;

          // Evaluate Loadcurve if defined. Put current load factor in parameter list
          const std::vector<int>*    curve  = cond.Get<std::vector<int> >("curve");
          int curvenum = -1;
          if (curve) curvenum = (*curve)[0];
          double curvefac = 1.0;
          if (curvenum>=0 && usetime)
            curvefac = Problem::Instance()->Funct(curvenum).EvaluateTime(time);

          // Get ConditionID of current condition if defined and write value in parameter list
          const std::vector<int>*    CondIDVec  = cond.Get<std::vector<int> >("ConditionID");
          if (CondIDVec)
          {
            params.set("ConditionID",(*CondIDVec)[0]);
            char factorname[30];
            sprintf(factorname,"LoadCurveFactor %d",(*CondIDVec)[0]);
            params.set(factorname,curvefac);
          }
          else
          {
            params.set("LoadCurveFactor",curvefac);
          }
          params.set<Teuchos::RCP<DRT::Condition> >("condition", fool->second);

          for (curr=geom.begin(); curr!=geom.end(); ++curr)
          {
            selectiveassembly = false;

            // get element location vector and ownerships
            // the LocationVector method will return the the location vector
            // of the dofs this condition is meant to assemble into.
            // These dofs do not need to be the same as the dofs of the element
            // (this is the standard case, though). Special boundary conditions,
            // like weak Dirichlet conditions, assemble into the dofs of the parent element.
            curr->second->LocationVector(*this,la,false,condstring,params);

            // get dimension of element matrices and vectors
            // Reshape element matrices and vectors and initialize to zero
            strategy.ClearElementStorage( la[row].Size(), la[col].Size() );

            // call the element specific evaluate method
            int err = curr->second->Evaluate(params,*this,la,
                strategy.Elematrix1(),
                strategy.Elematrix2(),
                strategy.Elevector1(),
                strategy.Elevector2(),
                strategy.Elevector3());
            if (err) dserror("error while evaluating elements");

            if( exchangemanager->GetPseudoBoundarySwitch() )
            {
              if( immersedeles->size()>0 )
              {
                for(int iele=0;iele<(int)immersedeles->size();iele++)
                {
                  if(curr->second->Id()==immersedeles->at(iele))
                  {
                    selectiveassembly = true;
                    break;
                  }
                }

                if(selectiveassembly)
                {
                  LocationVectorSelectiveAssembly(curr->second.getRawPtr(),*this,la_selective);
                }
              }
            } // if selective assembly is switched on externally

            // assembly
            /* If BlockMatrixes are used, the decision which assemble strategy is
             * used, is based on the element id. As this id is compared to a list
             * of conditioned volume elements, always the volume element id should
             * be given to the Assembling! (comment: eid is not used by
             * sysmat.assemble(...,eid,...))*/
            int eid;
            if (DRT::FaceElement* faceele = dynamic_cast<DRT::FaceElement*>(curr->second.get()))
              eid = faceele->ParentElement()->Id();
            else
              eid = curr->second->Id();

            if(selectiveassembly == false)
            {
              strategy.AssembleMatrix1( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
              strategy.AssembleMatrix2( eid, la[row].lm_, la[col].lm_, la[row].lmowner_, la[col].stride_ );
              strategy.AssembleVector1( la[row].lm_, la[row].lmowner_ );
              strategy.AssembleVector2( la[row].lm_, la[row].lmowner_ );
              strategy.AssembleVector3( la[row].lm_, la[row].lmowner_ );
            }
            else
            {
              strategy.AssembleMatrix1( eid, la_selective[row].lm_, la[col].lm_, la_selective[row].lmowner_, la[col].stride_ );
              strategy.AssembleMatrix2( eid, la_selective[row].lm_, la[col].lm_, la_selective[row].lmowner_, la[col].stride_ );
              strategy.AssembleVector1( la_selective[row].lm_, la_selective[row].lmowner_ );
              strategy.AssembleVector2( la_selective[row].lm_, la_selective[row].lmowner_ );
              strategy.AssembleVector3( la_selective[row].lm_, la_selective[row].lmowner_ );
            }
          }
        }
      }
    } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  }
  return;
}// DRT::DiscretizationImmersed::EvaluateCondition


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::DiscretizationImmersed::LocationVectorSelectiveAssembly(
    const Element* ele,
    const Discretization& dis,
    DRT::Element::LocationArray& la) const
{
  const int numnode = ele->NumNode();
  const DRT::Node*const* nodes = ele->Nodes();

  la.Clear();

  // we need to look at all DofSets of our Discretization
  for (int dofset=0; dofset<la.Size(); ++dofset)
  {
    std::vector<int>& lm  = la[dofset].lm_;
    std::vector<int>& lmowner  = la[dofset].lmowner_;
    std::vector<int>& lmstride = la[dofset].stride_;

    // fill the vector with nodal dofs
    if (nodes)
    {
      for (int i=0; i<numnode; ++i)
      {
        const DRT::Node* node = nodes[i];

        const int owner = node->Owner();
        std::vector<int> dof;
        dis.Dof(dof,node,dofset,0,ele);

        // if there are more dofs on the node than the element can handle, this cannot work
        dsassert(
            ele->NumDofPerNode(*node) <= (int) dof.size() or dofset != 0,
            "More dofs on node than element can handle! Internal error!"
        );

        // assume that the first dofs are the relevant ones
        const int size = dofset == 0 ? ele->NumDofPerNode(*node) : dof.size();

        IMMERSED::ImmersedNode* immersednode = NULL;
        DRT::ELEMENTS::FluidPoroImmersed* immersedele = NULL;

        Teuchos::RCP<DRT::Discretization> immersedinfodis =
            DRT::Problem::Instance()->GetDis("porofluid");
        DRT::Element* iele = immersedinfodis->gElement(ele->Id());
        immersedele =
            dynamic_cast<DRT::ELEMENTS::FluidPoroImmersed* >(iele);
        if(immersedele==NULL)
          dserror("dynamic cast to FluidPoroImmersed* failed");

        int structnodeid = node->Id();
        immersednode =
            dynamic_cast<IMMERSED::ImmersedNode* >(immersedinfodis->gNode(structnodeid));
        if(immersednode==NULL)
          dserror("dynamic cast to ImmersedNode* failed");

        if (size) lmstride.push_back(1);

        for (int j=0; j<size; ++j)
        {
          if(immersednode->IsBoundaryImmersed()==false)
          {
            // do not touch last dof at each node
            if(j<(size-1))
            {
              lmowner.push_back(owner);
              lm.push_back(dof[j]);
            }
          }
        }
      } // end loop over nodes of this element

      // fill the vector with element dofs
      const int owner = ele->Owner();
      std::vector<int> dof = dis.Dof(dofset,ele);
      if (dof.size()) lmstride.push_back(dof.size());
      for (unsigned j=0; j<dof.size(); ++j)
      {
        lmowner.push_back(owner);
        lm.push_back(dof[j]);
      }

    } // if the element has nodes
  } // for (int dofset=0; dofset<la.Size(); ++dofset)

  return;
}// DRT::DiscretizationImmersed::LocationVectorSelectiveAssembly
