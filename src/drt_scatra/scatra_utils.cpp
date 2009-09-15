/*----------------------------------------------------------------------*/
/*!
\file scatra_utils.cpp

\brief utility functions for scalar transport problems

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_mat/matpar_material.H"
#include "scatra_utils.H"

// we need to know all necessary element types for the mesh creation
#include "../drt_f2/fluid2.H"
#include "../drt_f3/fluid3.H"
//#include "../drt_combust/combust3.H"
#include "../drt_scatra/scatra_element.H"


/*----------------------------------------------------------------------*/
//! create scalar transport discretization parallel to the fluid one
/*----------------------------------------------------------------------*/
void SCATRA::CreateScaTraDiscretization(
    Teuchos::RefCountPtr<DRT::Discretization>& fluiddis,
    Teuchos::RefCountPtr<DRT::Discretization>& scatradis,
    const std::map<string,string>& conditions_to_copy,
    int matid,
    const bool makequadratic
    )
{
  // is the fluid discretization ready?
//  if (!fluiddis->Filled() or !fluiddis->HaveDofs()) fluiddis->FillComplete();
  if (!fluiddis->Filled()) fluiddis->FillComplete(false,false,false);

  // is the second discretization really empty?
  if (scatradis->NumGlobalElements() or scatradis->NumGlobalNodes())
  {
    dserror("There are %d elements and %d nodes in empty discretization. Panic.",
            scatradis->NumGlobalElements(), scatradis->NumGlobalNodes());
  }

  // prepare some variables we need
  int myrank = scatradis->Comm().MyPID();

  vector<int> egid;
  egid.reserve(fluiddis->NumMyRowElements());

  vector<string> eletype;
  eletype.reserve(fluiddis->NumMyRowElements());

  set<int> rownodeset;
  set<int> colnodeset;
  const Epetra_Map* noderowmap = fluiddis->NodeRowMap();

  // Loop all fluid elements. Here we do ugly castings.
  // Please do not take this for an example of how to code in baci!

  // We need to test for all elements (including ghosted ones) to
  // catch all nodes attached to fluid elements
  int numelements = fluiddis->NumMyColElements();

  for (int i=0; i<numelements; ++i)
  {
    DRT::Element* actele = fluiddis->lColElement(i);
    bool ismyele = fluiddis->ElementRowMap()->MyGID(actele->Id());

    // we only support transport elements here
    eletype.push_back("TRANSP");
    {
      if (ismyele)
        egid.push_back(actele->Id());

      // copy node ids of actele to rownodeset but leave those that do
      // not belong to this processor
      remove_copy_if(actele->NodeIds(), actele->NodeIds()+actele->NumNode(),
                     inserter(rownodeset, rownodeset.begin()),
                     not1(DRT::UTILS::MyGID(noderowmap)));

      copy(actele->NodeIds(), actele->NodeIds()+actele->NumNode(),
           inserter(colnodeset, colnodeset.begin()));
    }
  } // loop over my elements

  if (makequadratic == true) dserror("Automatic conversion to quadratic elements not implemented.");

  // construct nodes in the new discretization
  for (int i=0; i<noderowmap->NumMyElements(); ++i)
  {
    int gid = noderowmap->GID(i);
    if (rownodeset.find(gid)!=rownodeset.end())
    {
      DRT::Node* fluidnode = fluiddis->lRowNode(i);
      scatradis->AddNode(rcp(new DRT::Node(gid, fluidnode->X(), myrank)));
    }
  }

  // we get the node maps almost for free
  vector<int> scatranoderowvec(rownodeset.begin(), rownodeset.end());
  rownodeset.clear();
  RefCountPtr<Epetra_Map> scatranoderowmap = rcp(new Epetra_Map(-1,
                                                             scatranoderowvec.size(),
                                                             &scatranoderowvec[0],
                                                             0,
                                                             scatradis->Comm()));
  scatranoderowvec.clear();

  vector<int> scatranodecolvec(colnodeset.begin(), colnodeset.end());
  colnodeset.clear();
  RefCountPtr<Epetra_Map> scatranodecolmap = rcp(new Epetra_Map(-1,
                                                             scatranodecolvec.size(),
                                                             &scatranodecolvec[0],
                                                             0,
                                                             scatradis->Comm()));
  scatranodecolvec.clear();

  // now do the elements

  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_scatra) &&
      (mtype != INPAR::MAT::m_mixfrac_scatra) &&
      (mtype != INPAR::MAT::m_sutherland_scatra) &&
      (mtype != INPAR::MAT::m_arrhenius_pv_scatra) &&
      (mtype != INPAR::MAT::m_matlist))
    dserror("Material with ID %d is not admissible for scalar transport elements",matid);

  // construct scalar transport elements
  // The order of the elements might be different from that of the
  // fluid elements. We don't care. There are no dofs to these elements.
  for (unsigned i=0; i<egid.size(); ++i)
  {
    DRT::Element* fluidele = fluiddis->gElement(egid[i]);

    // create a scalar transport element with the same global element id
    RCP<DRT::Element> scatraele = DRT::UTILS::Factory(eletype[i],"Polynomial",egid[i],myrank);

    // get global node ids of fluid element
    vector<int> nids;
    nids.reserve(fluidele->NumNode());
    transform(fluidele->Nodes(), fluidele->Nodes()+fluidele->NumNode(),
              back_inserter(nids), mem_fun(&DRT::Node::Id));

    // set the same global node ids to the new scalar transport element
    scatraele->SetNodeIds(nids.size(), &nids[0]);

    // We need to set material and gauss points to complete element setup.
    // This is again really ugly as we have to extract the actual
    // element type in order to access the material property
    // note: SetMaterial() was reimplemented by the transport element!
#if defined(D_FLUID2) || defined(D_FLUID3)
        DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(scatraele.get());
        if (trans!=NULL)
        {
          trans->SetMaterial(matid);
          trans->SetDisType(fluidele->Shape());
        }
        else
#endif
      {
        dserror("unsupported element type '%s'", typeid(*scatraele).name());
      }

    // add new scalar transport element to discretization
    scatradis->AddElement(scatraele);
  }

  // copy selected conditions to the new discretization (and rename them if desired)
  for (map<string,string>::const_iterator conditername = conditions_to_copy.begin();
       conditername != conditions_to_copy.end();
       ++conditername)
  {
    vector<DRT::Condition*> conds;
    fluiddis->GetCondition((*conditername).first, conds);
    for (unsigned i=0; i<conds.size(); ++i)
    {
      // We use the same nodal ids and therefore we can just copy the conditions.
      // The map gives the new condition names (e.g. renaming from TransportDirichlet to Dirichlet)
      scatradis->SetCondition((*conditername).second, rcp(new DRT::Condition(*conds[i])));
    }
    conds.clear();
  }

  // redistribute nodes to column (ghost) map
  DRT::UTILS::RedistributeWithNewNodalDistribution(*scatradis, *scatranoderowmap, *scatranodecolmap);
  scatradis->FillComplete();
}


#endif  // CCADISCRET
