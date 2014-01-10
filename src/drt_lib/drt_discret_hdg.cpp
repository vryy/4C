/*!----------------------------------------------------------------------
\file drt_discret_hdg.cpp

\brief a class to manage an enhanced discretization including all faces for HDG


</pre>

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "drt_discret_hdg.H"

#include "drt_exporter.H"
#include "drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

DRT::DiscretizationHDG::DiscretizationHDG(const std::string name,
                                          Teuchos::RCP<Epetra_Comm> comm)
  :
  DiscretizationFaces(name, comm)
{
  this->doboundaryfaces_ = true;
}

/*----------------------------------------------------------------------*
 |  Finalize construction (public)                     kronbichler 12/13|
 *----------------------------------------------------------------------*/
int DRT::DiscretizationHDG::FillComplete(bool assigndegreesoffreedom,
                                         bool initelements,
                                         bool doboundaryconditions)

{
  // call FillComleteFaces of base class with create_faces set to true
  this->FillCompleteFaces(Teuchos::null, assigndegreesoffreedom, initelements,
                          doboundaryconditions, true);

  // get the correct face orientation from the owner. since the elements in general do not allow
  // packing, extract the node ids, communicate them, and change the node ids in the element
  Exporter nodeexporter( *facerowmap_, *facecolmap_, Comm() );
  std::map<int, std::vector<int> > nodeIds;
  for (typename std::map<int, Teuchos::RCP<DRT::Element> >::const_iterator f=faces_.begin();
       f!= faces_.end(); ++f) {
    std::vector<int> ids(f->second->NumNode());
    for (int i=0; i<f->second->NumNode(); ++i)
      ids[i] = f->second->NodeIds()[i];
    nodeIds[f->first] = ids;
  }

  nodeexporter.Export( nodeIds );

  for (typename std::map<int, Teuchos::RCP<DRT::Element> >::iterator f=faces_.begin();
       f!= faces_.end(); ++f) {
    if ( f->second->Owner() == Comm().MyPID() )
      continue;
    std::vector<int> &ids = nodeIds[f->first];
    dsassert(ids.size() > 0, "Lost a face during communication");
    f->second->SetNodeIds(ids.size(), &ids[0]);

    // refresh node pointers if they have been set up
    DRT::Node** oldnodes = f->second->Nodes();
    if (oldnodes != 0) {
      std::vector<DRT::Node*> nodes(ids.size(), 0);

      for (unsigned int i=0; i<ids.size(); ++i) {
        for (unsigned int j=0; j<ids.size(); ++j)
          if (oldnodes[j]->Id() == ids[i]) {
            nodes[i] = oldnodes[j];
          }
        dsassert(nodes[i] != 0, "Could not find node.");
      }
      f->second->BuildNodalPointers(&nodes[0]);
    }
  }

  return 0;
}



void DRT::DiscretizationHDG::DoDirichletCondition(DRT::Condition&             cond,
                                                  const bool                  usetime,
                                                  const double                time,
                                                  Teuchos::RCP<Epetra_Vector> systemvector,
                                                  Teuchos::RCP<Epetra_Vector> systemvectord,
                                                  Teuchos::RCP<Epetra_Vector> systemvectordd,
                                                  Teuchos::RCP<Epetra_Vector> toggle,
                                                  Teuchos::RCP<std::set<int> > dbcgids)
{
  Discretization::DoDirichletCondition(cond, usetime,time,systemvector,systemvectord,systemvectordd,toggle,dbcgids);
  if (FaceRowMap() == NULL)
    return;

  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const std::vector<int>*    curve  = cond.Get<std::vector<int> >("curve");
  const std::vector<int>*    funct  = cond.Get<std::vector<int> >("funct");
  const std::vector<int>*    onoff  = cond.Get<std::vector<int> >("onoff");
  const std::vector<double>* val    = cond.Get<std::vector<double> >("val");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvector != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvector;
  }
  if (systemvectord != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectord;
  }
  if (systemvectordd != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null)
      systemvectoraux = systemvectordd;
  }
  dsassert(systemvectoraux!=Teuchos::null, "At least one vector must be unequal to null");

  // factor given by time curve
  std::vector<std::vector<double> > curvefacs(onoff->size());
  for (unsigned int j=0; j<onoff->size(); ++j)
  {
    int curvenum = -1;
    if (curve) curvenum = (*curve)[j];
    if (curvenum>=0 && usetime)
      curvefacs[j] = DRT::Problem::Instance()->Curve(curvenum).FctDer(time,deg);
    else
    {
      curvefacs[j].resize(deg+1, 1.0);
      for (unsigned i=1; i<(deg+1); ++i) curvefacs[j][i] = 0.0;
    }
  }

  if (NumMyRowFaces() > 0)
  {
    Epetra_SerialDenseVector elevec1, elevec2, elevec3;
    Epetra_SerialDenseMatrix elemat1, elemat2;
    std::vector<int> dummy;
    Teuchos::ParameterList initParams;
    //initParams.set<int>("action", 500); // TODO: Introduce a general action type that is valid for all problems
    if (funct != NULL) {
      Teuchos::Array<int> functarray(*funct);
      initParams.set("funct",functarray);
    }
    Teuchos::Array<int> onoffarray(*onoff);
    initParams.set("onoff",onoffarray);
    initParams.set("time",time);

    bool pressureDone = this->Comm().MyPID() != 0;

    for (int i=0; i<NumMyRowFaces(); ++i)
    {
      // do not consider internal faces
      if (lRowFace(i)->ParentSlaveElement() != NULL) continue;

      const unsigned int dimension = DRT::UTILS::getDimension(lRowFace(i)->ParentMasterElement()->Shape());
      if (onoff->size() <= dimension || (*onoff)[dimension] == 0)
        pressureDone = true;
      if (!pressureDone) {
        if (this->NumMyRowElements() > 0 && this->Comm().MyPID()==0) {
          std::vector<int> predof = this->Dof(0, lRowElement(0));
          const int gid = predof[0];
          const int lid = this->DofRowMap(0)->LID(gid);
          // amend vector of DOF-IDs which are Dirichlet BCs
          if (systemvector != Teuchos::null)
            (*systemvector)[lid] = 0;
          if (systemvectord != Teuchos::null)
            (*systemvectord)[lid] = 0;
          if (systemvectordd != Teuchos::null)
            (*systemvectordd)[lid] = 0;
          // set toggle vector
          if (toggle != Teuchos::null)
            (*toggle)[lid] = 1.0;
          // amend vector of DOF-IDs which are Dirichlet BCs
          if (dbcgids != Teuchos::null)
            (*dbcgids).insert(gid);
          pressureDone = true;
        }
      }

      int nummynodes = lRowFace(i)->NumNode();
      const int * mynodes = lRowFace(i)->NodeIds();

      // do only faces where all nodes are present in the node list
      bool faceRelevant = true;
      for (int j=0; j<nummynodes; ++j)
        if (!cond.ContainsNode(mynodes[j]))
        {
          faceRelevant = false;
          break;
        }
      if (!faceRelevant) continue;

      const unsigned int dofperface = lRowFace(i)->ParentMasterElement()->NumDofPerFace(lRowFace(i)->FaceMasterNumber());
      initParams.set<unsigned int>("faceconsider",
          static_cast<unsigned int>(lRowFace(i)->FaceMasterNumber()));
      if (static_cast<unsigned int>(elevec1.M()) != dofperface)
        elevec1.Shape(dofperface, 1);
      std::vector<int> dofs = this->Dof(0,lRowFace(i));

      bool do_evaluate = false;
      if (funct != NULL)
        for (unsigned int i=0; i<dimension; ++i)
          if ((*funct)[i] > 0)
            do_evaluate = true;

      if (do_evaluate)
        lRowFace(i)->ParentMasterElement()->Evaluate(initParams,*this,dummy,elemat1,elemat2,elevec1,elevec2,elevec3);
      else
        for (unsigned int i=0; i<dofperface; ++i)
          elevec1(i) = 1.;

      const unsigned int dofpercomponent = dofperface/dimension;
      for (unsigned int i=0; i<dofperface; ++i) {

        // TODO make general, not only for fluid with dimension components
        int onesetj = i / dofpercomponent;

        if ((*onoff)[onesetj]==0)
        {
          const int lid = (*systemvectoraux).Map().LID(dofs[i]);
          if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[i]);
          if (toggle!=Teuchos::null)
            (*toggle)[lid] = 0.0;
          // get rid of entry in DBC map - if it exists
          if (dbcgids != Teuchos::null)
            (*dbcgids).erase(dofs[i]);
          continue;
        }

        // get global id
        const int gid = dofs[i];

        // assign value
        const int lid = (*systemvectoraux).Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        if (systemvector != Teuchos::null)
          (*systemvector)[lid] = (*val)[onesetj] * elevec1(i) * curvefacs[onesetj][0];
        if (systemvectord != Teuchos::null)
          (*systemvectord)[lid] = (*val)[onesetj] * elevec1(i) * curvefacs[onesetj][1];
        if (systemvectordd != Teuchos::null)
          (*systemvectordd)[lid] = (*val)[onesetj] * elevec1(i) * curvefacs[onesetj][2];
        // set toggle vector
        if (toggle != Teuchos::null)
          (*toggle)[lid] = 1.0;
        // amend vector of DOF-IDs which are Dirichlet BCs
        if (dbcgids != Teuchos::null)
          (*dbcgids).insert(gid);
      }
    }  // loop over faces
  }
}



/*----------------------------------------------------------------------*
 |  << operator                                        kronbichler 12/13|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const DRT::DiscretizationHDG& dis)
{
  // print standard discretization info
  dis.Print(os);
  // print additional info about internal faces
  dis.PrintFaces(os);

  return os;
}

