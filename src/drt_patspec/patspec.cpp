/*!----------------------------------------------------------------------
\file patspec.cpp

\brief A collection of methods to modify patient specific geometries

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "patspec.H"
#include "../drt_mat/material.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/global_inp_control2.H"
#include <iostream>


using namespace std;


/*----------------------------------------------------------------------*
 |                                                             gee 03/10|
 *----------------------------------------------------------------------*/
void PATSPEC::PatientSpecificGeometry(DRT::Discretization& dis)
{
  if (!dis.Comm().MyPID())
  {
    cout << "____________________________________________________________\n";
    cout << "Entering patient specific structural preprocessing (PATSPEC)\n";
    cout << "\n";
  }




  //------------- test discretization for presence of the Gasser ILT material
  int lfoundit = 0;
  for (int i=0; i<dis.ElementRowMap()->NumMyElements(); ++i)
  {
    INPAR::MAT::MaterialType type = dis.lRowElement(i)->Material()->MaterialType();
    if (type == INPAR::MAT::m_aaagasser)
    {
      lfoundit = 1;
      break;
    }


    if (type == INPAR::MAT::m_aaa_mixedeffects)
    {
      lfoundit = 1;
      break;
    }
  }
  int gfoundit = 0;
  dis.Comm().SumAll(&lfoundit,&gfoundit,1);
  if (gfoundit)
  {
    if (!dis.Comm().MyPID())
      cout << "Computing distance functions...\n";
    PATSPEC::ComputeEleNormalizedLumenDistance(dis);
    PATSPEC::ComputeEleLocalRadius(dis);

  }
  //-------------------------------------------------------------------------
  
  
  
  
  if (!dis.Comm().MyPID())
  {
    cout << "\n";
    cout << "Leaving patient specific structural preprocessing (PATSPEC)\n";
    cout << "____________________________________________________________\n";
  }

  return;
}



/*----------------------------------------------------------------------*
 |                                                             gee 03/10|
 *----------------------------------------------------------------------*/
void PATSPEC::ComputeEleNormalizedLumenDistance(DRT::Discretization& dis)
{
  // find out whether we have a orthopressure or FSI condition
  vector<DRT::Condition*> conds;
  vector<DRT::Condition*> ortho(0);
  vector<DRT::Condition*> fsi(0);
  dis.GetCondition("SurfaceNeumann",ortho);
  for (int i=0; i<(int)ortho.size(); ++i)
  {
    if (ortho[i]->GType() != DRT::Condition::Surface) continue;
    const string* type = ortho[i]->Get<string>("type");
    if (*type == "neum_orthopressure")
      conds.push_back(ortho[i]);
  }

  dis.GetCondition("FSICoupling",fsi);
  for (int i=0; i<(int)fsi.size(); ++i)
  {
    if (fsi[i]->GType() != DRT::Condition::Surface) continue;
    conds.push_back(fsi[i]);
  }
  
  if (!conds.size())
    dserror("There is no orthopressure nor FSI condition in this discretization");
  
  // measure time as there is a brute force search in here
  Epetra_Time timer(dis.Comm());
  timer.ResetStartTime();
  
  // start building the distance function
  set<int> allnodes;
  for (int i=0; i<(int)conds.size(); ++i)
  {
    const vector<int>* nodes = conds[i]->Nodes();
    if (!nodes) dserror("Cannot find node ids in condition");
    for (int j=0; j<(int)nodes->size(); ++j)
      allnodes.insert((*nodes)[j]);
  }
  
/*
  for (int i=0; i<dis.Comm().NumProc(); ++i)
  {
    if (i == dis.Comm().MyPID())
    {
      set<int>::iterator fool;
      for (fool=allnodes.begin(); fool != allnodes.end(); ++fool)
        cout << "Proc " << i << " Node " << *fool << endl;
      cout << "Proc " << i << " sizeof allnodes " << (int)allnodes.size() << endl;
      fflush(stdout);
    }
    dis.Comm().Barrier();
  }
*/

  // create coordinates for all these nodes
  const int nnodes = (int)allnodes.size();
  vector<double> lcoords(nnodes*3,0.0);
  vector<double> gcoords(nnodes*3,0.0);
  set<int>::iterator fool;
  int count=0;
  for (fool=allnodes.begin(); fool != allnodes.end(); ++fool)
  {
    if (!dis.NodeRowMap()->MyGID(*fool)) 
    {
      count++;
      continue;
    }
    DRT::Node* node = dis.gNode(*fool);
    lcoords[count*3]   = node->X()[0];
    lcoords[count*3+1] = node->X()[1];
    lcoords[count*3+2] = node->X()[2];
    count++;
  }
  dis.Comm().SumAll(&lcoords[0],&gcoords[0],nnodes*3);
  lcoords.clear();
  allnodes.clear();

/*  
  for (int i=0; i<dis.Comm().NumProc(); ++i)
  {
    if (i == dis.Comm().MyPID())
    {
      for (int j=0; j<gcoords.size(); ++j)
        printf("Proc %d Coord %15.10e\n",i,gcoords[j]);
      fflush(stdout);
    }
    dis.Comm().Barrier();
  }
*/

  // compute distance of all of my nodes to these nodes
  // create vector for nodal values of ilt thickness
  // WARNING: This is a brute force expensive minimum distance search!
  const Epetra_Map* nrowmap = dis.NodeRowMap();
  RCP<Epetra_Vector> iltthick = LINALG::CreateVector(*nrowmap,true);
  for (int i=0; i<nrowmap->NumMyElements(); ++i)
  {
    const double* x = dis.gNode(nrowmap->GID(i))->X();
    // loop nodes from the condition and find minimum distance
    double mindist = 1.0e12;
    for (int j=0; j<nnodes; ++j)
    {
      double* xorth = &gcoords[j*3];
      double dist = sqrt( (x[0]-xorth[0])*(x[0]-xorth[0]) + 
                          (x[1]-xorth[1])*(x[1]-xorth[1]) + 
                          (x[2]-xorth[2])*(x[2]-xorth[2]) );
      if (dist<mindist) mindist = dist;
    }
    (*iltthick)[i] = mindist;
  }
  gcoords.clear();
  
  double maxiltthick;
  iltthick->MaxValue(&maxiltthick);
  maxiltthick -= 1.0; // subtract an approximate arterial wall thickness
  iltthick->Scale(1.0/maxiltthick);
  if (!dis.Comm().MyPID())
    printf("Max ILT thickness %10.5e\n",maxiltthick);
  
  // export nodal distances to column map
  RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(dis.NodeColMap()),true);
  LINALG::Export(*iltthick,*tmp);
  iltthick = tmp;
  
  // compute element-wise mean distance from nodal distance
  RCP<Epetra_Vector> iltele = LINALG::CreateVector(*(dis.ElementRowMap()),true);
  for (int i=0; i<dis.ElementRowMap()->NumMyElements(); ++i)
  {
    DRT::Element* actele = dis.gElement(dis.ElementRowMap()->GID(i));
    double mean = 0.0;
    for (int j=0; j<actele->NumNode(); ++j)
      mean += (*iltthick)[iltthick->Map().LID(actele->Nodes()[j]->Id())];
    mean /= actele->NumNode();
    (*iltele)[i] = mean;
  }
  tmp = LINALG::CreateVector(*(dis.ElementColMap()),true);
  LINALG::Export(*iltele,*tmp);
  iltele = tmp;

  // put this column vector of element ilt thickness in a condition and store in
  // discretization
  Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp(
    new DRT::Condition(0,DRT::Condition::PatientSpecificData,false,DRT::Condition::Volume));
  cond->Add("normalized ilt thickness",*iltele);

  //const Epetra_Vector* bla = cond->Get<Epetra_Vector>("ilt thickness");
  //cout << *bla;

  // check whether discretization has been filled before putting ocndition to dis
  bool filled = dis.Filled();
  dis.SetCondition("PatientSpecificData",cond);
  if (filled && !dis.Filled()) dis.FillComplete();
    
  
  if (!dis.Comm().MyPID())
    printf("Normalized ILT thickness computed in %10.5e sec\n",timer.ElapsedTime());
  
  return;
}



/*----------------------------------------------------------------------*
 |                                                           maier 11/10|
 *----------------------------------------------------------------------*/
void PATSPEC::ComputeEleLocalRadius(DRT::Discretization& dis)
{
  const ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  string filename = pslist.get<string>("CENTERLINEFILE");
  //cout << filename << endl;
  if (filename=="name.txt") 
  {
    cout << "No centerline file provided" << endl;
    // set element-wise mean distance to zero
    RCP<Epetra_Vector> locradele = LINALG::CreateVector(*(dis.ElementRowMap()),true);
     RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(dis.ElementColMap()),true);
    LINALG::Export(*locradele,*tmp);
    locradele = tmp;

    // put this column vector of element local radius in a condition and store in
    // discretization
    Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp( new DRT::Condition(0,DRT::Condition::PatientSpecificData,false,DRT::Condition::Volume));
    cond->Add("local radius",*locradele);

    // check whether discretization has been filled before putting ocndition to dis
    bool filled = dis.Filled();
    dis.SetCondition("PatientSpecificData",cond);
    if (filled && !dis.Filled()) dis.FillComplete();

    if (!dis.Comm().MyPID())
      printf("No local radii computed\n");
  
    return;
  }

  ifstream file (filename.c_str());
  if (file == NULL) dserror ("Error opening centerline file");
  string sLine;
  char buffer[1000];
  vector<double> clcoords;
  while (getline(file, sLine))
  {
    if (sLine.size() != 0) 
    {
      strcpy(buffer,sLine.c_str());
      char * pEnd;
    //only get values if first occurence of line is double number
      if (strtod(buffer, &pEnd) != 0.0)
      {
	//store x,y,z corrdinates
	clcoords.push_back(strtod(buffer, &pEnd));
	clcoords.push_back(strtod(pEnd, &pEnd));
	clcoords.push_back(strtod(pEnd, NULL));
      }
    }
  }
  //  for (int i=0; i<(clcoords.size()/3); ++i) cout << clcoords[i*3] << " " << clcoords[i*3+1] << " " << clcoords[i*3+2] << " " << endl;

  // measure time as there is a brute force search in here
  Epetra_Time timer(dis.Comm());
  timer.ResetStartTime();

  // compute distance of all of my nodes to the clpoints
  // create vector for nodal values of ilt thickness
  // WARNING: This is a brute force expensive minimum distance search!
  const Epetra_Map* nrowmap = dis.NodeRowMap();
  RCP<Epetra_Vector> localrad = LINALG::CreateVector(*nrowmap,true);
  for (int i=0; i<nrowmap->NumMyElements(); ++i)
  {
    const double* x = dis.gNode(nrowmap->GID(i))->X();
    // loop nodes from the condition and find minimum distance
    double mindist = 1.0e12;
    for (int j=0; j<(int)(clcoords.size()/3); ++j)
    {
      double* xcl = &clcoords[j*3];
      double dist = sqrt( (x[0]-xcl[0])*(x[0]-xcl[0]) + 
                          (x[1]-xcl[1])*(x[1]-xcl[1]) + 
                          (x[2]-xcl[2])*(x[2]-xcl[2]) );
      if (dist<mindist) mindist = dist;
    }
    (*localrad)[i] = mindist;
  }
  clcoords.clear();
  // max local rad just for information purpose:
  double maxlocalrad;
  localrad->MaxValue(&maxlocalrad);
  if (!dis.Comm().MyPID())
    printf("Max local radius %10.5e\n",maxlocalrad);
  
  // export nodal distances to column map
  RCP<Epetra_Vector> tmp = LINALG::CreateVector(*(dis.NodeColMap()),true);
  LINALG::Export(*localrad,*tmp);
  localrad = tmp;
  
  // compute element-wise mean distance from nodal distance
  RCP<Epetra_Vector> locradele = LINALG::CreateVector(*(dis.ElementRowMap()),true);
  for (int i=0; i<dis.ElementRowMap()->NumMyElements(); ++i)
  {
    DRT::Element* actele = dis.gElement(dis.ElementRowMap()->GID(i));
    double mean = 0.0;
    for (int j=0; j<actele->NumNode(); ++j)
      mean += (*localrad)[localrad->Map().LID(actele->Nodes()[j]->Id())];
    mean /= actele->NumNode();
    (*locradele)[i] = mean;
  }

  tmp = LINALG::CreateVector(*(dis.ElementColMap()),true);
  LINALG::Export(*locradele,*tmp);
  locradele = tmp;

  // put this column vector of element local radius in a condition and store in
  // discretization
  Teuchos::RCP<DRT::Condition> cond = Teuchos::rcp(
    new DRT::Condition(0,DRT::Condition::PatientSpecificData,false,DRT::Condition::Volume));
  cond->Add("local radius",*locradele);

  // check whether discretization has been filled before putting ocndition to dis
  bool filled = dis.Filled();
  dis.SetCondition("PatientSpecificData",cond);
  if (filled && !dis.Filled()) dis.FillComplete();

  if (!dis.Comm().MyPID())
    printf("Local radii computed in %10.5e sec\n",timer.ElapsedTime());
  
  return;
}


/*----------------------------------------------------------------------*
 |                                                             gee 03/10|
 *----------------------------------------------------------------------*/
void PATSPEC::GetILTDistance(const int eleid, 
                             Teuchos::ParameterList& params, 
                             DRT::Discretization& dis)
{
  vector<DRT::Condition*> mypatspeccond;
  dis.GetCondition("PatientSpecificData", mypatspeccond);
  if (!mypatspeccond.size()) return;

  const Epetra_Vector* fool = mypatspeccond[0]->Get<Epetra_Vector>("normalized ilt thickness");
  if (!fool) return;
  
  if (!fool->Map().MyGID(eleid)) dserror("I do not have this element");
  
  double meaniltthick = (*fool)[fool->Map().LID(eleid)];
  params.set("iltthick meanvalue",meaniltthick);

  return;
}


/*----------------------------------------------------------------------*
 |                                                           maier 11/10|
 *----------------------------------------------------------------------*/
void PATSPEC::GetLocalRadius(const int eleid, 
                             Teuchos::ParameterList& params, 
                             DRT::Discretization& dis)
{
  vector<DRT::Condition*> mypatspeccond;
  dis.GetCondition("PatientSpecificData", mypatspeccond);
  if (!mypatspeccond.size()) return;

  const Epetra_Vector* fool = mypatspeccond[1]->Get<Epetra_Vector>("local radius");
  if (!fool) return;
  
  if (!fool->Map().MyGID(eleid)) dserror("I do not have this element");
  
  double meanlocalrad = (*fool)[fool->Map().LID(eleid)];
  params.set("localrad meanvalue",meanlocalrad);

  return;
}




#endif
