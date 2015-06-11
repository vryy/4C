/*----------------------------------------------------------------------*/
/*!
\file pat_matpar_manager.cpp
\brief manage material parameters during optimization

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            089 - 289-15271
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include "pat_matpar_manager.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*/
ACOU::PatMatParManagerUniform::PatMatParManagerUniform(Teuchos::RCP<DRT::Discretization> discret)
:MatParManagerUniform(discret)
{
// TODO think about different inheritance
}

/*----------------------------------------------------------------------*/
ACOU::PatMatParManagerPerElement::PatMatParManagerPerElement(Teuchos::RCP<DRT::Discretization> discret, bool scaleele)
:MatParManagerPerElement(discret)
{
  scalegradele_ = scaleele;
}

/*----------------------------------------------------------------------*/
void ACOU::PatMatParManagerUniform::AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  Teuchos::RCP<const Epetra_Vector> phi = Discret()->GetState("phi");

  for (int i=0; i<Discret()->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = Discret()->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (ParaMap().find(elematid) == ParaMap().end() )
      continue;

    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_integr_grad_reac);

    // this works if we optimize only with respect to reac
    p.set<bool>("signum_mu",actele->Material()->Parameter()->GetParameter(1,Discret()->ElementColMap()->LID(actele->Id()))<0.0);
    p.set<bool>("scaleele",false);

    std::vector<int> actparams = ParaMap().at( elematid );
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      //initialize element vectors
      int ndof = actele->NumNode();
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);

      DRT::Element::LocationArray la(Discret()->NumDofSets());
      actele->LocationVector(*Discret(),la,false);
      actele->Evaluate(p,*Discret(),la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

      //reuse elevector2
      for (int l=0; l<(int)la[0].lm_.size(); l++)
      {
        int lid=phi->Map().LID(la[0].lm_.at(l));
        if (lid==-1) dserror("not found on this processor");
        elevector2[l] = (*phi)[lid];
      }
      double val2 = elevector2.Dot(elevector1);

      // Assemble the final gradient; this is parametrization class business
      // (i.e contraction to (optimization)-parameter space:
      ContractGradient(dfint,val2,actele->Id(),ParaPos().at(elematid).at(it-actparams.begin()), it-actparams.begin());

    }//loop this elements material parameters (only the ones to be optimized)

  }//loop elements

  return;
}
/*----------------------------------------------------------------------*/
void ACOU::PatMatParManagerPerElement::AddEvaluate(double time, Teuchos::RCP<Epetra_MultiVector> dfint)
{
  Teuchos::RCP<const Epetra_Vector> phi = Discret()->GetState("phi");

  for (int i=0; i<Discret()->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = Discret()->lRowElement(i);
    int elematid = actele->Material()->Parameter()->Id();

    if (ParaMap().find(elematid) == ParaMap().end() )
    {
      std::cout<<"Warning, skipping elematid "<<elematid<<" in ele "<<actele->Id()<<std::endl;
      continue;
    }
    // list to define routines at elementlevel
    Teuchos::ParameterList p;
    p.set<int>("action",SCATRA::calc_integr_grad_reac);

    // this works if we optimize only with respect to reac
    p.set<bool>("signum_mu",actele->Material()->Parameter()->GetParameter(1,Discret()->ElementColMap()->LID(actele->Id()))<0.0);
    p.set<bool>("scaleele",scalegradele_);

    std::vector<int> actparams = ParaMap().at( elematid );
    std::vector<int>::const_iterator it;
    for ( it=actparams.begin(); it!=actparams.end(); it++)
    {
      //initialize element vectors
      int ndof = actele->NumNode();
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);

      DRT::Element::LocationArray la(Discret()->NumDofSets());
      actele->LocationVector(*Discret(),la,false);
      actele->Evaluate(p,*Discret(),la,elematrix1,elematrix2,elevector1,elevector2,elevector3);

      //reuse elevector2
      for (int l=0; l<(int)la[0].lm_.size(); l++)
      {
        int lid=phi->Map().LID(la[0].lm_.at(l));
        if (lid==-1) dserror("not found on this processor");
        elevector2[l] = (*phi)[lid];
      }
      double val2 = elevector2.Dot(elevector1);

      // Assemble the final gradient; this is parametrization class business
      // (i.e contraction to (optimization)-parameter space:
      ContractGradient(dfint,val2,actele->Id(),ParaPos().at(elematid).at(it-actparams.begin()), it-actparams.begin());

    }//loop this elements material parameters (only the ones to be optimized)

  }//loop elements

  return;
}

/*----------------------------------------------------------------------*/
/* build blockwise connectivity graphs                   schoeder 03/15 */
/*----------------------------------------------------------------------*/
void ACOU::PatMatParManagerPerElement::FillAdjacencyMatrix(const Epetra_Map& paramrowmap, Teuchos::RCP<Epetra_CrsMatrix> graph)
{
  // this is basically a copy of sebastian's method in the father class, only the call for the element,
  // e.g. the action is different

  /*------------------------------------------------------------------- */
  // STEP 1: loop elements in elerowmap and store map of faces(vector
  // of NodeIds in sorted order) with corresponding gids of the elements

  std::map< std::vector<int>, std::vector<int> > facemap; //map faces to corresponding elements/parameters
  std::map< std::vector<int>, double > faceweight;        //map of faces to its area
  for (int i=0; i<paramrowmap.NumMyElements(); i++)
  {
    // the current element
    int pgid = paramrowmap.GID(i); // !! the local id of the partial map is not the local parameter id!!
    int plid = paramlayoutmap_->LID(pgid);
    int elegid = ParamsLIDtoeleGID()[plid];
    DRT::Element* ele=Discret()->gElement(elegid);
    if (ele == NULL) dserror("element not found here");

    // decide whether we deal with a 2D or 3D discretization
    unsigned int nele=0;
    const DRT::Element::DiscretizationType distype = ele->Shape();
    std::vector< std::vector<int> > faces;
    if (ele->NumSurface() > 1)   // 2D boundary element and 3D parent element
    {
      nele = ele->NumSurface();
      faces = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
    }
    else if (ele->NumSurface() == 1) // 1D boundary element and 2D parent element
    {
      nele = ele->NumLine();
      faces = DRT::UTILS::getEleNodeNumberingLines(distype);
    }
    else dserror("creating internal faces for 1D elements (would be points) not implemented yet");

    if (nele != faces.size()) dserror("number of surfaces or lines does not match!");

    // get the surface/line elements for area computation
    std::vector<Teuchos::RCP<DRT::Element> >  surfs;
    if (ele->NumSurface() > 1)
      surfs = ele->Surfaces();
    else if (ele->NumSurface() == 1)
      surfs = ele->Lines();
    else dserror("creating internal faces for 1D elements (would be points) not implemented yet");

    // get nodes of each of this element's face
    for (unsigned int iele = 0; iele < nele; iele++)
    {
      // allocate node vectors
      unsigned int nnode = faces[iele].size();
      std::vector<int> nodeids(nnode);

      // get connectivity info
      for (unsigned int inode=0;inode<nnode;inode++)
      {
        nodeids[inode] = ele->NodeIds()[faces[iele][inode]];
      }

      //get the area of this face
      Teuchos::ParameterList p;
      p.set<int>("action", SCATRA::bd_integrate_shape_functions);
      p.set<double>("boundaryint",0.0);

      DRT::Element::LocationArray la(Discret()->NumDofSets());
      surfs[iele]->LocationVector(*Discret(),la,false);
      //initialize element vectors
      int ndof = ele->NumNode()*3;
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);
      surfs[iele]->Evaluate(p,*Discret(),la,elematrix1,elematrix2,elevector1,elevector2,elevector3);
      double area = p.get("boundaryint",-1.0);
      if (area < 0.0) dserror("area computation of surface failed");

      // sort the nodes in faces to make them to be used for keys in the facemap.
      // they are not unique (in a parallel layout sense) though since a face can
      // be on multiple processors
      std::sort( nodeids.begin(), nodeids.end() );

      // store parameter global id corresponding to this face
      facemap[nodeids].push_back(pgid);
      faceweight[nodeids]=area;
    }
  } //loop local elements


  /*------------------------------------------------------------------- */
  // STEP 2: gather for each face on each proc the same set of corresponding
  // elements; ntargetprocs is equal to the total number of processors to make
  // data redundant on all procs
  // the face weight dont need to be communicated since they should be the
  // same on every proc having a face
  const int numprocs = Discret()->Comm().NumProc();
  int allproc[numprocs];
  for (int i = 0; i < numprocs; ++i)
    allproc[i] = i;

  for( int i=0; i<Discret()->Comm().NumProc(); i++ )
  {
    // by looping the procs one by one we need to make sure that every proc
    // does the same number of loops regardless of its own lenght of keys
    // in the facemap. so the actual number of loops is defined by the current
    // proc and distributed in localmapsize
    int localmapsize;
    if( Discret()->Comm().MyPID()==i ) localmapsize = facemap.size();
    Discret()->Comm().Broadcast(&localmapsize,1,i);

    // now iterate as often as there are faces on proc i
    std::map< std::vector<int>,std::vector<int> >::iterator face_it(facemap.begin());
    for( int j=0; j<localmapsize; j++ )
    {
      //get length of current face-key on all procs
      int keylength=0;
      if( Discret()->Comm().MyPID()==i ) keylength=face_it->first.size();
      Discret()->Comm().Broadcast(&(keylength),1,i);

      //send current face-key to all procs
      std::vector<int> facekey(keylength,0);
      if( Discret()->Comm().MyPID()==i ) facekey=face_it->first;
      Discret()->Comm().Broadcast(&(facekey[0]),keylength,i);

      //check whether one of the other procs also has this key and write IDs of
      // procs who own this face in "sowningprocs" and distribute this knowledge
      // among all procs to rowningprocs
      std::map< std::vector<int>,std::vector<int> >::iterator face_abroad(facemap.find(facekey));
      std::vector<int> sowningprocs;
      std::vector<int> rowningprocs;
      if( face_abroad!=facemap.end() && Discret()->Comm().MyPID() != i)
        sowningprocs.push_back(Discret()->Comm().MyPID());
      LINALG::Gather(sowningprocs,rowningprocs,numprocs,allproc,Discret()->Comm());

      // now bring parameters corresponding to this face on the other procs to proc i
      // (they are send to all procs but only proc i stores them in the map with
      // the current face-key)
      std::vector<int> sparams;
      std::vector<int> rparams;
      if( std::find(rowningprocs.begin(),rowningprocs.end(),Discret()->Comm().MyPID()) != rowningprocs.end())
        sparams=facemap[facekey];
      LINALG::Gather(sparams,rparams,numprocs,allproc,Discret()->Comm());

      // store additional elements on proc i
      if( Discret()->Comm().MyPID() == i )
      {
        for (int iele=0; iele<(int)rparams.size(); iele++)
        {
          // depending on which proc comes first it might be that parameters
          // are added redundantly to a face which leads to undesired summation
          // upon inserting weights into the matrix. So only add if not existent yet
          if ( std::find(facemap[facekey].begin(),facemap[facekey].end(),rparams[iele]) == facemap[facekey].end() )
            facemap[facekey].push_back(rparams[iele]);
        }
      }
      // increase face pointer on this proc
      if( Discret()->Comm().MyPID()==i ) face_it++;

    } // faces on each proc
  } // procs

//  //Debug print out
//  std::map<std::vector<int>, std::vector<int> >::iterator tmp;
//  for( tmp=facemap.begin(); tmp!=facemap.end(); tmp++)
//  {
//    std::cout << "Proc : " << Discret()->Comm().MyPID() << " FACE: ";
//    for (int i=0; i<(int)tmp->first.size(); i++)
//      std::cout << tmp->first[i] << " ";
//
//    std::cout << "PARAMS: ";
//    for (int i=0; i<(int)tmp->second.size(); i++)
//      std::cout << tmp->second[i] << " ";
//
//    std::cout << " " << std::endl;
//  }
//
//  Discret()->Comm().Barrier();
//  if( Discret()->Comm().MyPID()==0 )std::cout << " " << std::endl;
//  if( Discret()->Comm().MyPID()==0 )std::cout << " " << std::endl;

  /*------------------------------------------------------------------- */
  // STEP 3: sort elements into graph according to neighbour information
  // in facemap
  std::map<std::vector<int>, std::vector<int> >::iterator faces;
  for( faces=facemap.begin(); faces!=facemap.end(); faces++)
  {
    // all parameter ids connected to this face
    std::vector<int> parameters = faces->second;

    // weight corresponding to these parameters
    std::vector<double> weights(parameters.size(),faceweight[faces->first]);

    for(int iele=0; iele<(int)parameters.size(); iele++)
    {
      int globalrow=parameters[iele];
      if (paramrowmap.MyGID(globalrow))
      {
        // like this the diagonal entries are inserted redundantly and summed up
        // after FillComplete() is called; they are more or less useless anyways
        graph->InsertGlobalValues(globalrow,parameters.size(),&weights[0],&parameters[0]);
      }
    }
  }

  return;
}
