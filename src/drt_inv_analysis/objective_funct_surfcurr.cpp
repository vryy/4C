/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct_surfcurr.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/
/*----------------------------------------------------------------------*/


#include "objective_funct_surfcurr.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "Epetra_MultiVector.h"
#include "../linalg/linalg_utils.H"

#include "matpar_manager.H"


/*----------------------------------------------------------------------*/
/* standard constructor of current representation                       */
/*----------------------------------------------------------------------*/
STR::INVANA::ObjectiveFunctSurfCurrRepresentation::ObjectiveFunctSurfCurrRepresentation(Teuchos::RCP<DRT::Discretization> discret,
                                                                                        int steps,
                                                                                        Teuchos::RCP<std::vector<double> > timesteps):
sourcedis_(discret),
targetdis_(Teuchos::null),
timesteps_(timesteps),
msteps_(steps)
{
  // set up source discretization
  if (not sourcedis_->Filled() || not sourcedis_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");

  // set up target discretization we don't need dofmaps for the target
  targetdis_ = DRT::Problem::Instance(1)->GetDis("structure");
  if (not targetdis_->Filled()) targetdis_->FillComplete();

  // get the conditions for the current evaluation
  std::vector<DRT::Condition* > scc_source;
  std::vector<DRT::Condition* > scc_target;
  sourcedis_->GetCondition("SurfaceCurrent",scc_source);
  targetdis_->GetCondition("SurfaceCurrent",scc_target);

  // check wether length of both conditions is equal:
  if (scc_source.size()!=scc_target.size())
    dserror("Size of conditions in source and target differs. We need the size to be equal for meaningful computation");

  // find corresponding conditions in source and target discretization
  // and initialize the associated objective functions
  for (int i=0; i<(int)scc_source.size(); i++)
  {
    bool foundit=false;
    int ids = scc_source[i]->GetInt("matching id");
    for (int j=0; j<(int)scc_target.size(); j++)
    {
      int idt = scc_target[j]->GetInt("matching id");
      if (idt==ids)
      {
        currents_.push_back(Teuchos::rcp(new ObjectiveFunctSurfCurr(sourcedis_,targetdis_,scc_source[i],scc_target[j])));
        foundit=true;
        break;
      }
    }

    if (foundit)
      continue;
    else
      dserror("corresponding condition in target not found");
  }

  if (currents_.size()!=scc_source.size())
    dserror("problem in finding corresponding current conditions in source and target discretization");

}

/*----------------------------------------------------------------------*/
/* Evaluate value of the objective function                  keh 11/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::ObjectiveFunctSurfCurrRepresentation::Evaluate(Teuchos::RCP<Epetra_MultiVector> disp,
                                                                 double& val)
{
  val = 0.0;

  //so far only the case were one single "measurement" for the final timestep of the simulation is considered
  int step=msteps_-1; // so this is the last one
  sourcedis_->SetState("displacements", Teuchos::rcp((*disp)(step), false));

  for(int i=0; i<(int)currents_.size(); ++i)
  {
    // evaluate every single surface combination
    val+=currents_[i]->WSpaceNorm();
  }
}

/*----------------------------------------------------------------------*/
/* Evaluate the gradient of the objective function                      */
/* w.r.t the displacements                                   keh 11/13  */
/*----------------------------------------------------------------------*/
void STR::INVANA::ObjectiveFunctSurfCurrRepresentation::EvaluateGradient(Teuchos::RCP<Epetra_MultiVector> disp,
                                                                         Teuchos::RCP<Epetra_MultiVector> gradient)
{
  //so far only the case were one single "measurement" for the final timestep of the simulation is considered
  int step=msteps_-1; // so this is the last one
  gradient->PutScalar(0.0);

  sourcedis_->SetState("displacements", Teuchos::rcp((*disp)(step), false));

  for(int i=0; i<(int)currents_.size(); ++i)
  {
    // evaluate every single surface combination
    currents_[i]->GradientWSpaceNorm(gradient,step);
  }
}

void STR::INVANA::ObjectiveFunctSurfCurrRepresentation::SetScale(double sigmaW)
{
  for(int i=0; i<(int)currents_.size(); ++i)
  {
    currents_[i]->SetScale(sigmaW);
  }
}


/*----------------------------------------------------------------------*/
/* standard constructor for a surface current                keh 11/13  */
/*----------------------------------------------------------------------*/
STR::INVANA::ObjectiveFunctSurfCurr::ObjectiveFunctSurfCurr(Teuchos::RCP<DRT::Discretization> sourcedis,
                                                            Teuchos::RCP<DRT::Discretization> targetdis,
                                                            DRT::Condition* sourcecond,
                                                            DRT::Condition* targetcond):
sourcedis_(sourcedis),
targetdis_(targetdis),
sourcecond_(sourcecond),
targetcond_(targetcond),
sourcemap_(Teuchos::null),
targetmap_(Teuchos::null),
sourcemapred_(Teuchos::null),
targetmapred_(Teuchos::null)
{
  // get the scale of the kernel
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  double sigmaW = statinvp.get<double>("KERNELSCALE");
  if (sigmaW<0.0) dserror("you need to choose a proper scale (KERNELSCALE) at which to evaluate the current");
  sigmaW2_=sigmaW*sigmaW;

  // generate map of boundary elements of source and target
  sourcemap_ = Teuchos::rcp(new Epetra_Map(SetupConditionMap(sourcecond_)));
  targetmap_ = Teuchos::rcp(new Epetra_Map(SetupConditionMap(targetcond_)));

  //cout << "source map: " << endl;
  //cout << *sourcemap_ << endl;

  // compute the all reduced maps
  sourcemapred_ = LINALG::AllreduceEMap(*sourcemap_,0);
  targetmapred_ = LINALG::AllreduceEMap(*targetmap_,0);

  //cout << "source ele col map " << *(sourcedis_->ElementColMap()) << endl;
  //cout << "target ele col map " << *(targetdis_->ElementColMap()) << endl;
  //cout << "full map" << *sourcemap_ << endl;
  //cout << "reduced map" << *sourcemapred_ << endl;

  //complete evaluation of the target data can be done here once
  ComputeNormalCenterMaterialConfig(targetcond_,targetdis_,&normal_t_,&center_t_);

  // all reduce target data to proc0 here already, since it is the same throughout
  DRT::Exporter ex(*targetmap_,*targetmapred_,targetdis_->Comm());
  ex.Export(normal_t_);
  ex.Export(center_t_);

  // precompute structural integrity component of the functional
  structinttarget_=Convolute(normal_t_,center_t_,normal_t_,center_t_);

  //cout << "target stuff" << endl;
  //PrintDataToScreen(center_t_);
  //PrintDataToScreen(normal_t_);

}



const Epetra_Map STR::INVANA::ObjectiveFunctSurfCurr::SetupConditionMap(DRT::Condition* cond)
{
  // get total number of elements in this condition
  int lnumele = cond->Geometry().size();
  int gnumele;
  cond->Comm()->SumAll(&lnumele,&gnumele,1);


  // a vector to hold the gids
  std::vector<int> gids;

  // get the geometry and loop the elements
  std::map<int,RCP<DRT::Element> >& geom = cond->Geometry();
  std::map<int,RCP<DRT::Element> >::iterator ele;
  for (ele=geom.begin(); ele != geom.end(); ++ele)
  {
    DRT::Element* element = ele->second.get();
    //cout << "pid: " << cond->Comm()->MyPID() << " eleid: " << element->Id() << endl;
    switch(element->NumNode())
    {
    case 3:
    {
      gids.push_back(element->Id());
    }
    break;
    case 4:
    {
      gids.push_back(element->Id());
      gids.push_back(element->Id()+gnumele);
    }
    break;
    default:
        dserror("shape type not supported in here!\n");
    break;
    }
  }

  // maybe use the LINALG map creator here?
  Epetra_Map elemap = Epetra_Map(-1,(int)gids.size(),gids.data(),0,*(cond->Comm()));
  return elemap;
}

void  STR::INVANA::ObjectiveFunctSurfCurr::ComputeNormalCenterMaterialConfig(DRT::Condition* cond,
                                                                             Teuchos::RCP<DRT::Discretization> discret,
                                                                             std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* normals,
                                                                             std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* centers)
{
  //clear data maps
  normals->clear();
  centers->clear();

  // get total number of elements in this condition
  int lnumele = cond->Geometry().size();
  int gnumele;
  cond->Comm()->SumAll(&lnumele,&gnumele,1);

  // get the geometry and loop the elements
  std::map<int,RCP<DRT::Element> >& geom = cond->Geometry();
  std::map<int,RCP<DRT::Element> >::iterator ele;
  for (ele=geom.begin(); ele != geom.end(); ++ele)
  {
    DRT::Element* element = ele->second.get();
    int numnodes = element->NumNode();
    int egid = element->Id();

    // fill maps of normals and centers with data:
    // case tri  -> no special treatment
    // case quad -> split into 2 tris and associate to the second tri a gid as gid=egid+maxele
    LINALG::SerialDenseMatrix x(3,3); //node coordinates
    Epetra_SerialDenseMatrix c(3,1); // facet center
    Epetra_SerialDenseMatrix n(3,1); // facet normal
    if (numnodes==3)
    {
      for (int i=0; i<3; ++i)
      {
        x(i,0) = element->Nodes()[i]->X()[0];
        x(i,1) = element->Nodes()[i]->X()[1];
        x(i,2) = element->Nodes()[i]->X()[2];
      }
      ComputeNormalCenter(x,n,c);

      // plug everything into the data maps
      normals->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(n)) ));
      centers->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(c)) ));
    }
    else if (numnodes==4)
    {
      //permutations to split the quad
      int perm0[]={0,1,2}; // nodes of the first tri
      int perm1[]={0,2,3}; // nodes of the second tri

      // first tri
      for (int i=0; i<3; ++i)
      {
        x(i,0) = element->Nodes()[perm0[i]]->X()[0];
        x(i,1) = element->Nodes()[perm0[i]]->X()[1];
        x(i,2) = element->Nodes()[perm0[i]]->X()[2];
      }
      ComputeNormalCenter(x,n,c);

      // plug everything into the data maps
      normals->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(n)) ));
      centers->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(c)) ));

      // second tri
      for (int i=0; i<3; ++i)
      {
        x(i,0) = element->Nodes()[perm1[i]]->X()[0];
        x(i,1) = element->Nodes()[perm1[i]]->X()[1];
        x(i,2) = element->Nodes()[perm1[i]]->X()[2];
      }
      ComputeNormalCenter(x,n,c);

      // plug everything into the data maps
      normals->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid+gnumele,Teuchos::rcp(new Epetra_SerialDenseMatrix(n)) ));
      centers->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid+gnumele,Teuchos::rcp(new Epetra_SerialDenseMatrix(c)) ));
    }
    else
      dserror("only linear tris or bilinear quads are processed in here");
  }
}

void  STR::INVANA::ObjectiveFunctSurfCurr::ComputeNormalCenterSpatialConfig(DRT::Condition* cond,
                                                                            Teuchos::RCP<DRT::Discretization> discret,
                                                                            std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* normals,
                                                                            std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* centers,
                                                                            std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >* derivnormal,
                                                                            std::map<int, std::vector<int> >* facetdofmap,
                                                                            bool wantderiv)
{
  // clear data maps
  normals->clear();
  centers->clear();
  derivnormal->clear();
  facetdofmap->clear();


  // get state of displacement
  if (!discret->HasState("displacements"))
    dserror("state needs to be set in advance");

  Teuchos::RCP<const Epetra_Vector> disp = discret->GetState("displacements");
  disp = discret->GetState("displacements");

  // get total number of elements in this condition
  int lnumele = cond->Geometry().size();
  int gnumele;
  cond->Comm()->SumAll(&lnumele,&gnumele,1);

  // get the geometry and loop the elements
  std::map<int,RCP<DRT::Element> >& geom = cond->Geometry();
  std::map<int,RCP<DRT::Element> >::iterator ele;
  for (ele=geom.begin(); ele != geom.end(); ++ele)
  {
    DRT::Element* element = ele->second.get();
    int numnodes = element->NumNode();
    int egid = element->Id();

    //get this element's dofs
    DRT::Element::LocationArray la(discret->NumDofSets());
    element->LocationVector(*discret,la,false);

    // this element's displacements
    std::vector<double> mydisp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,la[0].lm_);

    // fill maps of normals and centers with data:
    // case tri  -> no special treatment
    // case quad -> split into 2 tris and associate to the second tri a gid as gid=egid+maxele
    LINALG::SerialDenseMatrix x(3,3); // node coordinates
    Epetra_SerialDenseMatrix c(3,1); // facet center
    Epetra_SerialDenseMatrix n(3,1); // facet normal
    Epetra_SerialDenseMatrix dn(3,9); // facet normal derivative
    std::vector<int> dofs; // each facets dof gids
    if (numnodes==3)
    {
      for (int i=0; i<3; ++i)
      {
        x(i,0) = element->Nodes()[i]->X()[0]+mydisp[i*3+0];
        x(i,1) = element->Nodes()[i]->X()[1]+mydisp[i*3+1];
        x(i,2) = element->Nodes()[i]->X()[2]+mydisp[i*3+2];

        if (wantderiv)
        {
          dofs.push_back(la[0].lm_[i*3+0]);
          dofs.push_back(la[0].lm_[i*3+1]);
          dofs.push_back(la[0].lm_[i*3+2]);
        }

      }
      ComputeNormalCenter(x,n,c);
      if (wantderiv) ComputeNormalDeriv(x,dn);

      normals->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(n)) ));
      centers->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(c)) ));
      if (wantderiv)
      {
        derivnormal->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(dn)) ));
        facetdofmap->insert(std::pair<int,std::vector<int> >(egid,dofs));
      }

    }
    else if (numnodes==4)
    {
      //permutations to split the quad
      int perm0[]={0,1,2}; // nodes of the first tri
      int perm1[]={0,2,3}; // nodes of the second tri

      // first tri
      for (int i=0; i<3; ++i)
      {
        x(i,0) = element->Nodes()[perm0[i]]->X()[0]+mydisp[perm0[i]*3+0];
        x(i,1) = element->Nodes()[perm0[i]]->X()[1]+mydisp[perm0[i]*3+1];
        x(i,2) = element->Nodes()[perm0[i]]->X()[2]+mydisp[perm0[i]*3+2];

        if (wantderiv)
        {
          dofs.push_back(la[0].lm_[perm0[i]*3+0]);
          dofs.push_back(la[0].lm_[perm0[i]*3+1]);
          dofs.push_back(la[0].lm_[perm0[i]*3+2]);
        }
      }
      ComputeNormalCenter(x,n,c);
      if (wantderiv) ComputeNormalDeriv(x,dn);

      normals->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(n)) ));
      centers->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(c)) ));
      if (wantderiv)
      {
        derivnormal->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid,Teuchos::rcp(new Epetra_SerialDenseMatrix(dn)) ));
        facetdofmap->insert(std::pair<int,std::vector<int> >(egid,dofs));
      }

      dofs.clear();
      // second tri
      for (int i=0; i<3; ++i)
      {
        x(i,0) = element->Nodes()[perm1[i]]->X()[0]+mydisp[perm1[i]*3+0];
        x(i,1) = element->Nodes()[perm1[i]]->X()[1]+mydisp[perm1[i]*3+1];
        x(i,2) = element->Nodes()[perm1[i]]->X()[2]+mydisp[perm1[i]*3+2];

        if (wantderiv)
        {
          dofs.push_back(la[0].lm_[perm1[i]*3+0]);
          dofs.push_back(la[0].lm_[perm1[i]*3+1]);
          dofs.push_back(la[0].lm_[perm1[i]*3+2]);
        }
      }
      ComputeNormalCenter(x,n,c);
      if (wantderiv) ComputeNormalDeriv(x,dn);

      // plug everything in to the data maps
      normals->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid+gnumele,Teuchos::rcp(new Epetra_SerialDenseMatrix(n)) ));
      centers->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid+gnumele,Teuchos::rcp(new Epetra_SerialDenseMatrix(c)) ));
      if (wantderiv)
      {
        derivnormal->insert(std::pair<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >(egid+gnumele,Teuchos::rcp(new Epetra_SerialDenseMatrix(dn)) ));
        facetdofmap->insert(std::pair<int,std::vector<int> >(egid+gnumele,dofs));
      }

    }
    else
      dserror("only linear tris or bilinear quads are processed in here");
  }
}

void STR::INVANA::ObjectiveFunctSurfCurr::ComputeNormalCenter(LINALG::SerialDenseMatrix& x,
                                                              Epetra_SerialDenseMatrix& n,
                                                              Epetra_SerialDenseMatrix& c)
{
  // normal
  n(0,0)=0.5*((x(2,1)-x(0,1))*(x(0,2)-x(1,2))-(x(2,2)-x(0,2))*(x(0,1)-x(1,1)));
  n(1,0)=0.5*((x(2,2)-x(0,2))*(x(0,0)-x(1,0))-(x(2,0)-x(0,0))*(x(0,2)-x(1,2)));
  n(2,0)=0.5*((x(2,0)-x(0,0))*(x(0,1)-x(1,1))-(x(2,1)-x(0,1))*(x(0,0)-x(1,0)));

  // center
  c(0,0)=(x(0,0)+x(1,0)+x(2,0))/3;
  c(1,0)=(x(0,1)+x(1,1)+x(2,1))/3;
  c(2,0)=(x(0,2)+x(1,2)+x(2,2))/3;

}

void STR::INVANA::ObjectiveFunctSurfCurr::ComputeNormalDeriv(LINALG::SerialDenseMatrix& x,
                                                             Epetra_SerialDenseMatrix& dn)
{
  // derivative of normal (wich is the cross product of vector made up from point coordinates in x, see function
  // ComputeNormalCenter) wrt the single dofs
  // -> a 3by9 matrix [dN/du1 dN/du2 ... dN/du9]

  // dofs of the first node
  // dN/du1 -> dN/x(0,0)
  dn(0,0)=0.0;
  dn(1,0)=(x(2,2)-x(0,2))+(x(0,2)-x(1,2));
  dn(2,0)=-(x(0,1)-x(1,1))-(x(2,1)-x(0,1));

  // dN/du2 -> dN/x(0,1)
  dn(0,1)=-(x(0,2)-x(1,2))-(x(2,2)-x(0,2));
  dn(1,1)=0.0;
  dn(2,1)=(x(2,0)-x(0,0))+(x(0,0)-x(1,0));

  // dN/du3 -> dN/dx(0,2)
  dn(0,2)=(x(2,1)-x(0,1))+(x(0,1)-x(1,1));
  dn(1,2)=-(x(0,0)-x(1,0))-(x(2,0)-x(0,0));
  dn(2,2)=0.0;

  // dofs of the second node
  // dN/du4 -> dN/dx(1,0)
  dn(0,3)=-0.0;
  dn(1,3)=-(x(2,2)-x(0,2));
  dn(2,3)=(x(2,1)-x(0,1));

  // dN/du5 -> dN/dx(1,1)
  dn(0,4)=(x(2,2)-x(0,2));
  dn(1,4)=0.0;
  dn(2,4)=-(x(2,0)-x(0,0));

  // dN/du6 -> dN/dx(1,2)
  dn(0,5)=-(x(2,1)-x(0,1));
  dn(1,5)=(x(2,0)-x(0,0));
  dn(2,5)=0.0;

  // dofs of the third node
  // dN/du7 -> dN/dx(2,0)
  dn(0,6)=0.0;
  dn(1,6)=-(x(0,2)-x(1,2));
  dn(2,6)=(x(0,1)-x(1,1));

  // dN/du8 -> dN/dx(2,1)
  dn(0,7)=(x(0,2)-x(1,2));
  dn(1,7)=0.0;
  dn(2,7)=-(x(0,0)-x(1,0));

  // dN/du9 -> dN/dx(2,2)
  dn(0,8)=-(x(0,1)-x(1,1));
  dn(1,8)=(x(0,0)-x(1,0));
  dn(2,8)=0.0;

  dn.Scale(0.5);

}

void STR::INVANA::ObjectiveFunctSurfCurr::PrintDataToScreen(std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > data)
{
  std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator it;

  std::cout << "id data x y z" << std::endl;
  for (it=data.begin(); it!=data.end(); ++it)
  {
    cout << "gid: " << it->first << " data: " << (*it->second)(0,0) << " " << (*it->second)(1,0) << " " << (*it->second)(2,0) << endl;
  }
}

void STR::INVANA::ObjectiveFunctSurfCurr::PrintDataToScreen(std::map<int,std::vector<int> > data)
{
  std::map<int,std::vector<int> >::iterator it;

  for (it=data.begin(); it!=data.end(); ++it)
  {
    cout << "gid: " << it->first << " data: ";
    for (int j=0; j<(int)it->second.size(); j++)
      cout << (it->second)[j] << " ";

    cout << " " << endl;
  }
}

void STR::INVANA::ObjectiveFunctSurfCurr::PrintDataToScreen(Epetra_Vector& vec)
{
  cout << "print values of EpetraVector " << endl;
  for (int i=0; i<vec.MyLength(); i++)
  {
    printf("mypid: %2d %10.8e \n", vec.Comm().MyPID(), vec[i]);
  }
}

double STR::INVANA::ObjectiveFunctSurfCurr::Convolute(std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n1,
                                                      std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c1,
                                                      std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n2,
                                                      std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c2)
{
  double val = 0.0;

  if (sourcedis_->Comm().MyPID()==0)
  {
    std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator ito;
    std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator iti;
    //Epetra_SerialDenseMatrix dum(3,1);
    for (ito=n1.begin(); ito!=n1.end(); ++ito)
    {
      for (iti=n2.begin(); iti!=n2.end(); ++iti)
      {
        //double dot = dum.DOT(3,(*(ito->second))[0],(*(iti->second))[0],1,1);
        double dot = dotp(*(ito->second),*(iti->second));
        val+=dot*kernel(*c2[iti->first], *c1[ito->first]);
      }
    }
  }

  sourcedis_->Comm().Broadcast(&val,1,0);

  return val;
}

void STR::INVANA::ObjectiveFunctSurfCurr::DerivComponent2(Teuchos::RCP<Epetra_Vector> lgradient,
                                                          std::map<int,std::vector<int> >& facetdofmap,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n1,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c1,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n2,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c2)
{
  if (sourcedis_->Comm().MyPID()==0)
  {
    Epetra_SerialDenseMatrix dum(9,1);

    std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator ito;
    std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator iti;
    for (ito=n1.begin(); ito!=n1.end(); ++ito)
    {
      for (iti=n2.begin(); iti!=n2.end(); ++iti)
      {
        double dot = dum.DOT(3,(*(ito->second))[0],(*(iti->second))[0],1,1);

        // for the dofs of the inner loop facet
        dum = kernelderiv(*c2[iti->first], *c1[ito->first],0);
        dum.Scale(dot);
        //add result to the corresponding entries in the global gradient vector
        int err = lgradient->SumIntoGlobalValues(9,dum[0],facetdofmap[iti->first].data());
        if (err)
          dserror("one of these gids not living on this proc");

        // for the dofs of the outer loop facet
        dum = kernelderiv(*c2[iti->first], *c1[ito->first],1);
        dum.Scale(dot);
        //add result to the corresponding entries in the global gradient vector
        err = lgradient->SumIntoGlobalValues(9,dum[0],facetdofmap[ito->first].data());
        if (err)
          dserror("one of these gids not living on this proc");
      }
    }
  }
}

void STR::INVANA::ObjectiveFunctSurfCurr::DerivComponent3(Teuchos::RCP<Epetra_Vector> lgradient,
                                                          std::map<int,std::vector<int> >& facetdofmap,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n1,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c1,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n2,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c2)
{
  if (sourcedis_->Comm().MyPID()==0)
  {
    Epetra_SerialDenseMatrix dum(9,1);

    std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator ito;
    std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator iti;
    for (ito=n1.begin(); ito!=n1.end(); ++ito)
    {
      for (iti=n2.begin(); iti!=n2.end(); ++iti)
      {
        double dot = dum.DOT(3,(*(ito->second))[0],(*(iti->second))[0],1,1);

        // for the dofs of the outer loop facet
        dum = kernelderiv(*c2[iti->first], *c1[ito->first],1);
        dum.Scale(dot);
        //add result to the corresponding entries in the global gradient vector
        int err = lgradient->SumIntoGlobalValues(9,dum[0],facetdofmap[ito->first].data());
        if (err)
          dserror("one of these gids not living on this proc");
      }
    }
  }
}

void STR::INVANA::ObjectiveFunctSurfCurr::DerivComponent1(Teuchos::RCP<Epetra_Vector> lgradient,
                                                          std::map<int,std::vector<int> >& facetdofmap,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& dn1,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c1,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& n2,
                                                          std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >& c2)
{
  if (sourcedis_->Comm().MyPID()==0)
  {
    Epetra_SerialDenseMatrix dum(9,1);

    std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator ito;
    std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> >::iterator iti;
    for (ito=dn1.begin(); ito!=dn1.end(); ++ito)
    {
      for (iti=n2.begin(); iti!=n2.end(); ++iti)
      {
        dum.Multiply('T','N',1.0,*(ito->second), *(iti->second),0.0);
        dum.Scale(kernel(*c2[iti->first], *c1[ito->first]));

        //add result to the corresponding entries in the global gradient vector
        int err = lgradient->SumIntoGlobalValues(9,dum[0],facetdofmap[ito->first].data());
        if (err)
          dserror("one of these gids not living on this proc");
      }
    }
  }
}

double STR::INVANA::ObjectiveFunctSurfCurr::WSpaceNorm()
{
  // the structural integrity component of the target is already precomputed
  double val=structinttarget_;

  //evaluate centers and normals in the spatial configuration (no derivatives at this point)
  ComputeNormalCenterSpatialConfig(sourcecond_,sourcedis_,&normal_s_,&center_s_,&derivnormal_s_,&facetdofmap_,false);

  //communicate the data maps to proc0
  DRT::Exporter ex(*sourcemap_,*sourcemapred_,sourcedis_->Comm());
  ex.Export(normal_s_);
  ex.Export(center_s_);

  //do the looping and plug result into val.
  val+=Convolute(normal_s_,center_s_,normal_s_,center_s_);
  val-=2*Convolute(normal_s_,center_s_,normal_t_,center_t_);

  return val;
}

void STR::INVANA::ObjectiveFunctSurfCurr::GradientWSpaceNorm(Teuchos::RCP<Epetra_MultiVector> gradient, int vind)
{
  // we need a dof row map reduced to proc0 to be able to fill gradient data into the correct
  // global dof and then distribute it afterwards
  Teuchos::RCP<Epetra_Map> dofrowmapred =  LINALG::AllreduceEMap(*(sourcedis_->DofRowMap()),0);

//  cout << "original target dofrowmap: " << *(targetdis_->DofRowMap()) << endl;
//  cout << "original source dofrowmap: " << *(sourcedis_->DofRowMap()) << endl;
//  cout << "reduced dofrowmap: " << *dofrowmapred << endl;

  // the gradient living on proc0 only
  Teuchos::RCP<Epetra_Vector> lgradient = Teuchos::rcp(new Epetra_Vector(*dofrowmapred,true));
  // temporal helpers
  Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*dofrowmapred,true));
  Teuchos::RCP<Epetra_Vector> tempdist = Teuchos::rcp(new Epetra_Vector(*(sourcedis_->DofRowMap()),false));

  // evaluate centers, normals and derivatives of the source in the spatial configuration
  ComputeNormalCenterSpatialConfig(sourcecond_,sourcedis_,&normal_s_,&center_s_,&derivnormal_s_,&facetdofmap_,true);

  //cout << "source stuff" << endl;
  //PrintDataToScreen(center_s_);
  //PrintDataToScreen(normal_s_);

  //PrintDataToScreen(facetdofmap_);

  //communicate the data maps to proc0
  DRT::Exporter ex(*sourcemap_,*sourcemapred_,sourcedis_->Comm());
  ex.Export(normal_s_);
  ex.Export(center_s_);
  ex.Export(derivnormal_s_);
  ex.Export(facetdofmap_);

  //maps of the target already live on proc0 since construction

  // compute the different summands making up the gradient
  DerivComponent1(lgradient,facetdofmap_,derivnormal_s_,center_s_,normal_s_,center_s_);
  //DerivComponent1(lgradient,facetdofmap_,derivnormal_s_,center_s_,normal_s_,center_s_);
  lgradient->Scale(2.0);
  DerivComponent2(lgradient,facetdofmap_,normal_s_,center_s_,normal_s_,center_s_);
  DerivComponent1(temp,facetdofmap_,derivnormal_s_,center_s_,normal_t_,center_t_);
  lgradient->Update(-2.0,*temp,1.0);
  temp->PutScalar(0.0);
  DerivComponent3(temp,facetdofmap_,normal_s_,center_s_,normal_t_,center_t_);
  lgradient->Update(-2.0,*temp,1.0);

  //PrintDataToScreen(*lgradient);

  LINALG::Export(*lgradient, *tempdist);
// !!!!!!!!!!! THE ADD COMBO METHOD DOES NOT WORK!!!!!!!!!!!!!!!
// SO ONE CAN USE THE LINALG VERSION BUT THAT NEED AOTHER TEMPORAL VECTOR SINCE ONLY "INSERT" IS PROVIDED
//  // import the gradient according to the global problems dofrowmap
//  Epetra_Import im = Epetra_Import(gradient->Map(),*dofrowmapred);
//  int err = (*gradient)(vind)->Import(*lgradient,im,Insert);
//  if (err)
//    dserror("Gradient distribution failed");

  // set into the gradient (update here, since there migth be other pairs of surfaces which contribute)
  (*gradient)(vind)->Update(1.0,*tempdist,1.0);

}

inline double STR::INVANA::ObjectiveFunctSurfCurr::kernel(Epetra_SerialDenseMatrix x, Epetra_SerialDenseMatrix y)
{
  return exp(-((x(0,0)-y(0,0))*(x(0,0)-y(0,0))+(x(1,0)-y(1,0))*(x(1,0)-y(1,0))+(x(2,0)-y(2,0))*(x(2,0)-y(2,0)))/sigmaW2_);
}

Epetra_SerialDenseMatrix STR::INVANA::ObjectiveFunctSurfCurr::kernelderiv(Epetra_SerialDenseMatrix x, Epetra_SerialDenseMatrix y, int fid)
{
  double val=-(2.0/(3.0*sigmaW2_))*kernel(x,y);

  Epetra_SerialDenseMatrix deriv(9,1);
  Epetra_SerialDenseMatrix diff = x;
  y.Scale(-1.0);
  diff+=y;

  for (int i=0; i<3; i++)
  {
    deriv(i*3+0,0) = diff(0,0);
    deriv(i*3+1,0) = diff(1,0);
    deriv(i*3+2,0) = diff(2,0);
  }

  if (fid==0)
    deriv.Scale(val);
  else if (fid==1)
    deriv.Scale(-val);
  else
    dserror("give fid as 0 or 1 only!");

  return deriv;
}

inline double STR::INVANA::ObjectiveFunctSurfCurr::dotp(Epetra_SerialDenseMatrix x, Epetra_SerialDenseMatrix y)
{
  return x(0,0)*y(0,0)+x(1,0)*y(1,0)+x(2,0)*y(2,0);
}
