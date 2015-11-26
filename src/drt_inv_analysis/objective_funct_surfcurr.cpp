#if defined( HAVE_Kokkos )
/*----------------------------------------------------------------------*/
/*!
 * \file objective_funct_surfcurr.cpp

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#include "objective_funct_surfcurr.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_io/io_pstream.H"

#include <sys/time.h>


/*----------------------------------------------------------------------*/
/* standard constructor of current representation                       */
/*----------------------------------------------------------------------*/
INVANA::SurfCurrentGroup::SurfCurrentGroup(Teuchos::RCP<DRT::Discretization> discret):
sourcedis_(discret),
targetdis_(Teuchos::null)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

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
  std::vector<Teuchos::RCP<SurfCurrentPair> > dump;
  currents_.push_back(dump);
  for (int i=0; i<(int)scc_source.size(); i++)
  {
    bool foundit=false;
    int ids = scc_source[i]->GetInt("matching id");
    for (int j=0; j<(int)scc_target.size(); j++)
    {
      int idt = scc_target[j]->GetInt("matching id");
      if (idt==ids)
      {
        // add to "0" since the loop of multiple measurements over timestill missing
        currents_[0].push_back(Teuchos::rcp(new SurfCurrentPair(sourcedis_,targetdis_,scc_source[i],scc_target[j])));
        foundit=true;
        break;
      }
    }

    if (foundit)
      continue;
    else
      dserror("corresponding condition in target not found");
  }

  // currents so far only for the last time step;
  timesteps_.push_back(sdyn.get<double>("MAXTIME"));

  if (currents_[0].size()!=scc_source.size())
    dserror("problem in finding corresponding current conditions in source and target discretization");

  //initialize parallel environment for the computation of each surface current pair
  Kokkos::initialize();
#if defined( KOKKOS_HAVE_OPENMP )
  std::cout << "Surface Currents in OpenMP-mode" << std::endl;
#else
  std::cout << "Surface Currents in Serial-mode" << std::endl;
#endif

}

/*----------------------------------------------------------------------*/
/* find step of measurement according to given time          keh 10/14  */
/*----------------------------------------------------------------------*/
int INVANA::SurfCurrentGroup::FindStep(double time)
{
  // find step of the evaluation according to time:
  int step=-1;
  for (int i=0; i<(int)timesteps_.size(); i++)
  {
    double dt=abs(timesteps_[i]-time);
    if (dt<1.0e-10)
      step=i;
  }

  return step;
}

/*----------------------------------------------------------------------*/
/* Evaluate value of the objective function                  keh 11/13  */
/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentGroup::Evaluate(Teuchos::RCP<Epetra_Vector> state,
    double time,
    double& val)
{
  int step = FindStep(time); // not needed so far
  if (step == -1)
    dserror("no measurements for the requested time step");

  val = 0.0;

  //so far only the case were one single "measurement" for the final timestep of the simulation is considered
  sourcedis_->SetState("displacements", state);

  struct timeval begin,end;
  gettimeofday(&begin,NULL);

  for(int i=0; i<(int)currents_[0].size(); ++i)
  {
    // evaluate every single surface combination
    val+=currents_[0][i]->WSpaceNorm();
  }
  gettimeofday(&end,NULL);
  double dtime = 1.0*(end.tv_sec-begin.tv_sec) + 1.0e-6*(end.tv_usec-begin.tv_usec);

  if (sourcedis_->Comm().MyPID()==0)
    IO::cout << "SURFACE CURRENT function evaluation took: " << dtime << " seconds" << IO::endl;
}

/*----------------------------------------------------------------------*/
/* Evaluate the gradient of the objective function                      */
/* w.r.t the displacements                                   keh 11/13  */
/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentGroup::EvaluateGradient(Teuchos::RCP<Epetra_Vector> state,
    double time,
    Teuchos::RCP<Epetra_Vector> gradient)
{
  int step = FindStep(time); // not needed so far
  if (step == -1)
    dserror("no measurements for the requested time step");

  gradient->PutScalar(0.0);

  sourcedis_->SetState("displacements", state);

  struct timeval begin,end;
  gettimeofday(&begin,NULL);

  for(int i=0; i<(int)currents_[0].size(); ++i)
  {
    // evaluate every single surface combination
    currents_[0][i]->GradientWSpaceNorm(gradient);
  }
  gettimeofday(&end,NULL);
  double dtime = 1.0*(end.tv_sec-begin.tv_sec) + 1.0e-6*(end.tv_usec-begin.tv_usec);

  if (sourcedis_->Comm().MyPID()==0)
    IO::cout << "SURFACE CURRENT gradient evaluation took: " << dtime << " seconds" << IO::endl;
}

void INVANA::SurfCurrentGroup::SetScale(double sigmaW)
{
  for(int i=0; i<(int)currents_[0].size(); ++i)
  {
    currents_[0][i]->SetScale(sigmaW);
  }
}


/*----------------------------------------------------------------------*/
/* standard constructor for a surface current                keh 11/13  */
/*----------------------------------------------------------------------*/
INVANA::SurfCurrentPair::SurfCurrentPair(
    Teuchos::RCP<DRT::Discretization> sourcedis,
    Teuchos::RCP<DRT::Discretization> targetdis,
    DRT::Condition* sourcecond,
    DRT::Condition* targetcond)
{
  // get the scale of the kernel
  const Teuchos::ParameterList& statinvp = DRT::Problem::Instance()->StatInverseAnalysisParams();
  double sigmaW = statinvp.get<double>("KERNELSCALE");
  if (sigmaW<0.0) dserror("you need to choose a proper scale (KERNELSCALE) at which to evaluate the current");
  sigmaW2_=sigmaW*sigmaW;

  scaling_ = DRT::INPUT::IntegralValue<bool>(statinvp, "OBJECTIVEFUNCTSCAL");

  // Setup the triangulations
  tri_target_ = Teuchos::rcp(new INVANA::Triangulation(targetdis,Teuchos::rcp(targetcond,false)));
  tri_source_ = Teuchos::rcp(new INVANA::Triangulation(sourcedis,Teuchos::rcp(sourcecond,false)));

  extract_type x_target = tri_target_->Points();
  int xsize=x_target.size();
  trimesh_type x_target_view("X_View_t", xsize);
  trimesh_host_type h_x_target_view = Kokkos::create_mirror_view (x_target_view);

  // Fill View and Copy to device
  ExtractToHView(x_target, h_x_target_view);
  Kokkos::deep_copy (x_target_view, h_x_target_view);

  // normals
  currents_type n_target("N_View_t",xsize);
  Kokkos::parallel_for(xsize,Normal(n_target, x_target_view));
  n_target_=n_target;

  // centers
  currents_type c_target("C_View_t",xsize);
  Kokkos::parallel_for(xsize,Centers(c_target, x_target_view));
  c_target_=c_target;

  double sum=0.0;
  Kokkos::parallel_reduce(xsize,Convolute<currents_type>(c_target,n_target,c_target,n_target,sigmaW2_),sum);

  // since the work was done redundantly by all mpi-ranks
  // this value should be the same for all
  preconvtarget_ = sum;
}

/*----------------------------------------------------------------------*/
double INVANA::SurfCurrentPair::WSpaceNorm()
{
  int myrank = tri_source_->Comm()->MyPID();
  int numrnk = tri_source_->Comm()->NumProc();

  // the structural integrity component of the target is already precomputed
  double val=preconvtarget_;

  // let the triangulation be warped according to the state of whatever, e.g the pde
  // solved by the 'discretization'
  tri_source_->EvaluateWarp();

  //-------------------------------------------------
  // All the points and the data reduced to myrank
  extract_type x_source = tri_source_->Points(myrank);
  extract_type disp = tri_source_->Data(myrank);
  unsigned xsize=x_source.size();

  // View on all points of the source
  trimesh_type x_source_view("X_View_s", xsize);
  trimesh_host_type h_x_source_view = Kokkos::create_mirror_view (x_source_view);
  ExtractToHView(x_source, h_x_source_view);
  Kokkos::deep_copy (x_source_view, h_x_source_view);

  // View on all the data associated to the points
  trimesh_type d_source_view("D_View",xsize);
  trimesh_host_type h_d_source_view = Kokkos::create_mirror_view (d_source_view);
  ExtractToHView(disp,h_d_source_view);
  Kokkos::deep_copy (d_source_view, h_d_source_view);

  //-------------------------------------------------
  // Normals and Centers
  // push forward points
  Kokkos::parallel_for(xsize,Update(x_source_view,d_source_view));

  // normals
  currents_type n_source("N_View_s",xsize);
  Kokkos::parallel_for(xsize,Normal(n_source, x_source_view));

  // centers
  currents_type c_source("C_View_s",xsize);
  Kokkos::parallel_for(xsize,Centers(c_source, x_source_view));

  //-------------------------------------------------
  // chunks of the data for myrank
  std::pair<currents_type::size_type, currents_type::size_type> mychunk = MyIndices(numrnk,myrank,xsize);
  auto my_c_source = Kokkos::subview(c_source, mychunk, Kokkos::ALL());
  auto my_n_source = Kokkos::subview(n_source, mychunk, Kokkos::ALL());
  int my_xsize = my_c_source.dimension_0();

  //-------------------------------------------------
  // Do the convolutions
  double sum=0.0;
  Kokkos::parallel_reduce(my_xsize,Convolute<ViewStride>(my_c_source,my_n_source,c_source,n_source,sigmaW2_),sum);
  double gsum=0.0;
  tri_source_->Comm()->SumAll(&sum,&gsum,1);
  val+=gsum;

  sum=0.0;
  Kokkos::parallel_reduce(my_xsize,Convolute<ViewStride>(my_c_source,my_n_source,c_target_,n_target_,sigmaW2_),sum);
  gsum=0.0;
  tri_source_->Comm()->SumAll(&sum,&gsum,1);
  val-=2*gsum;

  //scaling by number of surface elements in source and target
  if (scaling_)
  {
    double fac = xsize+c_target_.dimension_0();
    val=val/fac;
  }

  return val;
}

/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentPair::GradientWSpaceNorm(Teuchos::RCP<Epetra_MultiVector> gradient)
{
  int myrank = tri_source_->Comm()->MyPID();
  int numrnk = tri_source_->Comm()->NumProc();

  // let the triangulation be warped according to the state of whatever, e.g the pde
  // solved by the 'discretization'
  tri_source_->EvaluateWarp();

  //-------------------------------------------------
  // All the points and the data reduced to myrank
  extract_type x_source = tri_source_->Points(myrank);
  extract_type disp = tri_source_->Data(myrank);
  unsigned xsize=x_source.size();

  // View on all points of the source
  trimesh_type x_source_view("X_View_s", xsize);
  trimesh_host_type h_x_source_view = Kokkos::create_mirror_view (x_source_view);
  Teuchos::RCP<Epetra_Map> trimap=Teuchos::rcp(new Epetra_Map(0,0,*tri_source_->Comm()));
  ExtractToHView(x_source, h_x_source_view, trimap);
  Kokkos::deep_copy (x_source_view, h_x_source_view);

  // View on all the data associated to the points
  trimesh_type d_source_view("D_View",xsize);
  trimesh_host_type h_d_source_view = Kokkos::create_mirror_view (d_source_view);
  ExtractToHView(disp,h_d_source_view);
  Kokkos::deep_copy (d_source_view, h_d_source_view);

  //-------------------------------------------------
  // Normals and Centers
  // push forward points
  Kokkos::parallel_for(xsize,Update(x_source_view,d_source_view));

  // normals
  currents_type n_source("N_View_s",xsize);
  Kokkos::parallel_for(xsize,Normal(n_source, x_source_view));

  // derivative of the normals
  trimesh_vecdata_type dn_source("DN_View_s",xsize);
  Kokkos::parallel_for(xsize,DNormal(x_source_view,dn_source));

  // centers
  currents_type c_source("C_View_s",xsize);
  Kokkos::parallel_for(xsize,Centers(c_source, x_source_view));

  // Initialize the gradient
  trimesh_type grad("G_View",xsize);
  Kokkos::parallel_for(xsize,Init9(grad));

  //-------------------------------------------------
  // chunks of the data for myrank
  std::pair<currents_type::size_type, currents_type::size_type> mychunk = MyIndices(numrnk,myrank,xsize);
  auto my_c_source = Kokkos::subview(c_source, mychunk, Kokkos::ALL());
  auto my_n_source = Kokkos::subview(n_source, mychunk, Kokkos::ALL());
  auto my_dn_source = Kokkos::subview(dn_source,mychunk, Kokkos::ALL());
  int my_xsize = my_c_source.dimension_0();

  int off = mychunk.first;
  Kokkos::parallel_for(my_xsize,ConvoluteDN<ViewStride,ViewStride>(my_c_source,my_dn_source,c_source,n_source,grad,off,sigmaW2_,2.0));
  Kokkos::parallel_for(my_xsize,ConvoluteDN<ViewStride,ViewStride>(my_c_source,my_dn_source,c_target_,n_target_,grad,off,sigmaW2_,-2.0));
  Kokkos::parallel_for(my_xsize, ConvoluteDk<ViewStride>(my_c_source,my_n_source,c_source,n_source,grad,off,sigmaW2_,1.0,true));
  Kokkos::parallel_for(my_xsize, ConvoluteDk<ViewStride>(my_c_source,my_n_source,c_target_,n_target_,grad,off,sigmaW2_,-2.0,false));

  // Bring back to host
  trimesh_host_type h_grad = Kokkos::create_mirror_view(grad);
  Kokkos::deep_copy(h_grad,grad);

  // restore in extract_type format to pass it to the SetData routine
  extract_type grad_data;
  HViewToExtract(h_grad,*trimap,grad_data);

  // Set gradient data into the global gradient vector
  Teuchos::RCP<Epetra_MultiVector> lgradient = Teuchos::rcp(new Epetra_MultiVector(gradient->Map(),true));
  tri_source_->SetDataGlobally(grad_data,*trimap,lgradient);

  //scaling by number of surface elements in source and target
  if (scaling_)
  {
    double fac = xsize+c_target_.dimension_0();
    tri_source_->Comm()->Broadcast(&fac,1,0);
    gradient->Update(1.0/fac,*lgradient,1.0);
  }
  else
    gradient->Update(1.0,*lgradient,1.0);
}

/*----------------------------------------------------------------------*/
std::pair<int,int> INVANA::SurfCurrentPair::MyIndices(int nrnk, int myrnk, int size)
{
  int l = size%nrnk;
  int s = (int)(size-l)/nrnk;

  int lower = myrnk*s;
  int upper = (myrnk+1)*s;

  // for the "last" rank
  if (myrnk==(nrnk-1))
    upper = size;

  return std::make_pair(lower, upper);
}

/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentPair::ExtractToHView(
    const extract_type& in,
    trimesh_host_type& out)
{
  dsassert(in.size()!=out.dimension_0(),"dimension mismatch");

  // fill the data in the view
  int i=0;
  for (auto it=in.begin(); it!=in.end(); it++)
  {
    out(i,0) = it->second[0];
    out(i,1) = it->second[1];
    out(i,2) = it->second[2];
    out(i,3) = it->second[3];
    out(i,4) = it->second[4];
    out(i,5) = it->second[5];
    out(i,6) = it->second[6];
    out(i,7) = it->second[7];
    out(i,8) = it->second[8];
    i++;
  }
}

/*----------------------------------------------------------------------*/
void INVANA::SurfCurrentPair::ExtractToHView(
    const extract_type& in,
    trimesh_host_type& out,
    Teuchos::RCP<Epetra_Map> map)
{
  dsassert(in.size()!=out.dimension_0(),"dimension mismatch");

  std::vector<int> gids;

  // fill the data in the view
  int i=0;
  for (auto it=in.begin(); it!=in.end(); it++)
  {
    gids.push_back(it->first);
    out(i,0) = it->second[0];
    out(i,1) = it->second[1];
    out(i,2) = it->second[2];
    out(i,3) = it->second[3];
    out(i,4) = it->second[4];
    out(i,5) = it->second[5];
    out(i,6) = it->second[6];
    out(i,7) = it->second[7];
    out(i,8) = it->second[8];
    i++;
  }

  *map=Epetra_Map(-1,(int)gids.size(),gids.data(),0,map->Comm());
}

void INVANA::SurfCurrentPair::HViewToExtract(const trimesh_host_type& in, const Epetra_Map& inmap, extract_type& out)
{
  int size0=in.dimension_0();
  std::vector<double> dum;
  for (int i=0; i<size0; i++)
  {
    for(unsigned j=0; j<in.dimension_1(); j++)
      dum.push_back(in(i,j));

    out.insert(std::pair<int, std::vector<double> >(inmap.GID(i),dum));
    dum.clear();
  }
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::EvaluatePoints()
{
  //clear data
  points_.clear();
  std::vector<int> gids;

  int lnumele = surface_->Geometry().size();
  int gnumele;
  surface_->Comm()->SumAll(&lnumele,&gnumele,1);

  for (auto ele=surface_->Geometry().begin(); ele != surface_->Geometry().end(); ++ele)
  {
    Teuchos::RCP<DRT::Element> element = ele->second;

    // exclude non row map elements to be processed
    Teuchos::RCP<DRT::FaceElement> elef = Teuchos::rcp_dynamic_cast<DRT::FaceElement>(element);
    if (discret_->ElementRowMap()->LID(elef->ParentElementId()) == -1) continue;

    int numnodes = element->NumNode();
    int egid = element->Id();

    // fill maps of normals and centers with data:
    // case tri  -> no special treatment
    // case quad -> split into 2 tris and associate to the second tri a gid as gid=egid+maxele
    std::vector<double> x; //node coordinates

    if (numnodes==3)
    {
      for (int i=0; i<3; ++i)
      {
        x.push_back(element->Nodes()[i]->X()[0]);
        x.push_back(element->Nodes()[i]->X()[1]);
        x.push_back(element->Nodes()[i]->X()[2]);
      }
      points_.insert(std::pair<int, std::vector<double> >(egid,x));
      gids.push_back(egid);

    }
    else if (numnodes==4)
    {
      //permutations to split the quad
      int perm0[]={0,1,2}; // nodes of the first tri
      int perm1[]={0,2,3}; // nodes of the second tri

      // first tri
      for (int i=0; i<3; ++i)
      {
        x.push_back(element->Nodes()[perm0[i]]->X()[0]);
        x.push_back(element->Nodes()[perm0[i]]->X()[1]);
        x.push_back(element->Nodes()[perm0[i]]->X()[2]);
      }
      points_.insert(std::pair<int, std::vector<double> >(egid,x));
      gids.push_back(egid);

      x.clear();
      // second tri
      for (int i=0; i<3; ++i)
      {
        x.push_back(element->Nodes()[perm1[i]]->X()[0]);
        x.push_back(element->Nodes()[perm1[i]]->X()[1]);
        x.push_back(element->Nodes()[perm1[i]]->X()[2]);
      }
      points_.insert(std::pair<int, std::vector<double> >(egid+gnumele,x));
      gids.push_back(egid+gnumele);
    }
    else
      dserror("only linear tris or bilinear quads are processed in here");
  }
  // set up new map of triangles (which is unique)
  trimap_ = Teuchos::rcp(new Epetra_Map(-1,(int)gids.size(),gids.data(),0,*(surface_->Comm())));

  haspoints_=true;
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::EvaluateDofs()
{
  std::vector<int> gids;

  // get total number of face elements in this condition
  int lnumele = surface_->Geometry().size();
  int gnumele;
  surface_->Comm()->SumAll(&lnumele,&gnumele,1);

  // loop the surface elements and build dof map for the faces
  facetdofstype facetdofs;
  for (auto ele=surface_->Geometry().begin(); ele != surface_->Geometry().end(); ++ele)
  {
    Teuchos::RCP<DRT::Element> element = ele->second;
    int numnodes = element->NumNode();
    int egid = element->Id();

    //get this element's dofs
    DRT::Element::LocationArray la(discret_->NumDofSets());
    element->LocationVector(*discret_,la,false);

    if (numnodes==3)
    {
      facetdofs.insert(std::pair<int,std::vector<int> >(egid,la[0].lm_));
      gids.push_back(egid);
    }
    else if (numnodes==4)
    {
      //permutations to split the quad
      int perm0[]={0,1,2}; // nodes of the first tri
      int perm1[]={0,2,3}; // nodes of the second tri

      std::vector<int> dofs;

      // first tri
      for (int i=0; i<3; ++i)
      {
        dofs.push_back(la[0].lm_[perm0[i]*3+0]);
        dofs.push_back(la[0].lm_[perm0[i]*3+1]);
        dofs.push_back(la[0].lm_[perm0[i]*3+2]);
      }
      facetdofs.insert(std::pair<int,std::vector<int> >(egid,dofs));
      gids.push_back(egid);

      dofs.clear();
      // second tri
      for (int i=0; i<3; ++i)
      {
        dofs.push_back(la[0].lm_[perm1[i]*3+0]);
        dofs.push_back(la[0].lm_[perm1[i]*3+1]);
        dofs.push_back(la[0].lm_[perm1[i]*3+2]);
      }
      facetdofs.insert(std::pair<int,std::vector<int> >(egid+gnumele,dofs));
      gids.push_back(egid+gnumele);
    }
    else
      dserror("only linear tris or bilinear quads are processed in here");
  }
  // set up new map of triangles (which is not unique)
  Epetra_Map trimap(-1,(int)gids.size(),gids.data(),0,*(surface_->Comm()));

  // make it unique
  DRT::Exporter ex(trimap,*trimap_,trimap_->Comm());
  ex.Export(facetdofs);

  facetmap_=facetdofs;
  hasdofs_=true;
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::EvaluateWarp()
{
  //clear data
  data_.clear();
  std::vector<int> gids;

  // get state of displacement
  if (!discret_->HasState("displacements"))
    dserror("state needs to be set in advance to be able to gather data");

  Epetra_Vector disp(*(discret_->GetState("displacements")));

  // get total number of face elements in this condition
  int lnumele = surface_->Geometry().size();
  int gnumele;
  surface_->Comm()->SumAll(&lnumele,&gnumele,1);

  // triangulation
  std::map<int,Teuchos::RCP<DRT::Element> >& geom = surface_->Geometry();
  std::map<int,Teuchos::RCP<DRT::Element> >::iterator ele;
  for (ele=geom.begin(); ele != geom.end(); ++ele)
  {
    Teuchos::RCP<DRT::Element> element = ele->second;
    int numnodes = element->NumNode();
    int egid = element->Id();

    //get this element's dofs
    DRT::Element::LocationArray la(discret_->NumDofSets());
    element->LocationVector(*discret_,la,false);

    // this element's displacements
    std::vector<double> mydisp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(disp,mydisp,la[0].lm_);

    // fill maps of normals and centers with data:
    // case tri  -> no special treatment
    // case quad -> split into 2 tris and associate to the second tri a gid as gid=egid+maxele
    std::vector<double> x; //node data

    if (numnodes==3)
    {
      for (int i=0; i<3; ++i)
      {
        x.push_back(mydisp[i*3+0]);
        x.push_back(mydisp[i*3+1]);
        x.push_back(mydisp[i*3+2]);
      }
      data_.insert(std::pair<int, std::vector<double> >(egid,x));
      gids.push_back(egid);

    }
    else if (numnodes==4)
    {
      //permutations to split the quad
      int perm0[]={0,1,2}; // nodes of the first tri
      int perm1[]={0,2,3}; // nodes of the second tri

      // first tri
      for (int i=0; i<3; ++i)
      {
        x.push_back(mydisp[perm0[i]*3+0]);
        x.push_back(mydisp[perm0[i]*3+1]);
        x.push_back(mydisp[perm0[i]*3+2]);
      }
      data_.insert(std::pair<int, std::vector<double> >(egid,x));
      gids.push_back(egid);

      x.clear();
      // second tri
      for (int i=0; i<3; ++i)
      {
        x.push_back(mydisp[perm1[i]*3+0]);
        x.push_back(mydisp[perm1[i]*3+1]);
        x.push_back(mydisp[perm1[i]*3+2]);
      }
      data_.insert(std::pair<int, std::vector<double> >(egid+gnumele,x));
      gids.push_back(egid+gnumele);
    }
    else
      dserror("only linear tris or bilinear quads are processed in here");
  }
  // set up new map of triangles (which is not unique)
  Epetra_Map trimap(-1,(int)gids.size(),gids.data(),0,*(surface_->Comm()));

  // make it unique
  DRT::Exporter ex(trimap,*trimap_,trimap_->Comm());
  ex.Export(data_);

  hasdata_=true;
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::SetDataGlobally(const extract_type& data,const Epetra_Map& datamap, Teuchos::RCP<Epetra_MultiVector> vector)
{
  int myrank = Comm()->MyPID();
  // reduce assemble vector to myrank
  Epetra_Map bla(-1,vector->Map().NumMyElements(),vector->Map().MyGlobalElements(),0,vector->Comm());
  Epetra_Map dofmap_red(*(LINALG::AllreduceEMap(bla,myrank)));

  //reduce facetmap to myrank
  facetdofstype facetmap_red = facetmap_;
  DRT::Exporter ex1(*trimap_,datamap,trimap_->Comm());
  ex1.Export(facetmap_red);

  //assemble everything on myrank
  Epetra_Vector vector_red(dofmap_red,true);
  std::vector<int> dofs;
  for (auto it=data.begin(); it!=data.end(); it++)
  {
    auto facet = facetmap_red.find(it->first);
    if (facet == facetmap_red.end())
      dserror("Facet %d could not be found", it->first);

    dofs = facet->second;
    int err = vector_red.SumIntoGlobalValues(9,it->second.data(),dofs.data());
    if (err) dserror("Something went wrong!");

  }

  // bring to the MPI parallel layout
  Epetra_Export ex2(dofmap_red,vector->Map());
  int err2 = vector->Export(vector_red, ex2, Add);
  if (err2)
    dserror("export into global gradien layout failed with code: %d",err2);
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::Points(int rank)
{
  if (not haspoints_)
    dserror("Points were not evaluated so far!");

  // copy since we want to keep the original points
  extract_type p(points_);
  CommunicateData(p,rank);
  return p;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::Points()
{
  if (not haspoints_)
    dserror("Points were not evaluated so far!");

  // copy since we want to keep the original points
  extract_type p(points_);
  CommunicateData(p);
  return p;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::Data(int rank)
{
  if (not hasdata_)
    dserror("Data was not evaluated so far!");

  // copy since we want to keep the original data
  extract_type p(data_);
  CommunicateData(p,rank);
  return p;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::Data()
{
  if (not hasdata_)
    dserror("Data was not evaluated so far!");

  // copy since we want to keep the original data
  extract_type p(data_);
  CommunicateData(p);
  return p;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::MyPoints()
{
  if (not haspoints_)
    dserror("Points were not evaluated so far!");

  return points_;
}

/*----------------------------------------------------------------------*/
INVANA::extract_type INVANA::Triangulation::MyData()
{
  if (not hasdata_)
    dserror("Data was not evaluated so far!");

  return data_;
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::CommunicateData(extract_type& data,int rank)
{

  Epetra_Map map_red(*LINALG::AllreduceEMap(*trimap_, rank));

  // build the exporter
  DRT::Exporter ex(*trimap_,map_red,trimap_->Comm());

  // export
  ex.Export(data);
}

/*----------------------------------------------------------------------*/
void INVANA::Triangulation::CommunicateData(extract_type& data)
{

  Epetra_Map map_red(*LINALG::AllreduceEMap(*trimap_));

  // build the exporter
  DRT::Exporter ex(*trimap_,map_red,trimap_->Comm());

  // export
  ex.Export(data);
}

/*----------------------------------------------------------------------*/
const Teuchos::RCP<Epetra_Comm> INVANA::Triangulation::Comm() {return surface_->Comm();}
#endif
