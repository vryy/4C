/*!----------------------------------------------------------------------
\file pat_utils.cpp

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/staff/svenja-schoeder/
            089 - 289-15271
</pre>
*----------------------------------------------------------------------*/

#include "pat_utils.H"
#include "pat_imagereconstruction.H"

#include "Epetra_CrsMatrix.h"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_timestepping/timintmstep.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.cpp"
#include "acou_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_action.H"



/*----------------------------------------------------------------------*/
ACOU::PATSearchDirection::PATSearchDirection(INPAR::ACOU::OptimizationType opti)
{
  opti_ = opti;
}

/*----------------------------------------------------------------------*/
void ACOU::PATSearchDirection::Setup(const Epetra_Map* map, const Epetra_Map* uniquemap)
{
  if(opti_==INPAR::ACOU::inv_lbfgs)
  {
    actsize_ = 0;
    ssize_ = 10;

    sstore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, map, true));
    ystore_ = Teuchos::rcp(new TIMINT::TimIntMStep<Epetra_Vector>(0, 0, map, true));

    map_ = Teuchos::rcp(new Epetra_Map(*map));
    uniquemap_ = Teuchos::rcp(new Epetra_Map(*uniquemap));

    oldparams_ = Teuchos::rcp(new Epetra_Vector(*map_,false));
    oldgradient_ = Teuchos::rcp(new Epetra_Vector(*map_,false));
  }
  // nothing to do for steepest descent

  return;
}

/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ACOU::PATSearchDirection::ComputeDirection(Teuchos::RCP<Epetra_Vector> gradient, Teuchos::RCP<Epetra_Vector> params, int iter)
{
  Teuchos::RCP<Epetra_Vector> direction = Teuchos::rcp(new Epetra_Vector(gradient->Map(),false));
  if(opti_==INPAR::ACOU::inv_lbfgs)
  {
    if(iter==0)
    {
      oldgradient_->Update(1.0,*gradient,0.0);
      oldparams_->Update(1.0,*params,0.0);
      direction->Update(-1.0,*gradient,0.0);
    }
    else
    {
      // store vectors
      if(iter<=ssize_)
      {
        actsize_+=1;

        sstore_->Resize(-actsize_+1,0,map_.get(),false);
        ystore_->Resize(-actsize_+1,0,map_.get(),false);

        Epetra_Vector s(*map_,false);
        s.Update(1.0,*params,-1.0,*oldparams_,0.0);
        sstore_->UpdateSteps(s);

        s.Update(1.0,*gradient,-1.0,*oldgradient_,false);
        ystore_->UpdateSteps(s);
      }

      // compute direction
      direction->Update(1.0,*gradient,0.0);
      std::vector<double> alpha;

      // loop steps
      for (int i=0; i>-actsize_; i--)
      {
        double aa=0.0;
        double bb=0.0;

        ((*ystore_)(i))->Dot(*(*sstore_)(i),&aa);
        ((*sstore_)(i))->Dot((*direction),&bb);

        alpha.push_back(1/aa*bb);

        direction->Update(-1.0*alpha.back(), *(*ystore_)(i),1.0 );
      }

      for (int i=-actsize_+1; i<=0; i++)
      {
        double aa=0.0;
        double bb=0.0;
        double beta=0.0;

        ((*ystore_)(i))->Dot(*(*sstore_)(i),&aa);
        ((*ystore_)(i))->Dot((*direction),&bb);

        beta=1/aa*bb;
        double alphac=alpha.back();
        alpha.pop_back();

        direction->Update(alphac-beta, *(*sstore_)(i),1.0 );
      }

      // minimization not maximization
      direction->Scale(-1.0);
    }
  }
  else
  {
    if(0)
      direction->Update(-1000.0,*gradient,0.0);
    else
    {
      double maxval = 0.0;
      gradient->MaxValue(&maxval);
      double minval = 0.0;
      gradient->MinValue(&minval);

      if(abs(minval)>maxval && (minval>1.0e-10 || minval<-1.0e-10))
        direction->Update(-0.1/abs(minval),*gradient,0.0);
      else if( abs(maxval) > 1.0e-10)
        direction->Update(-0.1/abs(maxval),*gradient,0.0);
      else
      {
        double mean = 0.0;
        gradient->MeanValue(&mean);
        direction->Update(-0.5/abs(mean),*gradient,0.0);
      }
    }
  }
  return direction;
}

/*----------------------------------------------------------------------*/
ACOU::PATLineSearch::PATLineSearch(Teuchos::RCP<PatImageReconstruction> imagereconstruction)
{
  itermax_ = imagereconstruction->acouparams_->sublist("PA IMAGE RECONSTRUCTION").get<int>("INV_LS_MAX_RUN");
  alpha_max_ = 15.0;

  c1_ =  1.0e-12;
  c2_ = 0.9;

  imagereconstruction_ =  imagereconstruction;
  myrank_ = imagereconstruction_->myrank_;
}

/*----------------------------------------------------------------------*/
void ACOU::PATLineSearch::Init(double J0, Teuchos::RCP<Epetra_Vector> gradient, Teuchos::RCP<Epetra_Vector> direction, Teuchos::RCP<Epetra_Vector> state, const Epetra_Map* uniquemap)
{
  // initialize objective function value
  J_0_ = J0;
  J_i_ = J0;
  J_im1_ = J0;

  // initialize vectorial quantities
  dir_ = direction;
  step_ = Teuchos::rcp(new Epetra_Vector(dir_->Map()));
  state_ = Teuchos::rcp(new Epetra_Vector(*state));

  // initialize unique map
  uniquemap_ = Teuchos::rcp(new Epetra_Map(*uniquemap));

  // initialize norm of derivative of line search function
  imagereconstruction_->CalculateGradDirNorm(*dir_,*uniquemap_,&normgradphi_0_);
  normgradphi_i_ = 0.0;

  // set step lengths
  alpha_i_ = 1.0;
  alpha_im1_ = 0.0;
  alpha_x_ = 0.0;

  return;
}

/*----------------------------------------------------------------------*/
bool ACOU::PATLineSearch::Run()
{
  if(!myrank_)
    std::cout<<"*************** RUN LINE SEARCH: J_0 "<<J_0_<<" normgradphi_0 "<<normgradphi_0_<<std::endl;

  for(int i=0; i<itermax_; ++i)
  {
    if(!myrank_)
      std::cout<<"*************** line search iteration "<<i<<" of maximal "<<itermax_<<" line search iterations, alpha_i_ "<<alpha_i_<<", J_0_ "<<J_0_<<", J_i_ "<< J_i_<<", normgradphi_0_ "<<normgradphi_0_<<", normgradphi_i_ "<<normgradphi_i_<<std::endl;

    // update parameters
    step_->Update(1.0,*state_,0.0);
    step_->Update(alpha_i_,*dir_,1.0);
    imagereconstruction_->ReplaceParams(step_);

    // solve forward problem
    imagereconstruction_->SolveStandardScatra();
    imagereconstruction_->SolveStandardAcou();

    // evaluate objective function
    J_i_ = imagereconstruction_->EvalulateObjectiveFunction();

    // check first condition
    if( J_i_ > J_0_ + c1_ * alpha_i_ * normgradphi_0_ || (J_i_ >= J_im1_ && i>0) )
    {
      if(!myrank_)
        std::cout<<"*************** line search condition 1 met, J_i "<<J_i_<<", J_0 "<<J_0_<<std::endl;
      alpha_x_ = Zoom(alpha_im1_,alpha_i_,J_im1_,false);
      break;
    }
    else if(!myrank_)
      std::cout<<"*************** line search condition 1 NOT met, J_i "<<J_i_<<", J_0 "<<J_0_<<std::endl;

    // solve adjoint problem
    imagereconstruction_->SolveAdjointAcou();
    imagereconstruction_->SolveAdjointScatra();

    // calculate gradient
    imagereconstruction_->EvaluateGradient();
    imagereconstruction_->CalculateGradDirNorm(*dir_,*uniquemap_,&normgradphi_i_);

    // check second condition
    if(std::abs(normgradphi_i_) <= -c2_*normgradphi_0_)
    {
      alpha_x_ = alpha_i_;
      if(!myrank_)
        std::cout<<"*************** line search condition 2 met, |\\/phi|_i "<<normgradphi_i_<<", |\\/phi|_0 "<<normgradphi_0_<<std::endl;
      break;
    }
    else if(!myrank_)
      std::cout<<"*************** line search condition 2 NOT met, |\\/phi|_i "<<normgradphi_i_<<", |\\/phi|_0 "<<normgradphi_0_<<std::endl;

    // check third condition
    if(normgradphi_i_>= 0)
    {
      if(!myrank_)
        std::cout<<"*************** line search condition 3 met"<<std::endl;
      alpha_x_ = Zoom(alpha_i_,alpha_im1_,J_im1_,true);
      break;
    }
    else if(!myrank_)
      std::cout<<"*************** line search condition 3 not met"<<std::endl;

    // update alphas
    alpha_im1_ = alpha_i_;
    alpha_i_ *= 2.0; // = PredictStepLength();
  }

  if(alpha_x_ != 0.0)
  {
    if(!myrank_)
      std::cout<<"*************** line search succeeded"<<std::endl;
    return true;
  }
  else
  {
    if(!myrank_)
      std::cout<<"*************** line search failed"<<std::endl;

    // replace params by original params
    imagereconstruction_->ReplaceParams(state_);

    return false;
  }
}

/*----------------------------------------------------------------------*/
double ACOU::PATLineSearch::PredictStepLength()
{
  // step length must become longer when we are here!
  double alpha = 0.0;

  // calculate the coefficients of the quartic polynomial
  // double d = J_0_;
  double c = normgradphi_0_;
  double b = -1.0/alpha_i_/alpha_i_*(alpha_i_*normgradphi_i_+2.0*normgradphi_0_*alpha_i_-3.0*J_i_+3*J_0_);
  double a = 1.0/alpha_i_/alpha_i_/alpha_i_*(J_i_-J_0_-normgradphi_0_*alpha_i_-b*alpha_i_*alpha_i_);

  // calculate the minima of the quartic polynomial
  double radi = b*b/9.0/a/a-c/3.0/a;
  if(radi>0.0)
  {
    double alpha1 = -b/3.0/a+sqrt(radi);
    double alpha2 = -b/3.0/a-sqrt(radi);
    // check if the results suit us
    if(alpha1>alpha_i_&&alpha1<10.0*alpha_i_)
      alpha = alpha1;
    else if(alpha2>alpha_i_&&alpha2<10.0*alpha_i_)
      alpha = alpha2;
    else if(alpha1>10.0*alpha_i_&&alpha2>10.0*alpha_i_)
      alpha = 10.0*alpha_i_;
    else if(alpha1<alpha_i_&&alpha2<alpha_i_)
      alpha = 2.0*alpha_i_;
    else
      alpha = 2.0*alpha_i_;
  }
  else
  {
    // quadratic interpolation
    alpha = -normgradphi_0_/2.0/(J_i_-J_0_-normgradphi_0_*alpha_i_)*alpha_i_;
    if(alpha<alpha_i_)
      alpha = 2.0*alpha_i_;
    else if(alpha>10.0*alpha_i_)
      alpha = 10.0*alpha_i_;
  }
  return alpha;
}

/*----------------------------------------------------------------------*/
double ACOU::PATLineSearch::Zoom(double alpha_lo, double alpha_hi, double J_alpha_lo, bool turn)
{
  double alpha_j = 0.0;
  double J_j = 0.0;
  double normgradphi_j = 0.0;

  // zoom iterations
  for(int j=0; j<itermax_; ++j)
  {
    // update alpha
    if(turn)
      alpha_j = alpha_hi+(alpha_lo-alpha_hi)/3.0;
    else
      alpha_j = (alpha_lo + alpha_hi) / 2.0;

    // output for user
    if(!myrank_)
      std::cout<<"*************** zoom iteration "<<j<<": alpha_lo "<<alpha_lo<<" alpha_hi "<<alpha_hi<<" alpha_j "<<alpha_j<<std::endl;

    // update parameters
    step_->Update(1.0,*state_,0.0);
    step_->Update(alpha_j,*dir_,1.0);
    imagereconstruction_->ReplaceParams(step_);

    // solve forward problem
    imagereconstruction_->SolveStandardScatra();
    imagereconstruction_->SolveStandardAcou();

    // evaluate objective function
    J_j = imagereconstruction_->EvalulateObjectiveFunction();

    // output for user
    if(!myrank_)
      std::cout<<"J_j "<<J_j<<" J_0_ "<<J_0_<<" J_0_+c1_... "<<J_0_ + c1_ * alpha_j * normgradphi_0_<<" J_alpha_lo "<<J_alpha_lo<<std::endl;

    if( J_j > J_0_ + c1_ * alpha_j * normgradphi_0_ || J_j >= J_alpha_lo )
    {
      if(!myrank_)
        std::cout<<"*************** zoom condition 1 met"<<", J_j "<<J_j<<", J_0 "<<J_0_<<std::endl;
      alpha_hi = alpha_j;
    }
    else
    {
      // solve adjoint problem
      imagereconstruction_->SolveAdjointAcou();
      imagereconstruction_->SolveAdjointScatra();

      // calculate gradient
      imagereconstruction_->EvaluateGradient();
      imagereconstruction_->CalculateGradDirNorm(*dir_,*uniquemap_,&normgradphi_j);

      // check second condition
      if(std::abs(normgradphi_j) <= -c2_*normgradphi_0_)
      {
        if(!myrank_)
          std::cout<<"*************** zoom condition 2 met, |\\/phi|_j "<<normgradphi_j<<", |\\/phi|_0 "<<normgradphi_0_<<std::endl;
        return alpha_j;
      }

      // check third condition
      if( normgradphi_j * (alpha_hi - alpha_lo) >= 0 )
      {
        if(!myrank_)
          std::cout<<"*************** zoom condition 3 met"<<std::endl;
        alpha_hi = alpha_lo;
      }

      alpha_lo = alpha_j;
    }
  }

  return 0.0;
}


/*----------------------------------------------------------------------*/
ACOU::PATRegula::PATRegula(INPAR::ACOU::RegulaType regulatype, double tikhweight, double tvdweight, double tvdeps, Teuchos::RCP<DRT::Discretization> discret)
:
type_(regulatype),
tikhweight_(tikhweight),
tvdweight_(tvdweight),
tvdeps_(tvdeps),
tvdmatrix_(Teuchos::null),
discret_(discret)
{
  if(type_==INPAR::ACOU::pat_regula_tvd || type_==INPAR::ACOU::pat_regula_tikhtvd)
  {
    // create matrix
    tvdmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*discret_->ElementRowMap(),6,false));

    // fill matrix
    FillTvdMatrix();
  }
}

/*----------------------------------------------------------------------*/
void ACOU::PATRegula::Evaluate(Teuchos::RCP<Epetra_Vector> params, double* val)
{
  if(type_==INPAR::ACOU::pat_regula_tikh || type_==INPAR::ACOU::pat_regula_tikhtvd)
  {
    double norm2 = 0.0;
    params->Norm2(&norm2);
    *val += 0.5 * tikhweight_ * norm2 * norm2;
  }
  if(type_==INPAR::ACOU::pat_regula_tvd || type_==INPAR::ACOU::pat_regula_tikhtvd)
  {
    // export parameters to column map
    Epetra_Vector paramscol(tvdmatrix_->ColMap(),false);
    LINALG::Export(*params,paramscol);

    double functvalue=0.0;
    for (int i=0; i<params->MyLength(); i++)
    {
      // get weights of neighbouring parameters
      int lenindices = tvdmatrix_->NumMyEntries(i);
      int numindex;
      std::vector<int> indices(lenindices,0);
      std::vector<double> weights(lenindices,0);
      tvdmatrix_->ExtractMyRowCopy(i,lenindices,numindex,&weights[0],&indices[0]);

      double rowsum=0.0;
      double rowval=(*params)[i];
      for (int j=0; j<lenindices; j++)
      {
        if (indices[j]!=i) // skip substracting from itself
        {
          double colval=(paramscol)[indices[j]];
          rowsum += weights[j]*(colval-rowval)*(colval-rowval);
        }
      }

      // sum over all the rows
      functvalue+=sqrt(rowsum+tvdeps_);
    }

    // collect contributions from all procs
    double gfunctvalue=0.0;
    discret_->Comm().SumAll(&functvalue,&gfunctvalue,1);

    *val +=tvdweight_*gfunctvalue;
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PATRegula::EvaluateGradient(Teuchos::RCP<Epetra_Vector> params, Teuchos::RCP<Epetra_Vector> gradient)
{
  if(type_==INPAR::ACOU::pat_regula_tikh || type_==INPAR::ACOU::pat_regula_tikhtvd)
  {
    gradient->Update(1.0,*params,1.0);
  }
  if(type_==INPAR::ACOU::pat_regula_tvd || type_==INPAR::ACOU::pat_regula_tikhtvd)
  {
    // export parameters to column map
    Epetra_Vector paramscol(tvdmatrix_->ColMap(),false);
    LINALG::Export(*params,paramscol);

    Epetra_Vector gradientcol(tvdmatrix_->ColMap(),true);
    for (int i=0; i<params->MyLength(); i++)
    {
      // get weights of neighbouring parameters
      int lenindices = tvdmatrix_->NumMyEntries(i);
      int numindex;
      std::vector<int> indices(lenindices,0);
      std::vector<double> weights(lenindices,0);
      tvdmatrix_->ExtractMyRowCopy(i,lenindices,numindex,&weights[0],&indices[0]);

      // value of the parameter in this row
      double rowval=(*params)[i];

      // denominator
      double denom=0.0;
      for (int j=0; j<lenindices; j++)
      {
        if (indices[j]!=i) // skip substracting from itself
        {
          double colval=(paramscol)[indices[j]];
          denom += weights[j]*(colval-rowval)*(colval-rowval);
        }
      }
      denom = sqrt(denom+tvdeps_);

      // nominator i
      double nomi=0.0;
      for (int j=0; j<lenindices; j++)
      {
        if (indices[j]!=i) // skip substracting from itself
        {
          double colval=(paramscol)[indices[j]];
          nomi += weights[j]*(colval-rowval);
        }
      }
      nomi=nomi*(-1);

      // insert contributions into dof i of the gradient
      gradientcol.SumIntoMyValue(i,0,nomi/denom);

      // nominator j
      double nomj=0.0;
      for (int j=0; j<lenindices; j++)
      {
        if (indices[j]!=i) // skip substracting from itself
        {
          double colval=(paramscol)[indices[j]];
          nomj = weights[j]*(colval-rowval);
          // insert contributions into dof j of the gradient
          gradientcol.SumIntoMyValue(indices[j],0,nomj/denom);
        }
      }
    }

    // bring back gradient to the unique layout he came here with.
    // we have to add up off proc components so we cannot use
    // the LINALG::Export since it only provides CombineMode::Insert
    Epetra_MultiVector tmp(gradient->Map(), gradient->NumVectors(),true);
    Epetra_Export exporter(gradientcol.Map(), tmp.Map());
    int err = tmp.Export(gradientcol, exporter, Add);
    if (err)
      dserror("Export using exporter returned err=%d", err);

    gradient->Update(tvdweight_,tmp,1.0);
  }

  return;
}

/*----------------------------------------------------------------------*/
void ACOU::PATRegula::FillTvdMatrix()
{
  std::map< std::vector<int>, std::vector<int> > facemap; //map faces to corresponding elements/parameters
  std::map< std::vector<int>, double > faceweight;        //map of faces to its area

  // loop all elements
  for(int e=0; e<discret_->NumMyRowElements(); ++e)
  {
    // get the element
    DRT::Element* ele = discret_->lRowElement(e);
    int pgid = ele->Id();

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
    else
      dserror("creating internal faces for 1D elements (would be points) not implemented");

    // safety check
    if (nele != faces.size()) dserror("number of surfaces or lines does not match!");

    // get the surface/line elements for area computation
    std::vector<Teuchos::RCP<DRT::Element> >  surfs;
    if (ele->NumSurface() > 1)
      surfs = ele->Surfaces();
    else if (ele->NumSurface() == 1)
      surfs = ele->Lines();
    else
      dserror("creating internal faces for 1D elements (would be points) not implemented");

    // get nodes of each of this element's face
    for (unsigned int iele = 0; iele < nele; iele++)
    {
      // allocate node vectors
      unsigned int nnode = faces[iele].size();
      std::vector<int> nodeids(nnode);

      // get connectivity info
      for (unsigned int inode=0;inode<nnode;inode++)
        nodeids[inode] = ele->NodeIds()[faces[iele][inode]];

      //get the area of this face
      Teuchos::ParameterList p;
      //SetAction(p);
      if(discret_->Name()=="scatra")
        p.set<int>("action",SCATRA::bd_integrate_shape_functions);
      else if(discret_->Name()=="acou")
        p.set<int>("action", ACOU::bd_integrate);
      else
        dserror("what happenend?");

      p.set("area",0.0);
      DRT::Element::LocationArray la(discret_->NumDofSets());
      surfs[iele]->LocationVector(*discret_,la,false);
      //initialize element vectors
      int ndof = ele->NumNode()*3;
      Epetra_SerialDenseMatrix elematrix1(ndof,ndof,false);
      Epetra_SerialDenseMatrix elematrix2(ndof,ndof,false);
      Epetra_SerialDenseVector elevector1(ndof);
      Epetra_SerialDenseVector elevector2(ndof);
      Epetra_SerialDenseVector elevector3(ndof);
      surfs[iele]->Evaluate(p,*discret_,la,elematrix1,elematrix2,elevector1,elevector2,elevector3);
      double area=p.get("area",-1.0);
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

  const int numprocs = discret_->Comm().NumProc();
  int allproc[numprocs];
  for (int i = 0; i < numprocs; ++i)
    allproc[i] = i;

  for( int i=0; i<discret_->Comm().NumProc(); i++ )
  {
    // by looping the procs one by one we need to make sure that every proc
    // does the same number of loops regardless of its own lenght of keys
    // in the facemap. so the actual number of loops is defined by the current
    // proc and distributed in localmapsize
    int localmapsize;
    if( discret_->Comm().MyPID()==i ) localmapsize = facemap.size();
    discret_->Comm().Broadcast(&localmapsize,1,i);

    // now iterate as often as there are faces on proc i
    std::map< std::vector<int>,std::vector<int> >::iterator face_it(facemap.begin());
    for( int j=0; j<localmapsize; j++ )
    {
      //get length of current face-key on all procs
      int keylength=0;
      if( discret_->Comm().MyPID()==i ) keylength=face_it->first.size();
      discret_->Comm().Broadcast(&(keylength),1,i);

      //send current face-key to all procs
      std::vector<int> facekey(keylength,0);
      if( discret_->Comm().MyPID()==i ) facekey=face_it->first;
      discret_->Comm().Broadcast(&(facekey[0]),keylength,i);

      //check whether one of the other procs also has this key and write IDs of
      // procs who own this face in "sowningprocs" and distribute this knowledge
      // among all procs to rowningprocs
      std::map< std::vector<int>,std::vector<int> >::iterator face_abroad(facemap.find(facekey));
      std::vector<int> sowningprocs;
      std::vector<int> rowningprocs;
      if( face_abroad!=facemap.end() && discret_->Comm().MyPID() != i)
        sowningprocs.push_back(discret_->Comm().MyPID());
      LINALG::Gather(sowningprocs,rowningprocs,numprocs,allproc,discret_->Comm());

      // now bring parameters corresponding to this face on the other procs to proc i
      // (they are send to all procs but only proc i stores them in the map with
      // the current face-key)
      std::vector<int> sparams;
      std::vector<int> rparams;
      if( std::find(rowningprocs.begin(),rowningprocs.end(),discret_->Comm().MyPID()) != rowningprocs.end())
        sparams=facemap[facekey];
      LINALG::Gather(sparams,rparams,numprocs,allproc,discret_->Comm());

      // store additional elements on proc i
      if( discret_->Comm().MyPID() == i )
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
      if( discret_->Comm().MyPID()==i ) face_it++;

    } // faces on each proc
  } // procs

  /*------------------------------------------------------------------- */
  // STEP 3: sort elements into graph according to neighbor information
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
      if (discret_->ElementRowMap()->LID(globalrow)>=0)
      {
        // like this the diagonal entries are inserted redundantly and summed up
        // after FillComplete() is called; they are more or less useless anyways
        tvdmatrix_->InsertGlobalValues(globalrow,parameters.size(),&weights[0],&parameters[0]);
      }
    }
  }

  // finalize the graph ...
  tvdmatrix_->FillComplete();
  tvdmatrix_->OptimizeStorage();

  // delete diagonal values
  Teuchos::RCP<Epetra_Vector> diagonal=Teuchos::rcp(new Epetra_Vector(*discret_->ElementRowMap(), true));
  tvdmatrix_->ReplaceDiagonalValues(*diagonal);

  return;
}

