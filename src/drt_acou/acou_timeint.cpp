/*!----------------------------------------------------------------------
\file acou_timeint.cpp
\brief Base class functions for implicit and explicit time integration

<pre>
\level 2

\maintainer Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include "acou_timeint.H"
#include "acou_ele_action.H"
#include "acou_ele.H"
#include "acou_resulttest.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_mat/scatra_mat.H"
#include "../linalg/linalg_utils.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 |  Constructor (public)                                 schoeder 03/15 |
 *----------------------------------------------------------------------*/
ACOU::AcouTimeInt::AcouTimeInt(
  const Teuchos::RCP<DRT::DiscretizationHDG>&   actdis,
  const Teuchos::RCP<LINALG::Solver>&           solver,
  const Teuchos::RCP<Teuchos::ParameterList>&   params,
  const Teuchos::RCP<IO::DiscretizationWriter>& output
  ):
  discret_        (actdis),
  solver_         (solver),
  params_         (params),
  output_         (output),
  dyna_           (DRT::INPUT::IntegralValue<INPAR::ACOU::DynamicType>(*params_,"TIMEINT")),
  phys_           (DRT::INPUT::IntegralValue<INPAR::ACOU::PhysicalType>(*params,"PHYSICAL_TYPE")),
  invana_         (params_->get<bool>("invana")),
  padaptivity_    (DRT::INPUT::IntegralValue<bool>(*params_,"P_ADAPTIVITY")),
  adjoint_        (params_->get<bool>("adjoint")),
  myrank_         (actdis->Comm().MyPID()),
  writemonitor_   (DRT::INPUT::IntegralValue<bool>(*params_,"WRITEMONITOR")),
  writestress_    (DRT::INPUT::IntegralValue<bool>(*params_,"WRITESTRESS")),
  time_           (0.0),
  step_           (0),
  restart_        (params_->get<int>("restart")),
  maxtime_        (params_->get<double>("MAXTIME")),
  stepmax_        (params_->get<int>("NUMSTEP")),
  uprestart_      (params_->get<int>("RESTARTEVRY", -1)),
  upres_          (params_->get<int>("RESULTSEVRY", -1)),
  numdim_         (DRT::Problem::Instance()->NDim()),
  dtp_            (params_->get<double>("TIMESTEP")),
  adjoint_rhs_    (Teuchos::null)
{
  // create the global trace vectors
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  velnp_ = LINALG::CreateVector(*dofrowmap,true);

  if(params_->isParameter("rhsvec"))
    adjoint_rhs_ = params_->get<Teuchos::RCP<Epetra_MultiVector> >("rhsvec");

  if(!invana_)
    params_->set<bool>("timereversal",false);

} // AcouTimeInt

/*----------------------------------------------------------------------*
 |  Desctructor (public)                                 schoeder 03/15 |
 *----------------------------------------------------------------------*/
ACOU::AcouTimeInt::~AcouTimeInt()
{}

/*----------------------------------------------------------------------*
 |  Print some information (public)                      schoeder 03/15 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::PrintInformationToScreen()
{
  if (!myrank_)
  {
    std::cout << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    if(!adjoint_)
      std::cout << "INTEGRATION OF AN ACOUSTIC PROBLEM USING HDG" << std::endl;
    else
      std::cout << "INTEGRATION OF AN ADJOINT ACOUSTIC PROBLEM USING HDG" << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << "number of nodes             " << discret_->NumGlobalNodes() << std::endl;
    std::cout << "number of elements          " << discret_->NumGlobalElements() << std::endl;
    std::cout << "number of faces             " << discret_->NumGlobalFaces() << std::endl;
    std::cout << "number of trace unknowns    " << discret_->DofRowMap(0)->NumGlobalElements() << std::endl;
    std::cout << "number of interior unknowns " << discret_->DofRowMap(1)->NumGlobalElements() << std::endl;
    std::cout << "time step size              " << dtp_ << std::endl;
    std::cout << "time integration scheme     " << Name() << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
  }
  return;
} // PrintInformationToScreen

/*----------------------------------------------------------------------*
 | Initialization by given scatra solution vector (pub)  schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::SetInitialZeroField()
{
  // we have to call an init for the elements (for inverse analysis, otherwise adjoint run starts from values of previous forward run)
  Teuchos::ParameterList initParams;
  initParams.set<int>("action",ACOU::ele_init);
  initParams.set<bool>("padaptivity",false);
  initParams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  initParams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  // discret_->Evaluate(initParams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  Epetra_SerialDenseVector elevec;
  Epetra_SerialDenseMatrix elemat;
  DRT::Element::LocationArray la(2);
  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    DRT::Element *ele = discret_->lColElement(el);
    ele->Evaluate(initParams,*discret_,la[0].lm_,elemat,elemat,elevec,elevec,elevec);
  }

  return;
}


/*----------------------------------------------------------------------*
 | Initialization by given scatra solution vector (pub)  schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::SetInitialPhotoAcousticField(Teuchos::RCP<Epetra_Vector> light,
                                                     Teuchos::RCP<DRT::Discretization> scatradis,
                                                     bool meshconform)
{
  // we have to call an init for the elements first!
  Teuchos::ParameterList initParams;
  initParams.set<int>("action",ACOU::ele_init);
  initParams.set<bool>("padaptivity",false);
  initParams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  initParams.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  // discret_->Evaluate(initParams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  Epetra_SerialDenseVector elevec;
  Epetra_SerialDenseMatrix elemat;
  DRT::Element::LocationArray la(2);
  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    DRT::Element *ele = discret_->lColElement(el);
    ele->Evaluate(initParams,*discret_,la[0].lm_,elemat,elemat,elevec,elevec,elevec);
  }

  // export light vector to column map, this is necessary
  Teuchos::RCP<Epetra_Vector> lightcol = Teuchos::rcp(new Epetra_Vector(*(scatradis->DofColMap())));
  LINALG::Export(*light,*lightcol);

  Epetra_SerialDenseVector elevec1, elevec2, elevec3;
  Epetra_SerialDenseMatrix elemat1, elemat2;

  initParams.set<int>("action",ACOU::project_optical_field);
  initParams.set<bool>("mesh conform",meshconform);

  int numoptele = scatradis->NumGlobalElements();

  if(meshconform)
  {
    int minoptelegid = scatradis->ElementRowMap()->MinAllGID();
    int minacouelegid = discret_->ElementRowMap()->MinAllGID();

    std::vector<int> localrelevantghostelements;
    // loop all optical elements (usually there are less optical elements than acoustical elements
    for(int optel=0; optel<numoptele; ++optel)
    {
      int myopteleowner = -1;
      int opteleowner = -1;
      int myacoueleowner = -1;
      int acoueleowner = -1;

      // determine owner of the optical element
      DRT::Element* optele = NULL;
      if ( scatradis->HaveGlobalElement(optel+minoptelegid) )
      {
        optele = scatradis->gElement(optel+minoptelegid);
        myopteleowner = optele->Owner();
        if ( myopteleowner != myrank_ ) myopteleowner = -1; // do not want to consider ghosted elements
      }
      // now, every proc knows the owner of this optical element
      scatradis->Comm().MaxAll(&myopteleowner,&opteleowner,1);

      // who is the owner of the associated acoustical element?
      DRT::Element* acouele = NULL;
      if( discret_->HaveGlobalElement(optel+minacouelegid) )
      {
        acouele = discret_->gElement(optel+minacouelegid);
        myacoueleowner = acouele->Owner();
        if ( myacoueleowner != myrank_ )
        {
          myacoueleowner = -1;
          localrelevantghostelements.push_back(optel+minacouelegid);
        } // do not want to consider ghosted elements
      }

      // now, every proc knows the owner of this acoustical element
      discret_->Comm().MaxAll(&myacoueleowner,&acoueleowner,1);

      // easiest case: both elements are on the same processor
      if ( acoueleowner == opteleowner )
      {
        // the owning proc can do all his business
        if ( opteleowner == myrank_ )
        {
          elevec1.Scale(0.0);elevec2.Scale(0.0);
          acouele->LocationVector(*discret_,la,false);
          if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
            elevec1.Shape(la[0].lm_.size(), 1);
          if (elevec2.M() != discret_->NumDof(1,acouele))
            elevec2.Shape(discret_->NumDof(1,acouele), 1);

          // we need the absorption coefficient of the light element
          double absorptioncoeff = static_cast <MAT::ScatraMat*>((optele->Material()).get())->ReaCoeff(scatradis->ElementColMap()->LID(optele->Id()));

          initParams.set<double>("absorption",absorptioncoeff);
          Teuchos::RCP<std::vector<double> > nodevals = Teuchos::rcp(new std::vector<double>);
          int numlightnode = optele->NumNode();
          (*nodevals).resize((numdim_+1)*numlightnode);

          // fill nodevals with node coords and nodebased solution values
          DRT::Node** lightnodes = optele->Nodes();
          for(int i=0; i<numlightnode; ++i)
          {
            for(int j=0; j<numdim_; ++j)
            {
              (*nodevals)[i*(numdim_+1)+j] = lightnodes[i]->X()[j];
            }

            int dof = scatradis->Dof(0,lightnodes[i],0);
            int lid = lightcol->Map().LID(dof);
            if ( lid < 0 )
              dserror("given dof is not stored on proc %d although map is colmap",myrank_);
            else
              (*nodevals)[i*(numdim_+1)+numdim_] = (*(lightcol.get()))[lid];
          }

          initParams.set<Teuchos::RCP<std::vector<double> > >("nodevals",nodevals);

          // evaluate the element
          acouele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec1,elevec2,elevec3);

//          // fill evelvec1 into the global vector
//          for (unsigned int i=0; i<la[0].lm_.size(); ++i)
//          {
//            const int lid = dofrowmap->LID(la[0].lm_[i]);
//            if ( lid >= 0 )
//              (*velnp_)[lid] = elevec1(i);
//          }

        }
        discret_->Comm().Barrier(); // other procs please wait for the one, who did all the work
      } // if ( acoueleowner == opteleowner )
      else // the other case: the acoustical element and the optical element are owned by different procs
      {
        // now, we perform several steps:
        // 1. get the information from the optical element we need
        // 2. communicate this information to the proc owning the acoustical element
        // 3. the acoustical element owning proc shall do his business

        // this is the 1. step:
        Teuchos::RCP<std::vector<double> > nodevals = Teuchos::rcp(new std::vector<double>);
        int size = 0;
        double absorptioncoeff = 0.0;
        if ( myrank_ == opteleowner )
        {
          // we need the absorption coefficient of the light element
          absorptioncoeff = static_cast <MAT::ScatraMat*>((optele->Material()).get())->ReaCoeff(scatradis->ElementColMap()->LID(optele->Id()));

          int numlightnode = optele->NumNode();
          size = (numdim_+1)*numlightnode;
          (*nodevals).resize(size);

          DRT::Node** lightnodes = optele->Nodes();
          for(int i=0; i<numlightnode; ++i)
          {
            for(int j=0; j<numdim_; ++j)
            {
              (*nodevals)[i*(numdim_+1)+j] = lightnodes[i]->X()[j];
            }

            int dof = scatradis->Dof(0,lightnodes[i],0);
            int lid = lightcol->Map().LID(dof);
            if ( lid < 0 )
              dserror("given dof is not stored on proc %d although map is colmap",myrank_);
            else
              (*nodevals)[i*(numdim_+1)+numdim_] = (*(lightcol.get()))[lid];
          }
        }

        // this is the 2. step:
        discret_->Comm().Broadcast(&size,1,opteleowner);
        discret_->Comm().Broadcast(&absorptioncoeff,1,opteleowner);
        if ( myrank_ != opteleowner)
          (*nodevals).resize(size);
        discret_->Comm().Broadcast(&(*nodevals)[0],size,opteleowner);

        // this is the 3. step:
        if ( myrank_ == acoueleowner )
        {
          elevec1.Scale(0.0);elevec2.Scale(0.0);
          acouele->LocationVector(*discret_,la,false);
          if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
            elevec1.Shape(la[0].lm_.size(), 1);
          if (elevec2.M() != discret_->NumDof(1,acouele))
            elevec2.Shape(discret_->NumDof(1,acouele), 1);

          initParams.set<double>("absorption",absorptioncoeff);
          initParams.set<Teuchos::RCP<std::vector<double> > >("nodevals",nodevals);

          acouele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec1,elevec2,elevec3);

          // fill evelvec1 into the global vector
          // for (unsigned int i=0; i<la[0].lm_.size(); ++i)
          // {
          //   const int lid = dofrowmap->LID(la[0].lm_[i]);
          //   if ( lid >= 0 )
          //     (*velnp_)[lid] = elevec1(i);
          // }
        }

      } // else ** if ( acoueleowner == opteleowner )
    } // for(int optel=0; optel<numoptele; ++optel)

    int sumghostele = -1;
    int localnumghostele = localrelevantghostelements.size();
    discret_->Comm().SumAll(&localnumghostele,&sumghostele,1);
    std::vector<int> relevantghostelements;

    for(int proc=0; proc<discret_->Comm().NumProc(); ++proc)
    {
      int vals = -1;
      int locsize = localnumghostele;
      discret_->Comm().Broadcast(&locsize,1,proc);
      for(int j=0; j<locsize; ++j)
      {
        if(myrank_==proc) vals=localrelevantghostelements[j];
        discret_->Comm().Broadcast(&vals,1,proc);
        relevantghostelements.push_back(vals);
      }
    }

    // communicate to ghosted elements
    for(unsigned int i=0; i<relevantghostelements.size(); ++i)
    {
      DRT::Element* ele = NULL;
      int owner = 0;
      std::vector<double> pressvals;
      if( discret_->HaveGlobalElement(relevantghostelements[i]) ) // true for owner and the proc with ghost
      {
        ele = discret_->gElement(relevantghostelements[i]);
        owner = ele->Owner();
        Epetra_SerialDenseVector tempvec;
        if(owner == myrank_)
        {
          DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(ele);
          pressvals.resize(acouele->eleinteriorPressnp_.M());
          for(int j=0; j<acouele->eleinteriorPressnp_.M(); ++j)
            pressvals[j] = acouele->eleinteriorPressnp_(j);
        }
      }
      int actualowner = 0;

      discret_->Comm().MaxAll(&owner,&actualowner,1);
      int size = pressvals.size();
      discret_->Comm().Broadcast(&size,1,actualowner);
      if(myrank_!=actualowner)
        pressvals.resize(size);
      discret_->Comm().Broadcast(&pressvals[0],size,actualowner);
      // each proc has the values, now the proc with the ghost element has to update the values
      if( discret_->HaveGlobalElement(relevantghostelements[i]) && actualowner!=myrank_ )
      {
        DRT::Element* ele = discret_->gElement(relevantghostelements[i]);
        DRT::ELEMENTS::Acou * acouele = dynamic_cast<DRT::ELEMENTS::Acou*>(ele);
        Epetra_SerialDenseVector tempvec(size);
        for(int j=0; j<size; ++j) tempvec(j) = pressvals[j];
        acouele->eleinteriorPressnp_=tempvec;
      }
    }

  } // if(meshconform)
  else
  {
    // TODO: have to set the initial field for column elements as well!
    if(!myrank_) std::cout<<"Welcome to nonconform mapping, this might take a while! Please have a little patience."<<std::endl;
    double tcpumap=Teuchos::Time::wallTime();
    // first of all, we want to get a nodebased vector of pressure values
    // in a second step, we evaluate the internal pressure field
    // this is necessary, otherwise we get problems with the interpolation, especially for acoustical
    // elements which are bigger than optical elements!

    /*************************** STEP 1 - COMPUTE NODE BASED PRESSURE FIELD ***************************/
    Teuchos::RCP<Epetra_Vector> pressurenode = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeRowMap())));

    // loop all acoustical nodes
    int minacounodegid = discret_->NodeRowMap()->MinAllGID();
    int percent_counter = 0;
    for(int acound = 0; acound<discret_->NumGlobalNodes(); ++acound)
    {
      if( (int)((acound*100)/discret_->NumGlobalNodes()) > (int)(10*percent_counter))
      {
        if(!myrank_)
        {
          std::cout << "---------------------------------------" << std::endl;
          std::cout << (int)((acound*100)/discret_->NumGlobalNodes()) -1 << "% of Coupling Evaluations are done!" << std::endl;
        }
        percent_counter++;
      }

      // get node
      DRT::Node* acounode = NULL;
      int myacounodeowner = -1;
      if(discret_->HaveGlobalNode(acound+minacounodegid))
      {
        acounode = discret_->gNode(acound+minacounodegid);
        myacounodeowner = acounode->Owner();
        if(myacounodeowner != myrank_) myacounodeowner = -1;
      }
      int acounodeowner = -1;
      discret_->Comm().MaxAll(&myacounodeowner,&acounodeowner,1);

      double acounodecoords[numdim_];
      if(myrank_==acounodeowner)
      {
        for(int d=0; d<numdim_; ++d)
          acounodecoords[d] = acounode->X()[d];
      }
      discret_->Comm().Broadcast(&acounodecoords[0],numdim_,acounodeowner);

      double p = 0.0;
      // now find the corresponding optiele!
      for(int optel=0; optel<scatradis->NumMyRowElements(); ++optel)
      {

        DRT::Element* ele = scatradis->lRowElement(optel);
        // get the nodes of this element, and then check if acoustical node is inside
        if(ele->Shape()==DRT::Element::quad4)
        {
          double optnodecoords[4][numdim_];
          double minmaxvals[2][numdim_];
          for(int j=0; j<numdim_; ++j)
          {
            minmaxvals[0][j] = 1.0e6; // minvals
            minmaxvals[1][j] = -1.0e6; // maxvals
          }
          for(int nd=0;nd<4;++nd) // quad4 has 4 nodes
            for(int d=0;d<numdim_;++d)
            {
              optnodecoords[nd][d] = ele->Nodes()[nd]->X()[d];
              if(optnodecoords[nd][d] < minmaxvals[0][d]) minmaxvals[0][d]=optnodecoords[nd][d];
              if(optnodecoords[nd][d] > minmaxvals[1][d]) minmaxvals[1][d]=optnodecoords[nd][d];
            }
          // check, if acoustical node is in bounding box
          bool inside = true;
          for(int d=0;d<numdim_;++d)
            if(acounodecoords[d]>minmaxvals[1][d]+1.0e-5 || acounodecoords[d]<minmaxvals[0][d]-1.0e-5)
              inside=false;

          if(inside)
          {
            // solve for xi by local Newton
            LINALG::Matrix<2,1> F(true);
            LINALG::Matrix<2,2> dFdxi(true);
            LINALG::Matrix<2,1> xi(true);
            LINALG::Matrix<2,1> deltaxi(true);
            double deltaxinorm = 0.0;
            int count = 0;
            do{
              count++;
              F(0) = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * optnodecoords[0][0]
                   + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * optnodecoords[1][0]
                   + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * optnodecoords[2][0]
                   + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * optnodecoords[3][0]  - acounodecoords[0];
              F(1) = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * optnodecoords[0][1]
                   + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * optnodecoords[1][1]
                   + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * optnodecoords[2][1]
                   + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * optnodecoords[3][1]  - acounodecoords[1] ;

              dFdxi(0,0) = - 0.25 * (1. - xi(1)) * optnodecoords[0][0]
                           + 0.25 * (1. - xi(1)) * optnodecoords[1][0]
                           + 0.25 * (1. + xi(1)) * optnodecoords[2][0]
                           - 0.25 * (1. + xi(1)) * optnodecoords[3][0] ;
              dFdxi(0,1) = - 0.25 * (1. - xi(0)) * optnodecoords[0][0]
                           - 0.25 * (1. + xi(0)) * optnodecoords[1][0]
                           + 0.25 * (1. + xi(0)) * optnodecoords[2][0]
                           + 0.25 * (1. - xi(0)) * optnodecoords[3][0] ;
              dFdxi(1,0) = - 0.25 * (1. - xi(1)) * optnodecoords[0][1]
                           + 0.25 * (1. - xi(1)) * optnodecoords[1][1]
                           + 0.25 * (1. + xi(1)) * optnodecoords[2][1]
                           - 0.25 * (1. + xi(1)) * optnodecoords[3][1] ;
              dFdxi(1,1) = - 0.25 * (1. - xi(1)) * optnodecoords[0][1]
                           - 0.25 * (1. + xi(1)) * optnodecoords[1][1]
                           + 0.25 * (1. + xi(1)) * optnodecoords[2][1]
                           + 0.25 * (1. - xi(1)) * optnodecoords[3][1] ;

              LINALG::FixedSizeSerialDenseSolver<2,2,1> inverser;
              inverser.SetMatrix(dFdxi);
              inverser.SetVectors(deltaxi,F);
              inverser.Solve();

              deltaxinorm = deltaxi.Norm2();
              xi.Update(-1.0,deltaxi,1.0);
            } while ( deltaxinorm > 1.0e-8 && count < 10 );

            if(!(count == 10 || xi.NormInf()>1.0+0.1))
            {
              double absorptioncoeff = static_cast <MAT::ScatraMat*>((ele->Material()).get())->ReaCoeff(scatradis->ElementColMap()->LID(ele->Id()));
              // get the values!
              double values[4] = {0};
              for(int nd=0;nd<4;++nd) // quad4 has 4 nodes
              {
                int dof = scatradis->Dof(0,ele->Nodes()[nd],0);
                int lid = lightcol->Map().LID(dof);
                if ( lid < 0 )
                  dserror("given dof is not stored on proc %d although map is colmap",myrank_);
                else
                  values[nd] = (*(lightcol.get()))[lid];
              }

              p = 0.25 * ( (1. - xi(0))*(1. - xi(1)) ) * values[0]
                + 0.25 * ( (1. + xi(0))*(1. - xi(1)) ) * values[1]
                + 0.25 * ( (1. + xi(0))*(1. + xi(1)) ) * values[2]
                + 0.25 * ( (1. - xi(0))*(1. + xi(1)) ) * values[3];
              p *= -absorptioncoeff;
            }
          } // if(inside)
        }
        else dserror("up to now only implemented for quad4");
      } // for(int optel=0; optel<scatradis->NumMyRowElements(); ++optel)
      // one processor might provide a value

      double glob_p_min = 0.0;
      discret_->Comm().MinAll(&p,&glob_p_min,1);
      double glob_p_max = 0.0;
      discret_->Comm().MaxAll(&p,&glob_p_max,1);
      // take higher absolut values
      double glob_p = 0.0;
      if(std::abs(glob_p_min)>std::abs(glob_p_max))
        glob_p = glob_p_min;
      else
        glob_p = glob_p_max;

      // set p value in node based vector
      if(myrank_==acounodeowner && glob_p != 0.0)
        pressurenode->ReplaceGlobalValue(acounode->Id(),0,glob_p);
    } // for(int acound = 0; acound<discret_->NumGlobalNodes(); ++acound)

    /*************************** STEP 2 - NODE BASED -> DOF BASED FIELD ***************************/

    Teuchos::RCP<Epetra_Vector> pressurenodecol = Teuchos::rcp(new Epetra_Vector(*(discret_->NodeColMap())));
    LINALG::Export(*pressurenode,*pressurenodecol);
    initParams.set<Teuchos::RCP<Epetra_Vector> >("pressurenode",pressurenodecol);

    for(int acouel=0; acouel<discret_->NumMyRowElements(); ++acouel)
    {
      DRT::Element* acouele = discret_->lRowElement(acouel);

      elevec1.Scale(0.0);elevec2.Scale(0.0);
      acouele->LocationVector(*discret_,la,false);
      if (static_cast<std::size_t>(elevec1.M()) != la[0].lm_.size())
        elevec1.Shape(la[0].lm_.size(), 1);
      if (elevec2.M() != discret_->NumDof(1,acouele))
        elevec2.Shape(discret_->NumDof(1,acouele), 1);

      acouele->LocationVector(*discret_,la,false);
      acouele->Evaluate(initParams,*discret_,la[0].lm_,elemat1,elemat2,elevec1,elevec2,elevec3);
    }
    if(!myrank_)
    {
      std::cout <<"---------------------------------------" << std::endl;
      std::cout <<"100% of Coupling Evaluations are done! It took "<<Teuchos::Time::wallTime()-tcpumap<<" seconds!"<<std::endl;
      std::cout <<"---------------------------------------" << std::endl;
    }
  } // else ** if(meshconform)

  return;
} // SetInitialPhotoAcousticField

/*----------------------------------------------------------------------*
 |  Write restart vectors (public)                       schoeder 06/15 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::WriteRestart()
{
  if (myrank_ == 0 && !invana_)
    std::cout<<"======= Restart written in step "<<step_<<std::endl;

  output_->WriteVector("velnps", velnp_);

  // write internal field for which we need to create and fill the corresponding vectors
  // since this requires some effort, the WriteRestart method should not be used excessively!
  Teuchos::RCP<Epetra_Vector> intvelnp = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",ACOU::fill_restart_vecs);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("padaptivity",padaptivity_);
  eleparams.set<int>("useacouoptvecs",-6);

  discret_->SetState(1,"intvelnp",intvelnp);
  discret_->Evaluate(eleparams);

  Teuchos::RCP<const Epetra_Vector> matrix_state = discret_->GetState(1,"intvelnp");
  output_->WriteVector("intvelnp",matrix_state);
  discret_->ClearState(true);
  return;
} // WriteRestart

/*----------------------------------------------------------------------*
 |  ReadRestart (public)                                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  Teuchos::RCP<Epetra_Vector> intvelnp = Teuchos::rcp(new Epetra_Vector(*(discret_->DofRowMap(1))));
  reader.ReadVector(intvelnp,"intvelnp");
  reader.ReadVector(velnp_,"velnps");

  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",ACOU::ele_init_from_restart);
  eleparams.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  eleparams.set<bool>("padaptivity",padaptivity_);
  discret_->SetState(1,"intvelnp",intvelnp);
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  discret_->ClearState(true);

  return;
} // ReadRestart

/*----------------------------------------------------------------------*
 |  OutputDensityAndSpeedOfSound()                       schoeder 06/15 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::OutputDensityAndSpeedOfSound()
{
  // build the two vectors
  Teuchos::RCP<Epetra_Vector> densvec = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap()),false));
  Teuchos::RCP<Epetra_Vector> cvec = Teuchos::rcp(new Epetra_Vector(*(discret_->ElementRowMap()),false));

  for (int i=0; i<discret_->NumMyRowElements(); i++)
  {
    DRT::Element* actele;
    actele = discret_->lRowElement(i);
    double dens = actele->Material()->Parameter()->GetParameter(0,discret_->ElementColMap()->LID(actele->Id()));
    double c = actele->Material()->Parameter()->GetParameter(1,discret_->ElementColMap()->LID(actele->Id()));
    densvec->operator [](i) = dens;
    cvec->operator [](i) = c;
  }
  output_->WriteVector("density",densvec);
  output_->WriteVector("speedofsound",cvec);

  return;
}

namespace
{
  void getNodeVectorsHDG (DRT::Discretization               &dis,
                          const Teuchos::RCP<Epetra_Vector> &traceValues,
                          const int                          ndim,
                          Teuchos::RCP<Epetra_MultiVector>  &velocity,
                          Teuchos::RCP<Epetra_Vector>       &pressure,
                          Teuchos::RCP<Epetra_Vector>       &tracevel,
                          Teuchos::RCP<Epetra_Vector>       &cellPres,
                          INPAR::ACOU::PhysicalType         phys,
                          bool                              padapt)
  {
    //if (pressure.get() == NULL || pressure->MyLength() != dis.NumMyRowNodes())
    {
      const Epetra_Map* nodemap = dis.NodeRowMap();
      pressure.reset(new Epetra_Vector(*nodemap));
      velocity.reset(new Epetra_MultiVector(*nodemap,3));
      tracevel.reset(new Epetra_Vector(pressure->Map()));
      cellPres.reset(new Epetra_Vector(*dis.ElementRowMap()));
    }

    // call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::interpolate_hdg_to_node);
    params.set<int>("useacouoptvecs",-1);
    dis.SetState(0,"trace",traceValues);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys);
    params.set<bool>("padaptivity",padapt);

    DRT::Element::LocationArray la(2);

    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(pressure->MyLength());
    velocity->PutScalar(0.);
    pressure->PutScalar(0.);

    for (int el=0; el<dis.NumMyColElements();++el)
    {
      DRT::Element *ele = dis.lColElement(el);
      ele->LocationVector(dis,la,false);
      if (interpolVec.M() == 0)
        interpolVec.Resize(ele->NumNode()*(ndim+2)+1);

      ele->Evaluate(params,dis,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i=0; i<ele->NumNode(); ++i)
      {
        DRT::Node* node = ele->Nodes()[i];
        const int localIndex = pressure->Map().LID(node->Id());

        if (localIndex < 0)
          continue;

        touchCount[localIndex]++;
        for (int d=0; d<ndim; ++d)
        {
          velocity->SumIntoMyValue(localIndex,d,interpolVec(i+d*ele->NumNode()));
        }
        (*pressure)[localIndex] += interpolVec(i+ndim*ele->NumNode());
        (*tracevel)[localIndex] += interpolVec(i+(ndim+1)*ele->NumNode());
      }

      const int eleIndex = dis.ElementRowMap()->LID(ele->Id());
      if (eleIndex >= 0)
        (*cellPres)[eleIndex] += interpolVec((ndim+2)*ele->NumNode());
    }

    for (int i=0; i<pressure->MyLength(); ++i)
    {
      (*pressure)[i] /= touchCount[i];
      for (int d=0; d<ndim; ++d)
        (*velocity)[d][i] /= touchCount[i];
      (*tracevel)[i] /= touchCount[i];
    }
    dis.ClearState(true);

    return;
  } // getNodeVectorsHDG

  void getNodeVectorsHDGSolid(DRT::Discretization              &dis,
                            const Teuchos::RCP<Epetra_Vector> &traceValues,
                            const int                          ndim,
                            Teuchos::RCP<Epetra_MultiVector>  &velocitygradient,
                            Teuchos::RCP<Epetra_MultiVector>  &velocity,
                            Teuchos::RCP<Epetra_Vector>       &pressure,
                            Teuchos::RCP<Epetra_MultiVector>  &tracevelocity,
                            Teuchos::RCP<Epetra_Vector>       &cellPres,
                            INPAR::ACOU::PhysicalType         phys,
                            bool                              writestress)
  {
    {
      const Epetra_Map* nodemap = dis.NodeRowMap();
      velocity.reset(new Epetra_MultiVector(*nodemap,3));
      velocitygradient.reset(new Epetra_MultiVector(*nodemap,ndim*ndim));
      pressure.reset(new Epetra_Vector(*nodemap));
      tracevelocity.reset(new Epetra_MultiVector(*nodemap,3));
      cellPres.reset(new Epetra_Vector(*dis.ElementRowMap()));
    }

    // call element routine for interpolate HDG to elements
    Teuchos::ParameterList params;
    params.set<int>("action",ACOU::interpolate_hdg_to_node);
    params.set<INPAR::ACOU::PhysicalType>("physical type",phys);
    params.set<bool>("writestress",writestress);
    dis.SetState(0,"trace",traceValues);

    std::vector<int> dummy;
    DRT::Element::LocationArray la(2);

    Epetra_SerialDenseMatrix dummyMat;
    Epetra_SerialDenseVector dummyVec;
    Epetra_SerialDenseVector interpolVec;
    std::vector<unsigned char> touchCount(pressure->MyLength());

    velocity->PutScalar(0.0);
    pressure->PutScalar(0.0);
    tracevelocity->PutScalar(0.0);
    cellPres->PutScalar(0.0);

    for (int el=0; el<dis.NumMyColElements();++el)
    {
      DRT::Element *ele = dis.lColElement(el);
      ele->LocationVector(dis,la,false);
      if (interpolVec.M() == 0)
        interpolVec.Resize(ele->NumNode()*(2*ndim+2+ndim*ndim)+2);

      ele->Evaluate(params,dis,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

      // sum values on nodes into vectors and record the touch count (build average of values)
      for (int i=0; i<ele->NumNode(); ++i)
      {
        DRT::Node* node = ele->Nodes()[i];
        const int localIndex = pressure->Map().LID(node->Id());

        if (localIndex < 0)
          continue;

        touchCount[localIndex]++;
        for (int d=0; d<ndim; ++d)
        {
          velocity->SumIntoMyValue(localIndex,d,interpolVec(d*ele->NumNode()+i));
          tracevelocity->SumIntoMyValue(localIndex,d,interpolVec((d+ndim)*ele->NumNode()+i));
        }
        for (int d=0; d<ndim*ndim; ++d)
          velocitygradient->SumIntoMyValue(localIndex,d,interpolVec(ele->NumNode()*(2*ndim+2+d)+i+2));
        (*pressure)[localIndex] += interpolVec(ele->NumNode()*(2*ndim)+i);
      }
      const int eleIndex = dis.ElementRowMap()->LID(ele->Id());
      if (eleIndex >= 0)
      {
        (*cellPres)[eleIndex] += interpolVec(ele->NumNode()*(2*ndim+2));
      }
    } // for (int el=0; el<dis.NumMyColElements();++el)

    for (int i=0; i<pressure->MyLength(); ++i)
    {
      (*pressure)[i] /= touchCount[i];
      for (int d=0; d<ndim; ++d)
      {
        (*velocity)[d][i] /= touchCount[i];
        (*tracevelocity)[d][i] /= touchCount[i];
      }
      for (int d=0; d<ndim*ndim; ++d)
        (*velocitygradient)[d][i] /= touchCount[i];
    }
    dis.ClearState(true);

    return;
  } // getNodeVectorsHDGSolid

} // namespace


/*----------------------------------------------------------------------*
 |  Output (public)                                      schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::Output(Teuchos::RCP<Epetra_MultiVector> history)
{
  TEUCHOS_FUNC_TIME_MONITOR("ACOU::AcouImplicitTimeInt::Output");

  // output of solution
  Teuchos::RCP<Epetra_Vector> interpolatedPressure, traceVel, cellPres;
  Teuchos::RCP<Epetra_MultiVector> interpolatedVelocity;
  Teuchos::RCP<Epetra_MultiVector> traceVelocity;
  Teuchos::RCP<Epetra_MultiVector> interpolatedVelocityGradient;
  if(phys_ == INPAR::ACOU::acou_lossless)
  {
    getNodeVectorsHDG(*discret_, velnp_, numdim_,
                      interpolatedVelocity, interpolatedPressure, traceVel, cellPres, phys_,padaptivity_);
  }
  else // if(phys_ == INPAR::ACOU::acou_solid)
  {
    getNodeVectorsHDGSolid(*discret_, velnp_, numdim_,
        interpolatedVelocityGradient,interpolatedVelocity,interpolatedPressure,
        traceVelocity,cellPres,phys_,writestress_);
  }
  // fill in pressure values into monitor file, if required
  FillMonitorFile(interpolatedPressure);

  if( history != Teuchos::null )
  {
    // monitor boundary condition
    std::string condname = "PressureMonitor";
    std::vector<DRT::Condition*> pressuremon;
    discret_->GetCondition(condname,pressuremon);
    const std::vector<int> pressuremonnodes = *(pressuremon[0]->Nodes());
    for(unsigned int i=0; i<pressuremonnodes.size(); ++i)
    {
      if(discret_->NodeRowMap()->LID(pressuremonnodes[i])>=0)
        history->ReplaceMyValue(history->Map().LID(pressuremonnodes[i]),step_,interpolatedPressure->operator [](discret_->NodeRowMap()->LID(pressuremonnodes[i])));
    }
  } // if( history != Teuchos::null )


  if (step_%upres_ == 0)
  {
    Teuchos::RCP<Epetra_Vector> dmap;
    if(padaptivity_)
    {
      dmap.reset(new Epetra_Vector(*discret_->ElementRowMap()));
      for(int i=0; i<discret_->NumMyRowElements(); ++i)
      {
        dmap->operator [](i) = double(discret_->lRowElement(i)->Degree());
      }
    }

    if (myrank_ == 0 && !invana_)
      std::cout<<"======= Output written in step "<<step_<<std::endl;
    // step number and time
    output_->NewStep(step_,time_);
    // write element data only once
    if (step_==0)
    {
      output_->WriteElementData(true);
      if(phys_!=INPAR::ACOU::acou_solid)
        OutputDensityAndSpeedOfSound();
    }

    output_->WriteVector("velnp",interpolatedVelocity);
    output_->WriteVector("pressure",interpolatedPressure);
    output_->WriteVector("pressure_avg",cellPres);
    if(phys_ == INPAR::ACOU::acou_lossless)
    {
      output_->WriteVector("par_vel",traceVel);
    }
    else // (phys_ == INPAR::ACOU::acou_solid)
    {
      output_->WriteVector("trace_velocity",traceVelocity);
      //output_->WriteVector("stress",interpolatedVelocityGradient,output_->nodevector);
      if(numdim_==2)
      {
        Teuchos::RCP<Epetra_Vector> stress = Teuchos::rcp(new Epetra_Vector(*discret_->NodeRowMap()));
        stress->Update(1.0,*interpolatedVelocityGradient->operator ()(0),0.0);
        output_->WriteVector("stress_xx",stress);
        stress->Update(1.0,*interpolatedVelocityGradient->operator ()(1),0.0);
        output_->WriteVector("stress_xy",stress);
        stress->Update(1.0,*interpolatedVelocityGradient->operator ()(2),0.0);
        output_->WriteVector("stress_yx",stress);
        stress->Update(1.0,*interpolatedVelocityGradient->operator ()(3),0.0);
        output_->WriteVector("stress_yy",stress);
      }
      if(numdim_==3)
      {
        output_->WriteVector("stress_xx",Teuchos::rcp(interpolatedVelocityGradient->operator ()(0)));
        output_->WriteVector("stress_xy",Teuchos::rcp(interpolatedVelocityGradient->operator ()(1)));
        output_->WriteVector("stress_xz",Teuchos::rcp(interpolatedVelocityGradient->operator ()(2)));
        output_->WriteVector("stress_yx",Teuchos::rcp(interpolatedVelocityGradient->operator ()(3)));
        output_->WriteVector("stress_yy",Teuchos::rcp(interpolatedVelocityGradient->operator ()(4)));
        output_->WriteVector("stress_yz",Teuchos::rcp(interpolatedVelocityGradient->operator ()(5)));
        output_->WriteVector("stress_zx",Teuchos::rcp(interpolatedVelocityGradient->operator ()(6)));
        output_->WriteVector("stress_zy",Teuchos::rcp(interpolatedVelocityGradient->operator ()(7)));
        output_->WriteVector("stress_zz",Teuchos::rcp(interpolatedVelocityGradient->operator ()(8)));
      }

    }

    //if(errormaps_) output_->WriteVector("error",error_);
    if(padaptivity_) output_->WriteVector("degree",dmap);

    // add restart data
    if (uprestart_ != 0 && step_%uprestart_ == 0)
    {
      WriteRestart();
    }
  }

  return;
} // Output

/*----------------------------------------------------------------------*
 |  Calculate node based values (public)                 schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::NodalPsiField(Teuchos::RCP<Epetra_Vector> outvec)
{

  // call element routine for interpolate HDG to elements
  Teuchos::ParameterList params;
  params.set<int>("action",ACOU::interpolate_psi_to_node);
  discret_->SetState(0,"trace",velnp_);
  params.set<INPAR::ACOU::PhysicalType>("physical type",phys_);
  params.set<INPAR::ACOU::DynamicType>("dynamic type",dyna_);
  params.set<double>("dt",dtp_);
  params.set<bool>("padaptivity",false);
  params.set<int>("useacouoptvecs",-6);

  std::vector<int> dummy;
  DRT::Element::LocationArray la(2);

  Epetra_SerialDenseMatrix dummyMat;
  Epetra_SerialDenseVector dummyVec;
  Epetra_SerialDenseVector interpolVec;
  std::vector<unsigned char> touchCount(outvec->MyLength());

  outvec->PutScalar(0.);

  for (int el=0; el<discret_->NumMyColElements();++el)
  {
    DRT::Element *ele = discret_->lColElement(el);
    ele->LocationVector(*discret_,la,false);
    if (interpolVec.M() == 0)
      interpolVec.Resize(ele->NumNode());

    ele->Evaluate(params,*discret_,la[0].lm_,dummyMat,dummyMat,interpolVec,dummyVec,dummyVec);

    // sum values on nodes into vectors and record the touch count (build average of values)
    for (int i=0; i<ele->NumNode(); ++i)
    {
      DRT::Node* node = ele->Nodes()[i];
      const int localIndex = outvec->Map().LID(node->Id());

      if (localIndex < 0)
        continue;

      touchCount[localIndex]++;
      (*outvec)[localIndex] += interpolVec(i);
    }
  }

  for (int i=0; i<outvec->MyLength(); ++i)
    (*outvec)[i] /= touchCount[i];

  discret_->ClearState(true);

  return;
} // NodalPsiField

/*----------------------------------------------------------------------*
 |  Fill touch count vec (needed for inverse analysis)   schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::FillTouchCountVec(Teuchos::RCP<Epetra_Vector> touchcount)
{
  // absorbing boundary conditions
  std::string condname = "PressureMonitor";
  std::vector<DRT::Condition*> pressuremon;
  discret_->GetCondition(condname,pressuremon);

  std::vector<unsigned char> touchCount(touchcount->MyLength());

  for(unsigned int i=0; i<pressuremon.size(); ++i)
  {
    std::map<int,Teuchos::RCP<DRT::Element> >& geom = pressuremon[i]->Geometry();
    std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;
    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      for(int j=0; j<curr->second->NumNode(); ++j)
      {
        DRT::Node* node = curr->second->Nodes()[j];
        const int localIndex = touchcount->Map().LID(node->Id());

        if (localIndex < 0)
          continue;

        touchCount[localIndex]++;
      }
    }
  }
  for (int i=0; i<touchcount->MyLength(); ++i)
    (*touchcount)[i] = 1.0/touchCount[i];

  return;
} // FillTouchCountVec


/*----------------------------------------------------------------------*
 |  InitMonitorFile                                      schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::InitMonitorFile()
{
  if(writemonitor_)
  {
    FILE *fp = NULL;
    if(myrank_==0)
    {
      std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
      name.append(".monitor");
      fp = fopen(name.c_str(), "w");
      if(fp == NULL)
        dserror("Couldn't open file.");
    }

    //get condition
    std::string condname="PressureMonitor";
    std::vector<DRT::Condition*>pressuremon;
    discret_->GetCondition(condname,pressuremon);
    if(pressuremon.size()>1) dserror("write of monitor file only implemented for one pressure monitor condition");
    const std::vector<int> pressuremonmics = *(pressuremon[0]->Nodes());

    int mics=pressuremonmics.size();
    int steps=0;
    if(dtp_*stepmax_<maxtime_)
      steps=stepmax_;
    else
      steps=maxtime_/dtp_+3; // first, last and int cut off

    if(myrank_ == 0)
    {
      fprintf(fp,"steps %d ",steps);
      fprintf(fp,"mics %d\n",mics);
    }

    int speakingproc=-1;
    int helptospeak=-1;
    const double* nod_coords;
    double coords[3];

    for(unsigned int n=0;n<pressuremonmics.size();++n)
    {

      if(discret_->HaveGlobalNode(pressuremonmics[n]))
      {
        helptospeak = myrank_;
        nod_coords = discret_->gNode(pressuremonmics[n])->X();
        coords[0]=nod_coords[0];
        coords[1]=nod_coords[1];
        coords[2]=nod_coords[2];
      }
      else
        helptospeak = 0;
      discret_->Comm().MaxAll(&helptospeak,&speakingproc,1);
      discret_->Comm().Broadcast(coords,3,speakingproc);

      if(myrank_==0)
        fprintf(fp,"%e %e %e\n",coords[0],coords[1],coords[2]);
    }
    if(myrank_ == 0)
    {
      fprintf(fp,"#\n#\n#\n");
      fclose(fp);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  FillMonitorFile                                      schoeder 04/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouTimeInt::FillMonitorFile(Teuchos::RCP<Epetra_Vector> ip)
{
  if(writemonitor_)
  {
    FILE *fp = NULL;
    if(myrank_==0)
    {
     std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
     name.append(".monitor");
     fp = fopen(name.c_str(), "a");
    }

    // get condition
    std::string condname="PressureMonitor";
    std::vector<DRT::Condition*>pressuremon;
    discret_->GetCondition(condname,pressuremon);
    const std::vector<int> pressuremonmics = *(pressuremon[0]->Nodes());
    int mics=pressuremonmics.size();

    if(myrank_==0) fprintf(fp,"%e ",time_);
    int helptospeak=-1;
    int speakingproc=-1;
    double pressure = 0.0;

    for(int n=0;n<mics;n++)
    {
      if(discret_->HaveGlobalNode(pressuremonmics[n]))
      {
        if(ip->Map().LID(pressuremonmics[n])>=0)
        {
          helptospeak=myrank_;
          pressure=ip->operator [](ip->Map().LID(pressuremonmics[n]));
        }
      }
      else
        helptospeak = -1;
      discret_->Comm().MaxAll(&helptospeak,&speakingproc,1);
      discret_->Comm().Broadcast(&pressure,1,speakingproc);

      if(myrank_==0) fprintf(fp,"%e ", pressure);
    }
    if(myrank_==0)
    {
      fprintf(fp,"\n");
      fclose(fp);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Return discretization (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> ACOU::AcouTimeInt::Discretization()
{
  return Teuchos::rcp_dynamic_cast<DRT::Discretization>(discret_);
} // Discretization


/*----------------------------------------------------------------------*
 |  Create test field (public)                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ACOU::AcouTimeInt::CreateFieldTest()
{
  return Teuchos::rcp(new AcouResultTest(*this));
} // CreateFieldTest
