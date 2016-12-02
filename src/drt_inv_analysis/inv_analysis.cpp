/*----------------------------------------------------------------------*/
/*!
 * \file inv_analysis.cpp

\brief inv analysis
\level 2

<pre>
\maintainer Sebastian Kehl
</pre>
*/
/*----------------------------------------------------------------------*/


#include "inv_analysis.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_structure/strtimint.H"
#include "../drt_structure/strtimint_create.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_matelast/elast_coupblatzko.H"
#include "../drt_matelast/elast_couplogneohooke.H"
#include "../drt_matelast/elast_couplogmixneohooke.H"
#include "../drt_matelast/elast_coupexppol.H"
#include "../drt_matelast/elast_coupneohooke.H"
#include "../drt_matelast/elast_isoexpopow.H"
#include "../drt_matelast/elast_coupmooneyrivlin.H"
#include "../drt_matelast/elast_isomooneyrivlin.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_isoyeoh.H"
#include "../drt_matelast/elast_iso1pow.H"
#include "../drt_matelast/elast_iso2pow.H"
#include "../drt_matelast/elast_coup1pow.H"
#include "../drt_matelast/elast_coup2pow.H"
#include "../drt_matelast/elast_volpenalty.H"
#include "../drt_matelast/elast_vologden.H"
#include "../drt_matelast/elast_volsussmanbathe.H"
#include "../drt_matelast/visco_isoratedep.H"
#include "../drt_mat/material.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/viscoelasthyper.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::InvAnalysis::InvAnalysis(Teuchos::RCP<DRT::Discretization> dis,
                                Teuchos::RCP<LINALG::Solver> solver,
                                Teuchos::RCP<IO::DiscretizationWriter> output)
  : discret_(dis),
    solver_(solver),
    output_(output),
    sti_(Teuchos::null)
{

  // Getting boundary conditions
  discret_->GetCondition("SurfaceNeumann",surfneum_ );
  discret_->GetCondition("Dirichlet",surfdir_ );
  reset_out_count_=0;

  if (surfneum_.size()>1)
    dserror("The inverse analysis only works for 1 NBC with 2 DBC or for 2 DBC!");

  // input parameters structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // measured points, gives the number how many displacment steps are
  // measured
  nmp_   = 2*sdyn.get<int>("NUMSTEP");
  tstep_ = sdyn.get<double>("TIMESTEP");
  // get total timespan of simulation 0.5 due to factor 2 in nmp_
  double ttime_ = nmp_*tstep_*0.5;
  // input parameters inverse analysis
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();

  //  tolerance for the curve fitting
  tol_ = iap.get<double>("INV_ANA_TOL");

  // experimentally measured curve
  {
    mcurve_   = Epetra_SerialDenseVector(nmp_);  //
    double cpx0 = iap.get<double>("MC_X_0");
    double cpx1 = iap.get<double>("MC_X_1");
    double cpx2 = iap.get<double>("MC_X_2");
    double cpy0 = iap.get<double>("MC_Y_0");
    double cpy1 = iap.get<double>("MC_Y_1");
    double cpy2 = iap.get<double>("MC_Y_2");
    for (int i=0; i<nmp_; i++)
    {
      mcurve_[i] = cpx0*(1-exp(-pow(((1000./ttime_)*cpx1*(i)*ttime_/nmp_), cpx2)));
      i=i+1;
      mcurve_[i] = cpy0*(1-exp(-pow(((1000./ttime_)*cpy1*(i-1)*ttime_/nmp_), cpy2)));
     }
  }
  //dserror("Halt");

  // error: diference of the measured to the calculated curve
  error_  = 1.0E6;
  error_o_= 1.5E6;


  switch(DRT::INPUT::IntegralValue<INPAR::STR::RegStratUpdate>(iap,"UPDATE_REG"))
   {
     case INPAR::STR::reg_update_res:
     {
      // dummy different strategies are only implemented in generalized inv ana  reg_update_ = res_based;
     }
     break;
     default:
       dserror("Only RES based update strategy implemented for old inverse analysis switch to inv_generalized to use other update strategies");
       break;
   }
  if (DRT::INPUT::IntegralValue<int>(iap,"PARAM_BOUNDS"))
    dserror("Bounds Check only implemented for gen_inv_ana");

  // trainings parameter
  mu_  = 1.;
  tol_mu_ = tol_;

  // list of materials for each problem instance that should be fitted
  for (unsigned prob=0; prob<DRT::Problem::NumInstances(); ++prob)
  {
    const Teuchos::ParameterList& myiap = DRT::Problem::Instance(prob)->InverseAnalysisParams();
    std::set<int> myset;

    int word1;
    std::istringstream matliststream(Teuchos::getNumericStringParameter(myiap,"INV_LIST"));
    while (matliststream >> word1)
    {
      if (word1!=-1) // this means there was no matlist specified in the input file
        myset.insert(word1);
    }
    matset_.push_back(myset);
  }

  // list of elasthyper summands for each problem instance that should be fitted
  for (unsigned prob=0; prob<DRT::Problem::NumInstances(); ++prob)
  {
    const Teuchos::ParameterList& myiap = DRT::Problem::Instance(prob)->InverseAnalysisParams();
    std::set<int> myset;

    int word1;
    std::istringstream matliststream(Teuchos::getNumericStringParameter(myiap,"INV_EH_LIST"));
    while (matliststream >> word1)
    {
      if (word1!=-1) // this means there was no matlist specified in the input file
        myset.insert(word1);
    }
    eh_matset_.push_back(myset);
  }

  // read material parameters from input file
  ReadInParameters();

  // Number of material parameters
  np_ = p_.Length();

  // controlling parameter
  numb_run_ =  0;         // counter of how many runs were made in the inverse analysis
}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::InvAnalysis::Integrate()
{
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  double alpha  = iap.get<double>("INV_ALPHA");
  double beta   = iap.get<double>("INV_BETA");
  int max_itter = iap.get<int>("INV_ANA_MAX_RUN");
  int newfiles = DRT::INPUT::IntegralValue<int>(iap,"NEW_FILES");

  // fitting loop
  do
  {
    std::vector<double> inc(np_, 0.0);
    for (int i=0; i<np_;i++)
      inc[i] = alpha + beta * p_[i];

    Epetra_SerialDenseMatrix cmatrix(nmp_, np_+1);

    if (discret_->Comm().MyPID()==0)
      std::cout << "-----------------------------making Jacobian matrix-------------------------" <<std::endl;
    for (int i=0; i<np_+1;i++)
    {
      bool outputtofile = false;
      // output only for last run
      if (i==np_) outputtofile = true;

      if (outputtofile)
      {
        if (newfiles)
          output_->NewResultFile((numb_run_));
        else
          output_->OverwriteResultFile();

        output_->WriteMesh(0,0.0);
      }

      // Multi-scale: if an inverse analysis is performed on the micro-level,
      // the time and step need to be reset now. Furthermore, the result file
      // needs to be opened.
      MultiInvAnaInit();

      if (discret_->Comm().MyPID()==0)
        std::cout << "------------------------------- run "<< i+1 << " of: " << np_+1 <<" ---------------------------------" <<std::endl;
      Epetra_SerialDenseVector p_cur = p_;
      if (i!= np_)
        p_cur[i]=p_[i] + inc[i];
      SetParameters(p_cur);

      // compute nonlinear problem and obtain computed displacements
      Epetra_SerialDenseVector cvector;
      cvector = CalcCvector(outputtofile);

      for (int j=0; j<nmp_;j++)
        cmatrix(j, i)=cvector[j];
    }

    discret_->Comm().Barrier();
    for (int proc=0; proc<discret_->Comm().NumProc(); ++proc)
    {
      if (proc==discret_->Comm().MyPID())
      {
        if (proc == 0)
        {
          CalcNewParameters(cmatrix,  inc);
        }
      }
    }

    discret_->Comm().Barrier();
    // set new material parameters
    SetParameters(p_);
    numb_run_++;
    discret_->Comm().Broadcast(&error_o_,1,0);
    discret_->Comm().Broadcast(&error_,1,0);
    discret_->Comm().Broadcast(&numb_run_,1,0);
  }while(numb_run_<max_itter && error_>tol_) ;      // while (abs(error_o_-error_)>0.001 && error_>tol_ && numb_run_<max_itter);

  discret_->Comm().Barrier();
  for (int proc=0; proc<discret_->Comm().NumProc(); ++proc)
    if (proc==discret_->Comm().MyPID())
      if (proc == 0)
        PrintFile();
  discret_->Comm().Barrier();
  return;
}

void STR::InvAnalysis::CalcNewParameters(Epetra_SerialDenseMatrix cmatrix,  std::vector<double> inc)
{
  // initalization of the Jacobi and storage matrix
  Epetra_SerialDenseMatrix J(nmp_, np_);
  Epetra_SerialDenseMatrix sto(np_,  np_);
  Epetra_SerialDenseMatrix sto2(np_, nmp_);
  Epetra_SerialDenseVector delta_p(np_);
  Epetra_SerialDenseVector rcurve(nmp_);

  for (int i=0; i<nmp_; i++)
    for (int j=0; j<np_; j++)
      J(i, j) = (cmatrix(i, j)-cmatrix(i, np_)) / inc[j];     //calculating J(p)

  sto.Multiply('T',  'N',  1,  J, J,  0);     //calculating J.T*J

  for (int i=0; i<np_; i++)
    sto(i, i) += (mu_*sto(i, i));

  LINALG::NonSymmetricInverse(sto,  np_);     //calculating (J.T*J+mu*I).I
  sto2.Multiply('N', 'T', 1,  sto, J, 0);     //calculating (J.T*J+mu*I).I*J.T

  for (int i = 0; i<nmp_; i++)
    rcurve[i]=mcurve_[i]-cmatrix(i, np_);
  delta_p.Multiply('N', 'N', 1., sto2, rcurve, 0.);

  for (int i=0;i<np_;i++)
    p_[i] += delta_p[i];//*p_[i];

  // dependent on the # of steps
  error_o_   = error_;
  error_   = rcurve.Norm2()/sqrt(nmp_);

  //Adjust training parameter
  mu_ *= (error_/error_o_);
  if (numb_run_==0)
    mu_=1.;

  PrintStorage(cmatrix,  delta_p);
  return;
}

Epetra_SerialDenseVector STR::InvAnalysis::CalcCvector(bool outputtofile)
{
  // get input parameter lists
  const Teuchos::ParameterList& ioflags
    = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();
  Teuchos::ParameterList xparams;
  xparams.set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  xparams.set<int>("REDUCED_OUTPUT",0);

  // create time integrator
  sti_ = TimIntCreate(sdyn, ioflags, sdyn, xparams, discret_, solver_, solver_, output_);
  sti_->Init(sdyn,sdyn,xparams,discret_,solver_);
  sti_->Setup();
  if (sti_ == Teuchos::null) dserror("Failed in creating integrator.");

  // initialize time loop
  const int stepmax = sti_->NumStep();
  Epetra_SerialDenseVector cvector(2*stepmax);

  // time loop
  while ( sti_->NotFinished() )
  {
    // call the predictor
    sti_->PrepareTimeStep();

    // integrate time step
    // after this step we hold disn_, etc
    sti_->Solve();

    // calculate stresses, strains, energies
    sti_->PrepareOutput();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    sti_->PreUpdate();
    sti_->UpdateStepState();

    // gets the displacements per timestep
    {
      // get current step n
      int step = sti_->StepOld();

      Epetra_SerialDenseVector cvector_arg = GetCalculatedCurve();
      cvector[2*step]   = cvector_arg[0];
      cvector[2*step+1] = cvector_arg[1];
    }

    //dserror("Halt");
    // update time and step
    sti_->UpdateStepTime();

    // Update Element
    sti_->UpdateStepElement();
    sti_->PostUpdate();

    // print info about finished time step
    sti_->PrintStep();

    // write output
    if (outputtofile) sti_->Output();

  }
  return cvector;
}


/*----------------------------------------------------------------------*/
/* */
Epetra_SerialDenseVector STR::InvAnalysis::GetCalculatedCurve()
{
  Epetra_SerialDenseVector cvector_arg(2);

  // current displacement vector
  Teuchos::RCP<Epetra_Vector> disp = sti_->DisNew();

  std::vector<DRT::Condition*> invanacond;
  discret_->GetCondition("SurfInvAna",  invanacond);

  //nodes of the pulling direction
  const std::vector<int>* ia_nd_ps  = invanacond[0]->Nodes();

  for (std::vector<int>::const_iterator inodegid = ia_nd_ps->begin();
       inodegid !=ia_nd_ps->end();
       ++inodegid)
  {
    if (discret_->HaveGlobalNode(*inodegid))
    {
      if (discret_->gNode(*inodegid)->Owner() == discret_->Comm().MyPID())
      {
        const DRT::Node* node = discret_->gNode(*inodegid);
        std::vector<int> lm = discret_->Dof(node);
        const double disp_x = (*disp)[disp->Map().LID(lm[0])];
        cvector_arg[0] += disp_x;
      }
    }
  }

  {
    double test;
    discret_->Comm().SumAll(&cvector_arg[0],&test,1);
    cvector_arg[0] = test/(*ia_nd_ps).size();
  }


  //nodes to determine the compression
  const std::vector<int>* ia_nd_fs_p = invanacond[1]->Nodes();
  const std::vector<int>* ia_nd_fs_n = invanacond[2]->Nodes();

  for (std::vector<int>::const_iterator inodegid = ia_nd_fs_p->begin();
       inodegid !=ia_nd_fs_p->end();
       ++inodegid)
  {
    if (discret_->HaveGlobalNode(*inodegid))
    {
      if (discret_->gNode(*inodegid)->Owner() == discret_->Comm().MyPID())
      {
        const DRT::Node* node = discret_->gNode(*inodegid);
        std::vector<int> lm = discret_->Dof(node);
        const double disp_y = (*disp)[disp->Map().LID(lm[1])];
        cvector_arg[1] += disp_y;
      }
    }
  }

  for (std::vector<int>::const_iterator inodegid = ia_nd_fs_n->begin();
       inodegid !=ia_nd_fs_n->end();
       ++inodegid)
  {
    if (discret_->HaveGlobalNode(*inodegid))
    {
      if (discret_->gNode(*inodegid)->Owner() == discret_->Comm().MyPID())
      {
        const DRT::Node* node = discret_->gNode(*inodegid);
        std::vector<int> lm = discret_->Dof(node);
        const double disp_y = (*disp)[disp->Map().LID(lm[1])];
        cvector_arg[1] -= disp_y;
      }
    }
  }

  {
    double test;
    discret_->Comm().SumAll(&cvector_arg[1],&test,1);
    cvector_arg[1] = test/(2.*(*ia_nd_fs_p).size());
  }
  return cvector_arg;
}

/*----------------------------------------------------------------------*/
/* */
void STR::InvAnalysis::PrintStorage(Epetra_SerialDenseMatrix cmatrix,  Epetra_SerialDenseVector delta_p)
{

  // store the error and mu_

  p_s_.Reshape(numb_run_+1,  np_);
  for (int i=0; i<np_; i++)
    p_s_(numb_run_, i)=p_(i);

  delta_p_s_.Reshape(numb_run_+1,  np_);
  for (int i=0; i<np_; i++)
    delta_p_s_(numb_run_, i)=delta_p(i);

  ccurve_s_.Reshape(nmp_,  numb_run_+1);
  for (int i=0; i<nmp_; i++)
    ccurve_s_(i, numb_run_)= cmatrix(i, cmatrix.ColDim()-1);

  mu_s_.Resize(numb_run_+1);
  mu_s_(numb_run_)=mu_;

  error_s_.Resize(numb_run_+1);
  error_s_(numb_run_) = error_;
  // print error and parameter

  std::cout << std::endl;
  printf("################################################");
  printf("##############################################\n");
  printf("############################ Inverse Analysis ##");
  printf("##############################################\n");
  printf("################################### run ########");
  printf("##############################################\n");
  printf("################################### %3i ########",  numb_run_);
  printf("##############################################\n");
  printf("################################################");
  printf("##############################################\n");

  for (int i=0; i < numb_run_+1; i++)
  {
    printf("Error: ");
    printf("%10.3f", error_s_(i));
    printf("\tParameter: ");
    for (int j=0; j < delta_p.Length(); j++)
      printf("%10.3f", p_s_(i, j));
    //printf("\tDelta_p: ");
    //for (int j=0; j < delta_p.Length(); j++)
    //  printf("%10.3f", delta_p_s_(i, j));
    printf("\tmu: ");
    printf("%10.3f", mu_s_(i));
    printf("\n");
  }

  printf("\n");
  for (int i=0; i < nmp_/2.; i++)
  {
    printf(" %10.2f ",  mcurve_(i*2));
    if (numb_run_<15)
    {
      for (int j=0; j<numb_run_+1; j++)
        printf(" %10.2f ",  ccurve_s_((i)*2, j));
    }
    else
    {
      for (int j=numb_run_-14; j<numb_run_+1; j++)
        printf(" %10.2f ",  ccurve_s_((i)*2, j));
    }
    printf("\n");
  }

  printf("\n");

  for (int i=0; i < nmp_/2.; i++)
  {
    printf(" %10.2f ",  mcurve_((i)*2+1));
    if (numb_run_<15)
    {
      for (int j=0; j<numb_run_+1; j++)
        printf(" %10.2f ",  ccurve_s_((i)*2+1, j));
    }
    else
    {
      for (int j=numb_run_-14; j<numb_run_+1; j++)
        printf(" %10.2f ",  ccurve_s_((i)*2+1, j));
    }
    printf("\n");
  }

  printf("################################################");
  printf("##############################################\n");
  std::cout << std::endl;
}


void STR::InvAnalysis::PrintFile()
{
  FILE * cxFile;
  FILE * cyFile;
  FILE * pFile;

  std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
  name.append(filename_);

  if (name.rfind("_run_")!=std::string::npos)
  {
    size_t pos = name.rfind("_run_");
    if (pos==std::string::npos)
      dserror("inconsistent file name");
    name = name.substr(0, pos);
  }

  std::string gp     = name+"_plot.gp";
  std::string xcurve = name+"_Curve_x.txt";
  std::string ycurve = name+"_Curve_y.txt";
  std::string para   = name+"_Para.txt";

  cxFile = fopen((xcurve).c_str(), "w");
  for (int i=0; i < nmp_/2.; i++)
  {
    fprintf(cxFile, " %10.5f ,",  mcurve_(i*2));
    for (int j=0; j<numb_run_; j++)
      fprintf(cxFile, " %10.5f ,",  ccurve_s_(i*2, j));
    fprintf(cxFile, "\n");
  }
  fclose(cxFile);

  cyFile = fopen((ycurve).c_str(), "w");
  for (int i=0; i < nmp_/2.; i++)
  {
    fprintf(cyFile, " %10.5f ,",  mcurve_((i)*2+1));
    for (int j=0; j<numb_run_; j++)
      fprintf(cyFile, " %10.5f ,",  ccurve_s_((i)*2+1, j));
    fprintf(cyFile, "\n");
  }
  fclose(cyFile);

  pFile  = fopen((para).c_str(), "w");
  fprintf(pFile, "#Error       Parameter    Delta_p      mu \n");
  for (int i=0; i < numb_run_; i++)
  {
    fprintf(pFile, "%10.3f,", error_s_(i));
    for (int j=0; j < np_; j++)
      fprintf(pFile, "%10.3f,", p_s_(i, j));
    for (int j=0; j < np_; j++)
      fprintf(pFile, "%10.3f,", delta_p_s_(i, j));
    fprintf(pFile, "%10.3f", mu_s_(i));
    fprintf(pFile, "\n");
  }
  fclose(pFile);

  numb_run_=numb_run_-1;
}


void STR::InvAnalysis::ReadInParameters()
{
  for (unsigned prob=0; prob<DRT::Problem::NumInstances(); ++prob)
  {
    const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance(prob)->Materials()->Map();
    std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator curr;

    for (curr=mats.begin(); curr != mats.end(); ++curr)
    {
      const Teuchos::RCP<MAT::PAR::Material> material = curr->second;
      std::set<int> mymatset = matset_[prob];

      if (mymatset.size()==0 or mymatset.find(material->Id())!=mymatset.end())
      {
        switch(material->Type())
        {
        case INPAR::MAT::m_elasthyper:
        {
          // Create a pointer on the Material
          const MAT::PAR::ElastHyper* params = dynamic_cast<const MAT::PAR::ElastHyper*>(material->Parameter());
          const int nummat               = params->nummat_;
          const std::vector<int>* matids = params->matids_;

          // For each of the summands of the hyperelastic material we need to add the
          // parameters to the inverse analysis

          // Problems with beta, is it the only negative parameter? Maybe
          // we should exclude it

          for (int i=0; i<nummat; ++i)
          {
            //get the material of the summand
            filename_=filename_+"_";

            const int id = (*matids)[i];

            std::set<int> myehmatset = eh_matset_[prob];

            if (myehmatset.size()==0 or myehmatset.find(id)!=myehmatset.end())
            {
              const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;

              switch (actelastmat->Type())
              {
              case INPAR::MAT::mes_couplogneohooke:
              {
                filename_=filename_+"_couplogneohooke";
                const MAT::ELASTIC::PAR::CoupLogNeoHooke* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupLogNeoHooke*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->mue_;
                p_[j+1] = params2->lambda_;
                break;
              }
              case INPAR::MAT::mes_coupexppol:
              {
                filename_=filename_+"_coupexppol";
                const MAT::ELASTIC::PAR::CoupExpPol* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupExpPol*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+3);
                p_[j]   = params2->a_;
                p_[j+1] = params2->b_;
                p_[j+2] = params2->c_;
                break;
              }
              case INPAR::MAT::mes_coupneohooke:
              {
                filename_=filename_+"_coupneohooke";
                const MAT::ELASTIC::PAR::CoupNeoHooke* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupNeoHooke*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j] = params2->c_;
                p_[j+1] = params2->beta_;
                break;
              }
              case INPAR::MAT::mes_coupblatzko:
              {
                filename_=filename_+"_coupblatzko";
                const MAT::ELASTIC::PAR::CoupBlatzKo* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupBlatzKo*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->mue_;
                p_[j+1] = (1./(1.-2.*params2->nue_))-1.;
                //p_[j+1] = params2->f_;
                break;
              }
              case INPAR::MAT::mes_isoneohooke:
              {
                filename_=filename_+"_isoneohooke";
                const MAT::ELASTIC::PAR::IsoNeoHooke* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoNeoHooke*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->mue_;
                break;
              }
              case INPAR::MAT::mes_isoyeoh:
              {
                filename_=filename_+"_isoyeoh";
                const MAT::ELASTIC::PAR::IsoYeoh* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+3);
                p_[j]   = params2->c1_;
                p_[j+1] = params2->c2_;
                p_[j+2] = params2->c3_;
                break;
              }
              case INPAR::MAT::mes_iso1pow:
              {
                filename_=filename_+"_iso1pow";
                const MAT::ELASTIC::PAR::Iso1Pow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->c_;
                break;
              }
              case INPAR::MAT::mes_iso2pow:
              {
                filename_=filename_+"_iso2pow";
                const MAT::ELASTIC::PAR::Iso2Pow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::Iso2Pow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->c_;
                break;
              }
              case INPAR::MAT::mes_coup1pow:
              {
                filename_=filename_+"_coup1pow";
                const MAT::ELASTIC::PAR::Coup1Pow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::Coup1Pow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->c_;
                break;
              }
              case INPAR::MAT::mes_coup2pow:
              {
                filename_=filename_+"_coup2pow";
                const MAT::ELASTIC::PAR::Coup2Pow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::Coup2Pow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->c_;
                break;
              }
              case INPAR::MAT::mes_coupmooneyrivlin:
              {
                filename_=filename_+"_coupmooneyrivlin";
                const MAT::ELASTIC::PAR::CoupMooneyRivlin* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupMooneyRivlin*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+3);
                p_[j]   = params2->c1_;
                p_[j+1] = params2->c2_;
                p_[j+2] = params2->c3_;
                break;
              }
              case INPAR::MAT::mes_isoexpopow:
              {
                filename_=filename_+"_isoexpopow";
                const MAT::ELASTIC::PAR::IsoExpoPow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoExpoPow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->k1_;
                p_[j+1] = params2->k2_;
                break;
              }
              case INPAR::MAT::mes_isomooneyrivlin:
              {
                filename_=filename_+"_isomooneyrivlin";
                const MAT::ELASTIC::PAR::IsoMooneyRivlin* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoMooneyRivlin*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->c1_;
                p_[j+1] = params2->c2_;
                break;
              }
              case INPAR::MAT::mes_volsussmanbathe:
              {
                filename_=filename_+"_volsussmanbathe";
                const MAT::ELASTIC::PAR::VolSussmanBathe* params2 = dynamic_cast<const MAT::ELASTIC::PAR::VolSussmanBathe*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->kappa_;
                break;
              }
              case INPAR::MAT::mes_volpenalty:
              {
                filename_=filename_+"_volpenalty";
                const MAT::ELASTIC::PAR::VolPenalty* params2 = dynamic_cast<const MAT::ELASTIC::PAR::VolPenalty*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->eps_;
                p_[j+1]   = params2->gam_;
                break;
              }
              case INPAR::MAT::mes_vologden:
              {
                filename_=filename_+"_vologden";
                const MAT::ELASTIC::PAR::VolOgden* params2 = dynamic_cast<const MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->kappa_;
                //p_[j+1] = params2->beta_;
                break;
              }
              default:
                dserror("cannot deal with this material");
                break;
              }
            }
          }
          break;
        }
        // viscoelasthyper
        case INPAR::MAT::m_viscoelasthyper:
        {
          // Create a pointer on the Material
          const MAT::PAR::ViscoElastHyper* params = dynamic_cast<const MAT::PAR::ViscoElastHyper*>(material->Parameter());
          const int nummat               = params->nummat_;
          const std::vector<int>* matids = params->matids_;

          // For each of the summands of the viscohyperelastic material we need to add the
          // parameters to the inverse analysis

          // Problems with beta, is it the only negative parameter? Maybe
          // we should exclude it

          for (int i=0; i<nummat; ++i)
          {
            //get the material of the summand
            filename_=filename_+"_";

            const int id = (*matids)[i];

            std::set<int> myehmatset = eh_matset_[prob];

            if (myehmatset.size()==0 or myehmatset.find(id)!=myehmatset.end())
            {
              const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;

              switch (actelastmat->Type())
              {
              case INPAR::MAT::mes_couplogneohooke:
              {
                filename_=filename_+"_couplogneohooke";
                const MAT::ELASTIC::PAR::CoupLogNeoHooke* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupLogNeoHooke*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->mue_;
                p_[j+1] = params2->lambda_;
                break;
              }
              case INPAR::MAT::mes_coupexppol:
              {
                filename_=filename_+"_coupexppol";
                const MAT::ELASTIC::PAR::CoupExpPol* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupExpPol*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+3);
                p_[j]   = params2->a_;
                p_[j+1] = params2->b_;
                p_[j+2] = params2->c_;
                break;
              }
              case INPAR::MAT::mes_coupneohooke:
              {
                filename_=filename_+"_coupneohooke";
                const MAT::ELASTIC::PAR::CoupNeoHooke* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupNeoHooke*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j] = params2->c_;
                p_[j+1] = params2->beta_;
                break;
              }
              case INPAR::MAT::mes_coupblatzko:
              {
                filename_=filename_+"_coupblatzko";
                const MAT::ELASTIC::PAR::CoupBlatzKo* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupBlatzKo*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->mue_;
                p_[j+1] = (1./(1.-2.*params2->nue_))-1.;
                //p_[j+1] = params2->f_;
                break;
              }
              case INPAR::MAT::mes_isoneohooke:
              {
                filename_=filename_+"_isoneohooke";
                const MAT::ELASTIC::PAR::IsoNeoHooke* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoNeoHooke*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->mue_;
                break;
              }
              case INPAR::MAT::mes_isoyeoh:
              {
                filename_=filename_+"_isoyeoh";
                const MAT::ELASTIC::PAR::IsoYeoh* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+3);
                p_[j]   = params2->c1_;
                p_[j+1] = params2->c2_;
                p_[j+2] = params2->c3_;
                break;
              }
              case INPAR::MAT::mes_iso1pow:
              {
                filename_=filename_+"_iso1pow";
                const MAT::ELASTIC::PAR::Iso1Pow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->c_;
                break;
              }
              case INPAR::MAT::mes_iso2pow:
              {
                filename_=filename_+"_iso2pow";
                const MAT::ELASTIC::PAR::Iso2Pow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::Iso2Pow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->c_;
                break;
              }
              case INPAR::MAT::mes_coup1pow:
              {
                filename_=filename_+"_coup1pow";
                const MAT::ELASTIC::PAR::Coup1Pow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::Coup1Pow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->c_;
                break;
              }
              case INPAR::MAT::mes_coup2pow:
              {
                filename_=filename_+"_coup2pow";
                const MAT::ELASTIC::PAR::Coup2Pow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::Coup2Pow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->c_;
                break;
              }
              case INPAR::MAT::mes_coupmooneyrivlin:
              {
                filename_=filename_+"_coupmooneyrivlin";
                const MAT::ELASTIC::PAR::CoupMooneyRivlin* params2 = dynamic_cast<const MAT::ELASTIC::PAR::CoupMooneyRivlin*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+3);
                p_[j]   = params2->c1_;
                p_[j+1] = params2->c2_;
                p_[j+2] = params2->c3_;
                break;
              }
              case INPAR::MAT::mes_isoexpopow:
              {
                filename_=filename_+"_isoexpopow";
                const MAT::ELASTIC::PAR::IsoExpoPow* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoExpoPow*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->k1_;
                p_[j+1] = params2->k2_;
                break;
              }
              case INPAR::MAT::mes_isomooneyrivlin:
              {
                filename_=filename_+"_isomooneyrivlin";
                const MAT::ELASTIC::PAR::IsoMooneyRivlin* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoMooneyRivlin*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->c1_;
                p_[j+1] = params2->c2_;
                break;
              }
              case INPAR::MAT::mes_volsussmanbathe:
              {
                filename_=filename_+"_volsussmanbathe";
                const MAT::ELASTIC::PAR::VolSussmanBathe* params2 = dynamic_cast<const MAT::ELASTIC::PAR::VolSussmanBathe*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->kappa_;
                break;
              }
              case INPAR::MAT::mes_volpenalty:
              {
                filename_=filename_+"_volpenalty";
                const MAT::ELASTIC::PAR::VolPenalty* params2 = dynamic_cast<const MAT::ELASTIC::PAR::VolPenalty*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+2);
                p_[j]   = params2->eps_;
                p_[j+1]   = params2->gam_;
                break;
              }
              case INPAR::MAT::mes_vologden:
              {
                filename_=filename_+"_vologden";
                const MAT::ELASTIC::PAR::VolOgden* params2 = dynamic_cast<const MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->kappa_;
                //p_[j+1] = params2->beta_;
                break;
              }
              case INPAR::MAT::mes_isoratedep:
              {
                filename_=filename_+"_isoratedep";
                const MAT::ELASTIC::PAR::IsoRateDep* params2 = dynamic_cast<const MAT::ELASTIC::PAR::IsoRateDep*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j+1);
                p_[j]   = params2->n_;
                break;
              }
              default:
                dserror("cannot deal with this material");
                break;
              }
            }
          }
          break;
        }

        // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block or an interface to a micro material
        case INPAR::MAT::mes_couplogneohooke:
        case INPAR::MAT::mes_coupexppol:
        case INPAR::MAT::mes_coupneohooke:
        case INPAR::MAT::mes_coupblatzko:
        case INPAR::MAT::mes_isoneohooke:
        case INPAR::MAT::mes_isoyeoh:
        case INPAR::MAT::mes_iso1pow:
        case INPAR::MAT::mes_iso2pow:
        case INPAR::MAT::mes_coup1pow:
        case INPAR::MAT::mes_coup2pow:
        case INPAR::MAT::mes_isoexpopow:
        case INPAR::MAT::mes_isomooneyrivlin:
        case INPAR::MAT::mes_coupmooneyrivlin:
        case INPAR::MAT::mes_volsussmanbathe:
        case INPAR::MAT::mes_volpenalty:
        case INPAR::MAT::mes_vologden:
        case INPAR::MAT::mes_isoratedep:
        case INPAR::MAT::m_struct_multiscale:
          break;
        default:
          dserror("The inverse analysis is not implemented for this material");
          break;
        }
      }
    }
  }
  return;
}

void STR::InvAnalysis::SetParameters(Epetra_SerialDenseVector p_cur)
{

  discret_->Comm().Broadcast(&p_cur[0],np_,0);

  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();
  if(subcomm != Teuchos::null)
  {
    // tell supporting procs that material parameters have to be set
    int task[2] = {6, np_};
    subcomm->Broadcast(task, 2, 0);
    // broadcast p_cur to the micro scale
    subcomm->Broadcast(&p_cur[0], np_, 0);

    if(DRT::Problem::Instance(0)->GetNPGroup()->SubComm()->NumProc() > 1 and ( matset_[0].size()!=0 or eh_matset_[0].size()!=0))
      dserror("No fitting of macro material in case of nested parallelism for multi scale inverse analysis.");
  }

  // write new material parameter

  for (unsigned prob=0; prob<DRT::Problem::NumInstances(); ++prob)
  {
    std::set<int> mymatset = matset_[prob];
    std::set<int> myehmatset = eh_matset_[prob];
    if(subcomm != Teuchos::null)
    {
      // broadcast sets within micro scale once per problem instance
      LINALG::GatherAll<int>(mymatset, *subcomm);
      LINALG::GatherAll<int>(myehmatset, *subcomm);
    }

    // material parameters are set for the current problem instance
    STR::SetMaterialParameters(prob, p_cur, mymatset, myehmatset);
  }
}


void STR::SetMaterialParameters(int prob, Epetra_SerialDenseVector& p_cur, std::set<int>& mymatset, std::set<int>& myehmatset)
{
  const std::map<int,Teuchos::RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance(prob)->Materials()->Map();
  std::map<int,Teuchos::RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    const Teuchos::RCP<MAT::PAR::Material> material = curr->second;
    if (mymatset.size()==0 or mymatset.find(material->Id())!=mymatset.end())
    {
      if (material->Type() == INPAR::MAT::m_elasthyper)
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(material->Parameter());
        const int nummat               = params->nummat_;
        const std::vector<int>* matids = params->matids_;

        // For each of the summands of the hyperelastic material we need to add the
        // parameters to the inverse analysis

        // Problems with beta, is it the only negative parameter? Maybe
        // we should exclude it

        //itterator to go through the parameters
        int j = 0;

        for (int i=0; i<nummat; i++)
        {
          //get the material of the summand
          const int id = (*matids)[i];

          if (myehmatset.size()==0 or myehmatset.find(id)!=myehmatset.end())
          {
            const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;

            switch (actelastmat->Type())
            {
            case INPAR::MAT::mes_couplogneohooke:
            {
              MAT::ELASTIC::PAR::CoupLogNeoHooke* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupLogNeoHooke*>(actelastmat->Parameter());
              params2->SetMue(abs(p_cur(j)));
              params2->SetLambda(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_coupexppol:
            {
              MAT::ELASTIC::PAR::CoupExpPol* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupExpPol*>(actelastmat->Parameter());
              params2->SetA(abs(p_cur(j)));
              params2->SetB(abs(p_cur(j+1)));
              params2->SetC(abs(p_cur(j+2)));
              j = j+3;
              break;
            }
            case INPAR::MAT::mes_coupneohooke:
            {
              MAT::ELASTIC::PAR::CoupNeoHooke* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupNeoHooke*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              params2->SetBeta(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_coupblatzko:
            {
              MAT::ELASTIC::PAR::CoupBlatzKo* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupBlatzKo*>(actelastmat->Parameter());
              params2->SetMue(abs(p_cur(j)));
              //params2->SetF(abs(p_cur(j+1)));
              params2->SetNue((abs(p_cur(j+1)))/(2.*(abs(p_cur(j+1))+1.)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_isoneohooke:
            {
              MAT::ELASTIC::PAR::IsoNeoHooke* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoNeoHooke*>(actelastmat->Parameter());
              params2->SetMue(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_isoyeoh:
            {
              MAT::ELASTIC::PAR::IsoYeoh* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
              params2->SetC1(abs(p_cur(j)));
              params2->SetC2(abs(p_cur(j+1)));
              params2->SetC3(abs(p_cur(j+2)));
              j = j+3;
              break;
            }
            case INPAR::MAT::mes_iso1pow:
            {
              MAT::ELASTIC::PAR::Iso1Pow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_iso2pow:
            {
              MAT::ELASTIC::PAR::Iso2Pow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::Iso2Pow*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_coup1pow:
            {
              MAT::ELASTIC::PAR::Coup1Pow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::Coup1Pow*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_coup2pow:
            {
              MAT::ELASTIC::PAR::Coup2Pow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::Coup2Pow*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_coupmooneyrivlin:
            {
              MAT::ELASTIC::PAR::CoupMooneyRivlin* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupMooneyRivlin*>(actelastmat->Parameter());
              params2->SetC1(abs(p_cur(j)));
              params2->SetC2(abs(p_cur(j+1)));
              params2->SetC3(abs(p_cur(j+2)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_isoexpopow:
            {
              MAT::ELASTIC::PAR::IsoExpoPow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoExpoPow*>(actelastmat->Parameter());
              params2->SetK1(abs(p_cur(j)));
              params2->SetK2(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_isomooneyrivlin:
            {
              MAT::ELASTIC::PAR::IsoMooneyRivlin* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoMooneyRivlin*>(actelastmat->Parameter());
              params2->SetC1(abs(p_cur(j)));
              params2->SetC2(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_volsussmanbathe:
            {
              MAT::ELASTIC::PAR::VolSussmanBathe* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::VolSussmanBathe*>(actelastmat->Parameter());
              params2->SetKappa(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_volpenalty:
            {
              MAT::ELASTIC::PAR::VolPenalty* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::VolPenalty*>(actelastmat->Parameter());
              params2->SetEpsilon(abs(p_cur(j)));
              params2->SetGamma(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_vologden:
            {
              MAT::ELASTIC::PAR::VolOgden* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
              params2->SetKappa(abs(p_cur(j)));
              //params2->SetBeta(abs(p_cur(j+1)));
              j = j+1;
              break;
            }
            default:
              dserror("cannot deal with this material");
              break;
            }
          }
        }
      }
      if (material->Type() == INPAR::MAT::m_viscoelasthyper)
      {
        MAT::PAR::ViscoElastHyper* params = dynamic_cast<MAT::PAR::ViscoElastHyper*>(material->Parameter());
        const int nummat               = params->nummat_;
        const std::vector<int>* matids = params->matids_;

        // For each of the summands of the viscohyperelastic material we need to add the
        // parameters to the inverse analysis

        // Problems with beta, is it the only negative parameter? Maybe
        // we should exclude it

        //itterator to go through the parameters
        int j = 0;

        for (int i=0; i<nummat; i++)
        {
          //get the material of the summand
          const int id = (*matids)[i];

          if (myehmatset.size()==0 or myehmatset.find(id)!=myehmatset.end())
          {
            const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;

            switch (actelastmat->Type())
            {
            case INPAR::MAT::mes_couplogneohooke:
            {
              MAT::ELASTIC::PAR::CoupLogNeoHooke* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupLogNeoHooke*>(actelastmat->Parameter());
              params2->SetMue(abs(p_cur(j)));
              params2->SetLambda(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_coupexppol:
            {
              MAT::ELASTIC::PAR::CoupExpPol* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupExpPol*>(actelastmat->Parameter());
              params2->SetA(abs(p_cur(j)));
              params2->SetB(abs(p_cur(j+1)));
              params2->SetC(abs(p_cur(j+2)));
              j = j+3;
              break;
            }
            case INPAR::MAT::mes_coupneohooke:
            {
              MAT::ELASTIC::PAR::CoupNeoHooke* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupNeoHooke*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              params2->SetBeta(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_coupblatzko:
            {
              MAT::ELASTIC::PAR::CoupBlatzKo* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupBlatzKo*>(actelastmat->Parameter());
              params2->SetMue(abs(p_cur(j)));
              //params2->SetF(abs(p_cur(j+1)));
              params2->SetNue((abs(p_cur(j+1)))/(2.*(abs(p_cur(j+1))+1.)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_isoneohooke:
            {
              MAT::ELASTIC::PAR::IsoNeoHooke* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoNeoHooke*>(actelastmat->Parameter());
              params2->SetMue(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_isoyeoh:
            {
              MAT::ELASTIC::PAR::IsoYeoh* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
              params2->SetC1(abs(p_cur(j)));
              params2->SetC2(abs(p_cur(j+1)));
              params2->SetC3(abs(p_cur(j+2)));
              j = j+3;
              break;
            }
            case INPAR::MAT::mes_iso1pow:
            {
              MAT::ELASTIC::PAR::Iso1Pow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_iso2pow:
            {
              MAT::ELASTIC::PAR::Iso2Pow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::Iso2Pow*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_coup1pow:
            {
              MAT::ELASTIC::PAR::Coup1Pow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::Coup1Pow*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_coup2pow:
            {
              MAT::ELASTIC::PAR::Coup2Pow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::Coup2Pow*>(actelastmat->Parameter());
              params2->SetC(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_coupmooneyrivlin:
            {
              MAT::ELASTIC::PAR::CoupMooneyRivlin* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::CoupMooneyRivlin*>(actelastmat->Parameter());
              params2->SetC1(abs(p_cur(j)));
              params2->SetC2(abs(p_cur(j+1)));
              params2->SetC3(abs(p_cur(j+2)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_isoexpopow:
            {
              MAT::ELASTIC::PAR::IsoExpoPow* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoExpoPow*>(actelastmat->Parameter());
              params2->SetK1(abs(p_cur(j)));
              params2->SetK2(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_isomooneyrivlin:
            {
              MAT::ELASTIC::PAR::IsoMooneyRivlin* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoMooneyRivlin*>(actelastmat->Parameter());
              params2->SetC1(abs(p_cur(j)));
              params2->SetC2(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_volsussmanbathe:
            {
              MAT::ELASTIC::PAR::VolSussmanBathe* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::VolSussmanBathe*>(actelastmat->Parameter());
              params2->SetKappa(abs(p_cur(j)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_volpenalty:
            {
              MAT::ELASTIC::PAR::VolPenalty* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::VolPenalty*>(actelastmat->Parameter());
              params2->SetEpsilon(abs(p_cur(j)));
              params2->SetGamma(abs(p_cur(j+1)));
              j = j+2;
              break;
            }
            case INPAR::MAT::mes_vologden:
            {
              MAT::ELASTIC::PAR::VolOgden* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
              params2->SetKappa(abs(p_cur(j)));
              //params2->SetBeta(abs(p_cur(j+1)));
              j = j+1;
              break;
            }
            case INPAR::MAT::mes_isoratedep:
            {
              MAT::ELASTIC::PAR::IsoRateDep* params2 =
                dynamic_cast<MAT::ELASTIC::PAR::IsoRateDep*>(actelastmat->Parameter());
              params2->SetN(abs(p_cur(j)));
              j = j+1;
              break;
            }
            default:
              dserror("cannot deal with this material");
            break;
            }
          }
        }
      }
    }
  }

  return;
}


void STR::InvAnalysis::MultiInvAnaInit()
{
  for (int i=0; i<discret_->NumMyColElements(); i++)
  {
    DRT::Element* actele = discret_->lColElement(i);
    Teuchos::RCP<MAT::Material> mat = actele->Material();
    if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
      bool eleowner = false;
      if (discret_->Comm().MyPID()==actele->Owner()) eleowner = true;
      micro->InvAnaInit(eleowner,actele->Id());
    }
  }
}

