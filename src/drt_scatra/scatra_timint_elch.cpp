/*!-----------------------------------------------------------------------------------------------*
\file scatra_timint_elch.cpp

  \brief scatra time integration for elch

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 *------------------------------------------------------------------------------------------------*/

#include "scatra_timint_elch.H"
#include "scatra_utils_splitstrategy.H" // for blockmatrix-splitstrategy

#include "../drt_inpar/inpar_elch.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/ion.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "../drt_fluid/fluid_utils.H" // for splitter
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_krylov_projector.H"

#include "../drt_fluid/fluid_meshtying.H"


/*----------------------------------------------------------------------*
 | constructor                                              ehrl  01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElch::ScaTraTimIntElch(
        Teuchos::RCP<DRT::Discretization>        dis,
        Teuchos::RCP<LINALG::Solver>             solver,
        Teuchos::RCP<Teuchos::ParameterList>     params,
        Teuchos::RCP<Teuchos::ParameterList>     sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList>     extraparams,
        Teuchos::RCP<IO::DiscretizationWriter>   output)
  : ScaTraTimIntImpl(dis,solver,sctratimintparams,extraparams,output),
    elchparams_     (params),
    elchtype_       (DRT::INPUT::IntegralValue<INPAR::ELCH::ElchType>(*elchparams_,"ELCHTYPE")),
    equpot_         (DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(*elchparams_,"EQUPOT")),
    frt_            (0.0),
    elchdensnm_     (Teuchos::null),
    elchdensn_      (Teuchos::null),
    elchdensnp_     (Teuchos::null),
    gstatnumite_    (0),
    gstatincrement_ (0.0),
    c0_             (0,0.0),
    sigma_          (Teuchos::null),
    densific_       (0),
    dlcapexists_    (false),
    ektoggle_       (Teuchos::null),
    dctoggle_       (Teuchos::null)
{
  return;
}


/*----------------------------------------------------------------------*
 | initialize algorithm                                 rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::Init()
{
  if(elchtype_==INPAR::ELCH::elchtype_undefined)
    dserror("The ElchType is not specified in the input file!");

  // initialize time-dependent electrode kinetics variables (galvanostatic mode or double layer contribution)
  ComputeTimeDerivPot0(true);

  // Important: this adds the required ConditionID's to the single conditions.
  // It is necessary to do this BEFORE ReadRestart() is called!
  // Output to screen and file is suppressed
  OutputElectrodeInfo(false,false);

  // initializes variables for natural convection (ELCH) if necessary
  densific_.resize(numscal_);
  SetupElchNatConv();

  // initialize dirichlet toggle:
  // for certain ELCH problem formulations we have to provide
  // additional flux terms / currents across Dirichlet boundaries for the standard element call
  Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()),false);
  dirichones->PutScalar(1.0);
  dctoggle_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);
  dbcmaps_->InsertCondVector(dirichones, dctoggle_);

  // screen output (has to come after SetInitialField)
  // a safety check for the solver type
  if ((numscal_ > 1) && (solvtype_!=INPAR::SCATRA::solvertype_nonlinear))
    dserror("Solver type has to be set to >>nonlinear<< for ion transport.");

  // check validity of material and element formulation
  // important for porous Diffusion-Conduction formulation using material ElchMat
  {
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::check_scatra_element_parameter);
    eleparams.set<int>("scatratype",scatratype_);

    // call standard loop over elements
    discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  }

  frt_ = INPAR::ELCH::faraday_const/(INPAR::ELCH::gas_const * elchparams_->get<double>("TEMPERATURE"));

  if (myrank_==0)
  {
    std::cout<<"\nSetup of splitter: numscal = "<<numscal_<<std::endl;
    std::cout<<"Temperature value T (Kelvin)     = "<<elchparams_->get<double>("TEMPERATURE")<<std::endl;
    std::cout<<"Constant F/RT                    = "<<frt_<<std::endl;
  }

  sigma_= Teuchos::rcp(new Epetra_SerialDenseVector(numdofpernode_));
  // conductivity must be stored for the galvanostatic condition in a global variable
  ComputeConductivity(); // every processor has to do this call
  if (myrank_==0)
  {
    for (int k=0;k < numscal_;k++)
    {
      std::cout<<"Electrolyte conductivity (species "<<k+1<<")    = "<<(*sigma_)[k]<<std::endl;
    }
    if (equpot_==INPAR::ELCH::equpot_enc_pde_elim)
    {
      double diff = (*sigma_)[0];
      for (int k=1;k < numscal_;k++)
      {
        diff += (*sigma_)[k];
      }
      std::cout<<"Electrolyte conductivity (species elim) = "<<(*sigma_)[numscal_]-diff<<std::endl;
    }
    std::cout<<"Electrolyte conductivity (all species)  = "<<(*sigma_)[numscal_]<<std::endl<<std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 | initialization of system matrix                           ehrl 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::InitSystemMatrix()
{
  if (numscal_ > 1) // we have at least two ion species + el. potential
  {
    // number of concentrations transported is numdof-1
    numscal_ -= 1;

    Teuchos::ParameterList& diffcondparams = elchparams_->sublist("DIFFCOND");

    // currrent is a solution variable
    if(DRT::INPUT::IntegralValue<int>(diffcondparams,"CURRENT_SOLUTION_VAR"))
    {
      // shape of local row element(0) -> number of space dimensions
      //int dim = DRT::Problem::Instance()->NDim();
      int dim = DRT::UTILS::getDimension(discret_->lRowElement(0)->Shape());
      // number of concentrations transported is numdof-1-nsd
      numscal_ -= dim;
    }

    // The diffusion-conduction formulation does not support all options of the Nernst-Plack formulation
    // Let's check for valid options
    if(elchtype_==INPAR::ELCH::elchtype_diffcond)
      ValidParameterDiffCond();
  }
  // set up the concentration-el.potential splitter
  splitter_ = Teuchos::rcp(new LINALG::MapExtractor);
  FLD::UTILS::SetupFluidSplit(*discret_,numscal_,*splitter_);

  if (DRT::INPUT::IntegralValue<int>(*params_,"BLOCKPRECOND")
      and msht_ == INPAR::FLUID::no_meshtying)
  {
    if ((equpot_!=INPAR::ELCH::equpot_enc))
      dserror("Special ELCH assemble strategy for block-matrix will not assemble A_11 block!");

    // initial guess for non-zeros per row: 27 neighboring nodes for hex8
    // this is enough! A higher guess would require too much memory!
    // A more precise guess for every submatrix would read:
    // A_00: 27*1,  A_01: 27*1,  A_10: 27*numscal_ due to electroneutrality, A_11: EMPTY matrix !!!!!
    // usage of a split strategy that makes use of the ELCH-specific sparsity pattern
    Teuchos::RCP<LINALG::BlockSparseMatrix<SCATRA::SplitStrategy> > blocksysmat =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<SCATRA::SplitStrategy>(*splitter_,*splitter_,27,false,true));

    blocksysmat->SetNumScal(numscal_);

    sysmat_ = blocksysmat;
  }
  else
    ScaTraTimIntImpl::InitSystemMatrix();

  return;
}


/*----------------------------------------------------------------------*
 | setup meshtying system                                    ehrl 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::SetupMeshtying()
{
  // Important:
  // Meshtying in elch is not tested well!!

  bool onlypotential = (DRT::INPUT::IntegralValue<int>(*params_,"ONLYPOTENTIAL"));

  if(msht_== INPAR::FLUID::condensed_bmat)
    dserror("The 2x2 block solver algorithm, which is necessary for a block matrix system,\n"
            "is not integrated into the adapter_scatra_base_algorithm. Just do it!!");

  if ((msht_== INPAR::FLUID::condensed_bmat or
      msht_== INPAR::FLUID::condensed_bmat_merged) and
      (equpot_== INPAR::ELCH::equpot_enc))
    dserror("In the context of mesh-tying, the ion-transport system inluding the electroneutrality condition \n"
        "cannot be solved in a block matrix");

  // meshtying: all dofs (all scalars + potential) are coupled
  int numdof = numscal_+1;

  // define coupling
  std::vector<int> coupleddof(numdof, 1);
  if(onlypotential==true)
  {
    // meshtying: only potential is coupled
    // coupleddof = [0, 0, ..., 0, 1]
    for(int ii=0;ii<numscal_;++ii)
      coupleddof[ii]=0;
  }

  meshtying_ = Teuchos::rcp(new FLD::Meshtying(discret_, *solver_, msht_, DRT::Problem::Instance()->NDim()));
  sysmat_ = meshtying_->Setup(coupleddof);

  return;
} // ScaTraTimIntImpl::SetupMeshtying()


/*----------------------------------------------------------------------*
 | set all general parameters for element                    ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::SetElementGeneralScaTraParameter()
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action",SCATRA::set_elch_scatra_parameter);

  // set type of scalar transport problem
  eleparams.set<int>("scatratype",scatratype_);

  // general scatra parameter
  eleparams.set<int>("form of convective term",convform_);
  eleparams.set("isale",isale_);

  // set flag for writing the flux vector fields
  eleparams.set<int>("writeflux",writeflux_);

  //! set vector containing ids of scalars for which flux vectors are calculated
  eleparams.set<Teuchos::RCP<std::vector<int> > >("writefluxids",writefluxids_);

  // parameters for stabilization
  eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

  // general elch parameter
  eleparams.set<double>("frt",INPAR::ELCH::faraday_const/(INPAR::ELCH::gas_const*(elchparams_->get<double>("TEMPERATURE"))));
  eleparams.set<int>("elchtype",elchtype_);
  eleparams.set<int>("equpot",equpot_);

  // parameters for diffusion-conduction formulation
  eleparams.sublist("DIFFCOND") = elchparams_->sublist("DIFFCOND");

  // call standard loop over elements
  discret_->Evaluate(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | add parameters depending on the problem              rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::AddProblemSpecificParametersAndVectors(
  Teuchos::ParameterList& params //!< parameter list
)
{
  // ELCH specific factor F/RT
  params.set("frt",frt_);

  // parameters for Elch/DiffCond formulation
  params.sublist("DIFFCOND") = elchparams_->sublist("DIFFCOND");

  discret_->SetState("dctoggle",dctoggle_);

  return;
}

/*--------------------------------------------------------------------------*
 | add parameters depending on the problem for inital phidt rasthofer 12/13 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::AddProblemSpecificParametersAndVectorsForCalcInitialPhiDt(
  Teuchos::ParameterList& params //!< parameter list
)
{
  discret_->SetState("dctoggle",dctoggle_);

  return;
}


/*----------------------------------------------------------------------*
 | contains the elch-specific nonlinear iteration loop       ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::NonlinearSolve()
{
  bool stopgalvanostat(false);
  gstatnumite_=1;

  // galvanostatic control (ELCH)
  while (!stopgalvanostat)
  {
    ScaTraTimIntImpl::NonlinearSolve();

    stopgalvanostat = ApplyGalvanostaticControl();
  }  // end galvanostatic control

  return;
}


/*----------------------------------------------------------------------*
 | Calculate problem specific norm                            ehrl 01/14|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::CalcProblemSpecificNorm(
    double& conresnorm,
    double& incconnorm_L2,
    double& connorm_L2,
    double& incpotnorm_L2,
    double& potnorm_L2,
    double& potresnorm,
    double& conresnorminf)
{
  Teuchos::RCP<Epetra_Vector> onlycon = splitter_->ExtractOtherVector(residual_);
  onlycon->Norm2(&conresnorm);

  splitter_->ExtractOtherVector(increment_,onlycon);
  onlycon->Norm2(&incconnorm_L2);

  splitter_->ExtractOtherVector(phinp_,onlycon);
  onlycon->Norm2(&connorm_L2);

  Teuchos::RCP<Epetra_Vector> onlypot = splitter_->ExtractCondVector(residual_);
  onlypot->Norm2(&potresnorm);

  splitter_->ExtractCondVector(increment_,onlypot);
  onlypot->Norm2(&incpotnorm_L2);

  splitter_->ExtractCondVector(phinp_,onlypot);
  onlypot->Norm2(&potnorm_L2);

  return;
}


/*----------------------------------------------------------------------*
 |  calculate error compared to analytical solution            gjb 10/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::EvaluateErrorComparedToAnalyticalSol()
{
  const INPAR::SCATRA::CalcError calcerr
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::CalcError>(*params_,"CALCERROR");

  switch (calcerr)
  {
  case INPAR::SCATRA::calcerror_no: // do nothing (the usual case)
    break;
  case INPAR::SCATRA::calcerror_Kwok_Wu:
  {
    //   References:

    //   Kwok, Yue-Kuen and Wu, Charles C. K.
    //   "Fractional step algorithm for solving a multi-dimensional
    //   diffusion-migration equation"
    //   Numerical Methods for Partial Differential Equations
    //   1995, Vol 11, 389-397

    //   G. Bauer, V. Gravemeier, W.A. Wall, A 3D finite element approach for the coupled
    //   numerical simulation of electrochemical systems and fluid flow,
    //   International Journal for Numerical Methods in Engineering, 86
    //   (2011) 1339–1359. DOI: 10.1002/nme.3107

    // create the parameters for the error calculation
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::calc_error);
    eleparams.set<int>("scatratype",scatratype_);
    eleparams.set("total time",time_);
    eleparams.set("frt",frt_);
    eleparams.set<int>("calcerrorflag",calcerr);
    //provide displacement field in case of ALE
    eleparams.set("isale",isale_);
    if (isale_)
      discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);
    // parameters for Elch/DiffCond formulation
    eleparams.sublist("DIFFCOND") = elchparams_->sublist("DIFFCOND");

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);

    // get (squared) error values
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(eleparams, errors);
    discret_->ClearState();

    double conerr1 = 0.0;
    double conerr2 = 0.0;
    // for the L2 norm, we need the square root
    if(numscal_==2)
    {
      conerr1 = sqrt((*errors)[0]);
      conerr2 = sqrt((*errors)[1]);
    }
    else if(numscal_==1)
    {
      conerr1 = sqrt((*errors)[0]);
      conerr2 = 0.0;
    }
    else
      dserror("The analytical solution of Kwok and Wu is only defined for two species");

    double poterr  = sqrt((*errors)[2]);

    if (myrank_ == 0)
    {
      printf("\nL2_err for Kwok and Wu (time = %f):\n", time_);
      printf(" concentration1 %15.8e\n concentration2 %15.8e\n potential      %15.8e\n\n",
             conerr1,conerr2,poterr);
    }
#if 0
    if (myrank_ == 0)
    {
      // append error of the last time step to the error file
      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        ostd::stringstream temp;
        const std::string simulation = DRT::Problem::Instance()->OutputControlFile()->FileName();
        //const std::string fname = simulation+".relerror";
        const std::string fname = "XXX_kwok_xele.relerror";

        double elelength=0.0;
        if(simulation.find("5x5")!=std::string::npos)
          elelength=0.2;
        else if(simulation.find("10x10")!=std::string::npos)
          elelength=0.1;
        else if(simulation.find("20x20")!=std::string::npos)
          elelength=0.05;
        else if(simulation.find("40x40")!=std::string::npos)
          elelength=0.025;
        else if(simulation.find("50x50")!=std::string::npos)
          elelength=0.02;
        else if(simulation.find("80x80")!=std::string::npos)
          elelength=0.0125;
        else std::cout << "Warning: file name did not allow a evaluation of the element size!!!" << std::endl;

        std::ofstream f;
        f.precision(12);
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        //f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
        f << "#| Step | Time | c1 abs. error L2 | c2 abs. error L2 | phi abs. error L2 |  element length  |\n";
        //f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f << step_ << " " << time_ << " " << conerr1 << " " << conerr2 << " " <<  poterr << " "<< elelength << "" << "\n";
        f.flush();
        f.close();
      }
    }
#endif
  }
  break;
  case INPAR::SCATRA::calcerror_cylinder:
  {
    //   Reference:
    //   G. Bauer, V. Gravemeier, W.A. Wall, A 3D finite element approach for the coupled
    //   numerical simulation of electrochemical systems and fluid flow,
    //   International Journal for Numerical Methods in Engineering, 2011

    // create the parameters for the error calculation
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::calc_error);
    eleparams.set<int>("scatratype",scatratype_);
    eleparams.set("total time",time_);
    eleparams.set("frt",frt_);
    eleparams.set<int>("calcerrorflag",calcerr);
    //provide displacement field in case of ALE
    eleparams.set("isale",isale_);
    if (isale_)
      discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);

    // get (squared) error values
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(3));
    discret_->EvaluateScalars(eleparams, errors);
    discret_->ClearState();

    // for the L2 norm, we need the square root
    double conerr1 = sqrt((*errors)[0]);
    double conerr2 = sqrt((*errors)[1]);
    double poterr  = sqrt((*errors)[2]);

    if (myrank_ == 0)
    {
      printf("\nL2_err for concentric cylinders (time = %f):\n", time_);
      printf(" concentration1 %15.8e\n concentration2 %15.8e\n potential      %15.8e\n\n",
             conerr1,conerr2,poterr);
    }
  }
  break;
  case INPAR::SCATRA::calcerror_electroneutrality:
  {
    // compute L2 norm of electroneutrality condition

    // create the parameters for the error calculation
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::calc_error);
    eleparams.set<int>("scatratype",scatratype_);
    eleparams.set("total time",time_);
    eleparams.set("frt",frt_);
    eleparams.set<int>("calcerrorflag",calcerr);
    //provide displacement field in case of ALE
    eleparams.set("isale",isale_);
    if (isale_)
      discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phinp",phinp_);

    // get (squared) error values
    Teuchos::RCP<Epetra_SerialDenseVector> errors
      = Teuchos::rcp(new Epetra_SerialDenseVector(1));
    discret_->EvaluateScalars(eleparams, errors);
    discret_->ClearState();

    // for the L2 norm, we need the square root
    double err = sqrt((*errors)[0]);

    if (myrank_ == 0)
    {
      printf("\nL2_err for electroneutrality (time = %f):\n", time_);
      printf(" Deviation from ENC: %15.8e\n\n",err);
    }
  }
  break;
  default:
    dserror("Cannot calculate error. Unknown type of analytical test problem"); break;
  }
  return;
} // SCATRA::ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol


/*==========================================================================*/
// ELCH
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | compute density from ion concentrations                    gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::ComputeDensity()
{
  double newdensity(0.0);
  int err(0);

  // loop over all local nodes
  for(int lnodeid=0; lnodeid<discret_->NumMyRowNodes(); lnodeid++)
  {
    // get the processor's local node
    DRT::Node* lnode = discret_->lRowNode(lnodeid);

    // get the degrees of freedom associated with this node
    std::vector<int> nodedofs;
    nodedofs = discret_->Dof(lnode);
    int numdof = nodedofs.size();

    newdensity= 1.0;
    // loop over all ionic species
    for(int k=0; k<numscal_; k++)
    {
      /*
        //                  k=numscal_-1
        //          /       ----                         \
        //         |        \                            |
        // rho_0 * | 1 +    /       alfa_k * (c_k - c_0) |
        //         |        ----                         |
        //          \       k=0                          /
        //
        // For use of molar mass M_k:  alfa_k = M_k/rho_0  !!
       */

      // global and processor's local DOF ID
      const int globaldofid = nodedofs[k];
      const int localdofid = phinp_->Map().LID(globaldofid);
      if (localdofid < 0)
        dserror("localdofid not found in map for given globaldofid");

      // compute contribution to density due to ionic species k
      newdensity += densific_[k]*((*phinp_)[localdofid]-c0_[k]);
    }

    // insert the current density value for this node
    // (has to be at position of el potential/ the position of the last dof!
    const int globaldofid = nodedofs[numdof-1];
    const int localdofid = phinp_->Map().LID(globaldofid);
    if (localdofid < 0)
      dserror("localdofid not found in map for given globaldofid");

    err = elchdensnp_->ReplaceMyValue(localdofid,0,newdensity);

    if (err != 0) dserror("error while inserting a value into elchdensnp_");
  } // loop over all local nodes
  return;
} // SCATRA::ScaTraTimIntImpl::ComputeDensity


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::Update(const int num)
{
  // perform update of time-dependent electrode variables
  ElectrodeKineticsTimeUpdate();

  return;
}


/*----------------------------------------------------------------------*
 |  output of electrode status information to screen         gjb  01/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP< std::vector<double> > SCATRA::ScaTraTimIntElch::OutputElectrodeInfo(
    bool printtoscreen,
    bool printtofile)
{
  // results of the first electrode kinetics BC
  Teuchos::RCP< std::vector<double> > firstelectkin = Teuchos::rcp(new std::vector<double>(3, 0.0));

  // evaluate the following type of boundary conditions:
  std::string condname("ElectrodeKinetics");
  std::vector<DRT::Condition*> cond;
  discret_->GetCondition(condname,cond);

  // leave method, if there's nothing to do!
  if (!cond.size()) return firstelectkin;

  double sum(0.0);

  if ((myrank_ == 0) and printtoscreen)
  {
    std::cout<<"Status of '"<<condname<<"':\n"
    <<"++----+---------------------+------------------+----------------------+--------------------+----------------+----------------+"<<std::endl;
    printf("|| ID |    Total current    | Area of boundary | Mean current density | Mean overpotential | Electrode pot. | Mean Concentr. |\n");
  }

  // first, add to all conditions of interest a ConditionID
  for (int condid = 0; condid < (int) cond.size(); condid++)
  {
    // is there already a ConditionID?
    const std::vector<int>*    CondIDVec  = cond[condid]->Get<std::vector<int> >("ConditionID");
    if (CondIDVec)
    {
      if ((*CondIDVec)[0] != condid)
        dserror("Condition %s has non-matching ConditionID",condname.c_str());
    }
    else
    {
      // let's add a ConditionID
      cond[condid]->Add("ConditionID",condid);
    }
  }
  // now we evaluate the conditions and separate via ConditionID
  for (int condid = 0; condid < (int) cond.size(); condid++)
  {
    double currtangent(0.0); // this value remains unused here!
    double currresidual(0.0); // this value remains unused here!
    double electrodesurface(0.0); // this value remains unused here!
    double electrodepot(0.0); // this value remains unused here!
    double meanoverpot(0.0); // this value remains unused here!

    Teuchos::RCP< std::vector<double> > electkin = OutputSingleElectrodeInfo(
        cond[condid],
        condid,
        printtoscreen,
        printtofile,
        sum,
        currtangent,
        currresidual,
        electrodesurface,
        electrodepot,
        meanoverpot);

    // only the first condition is tested
    if(condid==0)
      firstelectkin = electkin;
  } // loop over condid

  if ((myrank_==0) and printtoscreen)
  {
    std::cout<<"++----+---------------------+------------------+----------------------+--------------------+----------------+----------------+"<<std::endl;
    // print out the net total current for all indicated boundaries
    printf("Net total current over boundary: %10.3E\n\n",sum);
  }

  // clean up
  discret_->ClearState();

  return firstelectkin;
} // ScaTraImplicitTimeInt::OutputElectrodeInfo


/*----------------------------------------------------------------------*
 |  get electrode status for single boundary condition       gjb  11/09 |
 *----------------------------------------------------------------------*/
Teuchos::RCP< std::vector<double> > SCATRA::ScaTraTimIntElch::OutputSingleElectrodeInfo(
    DRT::Condition* condition,
    const int condid,
    const bool printtoscreen,
    const bool printtofile,
    double& currentsum,
    double& currtangent,
    double& currresidual,
    double& electrodesurface,
    double& electrodepot,
    double& meanoverpot)
{
  // return vector for nightly test cases
  Teuchos::RCP< std::vector<double> > singleelectkin = Teuchos::rcp(new std::vector<double>(3, 0.0));

  // safety check: is there already a ConditionID?
  const std::vector<int>* CondIDVec  = condition->Get<std::vector<int> >("ConditionID");
  if (not CondIDVec) dserror("Condition has not yet a ConditionID");

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  // needed for double-layer capacity!
  discret_->SetState("phidtnp",phidtnp_);

  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::bd_calc_elch_electrode_kinetics);
  eleparams.set<int>("scatratype",scatratype_);
  eleparams.set("calc_status",true); // just want to have a status ouput!

  // parameters for Elch/DiffCond formulation
  eleparams.sublist("DIFFCOND") = elchparams_->sublist("DIFFCOND");

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_)
    discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // Since we just want to have the status ouput for t_{n+1},
  // we have to take care for Gen.Alpha!
  // AddTimeIntegrationSpecificVectors cannot be used since we do not want
  // an evaluation for t_{n+\alpha_f} !!!

  // Warning:
  // Specific time integration parameter are set in the following function.
  // In the case of a genalpha-time integration scheme the solution vector phiaf_ at time n+af
  // is passed to the element evaluation routine. Therefore, the electrode status is evaluate at a
  // different time (n+af) than our output routine (n+1), resulting in slightly different values at the electrode.
  // A different approach is not possible (without major hacks) since the time-integration scheme is
  // necessary to perform galvanostatic simulations, for instance.
  // Think about: double layer effects for genalpha time-integratio scheme

  // add element parameters according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  // values to be computed
  eleparams.set("currentintegral",0.0);
  eleparams.set("currentdlintegral",0.0);
  eleparams.set("boundaryintegral",0.0);
  eleparams.set("electpotentialintegral",0.0);
  eleparams.set("overpotentialintegral",0.0);
  eleparams.set("electrodedifferencepotentialintegral",0.0);
  eleparams.set("opencircuitpotentialintegral",0.0);
  eleparams.set("concentrationintegral",0.0);
  eleparams.set("currentderiv",0.0);
  eleparams.set("currentderivDL",0.0);
  eleparams.set("currentresidual",0.0);

  // would be nice to have a EvaluateScalar for conditions!!!
  discret_->EvaluateCondition(eleparams,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,"ElectrodeKinetics",condid);

  // get integral of current on this proc
  double currentintegral = eleparams.get<double>("currentintegral");
  // get integral of current on this proc
  double currentdlintegral = eleparams.get<double>("currentdlintegral");
  // get area of the boundary on this proc
  double boundaryint = eleparams.get<double>("boundaryintegral");
  // get integral of overpotential on this proc
  double electpotentialint = eleparams.get<double>("electpotentialintegral");
  // get integral of overpotential on this proc
  double overpotentialint = eleparams.get<double>("overpotentialintegral");
  // get integral of overpotential on this proc
  double epdint = eleparams.get<double>("electrodedifferencepotentialintegral");
  // get integral of overpotential on this proc
  double ocpint = eleparams.get<double>("opencircuitpotentialintegral");
  // get integral of reactant concentration on this proc
  double cint = eleparams.get<double>("concentrationintegral");
  // tangent of current w.r.t. electrode potential on this proc
  double currderiv = eleparams.get<double>("currentderiv");
  // tangent of current w.r.t. electrode potential on this proc
  double currderivDL = eleparams.get<double>("currentderivDL");
  // get negative current residual (rhs of galvanostatic balance equation)
  double currentresidual = eleparams.get<double>("currentresidual");

  // care for the parallel case
  double parcurrentintegral = 0.0;
  discret_->Comm().SumAll(&currentintegral,&parcurrentintegral,1);
  double parcurrentdlintegral = 0.0;
  discret_->Comm().SumAll(&currentdlintegral,&parcurrentdlintegral,1);
  double parboundaryint = 0.0;
  discret_->Comm().SumAll(&boundaryint,&parboundaryint,1);
  double parelectpotentialint = 0.0;
  discret_->Comm().SumAll(&electpotentialint,&parelectpotentialint,1);
  double paroverpotentialint = 0.0;
  discret_->Comm().SumAll(&overpotentialint,&paroverpotentialint,1);
  double parepdint = 0.0;
  discret_->Comm().SumAll(&epdint,&parepdint,1);
  double parocpint = 0.0;
  discret_->Comm().SumAll(&ocpint,&parocpint,1);
  double parcint = 0.0;
  discret_->Comm().SumAll(&cint,&parcint,1);
  double parcurrderiv = 0.0;
  discret_->Comm().SumAll(&currderiv,&parcurrderiv ,1);
  double parcurrderivDL = 0.0;
  discret_->Comm().SumAll(&currderivDL,&parcurrderivDL ,1);
  double parcurrentresidual = 0.0;
  discret_->Comm().SumAll(&currentresidual,&parcurrentresidual ,1);

  // specify some return values
  currentsum += parcurrentintegral; // sum of currents
  currtangent  = parcurrderiv;      // tangent w.r.t. electrode potential on metal side
  currresidual = parcurrentresidual;
  electrodesurface = parboundaryint;
  electrodepot = parelectpotentialint/parboundaryint;
  meanoverpot = paroverpotentialint/parboundaryint;

  // clean up
  discret_->ClearState();

  // print out results to screen/file if desired
  if (myrank_ == 0)
  {
    if (printtoscreen) // print out results to screen
    {
      printf("|| %2d |     %10.3E      |    %10.3E    |      %10.3E      |     %10.3E     |   %10.3E   |   %10.3E   |\n",
          condid,parcurrentintegral+currentdlintegral,parboundaryint,parcurrentintegral/parboundaryint+currentdlintegral/parboundaryint,paroverpotentialint/parboundaryint, electrodepot, parcint/parboundaryint);
    }

    if (printtofile)// write results to file
    {
      std::ostringstream temp;
      temp << condid;
      const std::string fname
      = DRT::Problem::Instance()->OutputControlFile()->FileName()+".electrode_status_"+temp.str()+".txt";

      std::ofstream f;
      if (Step() <= 1)
      {
        f.open(fname.c_str(),std::fstream::trunc);
        f << "#ID,Step,Time,Total_current,Area_of_boundary,Mean_current_density_electrode_kinetics,Mean_current_density_dl,Mean_overpotential,Mean_electrode_pot_diff,Mean_opencircuit_pot,Electrode_pot,Mean_concentration\n";
      }
      else
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);

      f << condid << "," << Step() << "," << Time() << "," << parcurrentintegral+currentdlintegral  << "," << parboundaryint
      << "," << parcurrentintegral/parboundaryint << "," << parcurrentdlintegral/parboundaryint
      << "," << paroverpotentialint/parboundaryint << "," <<
      parepdint/parboundaryint << "," << parocpint/parboundaryint << "," << electrodepot << ","
      << parcint/parboundaryint << "\n";
      f.flush();
      f.close();
    }

    (*singleelectkin)[0]= parcint/parboundaryint;
    (*singleelectkin)[1]= paroverpotentialint/parboundaryint;
  } // if (myrank_ == 0)

  // galvanostatic simulations:
  // add the double layer current to the Butler-Volmer current
  currentsum += parcurrentdlintegral;
  (*singleelectkin)[2]=currentsum;

  return singleelectkin;
} // SCATRA::ScaTraTimIntImpl::OutputSingleElectrodeInfo


/*----------------------------------------------------------------------*
 | perform setup of natural convection applications (ELCH)    gjb 07/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::SetupElchNatConv()
{
  // loads densification coefficients and the initial mean concentration

  // only required for ELCH with natural convection
  if (DRT::INPUT::IntegralValue<int>(*elchparams_,"NATURAL_CONVECTION") == true)
  {
    // allocate denselch_ with *dofrowmap and initialize it
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    elchdensnp_ = LINALG::CreateVector(*dofrowmap,true);
    elchdensnp_->PutScalar(1.0);

    // Calculate the initial mean concentration value
    if (numscal_ < 1) dserror("Error since numscal = %d. Not allowed since < 1",numscal_);
    c0_.resize(numscal_);

    discret_->ClearState();
    discret_->SetState("phinp",phinp_);
    // set action for elements
    Teuchos::ParameterList eleparams;
    eleparams.set<int>("action",SCATRA::calc_mean_scalars);
    eleparams.set<int>("scatratype",scatratype_);
    eleparams.set("inverting",false);

    //provide displacement field in case of ALE
    if (isale_)
      discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

    // evaluate integrals of concentrations and domain
    Teuchos::RCP<Epetra_SerialDenseVector> scalars
    = Teuchos::rcp(new Epetra_SerialDenseVector(numscal_+1));
    discret_->EvaluateScalars(eleparams, scalars);
    discret_->ClearState();   // clean up

    // calculate mean_concentration
    const double domint  = (*scalars)[numscal_];
    for(int k=0;k<numscal_;k++)
    {
      c0_[k] = (*scalars)[k]/domint;
    }

    //initialization of the densification coefficient vector
    DRT::Element*   element = discret_->lRowElement(0);
    Teuchos::RCP<MAT::Material>  mat = element->Material();

    if (mat->MaterialType() == INPAR::MAT::m_matlist)
    {
      Teuchos::RCP<const MAT::MatList> actmat = Teuchos::rcp_static_cast<const MAT::MatList>(mat);

      for (int k = 0;k<numscal_;++k)
      {
        const int matid = actmat->MatID(k);
        Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);

        if (singlemat->MaterialType() == INPAR::MAT::m_ion)
        {
          Teuchos::RCP<const MAT::Ion> actsinglemat = Teuchos::rcp_static_cast<const MAT::Ion>(singlemat);

          densific_[k] = actsinglemat->Densification();

          if (densific_[k] < 0.0) dserror("received negative densification value");
        }
        else
          dserror("material type is not allowed");
      }
    }

    if (mat->MaterialType() == INPAR::MAT::m_ion) // for a single species calculation
    {
      Teuchos::RCP<const MAT::Ion> actmat = Teuchos::rcp_static_cast<const MAT::Ion>(mat);
      densific_[0] = actmat->Densification();
      if (densific_[0] < 0.0) dserror("received negative densification value");
      if (numscal_ > 1) dserror("Single species calculation but numscal = %d > 1",numscal_);
    }
  }

  return;
} // ScaTraTimIntImpl::SetupElchNatConv


/*-------------------------------------------------------------------------*
 | valid parameters for the diffusion-conduction formulation    ehrl 12/12 |
 *-------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::ValidParameterDiffCond()
{
  if(myrank_==0)
  {
    if(DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(*elchparams_,"MOVINGBOUNDARY")!=INPAR::ELCH::elch_mov_bndry_no)
      dserror("Moving boundaries are not supported in the ELCH diffusion-conduction framework!!");

    if(DRT::INPUT::IntegralValue<int>(*elchparams_,"NATURAL_CONVECTION"))
      dserror("Natural convection is not supported in the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(*params_,"SOLVERTYPE"))!= INPAR::SCATRA::solvertype_nonlinear)
      dserror("The only solvertype supported by the ELCH diffusion-conduction framework is the non-linar solver!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(*params_,"CONVFORM"))!= INPAR::SCATRA::convform_convective)
      dserror("Only the convective formulation is supported so far!!");

    if((DRT::INPUT::IntegralValue<int>(*params_,"NEUMANNINFLOW"))== true)
      dserror("Neuman inflow BC's are not supported by the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<int>(*params_,"CONV_HEAT_TRANS"))== true)
      dserror("Convective heat transfer BC's are not supported by the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::FSSUGRDIFF>(*params_,"FSSUGRDIFF"))!= INPAR::SCATRA::fssugrdiff_no)
      dserror("Subgrid diffusivity is not supported by the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<int>(*params_,"BLOCKPRECOND"))== true)
          dserror("Block preconditioner is not supported so far!!");

    // Parameters defined in "SCALAR TRANSPORT DYNAMIC"
    Teuchos::ParameterList& scatrastabparams = params_->sublist("STABILIZATION");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(scatrastabparams,"STABTYPE"))!= INPAR::SCATRA::stabtype_no_stabilization)
      dserror("No stabilization is necessary for solving the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(scatrastabparams,"DEFINITION_TAU"))!= INPAR::SCATRA::tau_zero)
      dserror("No stabilization is necessary for solving the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalTau>(scatrastabparams,"EVALUATION_TAU"))!= INPAR::SCATRA::evaltau_integration_point)
      dserror("Evaluation of stabilization parameter only at Gauss points!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(scatrastabparams,"EVALUATION_MAT"))!= INPAR::SCATRA::evalmat_integration_point)
      dserror("Evaluation of material only at Gauss points!!");

    if((DRT::INPUT::IntegralValue<INPAR::SCATRA::Consistency>(scatrastabparams,"CONSISTENCY"))!= INPAR::SCATRA::consistency_no)
          dserror("Consistence formulation is not in the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<int>(scatrastabparams,"SUGRVEL"))== true)
          dserror("Subgrid velocity is not incoperated in the ELCH diffusion-conduction framework!!");

    if((DRT::INPUT::IntegralValue<int>(scatrastabparams,"ASSUGRDIFF"))== true)
          dserror("Subgrid diffusivity is not incoperated in the ELCH diffusion-conduction framework!!");
  }

  return;
}


/*----------------------------------------------------------------------*
 | Initialize Nernst-BC                                      ehrl 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::InitNernstBC()
{
  // access electrode kinetics condition
  std::vector<DRT::Condition*> Elchcond;
  discret_->GetCondition("ElectrodeKinetics",Elchcond);
  int numcond = Elchcond.size();

  for(int icond = 0; icond < numcond; icond++)
  {
    // check if Nernst-BC is defined on electrode kinetics condition
    if (Elchcond[icond]->GetInt("kinetic model")==INPAR::SCATRA::nernst)
    {
      switch(elchtype_)
      {
        case INPAR::ELCH::elchtype_diffcond:
        {
          // this vector must not be defined more than once!!
          if(icond==0)
            ektoggle_ = LINALG::CreateVector(*(discret_->DofRowMap()),true);

          // 1.0 for electrode-kinetics toggle
          const double one = 1.0;

          // global node id's which are part of the Nernst-BC
          const std::vector<int>* nodegids = Elchcond[icond]->Nodes();

          // loop over all global nodes part of the Nernst-BC
          for (int ii = 0; ii< (int (nodegids->size())); ++ii)
          {
            if(discret_->NodeRowMap()->MyGID((*nodegids)[ii]))
            {
              // get node with global node id (*nodegids)[ii]
              DRT::Node* node=discret_->gNode((*nodegids)[ii]);

              // get global dof ids of all dof's with global node id (*nodegids)[ii]
              std::vector<int> nodedofs=discret_->Dof(node);

              // define electrode kinetics toggle
              // later on this toggle is used to blanck the sysmat and rhs
              ektoggle_->ReplaceGlobalValues(1,&one,&nodedofs[numscal_]);
            }
          }
          break;
        }
        default:
        {
          dserror("Nernst BC is only available with div i as closure equation");
          break;
        }
      }
    }
  }

  // At element level the Nernst condition has to be handled like a DC
  if(ektoggle_!=Teuchos::null)
    dctoggle_->Update(1.0,*ektoggle_,1.0);

  return;
} //SCATRA::ScaTraTimIntImpl::InitNernstBC



/*----------------------------------------------------------------------*
 | calculate initial electric potential field at t=t_0         gjb 04/10|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::CalcInitialPotentialField()
{
  if (DRT::INPUT::IntegralValue<int>(*params_,"INITPOTCALC")==INPAR::SCATRA::initpotcalc_yes)
  {
    // time measurement:
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + calc initial potential field");
    if (myrank_ == 0)
      std::cout<<"SCATRA: calculating initial field for electric potential"<<std::endl;

    // are we really at step 0?
    dsassert(step_==0,"Step counter is not 0");

    // construct intermediate vectors
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    Teuchos::RCP<Epetra_Vector> rhs = LINALG::CreateVector(*dofrowmap,true);
    Teuchos::RCP<Epetra_Vector> phi0 = LINALG::CreateVector(*dofrowmap,true);
    phi0->Update(1.0,*phinp_,0.0);
    Teuchos::RCP<Epetra_Vector> inc = LINALG::CreateVector(*dofrowmap,true);

    // zero out matrix entries
    sysmat_->Zero();

    // evaluate Dirichlet boundary conditions at time t=0
    // the values should match your initial field at the boundary!
    ApplyDirichletBC(time_,phin_,Teuchos::null);

    // contributions due to Neumann b.c. or ElectrodeKinetics b.c.
    // have to be summed up here, and applied
    // as a current flux condition at the potential field!

    // so far: fluxes resulting from Neuman and electrochemical boundary conditions are not considered in the framework!

    // Electrode kinetics:
    // If, e.g., the initial field does not match the applied boundary conditions
    // (e.g. relaxation process of a stationary concentration field),
    // the aforementioned strategy cannot be applied to our system
    // but nevertheless the potential level has to be fixed.

    // create the parameters for the discretization
    Teuchos::ParameterList eleparams;

    // action for elements
    eleparams.set<int>("action",SCATRA::calc_elch_initial_potential);

    // set type of scalar transport problem
    eleparams.set<int>("scatratype",scatratype_);

    // factor F/RT
    eleparams.set("frt",frt_);

    // parameters for stabilization
    eleparams.sublist("STABILIZATION") = params_->sublist("STABILIZATION");

    // parameters for Elch/DiffCond formulation
    eleparams.sublist("DIFFCOND") = elchparams_->sublist("DIFFCOND");

    //provide displacement field in case of ALE
    eleparams.set("isale",isale_);
    if (isale_)
      discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

    // set vector values needed by elements
    discret_->ClearState();
    discret_->SetState("phi0",phin_);

    // call loop over elements
    discret_->Evaluate(eleparams,sysmat_,rhs);
    discret_->ClearState();

    // finalize the complete matrix
    sysmat_->Complete();

    // project residual such that only part orthogonal to nullspace is considered
    if (projector_!=Teuchos::null)
      projector_->ApplyPT(*residual_);

    // apply Dirichlet boundary conditions to system matrix
    LINALG::ApplyDirichlettoSystem(sysmat_,phi0,rhs,phi0,*(dbcmaps_->CondMap()));

    // solve the system linear incrementally:
    // the system is linear and therefore it converges in a single step, but
    // an incremental solution procedure allows the solution for the potential field
    // to converge to an "defined" potential level due to initial conditions!
    solver_->Solve(sysmat_->EpetraOperator(),inc,rhs,true,true,projector_);

    // update the original initial field
    phi0->Update(1.0,*inc,1.0);

    // copy solution of initial potential field to the solution vectors
    Teuchos::RCP<Epetra_Vector> onlypot = splitter_->ExtractCondVector(phi0);
    // insert values into the whole solution vectors
    splitter_->InsertCondVector(onlypot,phinp_);
    splitter_->InsertCondVector(onlypot,phin_);

    // reset the matrix (and its graph!) since we solved
    // a very special problem here that has a different sparsity pattern
    if (DRT::INPUT::IntegralValue<int>(*params_,"BLOCKPRECOND"))
      BlockSystemMatrix()->Reset();
    else
      SystemMatrix()->Reset();
  }

  // go on!
  return;
} // ScaTraTimIntImpl::CalcInitialPotentialField


/*----------------------------------------------------------------------*
 |  calculate conductivity of electrolyte solution             gjb 07/09|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_SerialDenseVector> SCATRA::ScaTraTimIntElch::ComputeConductivity()
{
  // we perform the calculation on element level hiding the material access!
  // the initial concentration distribution has to be uniform to do so!!

  // create the parameters for the elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action",SCATRA::calc_elch_conductivity);
  eleparams.set<int>("scatratype",scatratype_);

  //provide displacement field in case of ALE
  eleparams.set("isale",isale_);
  if (isale_)
    discret_->AddMultiVectorToParameterList(eleparams,"dispnp",dispnp_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // pointer to current element
  DRT::Element* actele = discret_->lRowElement(0);

  // get element location vector, dirichlet flags and ownerships
  std::vector<int> lm;  // location vector
  std::vector<int> lmowner;  // processor which owns DOFs
  std::vector<int> lmstride;  // nodal block sizes in element matrices

  actele->LocationVector(*discret_,lm,lmowner,lmstride);

  // define element matrices and vectors
  // -- which are empty and unused, just to satisfy element Evaluate()
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  // call the element evaluate method of the first row element
  int err = actele->Evaluate(eleparams,*discret_,lm,elematrix1,elematrix2,*sigma_,elevector2,elevector3);
  if (err) dserror("error while computing conductivity");
  discret_->ClearState();

  return sigma_;
} // ScaTraTimIntImpl::ComputeConductivity


/*----------------------------------------------------------------------*
 | apply galvanostatic control                                gjb 11/09 |
 *----------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntElch::ApplyGalvanostaticControl()
{
  // for galvanostatic ELCH applications we have to adjust the
  // applied cell voltage and continue Newton-Raphson iterations until
  // we reach the desired value for the electric current.

  if (DRT::INPUT::IntegralValue<int>(*elchparams_,"GALVANOSTATIC"))
  {
    // set time derivate parameters of applied voltage for a double layer capacitance current density,
    if(dlcapexists_)
      ComputeTimeDerivPot0(false);

    std::vector<DRT::Condition*> cond;
    discret_->GetCondition("ElectrodeKinetics",cond);
    if (!cond.empty())
    {
      const unsigned condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");
      const unsigned condid_anode = elchparams_->get<int>("GSTATCONDID_ANODE");
      int gstatitemax = (elchparams_->get<int>("GSTATITEMAX"));
      double gstatcurrenttol = (elchparams_->get<double>("GSTATCURTOL"));
      const int curvenum = elchparams_->get<int>("GSTATCURVENO");
      const double tol = elchparams_->get<double>("GSTATCONVTOL");
      const double effective_length = elchparams_->get<double>("GSTAT_LENGTH_CURRENTPATH");

      // There are maximal two electrode conditions by definition
      // current flow i at electrodes
      Teuchos::RCP<std::vector<double> > actualcurrent = Teuchos::rcp(new std::vector<double>(2,0.0));
      // residual at electrodes = i*timefac
      Teuchos::RCP<std::vector<double> > currresidual= Teuchos::rcp(new std::vector<double>(2,0.0));
      Teuchos::RCP<std::vector<double> > currtangent = Teuchos::rcp(new std::vector<double>(2,0.0));
      Teuchos::RCP<std::vector<double> > electrodesurface = Teuchos::rcp(new std::vector<double>(2,0.0));
      Teuchos::RCP<std::vector<double> > electrodepot = Teuchos::rcp(new std::vector<double>(2,0.0));
      Teuchos::RCP<std::vector<double> > meanoverpot = Teuchos::rcp(new std::vector<double>(2,0.0));
      double meanelectrodesurface(0.0);
      //Assumption: Residual at BV1 is the negative of the value at BV2, therefore only the first residual is calculated
      double newtonrhs(0.0);

      // for all time integration schemes, compute the current value for phidtnp
      // this is needed for evaluating charging currents due to double-layer capacity
      // This may only be called here and not inside OutputSingleElectrodeInfo!!!!
      // Otherwise you modify your output to file called during Output()
      ComputeTimeDerivative();

      double targetcurrent = DRT::Problem::Instance()->Curve(curvenum-1).f(time_);
      double timefac = 1.0/ResidualScaling();

      double currtangent_anode(0.0);
      double currtangent_cathode(0.0);
      double potinc_ohm(0.0);
      double potdiffbulk(0.0);
      double potdiffcell(0.0);

      if(cond.size()>2)
        dserror("The framework may not work for geometrical setups containing more than two electrodes! \n"
                "If you need it, check the framework exactly!!");

      // loop over all BV
      // degenerated to a loop over 2 (user-specified) BV conditions
      for (unsigned int icond = 0; icond < cond.size(); icond++)
      {
        // note: only the potential at the boundary with id condid_cathode will be adjusted!
        OutputSingleElectrodeInfo(
            cond[icond],
            icond,
            false,
            false,
            (*actualcurrent)[icond],
            (*currtangent)[icond],
            (*currresidual)[icond],
            (*electrodesurface)[icond],
            (*electrodepot)[icond],
            (*meanoverpot)[icond]
            );

        if(cond.size()==2)
        {
          // In the case the actual current is zero, we assume that the first electrode is the cathode
          if((*actualcurrent)[icond]<0.0 and condid_cathode != icond)
            dserror("The defined GSTATCONDID_CATHODE does not match the actual current flow situation!!");
          else if((*actualcurrent)[icond]>0.0 and condid_anode != icond)
            dserror("The defined GSTATCONDID_ANODE does not match the actual current flow situation!!");
        }
      } // end loop over electrode kinetics

      if(cond.size()==1)
      {
        if(condid_cathode != 0 or condid_anode!=1)
          dserror("The defined GSTATCONDID_CATHODE and GSTATCONDID_ANODE is wrong for a setup with only one electrode!!\n"
                  "Choose: GSTATCONDID_CATHODE=0 and GSTATCONDID_ANODE=1");
      }

      // get the applied electrode potential of the cathode
      const double potold = cond[condid_cathode]->GetDouble("pot");
      double potnew = potold;

      // bulk voltage loss = V_A - eta_A - (V_C -eta_C)
      potdiffbulk = ((*electrodepot)[condid_anode]-(*meanoverpot)[condid_anode])-((*electrodepot)[condid_cathode]-(*meanoverpot)[condid_cathode]);
      // cell voltage loss = V_A - V_C
      potdiffcell=(*electrodepot)[condid_anode]-(*electrodepot)[condid_cathode];
      // tanget at anode and cathode
      currtangent_anode=(*currtangent)[condid_anode];
      currtangent_cathode=(*currtangent)[condid_cathode];

      // The linarization of potential increment is based on the cathode side!!
      {
        //Assumption: Residual at BV1 is the negative of the value at BV2, therefore only the first residual is calculated
        // newtonrhs = -residual, with the definition:  residual := timefac*(-I + I_target)
        newtonrhs = + (*currresidual)[condid_cathode] - (timefac*targetcurrent); // newtonrhs is stored only from cathode!
        if (myrank_==0)
        {
          // format output
          std::cout.precision(3);
          std::cout<<"\n  GALVANOSTATIC MODE:\n";
          std::cout<<"  +--------------------------------------------------------------------------" <<std::endl;
          std::cout<<"  | Convergence check: " <<std::endl;
          std::cout<<"  +--------------------------------------------------------------------------" <<std::endl;
          std::cout<<"  | iteration:                          "<<setw(7)<<std::right<<gstatnumite_<<" / "<<gstatitemax<<std::endl;
          std::cout<<"  | actual reaction current at cathode: "<<std::scientific<<setw(12)<<std::right<<(*actualcurrent)[condid_cathode]<<std::endl;
          std::cout<<"  | required total current at cathode:  "<<setw(12)<<std::right<<targetcurrent<<std::endl;
          std::cout<<"  | negative residual (rhs):            "<<setw(12)<<std::right<<newtonrhs<<std::endl;
          std::cout<<"  +--------------------------------------------------------------------------" <<std::endl;
        }

        if (gstatnumite_ > gstatitemax)
        {
          if (myrank_==0)
          {
            std::cout<<"  | --> maximum number iterations reached. Not yet converged!"<<std::endl;
            std::cout<<"  +--------------------------------------------------------------------------" <<std::endl <<std::endl;
          }
          return true; // we proceed to next time step
        }
        else if (abs(newtonrhs)< gstatcurrenttol)
        {
          if (myrank_==0)
          {
            std::cout<<"  | --> Newton-RHS-Residual is smaller than " << gstatcurrenttol<< "!" << std::endl;
            std::cout<<"  +--------------------------------------------------------------------------" <<std::endl <<std::endl;
          }
          return true; // we proceed to next time step
        }
        // electric potential increment of the last iteration
        else if ((gstatnumite_ > 1) and (abs(gstatincrement_)< (1+abs(potold))*tol)) // < ATOL + |pot|* RTOL
        {
          if (myrank_==0)
          {
            std::cout<<"  | --> converged: |"<<gstatincrement_<<"| < "<<(1+abs(potold))*tol<<std::endl;
            std::cout<<"  +--------------------------------------------------------------------------" <<std::endl <<std::endl;
          }
          return true; // galvanostatic control has converged
        }

        // update applied electric potential
        // potential drop ButlerVolmer conditions (surface ovepotential) and in the electrolyte (ohmic overpotential) are conected in parallel:

        // 2 different versions:
        // I_0 = I_BV1 = I_ohmic = I_BV2
        // R(I_target, I) = R_BV1(I_target, I) = R_ohmic(I_target, I) = -R_BV2(I_target, I)
        // delta E_0 = delta U_BV1 + delta U_ohmic - (delta U_BV2)
        // => delta E_0 = (R_BV1(I_target, I)/J) + (R_ohmic(I_target, I)/J) - (-R_BV2(I_target, I)/J)
        potinc_ohm=(-1.0*effective_length*newtonrhs)/((*sigma_)[numscal_]*timefac*(*electrodesurface)[condid_cathode]);

        // electrode surface of the cathode
        meanelectrodesurface=(*electrodesurface)[condid_cathode];

        // safety check
        if (abs((*currtangent)[condid_cathode])<EPS13)
          dserror("Tangent in galvanostatic control is near zero: %lf",(*currtangent)[condid_cathode]);
      }

      if(cond.size()==2)
      {
        // mean electrode surface of the cathode abd anode
        meanelectrodesurface=((*electrodesurface)[0]+(*electrodesurface)[1])/2;

        // actual potential difference is used to calculate the current path length
        // -> it is possible to compute the new ohmic potential step
        //    without the input parameter GSTAT_LENGTH_CURRENTPATH
        if (effective_length==-1.0)
        {
          potinc_ohm = (potdiffbulk*newtonrhs)/(timefac*(*actualcurrent)[condid_cathode]);
        }

        // Do not update the cell potential for small currents
        if(abs((*actualcurrent)[condid_cathode]) < EPS10)
          potinc_ohm = 0.0;

        // the current flow at both electrodes has to be the same within the solution tolerances
        if(abs((*actualcurrent)[condid_cathode]+(*actualcurrent)[condid_anode])>EPS8)
        {
          if (myrank_==0)
          {
            std::cout.precision(3);
            std::cout << "Warning!!!" << std::endl;
            std::cout << "The difference of the current flow at anode and cathode is " << abs((*actualcurrent)[condid_cathode]+(*actualcurrent)[condid_anode])
                      << " larger than " << EPS8 << std::endl;
          }
        }
      }

      // Newton step:  Jacobian * \Delta pot = - Residual
      const double potinc_cathode = newtonrhs/((-1)*currtangent_cathode);
      double potinc_anode = 0.0;
      if (abs(currtangent_anode)>EPS13) // anode surface overpotential is optional
      {
        potinc_anode = newtonrhs/((-1)*currtangent_anode);
      }
      gstatincrement_ = (potinc_cathode+potinc_anode+potinc_ohm);
      // update electric potential
      potnew += gstatincrement_;

      if(myrank_==0)
      {
        std::cout<<"  | ohmic potential increment is calculated based on" <<std::endl;
        if(effective_length != -1.0)
          std::cout<<"  | the GSTAT_LENGTH_CURRENTPATH as defined in the input file!" <<std::endl;
        else
          std::cout<<"  | the ohmic resistance calculated from applied potential and current flow!" <<std::endl;
        std::cout<<"  +--------------------------------------------------------------------------" <<std::endl;
        std::cout<<"  | Defined GSTAT_LENGTH_CURRENTPATH:               "<<setw(12)<<std::right<< effective_length << std::endl;
        if((*actualcurrent)[condid_cathode]!=0.0)
          std::cout<<"  | Approximated GSTAT_LENGTH_CURRENTPATH:          "<<setw(12)<<std::right<< potdiffbulk/(*actualcurrent)[condid_cathode]*(*sigma_)[numscal_]*meanelectrodesurface <<std::endl;
        std::cout<<"  | (not directly used to calculate the ohmic potential increment)"<<std::endl;
        std::cout<<"  | New guess for:                                  "<<std::endl;
        std::cout<<"  | - ohmic potential increment:                    "<<setw(12)<<std::right<< potinc_ohm <<std::endl;
        std::cout<<"  | - overpotential increment cathode (condid " << condid_cathode <<"):   " <<setw(12)<<std::right<< potinc_cathode << std::endl;
        std::cout<<"  | - overpotential increment anode (condid " << condid_anode <<"):     " <<setw(12)<<std::right<< potinc_anode << std::endl;
        std::cout<<"  | -> total increment for potential:               " <<setw(12)<<std::right<< gstatincrement_ << std::endl;
        std::cout<<"  +--------------------------------------------------------------------------" <<std::endl;
        std::cout<<"  | old potential at the cathode (condid "<<condid_cathode <<"):     "<<setw(12)<<std::right<<potold<<std::endl;
        std::cout<<"  | new potential at the cathode (condid "<<condid_cathode <<"):     "<<setw(12)<<std::right<<potnew<<std::endl;
        std::cout<<"  | new applied potential difference (condid "<<condid_cathode <<"): "<<setw(12)<<std::right<<potdiffcell<<std::endl;
        std::cout<<"  +--------------------------------------------------------------------------" <<std::endl<<std::endl;
      }

//      // print additional information
//      if (myrank_==0)
//      {
//        std::cout<< "  actualcurrent - targetcurrent = " << ((*actualcurrent)[condid_cathode]-targetcurrent) << std::endl;
//        std::cout<< "  conductivity                  = " << (*sigma_)[numscal_] << std::endl<< std::endl;
//        std::cout<< "  currtangent_cathode:  " << currtangent_cathode << std::endl;
//        std::cout<< "  currtangent_anode:    " << currtangent_anode << std::endl;
//        std::cout<< "  actualcurrent cathode:    " << (*actualcurrent)[condid_cathode] << std::endl;
//        std::cout<< "  actualcurrent anode:  " << (*actualcurrent)[condid_anode] << std::endl;
//      }

      // replace potential value of the boundary condition (on all processors)
      cond[condid_cathode]->Add("pot",potnew);
      gstatnumite_++;
      return false; // not yet converged -> continue Newton iteration with updated potential
    }
  }
  return true; //default

} // end ApplyGalvanostaticControl()

/*----------------------------------------------------------------------*
 | evaluate contribution of electrode kinetics to eq. system  gjb 02/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::EvaluateSolutionDependingBC(
    Teuchos::RCP<LINALG::SparseOperator> matrix,
    Teuchos::RCP<Epetra_Vector>          rhs
)
{
  // time measurement: evaluate condition 'ElectrodeKinetics'
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ElectrodeKinetics'");

  discret_->ClearState();

  // create an parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_elch_electrode_kinetics);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("frt",frt_); // factor F/RT
  condparams.set("isale",isale_);

  // parameters for Elch/DiffCond formulation
  condparams.sublist("DIFFCOND") = elchparams_->sublist("DIFFCOND");

  if (isale_)   //provide displacement field in case of ALE
    discret_->AddMultiVectorToParameterList(condparams,"dispnp",dispnp_);

  // add element parameters and set state vectors according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  std::string condstring("ElectrodeKinetics");
  // evaluate ElectrodeKinetics conditions at time t_{n+1} or t_{n+alpha_F}
  discret_->EvaluateCondition(condparams,matrix,Teuchos::null,rhs,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  //Add linearization of NernstCondition to system matrix
  if(ektoggle_!=Teuchos::null)
    LinearizationNernstCondition();

  return;
} // ScaTraTimIntElch::EvaluateSolutionDependingBC


/*----------------------------------------------------------------------*
 | Add Linearization for Nernst-BC                           ehrl 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::LinearizationNernstCondition()
{
  // Blank rows with Nernst -BC (inclusive diagonal entry)
  // Nernst-BC is a additional constraint coupled to the original system of equation
  sysmat_->ApplyDirichlet(ektoggle_,false);
  LINALG::ApplyDirichlettoSystem(increment_,residual_,zeros_,ektoggle_);

  discret_->ClearState();

  // create an parameter list
  Teuchos::ParameterList condparams;
  //update total time for time curve actions
  AddTimeIntegrationSpecificVectors();
  // action for elements
  condparams.set<int>("action",SCATRA::bd_calc_elch_linearize_nernst);
  condparams.set<int>("scatratype",scatratype_);
  condparams.set("isale",isale_);

  // parameters for Elch/DiffCond formulation
  condparams.sublist("DIFFCOND") = elchparams_->sublist("DIFFCOND");

  // add element parameters and set state vectors according to time-integration scheme
  // we need here concentration at t+np
  discret_->SetState("phinp",phinp_);

  std::string condstring("ElectrodeKinetics");
  // evaluate ElectrodeKinetics conditions at time t_{n+1} or t_{n+alpha_F}
  // phinp (view to phinp)
  discret_->EvaluateCondition(condparams,sysmat_,Teuchos::null,residual_,Teuchos::null,Teuchos::null,condstring);
  discret_->ClearState();

  return;
} //  SCATRA::ScaTraTimIntImpl::LinearizationNernstCondition()


/*----------------------------------------------------------------------*
 | check for zero/negative concentration values               gjb 01/10 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::CheckConcentrationValues(Teuchos::RCP<Epetra_Vector> vec)
{
  // action only for ELCH applications

  // for NURBS discretizations we skip the following check.
  // Control points (i.e., the "nodes" and its associated dofs can be located
  // outside the domain of interest. Thus, they can have negative
  // concentration values although the concentration solution is positive
  // in the whole computational domain!
  if(dynamic_cast<DRT::NURBS::NurbsDiscretization*>(discret_.get())!=NULL)
    return;

  // this option can be helpful in some rare situations
  bool makepositive(false);

  std::vector<int> numfound(numscal_,0);
#if 0
  std::stringstream myerrormessage;
#endif
  for (int i = 0; i < discret_->NumMyRowNodes(); i++)
  {
    DRT::Node* lnode = discret_->lRowNode(i);
    std::vector<int> dofs;
    dofs = discret_->Dof(lnode);

    for (int k = 0; k < numscal_; k++)
    {
      const int lid = discret_->DofRowMap()->LID(dofs[k]);
      if (((*vec)[lid]) < EPS13 )
      {
        numfound[k]++;
        if (makepositive)
          ((*vec)[lid]) = EPS13;
#if 0
        myerrormessage<<"PROC "<<myrank_<<" dof index: "<<k<<setprecision(7)<<scientific<<
            " val: "<<((*vec)[lid])<<" node gid: "<<lnode->Id()<<
            " coord: [x] "<< lnode->X()[0]<<" [y] "<< lnode->X()[1]<<" [z] "<< lnode->X()[2]<<std::endl;
#endif
      }
    }
  }

  // print warning to screen
  for (int k = 0; k < numscal_; k++)
  {
    if (numfound[k] > 0)
    {
      std::cout<<"WARNING: PROC "<<myrank_<<" has "<<numfound[k]<<
      " nodes with zero/neg. concentration values for species "<<k;
      if (makepositive)
        std::cout<<"-> were made positive (set to 1.0e-13)"<<std::endl;
      else
        std::cout<<std::endl;
    }
  }

#if 0
  // print detailed info to error file
  for(int p=0; p < discret_->Comm().NumProc(); p++)
  {
    if (p==myrank_) // is it my turn?
    {
      // finish error message
      myerrormessage.flush();

      // write info to error file
      if ((errfile_!=NULL) and (myerrormessage.str()!=""))
      {
        fprintf(errfile_,myerrormessage.str().c_str());
        // std::cout<<myerrormessage.str()<<std::endl;
      }
    }
    // give time to finish writing to file before going to next proc ?
    discret_->Comm().Barrier();
  }
#endif

  // so much code for a simple check!
  return;
} // ScaTraTimIntImpl::CheckConcentrationValues
