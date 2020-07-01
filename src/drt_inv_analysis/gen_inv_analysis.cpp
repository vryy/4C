/*----------------------------------------------------------------------*/
/*! \file
\brief gen inv analysis

\level 2

 */
/*----------------------------------------------------------------------*/

#include "gen_inv_analysis.H"
#include "../drt_cut/cut_clnwrapper.H"
#include "../drt_inpar/inpar_material.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_hdf.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/growth_law.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/matpar_parameter.H"
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/neohooke.H"
#include "../drt_matelast/elast_coup1pow.H"
#include "../drt_matelast/elast_coup2pow.H"
#include "../drt_matelast/elast_coup3pow.H"
#include "../drt_matelast/elast_coup13apow.H"
#include "../drt_mat/viscoelasthyper.H"
#include "../drt_matelast/elast_coupblatzko.H"
#include "../drt_matelast/elast_couplogneohooke.H"
#include "../drt_matelast/elast_couplogmixneohooke.H"
#include "../drt_matelast/elast_coupexppol.H"
#include "../drt_matelast/elast_coupmooneyrivlin.H"
#include "../drt_matelast/elast_coupneohooke.H"
#include "../drt_matelast/elast_coupsimopister.H"
#include "../drt_matelast/elast_coupSaintVenantKirchhoff.H"
#include "../drt_matelast/elast_iso1pow.H"
#include "../drt_matelast/elast_iso2pow.H"
#include "../drt_matelast/elast_isoexpopow.H"
#include "../drt_matelast/elast_isomooneyrivlin.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_isoyeoh.H"
#include "../drt_matelast/elast_vologden.H"
#include "../drt_matelast/elast_volpenalty.H"
#include "../drt_matelast/elast_volsussmanbathe.H"
#include "../drt_matelast/elast_volpow.H"
#include "../drt_matelast/visco_coupmyocard.H"
#include "../drt_matelast/visco_isoratedep.H"
#include "../drt_matelast/visco_genmax.H"
#include "../drt_matelast/visco_fract.H"
#include "../drt_mat/stvenantkirchhoff.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_structure/strtimint_create.H"
#include "../drt_structure/strtimint.H"
#include "../drt_stru_multi/microstatic.H"
#include "../linalg/linalg_multiply.H"
#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_solver.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_CrsMatrix.h"
#include <cstdlib>
#include <ctime>
#include <iostream>
#include "../drt_lib/drt_discret.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_inpar/inpar_invanalysis.H"
#include "../drt_inpar/inpar_fsi.H"

#include "../drt_adapter/ad_str_factory.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_fsi/fsi_monolithic.H"

#include "../drt_ale/ale_utils_clonestrategy.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_lungmonolithic.H"
#include "../drt_fsi/fsi_lungmonolithic_fluidsplit.H"
#include "../drt_fsi/fsi_lungmonolithic_structuresplit.H"
#include "../drt_fsi/fsi_constrmonolithic_fluidsplit.H"
#include "../drt_fsi/fsi_constrmonolithic_structuresplit.H"
#include "../drt_fsi/fsi_mortarmonolithic_fluidsplit.H"
#include "../drt_fsi/fsi_mortarmonolithic_structuresplit.H"
#include "../drt_fsi/fsi_fluidfluidmonolithic_structuresplit_nonox.H"
#include "../drt_fsi/fsi_fluidfluidmonolithic_fluidsplit_nonox.H"
#include "../drt_fsi/fsi_fluidfluidmonolithic_structuresplit.H"
#include "../drt_fsi/fsi_fluidfluidmonolithic_fluidsplit.H"
#include "../drt_fsi/fsi_slidingmonolithic_fluidsplit.H"
#include "../drt_fsi/fsi_slidingmonolithic_structuresplit.H"
#include "../drt_fsi/fsi_dirichletneumann.H"
#include "../drt_fsi/fsi_dirichletneumannslideale.H"
#include "../drt_fsi/fsi_dirichletneumann_volcoupl.H"

#include "../drt_binstrategy/binning_strategy.H"

#include "../drt_structure/stru_resulttest.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_cut/cut_kernel.H"

/*----------------------------------------------------------------------*/
/* standard constructor */
STR::GenInvAnalysis::GenInvAnalysis(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<IO::DiscretizationWriter> output)
    : discret_(dis), solver_(solver), output_(output)
{
  int myrank = dis->Comm().MyPID();

  nodescentdirection_ = 0;

  // input parameters inverse analysis
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  //  tolerance for the curve fitting
  tol_ = iap.get<double>("INV_ANA_TOL");

  ndim_ = DRT::Problem::Instance()->NDim();

  forward_prb_type_ = DRT::INPUT::IntegralValue<ProblemType>(iap, "FORWARD_PROBLEMTYP");
  materialhashistory_ = false;

  switch (DRT::INPUT::IntegralValue<INPAR::STR::MeasurementType>(iap, "MEASUREMENT_TYPE"))
  {
    case INPAR::STR::meas_dofs:
    {
      meas_type_ = dof_based;
    }
    break;
    case INPAR::STR::meas_points:
    {
      meas_type_ = point_based;
    }
    break;
    default:
      dserror("Unknown Type of Monitor File MEASUREMENT_TYPE. Provide proper Type and file.");
      break;
  }

  if (meas_type_ == dof_based)
    // open monitor file and read it
    ReadMonitorDofBased(myrank);
  else if (meas_type_ == point_based)
    // open point based monitor file and read it
    ReadMonitorPointBased(myrank);

  // error: difference of the measured to the calculated curve
  error_ = 1.0E6;
  error_o_ = 1.0E6;
  //
  error_grad_ = 1.0E6;
  error_grad_o_ = 1.0E6;

  // training parameter
  mu_ = iap.get<double>("INV_INITREG");
  mu_o_ = mu_;
  kappa_multi_ = 1.0;

  // update strategy for mu
  switch (DRT::INPUT::IntegralValue<INPAR::STR::RegStratUpdate>(iap, "UPDATE_REG"))
  {
    case INPAR::STR::reg_update_grad:
    {
      reg_update_ = grad_based;
    }
    break;
    case INPAR::STR::reg_update_res:
    {
      reg_update_ = res_based;
    }
    break;
    default:
      dserror("Unknown update strategy for regularization parameter! Fix your input file");
      break;
  }

  // special inverse analysis types
  spec_inv_ana_mult_ = false;
  spec_inv_ana_coup_ = false;
  switch (DRT::INPUT::IntegralValue<INPAR::STR::SpecInvAnalysisType>(iap, "SPECIAL_INV_ANALYSIS"))
  {
    case INPAR::STR::spec_inv_mult:
    {
      spec_inv_ana_mult_ = true;
    }
    break;
    case INPAR::STR::spec_inv_coup:
    {
      spec_inv_ana_coup_ = true;
    }
    break;
    default:  // normal inverse analysis --> do nothing
      break;
  }

  check_neg_params_ = DRT::INPUT::IntegralValue<int>(iap, "PARAM_BOUNDS");

  // list of materials for each problem instance that should be fitted
  for (unsigned prob = 0; prob < DRT::Problem::NumInstances(); ++prob)
  {
    const Teuchos::ParameterList& myiap = DRT::Problem::Instance(prob)->InverseAnalysisParams();
    std::set<int> myset;

    std::string word1;
    std::istringstream matliststream(Teuchos::getNumericStringParameter(myiap, "INV_LIST"));
    while (matliststream >> word1)
    {
      int word1_int = std::atoi(word1.c_str());
      if (word1_int != -1)  // this means there was no matlist specified in the input file
        myset.insert(word1_int);
    }
    matset_.push_back(myset);
  }

  // do patch stuff
  patches_ = DRT::INPUT::IntegralValue<bool>(iap, "PATCHES");
  // list of materials for each problem instance that should be fitted with patches
  if (patches_)
  {
    // some parameters
    numpatches_ = iap.get<int>("NUMPATCHES");
    smoothingsteps_ = iap.get<int>("SMOOTHINGSTEPSPATCHES");

    // read material list
    if (DRT::Problem::NumInstances() > 1)
      dserror("More than one problem instance %d not feasible with patches",
          DRT::Problem::NumInstances());
    const Teuchos::ParameterList& myiap = DRT::Problem::Instance(0)->InverseAnalysisParams();

    int word1;
    std::istringstream matliststream(Teuchos::getNumericStringParameter(myiap, "INV_LIST_PATCHES"));
    while (matliststream >> word1)
    {
      if (word1 != -1)  // this means there was no matlist specified in the input file
        matset_patches_.insert(word1);
    }

    // define how the patches are set
    std::string patchtype = iap.get<std::string>("DEFINEPATCHES");
    if (patchtype == "MaterialNumber")
    {
      which_patches_ = material;
      int size = matset_patches_.size();
      if (size != numpatches_)
        dserror(
            "Fatal failure in the case of patches defined by the material number:\n"
            "number of patches %d has to be identical to the number of materials %d !",
            numpatches_, size);
    }
    else if (patchtype == "Uniform")
      which_patches_ = uniform;
    else
      dserror("Type of patches not known");

    // init smoothing operator
    InitPatches();
  }


  // read material parameters from input file
  ReadInParameters();
  p_o_ = p_;
  p_print_ = p_;

  // Number of material parameters
  np_ = p_.Length();

  // controlling parameter
  numb_run_ = 0;  // counter of how many runs were made in the inverse analysis
}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::GenInvAnalysis::Integrate()
{
  const int myrank = discret_->Comm().MyPID();
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  const int max_iter = iap.get<int>("INV_ANA_MAX_RUN");
  int newfiles = DRT::INPUT::IntegralValue<int>(iap, "NEW_FILES");

  // multiple inverse analysis is just implemented for parallel use
  if (spec_inv_ana_mult_)
    dserror("Multiple inverse analysis is just implemented for parallel use.");

  // fitting loop
  do
  {
    if (!myrank)
    {
      std::cout << "#################################################################\n"
                << "########################### making Jacobian matrix ##############\n"
                << "#################################################################\n"
                << "No. of measured points: nmp_= " << nmp_
                << " # no. of parameters to fit: np_= " << np_ << std::endl;
    }

    // perturbation of material parameter (should be relative to the value that is perturbed)
    std::vector<double> perturb(np_, 0.0);
    const double alpha = iap.get<double>("INV_ALPHA");
    const double beta = iap.get<double>("INV_BETA");

    for (int i = 0; i < np_; ++i)
    {
      perturb[i] = alpha + beta * p_[i];
      if (!myrank) printf("perturbation[%d] %15.10e\n", i, perturb[i]);
    }

    // as the actual inverse analysis is on proc 0 only, do only provide storage on proc 0
    Epetra_SerialDenseMatrix cmatrix;
    if (!myrank) cmatrix.Shape(nmp_, np_ + 1);

    // loop over parameters to fit and build cmatrix
    for (int i = 0; i < np_ + 1; i++)
    {
      bool outputtofile = false;
      // output only for last run
      if (i == np_) outputtofile = true;

      if (outputtofile)
      {
        // no parameters for newfile yet
        if (newfiles)
          output_->NewResultFile((numb_run_));
        else
          output_->OverwriteResultFile();

        output_->WriteMesh(0, 0.0);
      }

      // Multi-scale: if an inverse analysis is performed on the micro-level,
      // the time and step need to be reset now. Furthermore, the result file
      // needs to be opened.
      MultiInvAnaInit();

      if (!myrank)
        std::cout << "--------------------------- run " << i + 1 << " of: " << np_ + 1
                  << " -------------------------" << std::endl;
      // make current set of material parameters

      Epetra_SerialDenseVector p_cur = p_;
      // perturb parameter i
      if (i != np_) p_cur[i] = p_[i] + perturb[i];


      // put perturbed material parameters to material laws
      discret_->Comm().Broadcast(&p_cur[0], p_cur.Length(), 0);

      SetParameters(p_cur);

      // Solve forward problem and obtain computed displacements output at the last step
      Epetra_SerialDenseVector cvector;

      switch (forward_prb_type_)
      {
        case prb_structure:
        {
          cvector = CalcCvector(outputtofile);
          break;
        }
        case prb_fsi:
        {
          cvector = CalcCvectorFSI(outputtofile, i);
          break;
        }
        default:
        {
          dserror("Chosen FORWARD_PROBLEMTYP not considered for inverse analysis");
          break;
        }
      }
      // copy displacements to sensitivity matrix
      if (!myrank)
        for (int j = 0; j < nmp_; j++) cmatrix(j, i) = cvector[j];

      if (materialhashistory_)
      {
        Teuchos::ParameterList p;
        p.set("action", "calc_struct_reset_all");
        discret_->Evaluate(
            p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
      }
    }

    discret_->Comm().Barrier();
    if (!myrank) CalcNewParameters(cmatrix, perturb);

    // set new material parameters
    discret_->Comm().Broadcast(&p_[0], p_.Length(), 0);
    SetParameters(p_);
    numb_run_++;
    discret_->Comm().Broadcast(&error_o_, 1, 0);
    discret_->Comm().Broadcast(&error_, 1, 0);
    discret_->Comm().Broadcast(&error_grad_o_, 1, 0);
    discret_->Comm().Broadcast(&error_grad_, 1, 0);
    discret_->Comm().Broadcast(&numb_run_, 1, 0);
    discret_->Comm().Broadcast(&nodescentdirection_, 1, 0);

    if (reg_update_ == grad_based)
      error_i_ = error_grad_;
    else if (reg_update_ == res_based)
      error_i_ = error_;


  } while (error_i_ > tol_ && numb_run_ < max_iter && !nodescentdirection_);


  // print results to file
  if (!myrank) PrintFile();

  // stop supporting processors in multi scale simulations
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();
  if (subcomm != Teuchos::null)
  {
    STRUMULTI::stop_np_multiscale();
  }

  return;
}

/*----------------------------------------------------------------------*/
/* analyse */
void STR::GenInvAnalysis::NPIntegrate()
{
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = DRT::Problem::Instance()->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  const int lmyrank = lcomm->MyPID();
  const int gmyrank = gcomm->MyPID();
  const int groupid = group->GroupId();
  const int ngroup = group->NumGroups();

  // make some plausibility checks:
  // number of groups must be divisible by the number of material parameters np_+1
  if (ngroup % (np_ + 1) != 0)
    dserror("# groups (%d) must divide by # material parameters + 1 (%d) ", ngroup, np_ + 1);


  // group 0 does the unpermuted solution
  // groups 1 to np_ do the permuted versions
  // proc 0 of group 0 (which is also gproc 0) does the optimization step

  if (!discret_->Filled() || !discret_->HaveDofs()) discret_->FillComplete(true, true, true);

  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  int max_itter = iap.get<int>("INV_ANA_MAX_RUN");
  int newfiles = DRT::INPUT::IntegralValue<int>(iap, "NEW_FILES");

  // do special stuff of multiple inverse analysis
  if (spec_inv_ana_mult_) NPIntegrateMult();

  // fitting loop
  do
  {
    gcomm->Barrier();

    if (!gmyrank)
    {
      std::cout << "#################################################################" << std::endl;
      std::cout << "######################## NP making Jacobian matrix ##############" << std::endl;
      std::cout << "#################################################################" << std::endl;
      printf("Measured parameters nmp_ %d # parameters to fit np_ %d NP groups %d\n", nmp_, np_,
          ngroup);
      fflush(stdout);
    }
    gcomm->Barrier();

    // perturbation of material parameter (should be relative to the value that is perturbed)
    std::vector<double> perturb(np_, 0.0);
    double alpha = iap.get<double>("INV_ALPHA");
    double beta = iap.get<double>("INV_BETA");
    for (int i = 0; i < np_; ++i)
    {
      perturb[i] = alpha + beta * p_[i];
      if (!gmyrank) printf("perturbation[%d] %15.10e\n", i, perturb[i]);
    }

    // set variables if multiple inverse analysis or normal inverse analysis
    int i;  // every group start with different material parameter perturbed
    double nmp;
    if (spec_inv_ana_mult_)
    {
      // in case of multiple inverse analysis: use measured dofs of all experiments
      nmp = nmp_mult_;

      // get experiment ID of each optimization
      int expid = iap.get<int>("EXP_ID_MULT_INV_ANA");
      i = groupid - expid * (np_ + 1);

      // to ensure i does not get smaller than zero
      // potentially, one can run a single experiment and the exp_id
      // can then be a number other than zero
      if (nexp_ == 1) i = groupid;
    }
    // normal inverse analysis
    else
    {
      nmp = nmp_;
      i = groupid;  // every group start with different material parameter perturbed
    }

    gcomm->Barrier();

    // as the actual inv analysis is on proc 0 only, do only provide storage on proc 0
    Epetra_SerialDenseMatrix cmatrix;
    if (!gmyrank) cmatrix.Shape(nmp, np_ + 1);

    //---------------------------- loop over parameters to fit and build cmatrix
    // i: every group start with different material parameter perturbed
    for (; i < np_ + 1; i += ngroup)  // loop with increments of ngroup
    {
      bool outputtofile = false;
      // output only for last run
      if (i == np_) outputtofile = true;

      if (outputtofile)
      {
        // no parameters for newfile yet
        if (newfiles)
          output_->NewResultFile((numb_run_));
        else
          output_->OverwriteResultFile();

        output_->WriteMesh(0, 0.0);
      }

      // Multi-scale: if an inverse analysis is performed on the micro-level,
      // the time and step need to be reset now. Furthermore, the result file
      // needs to be opened.
      MultiInvAnaInit();  // this might not work in here, I'm not sure

      if (!lmyrank)
        std::cout << "--------------------------- run " << i + 1 << " of: " << np_ + 1
                  << " -------------------------" << std::endl;
      // make current set of material parameters
      Epetra_SerialDenseVector p_cur = p_;
      // perturb parameter i
      if (i != np_) p_cur[i] = p_[i] + perturb[i];

      // put perturbed material parameters to material laws
      // these are different for every group
      lcomm->Broadcast(&p_cur[0], p_cur.Length(), 0);
      SetParameters(p_cur);

      // compute nonlinear problem and obtain computed displacements
      // output at the last step
      Epetra_SerialDenseVector cvector;
      cvector = CalcCvector(outputtofile, group);

      // output cvector is on lmyrank=0 of every group
      // all other procs have zero vectors of the same size here
      // communicate cvectors to gmyrank 0 and put into cmatrix
      gcomm->Barrier();
      if (spec_inv_ana_mult_)
      {
        int pos = 0;
        // go through all experiments
        for (int k = 0; k < nexp_; ++k)
        {
          // get writing start position of every experiment
          if (k > 0) pos += nmp_array_mult_[k - 1];

          for (int j = 0; j < np_ + 1; ++j)
          {
            int ssender = 0;
            int sender = 0;
            if (lmyrank == 0 && groupid == (j + k * (np_ + 1)))
            {
              ssender = gmyrank;
            }
            gcomm->SumAll(&ssender, &sender, 1);

            int n[2];
            n[0] = cvector.Length();
            n[1] = i;
            gcomm->Broadcast(n, 2, sender);
            Epetra_SerialDenseVector sendbuff(n[0]);
            if (sender == gmyrank)
            {
              // normalize to handle different experiments
              cvector.Scale(1. / facnorm_mult_[k]);
              sendbuff = cvector;
            }
            gcomm->Broadcast(sendbuff.Values(), n[0], sender);

            if (gmyrank == 0)  // gproc 0 puts received cvector in cmatrix
            {
              int ii = n[1];
              for (int l = pos; l < n[0] + pos; ++l)
              {
                cmatrix(l, ii) = sendbuff[l - pos];
              }
            }
          }
        }
      }
      // normal or coupled inverse analysis
      else
      {
        for (int j = 0; j < ngroup; ++j)
        {
          int ssender = 0;
          int sender = 0;
          if (lmyrank == 0 && groupid == j) ssender = gmyrank;
          gcomm->SumAll(&ssender, &sender, 1);

          int n[2];
          n[0] = cvector.Length();
          n[1] = i;
          gcomm->Broadcast(n, 2, sender);
          Epetra_SerialDenseVector sendbuff(n[0]);
          if (sender == gmyrank) sendbuff = cvector;
          gcomm->Broadcast(sendbuff.Values(), n[0], sender);
          if (gmyrank == 0)  // gproc 0 puts received cvector in cmatrix
          {
            if (nmp_ != n[0]) dserror("Size mismatch nmp_ %d != n %d", nmp_, n[0]);
            int ii = n[1];
            for (int k = 0; k < nmp_; ++k) cmatrix(k, ii) = sendbuff[k];
          }
        }
      }

      // reset discretization to blank
      Teuchos::ParameterList p;
      p.set("action", "calc_struct_reset_all");
      discret_->Evaluate(
          p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
    }  //--------------------------------------------------------------------------------

    gcomm->Barrier();
    if (!gmyrank) CalcNewParameters(cmatrix, perturb);

    // set new material parameters
    // these are the same for all groups
    gcomm->Broadcast(&p_[0], p_.Length(), 0);
    SetParameters(p_);

    numb_run_++;
    gcomm->Broadcast(&error_o_, 1, 0);
    gcomm->Broadcast(&error_, 1, 0);
    gcomm->Broadcast(&error_grad_o_, 1, 0);
    gcomm->Broadcast(&error_grad_, 1, 0);
    gcomm->Broadcast(&numb_run_, 1, 0);
    gcomm->Broadcast(&nodescentdirection_, 1, 0);

    if (reg_update_ == grad_based)
      error_i_ = error_grad_;
    else if (reg_update_ == res_based)
      error_i_ = error_;

  } while (error_i_ > tol_ && numb_run_ < max_itter && !nodescentdirection_);
  // printf("gmyrank %d reached this point\n",gmyrank); fflush(stdout); exit(0);


  // print results to file
  if (!gmyrank) PrintFile();

  return;
}

//---------------------------------------------------------------------------------------------
// Additional calculations of Integrate for Multiple Inverse Analysis
void STR::GenInvAnalysis::NPIntegrateMult()
{
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = DRT::Problem::Instance()->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  const int lmyrank = lcomm->MyPID();
  const int gmyrank = gcomm->MyPID();
  const int groupid = group->GroupId();
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();


  /* general idea: partition the NPGroup into subgroups, each subgroup belongs to one experiment;
   * then, distribute the single simulation runs (unpertubed and pertubed) in each subgroup to the
   * single procs this is currently achieved by the input syntax (script) of baci ngroup has to have
   * the size of (#experiments * (#params +1)) afterwards, send the results of the simulation runs
   * to one cmatrix on proc 0
   *
   * the cmatrix is arranged as followed
   * [u_pertubed_param1_exp1 --- u_pertubed_param2_exp1 --- ... --- u_unpertubed_exp1
   *  u_pertubed_param1_exp2 --- u_pertubed_param2_exp2 --- ... --- u_unpertubed_exp2
   *  ........
   *  u_pertubed_param1_expn --- u_pertubed_param2_expn --- ... --- u_unpertubed_expn]
   *
   *  this cmatrix is then forwarded to the CalcNewParametersMult() method
   *  and a new set of parameters is calculated
   */

  // get experiment ID of each optimization (starts from 0)
  int expid = iap.get<int>("EXP_ID_MULT_INV_ANA");

  // get number of experiments
  int expid_tmp_min = 0;
  int expid_tmp_max = 0;
  gcomm->Barrier();
  gcomm->MinAll(&expid, &expid_tmp_min, 1);
  gcomm->MaxAll(&expid, &expid_tmp_max, 1);

  if (expid_tmp_min == expid_tmp_max)
    // if only running one experiment, exp_id has just to be the same on all procs
    nexp_ = 1;
  else
    // if running multiple experiments
    nexp_ = expid_tmp_max + 1;  // calculate number of experiments

  // calculate number of measured dofs of all experiments
  // in case of multiple inverse analysis: nmp_mult_ stores the number of measured dofs of all
  // experiments
  gcomm->Barrier();
  gcomm->SumAll(&nmp_, &nmp_mult_, 1);
  nmp_mult_ = nmp_mult_ / (np_ + 1);  // Input is read for each identified parameter (+1)

  // This method checks for differences in the dat files which could disturb the run
  CheckDiffDat(expid);



  if (!gmyrank)
    // setting up a mcurve_mult_ vector with all measured curves of every experiment
    mcurve_mult_ = Epetra_SerialDenseVector(nmp_mult_);


  gcomm->Barrier();
  // this is an vector which stores all nmp_ values of every experiment:
  nmp_array_mult_ = Epetra_SerialDenseVector(nexp_);
  // this is an vector which stores all facnorms of every experiment:
  facnorm_mult_ = Epetra_SerialDenseVector(nexp_);
  int pos = 0;
  for (int k = 0; k < nexp_; ++k)
  {
    int ssender = 0;
    int sender = 0;
    if (lmyrank == 0 && groupid == k * (np_ + 1))
    {
      ssender = gmyrank;  // save first proc of each experiment

      nmp_array_mult_[k] = nmp_;  // save nmp of each experiment

      // calculate maximum and minimum measurement
      double min_zug = mcurve_(0);
      double max_zug = mcurve_(0);
      for (int n = 0; n < nmp_; n++)
      {
        min_zug = std::min(min_zug, mcurve_(n));
        max_zug = std::max(max_zug, mcurve_(n));
      }
      // norm of each experiment 1/sqrt(nmp)*(max-min)
      facnorm_mult_[k] = sqrt(nmp_) * (max_zug - min_zug);
    }
    gcomm->SumAll(&ssender, &sender, 1);
    gcomm->Broadcast(nmp_array_mult_.Values(), nexp_, sender);
    gcomm->Broadcast(facnorm_mult_.Values(), nexp_, sender);

    if (k > 0) pos += nmp_array_mult_[k - 1];  // save position of last experiment

    int n;
    n = mcurve_.Length();
    if (nmp_ != n) dserror("Size mismatch nmp_ %d != n %d", nmp_, n);
    gcomm->Broadcast(&n, 1, sender);
    Epetra_SerialDenseVector sendbuff(n);
    if (sender == gmyrank)
    {
      // normalize mcurve
      mcurve_.Scale(1. / facnorm_mult_[k]);

      sendbuff = mcurve_;
    }
    gcomm->Broadcast(sendbuff.Values(), n, sender);

    if (gmyrank ==
        0)  // gproc 0 puts received mcurve in mcurve_mult (one experiment after each other)
    {
      for (int l = pos; l < pos + n; ++l)
      {
        mcurve_mult_(l) = sendbuff[l - pos];
      }
    }
  }

  return;
}

//---------------------------------------------------------------------------------------------
void STR::GenInvAnalysis::CalcNewParameters(
    Epetra_SerialDenseMatrix& cmatrix, std::vector<double>& perturb)
{
  // in case of multiple inverse analysis: use measured dofs of all experiments
  int nmp = 0.0;
  Epetra_SerialDenseVector mcurve;
  if (spec_inv_ana_mult_)
  {
    nmp = nmp_mult_;  // overwrite nmp_
    mcurve = mcurve_mult_;
  }
  // normal or coupled inverse analysis
  else
  {
    nmp = nmp_;
    mcurve = mcurve_;
  }

  // --------------------------------------------------------------------
  // initialization of the Jacobian and other storage
  // --------------------------------------------------------------------
  // matrix to be inverted in each LVM iteration
  Epetra_SerialDenseMatrix sto(np_, np_);

  // parameter increment (to be solved for)
  Epetra_SerialDenseVector delta_p(np_);

  // right hand side for LVM iteration
  Epetra_SerialDenseVector rightHandSide(np_);

  // residuum (i.e. deviation of measurements and simulation results)
  Epetra_SerialDenseVector rcurve(nmp);

  // Simulation results
  Epetra_SerialDenseVector ccurve(nmp);

  // --------------------------------------------------------------------

  // copy column with unperturbed values to extra vector
  for (int i = 0; i < nmp; i++) ccurve[i] = cmatrix(i, np_);

  // remove the extra column np_+1 with unperturbed values
  cmatrix.Reshape(nmp, np_);

  // reuse the cmatrix array as J to save storage
  Epetra_SerialDenseMatrix& J = cmatrix;

  // in case of coupled inverse analysis
  if (spec_inv_ana_coup_)
  {
    // for information on coupled inverse analysis see
    // Anna M. Birzle, Christian Martin, Stefan Uhlig, Wolfgang A. Wall (2018):
    // A Nonlinear and Compressible Hyperelastic Material Model of Lung Parenchyma:
    // Experiments and Numerical Identification (in preperation)

    // get data
    // of volume-pressure-change experiment

    // number of data points of pressure-volume-change experiment
    nmp_volexp_ = 287;  // ToDo Why is this value hard-coded? This seems to be dangerous!
    // initialization of measured pressure and volume-change values
    std::vector<double> volexp_deltap(nmp_volexp_, 0.0);
    std::vector<double> volexp_deltaV(nmp_volexp_, 0.0);

    // get data
    pVData(volexp_deltap, volexp_deltaV);


    // compute pressure from volume-change and paramaters of respective material
    // this function dependents on the material, which is identified
    // this function is implemented for: Coupneohooke + Iso1Pow + Coup3Pow
    // extend new functions, if different materials should be identified

    // computed pressure values
    Epetra_SerialDenseMatrix volexp_p_comp(nmp_volexp_, np_ + 1);

    // compute pressure values
    ComputePressureNHIso1Coup3(perturb, volexp_deltaV, volexp_p_comp);


    // normalize data sets
    // with size of data set and range of data values

    // facnorm of pressure-volume-change experiment
    // get minimum and maximum of measured pressure values
    double min_volexp = volexp_deltap[0];
    double max_volexp = volexp_deltap[0];
    for (int n = 0; n < nmp_volexp_; n++)
    {
      min_volexp = std::min(min_volexp, volexp_deltap[n]);
      max_volexp = std::max(max_volexp, volexp_deltap[n]);
    }
    // calculate facnorm (N_vol)
    double volexp_fac_norm = sqrt(nmp_volexp_) * (max_volexp - min_volexp);

    // calculate facnorm (N_uniax) of tensile test
    double min_zug = mcurve(0);
    double max_zug = mcurve(0);
    for (int n = 0; n < nmp; n++)
    {
      min_zug = std::min(min_zug, mcurve(n));
      max_zug = std::max(max_zug, mcurve(n));
    }
    double fac_norm = sqrt((double)nmp) * (max_zug - min_zug);

    // calculating Jacobi-Matrix J(p)
    // tensile test data
    for (int i = 0; i < nmp; i++)
    {
      for (int j = 0; j < np_; j++)
      {
        J(i, j) -= ccurve[i];
        J(i, j) /= perturb[j];
        J(i, j) /= fac_norm;  // normalize J
      }
    }
    // add data of pressure-volume-change experiment
    J.Reshape(nmp + nmp_volexp_, np_);
    for (int i = 0; i < nmp_volexp_; i++)
    {
      for (int j = 0; j < np_; j++)
      {
        J(i + nmp, j) = volexp_p_comp(i, j);     // perturbated values of pressure
        J(i + nmp, j) -= volexp_p_comp(i, np_);  // - unperturbated values of pressure
        J(i + nmp, j) /= perturb[j];             // /perturbation
        J(i + nmp, j) /= volexp_fac_norm;        // normalize J
      }
    }

    // calculating residuum R
    // compute residual (measured vs. computed)
    // tensile test data (displacements)
    for (int i = 0; i < nmp; i++)
    {
      rcurve[i] = -mcurve[i] + ccurve[i];  // sign of residual
      rcurve[i] /= fac_norm;               // normalize residuum
    }
    // add data (pressure) of pressure-volume-change experiment
    rcurve.Reshape(nmp + nmp_volexp_, 1);
    for (int i = 0; i < nmp_volexp_; i++)
    {
      rcurve[i + nmp] = -volexp_deltap[i] + volexp_p_comp(i, np_);  // sign of residual
      rcurve[i + nmp] /= volexp_fac_norm;
    }
    // end of coupled inverse analysis
  }
  // normal or multiple inverse analysis
  else
  {
    nmp_volexp_ = 0;  // do not add any data
    // calculating J(p)
    for (int i = 0; i < nmp; i++)
      for (int j = 0; j < np_; j++)
      {
        J(i, j) -= ccurve[i];
        J(i, j) /= perturb[j];
        // normal inverse analysis: not normalized
        // multiple inverse analysis: already normalized
      }

    if (meas_type_ == dof_based)
    {
      // calculating R
      // compute residual displacement (measured vs. computed)
      for (int i = 0; i < nmp; i++) rcurve[i] = -mcurve[i] + ccurve[i];  // sign of residual
      // multiple inverse analysis: both parts already normalized
    }
    else if (meas_type_ == point_based)
    {
      // as the computed distances already are a residuum we use these
      for (int i = 0; i < nmp; i++) rcurve[i] = ccurve[i];
    }
    else
      dserror("This should not happen. You have to chose either point or dof based measurement.");
  }

  // calculating J.T*J
  sto.Multiply('T', 'N', 1.0, J, J, 0.0);

  // calculating  J.T*J + mu*diag(J.T*J)
  // do regularization by adding artificial mass on the main diagonal
  for (int i = 0; i < np_; i++) sto(i, i) += mu_ * sto(i, i);

  // delta_p = - (J.T*J+mu*diag(J.T*J))^{-1} * J.T*R
  rightHandSide.Multiply('T', 'N', -1.0, J, rcurve, 0.0);  // sign of update rule
  Epetra_SerialDenseSolver solver;
  solver.SetMatrix(sto);
  solver.SetVectors(delta_p, rightHandSide);
  solver.SolveToRefinedSolution(true);
  solver.Solve();

  // dependent on the # of steps
  error_grad_o_ = error_grad_;
  error_o_ = error_;
  mu_o_ = mu_;
  p_o_ = p_;

  // update res based error
  if ((spec_inv_ana_coup_) || (spec_inv_ana_mult_))
    error_ = rcurve.Norm2();  // already normalized
  else
    error_ = rcurve.Norm2() / sqrt((double)nmp);

  // get jacobian:
  Epetra_SerialDenseVector Ji(np_);
  for (int i = 0; i < np_; i++)
    for (int j = 0; j < nmp + nmp_volexp_; j++) Ji[i] += rcurve[j] * J(j, i);

  error_grad_ = Ji.Norm2();

  if (!numb_run_) error_grad_o_ = error_grad_;
  // Adjust training parameter based on gradient
  // check for a descent step in general

  // Gradient based update of mu based on
  // Kelley, C. T., Liao, L. Z., Qi, L., Chu, M. T., Reese, J. P., & Winton, C. (2009)
  // Projected pseudotransient continuation. SIAM Journal on Numerical Analysis, 46(6), 3071.
  /* ToDo Which of the updates formulas in this paper is used here? The paper presents three
   * distinct update formulas, denoted as SER-A, SER-B, and TTE.
   */
  if (reg_update_ == grad_based)
  {
    if (error_ < error_o_)
    {
      // update mu_ only if error in df/dp decreases
      if (error_grad_ < error_grad_o_) mu_ *= (error_grad_ / error_grad_o_);

      // update output of parameters
      p_print_ = p_;

      // update parameters
      for (int i = 0; i < np_; i++) p_[i] += delta_p[i];
    }
    else
    {
      std::cout << "WARNING: MAT Params not updated! No descent direction" << std::endl;
      nodescentdirection_ = 1;
    }
  }
  else
  // res_based update
  {
    // update output of parameters
    p_print_ = p_;
    // update params no matter what
    for (int i = 0; i < np_; i++) p_[i] += delta_p[i];
    if (!numb_run_) error_o_ = error_;

    // Adjust training parameter based on residuum of objective function ()
    mu_ *= (error_ / error_o_);
  }

  // return cmatrix to previous size and zero out
  cmatrix.Shape(nmp, np_ + 1);

  PrintStorage(delta_p);

  if (spec_inv_ana_coup_)
    CheckPhysiologicalPVRelationNHIso1Coup3();  // be careful, this function depends on the chosen
                                                // materials

  if (check_neg_params_) CheckOptStep();

  // the parameter alpha of the fractional viscoelasticity model needs to be smaller then 1
  // apply constrain just, if viscofract material is used (otherwise the position is by default
  // 9999)
  if (viscofract_alpha_position_ != 9999) ConstrainAlpha();

  return;
}

/*----------------------------------------------------------------------*/
/* */
Epetra_SerialDenseVector STR::GenInvAnalysis::CalcCvector(bool outputtofile)
{
  int myrank = discret_->Comm().MyPID();
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  Teuchos::ParameterList xparams;
  xparams.set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  xparams.set<int>("REDUCED_OUTPUT", 0);

  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  // use the same control file for every run since usually the last one is of interest
  structdis->Writer()->OverwriteResultFile();
  // create an adapterbase and adapter
  Teuchos::RCP<ADAPTER::Structure> structadapter = Teuchos::null;
  // FixMe The following switch is just a temporal hack, such we can jump between the new and the
  // old structure implementation. Has to be deleted after the clean-up has been finished!
  const enum INPAR::STR::IntegrationStrategy intstrat =
      DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");
  switch (intstrat)
  {
    // -------------------------------------------------------------------
    // old implementation
    // -------------------------------------------------------------------
    case INPAR::STR::int_old:
    {
      Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> adapterbase_old_ptr =
          Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
              sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
      structadapter = adapterbase_old_ptr->StructureField();
      structadapter->Setup();
      break;
    }
    // -------------------------------------------------------------------
    // new implementation
    // -------------------------------------------------------------------
    default:
    {
      Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase_ptr =
          ADAPTER::STR::BuildStructureAlgorithm(sdyn);
      adapterbase_ptr->Init(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
      adapterbase_ptr->Setup();
      structadapter = adapterbase_ptr->StructureField();
      break;
    }
  }


  int writestep = 0;
  Epetra_SerialDenseVector cvector(nmp_);

  // time loop
  while (structadapter->NotFinished())
  {
    // call the predictor
    structadapter->PrepareTimeStep();

    // integrate time step
    // after this step we hold disn_, etc
    structadapter->Solve();

    // calculate stresses, strains, energies
    structadapter->PrepareOutput();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    // update everything on the element level
    structadapter->Update();

    // print info about finished time step
    structadapter->PrintStep();

    // write output
    if (outputtofile) structadapter->Output();

    // get current time
    double time = structadapter->TimeOld();

    // get the displacements of the monitored timesteps
    {
      if (abs(time - timesteps_[writestep]) < 1.0e-5)
      {
        Epetra_SerialDenseVector cvector_arg = GetCalculatedCurve(*(structadapter->Dispnp()));
        if (!myrank)
          for (int j = 0; j < ndofs_; ++j) cvector[writestep * ndofs_ + j] = cvector_arg[j];
        writestep += 1;
      }

      // check if timestepsize is smaller than the tolerance above
      const double deltat = structadapter->Dt();
      if (deltat < 1.0e-5)
        dserror(
            "your time step size is too small, you will have problems with the monitored steps, "
            "thus adapt the tolerance");
    }
  }

  if (!(writestep * ndofs_ == nmp_))
  {
    std::cout << "# of timesteps extracted from the simulation: " << writestep << std::endl;
    dserror(
        "# of monitored timesteps does not match # of timesteps extracted from the simulation ");
  }

  return cvector;
}

/*----------------------------------------------------------------------*/
/* FSI Method                                               HaWi 06/18  */
/* mainly a reduced and adapted copy of fsi_ale_drt()                   */
/*----------------------------------------------------------------------*/
Epetra_SerialDenseVector STR::GenInvAnalysis::CalcCvectorFSI(bool outputtofile, int iLMstep)
{
  int myrank = discret_->Comm().MyPID();

  DRT::Problem* problem = DRT::Problem::Instance();

  problem->SetProblemType(forward_prb_type_);

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  const Epetra_Comm& comm = structdis->Comm();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       structure dof < fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  // look into problem to check for available discretizations
  //    std::cout << "///////XXXX/////// problem->GetDisNames()" << std::endl;
  //    const std::vector<std::string> disnames = problem->GetDisNames();
  //
  //    std::vector<std::string >::const_iterator iter;
  //    for (iter = disnames.begin(); iter != disnames.end(); ++iter)
  //    {
  //      std::cout << *iter << std::endl;
  //    }

  // use the same control file for every run since usually the last one is of interest
  structdis->Writer()->OverwriteResultFile();

  structdis->FillComplete();

  problem->GetDis("fluid")->FillComplete();

  problem->GetDis("ale")->FillComplete();

  // get discretizations
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

  // i guess it was not meant to be, but it works. suggestions for improvement are welcome
  // without these lines multiple evaluation of this FSI method does produce an error,
  // when you ask for binary output, which usually isn't of great interest for inverse analysis
  fluiddis->Writer()->OverwriteResultFile();
  aledis->Writer()->OverwriteResultFile();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes() == 0)  // empty ale discretization
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis, aledis);
    aledis->FillComplete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else if (iLMstep != 0 or numb_run_ != 0)
  {
    // here we do not need a new ale discretization as it was generated in the first forward run
  }
  else  // filled ale discretization (i.e. read from input file)
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis, aledis))
      dserror(
          "Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");

    if ((not DRT::INPUT::IntegralValue<bool>(problem->FSIDynamicParams(), "MATCHGRID_FLUIDALE")) or
        (not DRT::INPUT::IntegralValue<bool>(problem->FSIDynamicParams(), "MATCHGRID_STRUCTALE")))
    {
      // create vector of discr.
      std::vector<Teuchos::RCP<DRT::Discretization>> dis;
      dis.push_back(structdis);
      dis.push_back(fluiddis);
      dis.push_back(aledis);

      std::vector<Teuchos::RCP<Epetra_Map>> stdelecolmap;
      std::vector<Teuchos::RCP<Epetra_Map>> stdnodecolmap;

      // redistribute discr. with help of binning strategy
      if (structdis->Comm().NumProc() > 1)
      {
        // binning strategy is created and parallel redistribution is performed
        Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
            Teuchos::rcp(new BINSTRATEGY::BinningStrategy());
        binningstrategy->Init(dis);
        binningstrategy->DoWeightedPartitioningOfBinsAndExtendGhostingOfDiscretToOneBinLayer(
            dis, stdelecolmap, stdnodecolmap);
      }
    }
  }

  // invana vector
  int writestep = 0;
  Epetra_SerialDenseVector cvector(nmp_);

  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  FSI_COUPLING coupling = DRT::INPUT::IntegralValue<FSI_COUPLING>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_pseudo_structureale:
    case fsi_iter_fluidfluid_monolithicfluidsplit_nonox:
    case fsi_iter_fluidfluid_monolithicstructuresplit_nonox:
    {
      dserror("FSI coupling algorithm not covered for inverse analysis");
      break;
    }
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
    case fsi_iter_lung_monolithicstructuresplit:
    case fsi_iter_lung_monolithicfluidsplit:
    case fsi_iter_constr_monolithicfluidsplit:
    case fsi_iter_constr_monolithicstructuresplit:
    case fsi_iter_mortar_monolithicstructuresplit:
    case fsi_iter_mortar_monolithicfluidsplit:
    case fsi_iter_fluidfluid_monolithicfluidsplit:
    case fsi_iter_fluidfluid_monolithicstructuresplit:
    case fsi_iter_sliding_monolithicfluidsplit:
    case fsi_iter_sliding_monolithicstructuresplit:
    {
      // monolithic solver settings
      const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

      Teuchos::RCP<FSI::Monolithic> fsi;

      INPAR::FSI::LinearBlockSolver linearsolverstrategy =
          DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

      // call constructor to initialize the base class
      if (coupling == fsi_iter_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::MonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_lung_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::LungMonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_lung_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::LungMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_constr_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::ConstrMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_constr_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::ConstrMonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_mortar_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::MortarMonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_mortar_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::MortarMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicStructureSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_sliding_monolithicfluidsplit)
      {
        fsi = Teuchos::rcp(new FSI::SlidingMonolithicFluidSplit(comm, fsidyn));
      }
      else if (coupling == fsi_iter_sliding_monolithicstructuresplit)
      {
        fsi = Teuchos::rcp(new FSI::SlidingMonolithicStructureSplit(comm, fsidyn));
      }
      else
      {
        dserror("Cannot find appropriate monolithic solver for coupling %d and linear strategy %d",
            coupling, linearsolverstrategy);
      }

      // read the restart information, set vectors and variables ---
      // be careful, dofmaps might be changed here in a Redistribute call
      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        fsi->ReadRestart(restart);
      }

      // now do the coupling setup and create the combined dofmap
      fsi->SetupSystem();

      // possibly redistribute domain decomposition
      {
        const INPAR::FSI::Redistribute redistribute =
            DRT::INPUT::IntegralValue<INPAR::FSI::Redistribute>(fsimono, "REDISTRIBUTE");

        const double weight1 = fsimono.get<double>("REDIST_WEIGHT1");
        const double weight2 = fsimono.get<double>("REDIST_WEIGHT2");
        if (redistribute == INPAR::FSI::Redistribute_structure or
            redistribute == INPAR::FSI::Redistribute_fluid)
        {
          // redistribute either structure or fluid domain
          fsi->RedistributeDomainDecomposition(redistribute, coupling, weight1, weight2, comm, 0);

          // do setup again after redistribution
          fsi->SetupSystem();

          const double secweight1 =
              fsidyn.sublist("MONOLITHIC SOLVER").get<double>("REDIST_SECONDWEIGHT1");
          const double secweight2 =
              fsidyn.sublist("MONOLITHIC SOLVER").get<double>("REDIST_SECONDWEIGHT2");
          if (secweight1 != -1.0)
          {
            // redistribute either structure or fluid domain
            fsi->RedistributeDomainDecomposition(
                redistribute, coupling, secweight1, secweight2, comm, 0);

            // do setup again after redistribution
            fsi->SetupSystem();
          }
        }
        else if (redistribute == INPAR::FSI::Redistribute_both)
        {
          int numproc = comm.NumProc();

          // redistribute both structure and fluid domain
          fsi->RedistributeDomainDecomposition(
              INPAR::FSI::Redistribute_structure, coupling, weight1, weight2, comm, numproc / 2);

          // do setup again after redistribution (do this again here in between because the P matrix
          // changed!)
          fsi->SetupSystem();
          fsi->RedistributeDomainDecomposition(
              INPAR::FSI::Redistribute_fluid, coupling, weight1, weight2, comm, numproc / 2);

          // do setup again after redistribution
          fsi->SetupSystem();
        }
        else if (redistribute == INPAR::FSI::Redistribute_monolithic)
        {
          fsi->RedistributeMonolithicGraph(coupling, comm);

          // do setup again after redistribution
          fsi->SetupSystem();
        }
      }

      // here we go...
      //        fsi->Timeloop(fsi);, which is split up here to extract nodal displacements for
      //        inverse analysis
      fsi->PrepareTimeloop();

      while (fsi->NotFinished())
      {
        fsi->PrepareTimeStep();
        fsi->TimeStep(fsi);
        fsi->PrepareOutput();
        fsi->Update();
        if (outputtofile) fsi->Output();

        // get current time
        double time = fsi->Time();

        // get the displacements of the monitored timesteps
        {
          if (writestep >= nsteps_)
            dserror("you should not proceed in time after last monitored step");
          // because it is useless and could lead to a segfault here

          if (abs(time - timesteps_[writestep]) < 1.0e-5)
          {
            switch (meas_type_)
            {
              case dof_based:
              {
                Epetra_SerialDenseVector cvector_arg =
                    GetCalculatedCurve(*fsi->StructureField()->Dispnp());
                if (!myrank)
                  for (int j = 0; j < ndofs_; ++j) cvector[writestep * ndofs_ + j] = cvector_arg[j];
                break;
              }
              case point_based:
              {
                Epetra_SerialDenseVector dvector;
                // get colmap for parallel execution
                Teuchos::RCP<Epetra_Vector> dispnpcol =
                    Teuchos::rcp(new Epetra_Vector(*discret_->DofColMap(), true));
                LINALG::Export(*fsi->StructureField()->Dispnp(), *dispnpcol);

                dvector = GetDistancePointsToInterfaceContour(*dispnpcol, writestep);

                for (int j = 0; j < nnodes_; ++j) cvector[writestep * nnodes_ + j] = dvector[j];
                break;
              }
              default:
                dserror("This should never happen. Did you forget to implement a new case?");
                break;
            }
            writestep += 1;
          }

          // check if timestepsize is smaller than the tolerance above
          const double deltat = fsi->Dt();
          if (deltat < 1.0e-5)
            dserror(
                "your time step size is too small, you will have problems with the monitored "
                "steps, thus adapt the tolerance");
        }
      }

      break;
    }
    default:
    {
      dserror("FSI coupling algorithm not covered for inverse analysis");
      break;
    }
  }

  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);

  if ((meas_type_ == dof_based and !(writestep * ndofs_ == nmp_)) or
      (meas_type_ == point_based and !(writestep * nnodes_ == nmp_)))
  {
    std::cout << "# of timesteps extracted from the simulation: " << writestep << std::endl;
    dserror(
        "# of monitored timesteps does not match # of timesteps extracted from the simulation ");
  }
  problem->SetProblemType(prb_invana);
  return cvector;
}
/*----------------------------------------------------------------------*/
/* nested parallelity version of the same method            mwgee 05/12 */
Epetra_SerialDenseVector STR::GenInvAnalysis::CalcCvector(
    bool outputtofile, Teuchos::RCP<COMM_UTILS::NestedParGroup> group)
{
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  const int lmyrank = lcomm->MyPID();
  // const int gmyrank  = gcomm->MyPID();
  // const int groupid  = group->GroupId();
  // const int ngroup   = group->NumGroups();

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  Teuchos::ParameterList xparams;
  xparams.set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  xparams.set<int>("REDUCED_OUTPUT", 0);


  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");
  // use the same control file for every run since usually the last one is of interest
  structdis->Writer()->OverwriteResultFile();
  // create an adapterbase and adapter
  Teuchos::RCP<ADAPTER::Structure> structadapter = Teuchos::null;
  // FixMe The following switch is just a temporal hack, such we can jump between the new and the
  // old structure implementation. Has to be deleted after the clean-up has been finished!
  const enum INPAR::STR::IntegrationStrategy intstrat =
      DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");
  switch (intstrat)
  {
    // -------------------------------------------------------------------
    // old implementation
    // -------------------------------------------------------------------
    case INPAR::STR::int_old:
    {
      Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> adapterbase_old_ptr =
          Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(
              sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
      structadapter = adapterbase_old_ptr->StructureField();
      structadapter->Setup();
      break;
    }
    // -------------------------------------------------------------------
    // new implementation
    // -------------------------------------------------------------------
    default:
    {
      Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase_ptr =
          ADAPTER::STR::BuildStructureAlgorithm(sdyn);
      adapterbase_ptr->Init(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), structdis);
      adapterbase_ptr->Setup();
      structadapter = adapterbase_ptr->StructureField();
      break;
    }
  }

  if (structadapter == Teuchos::null) dserror("Failed in creating integrator.");
  fflush(stdout);
  gcomm->Barrier();

  int writestep = 0;
  Epetra_SerialDenseVector cvector(nmp_);

  // time loop
  while (structadapter->NotFinished())
  {
    // call the predictor
    structadapter->PrepareTimeStep();

    // integrate time step
    // after this step we hold disn_, etc
    structadapter->Solve();

    // calculate stresses, strains, energies
    structadapter->PrepareOutput();

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    // update everything on the element level
    structadapter->Update();

    // print info about finished time step
    structadapter->PrintStep();

    // write output
    if (outputtofile) structadapter->Output();

    // get current time
    double time = structadapter->TimeOld();

    // get the displacements of the monitored timesteps
    {
      // if (abs(time-timesteps_[writestep]) < 1.0e-5)
      if (abs(time - timesteps_[writestep]) < 1.0e-4)
      {
        Epetra_SerialDenseVector cvector_arg = GetCalculatedCurve(*(structadapter->Dispnp()));
        if (!lmyrank)
          for (int j = 0; j < ndofs_; ++j) cvector[writestep * ndofs_ + j] = cvector_arg[j];
        writestep += 1;
      }

      // check if timestepsize is smaller than the tolerance above
      const double deltat = structadapter->Dt();
      if (deltat < 1.0e-5)
        dserror(
            "your time step size is too small, you will have problems with the monitored steps, "
            "thus adapt the tolerance");
    }
  }

  if (!(writestep * ndofs_ == nmp_))
  {
    std::cout << "writestep: " << writestep << std::endl;
    dserror(
        "# of monitored timesteps does not match # of timesteps extracted from the simulation ");
  }

  return cvector;
}

/*----------------------------------------------------------------------*/
/* */
Epetra_SerialDenseVector STR::GenInvAnalysis::GetCalculatedCurve(const Epetra_Vector& disp)
{
  int myrank = discret_->Comm().MyPID();

  // values observed at current time step
  Epetra_SerialDenseVector lcvector_arg(ndofs_);
  Epetra_SerialDenseVector gcvector_arg(ndofs_);
  int count = 0;
  for (int i = 0; i < nnodes_; ++i)
  {
    int gnode = nodes_[i];
    if (!discret_->NodeRowMap()->MyGID(gnode))
    {
      count += (int)dofs_[i].size();
      continue;
    }
    DRT::Node* node = discret_->gNode(gnode);
    if (!node) dserror("Cannot find node on this proc");
    if (myrank != node->Owner())
    {
      count += (int)dofs_[i].size();
      continue;
    }
    for (int j = 0; j < (int)dofs_[i].size(); ++j)
    {
      int ldof = dofs_[i][j];
      int gdof = discret_->Dof(node, ldof);
      if (!disp.Map().MyGID(gdof)) dserror("Cannot find dof on this proc");
      lcvector_arg[count + j] = disp[disp.Map().LID(gdof)];
    }
    count += (int)dofs_[i].size();
  }
  discret_->Comm().SumAll(&lcvector_arg[0], &gcvector_arg[0], ndofs_);

  return gcvector_arg;
}

/*----------------------------------------------------------------------*/
/* measure distance between points and deformed interface    HaWi 08/18 */
/*----------------------------------------------------------------------*/

Epetra_SerialDenseVector STR::GenInvAnalysis::GetDistancePointsToInterfaceContour(
    const Epetra_Vector& disp, int writestep)
{
  //  Epetra_SerialDenseVector lcvector_arg(nnodes_);
  Epetra_SerialDenseVector dvector(nnodes_);

  std::vector<DRT::Condition*> fsicond;
  discret_->GetCondition("FSICoupling", fsicond);

  if (fsicond.size() != 1)
    dserror("You are not supposed to be here, FSI condition is not unique for structure");

  std::map<int, Teuchos::RCP<DRT::Element>>& interfacegeom = fsicond[0]->Geometry();

  std::map<int, Teuchos::RCP<DRT::Element>>::iterator curr;
  DRT::Element::LocationArray la(discret_->NumDofSets());
  if (discret_->NumDofSets() != 1) dserror("method was developed for one dofset (structure) only!");

  int idofs = 0;

  // nnodes_ is the number of measurement points here
  for (int i = 0; i < nnodes_; i++)
  {
    // store point based input data in a matrix
    LINALG::Matrix<3, 2> line;
    line.PutScalar(0);
    for (unsigned int j = 0; j < dofs_[i].size(); j++)
    {
      line(dofs_[i][j], 0) = mcurve_[writestep * ndofs_ + idofs * 2 + j];
      line(dofs_[i][j], 1) = mcurve_[writestep * ndofs_ + idofs * 2 + dofs_[i].size() + j];
    }
    idofs += dofs_[i].size();

    // look for intersection of line with interface
    std::vector<std::vector<double>> xsifound(0);

    for (curr = interfacegeom.begin(); curr != interfacegeom.end(); ++curr)
    {
      if (curr->second->Owner() != discret_->Comm().MyPID()) continue;

      curr->second->LocationVector(*discret_, la, false);

      DRT::Node** ifacenodes = curr->second->Nodes();
      switch (curr->second->Shape())
      {
        case DRT::Element::tri3:
        {
          LINALG::Matrix<3, 1> xsi;
          GEO::CUT::KERNEL::ComputeIntersection<3, DRT::Element::line2, DRT::Element::tri3> ci(xsi);
          break;
        }
        case DRT::Element::line2:
        {
          LINALG::Matrix<2, 2> linee;
          for (unsigned int k = 0; k < 2; k++)
          {
            for (unsigned int l = 0; l < 2; l++)
            {
              linee(k, l) = line(k, l);
            }
          }
          LINALG::Matrix<2, 1> xsi;
          GEO::CUT::KERNEL::ComputeIntersection<2, DRT::Element::line2, DRT::Element::line2> ci(
              xsi);

          LINALG::Matrix<2, 2> surf;
          surf.PutScalar(0);

          std::vector<double> eledisp(4);
          DRT::UTILS::ExtractMyValues(disp, eledisp, la[0].lm_);

          for (int i = 0; i < 2; i++)
          {
            const double* x = ifacenodes[i]->X();
            surf(0, i) = x[0];
            surf(1, i) = x[1];

            surf(0, i) += eledisp[i * 2 + 0];
            surf(1, i) += eledisp[i * 2 + 1];
          }

          std::vector<double> initxsifound(2);
          initxsifound[0] = std::numeric_limits<double>::max();
          initxsifound[1] = std::numeric_limits<double>::max();
          xsifound.push_back(initxsifound);

          if (ci(surf, linee))
          {
            // intersection converged
            if (ci.SurfaceWithinLimits())
            {
              std::vector<double> myxsi(2);
              for (int m = 0; m < 2; m++)
              {
                myxsi[m] = xsi(m);
              }
              // store all intersection points
              xsifound.push_back(myxsi);
            }
          }
          else
          {
            break;
          }
          break;
        }
        default:
          dserror("surface Element for contour condition not yet implemented for FSI");
      }
    }

    // decide which element is the desired one and compute distance
    // this interferes with the structure of your monitor file
    // large dummy number, so that another proc will find a lower xsi
    double minxsi = std::numeric_limits<double>::max();
    int lastxsi = xsifound[0].size() - 1;
    for (unsigned int j = 0; j < xsifound.size(); j++)
    {
      if (j == 0)
        minxsi = xsifound[j][lastxsi];
      else if (xsifound[j][lastxsi] < minxsi)
        minxsi = xsifound[j][lastxsi];
      else
      {
        // xsifound doesn't contain the smallest xsi, do nothing
      }
    }

    // find minimal xsi over all processors
    double lminxsi = minxsi;
    discret_->Comm().MinAll(&lminxsi, &minxsi, 1);

    double distance = 0.;
    for (unsigned int j = 0; j < 3; j++)
    {
      distance += pow((minxsi + 1) / 2 * (line(j, 1) - line(j, 0)), 2);
    }

    // build LM-residual vector from distances
    distance = sqrt(distance);
    distance *= (minxsi + 1 > 0) ? 1 : ((minxsi + 1 < 0) ? -1 : 0);
    dvector[i] = distance;
  }
  return dvector;
}

/*----------------------------------------------------------------------*/
/* only global proc 0 comes in here! mwgee 5/2012 */
void STR::GenInvAnalysis::PrintStorage(Epetra_SerialDenseVector delta_p)
{
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = DRT::Problem::Instance()->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  const int gmyrank = group->GlobalComm()->MyPID();
  if (gmyrank != 0)
    dserror("Only gmyrank=0 is supposed to be in here, but is not, gmyrank=%d", gmyrank);

  // storing current material parameters
  // storing current material parameter increments
  p_s_.Reshape(numb_run_ + 1, np_);
  delta_p_s_.Reshape(numb_run_ + 1, np_);
  for (int i = 0; i < np_; i++)
  {
    p_s_(numb_run_, i) = p_(i);
    delta_p_s_(numb_run_, i) = delta_p(i);
  }

  // This routine runs just for material coupexppol (birzle 08/2014)
  // Recalculates the material parameters b and c for (and just for) printout
  // Inverse Analysis is done with ln(b) and ln(c) --> print exp(b) and exp(c)
  /* ToDo is the above comment really true? I don't see any guard that prevents this routine to run
   * also in other cases than for material coupexppol.
   */
  for (unsigned prob = 0; prob < DRT::Problem::NumInstances(); ++prob)
  {
    const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
        *DRT::Problem::Instance(prob)->Materials()->Map();
    std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;
    // loop all materials in problem
    for (curr = mats.begin(); curr != mats.end(); ++curr)
    {
      const Teuchos::RCP<MAT::PAR::Material> actmat = curr->second;
      if (actmat->Type() == INPAR::MAT::mes_coupexppol)
      {
        int j = numb_run_;
        p_s_(j, coupexppol_parameter_position_) = exp(p_s_(j, coupexppol_parameter_position_));
        p_s_(j, coupexppol_parameter_position_ + 1) =
            exp(p_s_(j, coupexppol_parameter_position_ + 1));
      }
    }
  }

  // this memory is going to explode, do we really need this? mwgee
  // ccurve_s_.Reshape(nmp_,  numb_run_+1);
  // for (int i=0; i<nmp_; i++)
  // ccurve_s_(i, numb_run_)= cmatrix(i, cmatrix.ColDim()-1);

  mu_s_.Resize(numb_run_ + 1);
  mu_s_(numb_run_) = mu_;

  error_s_.Resize(numb_run_ + 1);
  error_s_(numb_run_) = error_;
  error_grad_s_.Resize(numb_run_ + 1);
  error_grad_s_(numb_run_) = error_grad_;

  // extended output of inverse analysis (birzle 12/2016)
  // in last line: error and corresponding material parameters of current best fit
  p_s_.Reshape(numb_run_ + 2, np_);
  error_s_.Resize(numb_run_ + 2);
  error_grad_s_.Resize(numb_run_ + 2);
  // if parameters are not updated, store parameters and error with lowest error
  // (=second last parameter) in extra last line
  if (nodescentdirection_)
  {
    for (int i = 0; i < np_; i++) p_s_(numb_run_ + 1, i) = p_print_(i);
    error_s_(numb_run_ + 1) = error_s_(numb_run_ - 1);
    error_grad_s_(numb_run_ + 1) = error_grad_s_(numb_run_ - 1);
  }
  // if graderror reached minimum or max number of iterations is reached
  // or during inverse analysis
  else
  {
    for (int i = 0; i < np_; i++) p_s_(numb_run_ + 1, i) = p_print_(i);
    error_s_(numb_run_ + 1) = error_s_(numb_run_);
    error_grad_s_(numb_run_ + 1) = error_grad_s_(numb_run_);
  }


  // print error and parameter
  if (gmyrank == 0)  // this if should actually not be necessary since there is only gproc 0 in here
  {
    std::cout << std::endl;
    printf("################################################");
    printf("##############################################\n");
    printf("############################ Inverse Analysis ##");
    printf("##############################################\n");
    printf("################################### run ########");
    printf("##############################################\n");
    printf("################################### %3i ########", numb_run_);
    printf("##############################################\n");
    printf("################################################");
    printf("##############################################\n");

    for (int i = 0; i < numb_run_ + 1; i++)
    {
      printf("Error: ");
      printf("%10.3e", error_s_(i));
      printf("\tGrad_error: ");
      printf("%10.3e", error_grad_s_(i));
      printf("\tParameter: ");
      for (int j = 0; j < delta_p.Length(); j++) printf("%10.3e\t", p_s_(i, j));
      // printf("\tDelta_p: ");
      // for (int j=0; j < delta_p.Length(); j++)
      //  printf("%10.3e", delta_p_s_(i, j));
      printf("\tmu: ");
      printf("%10.3e", mu_s_(i));
      printf("\n");
    }

    // print final parameters and error with lowest error
    // in extra last line
    int i = numb_run_ + 1;
    printf("Error and parameters with lowest error:\n");
    printf("Error: ");
    printf("%10.3e", error_s_(i));
    printf("\tGrad_error: ");
    printf("%10.3e", error_grad_s_(i));
    printf("\tParameter: ");
    for (int j = 0; j < delta_p.Length(); j++) printf("%10.3e\t", p_s_(i, j));
    printf("\n");

    printf("\n");
    fflush(stdout);
  }

  return;
}


/*----------------------------------------------------------------------*/
/* only gmyrank=0 comes in here */
void STR::GenInvAnalysis::PrintFile()
{
  // FILE * cxFile;
  // FILE * cyFile;
  FILE* pFile;

  std::string name = DRT::Problem::Instance()->OutputControlFile()->FileName();
  name.append(filename_);

  if (name.rfind("_run_") != std::string::npos)
  {
    size_t pos = name.rfind("_run_");
    if (pos == std::string::npos) dserror("inconsistent file name");
    name = name.substr(0, pos);
  }

  std::string gp = name + "_plot.gp";
  std::string xcurve = name + "_Curve_x.txt";
  std::string ycurve = name + "_Curve_y.txt";
  std::string para = name + "_Para.txt";

#if 0  // this only produces columns of zeros anyway
  // it will also burst memory with ccurve_s_ pretty quickly for larger problems
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
#endif

  pFile = fopen((para).c_str(), "w");
  fprintf(pFile, "#Error    #Grad_error       Parameter    Delta_p      mu \n");
  for (int i = 0; i < numb_run_; i++)
  {
    fprintf(pFile, "%10.6f, ", error_s_(i));
    fprintf(pFile, "%10.6f, ", error_grad_s_(i));
    for (int j = 0; j < np_; j++) fprintf(pFile, "%10.6f, ", p_s_(i, j));
    for (int j = 0; j < np_; j++) fprintf(pFile, "%10.6f, ", delta_p_s_(i, j));
    fprintf(pFile, "%10.6f", mu_s_(i));
    fprintf(pFile, "\n");
  }
  // print final parameters and error with lowest error
  // in extra last line
  int i = numb_run_;
  fprintf(pFile, "%10.6f, ", error_s_(i));
  fprintf(pFile, "%10.6f, ", error_grad_s_(i));
  for (int j = 0; j < np_; j++) fprintf(pFile, "%10.6f, ", p_s_(i, j));
  fprintf(pFile, "\n");

  fclose(pFile);

  numb_run_ = numb_run_ - 1;
}


//-----------------------------------------------------------------------------------
void STR::GenInvAnalysis::ReadInParameters()
{
  const int myrank = discret_->Comm().MyPID();

  for (unsigned prob = 0; prob < DRT::Problem::NumInstances(); ++prob)
  {
    const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
        *DRT::Problem::Instance(prob)->Materials()->Map();
    std::set<int> mymatset = matset_[prob];

    unsigned int overallnummat = mats.size();
    if (mymatset.size() > 0 && overallnummat > mymatset.size()) overallnummat = mymatset.size();
    if (myrank == 0) printf("No. material laws/summands considered : %d\n", overallnummat);

    std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;

    // default value of alpha position (to check, if alpha is fitted)
    viscofract_alpha_position_ = 9999;

    for (curr = mats.begin(); curr != mats.end(); ++curr)
    {
      const Teuchos::RCP<MAT::PAR::Material> actmat = curr->second;

      switch (actmat->Type())
      {
        case INPAR::MAT::m_aaaneohooke:
        {
          if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
          {
            MAT::PAR::AAAneohooke* params =
                dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
            if (!params) dserror("Cannot cast material parameters");
            const int j = p_.Length();
            p_.Resize(j + 2);
            p_(j) = params->GetParameter(params->young, -1);
            p_(j + 1) = params->GetParameter(params->beta, -1);
            // p_(j+2) = params->nue_; // need also change resize above to invoke
            // nue
          }
        }
        break;
        case INPAR::MAT::m_neohooke:
        {
          if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
          {
            MAT::PAR::NeoHooke* params = dynamic_cast<MAT::PAR::NeoHooke*>(actmat->Parameter());
            if (!params) dserror("Cannot cast material parameters");
            const int j = p_.Length();
            p_.Resize(j + 2);
            p_(j) = params->youngs_;
            p_(j + 1) = params->poissonratio_;
          }
        }
        break;
        case INPAR::MAT::m_stvenant:
        {
          if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
          {
            MAT::PAR::StVenantKirchhoff* params =
                dynamic_cast<MAT::PAR::StVenantKirchhoff*>(actmat->Parameter());
            if (!params) dserror("Cannot cast material parameters");
            const int j = p_.Length();
            p_.Resize(j + 2);
            p_(j) = params->youngs_;
            p_(j + 1) = params->poissonratio_;
          }
        }
        break;
        case INPAR::MAT::m_constraintmixture:
        {
          materialhashistory_ = true;
          if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
          {
            MAT::PAR::ConstraintMixture* params =
                dynamic_cast<MAT::PAR::ConstraintMixture*>(actmat->Parameter());
            if (!params) dserror("Cannot cast material parameters");
            const int j = p_.Length();
            p_.Resize(j + 1);
            // p_(j)   = params->mue_;
            // p_(j+1) = params->k1_;
            // p_(j+2) = params->k2_;
            // p_(j) = params->prestretchcollagen_;
            p_(j) = params->GetParameter(params->growthfactor, -1);
          }
        }
        break;
        case INPAR::MAT::m_growth_iso_stress:
        {
          materialhashistory_ = true;
          if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
          {
            MAT::PAR::GrowthLawIsoStress* params =
                dynamic_cast<MAT::PAR::GrowthLawIsoStress*>(actmat->Parameter());
            if (!params) dserror("Cannot cast material parameters");
            const int j = p_.Length();
            p_.Resize(j + 1);
            p_(j) = params->kthetaplus_;
          }
        }
        break;
        case INPAR::MAT::m_elasthyper:
        {
          MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
          if (!params) dserror("Cannot cast material parameters");
          const int nummat = params->nummat_;
          const std::vector<int>* matids = params->matids_;
          for (int i = 0; i < nummat; ++i)
          {
            const int id = (*matids)[i];

            if (mymatset.size() == 0 or mymatset.find(id) != mymatset.end())
            {
              const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
              switch (actelastmat->Type())
              {
                case INPAR::MAT::mes_couplogneohooke:
                {
                  filename_ = filename_ + "_couplogneohooke";
                  const MAT::ELASTIC::PAR::CoupLogNeoHooke* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupLogNeoHooke*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->mue_;
                  p_[j + 1] = params2->lambda_;
                  break;
                }
                case INPAR::MAT::mes_coupexppol:
                {
                  filename_ = filename_ + "_coupexppol";
                  const MAT::ELASTIC::PAR::CoupExpPol* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupExpPol*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 3);
                  p_[j] = params2->a_;
                  p_[j + 1] = params2->b_;
                  p_[j + 2] = params2->c_;

                  // do inverse analysis with ln(b) and ln(c)
                  p_[j + 1] = log(p_[j + 1]);
                  p_[j + 2] = log(p_[j + 2]);
                  // remind position of b in p_
                  coupexppol_parameter_position_ = j + 1;
                  break;
                }
                case INPAR::MAT::mes_coupneohooke:
                {
                  filename_ = filename_ + "_coupneohooke";
                  const MAT::ELASTIC::PAR::CoupNeoHooke* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupNeoHooke*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->youngs_;
                  p_[j + 1] = params2->nue_;
                  break;
                }
                case INPAR::MAT::mes_coupblatzko:
                {
                  filename_ = filename_ + "_coupblatzko";
                  const MAT::ELASTIC::PAR::CoupBlatzKo* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupBlatzKo*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->mue_;
                  p_[j + 1] = params2->nue_ / (1. - 2. * params2->nue_);
                  // p_[j+2] = params2->f_;
                  break;
                }
                case INPAR::MAT::mes_coupsimopister:
                {
                  filename_ = filename_ + "_coupsimopister";
                  const MAT::ELASTIC::PAR::CoupSimoPister* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupSimoPister*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->mue_;
                  break;
                }
                case INPAR::MAT::mes_coupSVK:
                {  // Formulation with lambda and mue
                  /*
                filename_ = filename_ + "_coupsaintvenantkirchhoff";
                const MAT::ELASTIC::PAR::CoupSVK* params2 =
                    dynamic_cast<const MAT::ELASTIC::PAR::CoupSVK*>(actelastmat->Parameter());
                int j = p_.Length();
                p_.Resize(j + 2);
                p_[j] = params2->lambda_;
                p_[j + 1] = params2->mue_;*/

                  // Formulation with youngs and nue
                  filename_ = filename_ + "_coupSVK";
                  const MAT::ELASTIC::PAR::CoupSVK* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupSVK*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->mue_ * (3 * params2->lambda_ + 2 * params2->mue_) /
                          (params2->lambda_ + params2->mue_);
                  p_[j + 1] = params2->lambda_ / (2 * params2->mue_ + 2 * params2->lambda_);
                  break;
                }

                case INPAR::MAT::mes_isoneohooke:
                {
                  filename_ = filename_ + "_isoneohooke";
                  const MAT::ELASTIC::PAR::IsoNeoHooke* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoNeoHooke*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->mue_;
                  break;
                }
                case INPAR::MAT::mes_isoyeoh:
                {
                  filename_ = filename_ + "_isoyeoh";
                  const MAT::ELASTIC::PAR::IsoYeoh* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 3);
                  p_[j] = params2->c1_;
                  p_[j + 1] = params2->c2_;
                  p_[j + 2] = params2->c3_;
                  break;
                }
                case INPAR::MAT::mes_iso1pow:
                {
                  filename_ = filename_ + "_iso1pow";
                  const MAT::ELASTIC::PAR::Iso1Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_iso2pow:
                {
                  filename_ = filename_ + "_iso2pow";
                  const MAT::ELASTIC::PAR::Iso2Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Iso2Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_coup1pow:
                {
                  filename_ = filename_ + "_coup1pow";
                  const MAT::ELASTIC::PAR::Coup1Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Coup1Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_coup2pow:
                {
                  filename_ = filename_ + "_coup2pow";
                  const MAT::ELASTIC::PAR::Coup2Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Coup2Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_coup3pow:
                {
                  filename_ = filename_ + "_coup3pow";
                  const MAT::ELASTIC::PAR::Coup3Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Coup3Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_coup13apow:
                {
                  filename_ = filename_ + "_coup13apow";
                  const MAT::ELASTIC::PAR::Coup13aPow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Coup13aPow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->c_;
                  p_[j + 1] = params2->a_;
                  break;
                }
                case INPAR::MAT::mes_coupmooneyrivlin:
                {
                  filename_ = filename_ + "_coupmooneyrivlin";
                  const MAT::ELASTIC::PAR::CoupMooneyRivlin* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupMooneyRivlin*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 3);
                  p_[j] = params2->c1_;
                  p_[j + 1] = params2->c2_;
                  p_[j + 2] = params2->c3_;
                  break;
                }
                case INPAR::MAT::mes_isoexpopow:
                {
                  filename_ = filename_ + "_isoexpopow";
                  const MAT::ELASTIC::PAR::IsoExpoPow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoExpoPow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->k1_;
                  p_[j + 1] = params2->k2_;
                  break;
                }
                case INPAR::MAT::mes_isomooneyrivlin:
                {
                  filename_ = filename_ + "_isomooneyrivlin";
                  const MAT::ELASTIC::PAR::IsoMooneyRivlin* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoMooneyRivlin*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->c1_;
                  p_[j + 1] = params2->c2_;
                  break;
                }
                case INPAR::MAT::mes_volsussmanbathe:
                {
                  filename_ = filename_ + "_volsussmanbathe";
                  const MAT::ELASTIC::PAR::VolSussmanBathe* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::VolSussmanBathe*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->kappa_;
                  break;
                }
                case INPAR::MAT::mes_volpenalty:
                {
                  filename_ = filename_ + "_volpenalty";
                  const MAT::ELASTIC::PAR::VolPenalty* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::VolPenalty*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->eps_;
                  p_[j + 1] = params2->gam_;
                  break;
                }
                case INPAR::MAT::mes_vologden:
                {
                  filename_ = filename_ + "_vologden";
                  const MAT::ELASTIC::PAR::VolOgden* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->kappa_;
                  // p_[j+1] = params2->beta_;
                  break;
                }
                case INPAR::MAT::mes_volpow:
                {
                  filename_ = filename_ + "_volpow";
                  const MAT::ELASTIC::PAR::VolPow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::VolPow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->a_;
                  // p_[j+1]   = params2->expon_;
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
        case INPAR::MAT::m_viscoelasthyper:
        {
          materialhashistory_ = true;
          MAT::PAR::ViscoElastHyper* params =
              dynamic_cast<MAT::PAR::ViscoElastHyper*>(actmat->Parameter());
          if (!params) dserror("Cannot cast material parameters");
          const int nummat = params->nummat_;
          const std::vector<int>* matids = params->matids_;
          for (int i = 0; i < nummat; ++i)
          {
            const int id = (*matids)[i];

            if (mymatset.size() == 0 or mymatset.find(id) != mymatset.end())
            {
              const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
              switch (actelastmat->Type())
              {
                case INPAR::MAT::mes_couplogneohooke:
                {
                  filename_ = filename_ + "_couplogneohooke";
                  const MAT::ELASTIC::PAR::CoupLogNeoHooke* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupLogNeoHooke*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->mue_;
                  p_[j + 1] = params2->lambda_;
                  break;
                }
                case INPAR::MAT::mes_coupneohooke:
                {
                  filename_ = filename_ + "_coupneohooke";
                  const MAT::ELASTIC::PAR::CoupNeoHooke* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupNeoHooke*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->c_;
                  p_[j + 1] = params2->beta_;
                  break;
                }
                case INPAR::MAT::mes_coupblatzko:
                {
                  filename_ = filename_ + "_coupblatzko";
                  const MAT::ELASTIC::PAR::CoupBlatzKo* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupBlatzKo*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->mue_;
                  p_[j + 1] = params2->nue_ / (1. - 2. * params2->nue_);
                  // p_[j+2] = params2->f_;
                  break;
                }
                case INPAR::MAT::mes_coupsimopister:
                {
                  filename_ = filename_ + "_coupsimopister";
                  const MAT::ELASTIC::PAR::CoupSimoPister* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupSimoPister*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->mue_;
                  break;
                }
                case INPAR::MAT::mes_coupSVK:
                {
                  filename_ = filename_ + "_coupsaintvenantkirchhoff";
                  const MAT::ELASTIC::PAR::CoupSVK* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupSVK*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->lambda_;
                  p_[j + 1] = params2->mue_;
                  break;
                }
                case INPAR::MAT::mes_isoneohooke:
                {
                  filename_ = filename_ + "_isoneohooke";
                  const MAT::ELASTIC::PAR::IsoNeoHooke* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoNeoHooke*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->mue_;
                  break;
                }
                case INPAR::MAT::mes_isoyeoh:
                {
                  filename_ = filename_ + "_isoyeoh";
                  const MAT::ELASTIC::PAR::IsoYeoh* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 3);
                  p_[j] = params2->c1_;
                  p_[j + 1] = params2->c2_;
                  p_[j + 2] = params2->c3_;
                  break;
                }
                case INPAR::MAT::mes_iso1pow:
                {
                  filename_ = filename_ + "_iso1pow";
                  const MAT::ELASTIC::PAR::Iso1Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_iso2pow:
                {
                  filename_ = filename_ + "_iso2pow";
                  const MAT::ELASTIC::PAR::Iso2Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Iso2Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_coup1pow:
                {
                  filename_ = filename_ + "_coup1pow";
                  const MAT::ELASTIC::PAR::Coup1Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Coup1Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_coup2pow:
                {
                  filename_ = filename_ + "_coup2pow";
                  const MAT::ELASTIC::PAR::Coup2Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Coup2Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_coup3pow:
                {
                  filename_ = filename_ + "_coup3pow";
                  const MAT::ELASTIC::PAR::Coup3Pow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Coup3Pow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->c_;
                  break;
                }
                case INPAR::MAT::mes_coup13apow:
                {
                  filename_ = filename_ + "_coup13apow";
                  const MAT::ELASTIC::PAR::Coup13aPow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Coup13aPow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->c_;
                  p_[j + 1] = params2->a_;
                  break;
                }
                case INPAR::MAT::mes_coupmooneyrivlin:
                {
                  filename_ = filename_ + "_coupmooneyrivlin";
                  const MAT::ELASTIC::PAR::CoupMooneyRivlin* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupMooneyRivlin*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 3);
                  p_[j] = params2->c1_;
                  p_[j + 1] = params2->c2_;
                  p_[j + 2] = params2->c3_;
                  break;
                }
                case INPAR::MAT::mes_isoexpopow:
                {
                  filename_ = filename_ + "_isoexpopow";
                  const MAT::ELASTIC::PAR::IsoExpoPow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoExpoPow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->k1_;
                  p_[j + 1] = params2->k2_;
                  break;
                }
                case INPAR::MAT::mes_isomooneyrivlin:
                {
                  filename_ = filename_ + "_isomooneyrivlin";
                  const MAT::ELASTIC::PAR::IsoMooneyRivlin* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoMooneyRivlin*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->c1_;
                  p_[j + 1] = params2->c2_;
                  break;
                }
                case INPAR::MAT::mes_volsussmanbathe:
                {
                  filename_ = filename_ + "_volsussmanbathe";
                  const MAT::ELASTIC::PAR::VolSussmanBathe* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::VolSussmanBathe*>(
                          actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->kappa_;
                  break;
                }
                case INPAR::MAT::mes_volpenalty:
                {
                  filename_ = filename_ + "_volpenalty";
                  const MAT::ELASTIC::PAR::VolPenalty* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::VolPenalty*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 2);
                  p_[j] = params2->eps_;
                  p_[j + 1] = params2->gam_;
                  break;
                }
                case INPAR::MAT::mes_vologden:
                {
                  filename_ = filename_ + "_vologden";
                  const MAT::ELASTIC::PAR::VolOgden* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->kappa_;
                  // p_[j+1] = params2->beta_;
                  break;
                }
                case INPAR::MAT::mes_volpow:
                {
                  filename_ = filename_ + "_volpow";
                  const MAT::ELASTIC::PAR::VolPow* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::VolPow*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->a_;
                  // p_[j+1]   = params2->expon_;
                  break;
                }
                case INPAR::MAT::mes_coupexppol:
                {
                  filename_ = filename_ + "_coupexppol";
                  const MAT::ELASTIC::PAR::CoupExpPol* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupExpPol*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 3);
                  p_[j] = params2->a_;
                  p_[j + 1] = params2->b_;
                  p_[j + 2] = params2->c_;
                  // do inverse analysis with ln(b) and ln(c)
                  p_[j + 1] = log(p_[j + 1]);
                  p_[j + 2] = log(p_[j + 2]);
                  // remind position of b in p_
                  coupexppol_parameter_position_ = j + 1;
                  break;
                }
                case INPAR::MAT::mes_coupmyocard:
                {
                  filename_ = filename_ + "_coupmyocard";
                  const MAT::ELASTIC::PAR::CoupMyocard* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::CoupMyocard*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->n_;
                  break;
                }
                case INPAR::MAT::mes_isoratedep:
                {
                  filename_ = filename_ + "_isoratedep";
                  const MAT::ELASTIC::PAR::IsoRateDep* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::IsoRateDep*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->n_;
                  break;
                }
                case INPAR::MAT::mes_genmax:
                {
                  filename_ = filename_ + "_genmax";
                  const MAT::ELASTIC::PAR::GenMax* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::GenMax*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 1);
                  p_[j] = params2->tau_;
                  // p_[j+1]   = params2->beta_;
                  break;
                }
                case INPAR::MAT::mes_fract:
                {
                  filename_ = filename_ + "_fract";
                  const MAT::ELASTIC::PAR::Fract* params2 =
                      dynamic_cast<const MAT::ELASTIC::PAR::Fract*>(actelastmat->Parameter());
                  int j = p_.Length();
                  p_.Resize(j + 3);
                  p_[j] = params2->tau_;
                  p_[j + 1] = params2->alpha_;
                  p_[j + 2] = params2->beta_;
                  viscofract_alpha_position_ = j + 1;  // comment out this line if alpha is not fit!
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

        // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block or an interface
        // to a micro material
        case INPAR::MAT::mes_couplogneohooke:
        case INPAR::MAT::mes_coupexppol:
        case INPAR::MAT::mes_coupneohooke:
        case INPAR::MAT::mes_coupblatzko:
        case INPAR::MAT::mes_coupsimopister:
        case INPAR::MAT::mes_coupSVK:
        case INPAR::MAT::mes_isoneohooke:
        case INPAR::MAT::mes_isoyeoh:
        case INPAR::MAT::mes_iso1pow:
        case INPAR::MAT::mes_iso2pow:
        case INPAR::MAT::mes_coup1pow:
        case INPAR::MAT::mes_coup2pow:
        case INPAR::MAT::mes_coup3pow:
        case INPAR::MAT::mes_coup13apow:
        case INPAR::MAT::mes_isoexpopow:
        case INPAR::MAT::mes_isomooneyrivlin:
        case INPAR::MAT::mes_coupmooneyrivlin:
        case INPAR::MAT::mes_volsussmanbathe:
        case INPAR::MAT::mes_volpenalty:
        case INPAR::MAT::mes_volpow:
        case INPAR::MAT::mes_vologden:
        case INPAR::MAT::mes_isoratedep:
        case INPAR::MAT::mes_genmax:
        case INPAR::MAT::mes_fract:
        case INPAR::MAT::m_struct_multiscale:
          break;
        default:
          break;
      }
    }
  }

  // Read patches
  if (patches_)
  {
    if (matset_patches_.size() == 0) dserror("You have to specify a material in INV_LIST_PATCHES");

    // prepare check if only one material law is considered
    int mattype = -1;
    int numparams = 1;

    // read in starting values for patches
    std::vector<double> startvalues;
    if (which_patches_ == material) startvalues.resize(numpatches_);

    const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
        *DRT::Problem::Instance(0)->Materials()->Map();
    std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;

    for (curr = mats.begin(); curr != mats.end(); ++curr)
    {
      const Teuchos::RCP<MAT::PAR::Material> actmat = curr->second;

      if (matset_patches_.find(actmat->Id()) != matset_patches_.end())
      {
        if (mattype == -1)
          mattype = actmat->Type();
        else if (actmat->Type() != mattype)
          dserror("patchwise inverse analysis can only be performed with one material law");

        switch (actmat->Type())
        {
          case INPAR::MAT::m_constraintmixture:
          {
            numparams = 2;
            // in the case of patches defined by materials we want to use the values specified in
            // the material
            if (which_patches_ == material)
            {
              int patchid = 0;
              int count = 0;
              for (std::set<int>::const_iterator mat = matset_patches_.begin();
                   mat != matset_patches_.end(); mat++)
              {
                if (mat == matset_patches_.find(actmat->Id())) patchid = count;
                count += 1;
              }
              MAT::PAR::ConstraintMixture* params =
                  dynamic_cast<MAT::PAR::ConstraintMixture*>(actmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              startvalues[patchid] = sqrt(params->GetParameter(params->growthfactor, -1));
              startvalues.resize(numparams * numpatches_);
              // startvalues[patchid+numpatches_] =
              // acos(sqrt(params->GetParameter(params->elastin_survival,-1)));
              startvalues[patchid + numpatches_] =
                  tan(params->GetParameter(params->elastin_survival, -1) * PI - PI / 2.0);
            }
          }
          break;
          default:
            dserror("Material type not implemented for patches");
            break;
        }
      }
    }
    int numpatchparams = numparams * numpatches_;

    if (which_patches_ != material)
    {
      // input parameters inverse analysis
      const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
      double word1;
      std::istringstream matliststream(
          Teuchos::getNumericStringParameter(iap, "STARTVALUESFORPATCHES"));
      while (matliststream >> word1)
      {
        startvalues.push_back(sqrt(word1));
      }

      // check if size of startvalues fits to number of patches
      int size = startvalues.size();
      if (size != numpatchparams)
        dserror(
            "Number of start values %d does fit to number of patches %d and number of parameters "
            "%d!",
            startvalues.size(), numpatches_, numparams);
    }

    // copy parameters into the full parameter vector
    startindexpatches_ = p_.Length();
    p_.Resize(startindexpatches_ + numpatchparams);
    for (int idpatches = 0; idpatches < numpatchparams; idpatches++)
    {
      p_[startindexpatches_ + idpatches] = startvalues[idpatches];
    }
  }

  return;
}
//--------------------------------------------------------------------------------------
void STR::GenInvAnalysis::SetParameters(Epetra_SerialDenseVector p_cur)
{
  // check whether there is a micro scale
  Teuchos::RCP<Epetra_Comm> subcomm = DRT::Problem::Instance(0)->GetNPGroup()->SubComm();
  if (subcomm != Teuchos::null)
  {
    // tell supporting procs that material parameters have to be set
    int task[2] = {6, np_};
    subcomm->Broadcast(task, 2, 0);
    // broadcast p_cur to the micro scale
    subcomm->Broadcast(&p_cur[0], np_, 0);
  }

  // write new material parameter
  for (unsigned prob = 0; prob < DRT::Problem::NumInstances(); ++prob)
  {
    std::set<int> mymatset = matset_[prob];

    // in case of micro scale: broadcast mymatset once per problem instance
    if (subcomm != Teuchos::null)
    {
      LINALG::GatherAll<int>(mymatset, *subcomm);
    }

    // material parameters are set for the current problem instance
    STR::SetMaterialParameters(prob, p_cur, mymatset);
  }

  if (patches_)
  {
    SetMaterialParametersPatches(p_cur);
  }

  return;
}

void STR::SetMaterialParameters(int prob, Epetra_SerialDenseVector& p_cur, std::set<int>& mymatset)
{
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = DRT::Problem::Instance()->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  const int lmyrank = lcomm->MyPID();
  // const int gmyrank  = gcomm->MyPID();
  const int groupid = group->GroupId();
  // const int ngroup   = group->NumGroups();

  // loop all materials in problem
  const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
      *DRT::Problem::Instance(prob)->Materials()->Map();
  int j = 0;
  std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;
  for (curr = mats.begin(); curr != mats.end(); ++curr)
  {
    const Teuchos::RCP<MAT::PAR::Material> actmat = curr->second;

    switch (actmat->Type())
    {
      case INPAR::MAT::m_aaaneohooke:
      {
        if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
        {
          MAT::PAR::AAAneohooke* params = dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
          if (!params) dserror("Cannot cast material parameters");
          // Epetra_Map dummy_map(1,1,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm()));
          Epetra_Map dummy_map(1, 1, 0, *gcomm);

          Teuchos::RCP<Epetra_Vector> parameter = Teuchos::rcp(new Epetra_Vector(dummy_map, true));
          parameter->PutScalar(double(p_cur[j]));
          params->SetParameter(params->young, parameter);
          parameter->PutScalar(p_cur[j + 1]);
          params->SetParameter(params->beta, parameter);
          if (lmyrank == 0) printf("NPGroup %3d: ", groupid);
          if (lmyrank == 0)
            printf("MAT::PAR::AAAneohooke %20.15e %20.15e\n", p_cur[j], p_cur[j + 1]);
          j += 2;
        }
      }
      break;
      case INPAR::MAT::m_neohooke:
      {
        if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
        {
          MAT::PAR::NeoHooke* params = dynamic_cast<MAT::PAR::NeoHooke*>(actmat->Parameter());
          if (!params) dserror("Cannot cast material parameters");
          // This is a tiny little bit brutal!!!
          const_cast<double&>(params->youngs_) = p_cur[j];
          const_cast<double&>(params->poissonratio_) = p_cur[j + 1];
          if (lmyrank == 0) printf("NPGroup %3d: ", groupid);
          if (lmyrank == 0)
            printf("MAT::PAR::NeoHooke %20.15e %20.15e\n", params->youngs_, params->poissonratio_);
          j += 2;
        }
      }
      break;
      case INPAR::MAT::m_stvenant:
      {
        if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
        {
          MAT::PAR::StVenantKirchhoff* params =
              dynamic_cast<MAT::PAR::StVenantKirchhoff*>(actmat->Parameter());
          if (!params) dserror("Cannot cast material parameters");
          // This is a tiny little bit brutal!!!
          const_cast<double&>(params->youngs_) = p_cur[j];
          const_cast<double&>(params->poissonratio_) = p_cur[j + 1];
          if (lmyrank == 0)
            printf("MAT::PAR::StVenantKirchhoff %20.15e %20.15e\n", params->youngs_,
                params->poissonratio_);
          j += 2;
        }
      }
      break;
      case INPAR::MAT::m_elasthyper:
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i = 0; i < nummat; ++i)
        {
          const int id = (*matids)[i];

          if (mymatset.size() == 0 or mymatset.find(id) != mymatset.end())
          {
            const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
            switch (actelastmat->Type())
            {
              case INPAR::MAT::mes_couplogneohooke:
              {
                MAT::ELASTIC::PAR::CoupLogNeoHooke* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupLogNeoHooke*>(actelastmat->Parameter());
                params2->SetMue(abs(p_cur(j)));
                params2->SetLambda(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_coupexppol:
              {
                MAT::ELASTIC::PAR::CoupExpPol* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupExpPol*>(actelastmat->Parameter());
                params2->SetA(abs(p_cur(j)));
                params2->SetB(abs(p_cur(j + 1)));
                params2->SetC(abs(p_cur(j + 2)));
                j = j + 3;
                break;
              }
              case INPAR::MAT::mes_coupneohooke:
              {
                MAT::ELASTIC::PAR::CoupNeoHooke* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupNeoHooke*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)) / (4 * (1 + abs(p_cur(j + 1)))));
                params2->SetBeta(abs(p_cur(j + 1)) / (1 - 2 * abs(p_cur(j + 1))));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_coupblatzko:
              {
                MAT::ELASTIC::PAR::CoupBlatzKo* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupBlatzKo*>(actelastmat->Parameter());
                params2->SetMue(abs(p_cur(j)));
                params2->SetNue((abs(p_cur(j + 1))) / (2. * (abs(p_cur(j + 1)) + 1.)));
                // params2->SetF(abs(p_cur(j+2)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_coupsimopister:
              {
                MAT::ELASTIC::PAR::CoupSimoPister* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupSimoPister*>(actelastmat->Parameter());
                params2->SetMue(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coupSVK:
              {  // formulation with lambda and mue
                /*
                  MAT::ELASTIC::PAR::CoupSVK* params2 =
                      dynamic_cast<MAT::ELASTIC::PAR::CoupSVK*>(actelastmat->Parameter());
                  params2->SetLambda(abs(p_cur(j)));
                  params2->SetMue(abs(p_cur(j + 1)));
                  j = j + 2;
                  */
                // Formulation with youngs and nue
                MAT::ELASTIC::PAR::CoupSVK* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupSVK*>(actelastmat->Parameter());
                params2->SetLambda(
                    (p_cur(j) * p_cur(j + 1)) / ((1 + p_cur(j + 1)) * (1 - 2 * p_cur(j + 1))));
                params2->SetMue((p_cur(j)) / (2 * (1 + p_cur(j + 1))));
                j = j + 2;

                break;
              }
              case INPAR::MAT::mes_isoneohooke:
              {
                MAT::ELASTIC::PAR::IsoNeoHooke* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoNeoHooke*>(actelastmat->Parameter());
                params2->SetMue(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_isoyeoh:
              {
                MAT::ELASTIC::PAR::IsoYeoh* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
                params2->SetC1(abs(p_cur(j)));
                params2->SetC2(abs(p_cur(j + 1)));
                params2->SetC3(abs(p_cur(j + 2)));
                j = j + 3;
                break;
              }
              case INPAR::MAT::mes_iso1pow:
              {
                MAT::ELASTIC::PAR::Iso1Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_iso2pow:
              {
                MAT::ELASTIC::PAR::Iso2Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Iso2Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coup1pow:
              {
                MAT::ELASTIC::PAR::Coup1Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Coup1Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coup2pow:
              {
                MAT::ELASTIC::PAR::Coup2Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Coup2Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coup3pow:
              {
                MAT::ELASTIC::PAR::Coup3Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Coup3Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coup13apow:
              {
                MAT::ELASTIC::PAR::Coup13aPow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Coup13aPow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                params2->SetA(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_coupmooneyrivlin:
              {
                MAT::ELASTIC::PAR::CoupMooneyRivlin* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupMooneyRivlin*>(actelastmat->Parameter());
                params2->SetC1(abs(p_cur(j)));
                params2->SetC2(abs(p_cur(j + 1)));
                params2->SetC3(abs(p_cur(j + 2)));
                j = j + 3;
                break;
              }
              case INPAR::MAT::mes_isoexpopow:
              {
                MAT::ELASTIC::PAR::IsoExpoPow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoExpoPow*>(actelastmat->Parameter());
                params2->SetK1(abs(p_cur(j)));
                params2->SetK2(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_isomooneyrivlin:
              {
                MAT::ELASTIC::PAR::IsoMooneyRivlin* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoMooneyRivlin*>(actelastmat->Parameter());
                params2->SetC1(abs(p_cur(j)));
                params2->SetC2(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_volsussmanbathe:
              {
                MAT::ELASTIC::PAR::VolSussmanBathe* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::VolSussmanBathe*>(actelastmat->Parameter());
                params2->SetKappa(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_volpenalty:
              {
                MAT::ELASTIC::PAR::VolPenalty* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::VolPenalty*>(actelastmat->Parameter());
                params2->SetEpsilon(abs(p_cur(j)));
                params2->SetGamma(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_vologden:
              {
                MAT::ELASTIC::PAR::VolOgden* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
                params2->SetKappa(abs(p_cur(j)));
                // params2->SetBeta(abs(p_cur(j+1)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_volpow:
              {
                MAT::ELASTIC::PAR::VolPow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::VolPow*>(actelastmat->Parameter());
                params2->SetA(abs(p_cur(j)));
                // params2->SetExpon(abs(p_cur(j+1)));
                j = j + 1;
                break;
              }
              default:
                dserror("cannot deal with this material");
            }
          }
        }
      }
      break;
      case INPAR::MAT::mes_couplogneohooke:
      case INPAR::MAT::mes_coupexppol:
      case INPAR::MAT::mes_coupneohooke:
      case INPAR::MAT::mes_coupblatzko:
      case INPAR::MAT::mes_coupsimopister:
      case INPAR::MAT::mes_coupSVK:
      case INPAR::MAT::mes_isoneohooke:
      case INPAR::MAT::mes_isoyeoh:
      case INPAR::MAT::mes_iso1pow:
      case INPAR::MAT::mes_iso2pow:
      case INPAR::MAT::mes_coup1pow:
      case INPAR::MAT::mes_coup2pow:
      case INPAR::MAT::mes_coup3pow:
      case INPAR::MAT::mes_coup13apow:
      case INPAR::MAT::mes_isoexpopow:
      case INPAR::MAT::mes_isomooneyrivlin:
      case INPAR::MAT::mes_coupmooneyrivlin:
      case INPAR::MAT::mes_volsussmanbathe:
      case INPAR::MAT::mes_volpenalty:
      case INPAR::MAT::mes_vologden:
      case INPAR::MAT::mes_volpow:
      case INPAR::MAT::mes_coupmyocard:
      case INPAR::MAT::mes_isoratedep:
      case INPAR::MAT::mes_genmax:
      case INPAR::MAT::mes_fract:
      case INPAR::MAT::m_struct_multiscale:
        break;

      case INPAR::MAT::m_constraintmixture:
      {
        if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
        {
          MAT::PAR::ConstraintMixture* params =
              dynamic_cast<MAT::PAR::ConstraintMixture*>(actmat->Parameter());
          if (!params) dserror("Cannot cast material parameters");
          // This is a tiny little bit brutal!!!
          // const_cast<double&>(params->mue_)  = p_cur[j];
          // const_cast<double&>(params->k1_)  = p_cur[j+1];
          // const_cast<double&>(params->k2_)  = p_cur[j+2];
          // const_cast<double&>(params->prestretchcollagen_) = p_cur[j];

          Epetra_Map dummy_map(1, 1, 0, *gcomm);
          Teuchos::RCP<Epetra_Vector> parameter = Teuchos::rcp(new Epetra_Vector(dummy_map, true));
          parameter->PutScalar(double(p_cur[j]));
          params->SetParameter(params->growthfactor, parameter);
          if (lmyrank == 0) printf("NPGroup %3d: ", groupid);
          if (lmyrank == 0) printf("MAT::PAR::ConstraintMixture %20.15e\n", p_cur[j]);
          j += 1;
        }
      }
      break;
      case INPAR::MAT::m_growth_iso_stress:
      {
        if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
        {
          MAT::PAR::GrowthLawIsoStress* params =
              dynamic_cast<MAT::PAR::GrowthLawIsoStress*>(actmat->Parameter());
          if (!params) dserror("Cannot cast material parameters");
          // This is a tiny little bit brutal!!!
          const_cast<double&>(params->kthetaplus_) = p_cur[j];
          if (lmyrank == 0) printf("NPGroup %3d: ", groupid);
          if (lmyrank == 0) printf("MAT::PAR::Growth %20.15e\n", params->kthetaplus_);
          j += 1;
        }
      }
      break;
      case INPAR::MAT::m_viscoelasthyper:
      {
        MAT::PAR::ViscoElastHyper* params =
            dynamic_cast<MAT::PAR::ViscoElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i = 0; i < nummat; ++i)
        {
          const int id = (*matids)[i];

          if (mymatset.size() == 0 or mymatset.find(id) != mymatset.end())
          {
            const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
            switch (actelastmat->Type())
            {
              case INPAR::MAT::mes_couplogneohooke:
              {
                MAT::ELASTIC::PAR::CoupLogNeoHooke* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupLogNeoHooke*>(actelastmat->Parameter());
                params2->SetMue(abs(p_cur(j)));
                params2->SetLambda(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_coupexppol:
              {
                MAT::ELASTIC::PAR::CoupExpPol* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupExpPol*>(actelastmat->Parameter());
                params2->SetA(abs(p_cur(j)));
                params2->SetB(abs(p_cur(j + 1)));
                params2->SetC(abs(p_cur(j + 2)));
                j = j + 3;
                break;
              }
              case INPAR::MAT::mes_coupneohooke:
              {
                MAT::ELASTIC::PAR::CoupNeoHooke* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupNeoHooke*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                params2->SetBeta(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_coupblatzko:
              {
                MAT::ELASTIC::PAR::CoupBlatzKo* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupBlatzKo*>(actelastmat->Parameter());
                params2->SetMue(abs(p_cur(j)));
                params2->SetNue((abs(p_cur(j + 1))) / (2. * (abs(p_cur(j + 1)) + 1.)));
                //  params2->SetF(abs(p_cur(j+2)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_coupsimopister:
              {
                MAT::ELASTIC::PAR::CoupSimoPister* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupSimoPister*>(actelastmat->Parameter());
                params2->SetMue(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coupSVK:
              {
                MAT::ELASTIC::PAR::CoupSVK* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupSVK*>(actelastmat->Parameter());
                params2->SetLambda(abs(p_cur(j)));
                params2->SetMue(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_isoneohooke:
              {
                MAT::ELASTIC::PAR::IsoNeoHooke* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoNeoHooke*>(actelastmat->Parameter());
                params2->SetMue(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_isoyeoh:
              {
                MAT::ELASTIC::PAR::IsoYeoh* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
                params2->SetC1(abs(p_cur(j)));
                params2->SetC2(abs(p_cur(j + 1)));
                params2->SetC3(abs(p_cur(j + 2)));
                j = j + 3;
                break;
              }
              case INPAR::MAT::mes_iso1pow:
              {
                MAT::ELASTIC::PAR::Iso1Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_iso2pow:
              {
                MAT::ELASTIC::PAR::Iso2Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Iso2Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coup1pow:
              {
                MAT::ELASTIC::PAR::Coup1Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Coup1Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coup2pow:
              {
                MAT::ELASTIC::PAR::Coup2Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Coup2Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coup3pow:
              {
                MAT::ELASTIC::PAR::Coup3Pow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Coup3Pow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coup13apow:
              {
                MAT::ELASTIC::PAR::Coup13aPow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Coup13aPow*>(actelastmat->Parameter());
                params2->SetC(abs(p_cur(j)));
                params2->SetA(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_coupmooneyrivlin:
              {
                MAT::ELASTIC::PAR::CoupMooneyRivlin* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupMooneyRivlin*>(actelastmat->Parameter());
                params2->SetC1(abs(p_cur(j)));
                params2->SetC2(abs(p_cur(j + 1)));
                params2->SetC3(abs(p_cur(j + 2)));
                j = j + 3;
                break;
              }
              case INPAR::MAT::mes_isoexpopow:
              {
                MAT::ELASTIC::PAR::IsoExpoPow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoExpoPow*>(actelastmat->Parameter());
                params2->SetK1(abs(p_cur(j)));
                params2->SetK2(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_isomooneyrivlin:
              {
                MAT::ELASTIC::PAR::IsoMooneyRivlin* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoMooneyRivlin*>(actelastmat->Parameter());
                params2->SetC1(abs(p_cur(j)));
                params2->SetC2(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_volsussmanbathe:
              {
                MAT::ELASTIC::PAR::VolSussmanBathe* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::VolSussmanBathe*>(actelastmat->Parameter());
                params2->SetKappa(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_volpenalty:
              {
                MAT::ELASTIC::PAR::VolPenalty* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::VolPenalty*>(actelastmat->Parameter());
                params2->SetEpsilon(abs(p_cur(j)));
                params2->SetGamma(abs(p_cur(j + 1)));
                j = j + 2;
                break;
              }
              case INPAR::MAT::mes_vologden:
              {
                MAT::ELASTIC::PAR::VolOgden* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
                params2->SetKappa(abs(p_cur(j)));
                // params2->SetBeta(abs(p_cur(j+1)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_volpow:
              {
                MAT::ELASTIC::PAR::VolPow* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::VolPow*>(actelastmat->Parameter());
                params2->SetA(abs(p_cur(j)));
                // params2->SetExpon(abs(p_cur(j+1)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_coupmyocard:
              {
                MAT::ELASTIC::PAR::CoupMyocard* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::CoupMyocard*>(actelastmat->Parameter());
                params2->SetN(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_isoratedep:
              {
                MAT::ELASTIC::PAR::IsoRateDep* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::IsoRateDep*>(actelastmat->Parameter());
                params2->SetN(abs(p_cur(j)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_genmax:
              {
                MAT::ELASTIC::PAR::GenMax* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::GenMax*>(actelastmat->Parameter());
                params2->SetTau(abs(p_cur(j)));
                // params2->SetBeta(abs(p_cur(j+1)));
                j = j + 1;
                break;
              }
              case INPAR::MAT::mes_fract:
              {
                MAT::ELASTIC::PAR::Fract* params2 =
                    dynamic_cast<MAT::ELASTIC::PAR::Fract*>(actelastmat->Parameter());
                params2->SetTau(abs(p_cur(j)));
                params2->SetAlpha(abs(p_cur(j + 1)));
                params2->SetBeta(abs(p_cur(j + 2)));
                j = j + 3;
                break;
              }
              default:
                dserror("cannot deal with this material");
            }
          }
        }
      }
      break;
      default:
      {
        // ignore unknown materials ?
        if (mymatset.size() == 0 or mymatset.find(actmat->Id()) != mymatset.end())
          dserror("Unknown type of material");
      }
      break;
    }
  }
  return;
}

void STR::GenInvAnalysis::SetMaterialParametersPatches(Epetra_SerialDenseVector& p_cur)
{
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = DRT::Problem::Instance()->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> lcomm = group->LocalComm();
  const int lmyrank = lcomm->MyPID();
  const int groupid = group->GroupId();

  int numparams = (p_cur.Length() - startindexpatches_) / numpatches_;

  for (int id_param = 0; id_param < numparams; id_param++)
  {
    Epetra_Vector p_patches(smoother_->DomainMap(), true);
    for (int idpatches = 0; idpatches < numpatches_; idpatches++)
    {
      if (lmyrank == 0) printf("NPGroup %3d: ", groupid);
      if (lmyrank == 0)
        printf("MAT::PAR::ConstraintMixture %20.15e\n",
            p_cur[startindexpatches_ + idpatches + numpatches_ * id_param]);
      if (id_param == 0)
        p_patches[idpatches] = p_cur[startindexpatches_ + idpatches + numpatches_ * id_param] *
                               p_cur[startindexpatches_ + idpatches + numpatches_ * id_param];
      else if (id_param == 1)
      {
        // p_patches[idpatches] =
        // cos(p_cur[startindexpatches_+idpatches+numpatches_*id_param])*cos(p_cur[startindexpatches_+idpatches+numpatches_*id_param]);
        p_patches[idpatches] =
            (atan(p_cur[startindexpatches_ + idpatches + numpatches_ * id_param]) + PI / 2.0) / PI;
      }
    }

    Teuchos::RCP<Epetra_Vector> parameter =
        Teuchos::rcp(new Epetra_Vector(*discret_->ElementColMap(), true));
    ComputeParametersFromPatches(p_patches, parameter);

    // loop all materials in problem
    const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
        *DRT::Problem::Instance(0)->Materials()->Map();
    std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;
    for (curr = mats.begin(); curr != mats.end(); ++curr)
    {
      const Teuchos::RCP<MAT::PAR::Material> actmat = curr->second;

      switch (actmat->Type())
      {
        case INPAR::MAT::m_constraintmixture:
        {
          if (matset_patches_.find(actmat->Id()) != matset_patches_.end())
          {
            MAT::PAR::ConstraintMixture* params =
                dynamic_cast<MAT::PAR::ConstraintMixture*>(actmat->Parameter());
            if (!params) dserror("Cannot cast material parameters");
            // parameter->Multiply(1.0,*parameter,*parameter,0.0);
            if (id_param == 0)
              params->SetParameter(params->growthfactor, parameter);
            else if (id_param == 1)
              params->SetParameter(params->elastin_survival, parameter);
            else
              dserror("MAT::ConstraintMixture has only 2 parameters not %d", id_param);
          }
        }
        break;
        default:
        {
          // ignore unknown materials ?
          if (matset_patches_.find(actmat->Id()) != matset_patches_.end())
            dserror("Unknown type of material");
        }
        break;
      }
    }
  }
  return;
}

void STR::GenInvAnalysis::InitPatches()
{
  // some information
  int myrank = discret_->Comm().MyPID();
  int length = (*discret_->ElementRowMap()).NumMyElements();

  // -----------------------------------------------------------------------------------
  // build graph of neighbouring elements
  // for DEFINEPATCHES MaterialNumber this graph contains only elements that are part of the patches
  // in Uniform all elements belong to a patch
  int maxband = 18;  // just an approximation
  // graph of neighbouring elements
  Teuchos::RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, *discret_->ElementRowMap(), maxband, false));
  // contains a list of gid and ghostelegid or ghostelegid1 and ghostelegid2
  std::vector<int> ghostelements;
  for (int i = 0; i < length; ++i)
  {
    std::vector<int> neighbouringelements;
    int gid = (*discret_->ElementRowMap()).GID(i);
    DRT::Element* myele = discret_->gElement(gid);
    // only include elements in the graph that are part of patches
    int matid = myele->Material()->Parameter()->Id();
    if (which_patches_ != material || matset_patches_.find(matid) != matset_patches_.end())
    {
      DRT::Node** mynodes = myele->Nodes();
      for (int idnodes = 0; idnodes < myele->NumNode(); idnodes++)
      {
        DRT::Node* locnode = mynodes[idnodes];
        DRT::Element** node_ele = locnode->Elements();
        int otherowner = -1;
        for (int id_ele = 0; id_ele < locnode->NumElement(); id_ele++)
        {
          int locmatid = node_ele[id_ele]->Material()->Parameter()->Id();
          if (which_patches_ != material || matset_patches_.find(locmatid) != matset_patches_.end())
          {
            neighbouringelements.push_back(node_ele[id_ele]->Id());
            if (node_ele[id_ele]->Owner() != myrank)
            {
              // store the two coordinates and send them to the different procs later
              ghostelements.push_back(gid);
              ghostelements.push_back(node_ele[id_ele]->Id());
              // catch case when one node has ghostelements from two different procs
              if (otherowner == -1)
              {
                otherowner = node_ele[id_ele]->Owner();
              }
              else if (otherowner != node_ele[id_ele]->Owner())
              {
                for (int id_ele2 = 0; id_ele2 < locnode->NumElement(); id_ele2++)
                {
                  int locmatid2 = node_ele[id_ele2]->Material()->Parameter()->Id();
                  if (which_patches_ != material ||
                      matset_patches_.find(locmatid2) != matset_patches_.end())
                  {
                    if (node_ele[id_ele2]->Owner() != node_ele[id_ele]->Owner())
                    {
                      ghostelements.push_back(node_ele[id_ele2]->Id());
                      ghostelements.push_back(node_ele[id_ele]->Id());
                      ghostelements.push_back(node_ele[id_ele]->Id());
                      ghostelements.push_back(node_ele[id_ele2]->Id());
                    }
                  }
                }
              }
            }  // ghostele
          }
        }  // loop over elements of one node
      }    // loop over nodes
    }

    int num = neighbouringelements.size();
    int err = graph->InsertGlobalIndices(gid, num, &neighbouringelements[0]);
    if (err < 0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d", err);
  }

  // some elements are not ghosted, get them from the other procs and add them to the colmap
  std::vector<int> addcolumns;

  // RoundRobin
  int numproc = discret_->Comm().NumProc();
  const int torank = (myrank + 1) % numproc;              // to
  const int fromrank = (myrank + numproc - 1) % numproc;  // from
  DRT::Exporter exporter(discret_->Comm());

  for (int irobin = 0; irobin < numproc; ++irobin)
  {
    std::vector<char> sdata;
    std::vector<char> rdata;
    // ---- pack data for sending -----
    {
      DRT::PackBuffer data;
      {
        DRT::PackBuffer::SizeMarker sm(data);
        sm.Insert();
        DRT::ParObject::AddtoPack(data, ghostelements);
      }
      data.StartPacking();
      {
        // DRT::PackBuffer::SizeMarker sm2( data );
        // sm2.Insert();
        DRT::ParObject::AddtoPack(data, ghostelements);
      }
      std::swap(sdata, data());
    }

    // ---- send ----
    MPI_Request request;
    exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);


    // ---- receive ----
    int length = rdata.size();
    int tag = -1;
    int from = -1;
    exporter.ReceiveAny(from, tag, rdata, length);
    if (tag != 1234 or from != fromrank)
      dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank, myrank,
          from, myrank);

    // ---- unpack ----
    {
      // Put received ghostelements either into graph or into list of ghostelements
      ghostelements.clear();
      std::vector<int> receivedelements;
      std::vector<char>::size_type index = 0;
      if (rdata.size() > 0)
      {
        DRT::ParObject::ExtractfromPack(index, rdata, receivedelements);
        for (unsigned int id_eles = 0; id_eles < receivedelements.size(); id_eles += 2)
        {
          int gid = receivedelements[id_eles];
          int gid2 = receivedelements[id_eles + 1];
          // check if gid2 is one of my elements, otherwise send it to the next proc
          if (graph->MyGlobalRow(gid2))
          {
            int err = graph->InsertGlobalIndices(gid2, 1, &gid);
            if (err < 0) dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d", err);
            addcolumns.push_back(gid);
          }
          else
          {
            ghostelements.push_back(gid);
            ghostelements.push_back(gid2);
          }
        }
      }
    }

    // wait for all communication to finish
    exporter.Wait(request);
    discret_->Comm().Barrier();  // I feel better this way ;-)
  }                              // end for irobin
  // check if list is empty
  if (ghostelements.size() > 0) dserror("Some ghost elements did not find their owner!");

  // give the graph the modified map
  int errgraph = graph->FillComplete();
  if (errgraph != 0) dserror("FillComplete of the graph did not work");

  // -----------------------------------------------------------------------------------
  // build matrix that does the smoothing
  // (1 - omega* D^-1 A )
  double omega = 2.0 / 3.0;
  Teuchos::RCP<Epetra_CrsMatrix> laplace =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *discret_->ElementRowMap(), maxband));
  for (int i = 0; i < length; ++i)
  {
    int gid = (*discret_->ElementRowMap()).GID(i);
    int lengthrow = graph->NumGlobalIndices(gid);
    std::vector<int> myrow(lengthrow, 0);
    int graphlength;
    graph->ExtractGlobalRowCopy(gid, lengthrow, graphlength, &myrow[0]);
    if (graphlength != lengthrow) dserror("something went wrong with row length");
    double offdiagonal = omega;
    if (lengthrow > 1) offdiagonal = offdiagonal / (lengthrow - 1.0);
    for (int idrow = 0; idrow < lengthrow; idrow++)
    {
      if (myrow[idrow] != gid)
      {
        int err = laplace->InsertGlobalValues(gid, 1, &offdiagonal, &myrow[idrow]);
        if (err < 0) dserror("Epetra_CrsMatrix::InsertGlobalIndices returned %d", err);
      }
    }
    double diagonal = 1.0 - omega;
    if (lengthrow <= 1) diagonal = 1.0;
    int err = laplace->InsertGlobalValues(gid, 1, &diagonal, &gid);
    if (err < 0) dserror("Epetra_CrsMatrix::InsertGlobalIndices returned %d", err);
  }
  int errlaplace = laplace->FillComplete();
  if (errlaplace != 0) dserror("FillComplete of the laplace matrix did not work");

  // -----------------------------------------------------------------------------------
  // build matrix with patch information
  Teuchos::RCP<Epetra_CrsMatrix> patches =
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *discret_->ElementRowMap(), 1));
  if (which_patches_ == uniform)
  {
    int maxele = discret_->NumGlobalElements();
    int numeleperpatch =
        maxele / numpatches_;  // if this is not an integer the last patch will be the biggest
    double one = 1.0;
    for (int i = 0; i < length; ++i)
    {
      int gid = (*discret_->ElementRowMap()).GID(i);
      int patchnumber = 0;
      for (int idpatches = 1; idpatches < numpatches_; idpatches++)
      {
        if (gid > (idpatches * numeleperpatch - 1))
        {
          patchnumber += 1;
        }
      }
      int err = patches->InsertGlobalValues(gid, 1, &one, &patchnumber);
      if (err != 0)
        dserror("InsertGlobalValues for entry (%d,%d) failed with err=%d", gid, patchnumber, err);
    }
  }
  else if (which_patches_ == material)
  {
    double one = 1.0;
    for (int i = 0; i < length; ++i)
    {
      int gid = (*discret_->ElementRowMap()).GID(i);
      int matid = discret_->gElement(gid)->Material()->Parameter()->Id();
      int patchid = -1;
      int count = 0;
      for (std::set<int>::const_iterator mat = matset_patches_.begin();
           mat != matset_patches_.end(); mat++)
      {
        if (mat == matset_patches_.find(matid)) patchid = count;
        count += 1;
      }
      if (patchid != -1)
      {
        int err = patches->InsertGlobalValues(gid, 1, &one, &patchid);
        if (err != 0)
          dserror("InsertGlobalValues for entry (%d,%d) failed with err=%d", gid, matid, err);
      }
    }
  }
  Epetra_Map dummy_map(numpatches_, numpatches_, 0, discret_->Comm());
  patches->FillComplete(dummy_map, *discret_->ElementRowMap());

  // -----------------------------------------------------------------------------------
  // multiply matrices to get final smoothing matrix
  if (smoothingsteps_ == 0)
  {
    smoother_ = Teuchos::rcp(new LINALG::SparseMatrix(patches, LINALG::Copy));
  }
  else
  {
    smoother_ = Teuchos::rcp(new LINALG::SparseMatrix(*discret_->ElementRowMap(), numpatches_));
    smoother_ = LINALG::MLMultiply(*laplace, *patches, false, false, true);
    LINALG::SparseMatrix linalg_laplace(laplace, LINALG::Copy);
    for (int idsmooth = 1; idsmooth < smoothingsteps_; idsmooth++)
      smoother_ = LINALG::MLMultiply(linalg_laplace, *smoother_, false, false, true);
  }
}

void STR::GenInvAnalysis::ComputeParametersFromPatches(
    Epetra_Vector p_patches, Teuchos::RCP<Epetra_Vector>& eleparameters)
{
  Epetra_Vector rowparameters(*discret_->ElementRowMap(), true);
  smoother_->Multiply(false, p_patches, rowparameters);
  LINALG::Export(rowparameters, *eleparameters);
}

void STR::GenInvAnalysis::CheckOptStep()
{
  for (int i = 0; i < np_; i++)
  {
    if (p_(i) < 0.0)
    {
      std::cout << "at least one negative material parameter detected: reverse update" << std::endl;
      // reset the update step in case it was too large and increase regularization;
      // storage is not touched, so screen print out does not reflect what really happens
      p_ = p_o_;
      error_ = error_o_;
      error_grad_ = error_grad_o_;
      mu_ = mu_o_ * 2.0;
      return;
    }
  }
}


//===========================================================================================
void STR::GenInvAnalysis::MultiInvAnaInit()
{
  for (int i = 0; i < discret_->NumMyColElements(); i++)
  {
    DRT::Element* actele = discret_->lColElement(i);
    Teuchos::RCP<MAT::Material> mat = actele->Material();
    if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      MAT::MicroMaterial* micro = static_cast<MAT::MicroMaterial*>(mat.get());
      bool eleowner = false;
      if (discret_->Comm().MyPID() == actele->Owner()) eleowner = true;
      micro->InvAnaInit(eleowner, actele->Id());
    }
  }
}


//===========================================================================================
/*----------------------------------------------------------------------*/
/* data of volume-pressure-change experiment
 * (just applicable for lung parenchyma)                 abirzle 12/17  */
void STR::GenInvAnalysis::pVData(
    std::vector<double>& volexp_deltap, std::vector<double>& volexp_deltaV)

{
  // for description of experiment and measurements see
  // Birzle A.M., Martin C., Yoshihara L., Uhlig S., Wall W.A. (2018):
  // Experimental characterization and model identification of the nonlinear compressible material
  // behavior of lung parenchyma, Journal of the Mechanical Behavior of Biomedical Materials, 77,
  // 754-763

  // pressure values
  int iend = 0;
  // p=2.5
  for (int i = 0; i < 10; i++) volexp_deltap[i] = {245.};
  iend += 10;
  // p=5
  for (int i = iend; i < iend + 20; i++) volexp_deltap[i] = {490.};
  iend += 20;
  // p=7.5
  for (int i = iend; i < iend + 36; i++) volexp_deltap[i] = {735.};
  iend += 36;
  // p=10
  for (int i = iend; i < iend + 35; i++) volexp_deltap[i] = {980.};
  iend += 35;
  // p=12.5
  for (int i = iend; i < iend + 36; i++) volexp_deltap[i] = {1225.};
  iend += 36;
  // p=15
  for (int i = iend; i < iend + 37; i++) volexp_deltap[i] = {1470.};
  iend += 37;
  // p=17.5
  for (int i = iend; i < iend + 36; i++) volexp_deltap[i] = {1715.};
  iend += 36;
  // p=20
  for (int i = iend; i < iend + 29; i++) volexp_deltap[i] = {1960.};
  iend += 29;
  // p=25
  for (int i = iend; i < iend + 25; i++) volexp_deltap[i] = {2450.};
  iend += 25;
  // p=30
  for (int i = iend; i < iend + 23; i++) volexp_deltap[i] = {2940.};

  // volume-change values
  volexp_deltaV = {// p=2.5
      1.6229, 1.9215, 1.622, 1.8658, 2.0495, 1.5765, 1.4661, 1.4233, 1.7377, 1.7867,
      // p=5
      3.7862, 3.8465, 3.3126, 3.3772, 3.2526, 2.6585, 3.4688, 3.3288, 3.2557, 3.1412, 2.6198,
      2.7568, 4.2779, 3.6789, 3.6505, 3.5996, 4.263, 3.2006, 3.8828, 3.64140,
      // p=7.5
      5.5336, 5.1114, 5.2201, 5.8763, 5.5818, 5.0764, 5.6554, 5.2189, 5.4851, 5.8187, 4.8171,
      5.2599, 4.983, 4.8974, 5.927, 5.416, 4.7944, 5.1813, 5.8164, 5.5941, 5.4742, 4.7867, 5.8561,
      5.321, 5.3231, 5.9927, 5.6989, 6.0252, 5.8827, 6.098, 5.2408, 4.9775, 5.3245, 5.799, 5.4763,
      5.5281,
      // p=10
      5.1607, 5.5248, 6.3745, 5.2035, 5.949, 6.5293, 5.7254, 5.4395, 5.5039, 6.2416, 5.4003, 6.0278,
      5.3075, 5.0073, 5.5973, 5.4907, 5.599, 6.1918, 5.9618, 5.5298, 6.1709, 6.3988, 6.1538, 6.0434,
      6.0947, 5.7312, 5.6098, 5.0424, 5.963, 5.7041, 6.0714, 5.7577, 5.4448, 5.3786, 6.0825,
      // p=12.5
      6.334, 6.5669, 5.5371, 6.5367, 5.9702, 6.0108, 6.1958, 5.5998, 6.143, 6.3129, 6.0195, 6.4948,
      6.3385, 6.8904, 6.2244, 6.2614, 5.7878, 6.4638, 6.4166, 5.9079, 6.4859, 6.2747, 6.9097,
      6.4257, 5.9645, 6.293, 6.3622, 6.95119426478467, 6.95259299289317, 6.29096065891005,
      6.43146272902885, 6.76639025833022, 5.81392007558196, 6.64797014073652, 6.78575903857122,
      6.77079145333657,
      // p=15
      6.2694, 6.8156, 6.4222, 6.6779, 6.4477, 6.946, 6.0315, 7.1796, 5.9834, 6.8852, 6.816, 6.562,
      6.9666, 6.8255, 6.6346, 7.0178, 6.7803, 6.9803, 6.0881, 6.5954, 6.8001, 6.0148, 6.2209,
      5.9416, 5.9046, 6.1117, 6.8738, 6.2452, 7.1326, 6.1593, 6.4043, 6.7863, 6.6326, 7.0466,
      6.3976, 6.3152, 7.1319,
      // p=17.5
      6.8341, 7.0194, 7.1327, 6.8067, 6.6793, 6.4689, 6.5056, 6.9019, 7.0894, 7.1208, 6.636, 7.0989,
      6.6911, 6.813, 6.943, 6.7678, 6.2031, 6.3, 6.7398, 6.4996, 6.0093, 6.726, 6.4589,
      6.46228453249002, 6.72269688082884, 7.26996172483992, 6.59701895651191, 6.48882023137021,
      6.64521277786977, 6.53037380729109, 6.84867336679295, 7.27482643080461, 7.28380705707665,
      6.59131173065240, 7.32076484497370, 6.50402799776321,
      // p=20
      6.6185, 7.2038, 7.0779, 7.3916, 6.8441, 6.9745, 6.4551, 6.6086, 6.2844, 6.7463, 6.5037,
      7.2111, 6.8724, 6.1147, 6.7871, 7.1206, 6.8298, 6.7577, 6.9406, 6.5804, 6.6009, 6.772, 7.1248,
      6.9741, 6.7305, 6.6312, 6.9349, 6.6075, 6.53,
      // p=25
      7.1142, 7.1994, 7.1104, 6.9684, 7.0248, 6.9105, 7.2144, 7.3821, 6.5336, 6.9874, 7.402, 6.9957,
      7.0756, 6.8509, 7.0822, 7.1881, 7.2234, 6.9572, 7.0065, 7.1812, 7.2786, 7.0574, 6.6902, 7.288,
      6.7679,
      // p=30
      7.4768, 7.63, 7.5284, 7.3955, 8.0542, 7.5079, 7.8074, 7.1117, 7.4324, 7.5605, 7.8052, 7.2476,
      7.0182, 7.6021, 7.3106, 8.1168, 7.1584, 7.4902, 7.8548, 7.4139, 7.129, 7.5809, 7.8425};

  return;
}


//===========================================================================================
/*----------------------------------------------------------------------*/
/*  calculation of pressure-value from volume-change and material paramters
 * the material paramters depend on the material, which is identified
 * this function is implemented for  Coupneohooke + Iso1Pow + Coup3Pow
 * (just applicable for lung parenchyma)                 abirzle 12/17  */
void STR::GenInvAnalysis::ComputePressureNHIso1Coup3(std::vector<double> perturb,
    std::vector<double> volexp_deltaV, Epetra_SerialDenseMatrix& volexp_p_comp)
{
  // this implementation works for the material combination:
  // Coupneohooke + Iso1Pow + Coup3Pow

  // Saftey check and get not fitted material parameters from dat-file
  // (correct material combination is not checked)
  int d1_const = 0;
  int d3_const = 0;
  double beta_const = 0.;
  for (unsigned prob = 0; prob < DRT::Problem::NumInstances(); ++prob)
  {
    const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
        *DRT::Problem::Instance(prob)->Materials()->Map();
    std::set<int> mymatset = matset_[prob];
    std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;
    for (curr = mats.begin(); curr != mats.end(); ++curr)
    {
      const Teuchos::RCP<MAT::PAR::Material> actmat = curr->second;
      if (actmat->Type() == INPAR::MAT::m_elasthyper)
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i = 0; i < nummat; ++i)
        {
          const int id = (*matids)[i];
          if (mymatset.size() == 0 or mymatset.find(id) != mymatset.end())
          {
            const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
            if (actelastmat->Type() == INPAR::MAT::mes_iso1pow)
            {
              const MAT::ELASTIC::PAR::Iso1Pow* params2 =
                  dynamic_cast<const MAT::ELASTIC::PAR::Iso1Pow*>(actelastmat->Parameter());
              d1_const = params2->d_;
            }
            else if (actelastmat->Type() == INPAR::MAT::mes_coupneohooke)
            {
              const MAT::ELASTIC::PAR::CoupNeoHooke* params2 =
                  dynamic_cast<const MAT::ELASTIC::PAR::CoupNeoHooke*>(actelastmat->Parameter());
              beta_const = params2->beta_;
            }
            else if (actelastmat->Type() == INPAR::MAT::mes_coup3pow)
            {
              const MAT::ELASTIC::PAR::Coup3Pow* params2 =
                  dynamic_cast<const MAT::ELASTIC::PAR::Coup3Pow*>(actelastmat->Parameter());
              d3_const = params2->d_;
            }
            else
              dserror("This material is not implemented for coupled inverse analysis, yet.");
          }
        }
      }
    }
  }


  // Get Parameters
  // of this material combination
  Epetra_SerialDenseMatrix volexp_params(
      np_ + 1, 6);  // 6 = number of parameters in this material combination

  // parameter order: cNH, beta, c1, d1, c3, d3
  // case cNH, c1 and c3 are fitted
  if (np_ == 3)
  {
    for (int j = 0; j < np_ + 1; j++)
    {
      volexp_params(j, 0) = p_(0);  // cNH
      volexp_params(j, 4) = p_(2);  // c3
    }
    volexp_params(0, 0) += perturb[0];
    volexp_params(2, 4) += perturb[2];
  }
  // case cNH, beta, c1 and c3 are fitted
  else if (np_ == 4)
  {
    for (int j = 0; j < np_ + 1; j++)
    {
      volexp_params(j, 0) = p_(0);  // cNH
      volexp_params(j, 1) = p_(1);  // beta
      volexp_params(j, 4) = p_(3);  // c3
    }
    volexp_params(0, 0) += perturb[0];
    volexp_params(1, 1) += perturb[1];
    volexp_params(3, 4) += perturb[3];
  }
  else
    dserror("This case is not implemented, yet.");

  // set fixed parameters
  if (np_ == 3)  // case cNH, c1 and c3 are fitted
  {
    for (int j = 0; j < np_ + 1; j++)
    {
      volexp_params(j, 1) = beta_const;
      volexp_params(j, 3) = d1_const;
      volexp_params(j, 5) = d3_const;
      volexp_params(j, 2) = p_(2);  // c1 -> constant for pV
    }
  }
  else if (np_ == 4)  // case cNH, beta, c1 and c3 are fitted
  {
    for (int j = 0; j < np_ + 1; j++)
    {
      volexp_params(j, 3) = d1_const;
      volexp_params(j, 5) = d3_const;
      volexp_params(j, 2) = p_(2);  // c1 -> constant for pV
    }
  }
  else
    dserror("Something went wrong in the material parameter set.");

  // compute pressure for coupneohooke + iso1pow + coup3pow
  // Epetra_SerialDenseMatrix volexp_p_comp(nmp_volexp,np_+1);
  for (int i = 0; i < nmp_volexp_; i++)  // data set
    for (int j = 0; j < np_ + 1; j++)    // parameter set
    {
      // NeoHooke
      volexp_p_comp(i, j) =
          2. * volexp_params(j, 0) * std::pow(volexp_deltaV[i], -1. / 3.) -
          2. * volexp_params(j, 0) * std::pow(volexp_deltaV[i], -2. * volexp_params(j, 1) - 1.);

      // Iso1Pow --> no contribution of isomaterials

      // Coup3Pow
      volexp_p_comp(i, j) +=
          2. / 3. * volexp_params(j, 4) * volexp_params(j, 5) *
          (std::pow(volexp_deltaV[i], -1. / 3.)) *
          (std::pow((std::pow(volexp_deltaV[i], (2. / 3.)) - 1.), (volexp_params(j, 5) - 1.)));
    }
}


//===========================================================================================
/*----------------------------------------------------------------------*/
/* check if pressure-volume-change relation is physiological
 * for identified material parameters
 * otherwise do not accept identification step           abirzle 12/17  */
void STR::GenInvAnalysis::CheckPhysiologicalPVRelationNHIso1Coup3()
{
  // checked volume change range J in (1,10)
  const int Jrange = 901;
  double JpV[Jrange] = {};
  int j = 0;
  for (double i = 1.00; i <= 10.00; i = i + 0.01)
  {
    JpV[j] = i;
    j++;
  }

  // get constant values from dat-files
  int d3_const = 0;
  for (unsigned prob = 0; prob < DRT::Problem::NumInstances(); ++prob)
  {
    const std::map<int, Teuchos::RCP<MAT::PAR::Material>>& mats =
        *DRT::Problem::Instance(prob)->Materials()->Map();
    std::set<int> mymatset = matset_[prob];
    std::map<int, Teuchos::RCP<MAT::PAR::Material>>::const_iterator curr;
    for (curr = mats.begin(); curr != mats.end(); ++curr)
    {
      const Teuchos::RCP<MAT::PAR::Material> actmat = curr->second;
      if (actmat->Type() == INPAR::MAT::m_elasthyper)
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i = 0; i < nummat; ++i)
        {
          const int id = (*matids)[i];
          if (mymatset.size() == 0 or mymatset.find(id) != mymatset.end())
          {
            const Teuchos::RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
            if (actelastmat->Type() == INPAR::MAT::mes_coup3pow)
            {
              const MAT::ELASTIC::PAR::Coup3Pow* params2 =
                  dynamic_cast<const MAT::ELASTIC::PAR::Coup3Pow*>(actelastmat->Parameter());
              d3_const = params2->d_;
            }
          }
        }
      }
    }
  }

  // calculate dp/dJ = d^2Psi/dJ^2
  Epetra_SerialDenseVector dpdJ(Jrange);

  // case cNH, beta, c1 and c3 are fitted
  if (np_ == 4)
  {
    for (int i = 0; i < Jrange; i++)
    {
      // NeoHooke
      dpdJ[i] += -2. / 3. * p_(0) * pow(JpV[i], -4. / 3.) +
                 2. * p_(0) * (2. * p_(1) + 1.) * pow(JpV[i], -2. * p_(1) - 2.);

      // Iso1Pow --> no contribution of isomaterials

      // Coup3Pow
      dpdJ[i] +=
          2. / 3. * p_(3) * d3_const *
          (-1. / 3. * pow(JpV[i], -4. / 3.) * pow((pow(JpV[i], 2. / 3.) - 1.), d3_const - 1.) +
              2. / 3. * (d3_const - 1.) * pow(JpV[i], -2. / 3.) *
                  pow((pow(JpV[i], 2. / 3.) - 1.), d3_const - 2.));
    }
  }
  else
    dserror("This case is not implemented, yet.");


  for (int i = 0; i < Jrange; i++)
  {
    if (dpdJ[i] < 0)
    {
      std::cout << "at least one negative dpdJ detected: reverse update" << std::endl;
      // reset the update step in case it was too large and increase regularization;
      // storage is not touched, so screen print out does not reflect what really happens
      p_ = p_o_;
      error_ = error_o_;
      error_grad_ = error_grad_o_;
      mu_ = mu_o_ * 2.0;
      return;
    }
  }
}


//===========================================================================================
/*----------------------------------------------------------------------*/
/* check if pressure-volume-change relation is physiological
 * of multiple inverse analysis                          abirzle 12/17  */
void STR::GenInvAnalysis::CheckDiffDat(int& expid)
{
  Teuchos::RCP<COMM_UTILS::NestedParGroup> group = DRT::Problem::Instance()->GetNPGroup();
  Teuchos::RCP<Epetra_Comm> gcomm = group->GlobalComm();
  const int groupid = group->GroupId();
  const int ngroup = group->NumGroups();
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();

  // get parameters of each optimization
  double alpha = iap.get<double>("INV_ALPHA");
  double beta = iap.get<double>("INV_BETA");
  int max_itter = iap.get<int>("INV_ANA_MAX_RUN");

  gcomm->Barrier();

  // plausibility check's
  if (expid < 0)
    dserror(
        "Wrong (negative) exp_id number in dat file! EXP_ID has to be greater or equal than zero");

  if (!nexp_) dserror("number of experiments is zero!");

  if (!nmp_mult_) dserror("number of measured dofs of all experiments is zero!");

  if ((nexp_ * (np_ + 1)) != ngroup)
  {
    std::cout << "Number of experiments (nexp): " << nexp_ << std::endl;
    std::cout << "Number of parameter to optimize (np) +1 : " << (np_ + 1) << std::endl;
    std::cout << "Number of specified groups (ngroup): " << ngroup << std::endl;
    dserror("ngroup has to be nexp*(np+1)!");
  }

  double tol_min = 0.0;
  double tol_max = 0.0;
  gcomm->MinAll(&tol_, &tol_min, 1);
  gcomm->MaxAll(&tol_, &tol_max, 1);
  if (tol_min != tol_max) dserror("Tolerance is not the same on every proc!");

  int max_itter_min = 0;
  int max_itter_max = 0;
  gcomm->MinAll(&max_itter, &max_itter_min, 1);
  gcomm->MaxAll(&max_itter, &max_itter_max, 1);
  if (max_itter_min != max_itter_max) dserror("MAY RUN is not the same on every proc!");

  double alpha_min = 0;
  double alpha_max = 0;
  gcomm->MinAll(&alpha, &alpha_min, 1);
  gcomm->MaxAll(&alpha, &alpha_max, 1);
  if (alpha_min != alpha_max) dserror("Alpha is not the same in every proc!");

  double beta_min = 0;
  double beta_max = 0;
  gcomm->MinAll(&beta, &beta_min, 1);
  gcomm->MaxAll(&beta, &beta_max, 1);
  if (beta_min != beta_max) dserror("Beta is not the same in every proc!");

  int check_neg_params_min;
  int check_neg_params_max;
  gcomm->MinAll(&check_neg_params_, &check_neg_params_min, 1);
  gcomm->MaxAll(&check_neg_params_, &check_neg_params_max, 1);
  if (check_neg_params_min != check_neg_params_max)
    dserror("check_neg_params_ is not the same in every proc!");


  if (nexp_ > 1)
  {
    // every proc gets a different expid, this is the ID that identifies the different experiments
    int expidtmp = -1;
    for (int i = 1; i < nexp_ + 1; ++i)
    {
      if (groupid < i * (np_ + 1) && groupid >= (i - 1) * (np_ + 1)) expidtmp = i - 1;
    }
    if (expidtmp != expid)
      dserror(
          "The files were not read in properly!\n"
          "The experiment with ID number 0 has to be read in first etc.\n"
          "Please use following scheme as parameters: Input_1 Output_1 Input_1 Output_2 Input_1 "
          "Output_3 Input_1 Output_4 Input_2 Output_1 Input_2 Output_2 etc.");
  }
  return;
}

/*----------------------------------------------------------------------*/
/* Constraint optimization:
 * In viscofractional model: alpha has to be always below 1
 * otherwise: do not accept optimization step               abirzle 12/17
 */
void STR::GenInvAnalysis::ConstrainAlpha()
{
  if (p_(viscofract_alpha_position_) > 1.0)
  {
    std::cout << "WARNING: alpha > 1: reverse update" << std::endl;
    // reset the update step in case it was too large and increase regularization;
    // storage is not touched, so screen print out does not reflect what really happens
    p_ = p_o_;
    error_ = error_o_;
    error_grad_ = error_grad_o_;
    mu_ = mu_o_ * 2.0;
    return;
  }
}

/*----------------------------------------------------------------------*/
/* Read Monitor file (moved into own method)                HaWi 07/18  */
/*----------------------------------------------------------------------*/

void STR::GenInvAnalysis::ReadMonitorDofBased(const int myrank)
{
  ndofs_ = 0;

  // input parameters inverse analysis
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  {
    char* foundit = NULL;
    std::string monitorfilename = iap.get<std::string>("MONITORFILE");
    if (monitorfilename == "none.monitor") dserror("No monitor file provided");
    // insert path to monitor file if necessary
    if (monitorfilename[0] != '/')
    {
      std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos != std::string::npos)
      {
        std::string path = filename.substr(0, pos + 1);
        monitorfilename.insert(monitorfilename.begin(), path.begin(), path.end());
      }
    }

    FILE* file = fopen(monitorfilename.c_str(), "rb");
    if (file == NULL) dserror("Could not open monitor file %s", monitorfilename.c_str());

    char buffer[150000];
    DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
    // read steps
    foundit = strstr(buffer, "steps");
    foundit += strlen("steps");
    nsteps_ = strtol(foundit, &foundit, 10);
    timesteps_.resize(nsteps_);
    // read nnodes
    foundit = strstr(buffer, "nnodes");
    foundit += strlen("nnodes");
    nnodes_ = strtol(foundit, &foundit, 10);
    // read nodes
    nodes_.resize(nnodes_);
    dofs_.resize(nnodes_);
    for (int i = 0; i < nnodes_; ++i)
    {
      DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
      foundit = buffer;
      nodes_[i] = strtol(foundit, &foundit, 10);
      int ndofs = strtol(foundit, &foundit, 10);
      ndofs_ += ndofs;
      dofs_[i].resize(ndofs, -1);
      if (!myrank) printf("Monitored node %d ndofs %d dofs ", nodes_[i], (int)dofs_[i].size());
      for (int j = 0; j < ndofs; ++j)
      {
        dofs_[i][j] = strtol(foundit, &foundit, 10);
        if (!myrank) printf("%d ", dofs_[i][j]);
      }
      if (!myrank) printf("\n");
    }

    // total number of measured dofs
    nmp_ = ndofs_ * nsteps_;

    // read in measured curve

    {
      mcurve_ = Epetra_SerialDenseVector(nmp_);

      if (!myrank) printf("nsteps %d ndofs %d\n", nsteps_, ndofs_);

      // read comment lines
      foundit = buffer;
      DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
      while (strstr(buffer, "#"))
        DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);

      // read in the values for each node in dirs directions
      int count = 0;
      for (int i = 0; i < nsteps_; ++i)
      {
        // read the time step
        timesteps_[i] = strtod(foundit, &foundit);
        for (int j = 0; j < ndofs_; ++j) mcurve_[count++] = strtod(foundit, &foundit);
        DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
        foundit = buffer;
      }
      if (count != nmp_) dserror("Number of measured disps wrong on input");
    }
  }
}

/*----------------------------------------------------------------------*/
/* Read Monitor file based on measured points               HaWi 07/18  */
/*----------------------------------------------------------------------*/

void STR::GenInvAnalysis::ReadMonitorPointBased(int myrank)
{
  /* this is how a point based monitor file should look like:

steps 2 npoints 5
2 0 1
2 0 1
2 0 1
2 0 1
2 0 1
# any number of comment lines, but only at this position
# x and y coords at 5 points, x and y coords at corresponding points (x' y') in measurement
# direction (towards interface)
(first time point) x1 y1 x1' y1' x2 y2 x2' y2' x3 y3 x3' y3' x4 y4 x4' y4' x5 y5 x5' y5'
(second time point) x1 y1 x1' y1' x2 y2 x2' y2' x3 y3 x3' y3' x4 y4 x4' y4' x5 y5 x5' y5'
*/

  // the following is a minimal example, mind that Levenberg Marquardt needs
  // steps*npoints > number of parameters to fit
  /*
  steps 1 npoints 1
  2 0 1
  # for example with measurement downwards in vertical (y)-direction for
  # one point (npoints 1) in 2D at time t=1.0 (steps 1)
  1.0 2.000000e-02 4.000000e-02 2.000000e-02 3.000000e-02
  */
  ndofs_ = 0;

  // input parameters inverse analysis
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();
  {
    char* foundit = NULL;
    std::string monitorfilename = iap.get<std::string>("MONITORFILE");
    if (monitorfilename == "none.monitor") dserror("No monitor file provided");
    // insert path to monitor file if necessary
    if (monitorfilename[0] != '/')
    {
      std::string filename = DRT::Problem::Instance()->OutputControlFile()->InputFileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos != std::string::npos)
      {
        std::string path = filename.substr(0, pos + 1);
        monitorfilename.insert(monitorfilename.begin(), path.begin(), path.end());
      }
    }

    FILE* file = fopen(monitorfilename.c_str(), "rb");
    if (file == NULL) dserror("Could not open monitor file %s", monitorfilename.c_str());

    char buffer[150000];
    DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
    // read steps
    foundit = strstr(buffer, "steps");
    foundit += strlen("steps");
    nsteps_ = strtol(foundit, &foundit, 10);
    timesteps_.resize(nsteps_);
    // read npoints, nnodes_ could be replicated to npoints_ f.e. but is reused here
    foundit = strstr(buffer, "npoints");
    foundit += strlen("npoints");
    nnodes_ = strtol(foundit, &foundit, 10);

    // no nodes involved for point based measurement
    nodes_.resize(0);

    dofs_.resize(nnodes_);
    for (int i = 0; i < nnodes_; ++i)
    {
      DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
      foundit = buffer;
      int ndofs = strtol(foundit, &foundit, 10);
      ndofs_ += ndofs;
      dofs_[i].resize(ndofs, -1);
      if (!myrank) printf("Monitored point %d ndofs %d dofs", i, (int)dofs_[i].size());
      for (int j = 0; j < ndofs; ++j)
      {
        dofs_[i][j] = strtol(foundit, &foundit, 10);
        if (!myrank) printf(" %d", dofs_[i][j]);
      }
      if (!myrank) printf("\n");
    }

    // total number of measured coordinates
    nmp_ = nnodes_ * nsteps_;

    // read in measured points

    {
      mcurve_ = Epetra_SerialDenseVector(2 * ndofs_ * nsteps_);

      if (!myrank) printf("nsteps %d ndofs %d\n", nsteps_, ndofs_);

      // read comment lines
      foundit = buffer;
      DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
      while (strstr(buffer, "#"))
        DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);

      // read in the values for each node in dirs directions
      int count = 0;
      for (int i = 0; i < nsteps_; ++i)
      {
        // read the time step
        timesteps_[i] = strtod(foundit, &foundit);
        for (int j = 0; j < (2 * ndofs_); ++j) mcurve_[count++] = strtod(foundit, &foundit);
        DRT::UTILS::Checkfgets(fgets(buffer, 150000, file), file, monitorfilename);
        foundit = buffer;
      }
      if (count != 2 * ndofs_ * nsteps_) dserror("Number of measured coords wrong on input");
    }
  }
}
