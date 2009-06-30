/*----------------------------------------------------------------------*/
/*!
 * \file inv_analysis.cpp

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/rausch
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "inv_analysis.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/global_inp_control2.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_io/io_hdf.H"
#include "../drt_lib/linalg_ana.H"
#include "../drt_mat/material.H"
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_structure/strtimint_create.H"
#include "../drt_matelast/elast_coupanisoexpotwo.H"
#include "../drt_matelast/elast_coupanisoneohooketwo.H"
#include "../drt_matelast/elast_coupblatzko.H"
#include "../drt_matelast/elast_couplogneohooke.H"
#include "../drt_matelast/elast_isomooneyrivlin.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_isoyeoh.H"
#include "../drt_matelast/elast_summand.H"
#include "../drt_matelast/elast_vologden.H"
#include "../drt_matelast/elast_volsussmanbathe.H"


//using namespace LINALG::ANA;
using namespace std;
using namespace DRT;
using namespace MAT;


#include "../drt_structure/stru_resulttest.H"


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

  if (surfneum_.size()==1)
  {
    cout << "Load controlled inverse analysis." << endl;
    problem_type_ = 0;
    bc_nodes_ = *(surfneum_[0]->Nodes());
  }
  else if (surfneum_.size()==0)
  {
    cout << "Displacment controlled inverse analysis." << endl;
    problem_type_ = 1;
    bc_nodes_ = *(surfdir_[1]->Nodes());
  }
  else dserror("The inverse analysis only works for one Neumann or two Dirichlet Conditions!");

  // input parameters structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // measured points, gives the number how many displacment steps are
  // measured
  mp_ = sdyn.get<int>("NUMSTEP");
  tstep_ = sdyn.get<double>("TIMESTEP");

  // input parameters inverse analysis
  const Teuchos::ParameterList& iap = DRT::Problem::Instance()->InverseAnalysisParams();

  //  tolerance for the curve fitting
  tol_ = iap.get<double>("INV_ANA_TOL");

  // displacment vectors
  ccurve_o_.Resize(mp_);                // calculated displacment
  ccurve_.Resize(mp_);                  // calculated displacment of the previous run
  mcurve_.Resize(mp_);                  // measured displacment of the experiment
  rcurve_.Resize(mp_);                  // difference between the measured and the calculated displacment

  //set load/displacment curve
  double cp0 = iap.get<double>("MEASURED_CURVE0");
  double cp1 = iap.get<double>("MEASURED_CURVE1");
  double cp2 = iap.get<double>("MEASURED_CURVE2");
  // build measured curve
  for (int i=0; i<mp_; i++)
    mcurve_[i] = cp0*(1-exp(-pow((cp1*(500.0/mp_)*i), cp2)));

  // error: diference of the measured to the calculated curve
  error_  = 1.0E6;
  error_o_= 1.5E6;


  // trainings parameter
  mu_  = 1;                                     // start value
  mum_ = 0.1;
  mup_ = 10;

  // material parameters
  //p_.Resize(3);
  //p_o_.Resize(3);

  // Which material is used in the input file
  {
    Teuchos::RCP<const MAT::Material> material = discret_->lRowElement(0)->Material();
    if (material->MaterialType() == INPAR::MAT::m_lung_penalty)
    {
      const MAT::LungPenalty* actmat = static_cast<const MAT::LungPenalty*>(material.get());
      int j = p_.Length();
      p_.Resize(j+3);
      p_(j)   = actmat->C();
      p_(j+1) = actmat->K1();
      p_(j+2) = actmat->K2();
    }
    else if (material->MaterialType() == INPAR::MAT::m_lung_ogden)
    {
      const MAT::LungOgden* actmat = static_cast<const MAT::LungOgden*>(material.get());
      int j = p_.Length();
      p_.Resize(j+3);
      p_(j)   = actmat->C();
      p_(j+1) = actmat->K1();
      p_(j+2) = actmat->K2();
    }
    else if (material->MaterialType() == INPAR::MAT::m_elasthyper)
    {
      // Create a pointer on the Material
      const MAT::ElastHyper* actmat = static_cast<const MAT::ElastHyper*>(material.get());

      // For each of the summed up materials we need to add the
      // parameters to the inverse analysis
      cout << endl;
      cout << "Anzahl der Materialien: " << actmat->NumMat() << endl;

      for (int i=0; i< actmat->NumMat(); i++)
      {
        cout << " Run for the " << i << ". material" << endl;
        //get the material of the summand
        Teuchos::RCP< MAT::ELASTIC::Summand > summat = MAT::ELASTIC::Summand::Factory(actmat->MatID(i));

        switch (summat->MaterialType())
        {
          case INPAR::MAT::mes_couplogneohooke:
          {
            const MAT::ELASTIC::CoupLogNeoHooke* actmat2 = static_cast<const MAT::ELASTIC::CoupLogNeoHooke*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+5);
            p_[j]   = actmat2->Mue();
            p_[j+1] = actmat2->Lambda();
            p_[j+2] = actmat2->Parmode();
            p_[j+3] = actmat2->Youngs();
            p_[j+4] = actmat2->Nue();
            break;
          }
          case INPAR::MAT::mes_coupblatzko:
          {
            const MAT::ELASTIC::CoupBlatzKo* actmat2 = static_cast<const MAT::ELASTIC::CoupBlatzKo*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+3);
            p_[j]   = actmat2->Mue();
            p_[j+1] = actmat2->Nue();
            p_[j+2] = actmat2->F();
            break;
          }
          case INPAR::MAT::mes_isoneohooke:
          {
            const MAT::ELASTIC::IsoNeoHooke* actmat2 = static_cast<const MAT::ELASTIC::IsoNeoHooke*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+1);
            p_[j]   = actmat2->Mue();
            break;
          }
          case INPAR::MAT::mes_isoyeoh:
          {
            const MAT::ELASTIC::IsoYeoh* actmat2 = static_cast<const MAT::ELASTIC::IsoYeoh*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+3);
            p_[j]   = actmat2->C1();
            p_[j+1] = actmat2->C2();
            p_[j+2] = actmat2->C3();
            break;
          }
          case INPAR::MAT::mes_isomooneyrivlin:
          {
            const MAT::ELASTIC::IsoMooneyRivlin* actmat2 = static_cast<const MAT::ELASTIC::IsoMooneyRivlin*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+2);
            p_[j]   = actmat2->C1();
            p_[j+1] = actmat2->C2();
            break;
          }
          case INPAR::MAT::mes_volsussmanbathe:
          {
            const MAT::ELASTIC::VolSussmanBathe* actmat2 = static_cast<const MAT::ELASTIC::VolSussmanBathe*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+1);
            p_[j]   = actmat2->Kappa();
            break;
          }
          case INPAR::MAT::mes_vologden:
          {
            const MAT::ELASTIC::VolOgden* actmat2 = static_cast<const MAT::ELASTIC::VolOgden*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+2);
            p_[j]   = actmat2->Kappa();
            p_[j+1] = actmat2->Beta();
            break;
          }
          case INPAR::MAT::mes_coupanisoexpotwo:
          {
            const MAT::ELASTIC::CoupAnisoExpoTwo* actmat2 = static_cast<const MAT::ELASTIC::CoupAnisoExpoTwo*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+4);
            p_[j]   = actmat2->K1();
            p_[j+1] = actmat2->K2();
            p_[j+2] = actmat2->K3();
            p_[j+3] = actmat2->K4();
            break;
          }
          case INPAR::MAT::mes_coupanisoneohooketwo:
          {
            const MAT::ELASTIC::CoupAnisoNeoHookeTwo* actmat2 = static_cast<const MAT::ELASTIC::CoupAnisoNeoHookeTwo*>(summat.get());
            int j = p_.Length();
            p_.Resize(j+2);
            p_[j]   = actmat2->C1();
            p_[j+1] = actmat2->C2();
            break;
          }
          default:
            dserror("cannot deal with this material");
        }
      }

    }
    else dserror("The inverse analysis is not implemented for this material");
  }

  for (int i=0; i<p_.Length();i++)
  {
    cout << p_[i] << " ";
  }
  cout << endl;

  pp_ = p_.Length();
  p_neg_.Resize(pp_);
  delta_p_.Resize(pp_);
  for (int i=0;  i< pp_; i++)
  {
    p_neg_ [i] = p_[i];
    delta_p_[i] = 1;
  }

  // controlling parameter
  precond_ =   0;         // counter where in the preconditioning we are
  numb_run_ =  0;          // counter of how many runs were made in the inverse analysis
  direction_ = 0;       // gives the direction of the stretch default 0

  // storage
  ccurve_s_.Shape(mp_, 1);
  for (int i=0; i<mp_; i++)
    ccurve_s_(i, 0) = mcurve_[i];

  par_s_.Reshape(4,  numb_run_+2);
  par_s_(0, 0) = mcurve_[int(mp_*0.05)];
  par_s_(1, 0) = mcurve_[int(mp_*0.1)];
  par_s_(2, 0) = mcurve_[int(mp_*0.5)];
  par_s_(3, 0) = mcurve_[mp_-1];

}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::InvAnalysis::Integrate()
{
  int max_itter = 100000;

  // get input parameter lists
  //  = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags
    = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn
    = DRT::Problem::Instance()->StructuralDynamicParams();
  Teuchos::ParameterList xparams;
  xparams.set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
  // fitting loop
  do
  {
    // create time integrator
    sti_ = TimIntCreate(ioflags, sdyn, xparams, discret_, solver_, output_);
    if (sti_ == Teuchos::null) dserror("Failed in creating integrator.");

    // initialize time loop
    double timen = sti_->GetTime();
    const double timemax = sti_->GetTimeEnd();
    int stepn = sti_->GetStep();
    const int stepmax = sti_->GetTimeNumStep();

    // time loop
    while ( (timen <= timemax) && (stepn <= stepmax) )
    {
      // integrate time step
      // after this step we hold disn_, etc
      sti_->IntegrateStep();

      // update displacements, velocities, accelerations
      // after this call we will have disn_==dis_, etc
      sti_->UpdateStepState();

      // gets the displacments per timestep
      if (problem_type_==0)
      {
        // current displacements
        Teuchos::RCP<Epetra_Vector> dis = sti_->DisNew();
        // external forces
        Teuchos::RCP<Epetra_Vector> fext = sti_->Fext();
        Teuchos::RCP<Epetra_Vector> fextunit = Teuchos::rcp(new Epetra_Vector(*fext));
        for (int lid=0; lid<fextunit->MyLength(); ++lid)
        {
          if (fabs((*fextunit)[lid]) < EPS12)
            ; /* do nothing */
          else
            (*fextunit)[lid] = 1.0;
        }
        Teuchos::RCP<LINALG::MapExtractor> nbcmaps = LINALG::ConvertDirichletToggleVectorToMaps(fextunit);
        Teuchos::RCP<Epetra_Vector> nbcdis = nbcmaps->ExtractCondVector(dis);
        // Displacments at the end of the time step are saved anyway
        GetCalculatedCurve(nbcdis, stepn-1);
      }
      else
      {
        // current reaction forces
        Teuchos::RCP<Epetra_Vector> freact = sti_->Freact();
        // Loads have to be taken out of the predictor/correcter during the calculation
        GetCalculatedCurve(freact, stepn-1);
      }

      // update time and step
      sti_->UpdateStepTime();

      // print info about finished time step
      sti_->PrintStep();

      // write output
      sti_->OutputStep();

      // get current time ...
      timen = sti_->GetTime();
      // ... and step
      stepn = sti_->GetStep();
    }

    // shift curve
    for (int i=mp_-1; i>-1; i--)
      ccurve_[i]=abs(ccurve_[i-1]);
    ccurve_[0]=0.0;

    //
    Evaluate();

  } while (error_>tol_ && numb_run_<max_itter);

  // leave
  return;
}

/*----------------------------------------------------------------------*/
/* evaluate */
void STR::InvAnalysis::Evaluate()
{
  discret_->Comm().Barrier();
  for (int proc=0; proc<discret_->Comm().NumProc(); ++proc)
  {
    if (proc==discret_->Comm().MyPID())
    {
      if (proc == 0)
      {
        cout << "********************************************************************************************" << endl;
        cout << "*******************************\t Inverse Analysis \t************************************" << endl;
        cout << "*******************************\t\t run \t\t************************************" << endl;
        cout << "*******************************\t\t  " << numb_run_ <<"\t\t************************************" << endl;
        cout << "********************************************************************************************" << endl;

        if (precond_<2 && error_ > 3*tol_)
          PrecondParameters();
        else
          CalcNewParameters();
      }
    }
  }
  discret_->Comm().Barrier();

  ResetParameters();
  discret_->Comm().Broadcast(&p_[0],  3,  0);
  numb_run_++;

  // write new material parameter
  {
    Teuchos::RCP<MAT::Material> material = discret_->lRowElement(0)->Material();
    if (material->MaterialType() == INPAR::MAT::m_lung_penalty)
    {
      MAT::LungPenalty* actmat = static_cast<MAT::LungPenalty*>(material.get());
      actmat->SetC(p_(0)*p_(0));
      actmat->SetK1(p_(1)*p_(1));
      actmat->SetK2(p_(2)*p_(2));
    }
    else if (material->MaterialType() == INPAR::MAT::m_lung_ogden)
    {
      MAT::LungOgden* actmat = static_cast<MAT::LungOgden*>(material.get());
      actmat->SetC(p_(0)*p_(0));
      actmat->SetK1(p_(1)*p_(1));
      actmat->SetK2(p_(2)*p_(2));
    }
    else
      dserror("Cannot handle material with type %d", material->MaterialType());
  }

  // wonderful
  return;
}

/*----------------------------------------------------------------------*/
/*  */
void STR::InvAnalysis::PrecondParameters()
{

  cout << "--------------------------------------------------------------------" << endl;
  cout << "-------------------- preconditioning parameters --------------------" << endl;
  cout << "--------------------------------------------------------------------" << endl;

  // calculating the new parameter
  vector<int>    matpar(3);
  matpar[0] = 0;
  matpar[1] = 0;
  matpar[2] = 0;
  vector<int>    sto(4);
  sto[0]= int(mp_*0.25);
  sto[1]= int(mp_*0.5);
  sto[2]= int(mp_*0.75);

  switch(precond_) {
    case 0:                                             // precondition parameter c
      matpar[0] = 1;
      matpar[1] = 1;
      matpar[2] = 1;
      break;
    case 1:                                           // precondition parameter k2
       matpar[0]   = 0;
       matpar[1]   = 0;
       matpar[2]   = 1;
    break;
    default:
      dserror("Precond_ in Inv_Ana is out of range!");
    break;
  }

  if (ccurve_[mp_-1] == 0.0 && ccurve_o_[mp_-1]!=0.0)
  {
    cout << "New ccurve!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1" << endl;
    for (int i=0; i<mp_;i++)
    {
      ccurve_[i]=ccurve_o_[i];
      ccurve_o_[i] = 0.0;
    }
  }
  // calculate gradient of claculated curve
  double g_ccurve   = ccurve_[sto[1]+1]   - ccurve_[sto[1]-1];
  double g_ccurve_o = ccurve_o_[sto[1]+1] - ccurve_o_[sto[1]-1];
  double g_mcurve   = mcurve_[sto[1]+1]   - mcurve_[sto[1]-1];

  double gg_ccurve   = ccurve_[sto[2]+1]   - 2*ccurve_[sto[2]]   + ccurve_[sto[2]-1];
  double gg_ccurve_o = ccurve_o_[sto[2]+1] - 2*ccurve_o_[sto[2]] + ccurve_o_[sto[2]-1];
  double gg_mcurve   = mcurve_[sto[2]+1]   - 2*mcurve_[sto[2]]   + mcurve_[sto[2]-1];

  vector<double> par(3);
  par[0] = ((ccurve_[sto[0]] - mcurve_[sto[0]])/2)+((ccurve_[sto[2]] - mcurve_[sto[2]])/2);
  par[1] = g_ccurve        - g_mcurve;
  par[2] = gg_ccurve       - gg_mcurve;

  cout << "ccurve_[sto[0]]: " << ccurve_[sto[0]] << endl;
  cout << "mcurve_[sto[0]]: " << mcurve_[sto[0]] << endl;
  cout << "ccurve_[sto[2]]: " << ccurve_[sto[2]] << endl;
  cout << "mcurve_[sto[2]]: " << mcurve_[sto[2]] << endl;
  cout << "Parameter:     " << par[0] << endl;

  vector<double> par_o(3);
  par_o[0] = (ccurve_o_[sto[0]] - mcurve_[sto[0]])/2+(ccurve_o_[sto[2]] - mcurve_[sto[2]])/2;
  par_o[1] = g_ccurve_o        - g_mcurve;
  par_o[2] = gg_ccurve_o       - gg_mcurve;
  cout << "old Parameter: " << par_o[0] << endl;

  vector<double> acc(3);
  acc[0] = 0.025 * (mcurve_[sto[0]]+mcurve_[sto[2]]);
  acc[1] = 0.005 * g_mcurve;
  acc[2] = 0.1   * gg_mcurve;

  if (matpar[0]==1 && (par[0]>0 && par_o[0]<0 || par[0]<0 && par_o[0]>0))
      delta_p_[0] = sqrt(delta_p_[0]);

  if (matpar[1]==1 && (par[1]>0 && par_o[1]<0 || par[1]<0 && par_o[1]>0))
      delta_p_[1] = sqrt(delta_p_[1]);

  if (matpar[2]==1 && (par[2]>0 && par_o[2]<0 || par[2]<0 && par_o[2]>0))
      delta_p_[2] = sqrt(delta_p_[2]);

  // update variables
  PrintStorage();

  // update parameters
  for (int i=0; i<3; i++) {
    if (matpar[i]==1) {
      if (problem_type_==0) {
        if   (par[i]<0) p_[i] /= delta_p_[i];
        else            p_[i] *= delta_p_[i];}
      else{
        if   (par[i]<0) p_[i] *= delta_p_[i];
        else            p_[i] /= delta_p_[i];}}}


  cout << "new Parameter: " << par[0] << endl;

  // artivicial variation for the start of the LMA
  if (precond_==2) {
    for (unsigned int i=0; i<3; i++)
      p_o_(i) = p_o_(i)*1.5;
    for (int i=0; i<mp_; i++)
      ccurve_o_(i) = ccurve_o_(i)*1.1;}

  par_s_.Reshape(4,  numb_run_+2);
  par_s_(0, numb_run_+1) = ccurve_[sto[0]];
  par_s_(1, numb_run_+1) = ccurve_[sto[1]];
  par_s_(2, numb_run_+1) = ccurve_[sto[2]];
  par_s_(3, numb_run_+1) = ccurve_[mp_];
  cout << endl;
  for (int i=0; i<4; i++){
    for (int j=0; j<numb_run_+2; j++)
      printf(" %5.1f ",  par_s_(i, j));
    printf("\n"); }

//      p_=p_o_;
//      if (numb_run_==0)
//      {
//        p_[0] = 0.1 * p_[0];
//      }
//      if (numb_run_==1)
//      {
//        p_[0] = 100 * p_[0];
//      }
//      if (numb_run_==2)
//      {
//        p_[0] = 0.1 * p_[0];
//        p_[1] = 0.1 * p_[1];
//      }
//      if (numb_run_==3)
//      {
//        p_[1] = 100 * p_[1];
//      }
//      if (numb_run_==4)
//      {
//        p_[1] = 0.1 * p_[1];
//        p_[2] = 0.1 * p_[2];
//      }
//      if (numb_run_==5)
//      {
//        p_[2] = 100 * p_[2];
//      }

//      if (numb_run_==6)
//      {
//        precond_=6;
//      }

  return;
}

/*----------------------------------------------------------------------*/
/* */
void STR::InvAnalysis::CalcNewParameters()
{
  cout << "---------------- normal calculation -----------------------" << endl;
  // initalization of the Jacobi and storage matrix
  Epetra_SerialDenseMatrix J(mp_, pp_);
  Epetra_SerialDenseMatrix sto(pp_,  pp_);
  Epetra_SerialDenseMatrix sto2(pp_, mp_);
  Epetra_SerialDenseVector delta_p(pp_);

  //calculating J(p)
  for (int i=0; i<mp_; i++)
    for (int j=0; j<pp_; j++)
      J(i, j) = (ccurve_[i]-ccurve_o_[i])/(p_[j]-p_o_[j]);

  //calculating J.T*J
  sto.Multiply('T',  'N',  1,  J, J,  0);

  //calculating J.T*J+mu*I
  for (int i=0; i<pp_; i++)
    sto[i][i] +=  mu_;

  //calculating (J.T*J+mu*I).I
  LINALG::NonSymmetricInverse(sto,  pp_);

  //calculating (J.T*J+mu*I).I*J.T
  sto2.Multiply('N', 'T', 1,  sto, J, 0);

  //calculating (J.T*J+mu*I).I*J.T*rcurve_
  delta_p_.Multiply('N', 'N', 1, sto2, rcurve_, 0);

  PrintStorage();

  for (int i=0; i<pp_; i++)
  {
    if ( p_(i)-delta_p_(i) > 0)
    {
      p_(i)   = p_(i) - delta_p_(i);
      cout << "p_ - delta_p_ = " << p_(i) << " - " << delta_p_(i)
           << " = " << p_(i)-delta_p_(i) << endl;
    }
    else
    {
      p_(i)=abs(1/(p_(i) - delta_p_(i)-1));
      cout << "abs(1/(p_ - delta_p_-2)) = abs(1/(" << p_(i) << " - " << delta_p_(i)
           << " - 2)) = " << abs(1/(p_(i)-delta_p_(i)-2)) << endl;
    }
  }
  //p_(i)   = p_(i) - delta_p_(i);

  if (error_o_<error_) mu_ = mup_;
  else  mu_ = mum_;

  cout << "********************************************ENDE********************************************" << endl;
  return;
}

/*----------------------------------------------------------------------*/
/* */
void STR::InvAnalysis::GetCalculatedCurve(const Teuchos::RCP<Epetra_Vector> value, int numb)
{
//   // check in which direction the surface neumann conditions pulls
//   if (value->Map().MyGID(0) && numb==1)
//   {
//     if (abs((*value)[bc_nodes_[0]]) > abs((*value)[bc_nodes_[1]]) &&
//         abs((*value)[bc_nodes_[0]]) > abs((*value)[bc_nodes_[2]]))
//       direction_ = 0;
//     else if (abs((*value)[bc_nodes_[1]]) > abs((*value)[bc_nodes_[0]]) &&
//              abs((*value)[bc_nodes_[1]]) > abs((*value)[bc_nodes_[2]]))
//       direction_ = 1;
//     else if (abs((*value)[bc_nodes_[2]]) > abs((*value)[bc_nodes_[0]]) &&
//              abs((*value)[bc_nodes_[2]]) > abs((*value)[bc_nodes_[1]]))
//       direction_ = 2;
//     else direction_=0;//dserror("The displacment direction is not clear!");
//   }
//   discret_->Comm().Broadcast(&direction_, 1, 0);
//   cout << "Richtung: " << direction_ << endl;

  // blank residual at DOFs on Dirichlet BC
  if (problem_type_==0)                        // load controlled
  {
    value->MeanValue(&(ccurve_[numb]));
  }
  else                                         // disp controlled
  {
    for (unsigned int i=0; i<value->MyLength();i++)
    {
      ccurve_[numb]+=(*value)[i];
    }
    ccurve_[numb] /= 2;

    //this just works for symmetric problems because then the other two
    //directions neutralize themselfs
  }

  return;
}

/*----------------------------------------------------------------------*/
/* reset the values for next run */
void STR::InvAnalysis::ResetParameters()
{
  cout << "Resetting all the Parameters! " << endl;
  // reset calculated curve values
  ccurve_.Scale(0.0);

  //double iotime = params_.get<double>("total time",0.0);
  //int iostep = params_.get<int>("step",0);
  output_->NewResultFile((numb_run_));
  //output_.WriteMesh(iostep, iotime);

  return;
}

/*----------------------------------------------------------------------*/
/* */
void STR::InvAnalysis::PrintStorage()
{
   // update everything for the next run

  for (int i=0; i< mp_; i++)
    rcurve_[i] = ccurve_[i]-mcurve_[i];

  error_o_ = error_;
  error_ = rcurve_.Norm1();

  // store the error and mu_

  p_s_.Reshape(numb_run_+1,  pp_);
  for (int i=0; i<pp_; i++)
    p_s_(numb_run_, i)=p_(i);

  delta_p_s_.Reshape(numb_run_+1,  pp_);
  for (int i=0; i<pp_; i++)
    delta_p_s_(numb_run_, i)=delta_p_(i);

  ccurve_s_.Reshape(mp_,  numb_run_+2);
  for (int i=0; i<mp_; i++)
    ccurve_s_(i, numb_run_+1)=ccurve_[i];

  mu_s_.Resize(numb_run_+1);
  mu_s_(numb_run_)=mu_;

  error_s_.Resize(numb_run_+1);
  error_s_(numb_run_) = error_;

  ccurve_o_ = ccurve_;
  p_o_ = p_;

  // print error and parameter
  cout << endl;
  printf("error \t\t c \t\t\t k1 \t\t\t k2 \t\t\t delta_c \t\t delta_k1 \t\t delta_k2 \n");
  for (int i=0; i < numb_run_+1; i++)
    printf("%5.2f  %20.10f  %20.10f  %20.10f  %20.10f  %20.10f  %20.10f  %20.10f\n",
           error_s_(i), p_s_(i, 0),  p_s_(i, 1), p_s_(i, 2), p_s_(i, 1)/p_s_(i, 2),
           delta_p_s_(i, 0), delta_p_s_(i, 1), delta_p_s_(i, 2));
  for (int i=0; i<mp_; i++)
  {
    for (int j=0; j<numb_run_+2; j++)
    {
      printf(" %10.1f ",  ccurve_s_(i, j));
    }
    printf("\n");
  }
}



#endif
