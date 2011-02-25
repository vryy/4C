/*----------------------------------------------------------------------*/
/*!
 * \file mlmc.cpp

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de

</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "mlmc.H"
#include <ctime>
#include <cstdlib>
#include <iostream>
#include "Epetra_SerialDenseMatrix.h"
#include "../drt_lib/global_inp_control2.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_io/io_hdf.H"
#include "../drt_mat/material.H"
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/neohooke.H"
#include "../drt_structure/strtimint_create.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_matelast/elast_coupanisoexpotwo.H"
#include "../drt_matelast/elast_coupanisoneohooketwo.H"
#include "../drt_matelast/elast_coupblatzko.H"
#include "../drt_matelast/elast_couplogneohooke.H"
#include "../drt_matelast/elast_isoexpo.H"
#include "../drt_matelast/elast_isomooneyrivlin.H"
#include "../drt_matelast/elast_isoneohooke.H"
#include "../drt_matelast/elast_isoyeoh.H"
#include "../drt_matelast/elast_volpenalty.H"
#include "../drt_matelast/elast_vologden.H"
#include "../drt_matelast/elast_volsussmanbathe.H"
#include "../drt_mat/matpar_bundle.H"

using namespace std;
using namespace DRT;
using namespace MAT;


#include "../drt_structure/stru_resulttest.H"


/*----------------------------------------------------------------------*/
/* standard constructor */
STR::MLMC::MLMC(Teuchos::RCP<DRT::Discretization> dis,
                                    Teuchos::RCP<LINALG::Solver> solver,
                                    Teuchos::RCP<IO::DiscretizationWriter> output)
  : discret_(dis),
    solver_(solver),
    output_(output),
    sti_(Teuchos::null)
{
  //int myrank = dis->Comm().MyPID();

  reset_out_count_=0;

  // input parameters structural dynamics
  //const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // input parameters multi level monte carlo
  //const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // read material parameters from input file
  //ReadInParameters();
  // controlling parameter
  numb_run_ =  0;     // counter of how many runs were made monte carlo
}


/*----------------------------------------------------------------------*/
/* analyse */
void STR::MLMC::Integrate()
{
  //int myrank = discret_->Comm().MyPID();

  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();
  int numruns = mlmcp.get<int>("NUMRUNS");
  do
  {
    output_->NewResultFile((numb_run_));


     // get input lists
     const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

     // major switch to different time integrators
     switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn,"DYNAMICTYP"))
     {
     case INPAR::STR::dyna_centr_diff:
       dserror("no central differences in DRT");
       break;
     case INPAR::STR::dyna_gen_alfa:
     case INPAR::STR::dyna_gen_alfa_statics:
     case INPAR::STR::dyna_statics:
     case INPAR::STR::dyna_genalpha:
     case INPAR::STR::dyna_onesteptheta:
     case INPAR::STR::dyna_gemm:
     case INPAR::STR::dyna_ab2:
     case INPAR::STR::dyna_euma :
     case INPAR::STR::dyna_euimsto :
       dyn_nlnstructural_drt();
       break;
     case INPAR::STR::dyna_Gen_EMM:
       dserror("GEMM not supported");
       break;
     default:
       dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
     }
     numb_run_++;
  } while (numb_run_< numruns);
  MapResultsToFinestGrid();
  return;
}
//---------------------------------------------------------------------------------------------
//void STR::MLMC::MapResultsToFinestGrid(Teuchos::RCP<DRT::Discretization> dis_in,
 //                                               Teuchos::RCP<DRT::Discretization> dis_out,
 //                                               int run_id)
void STR::MLMC::MapResultsToFinestGrid()
{
  // Get finest Grid problem instance
  // access the discretization
   Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
   actdis = DRT::Problem::Instance(1)->Dis(genprob.numsf, 0);

   // set degrees of freedom in the discretization
   if (not actdis->Filled()) actdis->FillComplete();

   // context for output and restart
   Teuchos::RCP<IO::DiscretizationWriter> output2
     = Teuchos::rcp(new IO::DiscretizationWriter(actdis));
   // write mesh to file
   output2->NewResultFile((2321));
   output2->WriteMesh(1, 0.01);
   //output2->NewResultFile((2321));
   output2->NewStep( 1, 0.01);

   // Write some vector to file
   // init vector
   Teuchos::RCP<Epetra_MultiVector> test_vector = Teuchos::rcp(new Epetra_MultiVector(*(actdis->DofRowMap()),1,true));
   // fill Vector with arbitrary values for testing
   for (int i = 0; i< test_vector->MyLength();i++)
   {
     (*test_vector)[0][i] = 13.0 ; //->ReplaceMyValue(0,i,13.0);

   }
   // fill with some data
   output2->WriteVector("displacement",test_vector,output2->dofvector);

  cout << "Feine Disketisierung" << actdis<< endl;
  int num_ele_dis_1 = actdis->NumGlobalElements();
  cout << "num_ele_dis_1 "<< num_ele_dis_1 <<endl;

}


//-----------------------------------------------------------------------------------
void STR::MLMC::ReadInParameters()
{
  const int myrank = discret_->Comm().MyPID();

  // loop all materials in problem
  const map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();

  if (myrank == 0) printf("No. material laws considered : %d\n",mats.size());
  map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    const RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
      case INPAR::MAT::m_aaaneohooke:
      {
        MAT::PAR::AAAneohooke* params = dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int j = p_.Length();
        p_.Resize(j+2);
        p_(j)   = params->youngs_;
        p_(j+1) = params->beta_;
        //p_(j+2) = params->nue_; // need also change resize above to invoke nue
      }
      break;
      case INPAR::MAT::m_neohooke:
      {
        MAT::PAR::NeoHooke* params = dynamic_cast<MAT::PAR::NeoHooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int j = p_.Length();
        p_.Resize(j+2);
        p_(j)   = params->youngs_;
        p_(j+1) = params->poissonratio_;
      }
      break;
      case INPAR::MAT::m_elasthyper:
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat               = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i=0; i<nummat; ++i)
        {
          const int id = (*matids)[i];
          const RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
          switch (actelastmat->Type())
          {
            case INPAR::MAT::mes_isoyeoh:
            {
              MAT::ELASTIC::PAR::IsoYeoh* params = dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const int j = p_.Length();
              p_.Resize(j+3);
              p_(j)   = params->c1_;
              p_(j+1) = params->c2_;
              p_(j+2) = params->c3_;
            }
            break;
            case INPAR::MAT::mes_vologden:
            {
              MAT::ELASTIC::PAR::VolOgden* params = dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const int j = p_.Length();
              p_.Resize(j+1);
              p_(j)   = params->kappa_;
              // p_(j+1) = params->beta_; // need also change resize above to invoke beta_
            }
            break;
            default:
              dserror("Unknown type of elasthyper material");
            break;
          }
        }
      }
      case INPAR::MAT::mes_isoyeoh: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      case INPAR::MAT::mes_vologden: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      default:
        // ignore unknown materials ?
        dserror("Unknown type of material");
      break;
    }
  }
  return;
}
//--------------------------------------------------------------------------------------
void STR::MLMC::SetParameters(Epetra_SerialDenseVector p_cur)
{
  const int myrank = discret_->Comm().MyPID();

  // parameters are evaluated on proc 0 only? This should not be neccessary....
  discret_->Comm().Broadcast(&p_cur[0],np_,0);

  // loop all materials in problem
  const map<int,RCP<MAT::PAR::Material> >& mats = *DRT::Problem::Instance()->Materials()->Map();
  int count=0;
  map<int,RCP<MAT::PAR::Material> >::const_iterator curr;
  for (curr=mats.begin(); curr != mats.end(); ++curr)
  {
    const RCP<MAT::PAR::Material> actmat = curr->second;
    switch(actmat->Type())
    {
      case INPAR::MAT::m_aaaneohooke:
      {
        MAT::PAR::AAAneohooke* params = dynamic_cast<MAT::PAR::AAAneohooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        // This is a tiny little bit brutal!!!
        const_cast<double&>(params->youngs_) = p_cur[count];
        const_cast<double&>(params->beta_)   = p_cur[count+1];
        //const_cast<double&>(params->nue_)    = p_cur[count+2];
        if (myrank == 0) printf("MAT::PAR::AAAneohooke %20.15e %20.15e\n",p_cur[count],p_cur[count+1]);
        count += 2;
      }
      break;
      case INPAR::MAT::m_neohooke:
      {
        MAT::PAR::NeoHooke* params = dynamic_cast<MAT::PAR::NeoHooke*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        // This is a tiny little bit brutal!!!
        const_cast<double&>(params->youngs_)       = p_cur[count];
        const_cast<double&>(params->poissonratio_) = p_cur[count+1];
        if (myrank == 0) printf("MAT::PAR::NeoHooke %20.15e %20.15e\n",params->youngs_,params->poissonratio_);
        count += 2;
      }
      break;
      case INPAR::MAT::m_elasthyper:
      {
        MAT::PAR::ElastHyper* params = dynamic_cast<MAT::PAR::ElastHyper*>(actmat->Parameter());
        if (!params) dserror("Cannot cast material parameters");
        const int nummat               = params->nummat_;
        const std::vector<int>* matids = params->matids_;
        for (int i=0; i<nummat; ++i)
        {
          const int id = (*matids)[i];
          const RCP<MAT::PAR::Material> actelastmat = mats.find(id)->second;
          switch (actelastmat->Type())
          {
            case INPAR::MAT::mes_isoyeoh:
            {
              MAT::ELASTIC::PAR::IsoYeoh* params = dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const_cast<double&>(params->c1_) = p_cur[count];
              const_cast<double&>(params->c2_) = p_cur[count+1];
              const_cast<double&>(params->c3_) = p_cur[count+2];
              if (myrank == 0) printf("MAT::ELASTIC::PAR::IsoYeoh %20.15e %20.15e %20.15e\n",params->c1_,params->c2_,params->c3_);
              count += 3;
            }
            break;
            case INPAR::MAT::mes_vologden:
            {
              MAT::ELASTIC::PAR::VolOgden* params = dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(actelastmat->Parameter());
              if (!params) dserror("Cannot cast material parameters");
              const_cast<double&>(params->kappa_) = p_cur[count];
              //const_cast<double&>(params->beta_) = p_cur[count+1];
              if (myrank == 0) printf("MAT::ELASTIC::PAR::VolOgden %20.15e %20.15e\n",params->kappa_,params->beta_);
              count += 1;
            }
            break;
            default:
              dserror("Unknown type of elasthyper material");
            break;
          }
        }
      }
      break;
      case INPAR::MAT::mes_isoyeoh: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      case INPAR::MAT::mes_vologden: // at this level do nothing, its inside the INPAR::MAT::m_elasthyper block
      break;
      default:
        // ignore unknown materials ?
        dserror("Unknown type of material");
      break;
    }
  }


  return;
}



#endif
