/*----------------------------------------------------------------------*/
/*!
\file viscogenmax.cpp

\brief
This file contains the Generalized Maxwell model consisting of one spring
in parallel to a sequential branch of a spring and a dashpot. This model can
be applied to any hyperelastic law of the Elasthyper toolbox. The viscous
effect can be applied to any part of the SEF (isotropic,isotropic isochoric,
isotropic volumetric, anisotropic, anisotropic isochoric, anisotropic
volumetric).

The input line should read
MAT 0   MAT_ViscoGenMax   NUMMAT 0 MATIDS  DENS 0 GAMMA 0 RELAX_ISOT_PRINC 0 BETA_ISOT_PRINC 0 RELAX_ISOT_MOD_VOL 0 BETA_ISOT_MOD_VOL 0 RELAX_ISOT_MOD_ISOC 0 BETA_ISOT_MOD_ISOC 0 RELAX_ANISOT_PRINC 0 BETA_ANISOT_PRINC 0 RELAX_ANISOT_MOD_VOL 0 BETA_ANISOT_MOD_VOL 0 RELAX_ANISOT_MOD_ISOC 0 BETA_ANISOT_MOD_ISOC 0 INIT_MODE -1

<pre>
Maintainer: Aline Bel
            brunon@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/


// unnecessary
//#include <vector>
#include "viscogenmax.H"
#include "elasthyper.H"
#include "../drt_matelast/elast_summand.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ViscoGenMax::ViscoGenMax(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: ElastHyper(matdata),
  relax_isot_princ_(matdata->GetDouble("RELAX_ISOT_PRINC")),
  beta_isot_princ_(matdata->GetDouble("BETA_ISOT_PRINC")),
  relax_isot_mod_vol_(matdata->GetDouble("RELAX_ISOT_MOD_VOL")),
  beta_isot_mod_vol_(matdata->GetDouble("BETA_ISOT_MOD_VOL")),
  relax_isot_mod_isoc_(matdata->GetDouble("RELAX_ISOT_MOD_ISOC")),
  beta_isot_mod_isoc_(matdata->GetDouble("BETA_ISOT_MOD_ISOC")),
  relax_anisot_princ_(matdata->GetDouble("RELAX_ANISOT_PRINC")),
  beta_anisot_princ_(matdata->GetDouble("BETA_ANISOT_PRINC")),
  relax_anisot_mod_vol_(matdata->GetDouble("RELAX_ANISOT_MOD_VOL")),
  beta_anisot_mod_vol_(matdata->GetDouble("BETA_ANISOT_MOD_VOL")),
  relax_anisot_mod_isoc_(matdata->GetDouble("RELAX_ANISOT_MOD_ISOC")),
  beta_anisot_mod_isoc_(matdata->GetDouble("BETA_ANISOT_MOD_ISOC"))

{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ViscoGenMax::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ViscoGenMax(this));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

MAT::ViscoGenMaxType MAT::ViscoGenMaxType::instance_;


DRT::ParObject* MAT::ViscoGenMaxType::Create( const std::vector<char> & data )
{
  MAT::ViscoGenMax* viscogenmax = new MAT::ViscoGenMax();
  viscogenmax->Unpack(data);
  return viscogenmax;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ViscoGenMax::ViscoGenMax()
  : viscoparams_(NULL)
{
  isinit_=false;
  histstressisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodisocurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodisolast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisocurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisolast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodvolcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodvollast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvolcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvollast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressanisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressanisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ViscoGenMax::ViscoGenMax(MAT::PAR::ViscoGenMax* params)
  :   ElastHyper(params),
      viscoparams_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  MAT::ElastHyper::Pack(data);

  //  pack history data
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    histsize = histstressisoprinclast_->size();
  }
  AddtoPack(data,histsize);  // Length of history vector(s)
  for (int var = 0; var < histsize; ++var)
  {
    AddtoPack(data,histstressisoprinclast_->at(var));
    AddtoPack(data,artstressisoprinclast_->at(var));
    AddtoPack(data,histstressisomodisolast_->at(var));
    AddtoPack(data,artstressisomodisolast_->at(var));
    AddtoPack(data,histstressisomodvollast_->at(var));
    AddtoPack(data,artstressisomodvollast_->at(var));
    AddtoPack(data,histstressanisoprinclast_->at(var));
    AddtoPack(data,artstressanisoprinclast_->at(var));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  MAT::ElastHyper::Unpack(basedata);

  // history data
  int histsize;
  ExtractfromPack(position,data,histsize);

  if (histsize == 0) isinit_=false;

  histstressisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodisocurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodisolast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisocurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisolast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodvolcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodvollast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvolcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvollast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressanisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressanisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  for (int var=0; var<histsize; var+=1)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> tmp(true);
    histstressisoprinccurr_->push_back(tmp);
    artstressisoprinccurr_->push_back(tmp);
    histstressisomodisocurr_->push_back(tmp);
    artstressisomodisocurr_->push_back(tmp);
    histstressisomodvolcurr_->push_back(tmp);
    artstressisomodvolcurr_->push_back(tmp);
    histstressanisoprinccurr_->push_back(tmp);
    artstressanisoprinccurr_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstressisoprinclast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstressisoprinclast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstressisomodisolast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstressisomodisolast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstressisomodvollast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstressisomodvollast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    histstressanisoprinclast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    artstressanisoprinclast_->push_back(tmp);

  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Setup(const int numgp,DRT::INPUT::LineDefinition* linedef)
{
  MAT::ElastHyper::Setup(linedef);

  // flags for the viscous contribution
  vis_isoprinc_ = false ;
  vis_isomod_vol_ = false ;
  vis_isomod_iso_ = false ;
  vis_anisoprinc_ = false ;
  vis_anisomod_vol_ = false ;
  vis_anisomod_iso_ = false ;

  if (viscoparams_->relax_isot_princ_!=0)
    {
      vis_isoprinc_ = true ;
      tau_isoprinc_  = viscoparams_->relax_isot_princ_;
      beta_isoprinc_ = viscoparams_->beta_isot_princ_;
    }
  if (viscoparams_->relax_isot_mod_vol_!=0)
    {
      vis_isomod_vol_ = true ;
      tau_isomod_vol_  = viscoparams_->relax_isot_mod_vol_;
      beta_isomod_vol_ = viscoparams_->beta_isot_mod_vol_;
    }
  if (viscoparams_->relax_isot_mod_isoc_!=0)
    {
      vis_isomod_iso_ = true ;
      tau_isomod_iso_  = viscoparams_->relax_isot_mod_isoc_;
      beta_isomod_iso_ = viscoparams_->beta_isot_mod_isoc_;
    }
  if (viscoparams_->relax_anisot_princ_!=0)
    {
      vis_anisoprinc_ = true ;
      tau_anisoprinc_  = viscoparams_->relax_anisot_princ_;
      beta_anisoprinc_ = viscoparams_->beta_anisot_princ_;
    }
  if (viscoparams_->relax_anisot_mod_vol_!=0)
    {
      vis_anisomod_vol_ = true ;
      tau_anisomod_vol_  = viscoparams_->relax_anisot_mod_vol_;
      beta_anisomod_vol_ = viscoparams_->beta_anisot_mod_vol_;
    }
  if (viscoparams_->relax_anisot_mod_isoc_!=0)
    {
      vis_anisomod_iso_ = true ;
      tau_anisomod_iso_  = viscoparams_->relax_anisot_mod_isoc_;
      beta_anisomod_iso_ = viscoparams_->beta_anisot_mod_isoc_;
    }

  // Initialise/allocate internal stress variables

  histstressisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodisocurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodisolast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisocurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisolast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressisomodvolcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodvollast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvolcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvollast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  histstressanisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressanisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinclast_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  histstressisoprinccurr_->resize(numgp);
  histstressisoprinclast_->resize(numgp);
  artstressisoprinccurr_->resize(numgp);
  artstressisoprinclast_->resize(numgp);

  histstressisomodisocurr_->resize(numgp);
  histstressisomodisolast_->resize(numgp);
  artstressisomodisocurr_->resize(numgp);
  artstressisomodisolast_->resize(numgp);

  histstressisomodvolcurr_->resize(numgp);
  histstressisomodvollast_->resize(numgp);
  artstressisomodvolcurr_->resize(numgp);
  artstressisomodvollast_->resize(numgp);

  histstressanisoprinccurr_->resize(numgp);
  histstressanisoprinclast_->resize(numgp);
  artstressanisoprinccurr_->resize(numgp);
  artstressanisoprinclast_->resize(numgp);

  for (int j=0; j<numgp; ++j)
  {
    histstressisoprinccurr_->at(j) = emptyvec;
    histstressisoprinclast_->at(j) = emptyvec;
    artstressisoprinccurr_->at(j) = emptyvec;
    artstressisoprinclast_->at(j) = emptyvec;
    histstressisomodisocurr_->at(j) = emptyvec;
    histstressisomodisolast_->at(j) = emptyvec;
    artstressisomodisocurr_->at(j) = emptyvec;
    artstressisomodisolast_->at(j) = emptyvec;
    histstressisomodvolcurr_->at(j) = emptyvec;
    histstressisomodvollast_->at(j) = emptyvec;
    artstressisomodvolcurr_->at(j) = emptyvec;
    artstressisomodvollast_->at(j) = emptyvec;
    histstressanisoprinccurr_->at(j) = emptyvec;
    histstressanisoprinclast_->at(j) = emptyvec;
    artstressanisoprinccurr_->at(j) = emptyvec;
    artstressanisoprinclast_->at(j) = emptyvec;
  }
  isinit_ = true;
  return ;
}

/*----------------------------------------------------------------------*
 |  Update internal stress variables              (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Update()
{
  histstressisoprinclast_=histstressisoprinccurr_;
  artstressisoprinclast_=artstressisoprinccurr_;
  histstressisomodisolast_=histstressisomodisocurr_;
  artstressisomodisolast_=artstressisomodisocurr_;
  histstressisomodvollast_=histstressisomodvolcurr_;
  artstressisomodvollast_=artstressisomodvolcurr_;
  histstressanisoprinclast_=histstressanisoprinccurr_;
  artstressanisoprinclast_=artstressanisoprinccurr_;
  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  histstressisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodisocurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodisocurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressisomodvolcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressisomodvolcurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  histstressanisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  artstressanisoprinccurr_=Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  const int numgp=histstressisoprinclast_->size();
  histstressisoprinccurr_->resize(numgp);
  artstressisoprinccurr_->resize(numgp);
  histstressisomodisocurr_->resize(numgp);
  artstressisomodisocurr_->resize(numgp);
  histstressisomodvolcurr_->resize(numgp);
  artstressisomodvolcurr_->resize(numgp);
  histstressanisoprinccurr_->resize(numgp);
  artstressanisoprinccurr_->resize(numgp);
  for (int j=0; j<numgp; ++j)
  {
    histstressisoprinccurr_->at(j) = emptyvec;
    artstressisoprinccurr_->at(j) = emptyvec;
    histstressisomodisocurr_->at(j) = emptyvec;
    artstressisomodisocurr_->at(j) = emptyvec;
    histstressisomodvolcurr_->at(j) = emptyvec;
    artstressisomodvolcurr_->at(j) = emptyvec;
    histstressanisoprinccurr_->at(j) = emptyvec;
    artstressanisoprinccurr_->at(j) = emptyvec;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Reset internal stress variables               (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Reset()
{
  // do nothing,
  // because #histstresscurr_ and #artstresscurr_ are recomputed anyway at every iteration
  // based upon #histstresslast_ and #artstresslast_ untouched within time step
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ViscoGenMax::Evaluate(
  const LINALG::Matrix<6,1>& glstrain,
  LINALG::Matrix<6,6>& cmat,
  LINALG::Matrix<6,1>& stress,
  const int gp,
  Teuchos::ParameterList& params
  )
{
  LINALG::Matrix<6,1> id2(true) ;
  LINALG::Matrix<6,1> rcg(true) ;
  LINALG::Matrix<6,1> scg(true) ;
  LINALG::Matrix<6,1> icg(true) ;
  LINALG::Matrix<6,6> id4(true) ;
  LINALG::Matrix<6,6> id4sharp(true) ;

  LINALG::Matrix<3,1> prinv(true);
  LINALG::Matrix<3,1> modinv(true);
  LINALG::Matrix<6,1> pranisoinv(true);

  LINALG::Matrix<3,1> gamma(true);
  LINALG::Matrix<8,1> delta(true);
  LINALG::Matrix<3,1> modgamma(true);
  LINALG::Matrix<5,1> moddelta(true);
  LINALG::Matrix<3,1> anisogamma(true);
  LINALG::Matrix<15,1> anisodelta(true);

  EvaluateKinQuant(glstrain,id2,scg,rcg,icg,id4,id4sharp,prinv,modinv,pranisoinv);
  EvaluateGammaDelta(rcg,prinv,modinv,pranisoinv,gamma,delta,modgamma,moddelta);

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress.Clear();
  cmat.Clear();

  if (isoprinc_)
    {
    LINALG::Matrix<NUM_STRESS_3D,1> stressisoprinc(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisoprinc(true) ;
    EvaluateIsotropicPrinc(stressisoprinc,cmatisoprinc,scg,id2,icg,id4sharp,gamma,delta);
    stress.Update(1.0, stressisoprinc, 1.0);
    cmat.Update(1.0,cmatisoprinc,1.0);

    // viscous contribution
    if (vis_isoprinc_)
      {
        //initialize scalars
        double artscalar1(true);
        double artscalar2(true);
        double scalarvisco(true);

        const double theta= 0.5 ;//params_->theta_;

        // get time algorithmic parameters
        // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
        // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
        double dt = params.get<double>("delta time"); // TIMESTEP in the .dat file

        // evaluate scalars to compute
        // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
        artscalar1=(tau_isoprinc_ - dt + theta*dt)/tau_isoprinc_;
        artscalar2=tau_isoprinc_/(tau_isoprinc_ + theta*dt);

        // factor to calculate visco stiffness matrix from elastic stiffness matrix
        scalarvisco = 1+beta_isoprinc_*exp(-dt/(2*tau_isoprinc_));//+alpha1*tau/(tau+theta*dt);

        // viscous part of the stress
        // read history
        LINALG::Matrix<NUM_STRESS_3D,1> S_n (histstressisoprinclast_->at(gp));
        S_n.Scale(-beta_isoprinc_);
        LINALG::Matrix<NUM_STRESS_3D,1> Q_n (artstressisoprinclast_->at(gp));

        // artificial visco stresses
        LINALG::Matrix<NUM_STRESS_3D,1> Q(Q_n);
        Q.Scale(artscalar1);
        stressisoprinc.Scale(beta_isoprinc_);
        Q += stressisoprinc;
        Q += S_n;
        Q.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + beta*(S^(n+1) - S^n)]

        // update history
        histstressisoprinccurr_->at(gp) = stressisoprinc;
        artstressisoprinccurr_->at(gp) = Q;

        // add the viscous to the elastic contribution
        stress.Update(1.0,Q,1.0);

        // constitutive tensor
        // viscous contribution : Cmat = Cmatx(1+delta) with (1+delta)=scalarvisco
        cmat.Scale(scalarvisco);
      }
    }


  if (isomod_)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> stressisomodiso(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodiso(true) ;
    LINALG::Matrix<NUM_STRESS_3D,1> stressisomodvol(true) ;
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatisomodvol(true) ;
    EvaluateIsotropicMod(stressisomodiso,stressisomodvol,cmatisomodiso,cmatisomodvol,rcg,id2,icg,id4,id4sharp,modinv,prinv,modgamma,moddelta);
    stress.Update(1.0, stressisomodiso, 1.0);
    stress.Update(1.0, stressisomodvol, 1.0);
    cmat.Update(1.0,cmatisomodiso,1.0);
    cmat.Update(1.0,cmatisomodvol,1.0);

    // viscous contribution

    // viscous part of the isochoric contribution
    //-----------------------------------------------------------------------
    if (vis_isomod_iso_)
      {
        //initialize scalars
        double artscalar1(true);
        double artscalar2(true);
        //double scalarvisco(true);

        const double theta= 0.5 ;//params_->theta_;

        // get time algorithmic parameters
        // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
        // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
        double dt = params.get<double>("delta time"); // TIMESTEP in the .dat file

        // evaluate scalars to compute
        // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
        artscalar1=(tau_isomod_iso_ - dt + theta*dt)/tau_isomod_iso_;
        artscalar2=tau_isomod_iso_/(tau_isomod_iso_ + theta*dt);

        // factor to calculate visco stiffness matrix from elastic stiffness matrix
        //scalarvisco = 1+beta_isomod_iso_*exp(-dt/(2*tau_isomod_iso_));//+alpha1*tau/(tau+theta*dt);

        // read history
        LINALG::Matrix<NUM_STRESS_3D,1> Siso_n (histstressisomodisolast_->at(gp));
        Siso_n.Scale(-beta_isomod_iso_);
        LINALG::Matrix<NUM_STRESS_3D,1> Qiso_n (artstressisomodisolast_->at(gp));

        // artificial visco stresses
        LINALG::Matrix<NUM_STRESS_3D,1> Qiso(Qiso_n);
        Qiso.Scale(artscalar1);
        stressisomodiso.Scale(beta_isomod_iso_);
        Qiso += stressisomodiso;
        Qiso += Siso_n;
        Qiso.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + beta*(S^(n+1) - S^n)]

        // update history
        histstressisomodisocurr_->at(gp) = stressisomodiso;
        artstressisomodisocurr_->at(gp) = Qiso;

        // add the viscous to the elastic contribution
        stress.Update(1.0,Qiso,1.0); //--> at this point, stress includes the isochoric viscoelastic contributions
        //----------------------------------------------------------------------
      }

    // viscous part of the volumetric contribution
    //-----------------------------------------------------------------------
    if (vis_isomod_vol_)
      {
        //initialize scalars
        double artscalar1(true);
        double artscalar2(true);
        double scalarvisco(true);

        const double theta= 0.5 ;//params_->theta_;

        // get time algorithmic parameters
        // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
        // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
        double dt = params.get<double>("delta time"); // TIMESTEP in the .dat file

        // evaluate scalars to compute
        // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
        artscalar1=(tau_isomod_vol_ - dt + theta*dt)/tau_isomod_vol_;
        artscalar2=tau_isomod_vol_/(tau_isomod_vol_ + theta*dt);

        // factor to calculate visco stiffness matrix from elastic stiffness matrix
        scalarvisco = 1+beta_isomod_vol_*exp(-dt/(2*tau_isomod_vol_));//+alpha1*tau/(tau+theta*dt);
        // read history
        LINALG::Matrix<NUM_STRESS_3D,1> Svol_n (histstressisomodvollast_->at(gp));
        Svol_n.Scale(-beta_isomod_vol_);
        LINALG::Matrix<NUM_STRESS_3D,1> Qvol_n (artstressisomodvollast_->at(gp));

        // artificial visco stresses
        LINALG::Matrix<NUM_STRESS_3D,1> Qvol(Qvol_n);
        Qvol.Scale(artscalar1);
        stressisomodvol.Scale(beta_isomod_vol_);
        Qvol += stressisomodvol;
        Qvol += Svol_n;
        Qvol.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + S^(n+1) - S^n]

        // update history
        histstressisomodvolcurr_->at(gp) = stressisomodvol;
        artstressisomodvolcurr_->at(gp) = Qvol;

        // add the viscous contribution
        stress.Update(1.0,Qvol,1.0);//--> at this point, stress includes the isochoric and volumetric viscoelastic contributions
        //----------------------------------------------------------------------

        // constitutive tensor

        // Cmat = Cmat_iso_inf(1+delta_iso) + Cmat_vol_inf(1+delta_vol)
        // 1+delta = scalarvisco ; could be two different scalars for iso and vol contributions
        // overall contributions by taking into account the viscous contribution
        cmat.Scale(scalarvisco);
      }
  }

  /*----------------------------------------------------------------------*/
  // coefficients in principal stretches
  const bool havecoeffstrpr = HaveCoefficientsStretchesPrincipal();
  const bool havecoeffstrmod = HaveCoefficientsStretchesModified();
  if (havecoeffstrpr or havecoeffstrmod) {
    ResponseStretches(cmat,stress,rcg,havecoeffstrpr,havecoeffstrmod);
  }

  /*----------------------------------------------------------------------*/
  //Do all the anisotropic stuff!
  if (anisoprinc_)
  {
      LINALG::Matrix<NUM_STRESS_3D,1> stressanisoprinc(true) ;
      LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> cmatanisoprinc(true) ;
      EvaluateAnisotropicPrinc(stressanisoprinc,cmatanisoprinc,rcg,params);
      stress.Update(1.0, stressanisoprinc, 1.0);
      cmat.Update(1.0, cmatanisoprinc, 1.0);

      // viscous contribution
      if (vis_anisoprinc_)
        {
          //initialize scalars
          double artscalar1(true);
          double artscalar2(true);
          double scalarvisco(true);

          const double theta= 0.5 ;//params_->theta_;

          // get time algorithmic parameters
          // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
          // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
          double dt = params.get<double>("delta time"); // TIMESTEP in the .dat file

          // evaluate scalars to compute
          // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
          artscalar1=(tau_anisoprinc_ - dt + theta*dt)/tau_anisoprinc_;
          artscalar2=tau_anisoprinc_/(tau_anisoprinc_ + theta*dt);

          // factor to calculate visco stiffness matrix from elastic stiffness matrix
          scalarvisco = 1+beta_anisoprinc_*exp(-dt/(2*tau_anisoprinc_));//+alpha1*tau/(tau+theta*dt);

          // viscous part of the stress
          // read history
          LINALG::Matrix<NUM_STRESS_3D,1> S_n (histstressanisoprinclast_->at(gp));
          S_n.Scale(-beta_anisoprinc_);
          LINALG::Matrix<NUM_STRESS_3D,1> Q_n (artstressanisoprinclast_->at(gp));

          // artificial visco stresses
          LINALG::Matrix<NUM_STRESS_3D,1> Q(Q_n);
          Q.Scale(artscalar1);
          stressanisoprinc.Scale(beta_anisoprinc_);
          Q += stressanisoprinc;
          Q += S_n;
          Q.Scale(artscalar2);  // Q^(n+1) = artscalar2* [artscalar1* Q + S^(n+1) - S^n]

          // update history
          histstressanisoprinccurr_->at(gp) = stressanisoprinc;
          artstressanisoprinccurr_->at(gp) = Q;

          // add the viscous to the elastic contribution
          stress.Update(1.0,Q,1.0);

          // constitutive tensor
          // viscous contribution : Cmat = Cmatx(1+delta) with (1+delta)=scalarvisco
          cmat.Scale(scalarvisco);
        }
  }

}

