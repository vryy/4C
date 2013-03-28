/*----------------------------------------------------------------------*/
/*!
\file thermoplastichyperelast.cpp
\brief Contains the functions to establish local material law /
       stress-strain law for isotropic material following finite strain
       von-Mises plasticity with linear isotropic hardening and general
       hyperelasticity (for the time being: NeoHooke).

       implementation is based on
       Simo and Hughes: "Computational Inelasticity."
       Springer, 2000.

       geometrically nonlinear, for finite strains, rate-independent

       example input line:
       MAT 1 MAT_Struct_ThrPlasticHyperElast YOUNG 206.9 NUE 0.29 DENS 0.0
         YIELD 0.45 ISOHARD 0.0 KINHARD 0.0

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "thermoplastichyperelast.H"
#include "matpar_bundle.H"
#include "material_service.H"

#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 03/13 |
 *----------------------------------------------------------------------*/
MAT::PAR::ThermoPlasticHyperElast::ThermoPlasticHyperElast(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  yield_(matdata->GetDouble("YIELD")),
  isohard_(matdata->GetDouble("ISOHARD")),
  kinhard_(matdata->GetDouble("KINHARD"))
{
}


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 03/13 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ThermoPlasticHyperElast::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ThermoPlasticHyperElast(this));
}

MAT::ThermoPlasticHyperElastType MAT::ThermoPlasticHyperElastType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 03/13 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::ThermoPlasticHyperElastType::Create(const std::vector<char>& data )
{
  MAT::ThermoPlasticHyperElast* thrplhyper = new MAT::ThermoPlasticHyperElast();
  thrplhyper->Unpack(data);
  return thrplhyper;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 03/13 |
 *----------------------------------------------------------------------*/
MAT::ThermoPlasticHyperElast::ThermoPlasticHyperElast()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                                 dano 03/13 |
 *----------------------------------------------------------------------*/
MAT::ThermoPlasticHyperElast::ThermoPlasticHyperElast(MAT::PAR::ThermoPlasticHyperElast* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                             dano 03/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  // pack history data
  int histsize;
  // if material is not initialized, i.e. start simulation, nothing to pack
  if (!Initialized())
  {
    histsize=0;
  }
  else
  {
    // if material is initialized (restart): size equates number of gausspoints
    histsize = defgrdlast_->size();
  }

  AddtoPack(data,histsize); // Length of history vector(s)
  for (int var=0; var<histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data,defgrdlast_->at(var));
    AddtoPack(data,bebarlast_->at(var));
    AddtoPack(data,accplstrainlast_->at(var));
  }
  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                           dano 03/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ThermoPlasticHyperElast*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position,data,histsize);

  // if system is not yet initialized, the history vectors have to be intialized
  if (histsize == 0)
    isinit_ = false;

  defgrdlast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  defgrdcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );

  bebarlast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );
  bebarcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> > );

  accplstrainlast_ = Teuchos::rcp( new std::vector<double> );
  accplstraincurr_ = Teuchos::rcp( new std::vector<double> );

  for (int var=0; var<histsize; ++var)
  {
    LINALG::Matrix<3,3> tmp(true);
    ExtractfromPack(position,data,tmp);
    defgrdlast_->push_back(tmp);
    
    ExtractfromPack(position,data,tmp);
    bebarlast_->push_back(tmp);
    
    double tmp1=0.0;
    ExtractfromPack(position,data,tmp1);
    accplstrainlast_->push_back(tmp1);

    // current vectors have to be initialized
    defgrdcurr_->push_back(tmp);
    bebarcurr_->push_back(tmp);
    accplstraincurr_->push_back(tmp1);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

}  // Unpack


/*---------------------------------------------------------------------*
 | initialise / allocate internal variables (public)                   |
 *---------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // initialise hist variables
  defgrdlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  defgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);

  bebarlast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  bebarcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);

  accplstrainlast_ = Teuchos::rcp(new std::vector<double>);
  accplstraincurr_ = Teuchos::rcp(new std::vector<double>);

  defgrdlast_->resize(numgp);
  defgrdcurr_->resize(numgp);

  bebarlast_->resize(numgp);
  bebarcurr_->resize(numgp);

  accplstrainlast_->resize(numgp);
  accplstraincurr_->resize(numgp);

  LINALG::Matrix<3,3> emptymat(true);
  for (int i=0; i<3; i++)
    emptymat(i,i) = 1.0;

  for (int i=0; i<numgp; i++)
  {
    defgrdlast_->at(i) = emptymat;
    defgrdcurr_->at(i) = emptymat;

    bebarlast_->at(i) = emptymat;
    bebarcurr_->at(i) = emptymat;

    accplstrainlast_->at(i) = 0.0;
    accplstraincurr_->at(i) = 0.0;
  }

  isinit_ = true;
  return;

}  // Setup()


/*----------------------------------------------------------------------*
 | update internal variables                                 dano 03/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Update()
{
  // make current values at time step t_n+1 to values of last step t_n
  defgrdlast_ = defgrdcurr_;
  bebarlast_ = bebarcurr_;
  accplstrainlast_ = accplstraincurr_;

  // empty vectors of current data
  defgrdcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  bebarcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3,3> >);
  accplstraincurr_ = Teuchos::rcp(new std::vector<double>);

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = defgrdlast_->size();
  defgrdcurr_->resize(histsize);
  bebarcurr_ ->resize(histsize);
  accplstraincurr_ ->resize(histsize);

  LINALG::Matrix<3,3> emptymat(true);
  for (int i=0; i<histsize; i++)
  {
    defgrdcurr_->at(i) = emptymat;
    bebarcurr_ ->at(i) = emptymat;
    accplstraincurr_ ->at(i) = 0.0;
  }

  return;
}  // Update()


/*----------------------------------------------------------------------*
 | calculate stress and constitutive tensor                  dano 03/13 |
 *----------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::Evaluate(
  const LINALG::Matrix<3,3>* defgrd,
  const LINALG::Matrix<6,1>* glstrain,
  Teuchos::ParameterList& params,
  LINALG::Matrix<6,1>* stress,
  LINALG::Matrix<6,6>* cmat
  )
{
  // extract the gauss points from the parameter list
  const int gp = params.get<int>("gp",-1);
  if (gp == -1) dserror("no Gauss point number provided in material");

  // elastic material data
  // get material parameters
  const double ym = params_->youngs_;  // Young's modulus
  const double nu = params_->poissonratio_;  // Poisson's ratio
  const double G = ym / (2.0 * (1.0 + nu));  // shear modulus, mu=G
  const double kappa = ym / (3.0 * (1.0 - 2.0 * nu));  // bulk modulus
  double Hiso = params_->isohard_;  // linear isotropic hardening modulus
  double sigma_y0 = params_->yield_;  // initial yield stress

  // no kinematic hardening yet
  if (params_->kinhard_!=0.) 
    dserror("Sorry, but kinematic hardening not yet implemented.");

  // 3x3 2nd-order identity matrix
  LINALG::Matrix<3,3> id2(true);
  for (int i=0; i<3; i++) id2(i,i) = 1.0;

  // start with current deformation
  defgrdcurr_->at(gp) = *defgrd;
  // calculate the Jacobi-determinant J = det(F_{n+1})
  double J = defgrd->Determinant();

  // ------------------------------------------ multiplicative split of F
  // split deformation gradient in elastic and plastic part
  // F_{n+1} = F_{n+1}^e . F_{n+1}^p

  // relative deformation gradient
  // f_{n+1} = F_{n+1} . (F_n)^-1
  LINALG::Matrix<3,3> defgrddelta(false);
  LINALG::Matrix<3,3> invdefgrdlast(defgrdlast_->at(gp));
  invdefgrdlast.Invert();
  defgrddelta.Multiply(*defgrd,invdefgrdlast);

  // isochoric part of relative deformation gradient
  // fbar_{n+1} = Fbar_{n+1} . Fbar_n^{-1} = (J_{n+1}/J_n)^{1/3} . f_{n+1}
  // with J_{n+1}/J_n = det(fbar_)
  LINALG::Matrix<3,3> defgrddeltabar(defgrddelta);
  defgrddeltabar.Scale(pow(defgrddelta.Determinant(),-1.0/3.0));

  // -------------------------------------------------- elastic predictor

  // elastic left Cauchy-Green (LCG) trial state (isochoric) (9.3.13)
  // bbar_{n+1}^{e,trial} = Fbar_{n+1} (Cbar_{n+1}^{p-1})^{trial} . Fbar_{n+1}^T
  //                      = fbar_{n+1} (bbar_{n}^{e,trial} . fbar_{n+1}^T
  LINALG::Matrix<3,3> bebar_trial(false);
  LINALG::Matrix<3,3> tmp;
  tmp.Multiply(defgrddeltabar,bebarlast_->at(gp));
  bebar_trial.MultiplyNT(tmp,defgrddeltabar);
  // trace of strain vector
  double tracebebar = ( bebar_trial(0,0) + bebar_trial(1,1) + bebar_trial(2,2) );

  // trial Kirchhoff stress deviator (9.3.9)
  // s_{n+1)^{trial} = G . dev(bbar_{n+1}^{e,trial}
  LINALG::Matrix<3,3> devtau_trial(bebar_trial);
  for (int i=0; i<3; i++)
    devtau_trial(i,i) -= 1.0/3.0 * tracebebar;
  devtau_trial.Scale(G);

  // trial equivalent von Mises stress
  // q^{trial} = sqrt(s^{trial}_ij * s^{trial}_ij)
  double q_trial = 0.0;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      q_trial += devtau_trial(i,j) * devtau_trial(i,j);
  q_trial = sqrt(q_trial);

  // calculate yield function at trial state
  // Phi = || s_{n+1}^{trial} || - sqrt(2/3) . sigma_y
  double Phi_trial = 0.0;
  Phi_trial = q_trial - sqrt(2.0/3.0) * (sigma_y0 + Hiso * (accplstrainlast_->at(gp)));

  // stress variables tau = J_{n+1} . p_{n+1} . I + s_{n+1}
  LINALG::Matrix<3,3> devtau;

  // some computations
  // mubar = 1/3 mu tr(bebar_{n+1}^{e,trial})
  double mubar = G * 1.0/3.0 * tracebebar;

  // initialise incremental plastic multiplier Delta gamma
  double Dgamma = 0.0;

  // unit spatial flow vector
  LINALG::Matrix<3,3> n(devtau_trial);
  if (q_trial != 0.0) n.Scale(1.0/q_trial);

  //-------------------------------------------------------------------
  // IF: elastic step (Phi_trial <= 0.0, Dgamma = 0.0)
  //-------------------------------------------------------------------
  if (Phi_trial <= 0.0)
  {
    // trial state vectors = result vectors of time step n+1

    // _n^{trial} --> _{n+1}
    bebarcurr_->at(gp) = bebar_trial;
    accplstraincurr_->at(gp) = accplstrainlast_->at(gp);
    devtau = devtau_trial;

  } // end if (Phi_trial <= 0.0), i.e. elastic step

  //-------------------------------------------------------------------
  // ELSE consistency condition is violated, i.e. plastic load step
  // ( Phi_trial > 0.0, Dgamma >= 0.0 )
  //-------------------------------------------------------------------
  else
  {
    // actual return map
    Dgamma = (Phi_trial / (2.0 * mubar)) / ( 1 + (IsoHard() / (3.0 * mubar)) );

    // update stress
    // s_{n+1} = s_{n+1}^{trial} - 2 . G . Delta gamma . n
    devtau.Update(1.0, devtau_trial, 0.0);
    devtau.Update(-2.0 * mubar * Dgamma, n, 1.0);

    // --------------------------------------------------- update history

    // update accumulated plastic strain
    accplstraincurr_->at(gp) = accplstrainlast_->at(gp) + (sqrt(2.0/3.0) * Dgamma);

    // update elastic LCG
    // bbar_{n+1}^e = bbar_{n+1}^{e,trial} - 2/3 . Dgamma . tr(bbar_{n+1}^{e,trial}) . n
    bebarcurr_->at(gp) = (bebar_trial);
    bebarcurr_->at(gp).Update((-2.0/3.0 * Dgamma * tracebebar), n, 1.0);
  }

  // addition of mean stress to gain Kirchhoff stress tau (9.2.6)
  // tau = J^e . p . I + devtau
  // with p := U'(J^e) = kappa/2 ( (J^e)^2 -1 ) / J^e
  // --> tau = kappa/2 ( (J^e)^2 -1 ) . I + devtau
  // different to Miehe (2.37): p = kappa (J^2 - 1) / J
  double p = kappa/2.0 * (J * J - 1.0) / J;
  LINALG::Matrix<3,3> tau(devtau);
  for (int i=0; i<3; i++)
    tau(i,i) += J * p;

  // transform Kirchhoff stress to 2.PK-stress
  // PK2 = F^{-1} . tau . F^{-T}
  LINALG::Matrix<3,3> PK2;
  LINALG::Matrix<3,3> invdefgrdcurr(*defgrd);
  invdefgrdcurr.Invert();
  tmp.Multiply(invdefgrdcurr,tau);
  PK2.MultiplyNT(tmp,invdefgrdcurr);

  // output 2.PK-stress in Voigt-notation
  // stresses are correct compared with other
  // constitutive laws in domain of linear elasticity
  (*stress)(0) = PK2(0,0);
  (*stress)(1) = PK2(1,1);
  (*stress)(2) = PK2(2,2);
  (*stress)(3) = 0.5 * (PK2(0,1) + PK2(1,0));
  (*stress)(4) = 0.5 * (PK2(1,2) + PK2(2,1));
  (*stress)(5) = 0.5 * (PK2(0,2) + PK2(2,0));

  //---------------------------------------------------------------------------
  // ----------------------- consistent elastoplastic tangent modulus (Box 9.2)
  //---------------------------------------------------------------------------

  // C_ep = C_{n+1}^{trial} - beta1 Cbar_{n+1}^{trial} - 2 mubar beta3 n \otimes n
  //        - 2 mubar beta4 sym[ n \otimes dev[n^2] ]^s
  // with C_{n+1}^{trial} =  TODO

  // dev (n^2)
  LINALG::Matrix<3,3> devnsquare(true);
  devnsquare.Multiply(n,n);
  double tracensquare= (devnsquare(0,0) + devnsquare(1,1) + devnsquare(2,2));
  for (int i=0; i<3; i++)
    devnsquare(i,i) -= 1.0/3.0 * tracensquare;

  // scaling factors for spatial tangent
  double beta0, beta1, beta2, beta3, beta4;
  beta0 = 1 + Hiso / 3.0 / mubar;
  beta1 = 2.0 * mubar * Dgamma / q_trial;
  beta2 = (1.0 - 1.0 / beta0) * 2.0/3.0 * q_trial / mubar * Dgamma;
  beta3 = 1.0 / beta0 - beta1 + beta2;
  beta4 = (1.0 / beta0 - beta1) * q_trial / mubar;

  // calculate material deformation gradient right away
  LINALG::Matrix<6,6> Cmat(true);

  LINALG::Matrix<6,6> Cbar_trialMaterial(true);

  // calculate the right Cauchy Green (RCG) deformation tensor and its inverse
  LINALG::Matrix<3,3> RCG;
  RCG.MultiplyTN(*defgrd,*defgrd);
  LINALG::Matrix<3,3> invRCG;
  invRCG.Invert(RCG);

  // pull back of the spatial flow vector n to N
  // N = F^{-1} n F^{-T}
  LINALG::Matrix<3,3> N;
  tmp.Multiply(invdefgrdcurr,n);
  N.MultiplyNT(tmp,invdefgrdcurr);

  // pullback of dev(n^2)
  LINALG::Matrix<3,3> devNsquare;
  tmp.Multiply(invdefgrdcurr,devnsquare);
  devNsquare.MultiplyNT(tmp,invdefgrdcurr);

  // Cbar_trial
  // spatial: cbar_trial = 2 . mubar I_d - 2/3 qbar [n \otimes I + I \otimes n]
  // with I_d = I_s - 1/3 . I \otimes I
  // pullback of I --> invRCG
  ElastSymTensor_o_Multiply(Cbar_trialMaterial, 2.0 * mubar, invRCG, invRCG, 1.0);
  ElastSymTensorMultiply(Cbar_trialMaterial, -2.0/3.0 * mubar, invRCG, invRCG, 1.0);
  ElastSymTensorMultiplyAddSym(Cbar_trialMaterial, -2.0/3.0 * q_trial, N, invRCG, 1.0);

  // C_trial
  // spatial c_trial = (J . U')' . J . I \otimes I - 2 J U' I4
  // with U'(J^e) = kappa/2 ( (J^e)^2 -1 ) / J^e
  // C_trial = kappa/2 (J^e)^2 C^{-1} \otimes C^{-1} - kappa ( (J^e)^2 -1 ) C^{-1} \otimes C^{-1}
  ElastSymTensorMultiply(Cmat, kappa * J * J, invRCG, invRCG, 1.0);
  ElastSymTensor_o_Multiply(Cmat,-1.0 * kappa * (J * J - 1.0), invRCG, invRCG, 1.0);
  Cmat.Update(1.0, Cbar_trialMaterial, 1.0);

  // plastic step update
  if (Dgamma != 0.0)
  {
    // this is nonlinear mechanics
    Cmat.Update((-1.0 * beta1), Cbar_trialMaterial, 1.0);
    ElastSymTensorMultiply(Cmat, (-2.0 * mubar * beta3), N, N, 1.0);
    ElastSymTensorMultiply(Cmat, (-2.0 * mubar * beta4), N, devNsquare, 1.0);
  }

  // update material tangent
  // cmat = C_ep = C_trial + Cbar_trial + C_p
  *cmat = Cmat;

  return;

}  // Evaluate()


/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::ThermoPlasticHyperElast::VisNames(std::map<string,int>& names)
{
  string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1; // scalar

}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool MAT::ThermoPlasticHyperElast::VisData(
  const string& name,
  std::vector<double>& data,
  int numgp
  )
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    double temp = 0.0;
    for (int iter=0; iter<numgp; iter++)
      temp += AccumulatedStrain(iter);
    data[0] = temp/numgp;
  }
  return true;
}  // VisData()


/*----------------------------------------------------------------------*/

