/*----------------------------------------------------------------------*/
/*!
\file damage.cpp
\brief TODO Contains the functions to establish local material law /
       stress-strain law for isotropic material for a 3D hex element
       following perfectly von Mises plasticity and a linear elastic material law
       (St.Venant Kirchhoff).

       perfect plasticity:
        - no hardening allowed
        - independent yield stress level of degree of plastification
        - constant uniaxial yield stress \f$ \sigma_y \,=\, const.\f$

       extend to linear kinematic hardening
        - describing the Bauschinger effect via \f$ Hkin \,=\, const.\f$
        - constant uniaxial yield stress \f$ \sigma_y \,=\, const.\f$

       geometric linear, for small strains

       example input line:
       MAT 1 MAT_Struct_Damage YOUNG 206.9 NUE 0.29 DENS 0.0
         YIELD 0.45 KINHARD 0.0 TOL 1.0e-6

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/

/*----------------------------------------------------------------------*
 | Headers                                                   dano 04/11 |
 *----------------------------------------------------------------------*/
#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "damage.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 04/11 |
 *----------------------------------------------------------------------*/
MAT::PAR::Damage::Damage(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  yield_(matdata->GetDouble("YIELD")),
  kinhard_(matdata->GetDouble("KINHARD")),
  abstol_(matdata->GetDouble("TOL"))
{
}

/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Damage::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Damage(this));
}


MAT::DamageType MAT::DamageType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::DamageType::Create( const std::vector<char> & data )
{
  MAT::Damage* plastic = new MAT::Damage();
  plastic->Unpack(data);
  return plastic;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 04/11 |
 *----------------------------------------------------------------------*/
MAT::Damage::Damage()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                                 dano 04/11 |
 *----------------------------------------------------------------------*/
MAT::Damage::Damage(MAT::PAR::Damage* params)
: plastic_step(false),
  params_(params)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                             dano 04/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  // in case we are in post-process mode
  if (params_ != NULL) matid = params_->Id();
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
    histsize = strainpllast_->size();
  }
  AddtoPack(data,histsize); // Length of history vector(s)
  for (int var=0; var<histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data,strainpllast_->at(var));
    AddtoPack(data,backstresslast_->at(var));

    AddtoPack(data,strainbarpllast_->at(var));
  }

  return;
} // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                           dano 04/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Unpack(const std::vector<char>& data)
{
  isinit_=true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Damage*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int histsize;
  ExtractfromPack(position,data,histsize);

  // if system is not yet initialized, the history vectors have to be intialized
  if (histsize == 0)
    isinit_ = false;

  // unpack strain vectors
  strainpllast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  strainplcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );

  // unpack back stress vectors (for kinematic hardening)
  backstresslast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  backstresscurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> > );

  strainbarpllast_ = Teuchos::rcp( new std::vector<LINALG::Matrix<1,1> > );
  strainbarplcurr_ = Teuchos::rcp( new std::vector<LINALG::Matrix<1,1> > );

  for (int var=0; var<histsize; ++var)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> tmp(true);
    // vectors of last converged state are unpacked
    ExtractfromPack(position,data,tmp);
    strainpllast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    backstresslast_->push_back(tmp);

    LINALG::Matrix<1,1> tmp1(true);
    // scalar-valued vector of last converged state are unpacked
    ExtractfromPack(position,data,tmp1);
    strainbarpllast_->push_back(tmp1);

    // current vectors have to be initialized
    strainplcurr_->push_back(tmp);
    backstresscurr_->push_back(tmp);

    strainbarplcurr_->push_back(tmp1);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

} // Unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal stress variables (public)      04/11 |
 *---------------------------------------------------------------------*/
void MAT::Damage::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // initialize hist variables
  strainpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  backstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  strainbarpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1,1> >);
  strainbarplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1,1> >);

  LINALG::Matrix<NUM_STRESS_3D,1> emptymat(true);
  strainpllast_->resize(numgp);
  strainplcurr_->resize(numgp);

  backstresslast_->resize(numgp);
  backstresscurr_->resize(numgp);

  LINALG::Matrix<1,1> emptymat1(true);
  strainbarpllast_->resize(numgp);
  strainbarplcurr_->resize(numgp);

  for (int i=0; i<numgp; i++)
  {
    strainpllast_->at(i) = emptymat;
    strainplcurr_->at(i) = emptymat;

    backstresslast_->at(i) = emptymat;
    backstresscurr_->at(i) = emptymat;

    strainbarpllast_->at(i) = emptymat1;
    strainbarplcurr_->at(i) = emptymat1;
  }

  isinit_=true;
  return;

}  // Setup()


/*---------------------------------------------------------------------*
 | update internal stress variables (public)                dano 04/11 |
 *---------------------------------------------------------------------*/
void MAT::Damage::Update()
{
  // make current values at time step tlast+1 to values of last step tlast
  strainpllast_ = strainplcurr_;
  backstresslast_ = backstresscurr_;

  strainbarpllast_ = strainbarplcurr_;

  // empty vectors of current data
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  strainbarplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<1,1> >);

  // get the size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int histsize = strainpllast_->size();
  strainplcurr_->resize(histsize);
  backstresscurr_->resize(histsize);

  strainbarplcurr_->resize(histsize);

  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  const LINALG::Matrix<1,1> emptyvec1(true);
  for (int i=0; i<histsize; i++)
  {
    strainplcurr_->at(i) = emptyvec;
    backstresscurr_->at(i) = emptyvec;

    strainbarplcurr_->at(i) = emptyvec1;
  }

  return;
}  // Update()


/*----------------------------------------------------------------------*
 | evaluate material (public)                                dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Evaluate(
  const LINALG::Matrix<3,3>* defgrd,
  const LINALG::Matrix<NUM_STRESS_3D,1>* linstrain,  // linear strain vector
  Teuchos::ParameterList& params,  // parameter list for communication & HISTORY
  LINALG::Matrix<NUM_STRESS_3D,1>* stress, // 2nd PK-stress
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat // material stiffness matrix
  )
{
  const int gp = params.get<int>("gp",-1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  LINALG::Matrix<MAT::NUM_STRESS_3D,1> plstrain(true);

  // get material parameters
  // Young's modulus
  double young = params_->youngs_;
  // Poisson's ratio
  double nu = params_->poissonratio_;
  // initial yield stress
  double sigma_y0 = params_->yield_;
  // linear kinematic hardening modulus
  double Hkin = params_->kinhard_;

  // initialize scalars
  // lame constant
  // shear modulus parameter mu == G
  double G = 0.0;
  G = young / ( 2.0 * (1.0 + nu) );
  // bulk modulus kappa = E /( 3 ( 1 - 2 nu) )= lambda + 2/3 * mu
  double kappa = 0.0;
  kappa = young /( 3 * (1 - 2 * nu) );

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<NUM_STRESS_3D,1> id2(true);
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // linstrain (in): independent variable passed from the element
  //  strain^p: evolution is determined by the flow rule, history varible
  //  strain^e: definition of additive decomposition:
  //  strain^e = strain - strain^p
  // REMARK: stress-like 6-Voigt vector
  LINALG::Matrix<NUM_STRESS_3D,1> strain(*linstrain);

  //---------------------------------------------------------------------------
  // elastic predictor (trial values)
  //---------------------------------------------------------------------------

  // -------------------------------------------- old plastic strains
  // strain^{p,trial}_{n+1} = strain^p_n
  // equivalent plastic strain
  LINALG::Matrix<NUM_STRESS_3D,1> strain_p(true);
  for (int i=0; i<6; i++)
    strain_p(i,0) = strainpllast_->at(gp)(i,0);

  // get old equivalent plastic strain only in case of plastic step
  double strainbar_p;
  // accumulated or equivalent plastic strain (scalar-valued)
  // astrain^{p,trial}_{n+1} = astrain^p_n
  strainbar_p = (strainbarpllast_->at(gp))(0,0);

  // ------------------------------------------------ old back stress
  // beta^{trial}_{n+1} = beta_n
  LINALG::Matrix<NUM_STRESS_3D,1> beta(true);
  for (int i=0; i<6; i++)
    beta(i,0) = backstresslast_->at(gp)(i,0);

  // ----------------------------------------------- physical strains
  // convert engineering shear components into physical components
  // input strain is given in Voigt-notation

  // convert engineering shear component (in) into physical component
  for (int i=3; i<6; ++i) strain(i) /= 2.0;
  for (int i=3; i<6; ++i) strain_p(i) /= 2.0;

  // ----------------------------------------------- elastic trial strain
  // assume load step is elastic
  // strain^{e}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> strain_e(true);

  // strain^{e,trial}_{n+1} = strain_n+1 - strain^p_n
  LINALG::Matrix<NUM_STRESS_3D,1> trialstrain_e(true);
  trialstrain_e.Update( 1.0, strain, 0.0 );
  trialstrain_e.Update( (-1.0), strain_p, 1.0 );

  // volumetric strain
  // trace of strain vector
  double tracestrain = ( trialstrain_e(0)+trialstrain_e(1)+trialstrain_e(2) );
  // volstrain = 1/3 . tr( strain ) . Id
  LINALG::Matrix<NUM_STRESS_3D,1> volumetricstrain(true);
  volumetricstrain.Update( (tracestrain/3.0), id2, 0.0 );

  // deviatoric strain
  // dev^e = strain^e - volstrain^e
  LINALG::Matrix<NUM_STRESS_3D,1> devstrain(true);
  devstrain.Update( 1.0, trialstrain_e, 0.0 );
  devstrain.Update( -1.0, volumetricstrain, 1.0 );

  // ------------------------------------------------------- trial stress
  // pressure = kappa . tr( strain ): saved as scalar
  double p = kappa * tracestrain;

  // deviatoric stress = 2 . G . devstrain
  LINALG::Matrix<NUM_STRESS_3D,1> devstress(true);
  devstress.Update( 2*G, devstrain, 0.0 );
  // be careful for shear stresses (sigma_12)
  // in Voigt-notation the shear strains have to be scaled with 1/2
  // normally done in the material tangent (cf. id4sharp)

  // --------------------------------------  relative effective stress
  // eta^(trial)_{n+1} = s^(trial)_{n+1} - beta^(trial)_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> eta(true);
  RelStress( devstress, beta, eta);

  // --------------------------  elastic trial von Mises effective stress
  // q^(trial)_{n+1} := q(s^(trial)_{n+1}) = \sqrt{ 3 . J2 }
  // J2 = 1/2 (s11^2 + s22^2 + s33^2 + 2 . s12^2 + 2 . s23^2 + 2 . s13^2)
  double J2 = 0.0;
  J2 = 2*G * 2*G * ( 1/2.0 * ( devstrain(0)*devstrain(0)   +
                               devstrain(1)*devstrain(1)   +
                               devstrain(2)*devstrain(2) ) +
                   devstrain(3)*devstrain(3) +
                   devstrain(4)*devstrain(4) +
                   devstrain(5)*devstrain(5)   );
  double q = 0.0;
  q = sqrt( 3.0 * J2 );

  //---------------------------------------------------------------------------
  // check plastic admissibility
  //---------------------------------------------------------------------------

  // ----------------------------------------------- trial yield function

  // calculate the uniaxial yield stress
  // sigma_y = sigma_y0
  double sigma_y = sigma_y0;

  // calculate the yield function
  // Phi = \sqrt{ 3.0 . J2 } - sigma_y = q - sigma_y
  // with trial values: Phi_trial = q_trial - sigma_y
  double Phi_trial = 0.0;
  Phi_trial = q - sigma_y;

  // initialise
  // if trial state is violated, there are 2 possible states:
  double heaviside = 0.0;
  // incremental plastic multiplier Delta gamma
  double Dgamma = 0.0;
  // kinematic hardening curve of current time step and old time step
  // betabar = Hkin * strainbar_p
  // linear kinematic hardening: Hkin = const., else: Hkin = Hkin(strainnbar_p)
  double betabarold = 0.0;
  double betabar = 0.0;

  // flow vector N (Prandtl-Reuss)
  // (using the updated deviatoric stress s_n+1, no longer s_n+1^trial)
  // N = sqrt{3/2} . ( s_{n+1} / || s_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D,1> N(true);

  // flow vector flovec (Prandtl-Reuss)
  // (using the updated relative stress eta_n+1, no longer eta_n+1^trial)
  // bN = sqrt{3/2} . ( eta_{n+1} / || eta_{n+1} || )
  LINALG::Matrix<NUM_STRESS_3D,1> flovec(true);

  //-------------------------------------------------------------------
  // IF elastic step (Phi_trial <= 0.0, Dgamma = 0.0)
  //-------------------------------------------------------------------
  if (Phi_trial < 0.0)
  {
    // trial state vectors = result vectors of time step n+1
    // sigma^e_n+1 = sigma^(e,trial)_n+1 = s^(trial)_{n+1} + p. I
    Stress( p, devstress, *stress );

    // total strains
    // strain^e_{n+1} = strain^(e,trial)_{n+1}
    // compute converged engineering strain components (Voigt-notation)
    strain_e.Update( 1.0, trialstrain_e, 0.0);
    for (int i=3; i<6; ++i) strain_e(i) *= 2.0;

    // no plastic yielding
    Dgamma = 0.0;
    // kinematic hardening curve of current time step and old time step
    // betabar = Hkin * strainbar_p
    // linear kinematic hardening: Hkin = const., else: Hkin = Hkin(strainnbar_p)
    betabarold = 0.0;
    betabar = 0.0;

    // pass the current plastic strains to the element (for visualisation)
    // compute converged engineering strain components (Voigt-notation)
    for (int i=3; i<6; ++i) strain_p(i) *= 2.0;
    plstrain.Update(1.0, strain_p, 0.0);

    // constant values for
    //  - plastic strains
    //  - accumulated plastic strains
    //  - relative stress

  }  // elastic step

  //-------------------------------------------------------------------
  // ELSE IF plastic step ( Phi_trial > 0.0, Dgamma >= 0.0 )
  // violated consistency condition
  //-------------------------------------------------------------------
  else //if (Phi_trial > 0.0)
  {
    // only first plastic call is output at screen for every processor
    // visualisation of whole plastic behaviour via PLASTIC_STRAIN in postprocessing
    if (plastic_step == false)
    {
      if (gp == 0)
        cout << "first plastic step: plastic strains unequal zero!" << endl;
      plastic_step = true;
    }

    // calculate kinematic hardening stress of old time step
    // beta_{n} = Hkin * astrain^p_{n} = Hkin * astrain^{p, trial}_{n+1}
    betabarold = Hkin * strainbar_p;

    // -------------------------------------------------- return-mapping

    // local Newton-Raphson

    // initialize
    const int itermax = 50;  // max. number of iterations
    int itnum = 0;  // iteration counter

    // Res:= residual of Newton iteration == yield function
    double Res = 0.0;
    // calculate residual derivative/tangent
    // ResTan = Phi' = d(Phi)/d(Dgamma)
    double ResTan = 0.0;
    // initialise to zero
    Dgamma = 0.0;

    // start iteration
    while (true)
    {
      itnum++;
      // check for convergence

      // if not converged
      if ( itnum>itermax )
      {
        dserror("local Newton iteration did not converge after iteration %3d/%3d",itnum,itermax);
      }
      // continue loop

      // Res:= Phi = q^(trial)_{n+1} - 3 . G . Delta gamma - sigma_y
      // Delta gamma . Hkin := betabar - betabarn
      // Res = q - 3 * G * Dgamma - sigma_y;
      Res = q - 3 * G * Dgamma - betabar + betabarold - sigma_y;

      // check for convergence
      double norm = abs(Res);
      if ( norm < params_->abstol_)
      {
#ifdef DEBUGMATERIAL
        if (gp == 0)
          printf("Newton method converged after %i iterations; abs(Res)=  %-14.8E\n", itnum, abs(Res));
#endif  // #ifdef DEBUGMATERIAL
        break;
      }

      // plasticity
      // ResTan = -3G = const.
      ResTan = - 3 * G;

      // incremental plastic multiplier Dgamma
      // Dgamma = Dgamma - Phi / Phi'
      Dgamma += ( -Res )/ResTan;

      // -------------------------------------- update of plastic values
      // compute new residual of accumulatd plastic strains
      strainbar_p += Dgamma;

      // sigma = sigma_y0
      sigma_y = sigma_y0;

#ifdef DEBUGMATERIAL
      if (gp == 0)
      {
        cout << "am 1.GP: local Newton: Res " << Res << endl;
        cout << "local Newton: ResTan " << ResTan << endl;
        cout << "local Newton: Dgamma " << Dgamma << endl;
        cout << "local Newton: betabarold " << betabarold << endl;
        cout << "local Newton: betabar " << betabar << "\n"<< endl;
      }
#endif  // #ifdef DEBUGMATERIAL

    }  // end of local Newton iteration

    // ----------------------------------------------------------- update

    flovec.Update( 1.0, devstress, 0.0);
    // relative stress norm || eta_{n+1}^{trial} ||
    double devstressnorm = 0.0;
    devstressnorm = sqrt( devstress(0)*devstress(0) + devstress(1)*devstress(1) + devstress(2)*devstress(2) +
                    2 * ( devstress(3)*devstress(3) + devstress(4)*devstress(4) + devstress(5)*devstress(5) ) );
    const double facflovec = (sqrt( 3.0/2.0 )) /  devstressnorm ;
    // flow vector flovec = sqrt(3/2) devstress_{n+1}^{trial}  / || devstress_{n+1}^{trial} ||
    flovec.Scale(facflovec);

    // deviatoric stress
    // s = s_{n+1}^{trial} - 2 . G . Delta gamma . flovec
    const double facdevstress = (-2.0) * G * Dgamma;
    devstress.Update( facdevstress, flovec, 1.0);

    // total stress
    // sigma_{n+1} = s_{n+1} + p_{n+1} . Id
    // pressure/volumetric stress no influence due to plasticity
    Stress( p, devstress, *stress );

    // total strains
    // strain^e_{n+1} = strain^(e,trial)_{n+1} - Dgamma . flovec
    // compute converged engineering strain components (Voigt-notation)
    strain_e.Update( 1.0, trialstrain_e, 0.0);
    strain_e.Update( (-Dgamma), flovec, 1.0 );

    // strain^p_{n+1} = strain^p_n + Dgamma . flovec
    strain_p.Update( Dgamma, flovec, 1.0 );

    // compute converged engineering strain components (Voigt-notation)
    for (int i=3; i<6; ++i) strain_e(i) *= 2.0;
    for (int i=3; i<6; ++i) strain_p(i) *= 2.0;

    // pass the current plastic strains to the element (for visualisation)
    plstrain.Update(1.0, strain_p, 0.0);

    // --------------------------------------------------- update history
    // plastic strain
    strainplcurr_->at(gp) = strain_p;

    // accumulated plastic strain
    strainbarplcurr_->at(gp) = strainbar_p;

    // back stress
    backstresscurr_->at(gp) = beta;

#ifdef DEBUGMATERIAL
    cout << "end strain_p\n " << strain_p << endl;
    cout << "end strainplcurr_->at(gp)\n " << strainplcurr_->at(gp) << endl;
#endif //ifdef DEBUGMATERIAL

  }  // plastic corrector

  // if Phi^trial=0: two tangent stress-strain relations exist
  // plastic loading --> C == C_ep
  if (Dgamma > 0.0)
    heaviside = 1.0;
  // elastic unloading --> C == C_e
  else heaviside = 0.0;

  //---------------------------------------------------------------------------
  // --------------------------------- consistent elastoplastic tangent modulus
  //---------------------------------------------------------------------------

  // using an associative flow rule: C_ep is symmetric
  // ( generally C_ep is nonsymmetric )
  SetupCmatElastoPlastic(
    *cmat,
    Dgamma,
    G,
    q,
    flovec,
    devstress,
    heaviside
    );

#ifdef DEBUGMATERIAL
  cout << "Nach Setup Cep\n" << endl;
  cout << " Dgamma " << Dgamma << endl;
  cout << " G " << G << endl;
  cout << " qbar " << qbar << endl;
  cout << " flow vector " << flovec << endl;
  cout << " heaviside " << heaviside << endl;

  cout << "cmat " << cmat << endl;
#endif // #ifdef DEBUGMATERIAL

  params.set<LINALG::Matrix<MAT::NUM_STRESS_3D,1> >("plglstrain",plstrain);

  return;

}  // Evaluate()


/*----------------------------------------------------------------------*
 | computes linear stress tensor                             dano 05/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::Stress(
  const double p,  // volumetric stress
  const LINALG::Matrix<NUM_STRESS_3D,1>& devstress,  // deviatoric stress tensor
  LINALG::Matrix<NUM_STRESS_3D,1>& stress // 2nd PK-stress
  )
{
  // total stress = deviatoric + hydrostatic pressure . I
  // sigma = s + p . I
  stress.Update(1.0, devstress, 0.0);
  for (int i=0; i<3; ++i) stress(i) += p;

}  // Stress()


/*----------------------------------------------------------------------*
 | compute relative deviatoric stress tensor                 dano 08/11 |
 *----------------------------------------------------------------------*/
void MAT::Damage::RelStress(
  const LINALG::Matrix<NUM_STRESS_3D,1>& devstress,  // deviatoric stress tensor
  const LINALG::Matrix<NUM_STRESS_3D,1>& beta,  // back stress tensor
  LINALG::Matrix<NUM_STRESS_3D,1>& eta // relative stress
  )
{
  // relative stress = deviatoric - back stress
  // eta = s - beta
  eta.Update( 1.0, devstress, 0.0 );
  eta.Update( (-1.0), beta, 1.0);

}  // RelStress()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 04/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::Damage::SetupCmat(LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat)
{
  // get material parameters
  // Young's modulus (modulus of elasticity)
  double young = params_->youngs_;
  // Poisson's ratio
  double nu = params_->poissonratio_;

  // isotropic elasticity tensor C in Voigt matrix notation, cf. FEscript p.29
  //                       [ 1-nu     nu     nu |          0    0    0 ]
  //                       [        1-nu     nu |          0    0    0 ]
  //           E           [               1-nu |          0    0    0 ]
  //   C = --------------- [ ~~~~   ~~~~   ~~~~   ~~~~~~~~~~  ~~~  ~~~ ]
  //       (1+nu)*(1-2*nu) [                    | (1-2*nu)/2    0    0 ]
  //                       [                    |      (1-2*nu)/2    0 ]
  //                       [ symmetric          |           (1-2*nu)/2 ]
  //
  const double mfac = young/((1.0+nu)*(1.0-2.0*nu));  // factor

  // clear the material tangent
  cmat.Clear();
  // write non-zero components
  cmat(0,0) = mfac*(1.0-nu);
  cmat(0,1) = mfac*nu;
  cmat(0,2) = mfac*nu;
  cmat(1,0) = mfac*nu;
  cmat(1,1) = mfac*(1.0-nu);
  cmat(1,2) = mfac*nu;
  cmat(2,0) = mfac*nu;
  cmat(2,1) = mfac*nu;
  cmat(2,2) = mfac*(1.0-nu);
  // ~~~
  cmat(3,3) = mfac*0.5*(1.0-2.0*nu);
  cmat(4,4) = mfac*0.5*(1.0-2.0*nu);
  cmat(5,5) = mfac*0.5*(1.0-2.0*nu);

}  // SetupCmat()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 05/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::Damage::SetupCmatElastoPlastic(
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,  // elasto-plastic tangent modulus (out)
  double Dgamma,  // plastic multiplier
  double G,  // shear modulus
  double q,  // elastic trial von Mises effective stress
  LINALG::Matrix<NUM_STRESS_3D,1> flowvector,  // unit flow vector
  LINALG::Matrix<NUM_STRESS_3D,1> devstress,  // deviatorix stress s
  double heaviside  // Heaviside function
  )
{
  // incremental constitutive function for the stress tensor
  // sigma_n+1 = [ cmat - (Dgamma 6 G^2/q) I_d ] : strain_n+1^{e,trial}
  // consistent tangent operator
  // D^{ep} := dsigma_n+1 / dstrain_n+1^{e,trial}

  // depending on the flow vector Cmat_ep can be a fully-occupied matrix

  // C_ep = C_e - ( H^ . Dgamma . 6 . G^2 ) / q . I_d +
  //        +  H^ . 6 . G^2 ( Dgamma/q - 1/(3 G) ) N* \otimes N*
  //
  // I_d = I_s - 1/3 I . I
  // I_d in Voigt-notation applied to symmetric problem, like stress calculation
  //         [ 2/3   -1/3  -1/3 | 0    0    0  ]
  //         [-1/3    2/3  -1/3 | 0    0    0  ]
  //         [-1/3   -1/3   2/3 | 0    0    0  ]
  //   I_d = [ ~~~~  ~~~~  ~~~~  ~~~  ~~~  ~~~ ]
  //         [                  | 1/2   0   0  ]
  //         [    symmetric     |      1/2  0  ]
  //         [                  |          1/2 ]
  //

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<NUM_STRESS_3D,1> id2(true);
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> id4sharp(true);
  for (int i=0; i<3; i++) id4sharp(i,i) = 1.0;
  for (int i=3; i<6; i++) id4sharp(i,i) = 0.5;

  // unit flow vector N_ (cf. de Souza Neto (7.117)/(7.210) )
  // N_ = sqrt(2/3) . N
  double flowfac = 0.0;
  flowfac = sqrt(2.0/3.0);
  flowvector.Scale(flowfac);

  // add standard isotropic elasticity tensor C_e first
  SetupCmat(cmat);

  // if plastic loading:   heaviside = 1.0 --> use C_ep
  // if elastic unloading: heaviside = 0.0 --> use C_e
  double epfac = 0.0;
  double epfac3 = 0.0;
  // elastic trial von Mises effective stress
  if (q != 0.0)
  {
    epfac = (-1.0) * heaviside * Dgamma * 6 * G * G / q;
  }
  // constitutive tensor
  // I_d = id4sharp - 1/3 Id \otimes Id
  // contribution: Id4^#
  cmat.Update(epfac, id4sharp, 1.0);
  // contribution: Id \otimes Id
  double epfac1 = 0.0;
  epfac1 = epfac / (-3.0);
  cmat.MultiplyNT(epfac1, id2, id2, 1.0);

  // unit flow vector: N^{bar} = eta / || eta || = eta^{trial} / || eta^{trial} ||
  // factoring of shear entries not necessary for this term (devstress is used for N)
  double norm = 0.0;
  norm = sqrt( devstress(0)*devstress(0) + devstress(1)*devstress(1) + devstress(2)*devstress(2) +
                2 * ( devstress(3)*devstress(3) + devstress(4)*devstress(4) + devstress(5)*devstress(5) )
                );
  // loop strains (columns)
  for (int k=0; k<6; ++k)
  {
    // ---------------------------------------------------------- tangent
    // loop stresses (rows)
    for (int i=0; i<6; ++i)
    {
      if (q != 0.0)
      {
        epfac3 =  heaviside * 6 * G * G * ( Dgamma / q - 1.0 / (3*G) )/(norm * norm);
        cmat(i,k) += epfac3 * devstress(i) * devstress(k);
      }
    }
  }

#ifdef DEBUGMATERIAL
//  cout << "Ende SetupCmatElastPlast" << endl;
//  cout << "Cep\n" << " Dgamma " << Dgamma << endl;
//  cout << " G " << G << endl;
//  cout << " q " << q << endl;
//  cout << " flowvector " << flowvector << endl;
//  cout << " heaviside " << heaviside << endl;
//  cout << " epfac " << epfac << endl;
//  cout << " epfac1 " << epfac1 << endl;
//  cout << " cmat " << cmat << endl;
#endif // #ifdef DEBUGMATERIAL

}  // SetupCmatElastoPlastic()


void MAT::Damage::VisNames(std::map<string,int>& names)
{
  string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1; // scalar
}


bool MAT::Damage::VisData(const string& name, std::vector<double>& data, int numgp)
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    LINALG::Matrix<1,1> temp(true);
    for (int iter=0; iter<numgp; iter++)
      temp.Update(1.0,AccumulatedStrain(iter),1.0);
    data[0] = temp(0)/numgp;
  }
  return true;
}
