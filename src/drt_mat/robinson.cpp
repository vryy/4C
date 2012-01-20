/*----------------------------------------------------------------------*/
/*!
\file robinson.cpp
\brief Robinson's visco-plastic material

       example input line:
       MAT 1 MAT_Struct_Robinson  KIND Arya_NarloyZ  YOUNG POLY 2 1.47e9 -7.05e5
        NUE 0.34  DENS 8.89e-3  THEXPANS 0.0  HRDN_FACT 3.847e-12  HRDN_EXPO 4.0
        SHRTHRSHLD POLY 2 69.88e8 -0.067e8   RCVRY 6.083e-3  ACTV_ERGY 40000.0
        ACTV_TMPR 811.0  G0 0.04  M_EXPO 4.365  BETA POLY 3 0.8 0.0 0.533e-6
        H_FACT 1.67e16

  Robinson's visco-plastic material                        bborn 03/07
  material parameters
  [1] Butler, Aboudi and Pindera: "Role of the material constitutive
      model in simulating the reusable launch vehicle thrust cell
      liner response", J Aerospace Engrg, 18(1), 2005.
      --> kind = Butler
  [2] Arya: "Analytical and finite element solutions of some problems
      using a vsicoplastic model", Comput & Struct, 33(4), 1989.
      --> kind = Arya
      --> E  = 31,100 - 13.59 . T + 0.2505e-05 . T^2 - 0.2007e-13 . T^3
      --> nu = 0.254 + 0.154e-3 . T - 0.126e-06 . T^2
  [3] Arya: "Viscoplastic analysis of an experimental cylindrical
      thrust chamber liner", AIAA J, 30(3), 1992.
      --> kind = Arya_NarloyZ, Arya_CrMoSteel

      //  kind_  == "Butler"
      //         == "Arya")
      //         == "Arya_NarloyZ"
      //         == "Arya_CrMoSteel"

  this represents the backward Euler implementation established by Burkhard Bornemann

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*
 | Definitions                                               dano 11/11 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

/*----------------------------------------------------------------------*
 | Headers                                                   dano 11/11 |
 *----------------------------------------------------------------------*/
#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "robinson.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

#include "../drt_tsi/tsi_defines.H"

// include this header needed for KinematicType
#include "../drt_so3/so_hex8.H"

/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 11/11 |
 *----------------------------------------------------------------------*/
MAT::PAR::Robinson::Robinson(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  kind_((matdata->Get<std::string>("KIND"))),
  youngsnum_(matdata->GetInt("YOUNGNUM")),
  youngs_(*(matdata->Get<vector<double> >("YOUNG"))),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS")),
  thermexpans_(matdata->GetDouble("THEXPANS")),
  hrdn_fact_(matdata->GetDouble("HRDN_FACT")),
  hrdn_expo_(matdata->GetDouble("HRDN_EXPO")),
  shrthrshldnum_(matdata->GetInt("SHRTHRSHLDNUM")),
  shrthrshld_(*(matdata->Get<vector<double> >("SHRTHRSHLD"))),
  rcvry_(matdata->GetDouble("RCVRY")),
  actv_ergy_(matdata->GetDouble("ACTV_ERGY")),
  actv_tmpr_(matdata->GetDouble("ACTV_TMPR")),
  g0_(matdata->GetDouble("G0")),
  m_(matdata->GetDouble("M_EXPO")),
  betanum_(matdata->GetInt("BETANUM")),
  beta_(*(matdata->Get<vector<double> >("BETA"))),
  h_(matdata->GetDouble("H_FACT"))
{
//    cout << "matdata" << *matdata<< endl;
//    cout << "kind_ = " << *kind_ << endl;
//    cout << "youngs_ = " << (*youngs_)[0] << endl;
}


Teuchos::RCP<MAT::Material> MAT::PAR::Robinson::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Robinson(this));
}


MAT::RobinsonType MAT::RobinsonType::instance_;


DRT::ParObject* MAT::RobinsonType::Create( const std::vector<char> & data )
{
  MAT::Robinson* plastic = new MAT::Robinson();
  plastic->Unpack(data);
  return plastic;
}


/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 11/11 |
 *----------------------------------------------------------------------*/
MAT::Robinson::Robinson()
: params_(NULL)
{
}


/*----------------------------------------------------------------------*
 | copy-constructor (public)                                      11/11 |
 *----------------------------------------------------------------------*/
MAT::Robinson::Robinson(MAT::PAR::Robinson* params)
: plastic_step(false),
  params_(params)
{
}


/*----------------------------------------------------------------------*
 | pack (public)                                             dano 11/11 |
 *----------------------------------------------------------------------*/
void MAT::Robinson::Pack(DRT::PackBuffer& data) const
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
  int numgp;
  // if material is not initialized, i.e. start simulation, nothing to pack
  if (!Initialized())
  {
    numgp=0;
  }
  else
  {
    // if material is initialized (restart): size equates number of gausspoints
    numgp = strainpllast_->size();
  }
  AddtoPack(data,numgp); // Length of history vector(s)
  for (int var=0; var<numgp; ++var)
  {
    // pack history data
    AddtoPack(data,strainpllast_->at(var));
    AddtoPack(data,backstresslast_->at(var));

    AddtoPack(data,kvarva_->at(var));
    AddtoPack(data,kvakvae_->at(var));
  }

  return;

} // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                           dano 11/11 |
 *----------------------------------------------------------------------*/
void MAT::Robinson::Unpack(const vector<char>& data)
{
  isinit_=true;
  vector<char>::size_type position = 0;
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
        params_ = static_cast<MAT::PAR::Robinson*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // history data
  int numgp;
  ExtractfromPack(position,data,numgp);

  // if system is not yet initialized, the history vectors have to be intialized
  if (numgp == 0)
    isinit_ = false;

  // unpack strain vectors
  strainpllast_ = Teuchos::rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  strainplcurr_ = Teuchos::rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );

  // unpack back stress vectors (for kinematic hardening)
  backstresslast_ = Teuchos::rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );
  backstresscurr_ = Teuchos::rcp( new vector<LINALG::Matrix<NUM_STRESS_3D,1> > );

  // unpack matrices needed for condensed system
  kvarva_ = Teuchos::rcp( new vector<LINALG::Matrix<2*NUM_STRESS_3D,1> > );
  kvakvae_ = Teuchos::rcp( new vector<LINALG::Matrix<2*NUM_STRESS_3D,NUM_STRESS_3D> > );

  for (int var=0; var<numgp; ++var)
  {
    LINALG::Matrix<NUM_STRESS_3D,1> tmp(true);
    LINALG::Matrix<2*NUM_STRESS_3D,1> tmp1(true);
    LINALG::Matrix<2*NUM_STRESS_3D,NUM_STRESS_3D> tmp2(true);

    // unpack strain/stress vectors of last converged state
    ExtractfromPack(position,data,tmp);
    strainpllast_->push_back(tmp);
    ExtractfromPack(position,data,tmp);
    backstresslast_->push_back(tmp);

    // unpack matrices of last converged state
    ExtractfromPack(position,data,tmp1);
    kvarva_->push_back(tmp1);
    ExtractfromPack(position,data,tmp2);
    kvakvae_->push_back(tmp2);

    // current vectors have to be initialised
    strainplcurr_->push_back(tmp);
    backstresscurr_->push_back(tmp);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);

  return;

} // Unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal stress variables (public) dano 11/11 |
 *---------------------------------------------------------------------*/
void MAT::Robinson::Setup(
  const int numgp,
//  KinematicType kintype
  DRT::INPUT::LineDefinition* linedef
  )
{
  // temporary variable for read-in
   std::string buffer;
  // read kinematic flag (default: we expect kintype to be total lagrangian)
  linedef->ExtractString("KINTYP",buffer);
  if (buffer=="lin")
  {
   kintype_ = soh8_geolin;
  }
  else if (buffer=="nln")
  {
   kintype_ = soh8_totlag;
  }
  else dserror ("Reading SO_HEX8 element failed");

  // initialise history variables
  strainpllast_ = Teuchos::rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  strainplcurr_ = Teuchos::rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  backstresslast_ = Teuchos::rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  backstresscurr_ = Teuchos::rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  kvarva_ = Teuchos::rcp(new vector<LINALG::Matrix<2*NUM_STRESS_3D,1> >);
  kvakvae_ = Teuchos::rcp(new vector<LINALG::Matrix<2*NUM_STRESS_3D,NUM_STRESS_3D> >);

  LINALG::Matrix<NUM_STRESS_3D,1> emptymat(true);
  strainpllast_->resize(numgp);
  strainplcurr_->resize(numgp);

  backstresslast_->resize(numgp);
  backstresscurr_->resize(numgp);

  LINALG::Matrix<2*NUM_STRESS_3D,1> emptymat2(true);
  kvarva_->resize(numgp);
  LINALG::Matrix<2*NUM_STRESS_3D,NUM_STRESS_3D> emptymat3(true);
  kvakvae_->resize(numgp);

  for (int i=0; i<numgp; i++)
  {
    strainpllast_->at(i) = emptymat;
    strainplcurr_->at(i) = emptymat;

    backstresslast_->at(i) = emptymat;
    backstresscurr_->at(i) = emptymat;

    kvarva_->at(i) = emptymat2;
    kvakvae_->at(i) = emptymat3;
  }

  isinit_=true;
  return;

}  // Setup()


/*---------------------------------------------------------------------*
 | update internal stress variables (public)                dano 11/11 |
 *---------------------------------------------------------------------*/
void MAT::Robinson::Update()
{
  // make current values at time step t_n+1 to values of last step t_n
  strainpllast_ = strainplcurr_;
  backstresslast_ = backstresscurr_;

  // empty vectors of current data
  strainplcurr_ = Teuchos::rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);
  backstresscurr_ = Teuchos::rcp(new vector<LINALG::Matrix<NUM_STRESS_3D,1> >);

  kvarva_ = Teuchos::rcp(new vector<LINALG::Matrix<2*NUM_STRESS_3D,1> >);
  kvakvae_ = Teuchos::rcp(new vector<LINALG::Matrix<2*NUM_STRESS_3D,NUM_STRESS_3D> >);

  // get size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int numgp = strainpllast_->size();
  strainplcurr_->resize(numgp);
  backstresscurr_->resize(numgp);
  kvarva_->resize(numgp);
  kvakvae_->resize(numgp);

  const LINALG::Matrix<NUM_STRESS_3D,1> emptyvec(true);
  const LINALG::Matrix<1,1> emptyvec1(true);
  const LINALG::Matrix<2*NUM_STRESS_3D,1> emptyvec2(true);
  const LINALG::Matrix<2*NUM_STRESS_3D,NUM_STRESS_3D> emptyvec3(true);

  for (int i=0; i<numgp; i++)
  {
    strainplcurr_->at(i) = emptyvec;
    backstresscurr_->at(i) = emptyvec;
    kvarva_->at(i) = emptyvec2;
    kvakvae_->at(i) = emptyvec3;
  }

  return;

}  // Update()


/*---------------------------------------------------------------------*
 | reset internal stress variables (public)                 dano 11/11 |
 *---------------------------------------------------------------------*/
void MAT::Robinson::Reset()
{
  // do nothing,
  // because #histplasticrcgcurr_ and #histeplasticscurr_ are recomputed
  // anyway at every iteration based upon #histplasticrcglast_ and
  // #histeplasticslast_ untouched within time step

  return;

}  // Reset()


/*----------------------------------------------------------------------*
 | evaluate material (public)                                dano 11/11 |
 | Select Robinson's material and integrate internal variables          |
 | CCARAT: so3_mat_robinson_be_sel()                                    |
 *----------------------------------------------------------------------*/
void MAT::Robinson::Evaluate(
  const LINALG::Matrix<NUM_STRESS_3D,1>& strain,  //!< total strain vector
  LINALG::Matrix<NUM_STRESS_3D,1>& plstrain,  //!< plastic strain vector
  const int gp, //!< current Gauss point
  Teuchos::ParameterList& params,  //!< parameter list for communication & HISTORY
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat, //!< material stiffness matrix
  LINALG::Matrix<NUM_STRESS_3D,1>& stress //!< 2nd PK-stress
  )
{
  cout << "params" << params << endl;
  // name: Baci --> CCARAT:

  // STRAIN
  // total strain:   strain   --> stntotn
  // elastic strain: strain_e --> stnela
  // thermal strain: strain_t --> stnthr
  // viscous strain: strain_p --> stnvscn

  // STRESS
  // stress deviator:      devstress       --> devstsn
  // back stress:          backstress/beta --> bacsts/bacstsn
  // over/relative stress: eta             --> stsovrn

  // RESIDUALS
  // viscous strain residual: strain_pres     --> vscstnr;
  LINALG::Matrix<NUM_STRESS_3D,1> strain_pres(true);
  // back stress (residual):  beta/backstress --> bckstsr
  LINALG::Matrix<NUM_STRESS_3D,1> backstress_res(true);

  // TANGENT
  // declare single terms of elasto-plastic tangent Cmat_ep
  // tangent term resulting from linearisation \frac{\pd sig}{\pd eps^v}
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kev(true);
  // tangent term resulting from linearisation \frac{\pd sig}{\pd al}
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kea(true);
  // tangent term resulting from linearisation \frac{\pd res^v}{\pd eps}
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kve(true);
  // tangent term resulting from linearisation \frac{\pd res^v}{\pd eps^v}
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kvv(true);
  // tangent term resulting from linearisation \frac{\pd res^v}{\pd al}
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kva(true);
  // tangent term resulting from linearisation \frac{\pd res^al}{\pd eps}
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kae(true);
  // tangent term resulting from linearisation \frac{\pd res^al}{\pd eps^v}
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kav(true);
  // tangent term resulting from linearisation \frac{\pd res^al}{\pd al}
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kaa(true);

  // get material parameters
  // TODO pass time step size from element to the material
//   const double dt = params->delta time; dyn->dt;
  // Young's modulus
  // get the vector of the temperature dependent youngs modulus
  const vector<double> youngvec = params_->youngs_;
  double young = youngvec[0]; //0.0;
  // Poisson's ratio
  double nu = params_->poissonratio_;
  // initialise the thermal expansion coefficient
  const double thermexpans = params_->thermexpans_;

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
  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  LINALG::Matrix<6,6> id4(true);
  for (int i=0; i<6; i++) id4(i,i) = 1.0;

  // --------------------------------------------------- temperatures
  // get initial and current temperature
  // 12.01.12 TODO pass current, initial temperature from thermo field to the material
  // initial temperature at Gauss point
  double tempinit; // CCARAT: tempr0 --> so3_tsi_temper0()
  // current temperature at Gauss point
  double tempnp;  // CCARAT: tempr --> so3_tsi_temper()

  // ------------------------------------------------- thermal strain
  // thermal strain vector stnthr
  LINALG::Matrix<NUM_STRESS_3D,1> strain_t(true);
  for (int i=0; i<3; ++i) strain_t(i) = thermexpans * (tempnp - tempinit);
  //  strain_t(i) = 2*E_xy = 2*E_yz = 2*E_zx = 0.0
  // for (int i=3; i<6; ++i) strain_t = 0.0;

  // ------------------------------------------------- viscous strain
  // viscous strain strain_{n+1}^{v<i>} at t_{n+1}
  // strain_p --> stnvscn
  LINALG::Matrix<NUM_STRESS_3D,1> strain_p(true);
  for (int i=0; i<6; i++)
    strain_p(i,0) = strainpllast_->at(gp)(i,0);
  // OR: TODO 13.01.12 CHECK if this is correct here!!!
  // actually strain at new time, i.e., strain_p(i,0) = strainplcurr_->at(gp)(i,0);
  LINALG::Matrix<NUM_STRESS_3D,1> strain_pn(true);
  for (int i=0; i<6; i++)
    strain_pn(i,0) = strainplcurr_->at(gp)(i,0);

  // ------------------------------------------------- elastic strain
  // elastic strain at t_{n+1}
  // strain^{e}_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> strain_e(true);

  // strain^e_{n+1} = strain_n+1 - strain^p_n - strain^t
  strain_e.Update( 1.0, strain, 0.0 );
  strain_e.Update( (-1.0), strain_p, 1.0 );
  strain_e.Update( (-1.0), strain_t, 1.0 );

  // ---------------------------------------------- elasticity tensor
  // cmat = kee = pd(sig)/pd(eps)
  // CCARAT: so3_mat_robinson_elmat(mat_robin, tempr, cmat);
  // pass the current temperature to calculate the current youngs modulus
  SetupCmat(tempnp, cmat);

  // ------------------------------------ tangents of stress equation
  // kev = pd(sigma)/pd(eps^v)
  // CCARAT: so3_mv6_m_assscl(-1.0, cmat, kev);
  // assign vector by another vector and scale it
  // (i): scale (-1.0)
  // (i): input matrix cmat
  // (o): output matrix kev
  kev.Update((-1),cmat,0.0);

  // kea = pd(sigma)/pd(backstress)
  // initialise with 1 on the diagonals
  // CCARAT: so3_mv6_m_id(kea);
  kea.Update(1.0,id4,0.0);

  // ------------------------------------------------- elastic stress
  // stress sig_{n+1}^<i> at t_{n+1} */
  // (i): input matrix cmat
  // (i): input vector strain_e (CCARAT: stnela)
  // (o): output vector stress
  // CCARAT: so3_mv6_v_assmvp(cmat, stnela, stress);
  // \f \sigma = {\mathbf C} \,\varepsilon \f
  stress.MultiplyNN(cmat,strain_e);
//  // alternatively
//  LINALG::Matrix<NUM_STRESS_3D,1> stressadditive(true);
//  Stress(p, devstress, stressadditive);

  // ------------------------------------------------------ devstress
  // deviatoric stress s_{n+1}^<i> at t_{n+1}
  // calculate the deviator from stress
  // CAUTION: s = 2G . devstrain only in case of small strain
  // CCARAT: so3_mv6_v_dev(stress, devstsn);
  LINALG::Matrix<NUM_STRESS_3D,1> devstress(true);
  // trace of stress vector
  double tracestress = ( stress(0)+stress(1)+stress(2) );
  for (size_t i=0; i<3; i++)
    devstress(i) = stress(i) - tracestress/3.0;
  for (size_t i=3; i<NUM_STRESS_3D; i++)
    devstress(i) = stress(i);
  // be careful for shear stresses (e.g., sigma_12)
  // in Voigt-notation the shear strains(!) have to be scaled with 1/2
  // normally considered in material tangent (using id4sharp, instead of id4)

  // ---------------------------------------------------- back stress
  // TODO 17.01.12 not sure, if this is correct here: bckstsr represents residual of backstress
  // new back stress at beta_{n+1}^i
  // CCARAT: actso3->miv_rob->bacstsn.a.da[ip]
  LINALG::Matrix<NUM_STRESS_3D,1> backstress_np(true);
  for (size_t i=0; i<NUM_STRESS_3D; i++)
    backstress_np(i,0) = backstresscurr_->at(gp)(i,0);
  // old back stress at beta_{n}
  LINALG::Matrix<NUM_STRESS_3D,1> backstress_n(true);
  for (size_t i=0; i<NUM_STRESS_3D; i++)
    backstress_n(i,0) = backstresslast_->at(gp)(i,0);

  // ------------------------------------------ over/relativestress
  // overstress Sig_{n+1}^i = s_{n+1}^i - al_{n+1}^i
  // (i): input vector devstn
  // (i): input vector beta
  // (o): output vector stsovr: subtract 2 vectors
  // CCARAT: so3_mv6_v_sub(devstsn, actso3->miv_rob->bacstsn.a.da[ip], stsovrn);
  // eta_{n+1} = devstress_{n+1} - backstress_{n+1}
  LINALG::Matrix<NUM_STRESS_3D,1> eta(true);
  RelDevStress( devstress, backstress_np, eta);

  // ------------------------------------- residual of viscous strain
  // residual of visc. strain eps_{n+1}^<i> and its consistent tangent for <i>
  // CCARAT: so3_mat_robinson_be_rvscstn(ele, mat_robin, dt, tempr,
  //                                     actso3->miv_rob->vicstn.a.da[ip],
  //                                     actso3->miv_rob->vicstnn.a.da[ip],
  //                                     devstsn, stsovrn,
  //                                     &(actso3->miv_rob->vscstns.a.iv[ip]),
  //                                     vscstnr, kve, kvv, kva);
  // TODO 13.01.12 pass the time step from element to the material, or implement
  // possibility to ask element from the material
  const double dt = 0.1;
  CalcBEViscousStrainRate(
    dt,
    tempnp,
    strain_p,
    strain_pn,
    devstress,
    eta,
    strain_pres,
    kve,
    kvv,
    kva
    );

  // ---------------------------------------- residual of back stress
  // residual of back stress al_{n+1} and its consistent tangent
  // CCARAT: so3_mat_robinson_be_rbcksts(ele, mat_robin, dt, tempr,
  //                                     actso3->miv_rob->vicstn.a.da[ip],
  //                                     actso3->miv_rob->vicstnn.a.da[ip],
  //                                     devstsn,
  //                                     actso3->miv_rob->bacsts.a.da[ip],
  //                                     actso3->miv_rob->bacstsn.a.da[ip],
  //                                     &(actso3->miv_rob->bckstss.a.iv[ip]), --> no backstressmode in baci
  //                                     bckstsr, kae, kav, kaa);
  CalcBEBackStressFlow(
    dt,
    tempnp,
    strain_p,
    strain_pn,
    devstress,
    backstress_n,
    backstress_np,
    backstress_res,
    kae,
    kav,
    kaa
    );

  // ------------------------------------- reduced stress and tangent
  // build reduced stress and tangent
  // ==> static condensation
  // CCARAT: so3_mat_robinson_be_red(stress, cmat, kev, kea,
  //                          vscstnr, kve, kvv, kva,
  //                          bckstsr, kae, kav, kaa,
  //                          actso3->miv_rob->kvarva.a.da[ip],
  //                          actso3->miv_rob->kvakvae.a.da[ip]);
   CalculateCondensedSystem(
     stress,
     cmat,
     kev,
     kea,
     strain_pres,
     kve,
     kvv,
     kva,
     backstress_res,
     kae,
     kav,
     kaa,
     kvarva_->at(gp),
     kvakvae_->at(gp)
     );

}  // Evaluate()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 11/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::Robinson::SetupCmat(
  double tempnp,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat
  )
{
  // get material parameters
  // Young's modulus
  // get the vector of the temperature dependent youngs modulus
  // Young's modulus (modulus of elasticity)
  double emod;
  GetMatParameterAtTempnp(&(params_->youngs_), tempnp, emod);
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
  const double mfac = emod/((1.0+nu)*(1.0-2.0*nu));  // factor

  // clear the material tangent
  cmat.Clear();
  // write non-zero components --- axial
  cmat(0,0) = mfac*(1.0-nu);
  cmat(0,1) = mfac*nu;
  cmat(0,2) = mfac*nu;
  cmat(1,0) = mfac*nu;
  cmat(1,1) = mfac*(1.0-nu);
  cmat(1,2) = mfac*nu;
  cmat(2,0) = mfac*nu;
  cmat(2,1) = mfac*nu;
  cmat(2,2) = mfac*(1.0-nu);
  // write non-zero components --- shear
  cmat(3,3) = mfac*0.5*(1.0-2.0*nu);
  cmat(4,4) = mfac*0.5*(1.0-2.0*nu);
  cmat(5,5) = mfac*0.5*(1.0-2.0*nu);

}  // SetupCmat()


/*----------------------------------------------------------------------*
 | computes linear stress tensor                             dano 11/11 |
 *----------------------------------------------------------------------*/
void MAT::Robinson::Stress(
  const double p,  //!< volumetric stress
  const LINALG::Matrix<NUM_STRESS_3D,1>& devstress,  //!< deviatoric stress tensor
  LINALG::Matrix<NUM_STRESS_3D,1>& stress //!< 2nd PK-stress
  )
{
  // total stress = deviatoric + hydrostatic pressure . I
  // sigma = s + p . I
  stress.Update(1.0, devstress, 0.0);
  for (int i=0; i<3; ++i) stress(i) += p;

}  // Stress()


/*----------------------------------------------------------------------*
 | compute relative deviatoric stress tensor                 dano 11/11 |
 *----------------------------------------------------------------------*/
void MAT::Robinson::RelDevStress(
  const LINALG::Matrix<NUM_STRESS_3D,1>& devstress,  //!< deviatoric stress tensor
  const LINALG::Matrix<NUM_STRESS_3D,1>& beta,  //!< back stress tensor
  LINALG::Matrix<NUM_STRESS_3D,1>& eta //!< relative stress
  )
{
  // relative stress = deviatoric - back stress
  // eta = s - beta
  eta.Update( 1.0, devstress, 0.0 );
  eta.Update( (-1.0), beta, 1.0);

}  // RelDevStress()


/*----------------------------------------------------------------------*
 | residual of BE-discretised viscous strain rate            dano 11/11 |
 | at Gauss point                                                       |
 *----------------------------------------------------------------------*/
void MAT::Robinson::CalcBEViscousStrainRate(
//  \param   ele          ELEMENT*      (i)   element
//  \param   mat_robin    VP_ROBINSON*  (i)   element Robinson material
  const double dt,  // (i) time step size --> BE CAREFUL: dt=0.1;
  // TODO 12.01.12 check what temperature is needed here
  double tempnp,  //       DOUBLE        (i)   temperature
  const LINALG::Matrix<NUM_STRESS_3D,1>&  strain_p,  // DOUBLE[]      (i)   viscous strain at t_n
  const LINALG::Matrix<NUM_STRESS_3D,1>&  strain_pn,  // DOUBLE[]      (i)   viscous strain at t_{n+1}^<i>
  const LINALG::Matrix<NUM_STRESS_3D,1>&  devstress,  // DOUBLE[]      (i)   stress deviator at t_{n+1}^<i>
  const LINALG::Matrix<NUM_STRESS_3D,1>&  eta,  // DOUBLE[]      (i)   over stress at t_{n+1}^<i>
  LINALG::Matrix<NUM_STRESS_3D,1>& strain_pres,  // DOUBLE[]      (o)   viscous strain residual
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kve,  // (o)
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kvv,  // (o)
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kva  // (o)
  )
{
  // initialise
  // J2-invariant
  // CCARAT: 'J_2' = j2
  double J2 = 0.0;
  // shear threshold 'K^2'
  double kksq = 0.0;;
  // hardening factor 'A'
  double aa = 0.0;;
  // hardening exponent 'N'
  // TODO nn so definiert? 13.01.12
  double nn = 0.0;
  // 'F'
  double ff = 0.0;;
  // ss = 1/2 * eta : devstress
  // CCARAT: ss = 1/2 * s : Sig
  double ss;

  // viscous strain rate
  // CCARAT: dvscstn[NUMSTR_SOLID3]
  // strainrate_p --> dvscstn
  LINALG::Matrix<NUM_STRESS_3D,1> strainrate_p(true);
  // kvs = d strain_p / d eta
  // CCARAT: d eps^v/d Sig
  // kvs[NUMSTR_SOLID3][NUMSTR_SOLID3];
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kvs(true);

  // identity tensor in vector notation
  LINALG::Matrix<NUM_STRESS_3D,1> id2(true);
  for (int i=0; i<3; i++) id2(i) = 1.0;

  // -------------------------------------------------- preliminaries
  // J2 = 1/2 eta : eta
  // J2 = 1/2 (eta11^2 + eta22^2 + eta33^2 + 2 . eta12^2 + 2 . eta23^2 + 2 . eta13^2)
  // CCARAT: 1.) so3_mv6_v_dblctr(ovrstsn, ovrstsn, &(j2));
  //         2.) j2 *= 0.5;
  //         --> J_2 = 1/2 * Sig : Sig  with Sig...overstress
  J2 = 1/2.0 * ( eta(0)*eta(0) + eta(1)*eta(1) + eta(2)*eta(2) ) +
       + eta(3)*eta(3) + eta(4)*eta(4) + eta(5)*eta(5);

  // Bingham-Prager shear stress threshold at current temperature
  // calculate non-linear material parameter
  // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->shrthrshld, tmpr, &(kksq));
  GetMatParameterAtTempnp(&(params_->shrthrshld_), tempnp, kksq);

  // F = (J_2 - K^2)/K_2
  if (fabs(kksq) <= EPS10)  // shrthrshld = kksq
  {
    ff = -1.0;
    dserror("Division by zero: Shear threshold very close to zero");
  }
  else
  {
    // CCARAT: (J2 - kksq)/kksq;
    ff = (J2 - kksq)/kksq;
  }

  // ss = 1/2 * devstress : eta
  // with  eta: overstress, devstress:.deviat.stress
  // CCARAT: double contraction of 2vectors: so3_mv6_v_dblctr(ovrstsn, devstsn, &(ss));
  ss = 1/2.0 * ( devstress(0)*eta(0) + devstress(1)*eta(1) + devstress(2)*eta(2) ) +
         + devstress(3)*eta(3) + devstress(4)*eta(4) + devstress(5)*eta(5);

  // hardening factor aa
  // calculate the temperature dependent material constant \bar{\mu} := aa
  if (*(params_->kind_) == "Arya_CrMoSteel")
  {
    double mu = params_->hrdn_fact_;
    // calculate theta1 used for the material constant \bar{\mu}
    // \bar{\mu} = (23.8 . tempnp - 2635.0) . (1.0/811.0 - 1.0/tempnp)), vgl. (14)
    double th1 = (23.8*tempnp - 2635.0)*(1.0/811.0 - 1.0/tempnp);
    if (isinf(th1))
    {
      dserror("Infinite theta1");
    }
    // \bar{\mu} = \mu . exp(- theta1), vgl. (12)
    aa = 0.5/(mu*exp(-th1));
  }
  else
  {
    aa = params_->hrdn_fact_;
  }

  // hardening exponent nn
  // TODO 13.01.12 if possible define at the beginning of the function
  nn = params_->hrdn_expo_;

  //-------------------------------------------------------------------
  // check plastic admissibility
  //-------------------------------------------------------------------

  // ------------------------------------------------- determine mode
  // The mode is determined for every time step once. The iteration sticks to
  // the initially determined mode.
  // CCARAT: if (*vscstns == so3_mat_robinson_state_vague)
  //         if (*vscstns == so3_mat_robinson_state_elastic)

  // the string is set to so3_mat_robinson_state_vague, before we check if the
  // load step is inelastic or elastic
  // CCARAT: vscstns --> strainratemode
  //  if (strainratemode == so3_mat_robinson_state_vague)
  //  if (*vscstns == so3_mat_robinson_state_vague)

  //-------------------------------------------------------------------
  // IF plastic step ( Phi_trial > 0.0, Dgamma >= 0.0 )
  // violated consistency condition
  // TODO 13.01.12: ff = Phi_trial, ss = Dgamma???
  //-------------------------------------------------------------------
  // --------------------------------------------------- viscous mode
  // residual of viscous strain rate at t_{n+1}
  // CCARAT: else if (*vscstns == so3_mat_robinson_state_inelastic)
  if ( (ff > 0.0) && (ss > 0.0) )
  {
    // inelastic/viscous strain residual:
    // CCARAT: double fct = aa * pow(ff, nn) / sqrt(j2);
    // 10.11.11 TODO what does 'pow'???
    double fct = aa * pow(ff, nn) / sqrt(J2);
    // calculate the viscous strain rate respecting that strain vector component
    // has a doubled shear component, but stress vectors not!
    // --> scale shear components accordingly
    // CCARAT: so3_mv6_v2_assscl(fct, ovrstsn, dvscstn);
    // strainrate_p --> dvscstn
    for (size_t i=0; i<3; i++)
    {
      strainrate_p(i) = eta(i);
    }
    for (size_t i=3; i<NUM_STRESS_3D; i++)
    {
      strainrate_p(i) = 2.0 * eta(i);
    }
    strainrate_p.Scale(fct);
  }
  //-------------------------------------------------------------------
  // ELSE IF elastic step (Phi_trial <= 0.0, Dgamma = 0.0)
  //-------------------------------------------------------------------
  else   // CCARAT: if (*vscstns == so3_mat_robinson_state_elastic)
  {
    // elastic step, no inelastic strains: strain_n^v' == 0
    // CCARAT: so3_mv6_v_zero(dvscstn);
    // --> dvstcstn == strainrate_p
    strainrate_p.Scale(0.0);
  }

  //-------------------------------------------------------------------
  // ----------------------------------------------- viscous residual
  //-------------------------------------------------------------------
  // res_{n+1}^v = (strain_{n+1}^v - strain_n^v)/dt - deps_{n+1}^v
  {
    for (size_t i=0; i<NUM_STRESS_3D; i++)
    {
      // strain_pres = 1/dt * (strain_{n+1}^v - strain_{n}^v - dt * strain_n^v')
      // CCARAT: vscstnr[istr] = (vscstnn[istr] - vscstn[istr] - dt*dvscstn[istr])/dt;
      strain_pres(i) = (strain_pn(i) - strain_p(i) - dt*strainrate_p(i))/dt;
    }
  }

  //-------------------------------------------------------------------
  // --------------------------------- derivative of viscous residual
  // derivative of viscous residual with respect to over stress eta
  // CCARAT: kvs = pd(res^v)/pd(Sig) -- d eps^v/d Sig
  //  kvs = d(strain_pres) / d (eta)

  // IF plastic step ( Phi_trial > 0.0, Dgamma >= 0.0 )
  if ( (ff > 0.0) && (ss > 0.0) ) //  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    // add facu to all diagonal terms matrix kvs
    // CCARAT: so3_mv6_m_idscl(facu, kvs);
    double facu = -aa * pow(ff, nn) / sqrt(J2);
    for (size_t i=0; i<NUM_STRESS_3D; i++)
      kvs(i,i) = facu;

    // update matrix with scaled dyadic vector product
    // CCARAT: so3_mv6_m_upddydscl(faco, ovrstsn, ovrstsn, kvs); // Update matrix with scaled dyadic vector product
    // contribution: kvs = kvs + faco . (eta \otimes eta^T)
    double faco = -nn * aa * pow(ff, (nn-1.0)) / (kksq * sqrt(J2))
                  + aa * pow(ff, nn) / (2.0 * pow(J2, 1.5));
    kvs.MultiplyNT(faco,eta,eta,1.0);
    // TODO 13.01.12 CHECK if same result is obtained
//    for (size_t i=0; i<NUM_STRESS_3D; i++)
//    {
//      for (size_t j=0; j<NUM_STRESS_3D; j++)
//      {
//        kvs(i,j) += faco * eta(i) * eta(j);
//      }
//    }

    // multiply last 3 rows by 2 to conform with definition of strain vectors
    // CCARAT: so3_mv6_m2_updmtom2(kvs);
    // consider the difference between physical strains and Voigt notation
    for (size_t i=3; i<NUM_STRESS_3D; i++)
    {
      for (size_t j=0; j<NUM_STRESS_3D; j++)
      {
        kvs(i,j) *= 2.0;
      }
    }

  } // inelastic
  // ELSE IF elastic step (Phi_trial <= 0.0, Dgamma = 0.0)
  else //  (*vscstns == so3_mat_robinson_state_elastic)
  {
    //so3_mv6_m_zero(kvs);
    // in case of an elastic step, no contribution to kvs
    kvs.Scale(0.0);
  }

  //-------------------------------------------------------------------
  // derivative of viscous residual with respect to total strain eps
  // CCARAT: kve = ( pd(res_{n+1}^v)/pd(eps_{n+1}) )|^<i>
  // kve = ( pd strain_pres^{n+1} )/ (pd strain^{n+1})|^<i>

  // IF plastic step ( Phi_trial > 0.0, Dgamma >= 0.0 )
  if ( (ff > 0.0) && (ss > 0.0) ) //  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    // calculate elastic material tangent with temperature-dependent Young's modulus
    // CCARAT: so3_mat_robinson_elmat(mat_robin, tmpr, kse);
    // kse = (pd eta) / (pd strain)
    // CCARAT: kse = pd(Sig)/pd(eps) */
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> kse(true);
    // pass the current temperature to calculate the current youngs modulus
    SetupCmat(tempnp, kse);

//    // build Cartesian identity 2-tensor I_{AB}
//    // CCARAT: so3_mv6_v_id(iv);
//    LINALG::Matrix<NUM_STRESS_3D,1> id2(true);
//    for (int i=0; i<3; i++) id2(i) = 1.0;

    // Matrix vector product: cid2 = kse(i,j)*id2(j)
    // CCARAT: so3_mv6_v_assmvp(kse, iv, civ);
    LINALG::Matrix<NUM_STRESS_3D,1> cid2(true);
    cid2.Multiply(1.0, kse, id2, 0.0);
    // update matrix with scaled dyadic vector product
    // CCARAT: so3_mv6_m_upddydscl(-1.0/3.0, iv, civ, kse);
    // contribution: kse = kse + (-1/3) . (id2 \otimes cid2^T)
    kse.MultiplyNT((-1.0/3.0),id2,cid2,1.0);
//    // TODO 13.01.12 check if equivalent to MultiplyNT
//    for (i=0; i<NUM_STRESS_3D; i++)
//    {
//      for (j=0; j<NUM_STRESS_3D; j++)
//      {
//        kse(i,j) += (-1.0/3.0) * id2(i) * cid2(j);
//      }
//    }

    // Assign kve by matrix-matrix product kvs, kse (inner product)
    // kve = kvs . kse
    // CCARAT: so3_mv6_m_mprd(kvs, kse, kve);
// TODO 22.11.11: is this the correct option to multiply these 2 tensors??
    // kve[i][j] += kvs[i][k] . kse[k][j] kve(1,2) = kvs(1,k).kse(k,j)
    kve.Multiply(1.0, kvs, kse, 0.0);
//    for (i=0; i<NUM_STRESS_3D; i++)
//    {
//      for (j=0; j<NUM_STRESS_3D; j++)
//      {
//        rcsum = 0.0;
//        for (k=0; k<NUM_STRESS_3D; k++)
//        {
//          rcsum += kvs(i,k) * kse(k,j);
//        }
//        kve(i,j) = rcsum;
//      }
//    }
  } // inelastic

  // ELSE IF elastic step
  else //  (*vscstns == so3_mat_robinson_state_elastic)
  {
    // in case of an elastic step, no contribution to kvs
    // CCARAT: so3_mv6_m_zero(kve);
    kve.Scale(0.0);
  }

  //-------------------------------------------------------------------
  // derivative of viscous residual with respect to viscous strain strain_p
  // IF plastic step ( Phi_trial > 0.0, Dgamma >= 0.0 )
  if ( (ff > 0.0) && (ss > 0.0) ) //  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    // derivative ksv = (pd eta) / (pd strain_p)
    // CCARAT: ksv = pd(Sig)/pd(eps^v)
    LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D> ksv(true);
    // CCARAT: so3_mat_robinson_elmat(mat_robin, tmpr, ksv);
    // pass the current temperature to calculate the current youngs modulus
    SetupCmat(tempnp, ksv);

    // Matrix vector product: cid2 = kse(i,j)*id2(j)
    // CCARAT: so3_mv6_v_assmvp(ksv, iv, civ);
    LINALG::Matrix<NUM_STRESS_3D,1> cid2(true);
    cid2.Multiply(1.0, ksv, id2, 0.0);
    // update matrix with scaled dyadic vector product
    // CCARAT: so3_mv6_m_upddydscl(-1.0/3.0, iv, civ, ksv);
    // contribution: ksv = ksv + (-1/3) . (id2 \otimes cid2^T)
    ksv.MultiplyNT((-1.0/3.0),id2,cid2,1.0);
//    // TODO 13.01.12 check if equivalent to MultiplyNT
//    for (i=0; i<NUM_STRESS_3D; i++)
//    {
//      for (j=0; j<NUM_STRESS_3D; j++)
//      {
//        ksv(i,j) += (-1.0/3.0) * id2(i) * cid2(j);
//      }
//    }

    // kvv = d(res^v)/d(esp^v) + pd(res^v)/pd(Sig) . pd(Sig)/pd(eps^v)
    // which gives
    // kvv = 1/dt * Id  +  kvs . ksv
    // scale matrix kvv with 1.0/dt
    // CCARAT: so3_mv6_m_idscl(1.0/dt, kvv);
    for (size_t i=0; i<NUM_STRESS_3D; i++)
      kvv(i,i) = 1.0/dt;
    // assign matrix kvv by scaled matrix-matrix product kvs, ksv (inner product)
    // kvv = kvv + (-1.0) . kvs . ksv;
    // CCARAT: so3_mv6_m_updmprdscl(-1.0, kvs, ksv, kvv);
    kvv.Multiply((-1.0), kvs, ksv, 1.0);
    // TODO 13.01.12 check if Multiply equivalent to
//    for (size_t i=0; i<NUM_STRESS_3D; i++)
//    {
//      for (size_t j=0; j<NUM_STRESS_3D; j++)
//      {
//        double rcsum = 0.0;
//        for (size_t k=0; k<NUM_STRESS_3D; k++)
//        {
//          rcsum += (-1.0) * kvs(i,k) * ksv(k,j);
//        }
//        cm(i,j) += rcsum;
//      }
//    }
  }
  // ELSE IF elastic step
  else //  (*vscstns == so3_mat_robinson_state_elastic)
  {
    // so3_mv6_m_idscl(1.0/dt, kvv);  // scale a matrix kvv with 1.0/dt
    for (size_t i=0; i<NUM_STRESS_3D; i++)
      kvv(i,i) = 1.0/dt;
  }

  //-------------------------------------------------------------------
  // derivative of viscous residual with respect to back stress beta
  // beta --> alpha
  // kva = (pd res_{n+1}^v) / (pd beta)
  // CCARAT: kva = pd(res^v)/pd(al)
  // IF plastic step ( Phi_trial > 0.0, Dgamma >= 0.0 )
  if ( (ff > 0.0) && (ss > 0.0) ) //  else if (*vscstns == so3_mat_robinson_state_inelastic)
  {
    // assign vector by another vector and scale it
    // kva = kvs . ksa = kva . (-Id)
    // with ksa = (pd eta) / (pd beta) = - Id; eta = s - beta
    // CCARAT: ksa = pd(Sig)/pd(al) = -Id
    kva.Update(-1.0, kvs, 0.0);  // so3_mv6_m_assscl(-1.0, kvs, kva);
  }
  // ELSE IF elastic step
  else // (*vscstns == so3_mat_robinson_state_elastic)
  {
    for (int i=0; i<6; i++)
      kva(i,i) = 0.0;  // so3_mv6_m_zero(kva);
  }

}  // CalcBEViscousStrainRate()


/*----------------------------------------------------------------------*
 | residual of BE-discretised back stress and its            dano 11/11 |
 | consistent tangent                                                   |
 | residual of BE-discretised back stress rate flow rule at Gauss point |
 *----------------------------------------------------------------------*/
void MAT::Robinson::CalcBEBackStressFlow(
  const double dt,
  const double tempnp,
  const LINALG::Matrix<NUM_STRESS_3D,1>&  strain_p,  // double vscstn[NUMSTR_SOLID3],
  const LINALG::Matrix<NUM_STRESS_3D,1>&  strain_pn,  // double vscstnn[NUMSTR_SOLID3],
  const LINALG::Matrix<NUM_STRESS_3D,1>&  devstress,  // double devstsn[NUMSTR_SOLID3],
  const LINALG::Matrix<NUM_STRESS_3D,1>& backstresslast, //  const double bacsts[NUMSTR_SOLID3],
  const LINALG::Matrix<NUM_STRESS_3D,1>& backstresscurr, //const double bacstsn[NUMSTR_SOLID3],
  //TODO 22.11.11INT* bckstss --> mode is not needed in BACI,
  LINALG::Matrix<NUM_STRESS_3D,1>& backstress_res, // double bckstsr[NUMSTR_SOLID3],
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kae, // double kae[NUMSTR_SOLID3][NUMSTR_SOLID3],
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kav, // double kav[NUMSTR_SOLID3][NUMSTR_SOLID3],
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kaa // double kaa[NUMSTR_SOLID3][NUMSTR_SOLID3])
  )
{
  // get parameters
  // 'I_2'
  double i2;
  // activation temperature
  double tem0;
  // shear threshold 'K_0^2'
  double kk0sq;
  // 'H'
  double hh;
  // 'beta'
  double beta;
  // 'm'
  double mm;
  // activation energy
  double q0;
  // recovery factor 'R_0'
  double rr0;
  // recovery term 'R'
  double rr;
  // 'G_0'
  double gg0;
  // 'G'
  double gg;
  // double contraction of backstress and devstress, scale by 1/2
  double sa;

  // -------------------------------------------------- preliminaries
  // I_2 = 1/2 * Alpha : Alpha  with Alpha...back stress
  // CCARAT: so3_mv6_v_dblctr(backstresscurr, backstresscurr, &i2);
  // CCARAT: i2 *= 0.5;
  i2 = 1/2.0 * ( backstresscurr(0)*backstresscurr(0) +
                 backstresscurr(1)*backstresscurr(1) +
                 backstresscurr(2)*backstresscurr(2)
                 ) +
       + backstresscurr(3)*backstresscurr(3) +
       + backstresscurr(4)*backstresscurr(4) +
       + backstresscurr(5)*backstresscurr(5);

  // activation temperature
  tem0 = params_->actv_tmpr_;
  // Bingham-Prager shear stress threshold at activation temperature
  // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->shrthrshld, tem0, &(kk0sq));
  GetMatParameterAtTempnp(&(params_->shrthrshld_), tem0, kk0sq);

  // 'beta' at current temperature
  // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->beta, tmpr, &(beta));
  GetMatParameterAtTempnp(&(params_->beta_), tempnp, beta);

  // 'm'
  mm = params_->m_;
  // 'H' at current temperature
  // CCARAT: if (mat_robin->kind == vp_robinson_kind_arya_narloyz)
  if (*(params_->kind_) == "Arya_NarloyZ")
  {
    // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->h, tmpr, &(hh));
    GetMatParameterAtTempnp(params_->h_, tempnp, hh);
    hh *= pow(6.896,1.0+beta) / (3.0*kk0sq);
  }
  // CCARAT: else if (mat_robin->kind == vp_robinson_kind_arya_crmosteel)
  if (*(params_->kind_) == "Arya_CrMoSteel")
  {
    double mu = params_->hrdn_fact_;
    // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->h, tmpr, &(hh));
    GetMatParameterAtTempnp(params_->h_, tempnp, hh);
    hh *= 2.0 * mu;
  }
  else
  {
    // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->h, tmpr, &(hh));
    GetMatParameterAtTempnp(params_->h_, tempnp, hh);
  }

  // 'Q_0'
  q0 = params_->actv_ergy_;

  // recovery factor 'R_0'
  // CCARAT: if (mat_robin->kind == vp_robinson_kind_arya_narloyz)
  if (*(params_->kind_) == "Arya_NarloyZ")
  {
    // pressure unit scale : cN/cm^2 = 10^-4 MPa
    const double pus = 1.0e-4;
    // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->rcvry, tmpr, &(rr0));
    GetMatParameterAtTempnp(params_->rcvry_, tempnp, rr0);
    rr0 *= pow(6.896,1.0+beta+mm) * pow(3.0*kk0sq*pus*pus,mm-beta);  /* TODO '*' of '/' ? */
  }
  else
  {
    // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->rcvry, tmpr, &(rr0));
    GetMatParameterAtTempnp(params_->rcvry_, tempnp, rr0);
  }

  // 'R'
  // T_{n+1} . T_0 < (1.0E-12)
  if (fabs(tempnp*tem0) <= EPS12)
  {
    // T_0 < (1.0E-12)
    if (fabs(tem0) <= EPS12)
    {
      rr = rr0;
    }
    // T_0 > (1.0E-12)
    else
    {
      rr = rr0 * exp(q0/tem0);
    }
  }
  // T_{n+1} . T_0 > (1.0E-12)
  else
  {
    rr = rr0 * exp(q0*(tempnp-tem0)/(tempnp*tem0));
    if (isinf(rr))
    {
      rr = rr0;
    }
  }

  // 'G_0'
  gg0 = params_->g0_;

  // G = I_2/K_0^2
  // K_0^2 < 1.0E-10
  if (fabs(kk0sq) <= EPS10)
  {
    gg = 0.0;
    dserror("Division by zero: Shear threshold very close to zero");
  }
  // K_0^2 > 1.0E-10
  else
  {
    gg = sqrt(i2/kk0sq);
  }

  // sa = 1/2 * devstress : backstresscurr
  // CCARAT: sa = 1/2 . Alpha . s, with Alpha:bacstsn, s: devstsn
  // CCARAT: so3_mv6_v_dblctr(bacstsn, devstsn, &(sa));
  // CCARAT: sa *= 0.5;
  sa = 1/2.0 * ( backstresscurr(0)*devstress(0) + backstresscurr(1)*devstress(1)
                 + backstresscurr(2)*devstress(2) ) +
       + backstresscurr(3)*devstress(3) +
       + backstresscurr(4)*devstress(4) +
       + backstresscurr(5)*devstress(5);

  // ------------------------------------------------- determine mode
  // viscous mode
  // CCARAT: *bckstss = so3_mat_robinson_state_inelastic;
  if ( (gg > gg0) && (sa > 0.0) )
  {

  }
  // elastic mode
  // CCARAT: *bckstss = so3_mat_robinson_state_elastic;
  else
  {

  }

  // ----------------- difference of current and last viscous strains
  // Delta strain_p = strain_pn - strain_p
  //  \incr \eps^v = \eps_{n+1}^v - \eps_{n}^v
  //  with halved entries to conform with stress vectors */
  // CCARAT: so3_mv6_v_sub(vscstnn, vscstn, vscstnd05);
  LINALG::Matrix<NUM_STRESS_3D,1> strain_pd05(true);
  strain_pd05.Update(1.0, strain_pn, 0.0);
  strain_pd05.Update((-1.0), strain_p, 1.0);
  // Due to the fact that strain vector component have a doubled shear
  // components, i.e.
  //   strain = [ a11 a22 a33 | 2*a12 2*a23 2*a31 ]
  // but stress vectors have not, i.e.
  //   stress = [ a11 a22 a33 | a12 a23 a31 ]
  // we need to scale the last three entries.
  // CCARAT: so3_mv6_v05_updvtov05(vscstnd05);
  for (size_t i=3; i<NUM_STRESS_3D; i++)
  {
    strain_pd05(i) *= 0.5;
  }

  // ----------------------------------- residual of back stress rate
  // viscous mode
  // CCARAT: *bckstss = so3_mat_robinson_state_inelastic;
  if ( (gg > gg0) && (sa > 0.0) )
  {
    double fctv = hh / pow(gg, beta);
    double fcta = rr * pow(gg, (mm-beta)) / sqrt(i2);
    for (size_t i=0; i<NUM_STRESS_3D; i++)
    {
      backstress_res(i) = backstresscurr(i) - backstresslast(i)
                          - fctv * strain_pd05(i)
                          + dt * fcta * backstresscurr(i);
    }
    backstress_res.Scale(1.0 / dt);
  }
  // elastic mode
  // CCARAT: *bckstss = so3_mat_robinson_state_elastic;
  else
  {
    double fctv = hh / pow(gg0, beta);
    double fcta;
    if (sqrt(i2) < EPS10)
    {
      // sqrt(i2) := 1.0e6 assures units are OK
      fcta = rr * pow(gg0, (mm-beta)) / 1.0e6;
    }
    else  // (sqrt(i2) > EPS10)
    {
      fcta = rr * pow(gg0, (mm-beta)) / sqrt(i2);
    }

    for (size_t i=0; i<NUM_STRESS_3D; i++)
    {
      backstress_res(i) = backstresscurr(i) - backstresslast(i)
                          - fctv * strain_pd05(i)
                          + dt * fcta * backstresscurr(i);
    }
    backstress_res.Scale(1.0 / dt);
  }

  // ----------------------------- derivative of back stress residual
  // ----------------------- derivative with respect to total strains
  // kae = pd(res^al)/pd(eps)
  // CCARAT: if ( (*bckstss == so3_mat_robinson_state_elastic)
  //            || (*bckstss == so3_mat_robinson_state_inelastic) )
  {
    // CCARAT: so3_mv6_m_zero(kae);
    kae.Scale(0.0);
  }

  // --------------------- derivative with respect to viscous strains
  // kav = pd(res_{n+1}^al)/pd(eps_{n+1}^v)
  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  LINALG::Matrix<6,6> id4sharp(true);
  for (size_t i=0; i<3; i++) id4sharp(i,i) = 1.0;
  for (size_t i=3; i<6; i++) id4sharp(i,i) = 0.5;

  // CCARAT: else if (*bckstss == so3_mat_robinson_state_inelastic)
  if ( (gg > gg0) && (sa > 0.0) )
  {
    double fctv = -hh / (pow(gg, beta) * dt);
    // CCARAT: so3_mv6_m05_idscl(fctv, kav);
    id4sharp.Scale(fctv);
  }
  // CCARAT: if (*bckstss == so3_mat_robinson_state_elastic)
  else
  {
    double fctv = -hh / (pow(gg0, beta) * dt);
    // CCARAT: so3_mv6_m05_idscl(fctv, kav);
    id4sharp.Scale(fctv);
  }
  // copy the terms on kav
  kav.Update(1.0, id4sharp, 0.0);

  // ------------------------- derivative with respect to back stress
  // derivative of back stress residual with respect to back stress
  // kaa = pd(res_{n+1}^al)/pd(al_{n+1})
  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  LINALG::Matrix<6,6> id4(true);
  for (int i=0; i<6; i++) id4(i,i) = 1.0;

  // CCARAT: else if (*bckstss == so3_mat_robinson_state_inelastic)
  if ( (gg > gg0) && (sa > 0.0) )
  {
    double fctu = 1.0/dt  +  rr * pow(gg,(mm-beta)) / sqrt(i2);
    double fctv = beta * hh / ( pow(gg,(beta+1.0)) * dt * kk0sq );
    double fcta = rr * (mm-beta) * pow(gg,(mm-beta-1.0)) / (sqrt(i2) * kk0sq)
                - rr * pow(gg,(mm-beta)) / (2.0 * pow(i2, 1.5));
    // CCARAT: so3_mv6_m_idscl(fctu, kaa);
    id4.Scale(fctu);
    kaa.Update(1.0, id4, 0.0);
    // CCARAT: so3_mv6_m_upddydscl(fctv, vscstnd05, bacstsn, kaa);
    kaa.MultiplyNT(fctv,strain_pd05,backstresscurr,1.0);
    // CCARAT: so3_mv6_m_upddydscl(fcta, bacstsn, bacstsn, kaa);
    kaa.MultiplyNT(fcta,backstresscurr,backstresscurr,1.0);
  }
  else // CCARAT: (*bckstss == so3_mat_robinson_state_elastic)
  {
    double ii2;
    if (sqrt(i2) < EPS10)
    {
      ii2 = 1.0e12;  /* sqrt(i2) := 1.0e6 assures units are OK */
    }
    else
    {
      ii2 = i2;
    }
    double fctu = 1.0/dt + rr * pow(gg0,(mm-beta)) / sqrt(ii2);
    double fcta = -rr * pow(gg0,(mm-beta)) / (2.0 * pow(ii2, 1.5));
    // CCARAT: so3_mv6_m_idscl(fctu, kaa);
    id4.Scale(fctu);
    kaa.Update(1.0, id4, 0.0);
    // CCARAT: so3_mv6_m_upddydscl(fcta, bacstsn, bacstsn, kaa);
    kaa.MultiplyNT(fcta,backstresscurr,backstresscurr,1.0);
  }

}  // CalcBEBackStressFlow()


/*----------------------------------------------------------------------*
 | Get temperature-dependent material parameter at           dano 11/11 |
 | current temperature --> polynomial type                              |
 *----------------------------------------------------------------------*/
void MAT::Robinson::GetMatParameterAtTempnp(
  const vector<double>* paramvector, //  (i) given parameter is a vector
  double tempnp,  // tmpr (i) current temperature
  double parambytempnp // prmbytempr (o) parameter at current temperature
  )
{
  // polynomial type

  // initialise the temperature dependent material parameter
  parambytempnp = 0.0;
  double tempnp_pow = 1.0;

  // Param = a + b . T + c . T^2 + d . T^3 + ...
  // with T: current temperature
  for (unsigned i=0; i<(*paramvector).size(); ++i)
  {
    // calculate coefficient of variable T^i
    parambytempnp += (*paramvector)[i] * tempnp_pow;
    // for the higher polynom increase the exponent of the temperature
    tempnp_pow *= tempnp;
  }

}  // GetMatParameterAtTempnp()


/*----------------------------------------------------------------------*
 | Get temperature-dependent material parameter at           dano 11/11 |
 | current temperature                                                  |
 *----------------------------------------------------------------------*/
void MAT::Robinson::GetMatParameterAtTempnp(
  const double paramconst, // (i) given parameter is a constant
  double tempnp,  // tmpr (i) current temperature
  double parambytempnp // prmbytempr (o) parameter at current temperature
  )
{
  // initialise the temperature dependent material parameter
  parambytempnp = 0.0;

  // constant
  if (paramconst != 0.0)
  {
    // now calculate the parameter
    parambytempnp = paramconst;
  }

}  // GetMatParameterAtTempnp()


/*----------------------------------------------------------------------*
 | Reduce (statically condense) system in strain, strain_p,  dano 01/12 |
 | backstress to purely strain                                          |
 *----------------------------------------------------------------------*/
void MAT::Robinson::CalculateCondensedSystem(
  LINALG::Matrix<NUM_STRESS_3D,1>& stress,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kev,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kea,
  const LINALG::Matrix<NUM_STRESS_3D,1>& strain_pres,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kve,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kvv,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kva,
  const LINALG::Matrix<NUM_STRESS_3D,1>& backstress_res,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kae,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kav,
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& kaa,
  LINALG::Matrix<(2*NUM_STRESS_3D),1>& kvarva,  // (o)
  LINALG::Matrix<(2*NUM_STRESS_3D),NUM_STRESS_3D>& kvakvae  // (o)
  )
{
  // update vector for material internal variables (MIV) iterative increments
  // (stored as Fortranesque vector 2*NUMSTR_SOLID3x1)
  //             [ kvv  kva ]^{-1}   [ res^v  ]
  //    kvarva = [          ]      . [        ]
  //             [ kav  kaa ]      . [ res^al ]

  // update matrix for material internal variables (MIV)  iterative increments
  // (stored as Fortran-type (Fortranesque)vector (2*NUMSTR_SOLID3*NUMSTR_SOLID3)x1)
  //              [ kvv  kva ]^{-1}   [ kve ]
  //    kvakvae = [          ]      . [     ]
  //              [ kav  kaa ]        [ kae ]

  // build the matrix kvvkvakavkaa, consisting of the four submatrices, each
  // with size (6x6) --> kvvkvakavkaa: (12x12)
  //                [ kvv  kva ]
  // kvvkvakavkaa = [          ] and its inverse after factorisation
  //                [ kav  kaa ]
  // CCARAT: DOUBLE* kvvvaavaa = (DOUBLE*) CCACALLOC(numstr_2*numstr_2, sizeof(DOUBLE));
  LINALG::Matrix<2*NUM_STRESS_3D,2*NUM_STRESS_3D> kvvkvakavkaa(true);

  // build the matrix kevea (6x12)
  // kevea = [ kev  kea ] stored in column-major order
  // CCARAT: DOUBLE* kevea = (DOUBLE*) CCACALLOC(numstr*numstr_2, sizeof(DOUBLE));
  LINALG::Matrix<NUM_STRESS_3D,2*NUM_STRESS_3D> kevea(true);

  // ------------------ build tangent and right hand side to reduce
  {
    // first NUM_STRESS_3D rows (i=1--6)
    for (size_t i=0; i<NUM_STRESS_3D; i++)
    {
      // residual vector (i=1--6,j=1)
      kvarva(i) = strain_pres(i);

      // first NUM_STRESS_3D columns
      for(size_t j=0; j<NUM_STRESS_3D; j++)
      {
        // tangent (i=1--6,j=1--6)
        kvvkvakavkaa(i,j) = kvv(i,j);
        // RHS (i=1--6,j=1--6)
        kvakvae(i,j) = kve(i,j);
        // intermediate matrix - position in column-major vector matrix (i=1--6,j=1--6)
        kevea(i,j) = kev(i,j);
      }

      // second NUM_STRESS_3D columns
      for(size_t j=0; j<NUM_STRESS_3D; j++)
      {
        // tangent (i=1--6,j=6--12)
        kvvkvakavkaa(i,(j+NUM_STRESS_3D)) = kva(i,j);
        // position in column-major vector matrix (i=1--6,j=6--12)
        kevea(i,(j+NUM_STRESS_3D)) = kva(i,j);
      }
    } //  rows: i=1--6

    // second NUM_STRESS_3D rows (i=6--12)
    for (size_t i=NUM_STRESS_3D; i<(2*NUM_STRESS_3D); i++)
    {
      // residual vector (i=7--12,j=1)
      // add residual of backstresses
      kvarva(6+i) = backstress_res(i);

      // first NUM_STRESS_3D columns
      for(size_t j=0; j<NUM_STRESS_3D; j++)
      {
        // tangent (i=6--12,j=1--6)
        kvvkvakavkaa(i,j) = kav(i,j);
        // RHS (i=6--12,j=1--6)
        kvakvae(i,j) = kae(i,j);
      }
      // second NUM_STRESS_3D columns
      for(size_t j=0; j<NUM_STRESS_3D; j++)
      {
        // tangent (i=6--12,j=6--12)
        kvvkvakavkaa(i,(j+NUM_STRESS_3D)) = kaa(i,j);
      }
    } //  rows: i=6--12

  }  // end build tangent and rhs

  // ------------------------------- factorise kvvkvakavkaa and solve
  // adapt solving at drt_utils_gder2.H implemented by Georg Bauer
  // drt_utils_gder2: LINALG::FixedSizeSerialDenseSolver<numderiv2,numderiv2,numnode> solver;
  // solve x = A^{-1} . b

  // --------------------------------- back substitution of residuals
  // pass the size of the matrix and the number of columns to the solver
  // solve x = A^{-1} . b, with: A=kvvkvakavkaa(i), x=kvarva(o),b=kvarva (i)
  //           [ kvv  kva ]^{-1} [ res^v  ]^i
  // kvarva =  [          ]      [        ]
  //           [ kav  kaa ]      [ res^al ]
  LINALG::FixedSizeSerialDenseSolver<(2*NUM_STRESS_3D),(2*NUM_STRESS_3D),1> solver_res;
  solver_res.SetMatrix(kvvkvakavkaa);
  // No need for a separate rhs. We assemble the rhs to the solution vector.
  // The solver will destroy the rhs and return the solution.
  // solver_res.SetVectors(x,b); x: vector of unknowns, b: right hand side vector
  // kvarva = kvvkvakavkaa^{-1} . kvarva
  solver_res.SetVectors(kvarva,kvarva);
  solver_res.Solve();

  // ----------------------------------- back substitution of tangent
  // with: A=kvvkvakavkaa(i), x=kvakvae(o),b=kvakvae (i)
  //            [ kvv  kva ]^{-1} [ kve ]^i
  // kvakvae =  [          ]      [     ]
  //            [ kav  kaa ]      [ kae ]
  LINALG::FixedSizeSerialDenseSolver<(2*NUM_STRESS_3D),(2*NUM_STRESS_3D),NUM_STRESS_3D> solver_tang;
  solver_tang.SetMatrix(kvvkvakavkaa);
  // No need for a separate rhs. We assemble the rhs to the solution vector.
  // The solver will destroy the rhs and return the solution.
//  solver.SetVectors(derxy2,derxy2);
  // x: vector of unknowns, y: right hand side vector
  // kvakvae = kvvkvakavkaa^{-1} . kvakvae
  solver_tang.SetVectors(kvakvae,kvakvae);
  solver_tang.Solve();

  // TODO 18.01.12 CHECK if the result solver_* is correct

  // ----------------------------------------------------------------
  // reduce stress vector sigma_red
  {
    // stress row
    for (size_t i=0; i<NUM_STRESS_3D; i++)
    {
      double rcsum = 0.0;
      // column
      for (size_t j=0; j<(2*NUM_STRESS_3D); j++)
      {
        // kevea (6x12)
        // kvarva (12x1)
        rcsum += kevea(i,j) * kvarva(j);
      } // colums: j=1--12
      stress(i) -= rcsum;
    }  // stress rows: i=1--6
  }  // reduce stress vector sigma_red

  // ----------------------------------------------------------------
  // reduce tangent k_ee_red
  {
    for (size_t i=0; i<NUM_STRESS_3D; i++)
    {
      for (size_t j=0; j<NUM_STRESS_3D; j++)
      {
        double rcsum = 0.0;
        for (size_t k=0; k<(2*NUM_STRESS_3D); k++)
        {
          // kevea (6x12)
          // kvakvae (12x6)
          rcsum += kevea(i,k) * kvakvae(k,j);
        }
        cmat(i,j) -= rcsum;
      }
    }
  }  // reduce tangent k_ee_red

}  // CalculateCondensedSystem()


/*----------------------------------------------------------------------*
 | iterative update of material internal variables that      dano 01/12 |
 | are condensed out of the system within CalculateCondensedSystem()    |
 *----------------------------------------------------------------------*/
void MAT::Robinson::IterativeUpdateOfInternalVariables(
  const int gp,
  const LINALG::Matrix<NUM_STRESS_3D,1> epsii
  )
{
  // get the reduced residual
  //           [ kvv  kva ]^{-1} [ res^v  ]^i
  // kvarva =  [          ]      [        ]
  //           [ kav  kaa ]      [ res^al ]
  // CCARAT: DOUBLE* kvarva = actso3->miv_rob->kvarva.a.da[ip];
  LINALG::Matrix<2*NUM_STRESS_3D,1> kvarva(true);
  for (int i=0; i<12; i++)
    kvarva(i,0) = kvarva_->at(gp)(i,0);

  // get the scaled tangent
  //            [ kvv  kva ]^{-1} [ kve ]^i
  // kvakvae =  [          ]      [     ]
  //            [ kav  kaa ]      [ kae ]
  // CCARAT: DOUBLE* kvakvae = actso3->miv_rob->kvakvae.a.da[ip];
  LINALG::Matrix<2*NUM_STRESS_3D,NUM_STRESS_3D> kvakvae(true);
  for (size_t i=0; i<(2*NUM_STRESS_3D); i++)
    kvakvae(i,0) = kvakvae_->at(gp)(i,0);

  // get the new viscous strain
  // CCARAT: DOUBLE* vscstnn = actso3->miv_rob->vicstnn.a.da[ip];
  LINALG::Matrix<NUM_STRESS_3D,1> strain_pn(true);
  for (size_t i=0; i<NUM_STRESS_3D; i++)
    strain_pn(i,0) = strainplcurr_->at(gp)(i,0);

  // get the new back stresses
  // CCARAT: DOUBLE* bckstsn = actso3->miv_rob->bacstsn.a.da[ip];
  LINALG::Matrix<NUM_STRESS_3D,1> beta(true);
  for (size_t i=0; i<NUM_STRESS_3D; i++)
    beta(i,0) = backstresscurr_->at(gp)(i,0);

  // TODO 18.01.12 check if strain_pn used here is really the correct variable,
  // because ONLY the increment of the total strain is used as proposed by Burkhard!! epsii

  // ----------------------------------------------------------------
  // update viscous strain
  // [ iinc eps^v ] = [ kvv  kva ]^{-1} (   [ res^v  ] - [ kve ] [ iinc eps ] )
  // Delta strain_pn(i) = kvarva(i) - kvakvae(i) . Delta strain
  // with kvarva (12x1), kvakvae (12x6), Delta strain (6x1)
  {
    for (size_t i=0; i<NUM_STRESS_3D; i++)
    {
      // viscous residual contribution
      // Delta strain_pn(i) = [ kvv  kva ]^{-1} [ res^v  ] for i=1--6
      double rcsum = kvarva(i);

      // tangent contribution
      // Delta strain_pn(i)+ = [ kvv  kva ]^{-1} ( - [ kve ] [ iinc eps ] ) for i=1--6
      for (size_t j=0; j<NUM_STRESS_3D; j++)
      {
        // Delta strain_pn(i) += kvakvae(i,j) * epsii(j)
        // tangent kvakvae (12x6) (i=1--6,j=1--6)
        // increment of strain epsii (j=1--6)
        rcsum += kvakvae(i,j) * epsii(j);
      }
      // put all terms on strain_pn
      strain_pn(i) -= rcsum;
    }
  } // end update viscous strains

  // ----------------------------------------------------------------
  // update back stress
  {
    for (size_t i=0; i<NUM_STRESS_3D; i++)
    {
      // back stress residual contribution
      // CCARAT: DOUBLE rcsum = kvarva[NUM_STRESS_3D+istr];
      // Delta beta(i) = [ kav  kaa ]^{-1} [ res^v  ] for i=7--12
      double rcsum = kvarva(NUM_STRESS_3D+i);

      // tangent contribution
      // Delta beta(i)+ = [ kav  kaa ]^{-1} ( - [ kae ] [ iinc eps ] )
      //   for i=7--12, j=1--6
      for (size_t j=0; j<NUM_STRESS_3D; j++)
      {
        // Delta strain_pn(i) += kvakvae(i,j) * epsii(j)
        // tangent kvakvae (12x6) (i=7--12,j=1--6)
        // increment of strain epsii (j=1--6)
        rcsum += kvakvae(NUM_STRESS_3D+i,j) * epsii(j);
      }
      beta(i) -= rcsum;
    }
  } // end update back stress

  // TODO 29.11.11 check if we wanna update the history vectors
  // strain_p^{n+1} = strain_p^{n} + Delta strain_p^{n+1}
  // with Delta strain_pn(i) = [ kvv  kva ]^{-1} [ res^v  ]
  //                         + = [ kvv  kva ]^{-1} ( - [ kve ] [ iinc eps ] )
  //   for i=1--6, j=1--6
  strainplcurr_->at(gp) = strain_pn;
  // beta^{n+1} = beta^{n} + Delta beta^{n+1}
  // with Delta beta^{n+1} = [ kav  kaa ]^{-1} ( [ res^v  ]
  //                       + = [ kav  kaa ]^{-1} ( - [ kae ] [ iinc eps ] )
  //   for i=7--12, j=1--6
  backstresscurr_->at(gp) = beta;

}  // IterativeUpdateOfInternalVariables()


/*----------------------------------------------------------------------*
 | incremental update of Robinson's internal material        dano 01/12 |
 | variables (strain_p and backstress)                                  |
 | originally: so3_mat_robinson_be_mivupdincr()                         |
 *----------------------------------------------------------------------*/
void MAT::Robinson::IncrementalUpdateOfInternalVariables(
  const int gp  //!< total number of GPs in domain
  )
{
  // ----------------------------------------------------------------
  // update/reset material internal variables (MIVs)
  for (int i=0; i<gp; i++)
  {
    // reset visous strain mode
    // CCARAT: actso3->miv_rob->vscstns.a.iv[ip] = so3_mat_robinson_state_vague;
    // update viscous strain
    //    eps_{n}^v := eps_{n+1}^v at every Gauss point <g>
    // CCARAT: so3_mv6_v_ass(actso3->miv_rob->vicstnn.a.da[ip],
    //                       actso3->miv_rob->vicstn.a.da[ip]);
    for (size_t i=0; i<6; i++)
    {
      // strain_p = strain_pn
      // vicstn[i] = vicstnn[i];
      strainpllast_->at(gp)(i,0) = strainplcurr_->at(gp)(i,0);
    }

    // reset back stress mode
    // CCARAT: actso3->miv_rob->bckstss.a.iv[ip] = so3_mat_robinson_state_vague;
    // update back stress 'Alpha'
    //    al_{n} := al_{n+1} at every Gauss point <g>
    // CCARAT: so3_mv6_v_ass(actso3->miv_rob->bacstsn.a.da[ip],
    //                       actso3->miv_rob->bacsts.a.da[ip]);
    for (size_t i=0; i<6; i++)
    {
      // bacsts[i] = bacstsn[i];
      backstresslast_->at(gp)(i,0) = backstresscurr_->at(gp)(i,0);
    }

  } // gp

  return;
}  // IncrementalUpdateOfInternalVariables()


/*----------------------------------------------------------------------*
 | select Robinson's material and return stress              dano 12/11 |
 | originally: so3_mat_robinson_be_stress()                             |
 *----------------------------------------------------------------------*/
void MAT::Robinson::SelectRobinsonsMaterialBEReturnStress(
  const int gp, //!< number of Gauss points
  const LINALG::Matrix<NUM_STRESS_3D,1>& strain,  //!< total strain vector
  LINALG::Matrix<NUM_STRESS_3D,1>& stress,  //!< stress
                                            //!< DOUBLE stress[NUMSTR_SOLID3],
  LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>& cmat //!< material stiffness matrix
                                                     //!< DOUBLE cmat[NUMSTR_SOLID3][NUMSTR_SOLID3]
  )
{
  // elastic strain (CCARAT: elastn)
  LINALG::Matrix<NUM_STRESS_3D,1> strain_e(true);
  // thermal strain (CCARAT: thrstn)
  LINALG::Matrix<NUM_STRESS_3D,1> strain_t(true);
  // viscous strain (CCARAT: vscstn)
  LINALG::Matrix<NUM_STRESS_3D,1> strain_p(true);

  // temperature (CCARAT: tem)
  double tempnp;

  // ----------------------------------------------------------------
  // TODO 18.01.12 check difference between
  // CCARAT: gds->stnengv, L1577 and
  // CCARAT: gds->stnglv, L1589
  // total strain at t_{n}
  if (kintype_ == soh8_geolin) // CCARAT: so3_geo_lin
  {
    cout << "kintype_ "<< kintype_ << endl;
    // linear (engineering) strain vector
    // CCARAT: so3_mv6_v_ass(gds->stnengv, totstn);
    // assign gds->stnengv (in) by assigned vector totstn (output)
    //  strain{n} := al_{n+1} at every Gauss point <gds>
    for (int i=0; i<6; i++)
    {
      // totstn[i] = gds->stnengv[i];
//      strain(i,0) = strain;
    }
  }
  // TODO 18.01.12 check difference between
  else if (kintype_ == soh8_totlag) // CCARAT: so3_total_lagr
  {
    // Green-Lagrange strain vector
    // CCARAT: so3_mv6_v_ass(gds->stnglv, totstn);
    // totstn[i] = gds->stnglv[i];
    for (int i=0; i<6; i++)
    {
      // totstn[i] = gds->stnengv[i];
//      strain(i,0) = strain;
    }
  }
  else
  {
    dserror("Cannot digest chosen type of spatial kinematic\n");
  }

  // --------------------------------- temperature and thermal strain

  {
    // temperature at Gauss point
    // TODO 09.12.11 calculate strain due to temperature evaluation:
    // return current temperature solution tempnp, formerly: tem
    // CCARAT_ so3_tsi_temper(container, ele,gds->gpc[0], gds->gpc[1],
    //                        gds->gpc[2],&tem);

    // coefficient of linear thermal expansion
    // initialise the thermal expansion coefficient
    // CCARAT: DOUBLE thermexpans = mat_robin->thermexpans;
    const double thermexpans = params_->thermexpans_;

    // thermal strain vector (only volumetric strains)
    for (int i=0; i<3; i++)
      strain_t(i,0) = thermexpans * tempnp;
  }

  // ---------------------------------------- viscous strain at t_{n}
  // CCARAT: so3_mv6_v_ass(actso3->miv_rob->vicstn.a.da[ip], vscstn);

  // ---------------------------------------- elastic strain at t_{n}
  {
    // strain_e = strain - strain_p - strain_t
    // CCARAT: elastn[i] = totstn[i] - vscstn[i] - thrstn[i];
    strain_e.Update(1.0,strain,0.0);
    strain_e.Update(1.0,strain_p,1.0);
    strain_e.Update(1.0,strain_t,1.0);
  }

  // ---------------------------------------------- elasticity tensor
  // CCARAT: so3_mat_robinson_elmat(mat_robin, tem, cmat);
  SetupCmat(tempnp,cmat);

  // ------------------------------------------------ stress at t_{n}
  // index i: stress
  for (size_t i=0; i<6; i++)
  {
    // intermediate row * column sum
    double stresssum = 0.0;
    // index j: strain
    for (size_t j=0; j<6; j++)
    {
      stresssum += cmat(i,j) * strain_e(j);
    }
    stress(i) = stresssum;
  }

  return;

}  // SelectRobinsonsMaterialBEReturnStress()


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
