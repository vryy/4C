/*----------------------------------------------------------------------*/
/*!
\brief Robinson's visco-plastic material

       example input line:
       MAT 1 MAT_Struct_Robinson  KIND Arya_NarloyZ  YOUNG POLY 2 1.47e9 -7.05e5
        NUE 0.34  DENS 8.89e-3  THEXPANS 0.0  INITTEMP 293.15
        HRDN_FACT 3.847e-12  HRDN_EXPO 4.0  SHRTHRSHLD POLY 2 69.88e8 -0.067e8
        RCVRY 6.083e-3  ACTV_ERGY 40000.0  ACTV_TMPR 811.0  G0 0.04  M_EXPO 4.365
        BETA POLY 3 0.8 0.0 0.533e-6  H_FACT 1.67e16

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

  In contrast to MAT_Struct_ThrStVenantK the calculation of a temperature-
  dependent stress does not occur. Here we have a thermal strain.

  As in the implementation of Burkhard, the displacement-dependent load term in
  the thermal equation is neglected so far.
  --> Future topic: extend to fully TSI.

\level 2

<pre>
\maintainer Matthias Mayr
</pre>
*/
/*----------------------------------------------------------------------*
 | Definitions                                               dano 11/11 |
 *----------------------------------------------------------------------*/

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
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/matpar_bundle.H"

#include "../drt_tsi/tsi_defines.H"

// include this header needed for KinematicType
#include "../drt_so3/so_hex8.H"

/*----------------------------------------------------------------------*
 | constructor (public)                                      dano 11/11 |
 *----------------------------------------------------------------------*/
MAT::PAR::Robinson::Robinson(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      kind_((matdata->Get<std::string>("KIND"))),
      youngs_(*(matdata->Get<std::vector<double>>("YOUNG"))),
      poissonratio_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS")),
      thermexpans_(matdata->GetDouble("THEXPANS")),
      inittemp_(matdata->GetDouble("INITTEMP")),
      hrdn_fact_(matdata->GetDouble("HRDN_FACT")),
      hrdn_expo_(matdata->GetDouble("HRDN_EXPO")),
      shrthrshld_(*(matdata->Get<std::vector<double>>("SHRTHRSHLD"))),
      rcvry_(matdata->GetDouble("RCVRY")),
      actv_ergy_(matdata->GetDouble("ACTV_ERGY")),
      actv_tmpr_(matdata->GetDouble("ACTV_TMPR")),
      g0_(matdata->GetDouble("G0")),
      m_(matdata->GetDouble("M_EXPO")),
      beta_(*(matdata->Get<std::vector<double>>("BETA"))),
      h_(matdata->GetDouble("H_FACT"))
{
}


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Robinson::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Robinson(this));
}


/*----------------------------------------------------------------------*
 |                                                           dano 02/12 |
 *----------------------------------------------------------------------*/
MAT::RobinsonType MAT::RobinsonType::instance_;


/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       dano 02/12 |
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::RobinsonType::Create(const std::vector<char>& data)
{
  MAT::Robinson* robinson = new MAT::Robinson();
  robinson->Unpack(data);
  return robinson;
}


/*----------------------------------------------------------------------*
 | constructor (public) --> called in Create()               dano 11/11 |
 *----------------------------------------------------------------------*/
MAT::Robinson::Robinson() : params_(NULL) {}


/*----------------------------------------------------------------------*
 | copy-constructor (public) --> called in CreateMaterial()  dano 11/11 |
 *----------------------------------------------------------------------*/
MAT::Robinson::Robinson(MAT::PAR::Robinson* params) : plastic_step(false), params_(params) {}


/*----------------------------------------------------------------------*
 | pack (public)                                             dano 11/11 |
 *----------------------------------------------------------------------*/
void MAT::Robinson::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  // in case we are in post-process mode
  if (params_ != NULL) matid = params_->Id();
  AddtoPack(data, matid);

  // pack history data
  int numgp;
  // if material is not initialised, i.e. start simulation, nothing to pack
  if (!Initialized())
  {
    numgp = 0;
  }
  else
  {
    // if material is initialised (restart): size equates number of gausspoints
    numgp = strainpllast_->size();
  }
  AddtoPack(data, numgp);  // Length of history vector(s)
  for (int var = 0; var < numgp; ++var)
  {
    // pack history data
    AddtoPack(data, strainpllast_->at(var));
    AddtoPack(data, backstresslast_->at(var));

    AddtoPack(data, kvarva_->at(var));
    AddtoPack(data, kvakvae_->at(var));
    AddtoPack(data, strain_last_.at(var));
  }

  return;

}  // Pack()


/*----------------------------------------------------------------------*
 | unpack (public)                                           dano 11/11 |
 *----------------------------------------------------------------------*/
void MAT::Robinson::Unpack(const std::vector<char>& data)
{
  isinit_ = true;
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Robinson*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // history data
  int numgp;
  ExtractfromPack(position, data, numgp);

  // if system is not yet initialised, the history vectors have to be intialized
  if (numgp == 0) isinit_ = false;

  // unpack strain vectors
  strainpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  // unpack back stress vectors (for kinematic hardening)
  backstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<NUM_STRESS_3D, 1>>);

  // unpack matrices needed for condensed system
  kvarva_ = Teuchos::rcp(new std::vector<LINALG::Matrix<2 * NUM_STRESS_3D, 1>>);
  kvakvae_ = Teuchos::rcp(new std::vector<LINALG::Matrix<2 * NUM_STRESS_3D, NUM_STRESS_3D>>);
  strain_last_.resize(0);

  for (int var = 0; var < numgp; ++var)
  {
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> tmp(true);
    LINALG::Matrix<2 * MAT::NUM_STRESS_3D, 1> tmp1(true);
    LINALG::Matrix<2 * MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> tmp2(true);

    // unpack strain/stress vectors of last converged state
    ExtractfromPack(position, data, tmp);
    strainpllast_->push_back(tmp);
    ExtractfromPack(position, data, tmp);
    backstresslast_->push_back(tmp);

    // unpack matrices of last converged state
    ExtractfromPack(position, data, tmp1);
    kvarva_->push_back(tmp1);
    ExtractfromPack(position, data, tmp2);
    kvakvae_->push_back(tmp2);

    ExtractfromPack(position, data, tmp);
    strain_last_.push_back(tmp);

    // current vectors have to be initialised
    strainplcurr_->push_back(tmp);
    backstresscurr_->push_back(tmp);
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;

}  // Unpack()


/*---------------------------------------------------------------------*
 | initialise / allocate internal stress variables (public) dano 11/11 |
 *---------------------------------------------------------------------*/
void MAT::Robinson::Setup(const int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // temporary variable for read-in
  std::string buffer;

  // initialise history variables
  strainpllast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>);
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>);

  backstresslast_ = Teuchos::rcp(new std::vector<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>);

  kvarva_ = Teuchos::rcp(new std::vector<LINALG::Matrix<2 * MAT::NUM_STRESS_3D, 1>>);
  kvakvae_ =
      Teuchos::rcp(new std::vector<LINALG::Matrix<2 * MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>>);

  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> emptymat(true);
  strainpllast_->resize(numgp);
  strainplcurr_->resize(numgp);
  strain_last_.resize(numgp, LINALG::Matrix<6, 1>(true));

  backstresslast_->resize(numgp);
  backstresscurr_->resize(numgp);

  LINALG::Matrix<2 * MAT::NUM_STRESS_3D, 1> emptymat2(true);
  kvarva_->resize(numgp);
  LINALG::Matrix<2 * MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> emptymat3(true);
  kvakvae_->resize(numgp);

  for (int i = 0; i < numgp; i++)
  {
    strainpllast_->at(i) = emptymat;
    strainplcurr_->at(i) = emptymat;

    backstresslast_->at(i) = emptymat;
    backstresscurr_->at(i) = emptymat;

    kvarva_->at(i) = emptymat2;
    kvakvae_->at(i) = emptymat3;
  }

  isinit_ = true;

  return;

}  // Setup()


/*---------------------------------------------------------------------*
 | update after time step                          (public) dano 11/11 |
 *---------------------------------------------------------------------*/
void MAT::Robinson::Update()
{
  // make current values at time step t_n+1 to values of last step t_n
  // x_n := x_n+1
  strainpllast_ = strainplcurr_;
  backstresslast_ = backstresscurr_;
  // the matrices do not have to be updated. They are reset after each time step

  // empty vectors of current data
  strainplcurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>);
  backstresscurr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>);

  // get size of the vector
  // (use the last vector, because it includes latest results, current is empty)
  const int numgp = strainpllast_->size();
  strainplcurr_->resize(numgp);
  backstresscurr_->resize(numgp);

  const LINALG::Matrix<MAT::NUM_STRESS_3D, 1> emptymat(true);
  const LINALG::Matrix<2 * MAT::NUM_STRESS_3D, 1> emptymat1(true);
  const LINALG::Matrix<2 * MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> emptymat2(true);

  for (int i = 0; i < numgp; i++)
  {
    strainplcurr_->at(i) = emptymat;
    backstresscurr_->at(i) = emptymat;
  }

  return;

}  // Update()


/*----------------------------------------------------------------------*
 | evaluate material (public)                                dano 11/11 |
 | select Robinson's material, integrate internal variables and return  |
 | stress and material tangent                                          |
 | CCARAT: so3_mat_robinson_be_sel() and so3_mat_robinson_be_stress     |
 *----------------------------------------------------------------------*/
void MAT::Robinson::Evaluate(const LINALG::Matrix<3, 3>* defgrd, const LINALG::Matrix<6, 1>* strain,
    Teuchos::ParameterList& params, LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat,
    const int eleGID)
{
  // extract from parameter list
  const int gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("No Gauss point number provided in Robinson material");
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> straininc(*strain);
  straininc.Update(-1., strain_last_[gp], 1.);
  strain_last_[gp] = *strain;
  const double scalartemp = params.get<double>("scalartemp", -1.0);
  if (scalartemp < 0.0) dserror("No temperature available in Robinson material");

  // update history of the condensed variables plastic strain and back stress
  // iterative update of the current history vectors at current Gauss point gp
  IterativeUpdateOfInternalVariables(gp, straininc);

  // name:                 BACI            --> CCARAT:
  // total strain:         strain          --> stntotn
  // elastic strain:       strain_e        --> stnela
  // thermal strain:       strain_t        --> stnthr
  // viscous strain:       strain_p        --> stnvscn
  // stress deviator:      devstress       --> devstsn
  // back stress:          backstress/beta --> bacsts/bacstsn
  // over/relative stress: eta             --> stsovrn

  // ----------------------------------------------------------------
  // implementation is identical for linear and Green-Lagrange strains
  // strains are calculated on element level and is passed to the material
  // --> material has no idea what strains are passed
  //     --> NO kintype is needed

  // get material parameters
  const double dt_ = params.get<double>("delta time");

  // build Cartesian identity 2-tensor I_{AB}
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;
  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  LINALG::Matrix<6, 6> id4(true);
  for (int i = 0; i < 6; i++) id4(i, i) = 1.0;

  // -------------------------------- temperatures and thermal strain
  // get initial and current temperature
  // initial temperature at Gauss point
  // CCARAT: tempr0 --> so3_tsi_temper0()
  const double tempinit = params_->inittemp_;
  // initialise the thermal expansion coefficient
  const double thermexpans = params_->thermexpans_;
  // thermal strain vector
  // CCARAT: stnthr
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strain_t(true);
  // update current temperature at Gauss point
  // CCARAT: tempr --> so3_tsi_temper()
  for (int i = 0; i < 3; ++i) strain_t(i) = thermexpans * (scalartemp - tempinit);
  // for (int i=3; i<6; ++i){ strain_t(i) = 2E_xy = 2E_yz = 2E_zx = 0.0; }

  // ------------------------------------------------- viscous strain
  // viscous strain strain_{n+1}^{v,i} at t_{n+1}
  // use the newest plastic strains here, i.e., from latest Newton iteration
  // (strainplcurr_), not necessarily equal to newest temporal strains (last_)
  // CCARAT: stnvscn
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strain_pn(false);
  strain_pn.Update(strainplcurr_->at(gp));
  // get history vector of old visco-plastic strain at t_n
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strain_p(false);
  strain_p.Update(strainpllast_->at(gp));

  // ------------------------------------------------- elastic strain
  // elastic strain at t_{n+1}
  // strain^{e}_{n+1}
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strain_e(false);

  // strain^e_{n+1} = strain_n+1 - strain^p_n - strain^t
  strain_e.Update(*strain);
  strain_e.Update((-1.0), strain_pn, (-1.0), strain_t, 1.0);

  // ---------------------------------------------- elasticity tensor
  // cmat = kee = pd(sig)/pd(eps)
  // CCARAT: so3_mat_robinson_elmat(mat_robin, tempr, cmat);
  // pass the current temperature to calculate the current youngs modulus
  SetupCmat(scalartemp, *cmat);

  // ------------------------------------ tangents of stress equation
  // declare single terms of elasto-plastic tangent Cmat_ep

  // kev = pd(sigma)/pd(eps^v)
  // tangent term resulting from linearisation \frac{\pd sig}{\pd eps^v}
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kev(true);
  // CCARAT: so3_mv6_m_assscl(-1.0, cmat, kev);
  // assign vector by another vector and scale it
  // (i): scale (-1.0)
  // (i): input matrix cmat
  // (o): output matrix kev
  kev.Update((-1.0), *cmat);

  // kea = pd(sigma)/pd(backstress)
  // tangent term resulting from linearisation \frac{\pd sig}{\pd al}
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kea(true);
  // initialise with 1 on the diagonals
  // CCARAT: so3_mv6_m_id(kea);
  kea.Update(id4);

  // ------------------------------------------------- elastic stress
  // stress sig_{n+1}^i at t_{n+1} */
  // (i): input matrix cmat
  // (i): input vector strain_e (CCARAT: stnela)
  // (o): output vector stress
  // CCARAT: so3_mv6_v_assmvp(cmat, stnela, stress);
  // stress_{n+1} = cmat . strain^e_{n+1}
  stress->MultiplyNN(*cmat, strain_e);

  // ------------------------------------------------------ devstress
  // deviatoric stress s_{n+1}^i at t_{n+1}
  // calculate the deviator from stress
  // CAUTION: s = 2G . devstrain only in case of small strain
  // CCARAT: so3_mv6_v_dev(stress, devstsn);
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> devstress(false);
  // trace of stress vector
  double tracestress = ((*stress)(0) + (*stress)(1) + (*stress)(2));
  for (int i = 0; i < 3; i++) devstress(i) = (*stress)(i)-tracestress / 3.0;
  for (int i = 3; i < MAT::NUM_STRESS_3D; i++) devstress(i) = (*stress)(i);
  // CAUTION: shear stresses (e.g., sigma_12)
  // in Voigt-notation the shear strains (e.g., strain_12) have to be scaled with
  // 1/2 normally considered in material tangent (using id4sharp, instead of id4)

  // ---------------------------------------------------- back stress
  // new back stress at t_{n+1} backstress_{n+1}^i
  // CCARAT: actso3->miv_rob->bacstsn.a.da[ip]
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> backstress_n(false);
  backstress_n.Update(backstresscurr_->at(gp));
  // old back stress at t_{n} backstress_{n}
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> backstress(false);
  backstress.Update(backstresslast_->at(gp));

  // ------------------------------------------ over/relativestress
  // overstress Sig_{n+1}^i = s_{n+1}^i - al_{n+1}^i
  // (i): input vector devstn
  // (i): input vector backstress_n
  // (o): output vector stsovr: subtract 2 vectors
  // CCARAT: so3_mv6_v_sub(devstsn, actso3->miv_rob->bacstsn.a.da[ip], stsovrn);
  // eta_{n+1} = devstress_{n+1} - backstress_{n+1}
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> eta(true);
  RelDevStress(devstress, backstress_n, eta);

  // to calculate the new history vectors (strainplcurr_, backstresscurr_), the
  // submatrices of the complete problems, that are condensed later, have to be calculated

  // ---------------------- residual of viscous strain, kve, kvv, kva
  // residual of visc. strain eps_{n+1}^<i> and its consistent tangent for <i>
  // CCARAT: so3_mat_robinson_be_rvscstn(ele, mat_robin, dt, tempr,
  //                                     actso3->miv_rob->vicstn.a.da[ip],
  //                                     actso3->miv_rob->vicstnn.a.da[ip],
  //                                     devstsn, stsovrn,
  //                                     &(actso3->miv_rob->vscstns.a.iv[ip]),
  //                                     vscstnr, kve, kvv, kva);
  // tangent term resulting from linearisation \frac{\pd res^v}{\pd eps}
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kve(true);
  // tangent term resulting from linearisation \frac{\pd res^v}{\pd eps^v}
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kvv(true);
  // tangent term resulting from linearisation \frac{\pd res^v}{\pd al}
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kva(true);
  // initialise the visco-plastic strain residual
  // CCARAT: strain_pres --> vscstnr;
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strain_pres(true);
  CalcBEViscousStrainRate(
      dt_, scalartemp, strain_p, strain_pn, devstress, eta, strain_pres, kve, kvv, kva);

  // ------------------------- residual of back stress, kae, kav, kaa
  // residual of back stress al_{n+1} and its consistent tangent
  // CCARAT: so3_mat_robinson_be_rbcksts(ele, mat_robin, dt, tempr,
  //                                     actso3->miv_rob->vicstn.a.da[ip],
  //                                     actso3->miv_rob->vicstnn.a.da[ip],
  //                                     devstsn,
  //                                     actso3->miv_rob->bacsts.a.da[ip],
  //                                     actso3->miv_rob->bacstsn.a.da[ip],
  //                                     &(actso3->miv_rob->bckstss.a.iv[ip]),
  //                                     bckstsr, kae, kav, kaa);
  // initialise the sub matrices needed for evaluation of the complete coupled
  // problem
  // tangent term resulting from linearisation \frac{\pd res^al}{\pd eps}
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kae(true);
  // tangent term resulting from linearisation \frac{\pd res^al}{\pd eps^v}
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kav(true);
  // tangent term resulting from linearisation \frac{\pd res^al}{\pd al}
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kaa(true);
  // initialise the back stress residual
  // back stress (residual): beta/backstress --> bckstsr
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> backstress_res(true);
  CalcBEBackStressFlow(dt_, scalartemp, strain_p, strain_pn, devstress, backstress, backstress_n,
      backstress_res, kae, kav, kaa);

  // build reduced system by condensing the evolution equations: only stress
  // equation remains
  // ------------------------------------- reduced stress and tangent
  // ==> static condensation
  // CCARAT: so3_mat_robinson_be_red(stress, cmat, kev, kea,
  //                          vscstnr, kve, kvv, kva,
  //                          bckstsr, kae, kav, kaa,
  //                          actso3->miv_rob->kvarva.a.da[ip],
  //                          actso3->miv_rob->kvakvae.a.da[ip]);
  CalculateCondensedSystem(*stress,  // (o): sigma_red --> used to calculate f_int^{e} at GP
      *cmat,                         // (o): C_red --> used to calculate stiffness k^{e} at GP
      kev, kea, strain_pres, kve, kvv, kva, backstress_res, kae, kav, kaa, (kvarva_->at(gp)),
      (kvakvae_->at(gp)));

  // pass the current plastic strains to the element (for visualisation)
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> plstrain(false);
  plstrain.Update(strainplcurr_->at(gp));
  // set in parameter list
  params.set<LINALG::Matrix<MAT::NUM_STRESS_3D, 1>>("plglstrain", plstrain);

}  // Evaluate()


/*----------------------------------------------------------------------*
 | computes isotropic elasticity tensor in matrix notion     dano 11/11 |
 | for 3d                                                               |
 *----------------------------------------------------------------------*/
void MAT::Robinson::SetupCmat(
    double tempnp, LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& cmat)
{
  // get material parameters
  // Young's modulus
  double emod = GetMatParameterAtTempnp(&(params_->youngs_), tempnp);
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
  const double mfac = emod / ((1.0 + nu) * (1.0 - 2.0 * nu));  // factor

  // clear the material tangent
  cmat.Clear();
  // write non-zero components --- axial
  cmat(0, 0) = mfac * (1.0 - nu);
  cmat(0, 1) = mfac * nu;
  cmat(0, 2) = mfac * nu;
  cmat(1, 0) = mfac * nu;
  cmat(1, 1) = mfac * (1.0 - nu);
  cmat(1, 2) = mfac * nu;
  cmat(2, 0) = mfac * nu;
  cmat(2, 1) = mfac * nu;
  cmat(2, 2) = mfac * (1.0 - nu);
  // write non-zero components --- shear
  cmat(3, 3) = mfac * 0.5 * (1.0 - 2.0 * nu);
  cmat(4, 4) = mfac * 0.5 * (1.0 - 2.0 * nu);
  cmat(5, 5) = mfac * 0.5 * (1.0 - 2.0 * nu);

}  // SetupCmat()


/*----------------------------------------------------------------------*
 | computes linear stress tensor                             dano 11/11 |
 *----------------------------------------------------------------------*/
void MAT::Robinson::Stress(const double p,                   //!< volumetric stress
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& devstress,  //!< deviatoric stress tensor
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress            //!< 2nd PK-stress
)
{
  // total stress = deviatoric + hydrostatic pressure . I
  // sigma = s + p . I
  stress.Update(devstress);
  for (int i = 0; i < 3; ++i) stress(i) += p;

}  // Stress()


/*----------------------------------------------------------------------*
 | compute relative deviatoric stress tensor                 dano 11/11 |
 *----------------------------------------------------------------------*/
void MAT::Robinson::RelDevStress(
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& devstress,     // (i) deviatoric stress tensor
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& backstress_n,  // (i) back stress tensor
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& eta                  // (o) relative stress
)
{
  // relative stress = deviatoric - back stress
  // eta = s - backstress
  eta.Update(1.0, devstress, (-1.0), backstress_n);

}  // RelDevStress()


/*----------------------------------------------------------------------*
 | residual of BE-discretised viscous strain rate            dano 11/11 |
 | at Gauss point                                                       |
 *----------------------------------------------------------------------*/
void MAT::Robinson::CalcBEViscousStrainRate(const double dt,  // (i) time step size
    double tempnp,                                            // (i) current temperature
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& strain_p,    // (i) viscous strain at t_n
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& strain_pn,   // (i) viscous strain at t_{n+1}^<i>
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& devstress,   // (i) stress deviator at t_{n+1}^<i>
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& eta,         // (i) over stress at t_{n+1}^<i>
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& strain_pres,       // (i,o) viscous strain residual
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kve,  // (i,o)
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kvv,  // (i,o)
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kva   // (i,o)
)
{
  // strain_pn' = f^v in a plastic load step
  // at the end of a time step, i.e. equilibrium, the residual has to vanish
  // r^v = ( strain_pn - strain_p ) / dt - f^v != 0
  //
  // discretising the flow equation strain_pn' with a backward Euler scheme
  // --> strain_pn' = ( strain_pn - strain_p ) / dt

  // initialise hardening exponent 'N' -->  nn
  double nn = params_->hrdn_expo_;

  // identity tensor in vector notation
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;

  // -------------------------------------------------- preliminaries
  // J2-invariant
  // J2 = 1/2 eta : eta
  // J2 = 1/2 (eta11^2 + eta22^2 + eta33^2 + 2 . eta12^2 + 2 . eta23^2 + 2 . eta13^2)
  // CCARAT: 1.) so3_mv6_v_dblctr(ovrstsn, ovrstsn, &(j2));
  //         2.) j2 *= 0.5;
  //         --> J_2 = 1/2 * Sig : Sig  with Sig...overstress
  double J2 = 0.0;
  J2 = 1.0 / 2.0 * (eta(0) * eta(0) + eta(1) * eta(1) + eta(2) * eta(2)) + +eta(3) * eta(3) +
       eta(4) * eta(4) + eta(5) * eta(5);

  // Bingham-Prager shear stress threshold at current temperature 'K^2'
  // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->shrthrshld, tmpr, &(kksq));
  double kksq = 0.0;
  kksq = GetMatParameterAtTempnp(&(params_->shrthrshld_), tempnp);

  // F = (J_2 - K^2) / K^2 = (J_2 / K^2) - 1
  double ff = 0.0;
  if (fabs(kksq) <= EPS10)  // shrthrshld = kksq
  {
    ff = -1.0;
    dserror("Division by zero: Shear threshold very close to zero");
  }
  else
  {
    // CCARAT: (J2 - kksq)/kksq;
    ff = (J2 - kksq) / kksq;
  }

  // hardening factor 'A' --> aa
  // calculate the temperature dependent material constant \bar{\mu} := aa
  double aa = 0.0;
  if (*(params_->kind_) == "Arya_CrMoSteel")
  {
    double mu = params_->hrdn_fact_;
    // calculate theta1 used for the material constant \bar{\mu}
    // \bar{\mu} = (23.8 . tempnp - 2635.0) . (1.0/811.0 - 1.0/tempnp)), vgl. (14)
    double th1 = (23.8 * tempnp - 2635.0) * (1.0 / 811.0 - 1.0 / tempnp);
    if (std::isinf(th1)) dserror("Infinite theta1");
    // theory is the same as in literature, but present implementation differs, e.g.,
    // here: A == \bar{\mu} = 0.5/(mu exp(theta1)) = 1/(2 mu exp(theta1)
    // cf. Arya: \bar{\mu} := \mu . exp(- theta1), vgl. (12), f(F) includes mu
    aa = 0.5 / (mu * std::exp(-th1));
  }
  else  // "Butler","Arya","Arya_NarloyZ"
  {
    aa = params_->hrdn_fact_;
  }

  // se = 1/2 * devstress : eta
  // with  eta: overstress, devstress:.deviat.stress
  // CCARAT: double contraction of 2vectors: so3_mv6_v_dblctr(ovrstsn, devstsn, &(ss));
  // CCARAT: ss = 1/2 * s : Sig
  double se = 0.0;
  se = 1.0 / 2.0 * (devstress(0) * eta(0) + devstress(1) * eta(1) + devstress(2) * eta(2)) +
       +devstress(3) * eta(3) + devstress(4) * eta(4) + devstress(5) * eta(5);

  // viscous strain rate
  // CCARAT: dvscstn[NUMSTR_SOLID3]
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strainrate_p(true);
  //-------------------------------------------------------------------
  // IF plastic step ( F > 0.0, (1/2 * devstress : eta) > 0.0 )
  //-------------------------------------------------------------------
  if ((ff > 0.0) and (se > 0.0))
  {
    //    std::cout << "plastic step: viscous strain rate strain_p'!= 0" << std::endl;

    // inelastic/viscous strain residual
    // epsilon_p' = F^n / sqrt(j2) . eta

    // fct = A . F^n / (J_2)^{1/2}
    // CCARAT: double fct = aa . ff^{nn} / sqrt(j2);
    double fct = aa * std::pow(ff, nn) / sqrt(J2);
    // calculate the viscous strain rate respecting that strain vector component
    // has a doubled shear component, but stress vectors not!
    // --> scale shear components accordingly
    // CCARAT: so3_mv6_v2_assscl(fct, ovrstsn, dvscstn);
    for (unsigned i = 0; i < 3; i++)
    {
      strainrate_p(i) = eta(i);
    }
    for (int i = 3; i < MAT::NUM_STRESS_3D; i++)
    {
      strainrate_p(i) = 2.0 * eta(i);
    }
    strainrate_p.Scale(fct);
  }  // inelastic
  //-------------------------------------------------------------------
  // ELSE IF elastic step (F <= 0.0, se <= 0.0)
  //-------------------------------------------------------------------
  else
  {
    // elastic step, no inelastic strains: strain_n^v' == 0
    // CCARAT: so3_mv6_v_zero(dvscstn);
    // --> dvstcstn == strainrate_p
    strainrate_p.Scale(0.0);
  }  // elastic

  //-------------------------------------------------------------------
  // --------------------- residual of viscous strain rate at t_{n+1}
  // res_{n+1}^v = (strain_{n+1}^v - strain_n^v)/dt - d_eps_{n+1}^v
  // strainrate_p corresponds to f^v(d_{n+1}, strain^p_{n+1}, alpha_{n+1}, T_{n+1})
  for (int i = 0; i < MAT::NUM_STRESS_3D; i++)
  {
    // strain_pres = 1/dt * (strain_{n+1}^v - strain_{n}^v - dt * strain_n^v')
    // CCARAT: vscstnr[istr] = (vscstnn[istr] - vscstn[istr] - dt*dvscstn[istr])/dt;

    // here is the essential!!!
    strain_pres(i) = (strain_pn(i) - strain_p(i) - dt * strainrate_p(i)) / dt;
  }


  //-------------------------------------------------------------------

  // --------------------------------- derivative of viscous residual
  // derivative of viscous residual with respect to over stress eta
  // kvs = d(strain_pres) / d (eta)
  // CCARAT: kvs = pd(res^v)/pd(Sig) -- d eps^v/d Sig
  // kvs[NUMSTR_SOLID3][NUMSTR_SOLID3];
  LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kvs(true);
  //-------------------------------------------------------------------
  // IF plastic step ( F > 0.0, (1/2 * devstress : eta) > 0.0 )
  //-------------------------------------------------------------------
  if ((ff > 0.0) and (se > 0.0))
  {
    // add facu to all diagonal terms matrix kvs
    // facu = - A . F^n / sqrt(J2)
    // CCARAT: so3_mv6_m_idscl(facu, kvs);
    double facu = -aa * std::pow(ff, nn) / sqrt(J2);
    for (int i = 0; i < MAT::NUM_STRESS_3D; i++) kvs(i, i) = facu;

    // update matrix with scaled dyadic vector product
    // CCARAT: so3_mv6_m_upddydscl(faco, ovrstsn, ovrstsn, kvs);
    // contribution: kvs = kvs + faco . (eta \otimes eta^T)
    // faco = -n . a . F^(n-1) / (kappa . sqrt(J2)) + a . F^n / (2. J2^{1.5})
    double faco = -nn * aa * std::pow(ff, (nn - 1.0)) / (kksq * sqrt(J2)) +
                  aa * std::pow(ff, nn) / (2.0 * std::pow(J2, 1.5));
    kvs.MultiplyNT(faco, eta, eta, 1.0);
    // multiply last 3 rows by 2 to conform with definition of strain vectors
    // CCARAT: so3_mv6_m2_updmtom2(kvs);
    // consider the difference between physical strains and Voigt notation
    for (int i = 3; i < MAT::NUM_STRESS_3D; i++)
    {
      for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
      {
        kvs(i, j) *= 2.0;
      }
    }  // rows 1--6
  }    // inelastic
  //-------------------------------------------------------------------
  // ELSE IF elastic step (F <= 0.0, 1/2 (devstress . eta) < 0.0)
  //-------------------------------------------------------------------
  else
  {
    // so3_mv6_m_zero(kvs);
    // in case of an elastic step, no contribution to kvs
    kvs.Scale(0.0);
  }

  //-------------------------------------------------------------------
  // - derivative of viscous residual with respect to total strain eps
  // CCARAT: kve = ( pd(res_{n+1}^v)/pd(eps_{n+1}) )|^<i>
  // kve = ( pd strain_pres^{n+1} )/ (pd strain^{n+1})|^<i>

  //-------------------------------------------------------------------
  // IF plastic step ( F > 0.0, (1/2 * devstress : eta) > 0.0 )
  //-------------------------------------------------------------------
  if ((ff > 0.0) and (se > 0.0))
  {
    // calculate elastic material tangent with temperature-dependent Young's modulus
    // CCARAT: so3_mat_robinson_elmat(mat_robin, tmpr, kse);
    // kse = (pd eta) / (pd strain)
    // CCARAT: kse = pd(Sig)/pd(eps) */
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> kse(true);
    // pass the current temperature to calculate the current youngs modulus
    SetupCmat(tempnp, kse);
    // Matrix vector product: cid2 = kse(i,j)*id2(j)
    // CCARAT: so3_mv6_v_assmvp(kse, iv, civ);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> cid2(false);
    cid2.Multiply(kse, id2);
    // update matrix with scaled dyadic vector product
    // CCARAT: so3_mv6_m_upddydscl(-1.0/3.0, iv, civ, kse);
    // contribution: kse = kse + (-1/3) . (id2 \otimes cid2^T)
    kse.MultiplyNT((-1.0 / 3.0), id2, cid2, 1.0);

    // assign kve by matrix-matrix product kvs . kse (inner product)
    // kve = kvs . kse
    // CCARAT: so3_mv6_m_mprd(kvs, kse, kve);
    // kve(i,j) = kvs(i,k).kse(k,j)
    kve.MultiplyNN(kvs, kse);
  }  // inelastic
  //-------------------------------------------------------------------
  // ELSE IF elastic step (F <= 0.0, (devstress . eta) <= 0.0)
  //-------------------------------------------------------------------
  else
  {
    // in case of an elastic step, no contribution to kvs
    // CCARAT: so3_mv6_m_zero(kve);
    kve.Scale(0.0);
  }  // elastic

  //-------------------------------------------------------------------
  // derivative of viscous residual with respect to viscous strain strain_p

  //-------------------------------------------------------------------
  // IF plastic step ( F > 0.0, (1/2 * devstress : eta) > 0.0 )
  //-------------------------------------------------------------------
  if ((ff > 0.0) and (se > 0.0))
  {
    // derivative ksv = (pd eta) / (pd strain_p)
    // CCARAT: ksv = pd(Sig)/pd(eps^v)
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D> ksv(true);
    // CCARAT: so3_mat_robinson_elmat(mat_robin, tmpr, ksv);
    // pass the current temperature to calculate the current youngs modulus
    SetupCmat(tempnp, ksv);

    // Matrix vector product: cid2 = kse(i,j)*id2(j)
    // CCARAT: so3_mv6_v_assmvp(ksv, iv, civ);
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1> cid2(false);
    cid2.Multiply(ksv, id2);
    // update matrix with scaled dyadic vector product
    // CCARAT: so3_mv6_m_upddydscl(-1.0/3.0, iv, civ, ksv);
    // contribution: ksv = ksv + (-1/3) . (id2 \otimes cid2^T)
    ksv.MultiplyNT((-1.0 / 3.0), id2, cid2, 1.0);

    // kvv = d(res^v)/d(esp^v) + pd(res^v)/pd(Sig) . pd(Sig)/pd(eps^v)
    // kvv = 1/dt * Id  +  kvs . ksv
    // scale matrix kvv with 1.0/dt
    // CCARAT: so3_mv6_m_idscl(1.0/dt, kvv);
    for (int i = 0; i < MAT::NUM_STRESS_3D; i++) kvv(i, i) = 1.0 / dt;
    // assign matrix kvv by scaled matrix-matrix product kvs, ksv (inner product)
    // kvv = kvv + (-1.0) . kvs . ksv;
    // CCARAT: so3_mv6_m_updmprdscl(-1.0, kvs, ksv, kvv);
    kvv.MultiplyNN((-1.0), kvs, ksv, 1.0);
  }  // inelastic
  //-------------------------------------------------------------------
  // ELSE IF elastic step (F <= 0.0, (1/2 * devstress : eta) <= 0.0)
  //-------------------------------------------------------------------
  else
  {
    // scale diagonal terms of kvv with 1.0/dt
    // CCARAT: so3_mv6_m_idscl(1.0/dt, kvv);
    for (int i = 0; i < MAT::NUM_STRESS_3D; i++) kvv(i, i) = 1.0 / dt;
  }  // elastic

  //-------------------------------------------------------------------
  // ----- derivative of viscous residual with respect to back stress
  // backstress --> alpha
  // kva = (pd res_{n+1}^v) / (pd back stress)
  // CCARAT: kva = pd(res^v)/pd(al)
  //-------------------------------------------------------------------
  // IF plastic step ( F > 0.0, (1/2 * devstress : eta) > 0.0 )
  //-------------------------------------------------------------------
  if ((ff > 0.0) and (se > 0.0))
  {
    // assign vector by another vector and scale it
    // kva = kvs . ksa = kva . (-Id)
    // with ksa = (pd eta) / (pd backstress) = - Id; eta = s - backstress
    // CCARAT: ksa = pd(Sig)/pd(al) = -Id
    kva.Update(-1.0, kvs);
  }  // inelastic
  //-------------------------------------------------------------------
  // ELSE IF elastic step (F <= 0.0, (1/2 * devstress : eta) <= 0.0)
  //-------------------------------------------------------------------
  else
  {
    for (int i = 0; i < 6; i++) kva(i, i) = 0.0;  // so3_mv6_m_zero(kva);
  }

}  // CalcBEViscousStrainRate()


/*----------------------------------------------------------------------*
 | residual of BE-discretised back stress and its            dano 11/11 |
 | consistent tangent according to the flow rule at Gauss point         |
 *----------------------------------------------------------------------*/
void MAT::Robinson::CalcBEBackStressFlow(const double dt, const double tempnp,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& strain_p,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& strain_pn,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& devstress,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& backstress,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& backstress_n,
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& backstress_res,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kae,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kav,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kaa)
{
  // backstress_n' = f^alpha in a plastic load step
  // at the end of a time step, i.e. equilibrium, the residual has to vanish
  // r^alpha = ( backstress_n - backstress ) / dt - f^alpha != 0
  //
  // discretising the flow equation backstress_n' with a backward Euler scheme
  // --> backstress_n' = ( backstress_n - backstress ) / dt

  // -------------------------------------------------- preliminaries
  // set Cartesian identity 4-tensor in 6-Voigt matrix notation
  // this is fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are stress-like 6-Voigt
  LINALG::Matrix<6, 6> id4sharp(true);
  for (unsigned i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (unsigned i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  // I_2 = 1/2 * Alpha : Alpha  with Alpha...back stress
  // CCARAT: so3_mv6_v_dblctr(backstresscurr, backstresscurr, &i2);
  // CCARAT: i2 *= 0.5;
  double i2 = 0.0;
  i2 = 1.0 / 2.0 *
           (backstress_n(0) * backstress_n(0) + backstress_n(1) * backstress_n(1) +
               backstress_n(2) * backstress_n(2)) +
       +backstress_n(3) * backstress_n(3) + +backstress_n(4) * backstress_n(4) +
       +backstress_n(5) * backstress_n(5);

  // Bingham-Prager shear stress threshold 'K_0^2' at activation temperature
  // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->shrthrshld, tem0, &(kk0sq));
  double kk0sq = 0.0;
  // activation temperature
  double tem0 = params_->actv_tmpr_;
  kk0sq = GetMatParameterAtTempnp(&(params_->shrthrshld_), tem0);

  // 'beta' at current temperature
  // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->beta, tmpr, &(beta));
  double beta = GetMatParameterAtTempnp(&(params_->beta_), tempnp);

  // 'H' at current temperature
  double hh = 0.0;
  // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->h, tmpr, &(hh));
  hh = GetMatParameterAtTempnp(params_->h_, tempnp);
  // CCARAT: if (mat_robin->kind == vp_robinson_kind_arya_narloyz)
  if (*(params_->kind_) == "Arya_NarloyZ")
  {
    hh *= std::pow(6.896, (1.0 + beta)) / (3.0 * kk0sq);
  }
  // CCARAT: else if (mat_robin->kind == vp_robinson_kind_arya_crmosteel)
  if (*(params_->kind_) == "Arya_CrMoSteel")
  {
    double mu = params_->hrdn_fact_;
    hh *= 2.0 * mu;
  }
  else  // (*(params_->kind_) == "Butler", "Arya", "Arya_NarloyZ")
  {
    // go on, no further changes for H required
  }

  // recovery/softening factor 'R_0'
  double rr0 = 0.0;
  // CCARAT: so3_mat_robinson_prmbytmpr(mat_robin->rcvry, tmpr, &(rr0));
  rr0 = GetMatParameterAtTempnp(params_->rcvry_, tempnp);
  // exponent 'm'
  double mm = params_->m_;
  // CCARAT: if (mat_robin->kind == vp_robinson_kind_arya_narloyz)
  if (*(params_->kind_) == "Arya_NarloyZ")
  {
    // pressure unit scale : cN/cm^2 = 10^-4 MPa
    const double pus = 1.0e-4;
    rr0 *= std::pow(6.896, (1.0 + beta + mm)) * std::pow(3.0 * kk0sq * pus * pus, (mm - beta));
  }
  else  // (*(params_->kind_) == "Butler", "Arya", "Arya_NarloyZ")
  {
    // go on, no further changes for R_0 required
  }

  // recovery/softening term 'R'
  // R = R_0 . exp[ Q_0 ( (T - Theta_0 )(T . Theta_0) )]
  // tem0: activation temperature
  double rr = 0.0;
  // activation energy 'Q_0'
  double q0 = params_->actv_ergy_;
  // T_{n+1} . T_0 < (1.0E-12)
  if (fabs(tempnp * tem0) <= EPS12)
  {
    // T_0 < (1.0E-12)
    if (fabs(tem0) <= EPS12)
    {
      rr = rr0;
    }
    // T_0 > (1.0E-12)
    else
    {
      rr = rr0 * exp(q0 / tem0);
    }
  }
  // T_{n+1} . T_0 > (1.0E-12)
  else
  {
    rr = rr0 * exp(q0 * (tempnp - tem0) / (tempnp * tem0));
    if (std::isinf(rr))
    {
      rr = rr0;
    }
  }

  // 'G_0'
  double gg0 = params_->g0_;

  // initialise 'G'
  // G = I_2/K_0^2
  double gg = 0.0;
  // K_0^2 < 1.0E-10
  if (fabs(kk0sq) <= EPS10)
  {
    gg = 0.0;
    dserror("Division by zero: Shear threshold very close to zero");
  }
  // K_0^2 > 1.0E-10
  else  // (fabs(kk0sq) > EPS10)
  {
    gg = sqrt(i2 / kk0sq);
  }

  // sa = 1/2 * devstress : backstresscurr
  // CCARAT: sa = 1/2 . Alpha . s, with Alpha:bacstsn, s: devstsn
  // CCARAT: so3_mv6_v_dblctr(bacstsn, devstsn, &(sa));
  // CCARAT: sa *= 0.5;
  double sa = 0.0;
  sa = 1.0 / 2.0 *
           (backstress_n(0) * devstress(0) + backstress_n(1) * devstress(1) +
               backstress_n(2) * devstress(2)) +
       +backstress_n(3) * devstress(3) + backstress_n(4) * devstress(4) +
       +backstress_n(5) * devstress(5);

  // ----------------- difference of current and last viscous strains
  // (Delta strain_p)_{n+1} = strain_pn - strain_p
  //  \incr \eps^v = \eps_{n+1}^v - \eps_{n}^v
  //  with halved entries to conform with stress vectors */
  // CCARAT: so3_mv6_v_sub(vscstnn, vscstn, vscstnd05);
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strain_pd05(false);
  strain_pd05.Update(1.0, strain_pn, (-1.0), strain_p);
  // Due to the fact that strain vector component have a doubled shear
  // components, i.e.
  //   strain = [ a11 a22 a33 | 2*a12 2*a23 2*a31 ]
  // but stress vectors have not, i.e.
  //   stress = [ a11 a22 a33 | a12 a23 a31 ]
  // we need to scale the last three entries.
  // CCARAT: so3_mv6_v05_updvtov05(vscstnd05);
  for (int i = 3; i < MAT::NUM_STRESS_3D; i++)
  {
    strain_pd05(i) *= 0.5;
  }


  // ------------------------------------------ residual of back stress
  //-------------------------------------------------------------------
  // IF plastic step (G > G_0, 1/2 (devstress . backstress) > 0.0)
  //-------------------------------------------------------------------
  if ((gg > gg0) and (sa > 0.0))
  {
    //     std::cout << "plastic step: back stress rate alpha'!= 0" << std::endl;

    double fctv = hh / std::pow(gg, beta);
    double fcta = rr * std::pow(gg, (mm - beta)) / sqrt(i2);
    for (int i = 0; i < MAT::NUM_STRESS_3D; i++)
    {
      // alpha_res = (alpha_{n+1} - alpha_{n})/dt - f^alpha)

      // here is the essential!!!
      backstress_res(i) =
          (backstress_n(i) - backstress(i) - fctv * strain_pd05(i) + dt * fcta * backstress_n(i)) /
          dt;
    }
  }  // inelastic
  //-------------------------------------------------------------------
  // ELSE IF elastic step (G <= G_0, 1/2 (devstress . backstress) <= 0.0)
  //-------------------------------------------------------------------
  else
  {
    double fctv = hh / std::pow(gg0, beta);
    double fcta = 0.0;
    if (sqrt(i2) < EPS10)
    {
      // sqrt(i2) := 1.0e6 assures units are OK
      fcta = rr * std::pow(gg0, (mm - beta)) / 1.0e6;
    }
    else  // (sqrt(i2) > EPS10)
    {
      fcta = rr * std::pow(gg0, (mm - beta)) / sqrt(i2);
    }

    for (int i = 0; i < MAT::NUM_STRESS_3D; i++)
    {
      // here is the essential!!!
      backstress_res(i) =
          (backstress_n(i) - backstress(i) - fctv * strain_pd05(i) + dt * fcta * backstress_n(i)) /
          dt;
    }
  }  // elastic


  // ----------------------------- derivative of back stress residual

  // ----------------------- derivative with respect to total strains
  // kae = pd(res^al)/pd(eps) == 0
  // CCARAT: if ( (*bckstss == so3_mat_robinson_state_elastic)
  //            || (*bckstss == so3_mat_robinson_state_inelastic) )
  //         so3_mv6_m_zero(kae);
  kae.Scale(0.0);

  // --------------------- derivative with respect to viscous strains
  // kav = pd(res_{n+1}^al)/pd(eps_{n+1}^v)
  //-------------------------------------------------------------------
  // IF plastic step (G > G_0, 1/2 (devstress . backstress) > 0.0)
  //-------------------------------------------------------------------
  if ((gg > gg0) and (sa > 0.0))
  {
    double fctv = -hh / (pow(gg, beta) * dt);
    // CCARAT: so3_mv6_m05_idscl(fctv, kav);
    kav.Update(fctv, id4sharp);
  }  // inelastic
  //-------------------------------------------------------------------
  // ELSE IF elastic step (G <= G_0, 1/2 (devstress . backstress) <= 0.0)
  //-------------------------------------------------------------------
  else
  {
    double fctv = -hh / (pow(gg0, beta) * dt);
    // CCARAT: so3_mv6_m05_idscl(fctv, kav);
    kav.Update(fctv, id4sharp);
  }  // elastic

  // ------------------------- derivative with respect to back stress
  // derivative of back stress residual with respect to back stress
  // kaa = pd(res_{n+1}^al)/pd(al_{n+1})
  // set Cartesian identity 4-tensor in 6x6-matrix notation (stress-like)
  // this is a 'mixed co- and contra-variant' identity 4-tensor, ie I^{AB}_{CD}
  // REMARK: rows are stress-like 6-Voigt
  //         columns are strain-like 6-Voigt
  LINALG::Matrix<6, 6> id4(true);
  for (int i = 0; i < 6; i++) id4(i, i) = 1.0;
  //-------------------------------------------------------------------
  // IF plastic step (G > G_0, 1/2 (devstress . backstress) > 0.0)
  //-------------------------------------------------------------------
  if ((gg > gg0) and (sa > 0.0))
  {
    double fctu = 1.0 / dt + rr * std::pow(gg, (mm - beta)) / sqrt(i2);
    double fctv = beta * hh / (std::pow(gg, (beta + 1.0)) * dt * kk0sq);
    double fcta = rr * (mm - beta) * std::pow(gg, (mm - beta - 1.0)) / (sqrt(i2) * kk0sq) -
                  rr * std::pow(gg, (mm - beta)) / (2.0 * std::pow(i2, 1.5));
    // CCARAT: so3_mv6_m_idscl(fctu, kaa);
    kaa.Update(fctu, id4);
    // CCARAT: so3_mv6_m_upddydscl(fctv, vscstnd05, bacstsn, kaa);
    kaa.MultiplyNT(fctv, strain_pd05, backstress_n, 1.0);
    // CCARAT: so3_mv6_m_upddydscl(fcta, bacstsn, bacstsn, kaa);
    kaa.MultiplyNT(fcta, backstress_n, backstress_n, 1.0);
  }  // inelastic
  //-------------------------------------------------------------------
  // ELSE IF elastic step (G <= G_0, 1/2 (devstress . backstress) <= 0.0)
  //-------------------------------------------------------------------
  else
  {
    double ii2;
    if (sqrt(i2) < EPS10)
    {
      ii2 = 1.0e12; /* sqrt(i2) := 1.0e6 assures units are OK */
    }
    else
    {
      ii2 = i2;
    }
    double fctu = 1.0 / dt + rr * std::pow(gg0, (mm - beta)) / sqrt(ii2);
    double fcta = -rr * std::pow(gg0, (mm - beta)) / (2.0 * std::pow(ii2, 1.5));
    // CCARAT: so3_mv6_m_idscl(fctu, kaa);
    kaa.Update(fctu, id4);
    // CCARAT: so3_mv6_m_upddydscl(fcta, bacstsn, bacstsn, kaa);
    kaa.MultiplyNT(fcta, backstress_n, backstress_n, 1.0);
  }  // elastic

}  // CalcBEBackStressFlow()


/*----------------------------------------------------------------------*
 | return temperature-dependent material parameter           dano 02/12 |
 | at current temperature --> polynomial type                           |
 *----------------------------------------------------------------------*/
double MAT::Robinson::GetMatParameterAtTempnp(
    const std::vector<double>* paramvector,  // (i) given parameter is a vector
    const double& tempnp                     // tmpr (i) current temperature
)
{
  // polynomial type

  // initialise the temperature dependent material parameter
  double parambytempnp = 0.0;
  double tempnp_pow = 1.0;

  // Param = a + b . T + c . T^2 + d . T^3 + ...
  // with T: current temperature
  for (unsigned int i = 0; i < (*paramvector).size(); ++i)
  {
    // calculate coefficient of variable T^i
    parambytempnp += (*paramvector)[i] * tempnp_pow;
    // for the higher polynom increase the exponent of the temperature
    tempnp_pow *= tempnp;
  }

  return parambytempnp;

}  // GetMatParameterAtTempnp()


/*----------------------------------------------------------------------*
 | Get temperature-dependent material parameter at           dano 11/11 |
 | current temperature                                                  |
 *----------------------------------------------------------------------*/
double MAT::Robinson::GetMatParameterAtTempnp(
    const double paramconst,  // (i) given parameter is a constant
    const double& tempnp      // tmpr (i) current temperature
)
{
  // initialise the temperature dependent material parameter
  double parambytempnp = 0.0;

  // constant
  if (paramconst != 0.0)
  {
    // now calculate the parameter
    parambytempnp = paramconst;
  }

  return parambytempnp;

}  // GetMatParameterAtTempnp()


/*----------------------------------------------------------------------*
 | Reduce (statically condense) system in strain, strain_p,  dano 01/12 |
 | backstress to purely strain                                          |
 *----------------------------------------------------------------------*/
void MAT::Robinson::CalculateCondensedSystem(
    LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress,  // (io): updated stress with condensed terms
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>&
        cmat,  // (io): updated tangent with condensed terms
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kev,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kea,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& strain_pres,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kve,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kvv,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kva,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& backstress_res,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kae,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kav,
    LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& kaa,
    LINALG::Matrix<(2 * MAT::NUM_STRESS_3D), 1>& kvarva,  // (io) solved/condensed residual
    LINALG::Matrix<(2 * MAT::NUM_STRESS_3D), MAT::NUM_STRESS_3D>&
        kvakvae  // (io) solved/condensed tangent
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
  LINALG::Matrix<2 * MAT::NUM_STRESS_3D, 2 * MAT::NUM_STRESS_3D> kvvkvakavkaa(true);

  // build the matrix kevea (6x12)
  // kevea = [ kev  kea ] stored in column-major order
  // CCARAT: DOUBLE* kevea = (DOUBLE*) CCACALLOC(numstr*numstr_2, sizeof(DOUBLE));
  LINALG::Matrix<MAT::NUM_STRESS_3D, 2 * MAT::NUM_STRESS_3D> kevea(true);

  // ------------------ build tangent and right hand side to reduce
  {
    // first MAT::NUM_STRESS_3D rows (i=1--6)
    for (int i = 0; i < MAT::NUM_STRESS_3D; i++)
    {
      // residual vector (i=1--6,j=1)
      kvarva(i) = strain_pres(i);

      // first MAT::NUM_STRESS_3D columns
      for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
      {
        // tangent (i=1--6,j=1--6)
        kvvkvakavkaa(i, j) = kvv(i, j);

        // RHS (i=1--6,j=1--6)
        kvakvae(i, j) = kve(i, j);
        // intermediate matrix - position in column-major vector matrix (i=1--6,j=1--6)
        kevea(i, j) = kev(i, j);
      }  // column: j=1--6

      // second MAT::NUM_STRESS_3D columns
      for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
      {
        // tangent (i=1--6,j=7--12)
        kvvkvakavkaa(i, (j + MAT::NUM_STRESS_3D)) = kva(i, j);
        // position in column-major vector matrix (i=1--6,j=6--12)
        kevea(i, (j + MAT::NUM_STRESS_3D)) = kva(i, j);
      }  // column: j=7--12
    }    //  rows: i=1--6

    // second MAT::NUM_STRESS_3D rows, i.e. (i=6--12)
    for (int i = 0; i < MAT::NUM_STRESS_3D; i++)
    {
      // residual vector (i=7--12,j=1)
      // add residual of backstresses
      kvarva(MAT::NUM_STRESS_3D + i, 0) = backstress_res(i);

      // first MAT::NUM_STRESS_3D columns
      for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
      {
        // tangent (i=6--12,j=1--6)
        kvvkvakavkaa(MAT::NUM_STRESS_3D + i, j) = kav(i, j);
        // RHS (i=6--12,j=1--6)
        kvakvae(MAT::NUM_STRESS_3D + i, j) = kae(i, j);
      }  // column: j=1--6
      // second MAT::NUM_STRESS_3D columns
      for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
      {
        // tangent (i=6--12,j=6--12)
        kvvkvakavkaa(MAT::NUM_STRESS_3D + i, (j + MAT::NUM_STRESS_3D)) = kaa(i, j);
      }  // column: j=7--12
    }    //  rows: i=6--12

  }  // end build tangent and rhs

  // ------------------------------- factorise kvvkvakavkaa and solve
  // adapt solving at drt_utils_gder2.H implemented by Georg Bauer
  // solve x = A^{-1} . b

  // --------------------------------- back substitution of residuals
  // pass the size of the matrix and the number of columns to the solver
  // solve x = A^{-1} . b
  // with: A=kvvkvakavkaa(i), x=kvarva(o), b=kvarva (i)
  //           [ kvv  kva ]^{-1} [ res^v  ]^i
  // kvarva =  [          ]      [        ]
  //           [ kav  kaa ]      [ res^al ]
  LINALG::FixedSizeSerialDenseSolver<(2 * MAT::NUM_STRESS_3D), (2 * MAT::NUM_STRESS_3D), 1>
      solver_res;
  solver_res.SetMatrix(kvvkvakavkaa);
  // No need for a separate rhs. We assemble the rhs to the solution vector.
  // The solver will destroy the rhs and return the solution.
  // x: vector of unknowns, b: right hand side vector
  // kvarva = kvvkvakavkaa^{-1} . kvarva
  solver_res.SetVectors(kvarva, kvarva);
  solver_res.Solve();

  // ----------------------------------- back substitution of tangent
  // solve x = A^{-1} . b
  // with: A=kvvkvakavkaa(i), x=kvakvae(o), b=kvakvae (i)
  //            [ kvv  kva ]^{-1} [ kve ]^i
  // kvakvae =  [          ]      [     ]
  //            [ kav  kaa ]      [ kae ]
  LINALG::FixedSizeSerialDenseSolver<(2 * MAT::NUM_STRESS_3D), (2 * MAT::NUM_STRESS_3D),
      MAT::NUM_STRESS_3D>
      solver_tang;
  solver_tang.SetMatrix(kvvkvakavkaa);
  // No need for a separate rhs. We assemble the rhs to the solution vector.
  // The solver will destroy the rhs and return the solution.
  // x: vector of unknowns, b: right hand side vector: here: x=b=kvakvae
  // kvakvae = kvvkvakavkaa^{-1} . kvakvae
  solver_tang.SetVectors(kvakvae, kvakvae);
  solver_tang.Solve();

  // final condensed system expressed only in stress, strain, cmat
  // sig_red^i = kee_red^i . iinc eps --> stress_red = cmat_red . Delta strain

  // --------------------------------- reduce stress vector sigma_red
  // stress (6x6) = kevea (6x12) . kvarva (12x1)
  stress.Multiply((-1.0), kevea, kvarva, 1.0);

  // ---------------------------------------- reduce tangent k_ee_red
  // cmat (6x6) = kevea (6x12) . kvakvae (12x6)
  cmat.MultiplyNN((-1.0), kevea, kvakvae, 1.0);

}  // CalculateCondensedSystem()


/*----------------------------------------------------------------------*
 | iterative update of material internal variables that      dano 01/12 |
 | are condensed out of the system within CalculateCondensedSystem()    |
 *----------------------------------------------------------------------*/
void MAT::Robinson::IterativeUpdateOfInternalVariables(
    const int gp, const LINALG::Matrix<MAT::NUM_STRESS_3D, 1> straininc  // total strain increment
)
{
  // get the condensed / reduced residual vector
  //           [ kvv  kva ]^{-1} [ res^v  ]^i
  // kvarva =  [          ]      [        ]
  //           [ kav  kaa ]      [ res^al ]
  // CCARAT: DOUBLE* kvarva = actso3->miv_rob->kvarva.a.da[ip];
  // condensed vector of residual
  LINALG::Matrix<(2 * MAT::NUM_STRESS_3D), 1> kvarva(false);
  kvarva = kvarva_->at(gp);

  // get the condensed / reduced scaled tangent
  //            [ kvv  kva ]^{-1} [ kve ]^i
  // kvakvae =  [          ]      [     ]
  //            [ kav  kaa ]      [ kae ]
  // CCARAT: DOUBLE* kvakvae = actso3->miv_rob->kvakvae.a.da[ip];
  LINALG::Matrix<(2 * MAT::NUM_STRESS_3D), MAT::NUM_STRESS_3D> kvakvae(false);
  kvakvae = kvakvae_->at(gp);

  // initialise updated vector with newest values strainp_np = strainp_n
  // strainp_np = strainp_n + Delta strain_np
  // viscous strain \f$\varepsilon_{n+1}\f$ at t_{n+1}^i
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> strain_pn(false);
  strain_pn.Update(strainplcurr_->at(gp));
  // back stress \f$alpha_{n+1}\f$ at t_{n+1}^i
  LINALG::Matrix<MAT::NUM_STRESS_3D, 1> backstress_n(false);
  backstress_n.Update(backstresscurr_->at(gp));

  // ---------------------------------  update current viscous strain
  // [ iinc eps^v ] = [ kvv  kva ]^{-1} (   [ res^v  ] - [ kve ] [ iinc eps ] )
  // Delta strain_pn(i) = kvarva(i) - kvakvae(i) . Delta strain
  // with kvarva (12x1), kvakvae (12x6), Delta strain (6x1)
  {
    for (int i = 0; i < MAT::NUM_STRESS_3D; i++)
    {
      double rcsum = 0.0;
      double rcsum1 = 0.0;

      // viscous residual contribution
      // Delta strain_pn(i) = [ kvv  kva ]^{-1} [ res^v  ] for i=1--6
      rcsum = kvarva(i);

      // tangent contribution
      // Delta strain_pn(i)+ = [ kvv  kva ]^{-1} ( - [ kve ] [ iinc eps ] ) for i=1--6
      for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
      {
        // Delta strain_pn(i) += kvakvae(i,j) * straininc(j)
        // tangent kvakvae (12x6) (i=1--6,j=1--6)
        // increment of strain (j=1--6)
        rcsum1 = kvakvae(i, j) * straininc(j);
      }
      // update all terms on strain_pn
      strain_pn(i) += -rcsum - rcsum1;
    }
  }  // end update viscous strains

  // ------------------------------------  update current back stress
  {
    for (int i = 0; i < MAT::NUM_STRESS_3D; i++)
    {
      double rcsum = 0.0;
      double rcsum1 = 0.0;

      // back stress residual contribution
      // CCARAT: DOUBLE rcsum = kvarva[MAT::NUM_STRESS_3D+istr];
      // Delta backstress(i) = [ kav  kaa ]^{-1} [ res^v  ] for i=7--12
      rcsum = kvarva(MAT::NUM_STRESS_3D + i);

      // tangent contribution
      // Delta backstress(i)+ = [ kav  kaa ]^{-1} ( - [ kae ] [ iinc eps ] )
      //   for i=7--12, j=1--6
      for (int j = 0; j < MAT::NUM_STRESS_3D; j++)
      {
        // Delta strain_pn(i) += kvakvae(i,j) * straininc(j)
        // tangent kvakvae (12x6) (i=7--12,j=1--6)
        // increment of strain (j=1--6)
        rcsum1 = kvakvae(MAT::NUM_STRESS_3D + i, j) * straininc(j);
      }

      backstress_n(i) += -rcsum - rcsum1;
    }
  }  // end update back stress

  // update the history vectors
  // strain_p^{n+1} = strain_p^{n} + Delta strain_p^{n+1}
  // with Delta strain_pn(i) = [ kvv  kva ]^{-1} [ res^v  ]
  //                         + = [ kvv  kva ]^{-1} ( - [ kve ] [ iinc eps ] )
  //   for i=1--6, j=1--6
  strainplcurr_->at(gp) = strain_pn;

  // backstress^{n+1} = backstress^{n} + Delta backstress^{n+1}
  // with Delta backstress^{n+1} = [ kav  kaa ]^{-1} ( [ res^v  ]
  //                               + [ kav  kaa ]^{-1} ( - [ kae ] [ iinc eps ] )
  //   for i=7--12, j=1--6
  backstresscurr_->at(gp) = backstress_n;

}  // IterativeUpdateOfInternalVariables()


/*----------------------------------------------------------------------*/
