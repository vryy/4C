/*! \file
\brief

This routine implements a standard, i.e. local, crystal plasticity model.

See the header file for a detailed description.

\level 3

\maintainer Jan Schnabel
            jan-eike.schnabel@hzg.de
*/
/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/
#include "crystal_plasticity.H"
#include "matpar_bundle.H"
#include "material_service.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | constructor (public)                                      			|
 *----------------------------------------------------------------------*/
MAT::PAR::CrystalPlasticity::CrystalPlasticity(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      tol(matdata->GetDouble("TOL")),
      youngs(matdata->GetDouble("YOUNG")),
      poissonratio(matdata->GetDouble("NUE")),
      density(matdata->GetDouble("DENS")),
      latticetype(*matdata->Get<std::string>("LAT")),
      c_to_a_ratio(matdata->GetDouble("CTOA")),
      lattice_constant(matdata->GetDouble("ABASE")),
      num_slip_sys(matdata->GetInt("NUMSLIPSYS")),
      num_sub_sets(matdata->GetInt("NUMSUBSETS")),
      sub_set_members(*matdata->Get<std::vector<int>>("SUBSETMEMBERS")),
      rate_exp(*matdata->Get<std::vector<int>>("RATEEXP")),
      ref_shear_rate(*matdata->Get<std::vector<double>>("GAMMADOTREF")),
      dislocation_gen_coeff(*matdata->Get<std::vector<double>>("DISGENCOEFF")),
      dislocation_dyn_rec_coeff(*matdata->Get<std::vector<double>>("DISDYNRECCOEFF")),
      lattice_resistance(*matdata->Get<std::vector<double>>("TAUY0")),
      microstructural_boundaries(*matdata->Get<std::vector<double>>("GBS")),
      Hall_Petch_coeff(*matdata->Get<std::vector<double>>("HPCOEFF"))
{
}

/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       			|
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::CrystalPlasticity::CreateMaterial()
{
  return Teuchos::rcp(new MAT::CrystalPlasticity(this));
}

/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       			|
 *----------------------------------------------------------------------*/

MAT::CrystalPlasticityType MAT::CrystalPlasticityType::instance_;

/*----------------------------------------------------------------------*
 | is called in Material::Factory from ReadMaterials()       			|
 *----------------------------------------------------------------------*/
DRT::ParObject* MAT::CrystalPlasticityType::Create(const std::vector<char>& data)
{
  MAT::CrystalPlasticity* cp = new MAT::CrystalPlasticity();
  cp->Unpack(data);
  return cp;
}

/*----------------------------------------------------------------------*
 | constructor (public)                                                 |
 *----------------------------------------------------------------------*/
MAT::CrystalPlasticity::CrystalPlasticity() : params_(NULL) {}

/*----------------------------------------------------------------------*
 | copy-constructor (public)                                            |
 *----------------------------------------------------------------------*/
MAT::CrystalPlasticity::CrystalPlasticity(MAT::PAR::CrystalPlasticity* params) : params_(params) {}

/*----------------------------------------------------------------------*
 | pack (public)                                                        |
 *----------------------------------------------------------------------*/
void MAT::CrystalPlasticity::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // pack matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack history data
  int histsize;
  // if material is not yet initialised, i.e. at start of simulation, nothing to pack
  if (!(isinit and (F_curr != Teuchos::null)))
  {
    histsize = 0;
  }
  else
  {
    // if material is initialized already (restart): size equates number of Gauss points
    histsize = F_last->size();
  }
  AddtoPack(data, histsize);  // length of history vector(s)

  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, F_last->at(var));
    AddtoPack(data, FP_last->at(var));
    AddtoPack(data, gamma_last->at(var));
    AddtoPack(data, def_dens_last->at(var));
  }

}  // Pack()

/*----------------------------------------------------------------------*
 | unpack (public)                                                      |
 *----------------------------------------------------------------------*/
void MAT::CrystalPlasticity::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // recover matid and params_
  int matid;
  ExtractfromPack(position, data, matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const unsigned int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::CrystalPlasticity*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  // history data
  int histsize;
  ExtractfromPack(position, data, histsize);

  // if system is not yet initialised, the history vectors have to be intialized
  if (histsize == 0) isinit = false;

  if (params_ != NULL)
  {
    F_last = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
    FP_last = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
    gamma_last = Teuchos::rcp(new std::vector<std::vector<double>>);
    def_dens_last = Teuchos::rcp(new std::vector<std::vector<double>>);

    for (int var = 0; var < histsize; ++var)
    {
      LINALG::Matrix<3, 3> tmp_matrix(true);
      std::vector<double> tmp_vect(slip_sys_count);

      ExtractfromPack(position, data, tmp_matrix);
      F_last->push_back(tmp_matrix);

      ExtractfromPack(position, data, tmp_matrix);
      FP_last->push_back(tmp_matrix);

      ExtractfromPack(position, data, tmp_vect);
      gamma_last->push_back(tmp_vect);

      ExtractfromPack(position, data, tmp_vect);
      def_dens_last->push_back(tmp_vect);
    }

    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  }
}  // Unpack

/*---------------------------------------------------------------------*
 | initialize / allocate internal variables (public)                   |
 *---------------------------------------------------------------------*/
void MAT::CrystalPlasticity::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // import material / model parameters and calculate derived values

  // General Properties:
  //-------------------
  // tolerance for local Newton Raphson iteration
  Tol = params_->tol;

  // Elastic Properties:
  //-------------------------
  // Young's Modulus
  E = params_->youngs;
  // Poisson's ratio
  Poisson = params_->poissonratio;
  // 1st lamé constant
  Lambda = (Poisson * E) / ((1.0 + Poisson) * (1.0 - 2.0 * Poisson));
  // Shear modulus / 2nd Lame constant
  Mu = E / (2.0 * (1.0 + Poisson));
  // Bulk modulus
  K = E / (3.0 - (6.0 * Poisson));

  // Crystal Properties:
  //-------------------
  // setup lattice vectors
  // set up slip plane normals and directions and initialize crystal lattice related model
  // parameters according to chosen lattice type
  this->SetupLatticeVectors();

  //  read Lattice orientation matrix from .dat file
  this->SetupLatticeOrientation(linedef);

  // rotate lattice vectors according to lattice orientation
  for (int i = 0; i < slip_sys_count; i++)
  {
    LINALG::Matrix<3, 1> UnrotatedNor = SlipNor->at(i);
    LINALG::Matrix<3, 1> UnrotatedDir = SlipDir->at(i);

    SlipNor->at(i).MultiplyNN(*LatticeOrientation, UnrotatedNor);
    SlipDir->at(i).MultiplyNN(*LatticeOrientation, UnrotatedDir);
  }
  // TODO ADD TEST WHETHER INPUT ROTATION MATRIX IS ORTHOGONAL Q*Q^T=I TO AVOID UNWANTED SCALING OF
  // LATTICE VECTORS!

  // Index to which subset a slip system belongs
  SubSetIndex = params_->sub_set_members;

  // Viscoplastic Properties:
  //------------------------
  // TODO RATE EXPONENT INPUT
  // reference shear rates
  Gamma0 = params_->ref_shear_rate;

  // Dislocation Generation/Recovery:
  //--------------------------------
  // dislocation generation coefficients
  Dislocation_Gen_Coeff = params_->dislocation_gen_coeff;
  // dynamic dislocation removal coefficients
  Dislocation_Dyn_Rec_Coeff = params_->dislocation_dyn_rec_coeff;

  // Initial Slip System Strengths:
  //------------------------------
  // lattice resistances to slip
  TauY_0 = params_->lattice_resistance;
  // microstructural parameters which are relevant for Hall-Petch strengthening, e.g., grain size
  Micro_Boundary_Dist = params_->microstructural_boundaries;
  // Hall-Petch coefficients corresponding to above microstructural boundaries
  Hall_Petch_Coeffs = params_->Hall_Petch_coeff;



  // initialize history variables

  F_last = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  F_curr = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  FP_last = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  FP_curr = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  gamma_last = Teuchos::rcp(new std::vector<std::vector<double>>);
  gamma_curr = Teuchos::rcp(new std::vector<std::vector<double>>);

  def_dens_last = Teuchos::rcp(new std::vector<std::vector<double>>);
  def_dens_curr = Teuchos::rcp(new std::vector<std::vector<double>>);

  // set size to number of  Gauss points so each Gauss Point has its own set of history data
  F_last->resize(numgp);
  F_curr->resize(numgp);

  FP_last->resize(numgp);
  FP_curr->resize(numgp);

  gamma_last->resize(numgp);
  gamma_curr->resize(numgp);

  def_dens_last->resize(numgp);
  def_dens_curr->resize(numgp);

  // set initial values
  std::vector<double> emptyvect(slip_sys_count);
  for (int i = 0; i < slip_sys_count; i++) emptyvect.at(i) = 0.0;

  for (int i = 0; i < numgp; i++)
  {
    F_last->at(i) = Identity3();
    F_curr->at(i) = Identity3();

    FP_last->at(i) = Identity3();
    FP_curr->at(i) = Identity3();

    gamma_last->at(i).resize(slip_sys_count);
    gamma_curr->at(i).resize(slip_sys_count);

    gamma_last->at(i) = emptyvect;
    gamma_curr->at(i) = emptyvect;

    def_dens_last->at(i).resize(slip_sys_count);
    def_dens_curr->at(i).resize(slip_sys_count);

    for (int j = 0; j < slip_sys_count; j++)
      def_dens_last->at(i).at(j) = 1e5;  // TODO GET INITIAL DIS DENS ROM INPUT LINE

    def_dens_curr->at(i) = emptyvect;
  }

  isinit = true;
  return;
}  // Setup()


/*----------------------------------------------------------------------*
 | update internal variables (public)                                   |
 *----------------------------------------------------------------------*/
void MAT::CrystalPlasticity::Update()
{
  // update values of last step t_n to current values at time step t_n+1
  F_last = F_curr;
  FP_last = FP_curr;
  gamma_last = gamma_curr;
  def_dens_last = def_dens_curr;

  // empty vectors of current data
  F_curr = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  FP_curr = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  gamma_curr = Teuchos::rcp(new std::vector<std::vector<double>>);
  def_dens_curr = Teuchos::rcp(new std::vector<std::vector<double>>);

  // get the size of the vector (use the last vector, because it includes latest results, current is
  // already empty)
  const int histsize = F_last->size();
  F_curr->resize(histsize);
  FP_curr->resize(histsize);
  gamma_curr->resize(histsize);
  def_dens_curr->resize(histsize);

  LINALG::Matrix<3, 3> emptymat(true);
  std::vector<double> emptyvect(slip_sys_count);
  for (int i = 0; i < slip_sys_count; i++)
  {
    emptyvect.at(i) = 0.0;
  }

  for (int i = 0; i < histsize; i++)
  {
    F_curr->at(i) = emptymat;
    FP_curr->at(i) = emptymat;
    gamma_curr->at(i) = emptyvect;
    def_dens_curr->at(i) = emptyvect;
  }


  return;
}  // Update()

/*-------------------------------------------------------------------------------*
 | calculate stress, material stiffness matrix and evolution internal variables  |
 *-------------------------------------------------------------------------------*/
void MAT::CrystalPlasticity::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int eleGID)
{
  // simulation parameters
  //---------------------------------------------------------------
  // extract Gauss point number
  gp = params.get<int>("gp", -1);
  if (gp == -1) dserror("no Gauss point number provided in material");
  if (eleGID == -1) dserror("no element provided in material");

  // extract time increment
  dt = params.get<double>("delta time");

  // set current deformation gradient
  F_curr->at(gp) = *defgrd;


  // local Newton-Raphson to determine plastic shears
  //--------------------------------------------------------------------------

  // define result output from NR
  LINALG::Matrix<3, 3> PK2_res;

  // call Newton-Raphson method with current deformation gradient
  this->NewtonRaphson(
      F_curr->at(gp), gamma_curr->at(gp), def_dens_curr->at(gp), PK2_res, FP_curr->at(gp));

  // update stress results
  (*stress)(0) = PK2_res(0, 0);
  (*stress)(1) = PK2_res(1, 1);
  (*stress)(2) = PK2_res(2, 2);
  (*stress)(3) = PK2_res(0, 1);
  (*stress)(4) = PK2_res(1, 2);
  (*stress)(5) = PK2_res(0, 2);

  return;
}  // Evaluate()

/*---------------------------------------------------------------------------------*
 | Setup slip directions and slip plane normals                                    |
 *---------------------------------------------------------------------------------*/
void MAT::CrystalPlasticity::SetupLatticeVectors()
{
  // extract lattice type from user input
  LatticeType = params_->latticetype;

  // assign number of slip systems that corresponds to the given lattice
  if (LatticeType == "L10" || LatticeType == "D019")
  {
    slip_sys_count = 12;
  }
  else if (LatticeType == "TEST")
  {
    slip_sys_count = 1;
  }

  // check whether the number of slip systems given by the user coincides with the lattice type
  if (params_->num_slip_sys != slip_sys_count)
  {
    dserror(
        "Given number of slip systems NUMSLIPSYS = %d does not match the expected number for "
        "%s lattices!",
        params_->num_slip_sys, LatticeType.c_str());
  }

  // initialize slip plane normals and directions
  SlipNor = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 1>>);
  SlipDir = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 1>>);

  // initialize slip system ID vector
  SlipSysID = Teuchos::rcp(new std::vector<std::string>);

  // resize to number of slip systems
  SlipNor->resize(slip_sys_count);
  SlipDir->resize(slip_sys_count);
  SlipSysID->resize(slip_sys_count);

  // import c to a ratio of unit cell as provided by the user
  CtoA = params_->c_to_a_ratio;

  // set up corresponding set of slip planes depending on the chosen crystallographic lattice
  // L10 super lattice
  if (LatticeType == "L10")
  {
    // ordinary slip systems
    SlipSysID->at(0) = "ordinary_1";

    SlipNor->at(0)(0, 0) = 1.0;
    SlipNor->at(0)(1, 0) = 1.0;
    SlipNor->at(0)(2, 0) = 1.0 / CtoA;

    SlipDir->at(0)(0, 0) = 1.0 / 2.0;
    SlipDir->at(0)(1, 0) = -1.0 / 2.0;
    SlipDir->at(0)(2, 0) = 0.0 * CtoA;

    SlipSysID->at(1) = "ordinary_2";

    SlipNor->at(1)(0, 0) = -1.0;
    SlipNor->at(1)(1, 0) = -1.0;
    SlipNor->at(1)(2, 0) = 1.0 / CtoA;

    SlipDir->at(1)(0, 0) = 1.0 / 2.0;
    SlipDir->at(1)(1, 0) = -1.0 / 2.0;
    SlipDir->at(1)(2, 0) = 0.0 * CtoA;

    SlipSysID->at(2) = "ordinary_3";

    SlipNor->at(2)(0, 0) = 1.0;
    SlipNor->at(2)(1, 0) = -1.0;
    SlipNor->at(2)(2, 0) = 1.0 / CtoA;

    SlipDir->at(2)(0, 0) = 1.0 / 2.0;
    SlipDir->at(2)(1, 0) = 1.0 / 2.0;
    SlipDir->at(2)(2, 0) = 0.0 * CtoA;

    SlipSysID->at(3) = "ordinary_4";

    SlipNor->at(3)(0, 0) = -1.0;
    SlipNor->at(3)(1, 0) = 1.0;
    SlipNor->at(3)(2, 0) = 1.0 / CtoA;

    SlipDir->at(3)(0, 0) = 1.0 / 2.0;
    SlipDir->at(3)(1, 0) = 1.0 / 2.0;
    SlipDir->at(3)(2, 0) = 0.0 * CtoA;

    // super slip systems
    SlipSysID->at(4) = "super_1";

    SlipNor->at(4)(0, 0) = 1.0;
    SlipNor->at(4)(1, 0) = 1.0;
    SlipNor->at(4)(2, 0) = 1.0 / CtoA;

    SlipDir->at(4)(0, 0) = 0.0;
    SlipDir->at(4)(1, 0) = 1.0;
    SlipDir->at(4)(2, 0) = -1.0 * CtoA;

    SlipSysID->at(5) = "super_2";

    SlipNor->at(5)(0, 0) = 1.0;
    SlipNor->at(5)(1, 0) = 1.0;
    SlipNor->at(5)(2, 0) = 1.0 / CtoA;

    SlipDir->at(5)(0, 0) = 1.0;
    SlipDir->at(5)(1, 0) = 0.0;
    SlipDir->at(5)(2, 0) = -1.0 * CtoA;

    SlipSysID->at(6) = "super_3";

    SlipNor->at(6)(0, 0) = 1.0;
    SlipNor->at(6)(1, 0) = -1.0;
    SlipNor->at(6)(2, 0) = -1.0 / CtoA;

    SlipDir->at(6)(0, 0) = 0.0;
    SlipDir->at(6)(1, 0) = 1.0;
    SlipDir->at(6)(2, 0) = -1.0 * CtoA;

    SlipSysID->at(7) = "super_4";

    SlipNor->at(7)(0, 0) = 1.0;
    SlipNor->at(7)(1, 0) = -1.0;
    SlipNor->at(7)(2, 0) = 1.0 / CtoA;

    SlipDir->at(7)(0, 0) = 1.0;
    SlipDir->at(7)(1, 0) = 0.0;
    SlipDir->at(7)(2, 0) = -1.0 * CtoA;

    SlipSysID->at(8) = "super_5";

    SlipNor->at(8)(0, 0) = 1.0;
    SlipNor->at(8)(1, 0) = 1.0;
    SlipNor->at(8)(2, 0) = -1.0 / CtoA;

    SlipDir->at(8)(0, 0) = 0.0;
    SlipDir->at(8)(1, 0) = -1.0;
    SlipDir->at(8)(2, 0) = -1.0 * CtoA;

    SlipSysID->at(9) = "super_6";

    SlipNor->at(9)(0, 0) = 1.0;
    SlipNor->at(9)(1, 0) = 1.0;
    SlipNor->at(9)(2, 0) = -1.0 / CtoA;

    SlipDir->at(9)(0, 0) = -1.0;
    SlipDir->at(9)(1, 0) = 0.0;
    SlipDir->at(9)(2, 0) = -1.0 * CtoA;

    SlipSysID->at(10) = "super_7";

    SlipNor->at(10)(0, 0) = 1.0;
    SlipNor->at(10)(1, 0) = -1.0;
    SlipNor->at(10)(2, 0) = 1.0 / CtoA;

    SlipDir->at(10)(0, 0) = 0.0;
    SlipDir->at(10)(1, 0) = -1.0;
    SlipDir->at(10)(2, 0) = -1.0 * CtoA;

    SlipSysID->at(11) = "super_8";

    SlipNor->at(11)(0, 0) = 1.0;
    SlipNor->at(11)(1, 0) = -1.0;
    SlipNor->at(11)(2, 0) = -1.0 / CtoA;

    SlipDir->at(11)(0, 0) = -1.0;
    SlipDir->at(11)(1, 0) = 0.0;
    SlipDir->at(11)(2, 0) = -1.0 * CtoA;
  }
  // D019 super lattice
  else if (LatticeType == "D019")
  {
    // initialize slip plane normals and directions for Miller-Bravais indices
    std::vector<LINALG::Matrix<4, 1>> SlipNorHex(slip_sys_count);
    std::vector<LINALG::Matrix<4, 1>> SlipDirHex(slip_sys_count);

    // NOTE: the c to a ratio is incorporated during transformation from Miller-Bravais to Miller

    // prismatic slip systems
    SlipSysID->at(0) = "prismatic_1";

    SlipNorHex.at(0)(0, 0) = 1.0;
    SlipNorHex.at(0)(1, 0) = -1.0;
    SlipNorHex.at(0)(2, 0) = 0.0;
    SlipNorHex.at(0)(3, 0) = 0.0;

    SlipDirHex.at(0)(0, 0) = 1.0 / 3.0;
    SlipDirHex.at(0)(1, 0) = 1.0 / 3.0;
    SlipDirHex.at(0)(2, 0) = -2.0 / 3.0;
    SlipDirHex.at(0)(3, 0) = 0.0;

    SlipSysID->at(1) = "prismatic_2";

    SlipNorHex.at(1)(0, 0) = 0.0;
    SlipNorHex.at(1)(1, 0) = 1.0;
    SlipNorHex.at(1)(2, 0) = -1.0;
    SlipNorHex.at(1)(3, 0) = 0.0;

    SlipDirHex.at(1)(0, 0) = -2.0 / 3.0;
    SlipDirHex.at(1)(1, 0) = 1.0 / 3.0;
    SlipDirHex.at(1)(2, 0) = 1.0 / 3.0;
    SlipDirHex.at(1)(3, 0) = 0.0;

    SlipSysID->at(2) = "prismatic_3";

    SlipNorHex.at(2)(0, 0) = 1.0;
    SlipNorHex.at(2)(1, 0) = 0.0;
    SlipNorHex.at(2)(2, 0) = -1.0;
    SlipNorHex.at(2)(3, 0) = 0.0;

    SlipDirHex.at(2)(0, 0) = 1.0 / 3.0;
    SlipDirHex.at(2)(1, 0) = -2.0 / 3.0;
    SlipDirHex.at(2)(2, 0) = 1.0 / 3.0;
    SlipDirHex.at(2)(3, 0) = 0.0;

    // basal slip systems

    SlipSysID->at(3) = "basal_1";

    SlipNorHex.at(3)(0, 0) = 0.0;
    SlipNorHex.at(3)(1, 0) = 0.0;
    SlipNorHex.at(3)(2, 0) = 0.0;
    SlipNorHex.at(3)(3, 0) = 1.0;

    SlipDirHex.at(3)(0, 0) = 1.0 / 3.0;
    SlipDirHex.at(3)(1, 0) = 1.0 / 3.0;
    SlipDirHex.at(3)(2, 0) = -2.0 / 3.0;
    SlipDirHex.at(3)(3, 0) = 0.0;

    SlipSysID->at(4) = "basal_2";

    SlipNorHex.at(4)(0, 0) = 0.0;
    SlipNorHex.at(4)(1, 0) = 0.0;
    SlipNorHex.at(4)(2, 0) = 0.0;
    SlipNorHex.at(4)(3, 0) = 1.0;

    SlipDirHex.at(4)(0, 0) = -2.0 / 3.0;
    SlipDirHex.at(4)(1, 0) = 1.0 / 3.0;
    SlipDirHex.at(4)(2, 0) = 1.0 / 3.0;
    SlipDirHex.at(4)(3, 0) = 0.0;

    SlipSysID->at(5) = "basal_3";

    SlipNorHex.at(5)(0, 0) = 0.0;
    SlipNorHex.at(5)(1, 0) = 0.0;
    SlipNorHex.at(5)(2, 0) = 0.0;
    SlipNorHex.at(5)(3, 0) = 1.0;

    SlipDirHex.at(5)(0, 0) = 1.0 / 3.0;
    SlipDirHex.at(5)(1, 0) = -2.0 / 3.0;
    SlipDirHex.at(5)(2, 0) = 1.0 / 3.0;
    SlipDirHex.at(5)(3, 0) = 0.0;

    // pyramidal slip systems

    SlipSysID->at(6) = "pyramidal_1";

    SlipNorHex.at(6)(0, 0) = 1.0;
    SlipNorHex.at(6)(1, 0) = 1.0;
    SlipNorHex.at(6)(2, 0) = -2.0;
    SlipNorHex.at(6)(3, 0) = 1.0;

    SlipDirHex.at(6)(0, 0) = -1.0 / 3.0;
    SlipDirHex.at(6)(1, 0) = -1.0 / 3.0;
    SlipDirHex.at(6)(2, 0) = 2.0 / 3.0;
    SlipDirHex.at(6)(3, 0) = 6.0 / 3.0;

    SlipSysID->at(7) = "pyramidal_2";

    SlipNorHex.at(7)(0, 0) = 1.0;
    SlipNorHex.at(7)(1, 0) = -2.0;
    SlipNorHex.at(7)(2, 0) = 1.0;
    SlipNorHex.at(7)(3, 0) = 1.0;

    SlipDirHex.at(7)(0, 0) = -1.0 / 3.0;
    SlipDirHex.at(7)(1, 0) = 2.0 / 3.0;
    SlipDirHex.at(7)(2, 0) = -1.0 / 3.0;
    SlipDirHex.at(7)(3, 0) = 6.0 / 3.0;

    SlipSysID->at(8) = "pyramidal_3";

    SlipNorHex.at(8)(0, 0) = -2.0;
    SlipNorHex.at(8)(1, 0) = 1.0;
    SlipNorHex.at(8)(2, 0) = 1.0;
    SlipNorHex.at(8)(3, 0) = 1.0;

    SlipDirHex.at(8)(0, 0) = 2.0 / 3.0;
    SlipDirHex.at(8)(1, 0) = -1.0 / 3.0;
    SlipDirHex.at(8)(2, 0) = -1.0 / 3.0;
    SlipDirHex.at(8)(3, 0) = 6.0 / 3.0;

    SlipSysID->at(9) = "pyramidal_4";

    SlipNorHex.at(9)(0, 0) = -1.0;
    SlipNorHex.at(9)(1, 0) = -1.0;
    SlipNorHex.at(9)(2, 0) = 2.0;
    SlipNorHex.at(9)(3, 0) = 1.0;

    SlipDirHex.at(9)(0, 0) = 1.0 / 3.0;
    SlipDirHex.at(9)(1, 0) = 1.0 / 3.0;
    SlipDirHex.at(9)(2, 0) = -2.0 / 3.0;
    SlipDirHex.at(9)(3, 0) = 6.0 / 3.0;

    SlipSysID->at(10) = "pyramidal_5";

    SlipNorHex.at(10)(0, 0) = -1.0;
    SlipNorHex.at(10)(1, 0) = 2.0;
    SlipNorHex.at(10)(2, 0) = -1.0;
    SlipNorHex.at(10)(3, 0) = 1.0;

    SlipDirHex.at(10)(0, 0) = 1.0 / 3.0;
    SlipDirHex.at(10)(1, 0) = -2.0 / 3.0;
    SlipDirHex.at(10)(2, 0) = 1.0 / 3.0;
    SlipDirHex.at(10)(3, 0) = 6.0 / 3.0;

    SlipSysID->at(11) = "pyramidal_6";

    SlipNorHex.at(11)(0, 0) = 2.0;
    SlipNorHex.at(11)(1, 0) = -1.0;
    SlipNorHex.at(11)(2, 0) = -1.0;
    SlipNorHex.at(11)(3, 0) = 1.0;

    SlipDirHex.at(11)(0, 0) = -2.0 / 3.0;
    SlipDirHex.at(11)(1, 0) = 1.0 / 3.0;
    SlipDirHex.at(11)(2, 0) = 1.0 / 3.0;
    SlipDirHex.at(11)(3, 0) = 6.0 / 3.0;

    // transform Miller-Bravais index notation of hexagonal lattices to Miller index notation of
    // cubic lattices
    MAT::CrystalPlasticity::MillerBravaisToMiller(SlipNorHex, SlipDirHex);
  }
  else if (LatticeType == "HCP" or LatticeType == "BCC" or LatticeType == "FCC")
  {
    dserror("Sorry, %s lattices are not yet supported", LatticeType.c_str());
  }
  // academic test lattice containing only one slip system
  else if (LatticeType == "TEST")
  {
    SlipSysID->at(0) = "test_1";

    SlipNor->at(0)(0, 0) = 1.0;
    SlipNor->at(0)(1, 0) = 0.0;
    SlipNor->at(0)(2, 0) = 1.0;

    SlipDir->at(0)(0, 0) = 1.0;
    SlipDir->at(0)(1, 0) = 0.0;
    SlipDir->at(0)(2, 0) = -1.0;
  }
  else
  {
    dserror(
        "Lattice type not known. Please check the LAT parameter in the input file. Currently it "
        "has to be FCC, BCC, HCP, D019 or L10");
  }


  // test whether directions and normals are perpendicular
  std::vector<LINALG::Matrix<1, 1>> PerpTest;
  PerpTest.resize(slip_sys_count);

  for (int i = 0; i < slip_sys_count; i++)
  {
    PerpTest.at(i).MultiplyTN(SlipDir->at(i), SlipNor->at(i));
    if (PerpTest.at(i)(0, 0) != 0.0)
    {
      dserror(
          "Warning, slip direction and slip plane normal of slip system "
          "%d are not perpendicular!",
          i + 1);
    }
  }

  // determine magnitude of Burgers vector
  Burgers_Mag = Teuchos::rcp(new std::vector<double>);
  Burgers_Mag->resize(slip_sys_count);
  LINALG::Matrix<3, 1> tmp_scale;

  // import lattice constant from user input
  LatticeConstant = params_->lattice_constant;

  for (int i = 0; i < slip_sys_count; i++)
  {
    tmp_scale = SlipDir->at(i);
    tmp_scale.Scale(LatticeConstant);
    Burgers_Mag->at(i) = tmp_scale.Norm2();
  }

  // normalize slip plane normals
  for (int i = 0; i < slip_sys_count; i++)
  {
    SlipNor->at(i).Scale(1.0 / SlipNor->at(i).Norm2());
  }

  // normalize slip directions
  for (int i = 0; i < slip_sys_count; i++)
  {
    SlipDir->at(i).Scale(1.0 / SlipDir->at(i).Norm2());
  }
  return;
}


/*---------------------------------------------------------------------------------*
 | Read Lattice orientation matrix from .dat file                                    |
 *---------------------------------------------------------------------------------*/

void MAT::CrystalPlasticity::SetupLatticeOrientation(DRT::INPUT::LineDefinition* linedef)
{
  LatticeOrientation = Teuchos::rcp(new LINALG::Matrix<3, 3>);

  std::vector<double> fiber1;
  std::vector<double> fiber2;
  std::vector<double> fiber3;

  if (linedef->HaveNamed("FIBER1"))
  {
    // extract fiber vectors as columns of the rotation matrix
    linedef->ExtractDoubleVector("FIBER1", fiber1);
    linedef->ExtractDoubleVector("FIBER2", fiber2);
    linedef->ExtractDoubleVector("FIBER3", fiber3);

    // assemble rotation matrix
    for (int i = 0; i < 3; i++)
    {
      (*LatticeOrientation)(i, 0) = fiber1.at(i);
      (*LatticeOrientation)(i, 1) = fiber2.at(i);
      (*LatticeOrientation)(i, 2) = fiber3.at(i);
    }
  }
  // error
  else
  {
    dserror(
        "No lattice orientation matrix provided! Please add 'FIBER1', 'FIBER2' and 'FIBER3' as "
        "columns of the rotation matrix that relates the lattice orientation to the global "
        "coordinate system");
  }
  return;
}


/*---------------------------------------------------------------------------------*
 | Transform Miller Bravais to Miller index notation                               |
 *---------------------------------------------------------------------------------*/

void MAT::CrystalPlasticity::MillerBravaisToMiller(
    std::vector<LINALG::Matrix<4, 1>>& SlipNorHex, std::vector<LINALG::Matrix<4, 1>>& SlipDirHex)
{
  for (int unsigned i = 0; i < SlipNorHex.size(); i++)
  {
    SlipNor->at(i)(0, 0) = SlipNorHex.at(i)(0, 0);
    SlipNor->at(i)(1, 0) = (SlipNorHex.at(i)(0, 0) + 2.0 * SlipNorHex.at(i)(1, 0)) / sqrt(3.0);
    SlipNor->at(i)(2, 0) = SlipNorHex.at(i)(3, 0) / CtoA;

    SlipDir->at(i)(0, 0) = 1.5 * SlipDirHex.at(i)(0, 0);
    SlipDir->at(i)(1, 0) = sqrt(3.0) * (0.5 * SlipDirHex.at(i)(0, 0) + SlipDirHex.at(i)(1, 0));
    SlipDir->at(i)(2, 0) = SlipDirHex.at(i)(3, 0) * CtoA;
  }
  return;
}

/*---------------------------------------------------------------------------------*
 | local Newton-Raphson iteration                                                  |
 *---------------------------------------------------------------------------------*/
void MAT::CrystalPlasticity::NewtonRaphson(LINALG::Matrix<3, 3>& F, std::vector<double>& gamma_res,
    std::vector<double>& def_dens_res, LINALG::Matrix<3, 3>& PK2_res, LINALG::Matrix<3, 3>& FP_res)
{
  // initialize iteration
  // ----------------------------------
  // max. number of iterations
  const int itermax = 20;
  // iteration counter
  int itnum = 0;
  // convergence indicator
  bool converged = false;

  // elastic predictor
  //--------------------------------------------------------------------------
  // assume load step is elastic, i.e., there is no change in plastic shear
  std::vector<double> gamma_trial = gamma_last->at(gp);

  // trial values of internal variables
  std::vector<double> def_dens_trial(slip_sys_count);
  LINALG::Matrix<3, 3> PK2_trial;
  LINALG::Matrix<3, 3> FP_trial;

  // residuals
  std::vector<double> residuals_trial(slip_sys_count);
  double totalResidual;

  // empty matrix
  LINALG::Matrix<3, 3> emptymat(true);

  // iteration
  //--------------------------------------------------------------------------
  while (!converged && itnum++ < itermax)
  {
    // reset iteration relevant values
    totalResidual = 0.0;
    PK2_trial = emptymat;

    // set up flow rule with given deformation gradient F and trial plastic shears gamma_trial
    this->SetupFlowRule(F, gamma_trial, FP_trial, def_dens_trial, PK2_trial, residuals_trial);

    // determine total residuum as norm of the vector of slip system residuals
    for (int i = 0; i < slip_sys_count; i++)
    {
      totalResidual += pow(residuals_trial.at(i), 2.0);
    }
    totalResidual = sqrt(totalResidual);

    // convergence check
    if (totalResidual < Tol)
    {
      converged = true;

      // collect results
      //--------------------------------------------------------------------------
      gamma_res = gamma_trial;
      def_dens_res = def_dens_trial;
      PK2_res = PK2_trial;
      FP_res = FP_trial;
    }
    else
    {
      // update rule: gamma_trial(n+1) = gamma_trial(n) + d_gamma_trial(n)
      // therefore the Jacobean ResidualStiffness(j,i)=dr(j)/d_gamma_trial(i) is set up, where r is
      // the residuum solving ResidualStiffness * d_gamma_trial(n) = - r

      // for L10 or D019 lattices
      LINALG::Matrix<12, 12> ResidualStiffness_12x12;
      LINALG::Matrix<12, 1> d_gamma_trial_12(true);
      LINALG::Matrix<12, 1> residuals_trial_LIN_12;
      LINALG::FixedSizeSerialDenseSolver<12, 12, 1> NRSolver_12x12;
      // for TEST lattices
      LINALG::Matrix<1, 1> ResidualStiffness_1x1;
      LINALG::Matrix<1, 1> d_gamma_trial_1(true);
      LINALG::Matrix<1, 1> residuals_trial_LIN_1;
      LINALG::FixedSizeSerialDenseSolver<1, 1, 1> NRSolver_1x1;

      // differentiate Residuum with respect to gamma_trial by perturbation DEFAULT:  gammaEps =
      // 1.0e-9
      double gammaEps = 1.0e-9;

      // change of residuals.at(j) with perturbation of gamma_trial_.at(i)

      // perturb single shears
      for (int i = 0; i < slip_sys_count; i++)
      {
        // perturbed vector of plastic shears
        std::vector<double> gamma_pert(slip_sys_count);
        // resultant vector of residuals
        std::vector<double> residuals_pert(slip_sys_count);
        // resultant vector of defect densities
        std::vector<double> def_dens_pert(slip_sys_count);
        // resultant 2nd PK stress
        LINALG::Matrix<3, 3> PK2_pert(true);
        // resultant plastic part of deformation gradient
        LINALG::Matrix<3, 3> FP_pert;

        // initially unperturbed
        gamma_pert = gamma_trial;

        // perturbation of single gammas
        gamma_pert.at(i) += gammaEps;

        // evaluate flow rule for perturbed gamma
        this->SetupFlowRule(F, gamma_pert, FP_pert, def_dens_pert, PK2_pert, residuals_pert);

        // finite difference
        if (LatticeType == "L10" || LatticeType == "D019")
        {
          for (int j = 0; j < slip_sys_count; j++)
          {
            ResidualStiffness_12x12(j, i) =
                (residuals_pert.at(j) - residuals_trial.at(j)) / gammaEps;

            residuals_trial_LIN_12(i) = -residuals_trial.at(i);
          }
        }
        else if (LatticeType == "TEST")
        {
          for (int j = 0; j < slip_sys_count; j++)
          {
            ResidualStiffness_1x1(j, i) = (residuals_pert.at(j) - residuals_trial.at(j)) / gammaEps;

            residuals_trial_LIN_1(i) = -residuals_trial.at(i);
          }
        }
      }

      // solve resultant system of equations
      if (LatticeType == "L10" || LatticeType == "D019")
      {
        NRSolver_12x12.SetMatrix(ResidualStiffness_12x12);
        NRSolver_12x12.SetVectors(d_gamma_trial_12, residuals_trial_LIN_12);
        NRSolver_12x12.Solve();

        // update vector of plastic shears to new trial value
        for (int i = 0; i < slip_sys_count; i++)
        {
          gamma_trial.at(i) += d_gamma_trial_12(i);
        }
      }
      else if (LatticeType == "TEST")
      {
        NRSolver_1x1.SetMatrix(ResidualStiffness_1x1);
        NRSolver_1x1.SetVectors(d_gamma_trial_1, residuals_trial_LIN_1);
        NRSolver_1x1.Solve();

        // update vector of plastic shears to new trial value
        for (int i = 0; i < slip_sys_count; i++)
        {
          gamma_trial.at(i) += d_gamma_trial_1(i);
        }
      }
    }
  }  // while

  // check whether or not the result converged during itermax iterations.
  // TODO: If not: reduce time step and restart time increment!
  if (!converged)
  {
    dserror("Internal Newton Raphson did not converge within itermax iterations");
  }

  return;
}  // NewtonRaphson


/*---------------------------------------------------------------------------------*
 | Evaluate flow rule for a  given increment of plastic shear and set up residuum  |
 *---------------------------------------------------------------------------------*/

void MAT::CrystalPlasticity::SetupFlowRule(LINALG::Matrix<3, 3> F, std::vector<double> gamma_trial,
    LINALG::Matrix<3, 3>& FP_trial, std::vector<double>& def_dens_trial,
    LINALG::Matrix<3, 3>& PK2_trial, std::vector<double>& residuals_trial)
{
  // set up trial increment of plastic shear for each slip system
  std::vector<double> delta_gamma_trial(slip_sys_count);

  for (int i = 0; i < slip_sys_count; i++)
  {
    delta_gamma_trial.at(i) = gamma_trial.at(i) - gamma_last->at(gp).at(i);
  }

  // determine trial defect densities that would result from delta_gamma_trial
  //--------------------------------------------------------------------------

  double total_dislocation_density = 0.0;

  for (int i = 0; i < slip_sys_count; i++)
  {
    // dislocation generation
    //  TODO EQUATION
    def_dens_trial.at(i) =
        def_dens_last->at(gp).at(i) +
        (1.0 / (Burgers_Mag->at(i) * Dislocation_Gen_Coeff.at(SubSetIndex.at(i) - 1))) *
            sqrt(def_dens_last->at(gp).at(i)) * abs(delta_gamma_trial.at(i));

    // dynamic recovery
    def_dens_trial.at(i) -= Dislocation_Dyn_Rec_Coeff.at(SubSetIndex.at(i) - 1) *
                            def_dens_last->at(gp).at(i) * abs(delta_gamma_trial.at(i));

    // determine updated total defect densities
    total_dislocation_density += def_dens_trial.at(i);
  }

  // kinematics and plastic velocity gradient
  //--------------------------------------------------------------------------

  // determine trial plastic velocity gradient L_p
  // set up L_p_trial
  LINALG::Matrix<3, 3> LP_trial(true);
  LINALG::Matrix<3, 3> TempMat;

  for (int i = 0; i < slip_sys_count; i++)
  {
    TempMat.MultiplyNT(delta_gamma_trial.at(i), SlipDir->at(i), SlipNor->at(i));
    LP_trial.Update(LP_trial, TempMat);
  }

  // take unimodular part of I + L_p to ensure plastic incompressibility
  LINALG::Matrix<3, 3> Uni_IplusLP_trial(true);
  Uni_IplusLP_trial.Update(Identity3(), LP_trial);
  Uni_IplusLP_trial.Scale(pow(Uni_IplusLP_trial.Determinant(), -1.0 / 3.0));

  // determine trial plastic deformation gradient
  FP_trial.MultiplyNN(Uni_IplusLP_trial, FP_last->at(gp));

  // determine trial stress
  //--------------------------------------------------------------------------
  // determine trial elastic deformation gradient
  // get the inverse FP^{-1}
  LINALG::Matrix<3, 3> inv_FP_trial;
  inv_FP_trial.Invert(FP_trial);
  LINALG::Matrix<3, 3> FE_trial;
  FE_trial.MultiplyNN(F, inv_FP_trial);

  // calculate the Jacobi-determinant J = det(FE_{n+1}) and the logarithm of it
  double Jacobi_trial = FE_trial.Determinant();
  double ln_Jacobi_trial = log(Jacobi_trial);

  // set up elastic right cauchy green and its inverse
  LINALG::Matrix<3, 3> CE;
  CE.MultiplyTN(FE_trial, FE_trial);
  LINALG::Matrix<3, 3> inv_CE;
  inv_CE.Invert(CE);

  // 2nd Piola-Kirchhoff stress
  // S = lambda * ln_Jacobi_trial * inv_CE + mu * (Identity3() - inv_CE)

  PK2_trial.Update(Lambda * ln_Jacobi_trial, inv_CE, 1.0);
  PK2_trial.Update(Mu, Identity3(), 1.0);
  PK2_trial.Update(-Mu, inv_CE, 1.0);

  // Mandel stress
  LINALG::Matrix<3, 3> M;
  M.MultiplyNN(CE, PK2_trial);

  // setting up Residua with flow/creep rule
  //--------------------------------------------------------------------------

  for (int i = 0; i < slip_sys_count; i++)
  {
    // resolved shear stress/Schmid stress
    LINALG::Matrix<1, 1> RSS(true);

    LINALG::Matrix<3, 1> TempVec(true);

    TempVec.MultiplyNN(M, SlipNor->at(i));

    RSS.MultiplyTN(SlipDir->at(i), TempVec);

    // work hardening

    // lattice resistance
    double tau_Y0 = TauY_0.at(SubSetIndex.at(i) - 1);

    // Hall-Petch strengthening term
    double HP_strengthening = Hall_Petch_Coeffs.at(SubSetIndex.at(i) - 1) *
                              (1.0 / sqrt(Micro_Boundary_Dist.at(SubSetIndex.at(i) - 1)));

    // work hardening hardening increment due to accumulated dislocation density
    double delta_tau_Y = 0.5 * Mu * Burgers_Mag->at(i) * sqrt(total_dislocation_density);

    // slip system strength, i.e. tau_Y = tau_Y0 + HP_strengthening + delta_tau_Y
    double tau_Y = tau_Y0 + HP_strengthening + delta_tau_Y;
    // check for consitency
    if (tau_Y < 0.0) dserror("Negative slip systems strength! Please check your input");

    // shear rate as determined from the flow rule
    double gamma_dot = 0.0;

    gamma_dot = abs(RSS(0, 0)) / tau_Y;
    gamma_dot = pow(gamma_dot, 50.0);  // TODO EXPONENT AUS DAT FILE HOLEN
    gamma_dot = gamma_dot * Gamma0.at(SubSetIndex.at(i) - 1);
    gamma_dot = std::copysign(gamma_dot, RSS(0, 0));

    // set up corresponding residual
    residuals_trial.at(i) = delta_gamma_trial.at(i) / dt - gamma_dot;
  }
  return;
}  // FlowRuleSetup()


/*---------------------------------------------------------------------*
 | return names of visualization data (public)                         |
 *---------------------------------------------------------------------*/
void MAT::CrystalPlasticity::VisNames(std::map<std::string, int>& names)
{
  // temporary string variable for assembling the output names
  std::string ID;

  // plastic shears on single systems
  for (int sys = 0; sys < slip_sys_count; sys++)
  {
    ID = "plastic_shear_";
    ID += SlipSysID->at(sys);
    names[ID] = 1;  // scalar
  }

  // accumulated plastic shears
  ID = "accumulated_plastic_shear";
  names[ID] = 1;  // scalar

  // dislocation densities on single systems
  for (int sys = 0; sys < slip_sys_count; sys++)
  {
    ID = "dislocation_density_";
    ID += SlipSysID->at(sys);
    names[ID] = 1;  // scalar
  }

  // total dislocation density
  ID = "total_dislocation_density";
  names[ID] = 1;  // scalar
}  // VisNames()


/*---------------------------------------------------------------------*
 | return visualization data (public)                                  |
 *---------------------------------------------------------------------*/
bool MAT::CrystalPlasticity::VisData(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  // plastic shears gamma
  for (int sys = 0; sys < slip_sys_count; sys++)
  {
    if (name == "plastic_shear_" + SlipSysID->at(sys))
    {
      if ((int)data.size() != 1) dserror("size mismatch");
      for (int gp_index = 0; gp_index < numgp; gp_index++)
      {
        data[0] = gamma_last->at(gp_index).at(sys);
      }
    }
  }

  // accumulated plastic shear gamma_acc
  if (name == "accumulated_plastic_shear")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double gamma_acc;
    for (int gp_index = 0; gp_index < numgp; gp_index++)
    {
      gamma_acc = 0.0;
      for (int sys = 0; sys < slip_sys_count; sys++)
      {
        gamma_acc += abs(gamma_last->at(gp_index).at(sys));
      }
      data[0] = gamma_acc;
    }
  }

  // dislocation densities
  for (int sys = 0; sys < slip_sys_count; sys++)
  {
    if (name == "dislocation_density_" + SlipSysID->at(sys))
    {
      if ((int)data.size() != 1) dserror("size mismatch");
      for (int gp_index = 0; gp_index < numgp; gp_index++)
      {
        data[0] = def_dens_last->at(gp_index).at(sys);
      }
    }
  }

  // accumulated dislocation density
  if (name == "total_dislocation_density")
  {
    if ((int)data.size() != 1) dserror("size mismatch");
    double dis_dens_acc;
    for (int gp_index = 0; gp_index < numgp; gp_index++)
    {
      dis_dens_acc = 0.0;
      for (int sys = 0; sys < slip_sys_count; sys++)
      {
        dis_dens_acc += def_dens_last->at(gp_index).at(sys);
      }
      data[0] = dis_dens_acc;
    }
  }

  return true;
}  // VisData()
