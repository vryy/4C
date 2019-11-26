/*! \file
\brief

This routine implements a standard, i.e. local, crystal plasticity model.

See the header file for a detailed description.

\level 3

\maintainer Jan Schnabel
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
      tol_(matdata->GetDouble("TOL")),
      youngs_(matdata->GetDouble("YOUNG")),
      poissonratio_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS")),
      latticetype_(*matdata->Get<std::string>("LAT")),
      c_to_a_ratio_(matdata->GetDouble("CTOA")),
      lattice_constant_(matdata->GetDouble("ABASE")),
      num_slip_sys_(matdata->GetInt("NUMSLIPSYS")),
      num_sub_sets_(matdata->GetInt("NUMSUBSETS")),
      sub_set_members_(*matdata->Get<std::vector<int>>("SUBSETMEMBERS")),
      rate_exp_(*matdata->Get<std::vector<int>>("RATEEXP")),
      ref_shear_rate_(*matdata->Get<std::vector<double>>("GAMMADOTREF")),
      dislocation_gen_coeff_(*matdata->Get<std::vector<double>>("DISGENCOEFF")),
      dislocation_dyn_rec_coeff_(*matdata->Get<std::vector<double>>("DISDYNRECCOEFF")),
      lattice_resistance_(*matdata->Get<std::vector<double>>("TAUY0")),
      microstructural_boundaries_(*matdata->Get<std::vector<double>>("GBS")),
      Hall_Petch_coeff_(*matdata->Get<std::vector<double>>("HPCOEFF"))
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
MAT::CrystalPlasticity::CrystalPlasticity() : params_(nullptr) {}

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
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack history data
  int histsize;
  // if material is not yet initialised, i.e. at start of simulation, nothing to pack
  if (!(isinit_ and (F_curr_ != Teuchos::null)))
  {
    histsize = 0;
  }
  else
  {
    // if material is initialized already (restart): size equates number of Gauss points
    histsize = F_last_->size();
  }
  AddtoPack(data, histsize);  // length of history vector(s)

  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, (*F_last_)[var]);
    AddtoPack(data, (*FP_last_)[var]);
    AddtoPack(data, (*gamma_last_)[var]);
    AddtoPack(data, (*def_dens_last_)[var]);
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
  if (histsize == 0) isinit_ = false;

  if (params_ != nullptr)
  {
    F_last_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
    FP_last_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
    gamma_last_ = Teuchos::rcp(new std::vector<std::vector<double>>);
    def_dens_last_ = Teuchos::rcp(new std::vector<std::vector<double>>);

    for (int var = 0; var < histsize; ++var)
    {
      LINALG::Matrix<3, 3> tmp_matrix(true);
      std::vector<double> tmp_vect(slip_sys_count_);

      ExtractfromPack(position, data, tmp_matrix);
      F_last_->push_back(tmp_matrix);

      ExtractfromPack(position, data, tmp_matrix);
      FP_last_->push_back(tmp_matrix);

      ExtractfromPack(position, data, tmp_vect);
      gamma_last_->push_back(tmp_vect);

      ExtractfromPack(position, data, tmp_vect);
      def_dens_last_->push_back(tmp_vect);
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
  Tol_ = params_->tol_;

  // Elastic Properties:
  //-------------------------
  // Young's Modulus
  E_ = params_->youngs_;
  // Poisson's ratio
  Poisson_ = params_->poissonratio_;
  // 1st lamé constant
  Lambda_ = (Poisson_ * E_) / ((1.0 + Poisson_) * (1.0 - 2.0 * Poisson_));
  // Shear modulus / 2nd Lame constant
  Mu_ = E_ / (2.0 * (1.0 + Poisson_));
  // Bulk modulus
  K_ = E_ / (3.0 - (6.0 * Poisson_));

  // Crystal Properties:
  //-------------------
  // setup lattice vectors
  // set up slip plane normals and directions and initialize crystal lattice related model
  // parameters according to chosen lattice type
  this->SetupLatticeVectors();

  //  read Lattice orientation matrix from .dat file
  this->SetupLatticeOrientation(linedef);

  // rotate lattice vectors according to lattice orientation
  for (int i = 0; i < slip_sys_count_; i++)
  {
    LINALG::Matrix<3, 1> UnrotatedNor = SlipNor_[i];
    LINALG::Matrix<3, 1> UnrotatedDir = SlipDir_[i];

    SlipNor_[i].MultiplyNN(LatticeOrientation_, UnrotatedNor);
    SlipDir_[i].MultiplyNN(LatticeOrientation_, UnrotatedDir);
  }
  // TODO ADD TEST WHETHER INPUT ROTATION MATRIX IS ORTHOGONAL Q*Q^T=I TO AVOID UNWANTED SCALING OF
  // LATTICE VECTORS!

  // Index to which subset a slip system belongs
  SubSetIndex_ = params_->sub_set_members_;

  // Viscoplastic Properties:
  //------------------------
  // TODO RATE EXPONENT INPUT
  // reference shear rates
  Gamma0_ = params_->ref_shear_rate_;

  // Dislocation Generation/Recovery:
  //--------------------------------
  // dislocation generation coefficients
  Dislocation_Gen_Coeff_ = params_->dislocation_gen_coeff_;
  // dynamic dislocation removal coefficients
  Dislocation_Dyn_Rec_Coeff_ = params_->dislocation_dyn_rec_coeff_;

  // Initial Slip System Strengths:
  //------------------------------
  // lattice resistances to slip
  TauY_0_ = params_->lattice_resistance_;
  // microstructural parameters which are relevant for Hall-Petch strengthening, e.g., grain size
  Micro_Boundary_Dist_ = params_->microstructural_boundaries_;
  // Hall-Petch coefficients corresponding to above microstructural boundaries
  Hall_Petch_Coeffs_ = params_->Hall_Petch_coeff_;



  // initialize history variables

  F_last_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  F_curr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  FP_last_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  FP_curr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  gamma_last_ = Teuchos::rcp(new std::vector<std::vector<double>>);
  gamma_curr_ = Teuchos::rcp(new std::vector<std::vector<double>>);

  def_dens_last_ = Teuchos::rcp(new std::vector<std::vector<double>>);
  def_dens_curr_ = Teuchos::rcp(new std::vector<std::vector<double>>);

  // set size to number of  Gauss points so each Gauss Point has its own set of history data
  F_last_->resize(numgp);
  F_curr_->resize(numgp);

  FP_last_->resize(numgp);
  FP_curr_->resize(numgp);

  gamma_last_->resize(numgp);
  gamma_curr_->resize(numgp);

  def_dens_last_->resize(numgp);
  def_dens_curr_->resize(numgp);

  // set initial values
  std::vector<double> emptyvect(slip_sys_count_);
  for (int i = 0; i < slip_sys_count_; i++) emptyvect.at(i) = 0.0;

  for (int i = 0; i < numgp; i++)
  {
    (*F_last_)[i] = Identity3();
    (*F_curr_)[i] = Identity3();

    (*FP_last_)[i] = Identity3();
    (*FP_curr_)[i] = Identity3();

    (*gamma_last_)[i].resize(slip_sys_count_);
    (*gamma_curr_)[i].resize(slip_sys_count_);

    (*gamma_last_)[i] = emptyvect;
    (*gamma_curr_)[i] = emptyvect;

    (*def_dens_last_)[i].resize(slip_sys_count_);
    (*def_dens_curr_)[i].resize(slip_sys_count_);

    for (int j = 0; j < slip_sys_count_; j++)
      (*def_dens_last_)[i][j] = 1e5;  // TODO GET INITIAL DIS DENS ROM INPUT LINE

    (*def_dens_curr_)[i] = emptyvect;
  }

  isinit_ = true;
  return;
}  // Setup()


/*----------------------------------------------------------------------*
 | update internal variables (public)                                   |
 *----------------------------------------------------------------------*/
void MAT::CrystalPlasticity::Update()
{
  // update values of last step t_n to current values at time step t_n+1
  F_last_ = F_curr_;
  FP_last_ = FP_curr_;
  gamma_last_ = gamma_curr_;
  def_dens_last_ = def_dens_curr_;

  // empty vectors of current data
  F_curr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  FP_curr_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  gamma_curr_ = Teuchos::rcp(new std::vector<std::vector<double>>);
  def_dens_curr_ = Teuchos::rcp(new std::vector<std::vector<double>>);

  // get the size of the vector (use the last vector, because it includes latest results, current is
  // already empty)
  const int histsize = F_last_->size();
  F_curr_->resize(histsize);
  FP_curr_->resize(histsize);
  gamma_curr_->resize(histsize);
  def_dens_curr_->resize(histsize);

  LINALG::Matrix<3, 3> emptymat(true);
  std::vector<double> emptyvect(slip_sys_count_);
  for (int i = 0; i < slip_sys_count_; i++)
  {
    emptyvect.at(i) = 0.0;
  }

  for (int i = 0; i < histsize; i++)
  {
    (*F_curr_)[i] = emptymat;
    (*FP_curr_)[i] = emptymat;
    (*gamma_curr_)[i] = emptyvect;
    (*def_dens_curr_)[i] = emptyvect;
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
  (*F_curr_)[gp] = *defgrd;


  // local Newton-Raphson to determine plastic shears
  //--------------------------------------------------------------------------

  // define result output from NR
  LINALG::Matrix<3, 3> PK2_res;

  // call Newton-Raphson method with current deformation gradient
  this->NewtonRaphson(
      (*F_curr_)[gp], (*gamma_curr_)[gp], (*def_dens_curr_)[gp], PK2_res, (*FP_curr_)[gp]);

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
  LatticeType_ = params_->latticetype_;

  // assign number of slip systems that corresponds to the given lattice
  if (LatticeType_ == "L10" || LatticeType_ == "D019")
  {
    slip_sys_count_ = 12;
  }
  else if (LatticeType_ == "TEST")
  {
    slip_sys_count_ = 1;
  }

  // check whether the number of slip systems given by the user coincides with the lattice type
  if (params_->num_slip_sys_ != slip_sys_count_)
  {
    dserror(
        "Given number of slip systems NUMSLIPSYS = %d does not match the expected number for "
        "%s lattices!",
        params_->num_slip_sys_, LatticeType_.c_str());
  }

  // resize to number of slip systems
  SlipNor_.resize(slip_sys_count_);
  SlipDir_.resize(slip_sys_count_);
  SlipSysID_.resize(slip_sys_count_);

  // import c to a ratio of unit cell as provided by the user
  CtoA_ = params_->c_to_a_ratio_;

  // set up corresponding set of slip planes depending on the chosen crystallographic lattice
  // L10 super lattice
  if (LatticeType_ == "L10")
  {
    // ordinary slip systems
    SlipSysID_[0] = "ordinary_1";

    SlipNor_[0](0, 0) = 1.0;
    SlipNor_[0](1, 0) = 1.0;
    SlipNor_[0](2, 0) = 1.0 / CtoA_;

    SlipDir_[0](0, 0) = 1.0 / 2.0;
    SlipDir_[0](1, 0) = -1.0 / 2.0;
    SlipDir_[0](2, 0) = 0.0 * CtoA_;

    SlipSysID_[1] = "ordinary_2";

    SlipNor_[1](0, 0) = -1.0;
    SlipNor_[1](1, 0) = -1.0;
    SlipNor_[1](2, 0) = 1.0 / CtoA_;

    SlipDir_[1](0, 0) = 1.0 / 2.0;
    SlipDir_[1](1, 0) = -1.0 / 2.0;
    SlipDir_[1](2, 0) = 0.0 * CtoA_;

    SlipSysID_[2] = "ordinary_3";

    SlipNor_[2](0, 0) = 1.0;
    SlipNor_[2](1, 0) = -1.0;
    SlipNor_[2](2, 0) = 1.0 / CtoA_;

    SlipDir_[2](0, 0) = 1.0 / 2.0;
    SlipDir_[2](1, 0) = 1.0 / 2.0;
    SlipDir_[2](2, 0) = 0.0 * CtoA_;

    SlipSysID_[3] = "ordinary_4";

    SlipNor_[3](0, 0) = -1.0;
    SlipNor_[3](1, 0) = 1.0;
    SlipNor_[3](2, 0) = 1.0 / CtoA_;

    SlipDir_[3](0, 0) = 1.0 / 2.0;
    SlipDir_[3](1, 0) = 1.0 / 2.0;
    SlipDir_[3](2, 0) = 0.0 * CtoA_;

    // super slip systems
    SlipSysID_[4] = "super_1";

    SlipNor_[4](0, 0) = 1.0;
    SlipNor_[4](1, 0) = 1.0;
    SlipNor_[4](2, 0) = 1.0 / CtoA_;

    SlipDir_[4](0, 0) = 0.0;
    SlipDir_[4](1, 0) = 1.0;
    SlipDir_[4](2, 0) = -1.0 * CtoA_;

    SlipSysID_[5] = "super_2";

    SlipNor_[5](0, 0) = 1.0;
    SlipNor_[5](1, 0) = 1.0;
    SlipNor_[5](2, 0) = 1.0 / CtoA_;

    SlipDir_[5](0, 0) = 1.0;
    SlipDir_[5](1, 0) = 0.0;
    SlipDir_[5](2, 0) = -1.0 * CtoA_;

    SlipSysID_[6] = "super_3";

    SlipNor_[6](0, 0) = 1.0;
    SlipNor_[6](1, 0) = -1.0;
    SlipNor_[6](2, 0) = -1.0 / CtoA_;

    SlipDir_[6](0, 0) = 0.0;
    SlipDir_[6](1, 0) = 1.0;
    SlipDir_[6](2, 0) = -1.0 * CtoA_;

    SlipSysID_[7] = "super_4";

    SlipNor_[7](0, 0) = 1.0;
    SlipNor_[7](1, 0) = -1.0;
    SlipNor_[7](2, 0) = 1.0 / CtoA_;

    SlipDir_[7](0, 0) = 1.0;
    SlipDir_[7](1, 0) = 0.0;
    SlipDir_[7](2, 0) = -1.0 * CtoA_;

    SlipSysID_[8] = "super_5";

    SlipNor_[8](0, 0) = 1.0;
    SlipNor_[8](1, 0) = 1.0;
    SlipNor_[8](2, 0) = -1.0 / CtoA_;

    SlipDir_[8](0, 0) = 0.0;
    SlipDir_[8](1, 0) = -1.0;
    SlipDir_[8](2, 0) = -1.0 * CtoA_;

    SlipSysID_[9] = "super_6";

    SlipNor_[9](0, 0) = 1.0;
    SlipNor_[9](1, 0) = 1.0;
    SlipNor_[9](2, 0) = -1.0 / CtoA_;

    SlipDir_[9](0, 0) = -1.0;
    SlipDir_[9](1, 0) = 0.0;
    SlipDir_[9](2, 0) = -1.0 * CtoA_;

    SlipSysID_[10] = "super_7";

    SlipNor_[10](0, 0) = 1.0;
    SlipNor_[10](1, 0) = -1.0;
    SlipNor_[10](2, 0) = 1.0 / CtoA_;

    SlipDir_[10](0, 0) = 0.0;
    SlipDir_[10](1, 0) = -1.0;
    SlipDir_[10](2, 0) = -1.0 * CtoA_;

    SlipSysID_[11] = "super_8";

    SlipNor_[11](0, 0) = 1.0;
    SlipNor_[11](1, 0) = -1.0;
    SlipNor_[11](2, 0) = -1.0 / CtoA_;

    SlipDir_[11](0, 0) = -1.0;
    SlipDir_[11](1, 0) = 0.0;
    SlipDir_[11](2, 0) = -1.0 * CtoA_;
  }
  // D019 super lattice
  else if (LatticeType_ == "D019")
  {
    // initialize slip plane normals and directions for Miller-Bravais indices
    std::vector<LINALG::Matrix<4, 1>> SlipNorHex(slip_sys_count_);
    std::vector<LINALG::Matrix<4, 1>> SlipDirHex(slip_sys_count_);

    // NOTE: the c to a ratio is incorporated during transformation from Miller-Bravais to Miller

    // prismatic slip systems
    SlipSysID_[0] = "prismatic_1";

    SlipNorHex[0](0, 0) = 1.0;
    SlipNorHex[0](1, 0) = -1.0;
    SlipNorHex[0](2, 0) = 0.0;
    SlipNorHex[0](3, 0) = 0.0;

    SlipDirHex[0](0, 0) = 1.0 / 3.0;
    SlipDirHex[0](1, 0) = 1.0 / 3.0;
    SlipDirHex[0](2, 0) = -2.0 / 3.0;
    SlipDirHex[0](3, 0) = 0.0;

    SlipSysID_[1] = "prismatic_2";

    SlipNorHex[1](0, 0) = 0.0;
    SlipNorHex[1](1, 0) = 1.0;
    SlipNorHex[1](2, 0) = -1.0;
    SlipNorHex[1](3, 0) = 0.0;

    SlipDirHex[1](0, 0) = -2.0 / 3.0;
    SlipDirHex[1](1, 0) = 1.0 / 3.0;
    SlipDirHex[1](2, 0) = 1.0 / 3.0;
    SlipDirHex[1](3, 0) = 0.0;

    SlipSysID_[2] = "prismatic_3";

    SlipNorHex[2](0, 0) = 1.0;
    SlipNorHex[2](1, 0) = 0.0;
    SlipNorHex[2](2, 0) = -1.0;
    SlipNorHex[2](3, 0) = 0.0;

    SlipDirHex[2](0, 0) = 1.0 / 3.0;
    SlipDirHex[2](1, 0) = -2.0 / 3.0;
    SlipDirHex[2](2, 0) = 1.0 / 3.0;
    SlipDirHex[2](3, 0) = 0.0;

    // basal slip systems

    SlipSysID_[3] = "basal_1";

    SlipNorHex[3](0, 0) = 0.0;
    SlipNorHex[3](1, 0) = 0.0;
    SlipNorHex[3](2, 0) = 0.0;
    SlipNorHex[3](3, 0) = 1.0;

    SlipDirHex[3](0, 0) = 1.0 / 3.0;
    SlipDirHex[3](1, 0) = 1.0 / 3.0;
    SlipDirHex[3](2, 0) = -2.0 / 3.0;
    SlipDirHex[3](3, 0) = 0.0;

    SlipSysID_[4] = "basal_2";

    SlipNorHex[4](0, 0) = 0.0;
    SlipNorHex[4](1, 0) = 0.0;
    SlipNorHex[4](2, 0) = 0.0;
    SlipNorHex[4](3, 0) = 1.0;

    SlipDirHex[4](0, 0) = -2.0 / 3.0;
    SlipDirHex[4](1, 0) = 1.0 / 3.0;
    SlipDirHex[4](2, 0) = 1.0 / 3.0;
    SlipDirHex[4](3, 0) = 0.0;

    SlipSysID_[5] = "basal_3";

    SlipNorHex[5](0, 0) = 0.0;
    SlipNorHex[5](1, 0) = 0.0;
    SlipNorHex[5](2, 0) = 0.0;
    SlipNorHex[5](3, 0) = 1.0;

    SlipDirHex[5](0, 0) = 1.0 / 3.0;
    SlipDirHex[5](1, 0) = -2.0 / 3.0;
    SlipDirHex[5](2, 0) = 1.0 / 3.0;
    SlipDirHex[5](3, 0) = 0.0;

    // pyramidal slip systems

    SlipSysID_[6] = "pyramidal_1";

    SlipNorHex[6](0, 0) = 1.0;
    SlipNorHex[6](1, 0) = 1.0;
    SlipNorHex[6](2, 0) = -2.0;
    SlipNorHex[6](3, 0) = 1.0;

    SlipDirHex[6](0, 0) = -1.0 / 3.0;
    SlipDirHex[6](1, 0) = -1.0 / 3.0;
    SlipDirHex[6](2, 0) = 2.0 / 3.0;
    SlipDirHex[6](3, 0) = 6.0 / 3.0;

    SlipSysID_[7] = "pyramidal_2";

    SlipNorHex[7](0, 0) = 1.0;
    SlipNorHex[7](1, 0) = -2.0;
    SlipNorHex[7](2, 0) = 1.0;
    SlipNorHex[7](3, 0) = 1.0;

    SlipDirHex[7](0, 0) = -1.0 / 3.0;
    SlipDirHex[7](1, 0) = 2.0 / 3.0;
    SlipDirHex[7](2, 0) = -1.0 / 3.0;
    SlipDirHex[7](3, 0) = 6.0 / 3.0;

    SlipSysID_[8] = "pyramidal_3";

    SlipNorHex[8](0, 0) = -2.0;
    SlipNorHex[8](1, 0) = 1.0;
    SlipNorHex[8](2, 0) = 1.0;
    SlipNorHex[8](3, 0) = 1.0;

    SlipDirHex[8](0, 0) = 2.0 / 3.0;
    SlipDirHex[8](1, 0) = -1.0 / 3.0;
    SlipDirHex[8](2, 0) = -1.0 / 3.0;
    SlipDirHex[8](3, 0) = 6.0 / 3.0;

    SlipSysID_[9] = "pyramidal_4";

    SlipNorHex[9](0, 0) = -1.0;
    SlipNorHex[9](1, 0) = -1.0;
    SlipNorHex[9](2, 0) = 2.0;
    SlipNorHex[9](3, 0) = 1.0;

    SlipDirHex[9](0, 0) = 1.0 / 3.0;
    SlipDirHex[9](1, 0) = 1.0 / 3.0;
    SlipDirHex[9](2, 0) = -2.0 / 3.0;
    SlipDirHex[9](3, 0) = 6.0 / 3.0;

    SlipSysID_[10] = "pyramidal_5";

    SlipNorHex[10](0, 0) = -1.0;
    SlipNorHex[10](1, 0) = 2.0;
    SlipNorHex[10](2, 0) = -1.0;
    SlipNorHex[10](3, 0) = 1.0;

    SlipDirHex[10](0, 0) = 1.0 / 3.0;
    SlipDirHex[10](1, 0) = -2.0 / 3.0;
    SlipDirHex[10](2, 0) = 1.0 / 3.0;
    SlipDirHex[10](3, 0) = 6.0 / 3.0;

    SlipSysID_[11] = "pyramidal_6";

    SlipNorHex[11](0, 0) = 2.0;
    SlipNorHex[11](1, 0) = -1.0;
    SlipNorHex[11](2, 0) = -1.0;
    SlipNorHex[11](3, 0) = 1.0;

    SlipDirHex[11](0, 0) = -2.0 / 3.0;
    SlipDirHex[11](1, 0) = 1.0 / 3.0;
    SlipDirHex[11](2, 0) = 1.0 / 3.0;
    SlipDirHex[11](3, 0) = 6.0 / 3.0;

    // transform Miller-Bravais index notation of hexagonal lattices to Miller index notation of
    // cubic lattices
    MAT::CrystalPlasticity::MillerBravaisToMiller(SlipNorHex, SlipDirHex);
  }
  else if (LatticeType_ == "HCP" or LatticeType_ == "BCC" or LatticeType_ == "FCC")
  {
    dserror("Sorry, %s lattices are not yet supported", LatticeType_.c_str());
  }
  // academic test lattice containing only one slip system
  else if (LatticeType_ == "TEST")
  {
    SlipSysID_[0] = "test_1";

    SlipNor_[0](0, 0) = 1.0;
    SlipNor_[0](1, 0) = 0.0;
    SlipNor_[0](2, 0) = 1.0;

    SlipDir_[0](0, 0) = 1.0;
    SlipDir_[0](1, 0) = 0.0;
    SlipDir_[0](2, 0) = -1.0;
  }
  else
  {
    dserror(
        "Lattice type not known. Please check the LAT parameter in the input file. Currently it "
        "has to be FCC, BCC, HCP, D019 or L10");
  }


  // test whether directions and normals are perpendicular
  std::vector<LINALG::Matrix<1, 1>> PerpTest;
  PerpTest.resize(slip_sys_count_);

  for (int i = 0; i < slip_sys_count_; i++)
  {
    PerpTest[i].MultiplyTN(SlipDir_[i], SlipNor_[i]);
    if (PerpTest[i](0, 0) != 0.0)
    {
      dserror(
          "Warning, slip direction and slip plane normal of slip system "
          "%d are not perpendicular!",
          i + 1);
    }
  }

  // determine magnitude of Burgers vector
  Burgers_Mag_.resize(slip_sys_count_);
  LINALG::Matrix<3, 1> tmp_scale;

  // import lattice constant from user input
  LatticeConstant_ = params_->lattice_constant_;

  for (int i = 0; i < slip_sys_count_; i++)
  {
    tmp_scale = SlipDir_[i];
    tmp_scale.Scale(LatticeConstant_);
    Burgers_Mag_[i] = tmp_scale.Norm2();
  }

  // normalize slip plane normals
  for (int i = 0; i < slip_sys_count_; i++)
  {
    SlipNor_[i].Scale(1.0 / SlipNor_[i].Norm2());
  }

  // normalize slip directions
  for (int i = 0; i < slip_sys_count_; i++)
  {
    SlipDir_[i].Scale(1.0 / SlipDir_[i].Norm2());
  }
  return;
}


/*---------------------------------------------------------------------------------*
 | Read Lattice orientation matrix from .dat file                                    |
 *---------------------------------------------------------------------------------*/

void MAT::CrystalPlasticity::SetupLatticeOrientation(DRT::INPUT::LineDefinition* linedef)
{
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
      LatticeOrientation_(i, 0) = fiber1.at(i);
      LatticeOrientation_(i, 1) = fiber2.at(i);
      LatticeOrientation_(i, 2) = fiber3.at(i);
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
    SlipNor_[i](0, 0) = SlipNorHex.at(i)(0, 0);
    SlipNor_[i](1, 0) = (SlipNorHex.at(i)(0, 0) + 2.0 * SlipNorHex.at(i)(1, 0)) / sqrt(3.0);
    SlipNor_[i](2, 0) = SlipNorHex.at(i)(3, 0) / CtoA_;

    SlipDir_[i](0, 0) = 1.5 * SlipDirHex.at(i)(0, 0);
    SlipDir_[i](1, 0) = sqrt(3.0) * (0.5 * SlipDirHex.at(i)(0, 0) + SlipDirHex.at(i)(1, 0));
    SlipDir_[i](2, 0) = SlipDirHex.at(i)(3, 0) * CtoA_;
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
  std::vector<double> gamma_trial = gamma_last_->at(gp);

  // trial values of internal variables
  std::vector<double> def_dens_trial(slip_sys_count_);
  LINALG::Matrix<3, 3> PK2_trial;
  LINALG::Matrix<3, 3> FP_trial;

  // residuals
  std::vector<double> residuals_trial(slip_sys_count_);
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
    for (int i = 0; i < slip_sys_count_; i++)
    {
      totalResidual += pow(residuals_trial.at(i), 2.0);
    }
    totalResidual = sqrt(totalResidual);

    // convergence check
    if (totalResidual < Tol_)
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
      for (int i = 0; i < slip_sys_count_; i++)
      {
        // perturbed vector of plastic shears
        std::vector<double> gamma_pert(slip_sys_count_);
        // resultant vector of residuals
        std::vector<double> residuals_pert(slip_sys_count_);
        // resultant vector of defect densities
        std::vector<double> def_dens_pert(slip_sys_count_);
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
        if (LatticeType_ == "L10" || LatticeType_ == "D019")
        {
          for (int j = 0; j < slip_sys_count_; j++)
          {
            ResidualStiffness_12x12(j, i) =
                (residuals_pert.at(j) - residuals_trial.at(j)) / gammaEps;

            residuals_trial_LIN_12(i) = -residuals_trial.at(i);
          }
        }
        else if (LatticeType_ == "TEST")
        {
          for (int j = 0; j < slip_sys_count_; j++)
          {
            ResidualStiffness_1x1(j, i) = (residuals_pert.at(j) - residuals_trial.at(j)) / gammaEps;

            residuals_trial_LIN_1(i) = -residuals_trial.at(i);
          }
        }
      }

      // solve resultant system of equations
      if (LatticeType_ == "L10" || LatticeType_ == "D019")
      {
        NRSolver_12x12.SetMatrix(ResidualStiffness_12x12);
        NRSolver_12x12.SetVectors(d_gamma_trial_12, residuals_trial_LIN_12);
        NRSolver_12x12.Solve();

        // update vector of plastic shears to new trial value
        for (int i = 0; i < slip_sys_count_; i++)
        {
          gamma_trial.at(i) += d_gamma_trial_12(i);
        }
      }
      else if (LatticeType_ == "TEST")
      {
        NRSolver_1x1.SetMatrix(ResidualStiffness_1x1);
        NRSolver_1x1.SetVectors(d_gamma_trial_1, residuals_trial_LIN_1);
        NRSolver_1x1.Solve();

        // update vector of plastic shears to new trial value
        for (int i = 0; i < slip_sys_count_; i++)
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
  std::vector<double> delta_gamma_trial(slip_sys_count_);

  for (int i = 0; i < slip_sys_count_; i++)
  {
    delta_gamma_trial[i] = gamma_trial[i] - (*gamma_last_)[gp][i];
  }

  // determine trial defect densities that would result from delta_gamma_trial
  //--------------------------------------------------------------------------

  double total_dislocation_density = 0.0;

  for (int i = 0; i < slip_sys_count_; i++)
  {
    // dislocation generation
    //  TODO EQUATION
    def_dens_trial[i] = (*def_dens_last_)[gp][i] +
                        (1.0 / (Burgers_Mag_[i] * Dislocation_Gen_Coeff_[SubSetIndex_[i] - 1])) *
                            sqrt((*def_dens_last_)[gp][i]) * abs(delta_gamma_trial[i]);

    // dynamic recovery
    def_dens_trial[i] -= Dislocation_Dyn_Rec_Coeff_[SubSetIndex_[i] - 1] *
                         (*def_dens_last_)[gp][i] * abs(delta_gamma_trial[i]);

    // determine updated total defect densities
    total_dislocation_density += def_dens_trial[i];
  }

  // kinematics and plastic velocity gradient
  //--------------------------------------------------------------------------

  // determine trial plastic velocity gradient L_p
  // set up L_p_trial
  LINALG::Matrix<3, 3> LP_trial(true);
  LINALG::Matrix<3, 3> TempMat;

  for (int i = 0; i < slip_sys_count_; i++)
  {
    TempMat.MultiplyNT(delta_gamma_trial.at(i), SlipDir_[i], SlipNor_[i]);
    LP_trial.Update(LP_trial, TempMat);
  }

  // take unimodular part of I + L_p to ensure plastic incompressibility
  LINALG::Matrix<3, 3> Uni_IplusLP_trial(true);
  Uni_IplusLP_trial.Update(Identity3(), LP_trial);
  Uni_IplusLP_trial.Scale(pow(Uni_IplusLP_trial.Determinant(), -1.0 / 3.0));

  // determine trial plastic deformation gradient
  FP_trial.MultiplyNN(Uni_IplusLP_trial, (*FP_last_)[gp]);

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

  PK2_trial.Update(Lambda_ * ln_Jacobi_trial, inv_CE, 1.0);
  PK2_trial.Update(Mu_, Identity3(), 1.0);
  PK2_trial.Update(-Mu_, inv_CE, 1.0);

  // Mandel stress
  LINALG::Matrix<3, 3> M;
  M.MultiplyNN(CE, PK2_trial);

  // setting up Residua with flow/creep rule
  //--------------------------------------------------------------------------

  for (int i = 0; i < slip_sys_count_; i++)
  {
    // resolved shear stress/Schmid stress
    LINALG::Matrix<1, 1> RSS(true);

    LINALG::Matrix<3, 1> TempVec(true);

    TempVec.MultiplyNN(M, SlipNor_[i]);

    RSS.MultiplyTN(SlipDir_[i], TempVec);

    // work hardening

    // lattice resistance
    double tau_Y0 = TauY_0_[SubSetIndex_[i] - 1];

    // Hall-Petch strengthening term
    double HP_strengthening = Hall_Petch_Coeffs_[SubSetIndex_[i] - 1] *
                              (1.0 / sqrt(Micro_Boundary_Dist_[SubSetIndex_[i] - 1]));

    // work hardening hardening increment due to accumulated dislocation density
    double delta_tau_Y = 0.5 * Mu_ * Burgers_Mag_[i] * sqrt(total_dislocation_density);

    // slip system strength, i.e. tau_Y = tau_Y0 + HP_strengthening + delta_tau_Y
    double tau_Y = tau_Y0 + HP_strengthening + delta_tau_Y;
    // check for consitency
    if (tau_Y < 0.0) dserror("Negative slip systems strength! Please check your input");

    // shear rate as determined from the flow rule
    double gamma_dot = 0.0;

    gamma_dot = abs(RSS(0, 0)) / tau_Y;
    gamma_dot = pow(gamma_dot, 50.0);  // TODO EXPONENT AUS DAT FILE HOLEN
    gamma_dot = gamma_dot * Gamma0_[SubSetIndex_[i] - 1];
    gamma_dot = std::copysign(gamma_dot, RSS(0, 0));

    // set up corresponding residual
    residuals_trial[i] = delta_gamma_trial[i] / dt - gamma_dot;
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
  for (int sys = 0; sys < slip_sys_count_; sys++)
  {
    ID = "plastic_shear_";
    ID += SlipSysID_[sys];
    names[ID] = 1;  // scalar
  }

  // accumulated plastic shears
  ID = "accumulated_plastic_shear";
  names[ID] = 1;  // scalar

  // dislocation densities on single systems
  for (int sys = 0; sys < slip_sys_count_; sys++)
  {
    ID = "dislocation_density_";
    ID += SlipSysID_[sys];
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
  for (int sys = 0; sys < slip_sys_count_; sys++)
  {
    if (name == "plastic_shear_" + SlipSysID_[sys])
    {
      if ((int)data.size() != 1) dserror("size mismatch");
      for (int gp_index = 0; gp_index < numgp; gp_index++)
      {
        data[0] = (*gamma_last_)[gp_index][sys];
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
      for (int sys = 0; sys < slip_sys_count_; sys++)
      {
        gamma_acc += abs((*gamma_last_)[gp_index][sys]);
      }
      data[0] = gamma_acc;
    }
  }

  // dislocation densities
  for (int sys = 0; sys < slip_sys_count_; sys++)
  {
    if (name == "dislocation_density_" + SlipSysID_[sys])
    {
      if ((int)data.size() != 1) dserror("size mismatch");
      for (int gp_index = 0; gp_index < numgp; gp_index++)
      {
        data[0] = (*def_dens_last_)[gp_index][sys];
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
      for (int sys = 0; sys < slip_sys_count_; sys++)
      {
        dis_dens_acc += (*def_dens_last_)[gp_index][sys];
      }
      data[0] = dis_dens_acc;
    }
  }

  return true;
}  // VisData()
