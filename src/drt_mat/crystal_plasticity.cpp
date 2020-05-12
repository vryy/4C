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
      poisson_(matdata->GetDouble("NUE")),
      density_(matdata->GetDouble("DENS")),
      lattice_(*matdata->Get<std::string>("LAT")),
      ctoa_(matdata->GetDouble("CTOA")),
      lattice_const_(matdata->GetDouble("ABASE")),
      num_slip_sys_(matdata->GetInt("NUMSLIPSYS")),
      num_sub_sets_(matdata->GetInt("NUMSUBSETS")),
      sub_set_mem_(*matdata->Get<std::vector<int>>("SUBSETMEMBERS")),
      rate_exp_(*matdata->Get<std::vector<int>>("RATEEXP")),
      ref_shear_rate_(*matdata->Get<std::vector<double>>("GAMMADOTREF")),
      dis_gen_coeff_(*matdata->Get<std::vector<double>>("DISGENCOEFF")),
      dis_dyn_rec_coeff_(*matdata->Get<std::vector<double>>("DISDYNRECCOEFF")),
      lat_resist_(*matdata->Get<std::vector<double>>("TAUY0")),
      micro_bound_(*matdata->Get<std::vector<double>>("GBS")),
      hp_coeff_(*matdata->Get<std::vector<double>>("HPCOEFF"))
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
  if (!(isinit_ and (deform_grad_current_ != Teuchos::null)))
  {
    histsize = 0;
  }
  else
  {
    // if material is initialized already (restart): size equates number of Gauss points
    histsize = deform_grad_last_->size();
  }
  AddtoPack(data, histsize);  // length of history vector(s)

  for (int var = 0; var < histsize; ++var)
  {
    // insert history vectors to AddtoPack
    AddtoPack(data, (*deform_grad_last_)[var]);
    AddtoPack(data, (*plastic_deform_grad_last_)[var]);
    AddtoPack(data, (*gamma_last_)[var]);
    AddtoPack(data, (*defect_densities_last_)[var]);
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
    deform_grad_last_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
    plastic_deform_grad_last_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
    gamma_last_ = Teuchos::rcp(new std::vector<std::vector<double>>);
    defect_densities_last_ = Teuchos::rcp(new std::vector<std::vector<double>>);

    for (int var = 0; var < histsize; ++var)
    {
      LINALG::Matrix<3, 3> tmp_matrix(true);
      std::vector<double> tmp_vect(slip_system_count_);

      ExtractfromPack(position, data, tmp_matrix);
      deform_grad_last_->push_back(tmp_matrix);

      ExtractfromPack(position, data, tmp_matrix);
      plastic_deform_grad_last_->push_back(tmp_matrix);

      ExtractfromPack(position, data, tmp_vect);
      gamma_last_->push_back(tmp_vect);

      ExtractfromPack(position, data, tmp_vect);
      defect_densities_last_->push_back(tmp_vect);
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
  newton_tolerance_ = params_->tol_;

  // Elastic Properties:
  //-------------------------
  // Young's Modulus
  youngs_mod_ = params_->youngs_;
  // Poisson's ratio
  poisson_ratio_ = params_->poisson_;
  // 1st lamé constant
  lambda_ =
      (poisson_ratio_ * youngs_mod_) / ((1.0 + poisson_ratio_) * (1.0 - 2.0 * poisson_ratio_));
  // Shear modulus / 2nd Lame constant
  mue_ = youngs_mod_ / (2.0 * (1.0 + poisson_ratio_));
  // Bulk modulus
  bulk_mod_ = youngs_mod_ / (3.0 - (6.0 * poisson_ratio_));

  // Crystal Properties:
  //-------------------
  // setup lattice vectors
  // set up slip plane normals and directions and initialize crystal lattice related model
  // parameters according to chosen lattice type
  this->SetupLatticeVectors();

  //  read Lattice orientation matrix from .dat file
  this->SetupLatticeOrientation(linedef);

  // rotate lattice vectors according to lattice orientation
  for (int i = 0; i < slip_system_count_; i++)
  {
    LINALG::Matrix<3, 1> unrotated_normal = slip_plane_normal_[i];
    LINALG::Matrix<3, 1> unrotated_direction = slip_direction_[i];

    slip_plane_normal_[i].MultiplyNN(lattice_orientation_, unrotated_normal);
    slip_direction_[i].MultiplyNN(lattice_orientation_, unrotated_direction);
  }
  // TODO ADD TEST WHETHER INPUT ROTATION MATRIX IS ORTHOGONAL Q*Q^T=I TO AVOID UNWANTED SCALING OF
  // LATTICE VECTORS!

  // Index to which subset a slip system belongs
  subset_index_ = params_->sub_set_mem_;

  // Viscoplastic Properties:
  //------------------------
  // TODO RATE EXPONENT INPUT
  // reference shear rates
  gamma_dot_0_ = params_->ref_shear_rate_;

  // Dislocation Generation/Recovery:
  //--------------------------------
  // dislocation generation coefficients
  dislocation_generation_coeff_ = params_->dis_gen_coeff_;
  // dynamic dislocation removal coefficients
  dislocation_dyn_recovery_coeff_ = params_->dis_dyn_rec_coeff_;

  // Initial Slip System Strengths:
  //------------------------------
  // lattice resistances to slip
  tau_y_0_ = params_->lat_resist_;
  // microstructural parameters which are relevant for Hall-Petch strengthening, e.g., grain size
  micro_boundary_distance_ = params_->micro_bound_;
  // Hall-Petch coefficients corresponding to above microstructural boundaries
  hall_petch_coeffs_ = params_->hp_coeff_;

  // set up 3x3 identity matrix
  LINALG::IdentityMatrix(identity3);

  // initialize history variables

  deform_grad_last_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  deform_grad_current_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  plastic_deform_grad_last_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  plastic_deform_grad_current_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);

  gamma_last_ = Teuchos::rcp(new std::vector<std::vector<double>>);
  gamma_current_ = Teuchos::rcp(new std::vector<std::vector<double>>);

  defect_densities_last_ = Teuchos::rcp(new std::vector<std::vector<double>>);
  defect_density_current_ = Teuchos::rcp(new std::vector<std::vector<double>>);

  // set size to number of  Gauss points so each Gauss Point has its own set of history data
  deform_grad_last_->resize(numgp);
  deform_grad_current_->resize(numgp);

  plastic_deform_grad_last_->resize(numgp);
  plastic_deform_grad_current_->resize(numgp);

  gamma_last_->resize(numgp);
  gamma_current_->resize(numgp);

  defect_densities_last_->resize(numgp);
  defect_density_current_->resize(numgp);

  // set initial values
  std::vector<double> emptyvect(slip_system_count_);
  for (int i = 0; i < slip_system_count_; i++) emptyvect.at(i) = 0.0;

  for (int i = 0; i < numgp; i++)
  {
    (*deform_grad_last_)[i] = identity3;
    (*deform_grad_current_)[i] = identity3;

    (*plastic_deform_grad_last_)[i] = identity3;
    (*plastic_deform_grad_current_)[i] = identity3;

    (*gamma_last_)[i].resize(slip_system_count_);
    (*gamma_current_)[i].resize(slip_system_count_);

    (*gamma_last_)[i] = emptyvect;
    (*gamma_current_)[i] = emptyvect;

    (*defect_densities_last_)[i].resize(slip_system_count_);
    (*defect_density_current_)[i].resize(slip_system_count_);

    for (int j = 0; j < slip_system_count_; j++)
      (*defect_densities_last_)[i][j] = 1e5;  // TODO GET INITIAL DIS DENS ROM INPUT LINE

    (*defect_density_current_)[i] = emptyvect;
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
  deform_grad_last_ = deform_grad_current_;
  plastic_deform_grad_last_ = plastic_deform_grad_current_;
  gamma_last_ = gamma_current_;
  defect_densities_last_ = defect_density_current_;

  // empty vectors of current data
  deform_grad_current_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  plastic_deform_grad_current_ = Teuchos::rcp(new std::vector<LINALG::Matrix<3, 3>>);
  gamma_current_ = Teuchos::rcp(new std::vector<std::vector<double>>);
  defect_density_current_ = Teuchos::rcp(new std::vector<std::vector<double>>);

  // get the size of the vector (use the last vector, because it includes latest results, current is
  // already empty)
  const int histsize = deform_grad_last_->size();
  deform_grad_current_->resize(histsize);
  plastic_deform_grad_current_->resize(histsize);
  gamma_current_->resize(histsize);
  defect_density_current_->resize(histsize);

  LINALG::Matrix<3, 3> emptymat(true);
  std::vector<double> emptyvect(slip_system_count_);
  for (int i = 0; i < slip_system_count_; i++)
  {
    emptyvect.at(i) = 0.0;
  }

  for (int i = 0; i < histsize; i++)
  {
    (*deform_grad_current_)[i] = emptymat;
    (*plastic_deform_grad_current_)[i] = emptymat;
    (*gamma_current_)[i] = emptyvect;
    (*defect_density_current_)[i] = emptyvect;
  }


  return;
}  // Update()

/*-------------------------------------------------------------------------------*
 | calculate stress, material stiffness matrix and evolution internal variables  |
 *-------------------------------------------------------------------------------*/
void MAT::CrystalPlasticity::Evaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, const int gp, const int eleGID)
{
  // simulation parameters
  //---------------------------------------------------------------
  if (eleGID == -1) dserror("no element provided in material");

  // extract time increment
  dt = params.get<double>("delta time");

  // set current deformation gradient
  (*deform_grad_current_)[gp] = *defgrd;


  // local Newton-Raphson to determine plastic shears
  //--------------------------------------------------------------------------

  // define result output from NR
  LINALG::Matrix<3, 3> second_pk_stress_result;

  // call Newton-Raphson method with current deformation gradient
  this->NewtonRaphson((*deform_grad_current_)[gp], (*gamma_current_)[gp],
      (*defect_density_current_)[gp], second_pk_stress_result, (*plastic_deform_grad_current_)[gp]);

  // update stress results
  (*stress)(0) = second_pk_stress_result(0, 0);
  (*stress)(1) = second_pk_stress_result(1, 1);
  (*stress)(2) = second_pk_stress_result(2, 2);
  (*stress)(3) = second_pk_stress_result(0, 1);
  (*stress)(4) = second_pk_stress_result(1, 2);
  (*stress)(5) = second_pk_stress_result(0, 2);

  return;
}  // Evaluate()

/*---------------------------------------------------------------------------------*
 | Setup slip directions and slip plane normals                                    |
 *---------------------------------------------------------------------------------*/
void MAT::CrystalPlasticity::SetupLatticeVectors()
{
  // extract lattice type from user input
  lattice_type_ = params_->lattice_;

  // assign number of slip systems that corresponds to the given lattice
  if (lattice_type_ == "L10" || lattice_type_ == "D019")
  {
    slip_system_count_ = 12;
  }
  else if (lattice_type_ == "TEST")
  {
    slip_system_count_ = 1;
  }

  // check whether the number of slip systems given by the user coincides with the lattice type
  if (params_->num_slip_sys_ != slip_system_count_)
  {
    dserror(
        "Given number of slip systems NUMSLIPSYS = %d does not match the expected number for "
        "%s lattices!",
        params_->num_slip_sys_, lattice_type_.c_str());
  }

  // resize to number of slip systems
  slip_plane_normal_.resize(slip_system_count_);
  slip_direction_.resize(slip_system_count_);
  slip_system_id_.resize(slip_system_count_);

  // import c to a ratio of unit cell as provided by the user
  c_to_a_ratio_ = params_->ctoa_;

  // set up corresponding set of slip planes depending on the chosen crystallographic lattice
  // L10 super lattice
  if (lattice_type_ == "L10")
  {
    // ordinary slip systems
    slip_system_id_[0] = "ordinary_1";

    slip_plane_normal_[0](0, 0) = 1.0;
    slip_plane_normal_[0](1, 0) = 1.0;
    slip_plane_normal_[0](2, 0) = 1.0 / c_to_a_ratio_;

    slip_direction_[0](0, 0) = 1.0 / 2.0;
    slip_direction_[0](1, 0) = -1.0 / 2.0;
    slip_direction_[0](2, 0) = 0.0 * c_to_a_ratio_;

    slip_system_id_[1] = "ordinary_2";

    slip_plane_normal_[1](0, 0) = -1.0;
    slip_plane_normal_[1](1, 0) = -1.0;
    slip_plane_normal_[1](2, 0) = 1.0 / c_to_a_ratio_;

    slip_direction_[1](0, 0) = 1.0 / 2.0;
    slip_direction_[1](1, 0) = -1.0 / 2.0;
    slip_direction_[1](2, 0) = 0.0 * c_to_a_ratio_;

    slip_system_id_[2] = "ordinary_3";

    slip_plane_normal_[2](0, 0) = 1.0;
    slip_plane_normal_[2](1, 0) = -1.0;
    slip_plane_normal_[2](2, 0) = 1.0 / c_to_a_ratio_;

    slip_direction_[2](0, 0) = 1.0 / 2.0;
    slip_direction_[2](1, 0) = 1.0 / 2.0;
    slip_direction_[2](2, 0) = 0.0 * c_to_a_ratio_;

    slip_system_id_[3] = "ordinary_4";

    slip_plane_normal_[3](0, 0) = -1.0;
    slip_plane_normal_[3](1, 0) = 1.0;
    slip_plane_normal_[3](2, 0) = 1.0 / c_to_a_ratio_;

    slip_direction_[3](0, 0) = 1.0 / 2.0;
    slip_direction_[3](1, 0) = 1.0 / 2.0;
    slip_direction_[3](2, 0) = 0.0 * c_to_a_ratio_;

    // super slip systems
    slip_system_id_[4] = "super_1";

    slip_plane_normal_[4](0, 0) = 1.0;
    slip_plane_normal_[4](1, 0) = 1.0;
    slip_plane_normal_[4](2, 0) = 1.0 / c_to_a_ratio_;

    slip_direction_[4](0, 0) = 0.0;
    slip_direction_[4](1, 0) = 1.0;
    slip_direction_[4](2, 0) = -1.0 * c_to_a_ratio_;

    slip_system_id_[5] = "super_2";

    slip_plane_normal_[5](0, 0) = 1.0;
    slip_plane_normal_[5](1, 0) = 1.0;
    slip_plane_normal_[5](2, 0) = 1.0 / c_to_a_ratio_;

    slip_direction_[5](0, 0) = 1.0;
    slip_direction_[5](1, 0) = 0.0;
    slip_direction_[5](2, 0) = -1.0 * c_to_a_ratio_;

    slip_system_id_[6] = "super_3";

    slip_plane_normal_[6](0, 0) = 1.0;
    slip_plane_normal_[6](1, 0) = -1.0;
    slip_plane_normal_[6](2, 0) = -1.0 / c_to_a_ratio_;

    slip_direction_[6](0, 0) = 0.0;
    slip_direction_[6](1, 0) = 1.0;
    slip_direction_[6](2, 0) = -1.0 * c_to_a_ratio_;

    slip_system_id_[7] = "super_4";

    slip_plane_normal_[7](0, 0) = 1.0;
    slip_plane_normal_[7](1, 0) = -1.0;
    slip_plane_normal_[7](2, 0) = 1.0 / c_to_a_ratio_;

    slip_direction_[7](0, 0) = 1.0;
    slip_direction_[7](1, 0) = 0.0;
    slip_direction_[7](2, 0) = -1.0 * c_to_a_ratio_;

    slip_system_id_[8] = "super_5";

    slip_plane_normal_[8](0, 0) = 1.0;
    slip_plane_normal_[8](1, 0) = 1.0;
    slip_plane_normal_[8](2, 0) = -1.0 / c_to_a_ratio_;

    slip_direction_[8](0, 0) = 0.0;
    slip_direction_[8](1, 0) = -1.0;
    slip_direction_[8](2, 0) = -1.0 * c_to_a_ratio_;

    slip_system_id_[9] = "super_6";

    slip_plane_normal_[9](0, 0) = 1.0;
    slip_plane_normal_[9](1, 0) = 1.0;
    slip_plane_normal_[9](2, 0) = -1.0 / c_to_a_ratio_;

    slip_direction_[9](0, 0) = -1.0;
    slip_direction_[9](1, 0) = 0.0;
    slip_direction_[9](2, 0) = -1.0 * c_to_a_ratio_;

    slip_system_id_[10] = "super_7";

    slip_plane_normal_[10](0, 0) = 1.0;
    slip_plane_normal_[10](1, 0) = -1.0;
    slip_plane_normal_[10](2, 0) = 1.0 / c_to_a_ratio_;

    slip_direction_[10](0, 0) = 0.0;
    slip_direction_[10](1, 0) = -1.0;
    slip_direction_[10](2, 0) = -1.0 * c_to_a_ratio_;

    slip_system_id_[11] = "super_8";

    slip_plane_normal_[11](0, 0) = 1.0;
    slip_plane_normal_[11](1, 0) = -1.0;
    slip_plane_normal_[11](2, 0) = -1.0 / c_to_a_ratio_;

    slip_direction_[11](0, 0) = -1.0;
    slip_direction_[11](1, 0) = 0.0;
    slip_direction_[11](2, 0) = -1.0 * c_to_a_ratio_;
  }
  // D019 super lattice
  else if (lattice_type_ == "D019")
  {
    // initialize slip plane normals and directions for Miller-Bravais indices
    std::vector<LINALG::Matrix<4, 1>> slip_plane_normal_hex(slip_system_count_);
    std::vector<LINALG::Matrix<4, 1>> slip_direction_hex(slip_system_count_);

    // NOTE: the c to a ratio is incorporated during transformation from Miller-Bravais to Miller

    // prismatic slip systems
    slip_system_id_[0] = "prismatic_1";

    slip_plane_normal_hex[0](0, 0) = 1.0;
    slip_plane_normal_hex[0](1, 0) = -1.0;
    slip_plane_normal_hex[0](2, 0) = 0.0;
    slip_plane_normal_hex[0](3, 0) = 0.0;

    slip_direction_hex[0](0, 0) = 1.0 / 3.0;
    slip_direction_hex[0](1, 0) = 1.0 / 3.0;
    slip_direction_hex[0](2, 0) = -2.0 / 3.0;
    slip_direction_hex[0](3, 0) = 0.0;

    slip_system_id_[1] = "prismatic_2";

    slip_plane_normal_hex[1](0, 0) = 0.0;
    slip_plane_normal_hex[1](1, 0) = 1.0;
    slip_plane_normal_hex[1](2, 0) = -1.0;
    slip_plane_normal_hex[1](3, 0) = 0.0;

    slip_direction_hex[1](0, 0) = -2.0 / 3.0;
    slip_direction_hex[1](1, 0) = 1.0 / 3.0;
    slip_direction_hex[1](2, 0) = 1.0 / 3.0;
    slip_direction_hex[1](3, 0) = 0.0;

    slip_system_id_[2] = "prismatic_3";

    slip_plane_normal_hex[2](0, 0) = 1.0;
    slip_plane_normal_hex[2](1, 0) = 0.0;
    slip_plane_normal_hex[2](2, 0) = -1.0;
    slip_plane_normal_hex[2](3, 0) = 0.0;

    slip_direction_hex[2](0, 0) = 1.0 / 3.0;
    slip_direction_hex[2](1, 0) = -2.0 / 3.0;
    slip_direction_hex[2](2, 0) = 1.0 / 3.0;
    slip_direction_hex[2](3, 0) = 0.0;

    // basal slip systems

    slip_system_id_[3] = "basal_1";

    slip_plane_normal_hex[3](0, 0) = 0.0;
    slip_plane_normal_hex[3](1, 0) = 0.0;
    slip_plane_normal_hex[3](2, 0) = 0.0;
    slip_plane_normal_hex[3](3, 0) = 1.0;

    slip_direction_hex[3](0, 0) = 1.0 / 3.0;
    slip_direction_hex[3](1, 0) = 1.0 / 3.0;
    slip_direction_hex[3](2, 0) = -2.0 / 3.0;
    slip_direction_hex[3](3, 0) = 0.0;

    slip_system_id_[4] = "basal_2";

    slip_plane_normal_hex[4](0, 0) = 0.0;
    slip_plane_normal_hex[4](1, 0) = 0.0;
    slip_plane_normal_hex[4](2, 0) = 0.0;
    slip_plane_normal_hex[4](3, 0) = 1.0;

    slip_direction_hex[4](0, 0) = -2.0 / 3.0;
    slip_direction_hex[4](1, 0) = 1.0 / 3.0;
    slip_direction_hex[4](2, 0) = 1.0 / 3.0;
    slip_direction_hex[4](3, 0) = 0.0;

    slip_system_id_[5] = "basal_3";

    slip_plane_normal_hex[5](0, 0) = 0.0;
    slip_plane_normal_hex[5](1, 0) = 0.0;
    slip_plane_normal_hex[5](2, 0) = 0.0;
    slip_plane_normal_hex[5](3, 0) = 1.0;

    slip_direction_hex[5](0, 0) = 1.0 / 3.0;
    slip_direction_hex[5](1, 0) = -2.0 / 3.0;
    slip_direction_hex[5](2, 0) = 1.0 / 3.0;
    slip_direction_hex[5](3, 0) = 0.0;

    // pyramidal slip systems

    slip_system_id_[6] = "pyramidal_1";

    slip_plane_normal_hex[6](0, 0) = 1.0;
    slip_plane_normal_hex[6](1, 0) = 1.0;
    slip_plane_normal_hex[6](2, 0) = -2.0;
    slip_plane_normal_hex[6](3, 0) = 1.0;

    slip_direction_hex[6](0, 0) = -1.0 / 3.0;
    slip_direction_hex[6](1, 0) = -1.0 / 3.0;
    slip_direction_hex[6](2, 0) = 2.0 / 3.0;
    slip_direction_hex[6](3, 0) = 6.0 / 3.0;

    slip_system_id_[7] = "pyramidal_2";

    slip_plane_normal_hex[7](0, 0) = 1.0;
    slip_plane_normal_hex[7](1, 0) = -2.0;
    slip_plane_normal_hex[7](2, 0) = 1.0;
    slip_plane_normal_hex[7](3, 0) = 1.0;

    slip_direction_hex[7](0, 0) = -1.0 / 3.0;
    slip_direction_hex[7](1, 0) = 2.0 / 3.0;
    slip_direction_hex[7](2, 0) = -1.0 / 3.0;
    slip_direction_hex[7](3, 0) = 6.0 / 3.0;

    slip_system_id_[8] = "pyramidal_3";

    slip_plane_normal_hex[8](0, 0) = -2.0;
    slip_plane_normal_hex[8](1, 0) = 1.0;
    slip_plane_normal_hex[8](2, 0) = 1.0;
    slip_plane_normal_hex[8](3, 0) = 1.0;

    slip_direction_hex[8](0, 0) = 2.0 / 3.0;
    slip_direction_hex[8](1, 0) = -1.0 / 3.0;
    slip_direction_hex[8](2, 0) = -1.0 / 3.0;
    slip_direction_hex[8](3, 0) = 6.0 / 3.0;

    slip_system_id_[9] = "pyramidal_4";

    slip_plane_normal_hex[9](0, 0) = -1.0;
    slip_plane_normal_hex[9](1, 0) = -1.0;
    slip_plane_normal_hex[9](2, 0) = 2.0;
    slip_plane_normal_hex[9](3, 0) = 1.0;

    slip_direction_hex[9](0, 0) = 1.0 / 3.0;
    slip_direction_hex[9](1, 0) = 1.0 / 3.0;
    slip_direction_hex[9](2, 0) = -2.0 / 3.0;
    slip_direction_hex[9](3, 0) = 6.0 / 3.0;

    slip_system_id_[10] = "pyramidal_5";

    slip_plane_normal_hex[10](0, 0) = -1.0;
    slip_plane_normal_hex[10](1, 0) = 2.0;
    slip_plane_normal_hex[10](2, 0) = -1.0;
    slip_plane_normal_hex[10](3, 0) = 1.0;

    slip_direction_hex[10](0, 0) = 1.0 / 3.0;
    slip_direction_hex[10](1, 0) = -2.0 / 3.0;
    slip_direction_hex[10](2, 0) = 1.0 / 3.0;
    slip_direction_hex[10](3, 0) = 6.0 / 3.0;

    slip_system_id_[11] = "pyramidal_6";

    slip_plane_normal_hex[11](0, 0) = 2.0;
    slip_plane_normal_hex[11](1, 0) = -1.0;
    slip_plane_normal_hex[11](2, 0) = -1.0;
    slip_plane_normal_hex[11](3, 0) = 1.0;

    slip_direction_hex[11](0, 0) = -2.0 / 3.0;
    slip_direction_hex[11](1, 0) = 1.0 / 3.0;
    slip_direction_hex[11](2, 0) = 1.0 / 3.0;
    slip_direction_hex[11](3, 0) = 6.0 / 3.0;

    // transform Miller-Bravais index notation of hexagonal lattices to Miller index notation of
    // cubic lattices
    MAT::CrystalPlasticity::MillerBravaisToMiller(slip_plane_normal_hex, slip_direction_hex);
  }
  else if (lattice_type_ == "HCP" or lattice_type_ == "BCC" or lattice_type_ == "FCC")
  {
    dserror("Sorry, %s lattices are not yet supported", lattice_type_.c_str());
  }
  // academic test lattice containing only one slip system
  else if (lattice_type_ == "TEST")
  {
    slip_system_id_[0] = "test_1";

    slip_plane_normal_[0](0, 0) = 1.0;
    slip_plane_normal_[0](1, 0) = 0.0;
    slip_plane_normal_[0](2, 0) = 1.0;

    slip_direction_[0](0, 0) = 1.0;
    slip_direction_[0](1, 0) = 0.0;
    slip_direction_[0](2, 0) = -1.0;
  }
  else
  {
    dserror(
        "Lattice type not known. Please check the LAT parameter in the input file. Currently it "
        "has to be FCC, BCC, HCP, D019 or L10");
  }


  // test whether directions and normals are perpendicular
  std::vector<LINALG::Matrix<1, 1>> normality_test;
  normality_test.resize(slip_system_count_);

  for (int i = 0; i < slip_system_count_; i++)
  {
    normality_test[i].MultiplyTN(slip_direction_[i], slip_plane_normal_[i]);
    if (normality_test[i](0, 0) != 0.0)
    {
      dserror(
          "Warning, slip direction and slip plane normal of slip system "
          "%d are not perpendicular!",
          i + 1);
    }
  }

  // determine magnitude of Burgers vector
  burgers_vector_mag_.resize(slip_system_count_);
  LINALG::Matrix<3, 1> tmp_scale;

  // import lattice constant from user input
  lattice_constant_ = params_->lattice_const_;

  for (int i = 0; i < slip_system_count_; i++)
  {
    tmp_scale = slip_direction_[i];
    tmp_scale.Scale(lattice_constant_);
    burgers_vector_mag_[i] = tmp_scale.Norm2();
  }

  // normalize slip plane normals
  for (int i = 0; i < slip_system_count_; i++)
  {
    slip_plane_normal_[i].Scale(1.0 / slip_plane_normal_[i].Norm2());
  }

  // normalize slip directions
  for (int i = 0; i < slip_system_count_; i++)
  {
    slip_direction_[i].Scale(1.0 / slip_direction_[i].Norm2());
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
      lattice_orientation_(i, 0) = fiber1.at(i);
      lattice_orientation_(i, 1) = fiber2.at(i);
      lattice_orientation_(i, 2) = fiber3.at(i);
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
    std::vector<LINALG::Matrix<4, 1>>& slip_plane_normal_hex,
    std::vector<LINALG::Matrix<4, 1>>& slip_direction_hex)
{
  for (int unsigned i = 0; i < slip_plane_normal_hex.size(); i++)
  {
    slip_plane_normal_[i](0, 0) = slip_plane_normal_hex.at(i)(0, 0);
    slip_plane_normal_[i](1, 0) =
        (slip_plane_normal_hex.at(i)(0, 0) + 2.0 * slip_plane_normal_hex.at(i)(1, 0)) / sqrt(3.0);
    slip_plane_normal_[i](2, 0) = slip_plane_normal_hex.at(i)(3, 0) / c_to_a_ratio_;

    slip_direction_[i](0, 0) = 1.5 * slip_direction_hex.at(i)(0, 0);
    slip_direction_[i](1, 0) =
        sqrt(3.0) * (0.5 * slip_direction_hex.at(i)(0, 0) + slip_direction_hex.at(i)(1, 0));
    slip_direction_[i](2, 0) = slip_direction_hex.at(i)(3, 0) * c_to_a_ratio_;
  }
  return;
}

/*---------------------------------------------------------------------------------*
 | local Newton-Raphson iteration                                                  |
 *---------------------------------------------------------------------------------*/
void MAT::CrystalPlasticity::NewtonRaphson(LINALG::Matrix<3, 3>& deform_grad,
    std::vector<double>& gamma_result, std::vector<double>& defect_density_result,
    LINALG::Matrix<3, 3>& second_pk_stress_result, LINALG::Matrix<3, 3>& plastic_deform_grad_result)
{
  // initialize iteration
  // ----------------------------------
  // max. number of iterations
  const int max_iterations = 20;
  // iteration counter
  int iteration_number = 0;
  // convergence indicator
  bool converged = false;

  // elastic predictor
  //--------------------------------------------------------------------------
  // assume load step is elastic, i.e., there is no change in plastic shear
  std::vector<double> gamma_trial = gamma_last_->at(gp);

  // trial values of internal variables
  std::vector<double> defect_density_trial(slip_system_count_);
  LINALG::Matrix<3, 3> second_pk_stress_trial;
  LINALG::Matrix<3, 3> plastic_deform_grad_trial;

  // residuals
  std::vector<double> residuals_trial(slip_system_count_);
  double total_residual;

  // empty matrix
  LINALG::Matrix<3, 3> emptymat(true);

  // iteration
  //--------------------------------------------------------------------------
  while (!converged && iteration_number++ < max_iterations)
  {
    // reset iteration relevant values
    total_residual = 0.0;
    second_pk_stress_trial = emptymat;

    // set up flow rule with given deformation gradient F and trial plastic shears gamma_trial
    this->SetupFlowRule(deform_grad, gamma_trial, plastic_deform_grad_trial, defect_density_trial,
        second_pk_stress_trial, residuals_trial);

    // determine total residuum as norm of the vector of slip system residuals
    for (int i = 0; i < slip_system_count_; i++)
    {
      total_residual += pow(residuals_trial.at(i), 2.0);
    }
    total_residual = sqrt(total_residual);

    // convergence check
    if (total_residual < newton_tolerance_)
    {
      converged = true;

      // collect results
      //--------------------------------------------------------------------------
      gamma_result = gamma_trial;
      defect_density_result = defect_density_trial;
      second_pk_stress_result = second_pk_stress_trial;
      plastic_deform_grad_result = plastic_deform_grad_trial;
    }
    else
    {
      // update rule: gamma_trial(n+1) = gamma_trial(n) + d_gamma_trial(n)
      // therefore the Jacobean ResidualStiffness(j,i)=dr(j)/d_gamma_trial(i) is set up, where r is
      // the residuum solving ResidualStiffness * d_gamma_trial(n) = - r

      // for L10 or D019 lattices
      LINALG::Matrix<12, 12> residual_tangent_12x12;
      LINALG::Matrix<12, 1> d_gamma_trial_12(true);
      LINALG::Matrix<12, 1> residuals_trial_LIN_12;
      LINALG::FixedSizeSerialDenseSolver<12, 12, 1> newton_raphson_solver_12x12;
      // for TEST lattices
      LINALG::Matrix<1, 1> residual_tangent_1x1;
      LINALG::Matrix<1, 1> d_gamma_trial_1(true);
      LINALG::Matrix<1, 1> residuals_trial_LIN_1;
      LINALG::FixedSizeSerialDenseSolver<1, 1, 1> newton_raphson_solver_1x1;

      // differentiate Residuum with respect to gamma_trial by perturbation DEFAULT:  gamma_epsilon
      // = 1.0e-9
      double gamma_epsilon = 1.0e-9;

      // change of residuals.at(j) with perturbation of gamma_trial_.at(i)

      // perturb single shears
      for (int i = 0; i < slip_system_count_; i++)
      {
        // perturbed vector of plastic shears
        std::vector<double> gamma_perturbed(slip_system_count_);
        // resultant vector of residuals
        std::vector<double> residuals_perturbed(slip_system_count_);
        // resultant vector of defect densities
        std::vector<double> defect_densities_perturbed(slip_system_count_);
        // resultant 2nd PK stress
        LINALG::Matrix<3, 3> second_pk_stress_perturbed(true);
        // resultant plastic part of deformation gradient
        LINALG::Matrix<3, 3> plastic_deform_grad_perturbed;

        // initially unperturbed
        gamma_perturbed = gamma_trial;

        // perturbation of single gammas
        gamma_perturbed.at(i) += gamma_epsilon;

        // evaluate flow rule for perturbed gamma
        this->SetupFlowRule(deform_grad, gamma_perturbed, plastic_deform_grad_perturbed,
            defect_densities_perturbed, second_pk_stress_perturbed, residuals_perturbed);

        // finite difference
        if (lattice_type_ == "L10" || lattice_type_ == "D019")
        {
          for (int j = 0; j < slip_system_count_; j++)
          {
            residual_tangent_12x12(j, i) =
                (residuals_perturbed.at(j) - residuals_trial.at(j)) / gamma_epsilon;

            residuals_trial_LIN_12(i) = -residuals_trial.at(i);
          }
        }
        else if (lattice_type_ == "TEST")
        {
          for (int j = 0; j < slip_system_count_; j++)
          {
            residual_tangent_1x1(j, i) =
                (residuals_perturbed.at(j) - residuals_trial.at(j)) / gamma_epsilon;

            residuals_trial_LIN_1(i) = -residuals_trial.at(i);
          }
        }
      }

      // solve resultant system of equations
      if (lattice_type_ == "L10" || lattice_type_ == "D019")
      {
        newton_raphson_solver_12x12.SetMatrix(residual_tangent_12x12);
        newton_raphson_solver_12x12.SetVectors(d_gamma_trial_12, residuals_trial_LIN_12);
        newton_raphson_solver_12x12.Solve();

        // update vector of plastic shears to new trial value
        for (int i = 0; i < slip_system_count_; i++)
        {
          gamma_trial.at(i) += d_gamma_trial_12(i);
        }
      }
      else if (lattice_type_ == "TEST")
      {
        newton_raphson_solver_1x1.SetMatrix(residual_tangent_1x1);
        newton_raphson_solver_1x1.SetVectors(d_gamma_trial_1, residuals_trial_LIN_1);
        newton_raphson_solver_1x1.Solve();

        // update vector of plastic shears to new trial value
        for (int i = 0; i < slip_system_count_; i++)
        {
          gamma_trial.at(i) += d_gamma_trial_1(i);
        }
      }
    }
  }  // while

  // check whether or not the result converged during max_iterations iterations.
  // TODO: If not: reduce time step and restart time increment!
  if (!converged)
  {
    dserror("Internal Newton Raphson did not converge within max iterations");
  }

  return;
}  // NewtonRaphson


/*---------------------------------------------------------------------------------*
 | Evaluate flow rule for a  given increment of plastic shear and set up residuum  |
 *---------------------------------------------------------------------------------*/

void MAT::CrystalPlasticity::SetupFlowRule(LINALG::Matrix<3, 3> deform_grad,
    std::vector<double> gamma_trial, LINALG::Matrix<3, 3>& plastic_deform_grad_trial,
    std::vector<double>& defect_density_trial, LINALG::Matrix<3, 3>& second_pk_stress_trial,
    std::vector<double>& residuals_trial)
{
  // set up trial increment of plastic shear for each slip system
  std::vector<double> delta_gamma_trial(slip_system_count_);

  for (int i = 0; i < slip_system_count_; i++)
  {
    delta_gamma_trial[i] = gamma_trial[i] - (*gamma_last_)[gp][i];
  }

  // determine trial defect densities that would result from delta_gamma_trial
  //--------------------------------------------------------------------------

  double total_dislocation_density = 0.0;

  for (int i = 0; i < slip_system_count_; i++)
  {
    // dislocation generation
    //  TODO EQUATION
    defect_density_trial[i] =
        (*defect_densities_last_)[gp][i] +
        (1.0 / (burgers_vector_mag_[i] * dislocation_generation_coeff_[subset_index_[i] - 1])) *
            sqrt((*defect_densities_last_)[gp][i]) * abs(delta_gamma_trial[i]);

    // dynamic recovery
    defect_density_trial[i] -= dislocation_dyn_recovery_coeff_[subset_index_[i] - 1] *
                               (*defect_densities_last_)[gp][i] * abs(delta_gamma_trial[i]);

    // determine updated total defect densities
    total_dislocation_density += defect_density_trial[i];
  }

  // kinematics and plastic velocity gradient
  //--------------------------------------------------------------------------

  // determine trial plastic velocity gradient L_p
  // set up L_p_trial
  LINALG::Matrix<3, 3> plastic_velocity_grad_trial(true);
  LINALG::Matrix<3, 3> temp_mat;

  for (int i = 0; i < slip_system_count_; i++)
  {
    temp_mat.MultiplyNT(delta_gamma_trial.at(i), slip_direction_[i], slip_plane_normal_[i]);
    plastic_velocity_grad_trial.Update(plastic_velocity_grad_trial, temp_mat);
  }

  // take unimodular part of I + L_p to ensure plastic incompressibility
  LINALG::Matrix<3, 3> unimod_identity_plus_plastic_velocity_grad_trial(true);

  unimod_identity_plus_plastic_velocity_grad_trial.Update(identity3, plastic_velocity_grad_trial);
  unimod_identity_plus_plastic_velocity_grad_trial.Scale(
      pow(unimod_identity_plus_plastic_velocity_grad_trial.Determinant(), -1.0 / 3.0));

  // determine trial plastic deformation gradient
  plastic_deform_grad_trial.MultiplyNN(
      unimod_identity_plus_plastic_velocity_grad_trial, (*plastic_deform_grad_last_)[gp]);

  // determine trial stress
  //--------------------------------------------------------------------------
  // determine trial elastic deformation gradient
  // get the inverse FP^{-1}
  LINALG::Matrix<3, 3> inv_plastic_deform_grad_trial;
  inv_plastic_deform_grad_trial.Invert(plastic_deform_grad_trial);
  LINALG::Matrix<3, 3> elastic_deform_grad_trial;
  elastic_deform_grad_trial.MultiplyNN(deform_grad, inv_plastic_deform_grad_trial);

  // calculate the Jacobi-determinant J = det(FE_{n+1}) and the logarithm of it
  double jacobi_det_trial = elastic_deform_grad_trial.Determinant();
  double ln_jacobi_det_trial = log(jacobi_det_trial);

  // set up elastic right cauchy green and its inverse
  LINALG::Matrix<3, 3> elastic_right_cauchy_green;
  elastic_right_cauchy_green.MultiplyTN(elastic_deform_grad_trial, elastic_deform_grad_trial);
  LINALG::Matrix<3, 3> inv_elastic_right_cauchy_green;
  inv_elastic_right_cauchy_green.Invert(elastic_right_cauchy_green);

  // 2nd Piola-Kirchhoff stress
  // S = lambda * ln_jacobi_det_trial * inv_elastic_right_cauchy_green + mu * (identity3
  // - inv_elastic_right_cauchy_green)

  second_pk_stress_trial.Update(lambda_ * ln_jacobi_det_trial, inv_elastic_right_cauchy_green, 1.0);
  second_pk_stress_trial.Update(mue_, identity3, 1.0);
  second_pk_stress_trial.Update(-mue_, inv_elastic_right_cauchy_green, 1.0);

  // Mandel stress
  LINALG::Matrix<3, 3> mandel_stress_trial;
  mandel_stress_trial.MultiplyNN(elastic_right_cauchy_green, second_pk_stress_trial);

  // setting up Residua with flow/creep rule
  //--------------------------------------------------------------------------

  for (int i = 0; i < slip_system_count_; i++)
  {
    // resolved shear stress/Schmid stress
    LINALG::Matrix<1, 1> resolved_shear_stress(true);

    LINALG::Matrix<3, 1> TempVec(true);

    TempVec.MultiplyNN(mandel_stress_trial, slip_plane_normal_[i]);

    resolved_shear_stress.MultiplyTN(slip_direction_[i], TempVec);

    // work hardening

    // lattice resistance
    double tau_y0 = tau_y_0_[subset_index_[i] - 1];

    // Hall-Petch strengthening term
    double hall_petch_strengthening = hall_petch_coeffs_[subset_index_[i] - 1] *
                                      (1.0 / sqrt(micro_boundary_distance_[subset_index_[i] - 1]));

    // work hardening hardening increment due to accumulated dislocation density
    double delta_tau_y = 0.5 * mue_ * burgers_vector_mag_[i] * sqrt(total_dislocation_density);

    // slip system strength, i.e. tau_y = tau_y0 + hall_petch_strengthening + delta_tau_y
    double tau_y = tau_y0 + hall_petch_strengthening + delta_tau_y;
    // check for consitency
    if (tau_y < 0.0) dserror("Negative slip systems strength! Please check your input");

    // shear rate as determined from the flow rule
    double gamma_dot = 0.0;

    gamma_dot = abs(resolved_shear_stress(0, 0)) / tau_y;
    gamma_dot = pow(gamma_dot, 50.0);  // TODO EXPONENT AUS DAT FILE HOLEN
    gamma_dot = gamma_dot * gamma_dot_0_[subset_index_[i] - 1];
    gamma_dot = std::copysign(gamma_dot, resolved_shear_stress(0, 0));

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
  for (int sys = 0; sys < slip_system_count_; sys++)
  {
    ID = "plastic_shear_";
    ID += slip_system_id_[sys];
    names[ID] = 1;  // scalar
  }

  // accumulated plastic shears
  ID = "accumulated_plastic_shear";
  names[ID] = 1;  // scalar

  // dislocation densities on single systems
  for (int sys = 0; sys < slip_system_count_; sys++)
  {
    ID = "dislocation_density_";
    ID += slip_system_id_[sys];
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
  for (int sys = 0; sys < slip_system_count_; sys++)
  {
    if (name == "plastic_shear_" + slip_system_id_[sys])
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
      for (int sys = 0; sys < slip_system_count_; sys++)
      {
        gamma_acc += abs((*gamma_last_)[gp_index][sys]);
      }
      data[0] = gamma_acc;
    }
  }

  // dislocation densities
  for (int sys = 0; sys < slip_system_count_; sys++)
  {
    if (name == "dislocation_density_" + slip_system_id_[sys])
    {
      if ((int)data.size() != 1) dserror("size mismatch");
      for (int gp_index = 0; gp_index < numgp; gp_index++)
      {
        data[0] = (*defect_densities_last_)[gp_index][sys];
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
      for (int sys = 0; sys < slip_system_count_; sys++)
      {
        dis_dens_acc += (*defect_densities_last_)[gp_index][sys];
      }
      data[0] = dis_dens_acc;
    }
  }

  return true;
}  // VisData()
