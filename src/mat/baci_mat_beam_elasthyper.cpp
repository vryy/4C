/*-----------------------------------------------------------------------------------------------*/
/*! \file
\brief constitutive relations for beam cross-section resultants (hyperelastic stored energy
function)


\level 3
*/
/*-----------------------------------------------------------------------------------------------*/

#include "baci_mat_beam_elasthyper.H"

#include "baci_lib_globalproblem.H"
#include "baci_mat_beam_elasthyper_parameter.H"
#include "baci_mat_par_bundle.H"

#include <Sacado.hpp>


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
MAT::BeamElastHyperMaterialType<T> MAT::BeamElastHyperMaterialType<T>::instance_;



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
DRT::ParObject* MAT::BeamElastHyperMaterialType<T>::Create(const std::vector<char>& data)
{
  MAT::Material* matobject = new MAT::BeamElastHyperMaterial<T>();
  matobject->Unpack(data);
  return matobject;
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
MAT::BeamElastHyperMaterial<T>::BeamElastHyperMaterial() : params_(nullptr)
{
  // empty constructor
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
MAT::BeamElastHyperMaterial<T>::BeamElastHyperMaterial(
    MAT::PAR::BeamElastHyperMaterialParameterGeneric* params)
    : params_(params)
{
  // empty constructor
}



/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::EvaluateForceContributionsToStress(
    CORE::LINALG::Matrix<3, 1, T>& stressN, const CORE::LINALG::Matrix<3, 3, T>& CN,
    const CORE::LINALG::Matrix<3, 1, T>& Gamma, const unsigned int gp)
{
  // compute material stresses by multiplying strains with constitutive matrix
  stressN.Multiply(CN, Gamma);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::EvaluateMomentContributionsToStress(
    CORE::LINALG::Matrix<3, 1, T>& stressM, const CORE::LINALG::Matrix<3, 3, T>& CM,
    const CORE::LINALG::Matrix<3, 1, T>& Cur, const unsigned int gp)
{
  // compute material stresses by multiplying curvature with constitutive matrix
  stressM.Multiply(CM, Cur);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::ComputeConstitutiveParameter(
    CORE::LINALG::Matrix<3, 3, T>& C_N, CORE::LINALG::Matrix<3, 3, T>& C_M)
{
  // setup constitutive matrices
  MAT::BeamElastHyperMaterial<T>::GetConstitutiveMatrixOfForcesMaterialFrame(C_N);
  MAT::BeamElastHyperMaterial<T>::GetConstitutiveMatrixOfMomentsMaterialFrame(C_M);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  this->AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  this->AddtoPack(data, matid);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  this->ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror(
        "Wrong instance type data. The extracted type id is %d, while the UniqueParObjectId is %d",
        type, UniqueParObjectId());

  // matid and recover params_
  int matid;
  this->ExtractfromPack(position, data, matid);
  params_ = nullptr;

  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();

      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);

      /* the idea is that we have a generic type of material (this class), but various
       * possible sets of material parameters to 'feed' these very general constitutive relations */
      if (mat->Type() == INPAR::MAT::m_beam_reissner_elast_hyper or
          mat->Type() == INPAR::MAT::m_beam_reissner_elast_hyper_bymodes or
          mat->Type() == INPAR::MAT::m_beam_kirchhoff_elast_hyper or
          mat->Type() == INPAR::MAT::m_beam_reissner_elast_plastic or
          mat->Type() == INPAR::MAT::m_beam_kirchhoff_elast_hyper_bymodes or
          mat->Type() == INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper or
          mat->Type() == INPAR::MAT::m_beam_kirchhoff_torsionfree_elast_hyper_bymodes)
        params_ = static_cast<MAT::PAR::BeamElastHyperMaterialParameterGeneric*>(mat);
      else
        dserror("Type of material parameter %d does not fit to type of material law %d",
            mat->Type(), MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
MAT::PAR::Parameter* MAT::BeamElastHyperMaterial<T>::Parameter() const
{
  return params_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
const MAT::PAR::BeamElastHyperMaterialParameterGeneric& MAT::BeamElastHyperMaterial<T>::Params()
    const
{
  if (params_ == nullptr) dserror("pointer to parameter class is not set!");

  return *params_;
}


/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::GetConstitutiveMatrixOfForcesMaterialFrame(
    CORE::LINALG::Matrix<3, 3, T>& C_N) const
{
  // defining material constitutive matrix CN between Gamma and N
  // according to Jelenic 1999, section 2.4
  C_N.Clear();

  C_N(0, 0) = Params().GetAxialRigidity();
  C_N(1, 1) = Params().GetShearRigidity2();
  C_N(2, 2) = Params().GetShearRigidity3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::GetConstitutiveMatrixOfMomentsMaterialFrame(
    CORE::LINALG::Matrix<3, 3, T>& C_M) const
{
  // defining material constitutive matrix CM between curvature and moment
  // according to Jelenic 1999, section 2.4
  C_M.Clear();

  C_M(0, 0) = Params().GetTorsionalRigidity();
  C_M(1, 1) = Params().GetBendingRigidity2();
  C_M(2, 2) = Params().GetBendingRigidity3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
double MAT::BeamElastHyperMaterial<T>::GetTranslationalMassInertiaFactor() const
{
  return Params().GetTranslationalMassInertia();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::GetMassMomentOfInertiaTensorMaterialFrame(
    CORE::LINALG::Matrix<3, 3>& J) const
{
  J.Clear();

  J(0, 0) = Params().GetPolarMassMomentOfInertia();
  J(1, 1) = Params().GetMassMomentOfInertia2();
  J(2, 2) = Params().GetMassMomentOfInertia3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::GetMassMomentOfInertiaTensorMaterialFrame(
    CORE::LINALG::Matrix<3, 3, Sacado::Fad::DFad<double>>& J) const
{
  J.Clear();

  J(0, 0) = Params().GetPolarMassMomentOfInertia();
  J(1, 1) = Params().GetMassMomentOfInertia2();
  J(2, 2) = Params().GetMassMomentOfInertia3();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
double MAT::BeamElastHyperMaterial<T>::GetInteractionRadius() const
{
  return this->Params().GetInteractionRadius();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::GetStiffnessMatrixOfMoments(
    CORE::LINALG::Matrix<3, 3, T>& stiffness_matrix, const CORE::LINALG::Matrix<3, 3, T>& C_M,
    const int gp)
{
  stiffness_matrix = C_M;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
template <typename T>
void MAT::BeamElastHyperMaterial<T>::GetStiffnessMatrixOfForces(
    CORE::LINALG::Matrix<3, 3, T>& stiffness_matrix, const CORE::LINALG::Matrix<3, 3, T>& C_N,
    const int gp)
{
  stiffness_matrix = C_N;
}

// explicit template instantiations
template class MAT::BeamElastHyperMaterial<double>;
template class MAT::BeamElastHyperMaterial<Sacado::Fad::DFad<double>>;


template class MAT::BeamElastHyperMaterialType<double>;
template class MAT::BeamElastHyperMaterialType<Sacado::Fad::DFad<double>>;
