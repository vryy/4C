/*----------------------------------------------------------------------*/
/*! \file

\brief main file containing routines for calculation of solid element
       simple displacement based
\level 1
*/
/*----------------------------------------------------------------------*/

#include "solid_ele_calc.H"
#include <Teuchos_ParameterList.hpp>
#include <memory>
#include <optional>
#include "drt_utils_integration.H"
#include "drt_utils.H"
#include "voigt_notation.H"
#include "solid_ele.H"
#include "drt_discret.H"
#include "so3_material.H"
#include "solid_ele_calc_lib.H"
#include "solid_utils.H"
#include "drt_utils_local_connectivity_matrices.H"
#include "drt_fiber_node.H"
#include "drt_fiber_utils.H"
#include "nodal_fiber_holder.H"

#include "gauss_point_data_output_manager.H"
#include "so_element_service.H"

namespace
{
  template <DRT::Element::DiscretizationType distype>
  std::unique_ptr<DRT::UTILS::GaussIntegration> CreateDefaultGaussIntegration()
  {
    return CreateDefaultGaussIntegration<distype>(
        DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);
  }

  template <DRT::Element::DiscretizationType distype>
  std::unique_ptr<DRT::UTILS::GaussIntegration> CreateGaussIntegration(
      const DRT::UTILS::IntPointsAndWeights<DRT::UTILS::DisTypeToDim<distype>::dim>& intpoints)
  {
    // format as DRT::UTILS::GaussIntegration
    Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> gp =
        Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints);

    std::array<double, 3> xi = {0., 0., 0.};
    for (int i = 0; i < intpoints.IP().nquad; ++i)
    {
      for (int d = 0; d < DRT::UTILS::DisTypeToDim<distype>::dim; ++d)
        xi[d] = intpoints.IP().qxg[i][d];
      gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
    }

    return std::make_unique<DRT::UTILS::GaussIntegration>(gp);
  }

  std::vector<char>& GetMutableStressData(
      const DRT::ELEMENTS::Solid& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.ParamsInterface().MutableStressDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("stress");
    }
  }

  std::vector<char>& GetMutableStrainData(
      const DRT::ELEMENTS::Solid& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return *ele.ParamsInterface().MutableStrainDataPtr();
    }
    else
    {
      return *params.get<Teuchos::RCP<std::vector<char>>>("strain");
    }
  }

  INPAR::STR::StressType GetIOStressType(
      const DRT::ELEMENTS::Solid& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.ParamsInterface().GetStressOutputType();
    }
    else
    {
      return DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress");
    }
  }

  INPAR::STR::StrainType GetIOStrainType(
      const DRT::ELEMENTS::Solid& ele, const Teuchos::ParameterList& params)
  {
    if (ele.IsParamsInterface())
    {
      return ele.ParamsInterface().GetStrainOutputType();
    }
    else
    {
      return DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain");
    }
  }
}  // namespace

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleCalc<distype>* DRT::ELEMENTS::SolidEleCalc<distype>::Instance(bool create)
{
  static SolidEleCalc<distype>* instance;
  if (create)
  {
    if (!instance) instance = new SolidEleCalc<distype>();
  }
  else
  {
    if (instance) delete instance;
    instance = nullptr;
  }
  return instance;
}

template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::SolidEleCalc<distype>::SolidEleCalc()
    : DRT::ELEMENTS::SolidEleInterface::SolidEleInterface()
{
  InitializeDefaultQuadrature();
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::InitializeDefaultQuadrature()
{
  DRT::UTILS::GaussRule3D rule = DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule;
  if (distype == DRT::Element::DiscretizationType::tet10)
  {
    // HACK: Don't use the "OptGaussRule", but use only 4 Gauss points
    // TODO: for tet-elements we need two different Gauss rules (4 points for force and stiffness
    // matrix, 11 points for mass matrix)
    rule = DRT::UTILS::GaussRule3D::tet_4point;
  }
  else if (distype == DRT::Element::DiscretizationType::tet4)
  {
    rule = DRT::UTILS::GaussRule3D::tet_1point;
  }
  // setup default integration
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(rule);

  // format as DRT::UTILS::GaussIntegration
  Teuchos::RCP<DRT::UTILS::CollectedGaussPoints> gp =
      Teuchos::rcp(new DRT::UTILS::CollectedGaussPoints);

  std::array<double, 3> xi = {0., 0., 0.};
  for (int i = 0; i < intpoints.IP().nquad; ++i)
  {
    for (int d = 0; d < nsd_; ++d) xi[d] = intpoints.IP().qxg[i][d];
    gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
  }

  // save default integration rule
  default_integration_ = Teuchos::rcp(new DRT::UTILS::GaussIntegration(gp));
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateNonlinearForceStiffnessMass(
    const DRT::ELEMENTS::Solid& ele, const DRT::Discretization& discretization,
    const std::vector<int>& lm, Teuchos::ParameterList& params,
    Epetra_SerialDenseVector* force_vector, Epetra_SerialDenseMatrix* stiffness_matrix,
    Epetra_SerialDenseMatrix* mass_matrix)
{
  std::unique_ptr<LINALG::Matrix<nsd_ * nen_, nsd_* nen_>> stiff = nullptr;
  std::unique_ptr<LINALG::Matrix<nsd_ * nen_, nsd_* nen_>> mass = nullptr;
  std::unique_ptr<LINALG::Matrix<nsd_ * nen_, 1>> force = nullptr;
  if (stiffness_matrix != nullptr)
    stiff = std::make_unique<LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>>(*stiffness_matrix, true);
  if (mass_matrix != nullptr)
    mass = std::make_unique<LINALG::Matrix<nsd_ * nen_, nsd_ * nen_>>(*mass_matrix, true);
  if (force_vector != nullptr)
    force = std::make_unique<LINALG::Matrix<nsd_ * nen_, 1>>(*force_vector, true);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // Loop over all Gauss points
  for (int gp = 0; gp < default_integration_->NumPoints(); ++gp)
  {
    const LINALG::Matrix<nsd<distype>, 1> xi =
        EvaluateParameterCoordinate<distype>(*default_integration_, gp);

    const ShapeFunctionsAndDerivatives<distype> shape_functions =
        EvaluateShapeFunctionsAndDerivs<distype>(xi);

    const JacobianMapping<distype> jacobian_mapping =
        EvaluateJacobianMapping(shape_functions, nodal_coordinates, *default_integration_, gp);

    const Strains<distype> strains = EvaluateStrains<distype>(nodal_coordinates, jacobian_mapping);

    LINALG::Matrix<numstr<distype>, nsd<distype> * nen<distype>> Bop =
        EvaluateStrainGradient(jacobian_mapping, strains);

    const Stress<distype> stress =
        EvaluateMaterialStress(*ele.SolidMaterial(), strains, params, gp, ele.Id());

    if (force != nullptr)
      AddInternalForceVector(Bop, stress, jacobian_mapping.integration_factor_, *force);

    if (stiff != nullptr)
    {
      AddElasticStiffnessMatrix(Bop, stress, jacobian_mapping.integration_factor_, *stiff);
      AddGeometricStiffnessMatrix(jacobian_mapping, stress, *stiff);
    }

    if (mass != nullptr)
    {
      AddMassMatrix(shape_functions, jacobian_mapping.integration_factor_,
          ele.SolidMaterial()->Density(gp), *mass);
    }
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Recover(const DRT::ELEMENTS::Solid& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  // nothing to do for a standard element
  // except...
  // TODO: We need to recover history information of materials!
  // which was also not implemented in the old elements
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Update(const DRT::ELEMENTS::Solid& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // Loop over all Gauss points
  for (int gp = 0; gp < default_integration_->NumPoints(); ++gp)
  {
    const LINALG::Matrix<nsd<distype>, 1> xi =
        EvaluateParameterCoordinate<distype>(*default_integration_, gp);

    const ShapeFunctionsAndDerivatives<distype> shape_functions =
        EvaluateShapeFunctionsAndDerivs<distype>(xi);

    const JacobianMapping<distype> jacobian_mapping =
        EvaluateJacobianMapping(shape_functions, nodal_coordinates, *default_integration_, gp);

    const Strains<distype> strains = EvaluateStrains<distype>(nodal_coordinates, jacobian_mapping);

    ele.SolidMaterial()->Update(strains.defgrd_, gp, params, ele.Id());
  }  // gp loop
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::PostProcessStressStrain(const DRT::ELEMENTS::Solid& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  // TODO: This method would be much easier if we get rid of post_drt_*
  const Epetra_SerialDenseMatrix gpstress =
      *(*params.get<Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>>>(
          "gpstressmap"))[ele.Id()];

  std::string stresstype = params.get<std::string>("stresstype");
  Epetra_MultiVector& poststress = *params.get<Teuchos::RCP<Epetra_MultiVector>>("poststress");

  if (stresstype == "ndxyz")
  {
    ExtrapolateGPQuantityToNodesAndAssemble<distype>(
        ele, gpstress, poststress, true, *default_integration_);
  }
  else if (stresstype == "cxyz")
  {
    const Epetra_BlockMap& elemap = poststress.Map();
    int lid = elemap.LID(ele.Id());
    if (lid != -1)
    {
      for (unsigned i = 0; i < numstr_; ++i)
      {
        double& s = (*((poststress)(i)))[lid];  // resolve pointer for faster access
        s = 0.;
        for (int j = 0; j < gpstress.M(); ++j)
        {
          s += gpstress(j, i);
        }
        s *= 1.0 / gpstress.M();
      }
    }
  }
  else
    dserror("unknown type of stress/strain output on element level");
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::CalculateStress(const DRT::ELEMENTS::Solid& ele,
    const DRT::Discretization& discretization, const std::vector<int>& lm,
    Teuchos::ParameterList& params)
{
  // TODO: If we get rid of post_drt_*, we don't need this here anymore. We could directly use
  // InitializeGaussPointDataOutput and EvaluateGaussPointDataOutput and write the stresses there.
  if (discretization.Comm().MyPID() != ele.Owner()) return;

  std::vector<char>& serialized_stress_data = GetMutableStressData(ele, params);
  std::vector<char>& serialized_strain_data = GetMutableStrainData(ele, params);
  Epetra_SerialDenseMatrix stress_data(default_integration_->NumPoints(), numstr_);
  Epetra_SerialDenseMatrix strain_data(default_integration_->NumPoints(), numstr_);

  const NodalCoordinates<distype> nodal_coordinates =
      EvaluateNodalCoordinates<distype>(ele, discretization, lm);

  // Loop over all Gauss points
  for (int gp = 0; gp < default_integration_->NumPoints(); ++gp)
  {
    const LINALG::Matrix<nsd<distype>, 1> xi =
        EvaluateParameterCoordinate<distype>(*default_integration_, gp);

    const ShapeFunctionsAndDerivatives<distype> shape_functions =
        EvaluateShapeFunctionsAndDerivs<distype>(xi);

    const JacobianMapping<distype> jacobian_mapping =
        EvaluateJacobianMapping(shape_functions, nodal_coordinates, *default_integration_, gp);

    const Strains<distype> strains = EvaluateStrains<distype>(nodal_coordinates, jacobian_mapping);

    LINALG::Matrix<numstr<distype>, nsd<distype> * nen<distype>> Bop =
        EvaluateStrainGradient(jacobian_mapping, strains);

    const Stress<distype> stress =
        EvaluateMaterialStress(*ele.SolidMaterial(), strains, params, gp, ele.Id());

    AssembleStrainTypeToMatrixRow(strains, GetIOStrainType(ele, params), strain_data, gp);
    AssembleStressTypeToMatrixRow(strains, stress, GetIOStressType(ele, params), stress_data, gp);
  }  // gp loop

  Serialize(stress_data, serialized_stress_data);
  Serialize(strain_data, serialized_strain_data);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::Setup(
    const DRT::ELEMENTS::Solid& ele, DRT::INPUT::LineDefinition* linedef)
{
  ele.SolidMaterial()->Setup(default_integration_->NumPoints(), linedef);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::MaterialPostSetup(const DRT::ELEMENTS::Solid& ele)
{
  Teuchos::ParameterList params{};
  if (DRT::FIBER::UTILS::HaveNodalFibers<distype>(ele.Nodes()))
  {
    // This element has fiber nodes.
    // Interpolate fibers to the Gauss points and pass them to the material

    // Get shape functions
    const static std::vector<LINALG::Matrix<nen_, 1>> shapefcts = std::invoke(
        [&]
        {
          const DRT::UTILS::GaussIntegration& integration_rule = *default_integration_;

          std::vector<LINALG::Matrix<nen_, 1>> shapefcns(integration_rule.NumPoints());
          for (int gp = 0; gp < integration_rule.NumPoints(); ++gp)
          {
            LINALG::Matrix<nsd_, 1> xi(integration_rule.Point(gp), true);
            DRT::UTILS::shape_function<distype>(xi, shapefcns[gp]);
          }
          return shapefcns;
        });

    // add fibers to the ParameterList
    DRT::FIBER::NodalFiberHolder fiberHolder;

    // Do the interpolation
    DRT::FIBER::UTILS::ProjectFibersToGaussPoints<distype>(ele.Nodes(), shapefcts, fiberHolder);

    params.set("fiberholder", fiberHolder);
  }

  // Call PostSetup of material
  ele.SolidMaterial()->PostSetup(params, ele.Id());
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::InitializeGaussPointDataOutput(
    const DRT::ELEMENTS::Solid& ele) const
{
  const DRT::UTILS::GaussIntegration& integration = *default_integration_;

  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  // Save number of Gauss of the element for gauss point data output
  ele.ParamsInterface().MutableGaussPointDataOutputManagerPtr()->AddElementNumberOfGaussPoints(
      integration.NumPoints());

  // holder for output quantity names and their size
  std::unordered_map<std::string, int> quantities_map{};

  // Ask material for the output quantity names and sizes
  ele.SolidMaterial()->RegisterVtkOutputDataNames(quantities_map);

  // Add quantities to the Gauss point output data manager (if they do not already exist)
  ele.ParamsInterface().MutableGaussPointDataOutputManagerPtr()->MergeQuantities(quantities_map);
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::SolidEleCalc<distype>::EvaluateGaussPointDataOutput(
    const DRT::ELEMENTS::Solid& ele) const
{
  dsassert(ele.IsParamsInterface(),
      "This action type should only be called from the new time integration framework!");

  const DRT::UTILS::GaussIntegration& integration = *default_integration_;
  // Collection and assembly of gauss point data
  for (const auto& quantity :
      ele.ParamsInterface().MutableGaussPointDataOutputManagerPtr()->GetQuantities())
  {
    const std::string& quantity_name = quantity.first;
    const int quantity_size = quantity.second;

    // Step 1: Collect the data for each Gauss point for the material
    LINALG::SerialDenseMatrix gp_data(integration.NumPoints(), quantity_size, true);
    bool data_available = ele.SolidMaterial()->EvaluateVtkOutputData(quantity_name, gp_data);

    // Step 3: Assemble data based on output type (elecenter, postprocessed to nodes, Gauss
    // point)
    if (data_available)
    {
      switch (ele.ParamsInterface().MutableGaussPointDataOutputManagerPtr()->GetOutputType())
      {
        case INPAR::STR::GaussPointDataOutputType::element_center:
        {
          // compute average of the quantities
          Teuchos::RCP<Epetra_MultiVector> global_data =
              ele.ParamsInterface()
                  .MutableGaussPointDataOutputManagerPtr()
                  ->GetMutableElementCenterData()
                  .at(quantity_name);
          DRT::ELEMENTS::AssembleAveragedElementValues(*global_data, gp_data, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::nodes:
        {
          Teuchos::RCP<Epetra_MultiVector> global_data =
              ele.ParamsInterface()
                  .MutableGaussPointDataOutputManagerPtr()
                  ->GetMutableNodalData()
                  .at(quantity_name);

          Epetra_IntVector& global_nodal_element_count =
              *ele.ParamsInterface()
                   .MutableGaussPointDataOutputManagerPtr()
                   ->GetMutableNodalDataCount()
                   .at(quantity_name);

          ExtrapolateGPQuantityToNodesAndAssemble<distype>(
              ele, gp_data, *global_data, false, integration);
          DRT::ELEMENTS::AssembleNodalElementCount(global_nodal_element_count, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::gauss_points:
        {
          std::vector<Teuchos::RCP<Epetra_MultiVector>>& global_data =
              ele.ParamsInterface()
                  .MutableGaussPointDataOutputManagerPtr()
                  ->GetMutableGaussPointData()
                  .at(quantity_name);
          DRT::ELEMENTS::AssembleGaussPointValues(global_data, gp_data, ele);
          break;
        }
        case INPAR::STR::GaussPointDataOutputType::none:
          dserror(
              "You specified a Gauss point data output type of none, so you should not end up "
              "here.");
        default:
          dserror("Unknown Gauss point data output type.");
      }
    }
  }
}

// template classes
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex18>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::SolidEleCalc<DRT::Element::pyramid5>;
