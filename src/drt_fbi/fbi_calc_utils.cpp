/*----------------------------------------------------------------------*/
/*! \file
 *
 *\brief Utility functions for fluid beam interaction related calculations
 *
 *\level 3
 *
 */
/*----------------------------------------------------------------------*/
#include "fbi_calc_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_element.H"
#include "../drt_beam3/beam3_base.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/

void FBI::UTILS::GetFBIElementCenterlineDOFIndices(DRT::Discretization const& discret,
    const DRT::Element* ele, std::vector<unsigned int>& ele_centerline_dof_indices,
    unsigned int& num_dof)
{
  // Todo implement method in DRT::Element or find alternative way of doing this
  // find out the elements' number of Dofs (=dimension of element vector/matrices)
  std::vector<int> lmrow;
  std::vector<int> dummy1, dummy2;

  ele->LocationVector(discret, lmrow, dummy1, dummy2);
  num_dof = lmrow.size();

  const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

  if (beamele != NULL)
  {
    beamele->CenterlineDofIndicesOfElement(ele_centerline_dof_indices);
  }
  else
  {
    ele_centerline_dof_indices.resize(num_dof * 3 / 4);
    int j = 0;
    for (unsigned int i = 0; i < num_dof; ++i)
    {
      if (!((i + 1) % 4)) ++i;
      ele_centerline_dof_indices[j] = i;
      j++;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void FBI::UTILS::AssembleCenterlineDofForceStiffIntoFBIElementForceStiff(
    const DRT::Discretization& discretization1, const DRT::Discretization& discretization2,
    std::vector<int> const& elegid,
    std::vector<LINALG::SerialDenseVector> const& eleforce_centerlineDOFs,
    std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elestiff_centerlineDOFs,
    std::vector<LINALG::SerialDenseVector>* eleforce,
    std::vector<std::vector<LINALG::SerialDenseMatrix>>* elestiff)
{
  std::vector<unsigned int> numdof_ele(2);
  std::vector<std::vector<unsigned int>> ele_centerlinedofindices(2);

  // Get DOFs for beam element
  DRT::Element* ele = discretization1.gElement(elegid[0]);
  GetFBIElementCenterlineDOFIndices(
      discretization1, ele, ele_centerlinedofindices[0], numdof_ele[0]);

  // Get DOFs for fluid element
  ele = discretization2.gElement(elegid[1]);
  GetFBIElementCenterlineDOFIndices(
      discretization2, ele, ele_centerlinedofindices[1], numdof_ele[1]);


  // assemble centerline DOF values correctly into element DOFvec vectors/matrices
  if (eleforce != NULL)
  {
    for (unsigned int iele = 0; iele < 2; ++iele)
    {
      // resize and clear variable
      ((*eleforce)[iele]).Size(numdof_ele[iele]);

      // safety check: dimensions
      if ((unsigned int)eleforce_centerlineDOFs[iele].RowDim() !=
          ele_centerlinedofindices[iele].size())
        dserror(
            "size mismatch! need to assemble %d values of centerline-Dof based "
            "force vector into element vector but only got %d element-local Dof indices",
            eleforce_centerlineDOFs[iele].RowDim(), ele_centerlinedofindices[iele].size());

      // Todo maybe use a more general 'SerialDenseAssemble' method here
      for (unsigned int idof = 0; idof < ele_centerlinedofindices[iele].size(); ++idof)
        ((*eleforce)[iele])(ele_centerlinedofindices[iele][idof]) =
            eleforce_centerlineDOFs[iele](idof);
    }
  }

  if (elestiff != NULL)
  {
    for (unsigned int iele = 0; iele < 2; ++iele)
    {
      for (unsigned int jele = 0; jele < 2; ++jele)
      {
        // resize and clear variable
        ((*elestiff)[iele][jele]).Shape(numdof_ele[iele], numdof_ele[jele]);

        // safety check: dimensions
        if ((unsigned int)elestiff_centerlineDOFs[iele][jele].RowDim() !=
            ele_centerlinedofindices[iele].size())
          dserror(
              "size mismatch! need to assemble %d row values of centerline-Dof based "
              "stiffness matrix into element matrix but only got %d element-local Dof indices",
              elestiff_centerlineDOFs[iele][jele].RowDim(), ele_centerlinedofindices[iele].size());

        if ((unsigned int)elestiff_centerlineDOFs[iele][jele].ColDim() !=
            ele_centerlinedofindices[jele].size())
          dserror(
              "size mismatch! need to assemble %d column values of centerline-Dof based "
              "stiffness matrix into element matrix but only got %d element-local Dof indices",
              elestiff_centerlineDOFs[iele][jele].ColDim(), ele_centerlinedofindices[jele].size());

        for (unsigned int idof = 0; idof < ele_centerlinedofindices[iele].size(); ++idof)
          for (unsigned int jdof = 0; jdof < ele_centerlinedofindices[jele].size(); ++jdof)
            ((*elestiff)[iele][jele])(
                ele_centerlinedofindices[iele][idof], ele_centerlinedofindices[jele][jdof]) =
                elestiff_centerlineDOFs[iele][jele](idof, jdof);
      }
    }
  }
}
