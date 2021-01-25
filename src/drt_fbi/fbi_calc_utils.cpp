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
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_sparseoperator.H"

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

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void FBI::UTILS::FEAssembleEleForceStiffIntoSystemVectorMatrices(
    const DRT::Discretization& discretization1, const DRT::Discretization& discretization2,
    std::vector<int> const& elegid, std::vector<LINALG::SerialDenseVector> const& elevec,
    std::vector<std::vector<LINALG::SerialDenseMatrix>> const& elemat,
    Teuchos::RCP<Epetra_FEVector>& f1, Teuchos::RCP<Epetra_FEVector>& f2,
    Teuchos::RCP<LINALG::SparseMatrix>& c11, Teuchos::RCP<LINALG::SparseOperator> c22,
    Teuchos::RCP<LINALG::SparseMatrix>& c12, Teuchos::RCP<LINALG::SparseMatrix>& c21)
{
  // the entries of elevecX  belong to the Dofs of the element with GID elegidX
  // the rows    of elematXY belong to the Dofs of the element with GID elegidX
  // the columns of elematXY belong to the Dofs of the element with GID elegidY
  const DRT::Element* ele1 = discretization1.gElement(elegid[0]);
  const DRT::Element* ele2 = discretization2.gElement(elegid[1]);

  // get element location vector and ownerships
  std::vector<int> lmrow1;
  std::vector<int> lmrow2;
  std::vector<int> lmrowowner1;
  std::vector<int> lmrowowner2;
  std::vector<int> lmstride;

  ele1->LocationVector(discretization1, lmrow1, lmrowowner1, lmstride);
  ele2->LocationVector(discretization2, lmrow2, lmrowowner2, lmstride);

  // assemble both element vectors into global system vector
  if (f1 != Teuchos::null)
  {
    f1->SumIntoGlobalValues(elevec[0].Length(), &lmrow1[0], elevec[0].Values());
  }
  if (f2 != Teuchos::null)
  {
    f2->SumIntoGlobalValues(elevec[1].Length(), &lmrow2[0], elevec[1].Values());
  }

  // and finally also assemble stiffness contributions
  if (c11 != Teuchos::null)
  {
    c11->FEAssemble(elemat[0][0], lmrow1, lmrow1);
  }
  if (c12 != Teuchos::null)
  {
    c12->FEAssemble(elemat[0][1], lmrow1, lmrow2);
  }
  if (c21 != Teuchos::null)
  {
    c21->FEAssemble(elemat[1][0], lmrow2, lmrow1);
  }
  if (c22 != Teuchos::null)
  {
    /* TODO Find out how to deal with this in parallel! For now the coupling contributions are
     * computed on the beam elements. Maybe we have to tackle that before being able to compute
     * stuff in parallel.
     *
     */
    c22->Assemble(elegid[1], lmstride, elemat[1][1], lmrow2, lmrowowner2, lmrow2);
  }
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
