/*----------------------------------------------------------------------*/
/*!
\file drt_meshfree_utils.cpp

\brief service methods for a given meshfree discretisations

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/

#include "drt_meshfree_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_utils.H"

/*--------------------------------------------------------------------------*
 | output errror code as specified in inpar_meshfree.H            nis Feb14 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::OutputMeshfreeError(
  int error)
{
  errorcode err = (errorcode)error;
  switch (err)
  {
  case det_zero:
  {
    std::cout << "Determinant of Jacobian of dual problem too close to zero!" << std::endl;
    break;
  }
  default:
  {
    std::cout << "Error not specified!" << std::endl;
  }
  }

  dserror("Error in meshfree computation. See output above for further information!");
}

/*--------------------------------------------------------------------------*
 | reduces node position to dimension of face                     nis Jan14 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MovePointOfReferenceToFace(
  Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz
  )
{
  Teuchos::RCP<LINALG::SerialDenseMatrix> xyz2 = Teuchos::null;
  return MovePointOfReferenceToFace(xyz,xyz2);
}

void DRT::MESHFREE::MovePointOfReferenceToFace(
  Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz1,
  Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz2
  )
{
  // get auxiliary variables
  const int dim = xyz1->M();
  const int na  = xyz1->N();

  // initialize center
  double center[dim];
  for (int i=0; i<dim; ++i)
    center[i] = 0.0;

  // compute center
  for (int i=0; i<na; ++i)
  {
    const double * a = (*xyz1)[i];
    for (int j=0; j<dim; ++j)
      center[j] += a[j];
  }

  // scale center
  for (int i=0; i<dim; ++i)
    center[i] /= na;

  // translate to center
  for (int i=0; i<na; ++i)
  {
    double * a = (*xyz1)[i];
    for (int j=0; j<dim; ++j)
      a[j] -= center[j];
  }

  if (xyz2!=Teuchos::null)
  {
    const int na2  = xyz2->N();
    // translate to center
    for (int i=0; i<na2; ++i)
    {
      double * a = (*xyz2)[i];
      for (int j=0; j<dim; ++j)
        a[j] -= center[j];
    }
  }

  return;
}

int DRT::MESHFREE::ReduceCloudDimensionOnFaces(
  Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz
  )
{
  Teuchos::RCP<LINALG::SerialDenseMatrix> xyz2 = Teuchos::null;
  return ReduceCloudDimensionOnFaces(xyz,xyz2);
};

int DRT::MESHFREE::ReduceCloudDimensionOnFaces(
  Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz1,
  Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz2
  )
{
  // get auxiliary variables
  const int dim = xyz1->M();
  const int na  = xyz1->N();

  // compute covariance matrix
  LINALG::SerialDenseMatrix C(dim,dim);
  C.Multiply('N','T',1.0,*xyz1,*xyz1,0.0);

  // get eigenvlaues and eigenvectors of covariance matrix
  LINALG::SerialDenseVector E(dim);
  LINALG::SymmetricEigenProblem(C,E,true);

  // count zero eigenvalues (eigenvalues are sorted in ascending order)
  int num_reduced = 0;
  for (int i=0; i<dim; ++i)
    if (E(i)<1e-14)
    {
      if (E(i)>-1e-14)
        ++num_reduced;
      else
        dserror("Covariance matrix has negative eigenvalues.");
    }
    else
      break;

  // determinant of C might be negative for duplicate eigenvalue
  if (dim==3)
  {
    LINALG::Matrix<3,3> C_fsm(C,true);
    if (C_fsm.Determinant()<0)
    {
      if (E(0)<1e-14 and (std::abs(E(1)-E(2))<1e-14))
      {
        for(int i=0; i<3; ++i)
          for(int j=0; j<3; ++j)
            C_fsm(i,j) = C_fsm(i,(j+1)%3);
      }
    }
  }

  // reduce dimension of diffx
  LINALG::SerialDenseMatrix C_reduced(View,C[num_reduced],dim,dim,dim-num_reduced);
  Teuchos::RCP<LINALG::SerialDenseMatrix> xyz_reduced = Teuchos::rcp(new LINALG::SerialDenseMatrix(dim-num_reduced,na));
  xyz_reduced->Multiply('T','N',1.0,C_reduced,*xyz1,0.0);
  xyz1 = xyz_reduced;
  if (xyz2!=Teuchos::null)
  {
    const int na2  = xyz2->N();
    Teuchos::RCP<LINALG::SerialDenseMatrix> xyz2_reduced = Teuchos::rcp(new LINALG::SerialDenseMatrix(dim-num_reduced,na2));
    xyz2_reduced->Multiply('T','N',1.0,C_reduced,*xyz2,0.0);
    xyz2 = xyz2_reduced;
  }

  return num_reduced;
}

// hush! this is a secret place!
