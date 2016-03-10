/*!-----------------------------------------------------------------------------------------------------------
\file beam3tobeampotential.cpp
\brief One beam potential pair (two beam elements)

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*-----------------------------------------------------------------------------------------------------------*/

#include "beam3tobeampotential.H"
#include "beam3contact_defines.H"
#include "beam3contact_utils.H"
#include "../drt_inpar/inpar_beampotential.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_structure/strtimint_impl.H"
#include "../drt_beam3/beam3.H"
#include "../drt_beam3r/beam3r.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../headers/FAD_utils.H"

#include "Teuchos_TimeMonitor.hpp"

/*----------------------------------------------------------------------*
 |  constructor (public)                                     grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::Beam3tobeampotential(
                                                                     const DRT::Discretization& pdiscret,
                                                                     const DRT::Discretization& cdiscret,
                                                                     const std::map<int,int>& dofoffsetmap,
                                                                     DRT::Element* element1,
                                                                     DRT::Element* element2,
                                                                     Teuchos::ParameterList beampotparams,
                                                                     std::vector<DRT::Condition*> linechargeconds):
pdiscret_(pdiscret),
cdiscret_(cdiscret),
dofoffsetmap_(dofoffsetmap),
element1_(element1),
element2_(element2),
bpotparams_(beampotparams),
k_(0.0),
m_(0.0),
iter_(0),
numstep_(0),
dt_(0.0),
ele1length_(0.0),
ele2length_(0.0),
radius1_(0.0),
radius2_(0.0)
{
  ele1pos_.Clear();
  ele2pos_.Clear();

  fpot1_.Clear();
  fpot2_.Clear();
  stiffpot1_.Clear();
  stiffpot2_.Clear();

  const DRT::ElementType & eot1 = element1_->ElementType();

  //Calculate initial length of beam elements (approximation for initially curved elements!)
  LINALG::TMatrix<double,3,1> lvec1(true);
  LINALG::TMatrix<double,3,1> lvec2(true);
  for(int i=0;i<3;i++)
  {
    lvec1(i)=(element1_->Nodes())[0]->X()[i]-(element1_->Nodes())[1]->X()[i];
    lvec2(i)=(element2_->Nodes())[0]->X()[i]-(element2_->Nodes())[1]->X()[i];
  }
  ele1length_=lvec1.Norm2();
  ele2length_=lvec2.Norm2();

  if (element1->ElementType() != element2->ElementType())
    dserror("The class beam3tobeampotential only works for pairs of the same beam element type!");

  if (element1->Id() >= element2->Id())
    dserror("Element 1 has to have the smaller element-ID. Adapt your contact search!");

  double Iyy1 = 0.0;
  double Iyy2 = 0.0;
  if (eot1 == DRT::ELEMENTS::Beam3Type::Instance())
  {
    Iyy1 = (static_cast<DRT::ELEMENTS::Beam3*>(element1_))->Iyy();
    Iyy2 = (static_cast<DRT::ELEMENTS::Beam3*>(element2_))->Iyy();
  }
  else if (eot1 == DRT::ELEMENTS::Beam3rType::Instance())
  {
    Iyy1 = (static_cast<DRT::ELEMENTS::Beam3r*>(element1_))->Iyy();
    Iyy2 = (static_cast<DRT::ELEMENTS::Beam3r*>(element2_))->Iyy();
  }
  else if (eot1 == DRT::ELEMENTS::Beam3ebType::Instance())
  {
    Iyy1 = (static_cast<DRT::ELEMENTS::Beam3eb*>(element1_))->Iyy();
    Iyy2 = (static_cast<DRT::ELEMENTS::Beam3eb*>(element2_))->Iyy();
  }
  else
    dserror("Only Beam3, Beam3r and Beam3eb elements allowed in Beam3tobeampotential!");

  radius1_ = MANIPULATERADIUS * sqrt(sqrt(4 * Iyy1 / M_PI));
  radius2_ = MANIPULATERADIUS * sqrt(sqrt(4 * Iyy2 / M_PI));

  // initialize line charge conditions applied to element1 and element2
  linechargeconds_.clear();
  if (linechargeconds.size() == 2)
  {
    for (int i=0; i<2; ++i)
    {
      if (linechargeconds[i]->Type() == DRT::Condition::BeamPotential_LineChargeDensity)
        linechargeconds_.push_back(linechargeconds[i]);
      else
        dserror("Provided line charge condition is not of correct type BeamPotential_LineChargeDensity!");
    }
  }
  else
    dserror("Expected TWO dline charge conditions!");

  return;
}

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                                grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::Beam3tobeampotential(const Beam3tobeampotential& old):
pdiscret_(old.pdiscret_),
cdiscret_(old.cdiscret_),
dofoffsetmap_(old.dofoffsetmap_)
{
  dserror("ERROR: Copy constructor incomplete");
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate the element (public)                             grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::Evaluate( LINALG::SparseMatrix& stiffmatrix,
                                                                   Epetra_Vector& fint,
                                                                   const double k,
                                                                   const double m)
{
  // reset fpot and stiffpot class variables
  fpot1_.Clear();
  fpot2_.Clear();
  stiffpot1_.Clear();
  stiffpot2_.Clear();

  // set class variables
  k_=k;
  m_=m;

  // prepare FAD
#ifdef AUTOMATICDIFF
  // The 2*3*numnodes*numnodalvalues primary DoFs are the components of the nodal positions / tangents.
  for (int i=0;i<3*numnodes*numnodalvalues;i++)
    ele1pos_(i).diff(i,2*3*numnodes*numnodalvalues);

  for (int i=0;i<3*numnodes*numnodalvalues;i++)
    ele2pos_(i).diff(3*numnodes*numnodalvalues+i,2*3*numnodes*numnodalvalues);
#endif

  // compute the values for fpot1, fpot2, stiffpot1, stiffpot2
  // TODO specify approximation type in input file: FullInt, LargeSepApprox, SmallSepApprox (Langbein)
  EvaluateFpotandStiffpot_LargeSepApprox();

  // assemble fpot1 and fpot2 into fint
  AssembleFpot(fint);

  // assemble stiffpot1 and 2 into stiffmatrix
  AssembleStiffpot(stiffmatrix);

  return (true);
}

/*----------------------------------------------------------------------*
 |  Compute contact forces                                   grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::EvaluateFpotandStiffpot_LargeSepApprox()
{
  // Set gauss integration rule
  DRT::UTILS::GaussRule1D gaussrule = DRT::UTILS::intrule_line_5point;

  // Get gauss points (gp) for integration
  DRT::UTILS::IntegrationPoints1D gausspoints(gaussrule);
  // number of gps
  const int numgp=gausspoints.nquad;

  // vectors for shape functions and their derivatives
  // Attention: these are individual shape function values, NOT shape function matrices
  // values at all gauss points are stored in advance (more efficient due to double integral)
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N1_i(numgp);        // = N1_i
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N2_i(numgp);        // = N2_i
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N1_i_xi(numgp);     // = N1_i,xi
  std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> > N2_i_xi(numgp);     // = N2_i,eta

  // coords and derivatives of the two gauss points
  LINALG::TMatrix<TYPE, 3, 1> r1(true);                               // = r1
  LINALG::TMatrix<TYPE, 3, 1> r2(true);                               // = r2
  LINALG::TMatrix<TYPE, 3, 1> dist(true);                             // = r1-r2
  TYPE norm_dist= 0.0;                                                // = |r1-r2|

  // Evaluate shape functions at gauss points and store values
  GetShapeFunctions(N1_i,N2_i,N1_i_xi,N2_i_xi,gausspoints);

  // evaluate charge densities from DLINE charge condition specified in input file
  double q1 = linechargeconds_[0]->GetDouble("val");
  double q2 = linechargeconds_[1]->GetDouble("val");

  // TODO evaluate given functions in line charge conditions! for now: dserror
  if (linechargeconds_[0]->GetInt("funct") != 0 or linechargeconds_[1]->GetInt("funct") != 0)
    dserror("DLINE beam potential charge condition: No functions allowed yet!");

  // compute Jacobi factors
  std::vector<double> jacobifac1(numgp,1.0);
  std::vector<double> jacobifac2(numgp,1.0);

  ComputeJacobiFacs(jacobifac1, jacobifac2,N1_i_xi,N2_i_xi,gausspoints);

  // auxiliary variable
  LINALG::TMatrix<TYPE, 3, 1> fpot_tmp(true);

  // determine prefactor of the integral (depends on whether surface or volume potential is applied)
  double prefactor = k_ * m_;

  switch(DRT::INPUT::IntegralValue<INPAR::BEAMPOTENTIAL::BeamPotentialType>(bpotparams_,"BEAMPOTENTIAL_TYPE"))
  {
  case INPAR::BEAMPOTENTIAL::beampot_surf:
    prefactor *= 4 * radius1_ * radius2_ * pow(M_PI,2);
    break;
  case INPAR::BEAMPOTENTIAL::beampot_vol:
    prefactor *= pow(radius1_,2) * pow(radius2_,2) * pow(M_PI,2);
    break;
  default:
    dserror("No valid BEAMPOTENTIAL_TYPE specified. Choose either Surface or Volume in input file!");
  }

  // loop over gauss points on ele1
  for (int gp1=0; gp1 < numgp; ++gp1)
  {
    ComputeCoords(r1,N1_i[gp1],ele1pos_);

    // loop over gauss points on ele2
    for (int gp2=0; gp2 < numgp; ++gp2)
    {

      ComputeCoords(r2,N2_i[gp2],ele2pos_);

      dist = FADUTILS::DiffVector(r1,r2);

      norm_dist = FADUTILS::VectorNorm<3>(dist);

      // auxiliary variables to store pre-calculated common terms
      TYPE norm_dist_exp1 = 0.0;
      if(norm_dist !=0.0)
      {
        norm_dist_exp1 = pow(norm_dist,-m_-2);
      }
      else
      {
        dserror("\n|r1-r2|=0 ! Interacting points are identical! Potential law not defined in this case! Think about shifting nodes in unconverged state?!");
      }
      double q1q2_JacFac_GaussWeights = q1 * q2 * jacobifac1[gp1] * jacobifac2[gp2] * gausspoints.qwgt[gp1] * gausspoints.qwgt[gp2];

      // compute fpot_tmp here, same for both element forces
      for (int i=0; i<3; ++i) fpot_tmp(i) = q1q2_JacFac_GaussWeights * norm_dist_exp1 * dist(i);

      //********************************************************************
      // calculate fpot1: force on element 1
      //********************************************************************
      // sum up the contributions of all nodes (in all dimensions)
      for (int i=0; i<(numnodes*numnodalvalues); ++i)
      {
        // loop over dimensions
        for (int j=0; j<3; ++j)
        {
          fpot1_(3*i+j) -= N1_i[gp1](i)*fpot_tmp(j);
        }
      }

      //********************************************************************
      // calculate fpot2: force on element 2
      //********************************************************************
      // sum up the contributions of all nodes (in all dimensions)
      for (int i=0; i<(numnodes*numnodalvalues); ++i)
      {
        // loop over dimensions
        for (int j=0; j<3; ++j)
        {
          fpot2_(3*i+j) += N2_i[gp2](i)*fpot_tmp(j);
        }
      }

      //********************************************************************
      // calculate stiffpot1
      //********************************************************************
      // auxiliary variables (same for both elements)
      TYPE norm_dist_exp2 = (m_+2) * pow(norm_dist,-m_-4);

      LINALG::TMatrix<TYPE, 3, 3> dist_dist_T(true);

      for (int i=0; i<3; ++i)
      {
        for (int j=0; j<=i; ++j)
        {
          dist_dist_T(i,j) = dist(i) * dist(j);
          if(i!=j) dist_dist_T(j,i) = dist_dist_T(i,j);
        }
      }

      for (int i=0; i<(numnodes*numnodalvalues); ++i)
      {

        // d (Res_1) / d (d_1)
        for (int j=0; j<(numnodes*numnodalvalues); ++j)
        {

          for (int idim=0; idim<3; ++idim)
          {
            stiffpot1_(3*i+idim,3*j+idim) -= norm_dist_exp1 * N1_i[gp1](i)*N1_i[gp1](j) *  q1q2_JacFac_GaussWeights;

            for (int jdim=0; jdim<3; ++jdim)
            {
              stiffpot1_(3*i+idim,3*j+jdim) += norm_dist_exp2 * N1_i[gp1](i) * dist_dist_T(idim,jdim) * N1_i[gp1](j) *  q1q2_JacFac_GaussWeights;
            }
          }
        }

        // d (Res_1) / d (d_2)
        for (int j=0; j<(numnodes*numnodalvalues); ++j)
        {

          for (int idim=0; idim<3; ++idim)
          {
            stiffpot1_(3*i+idim,3*(numnodes*numnodalvalues+j)+idim) += norm_dist_exp1 * N1_i[gp1](i)*N2_i[gp2](j) *  q1q2_JacFac_GaussWeights;

            for (int jdim=0; jdim<3; ++jdim)
            {
              stiffpot1_(3*i+idim,3*(numnodes*numnodalvalues+j)+jdim) -= norm_dist_exp2 * N1_i[gp1](i) * dist_dist_T(idim,jdim) * N2_i[gp2](j) *  q1q2_JacFac_GaussWeights;
            }
          }
        }

      }

      //********************************************************************
      // calculate stiffpot2
      //********************************************************************
      for (int i=0; i<(numnodes*numnodalvalues); ++i)
      {

        // d (Res_2) / d (d_1)
        for (int j=0; j<(numnodes*numnodalvalues); ++j)
        {

          for (int idim=0; idim<3; ++idim)
          {
            stiffpot2_(3*i+idim,3*j+idim) += norm_dist_exp1 * N2_i[gp2](i)*N1_i[gp1](j) *  q1q2_JacFac_GaussWeights;

            for (int jdim=0; jdim<3; ++jdim)
            {
              stiffpot2_(3*i+idim,3*j+jdim) -= norm_dist_exp2 * N2_i[gp2](i) * dist_dist_T(idim,jdim) * N1_i[gp1](j) *  q1q2_JacFac_GaussWeights;
            }
          }
        }

        // d (Res_2) / d (d_2)
        for (int j=0; j<(numnodes*numnodalvalues); ++j)
        {

          for (int idim=0; idim<3; ++idim)
          {
            stiffpot2_(3*i+idim,3*(numnodes*numnodalvalues+j)+idim) -= norm_dist_exp1 * N2_i[gp2](i)*N2_i[gp2](j) *  q1q2_JacFac_GaussWeights;

            for (int jdim=0; jdim<3; ++jdim)
            {
              stiffpot2_(3*i+idim,3*(numnodes*numnodalvalues+j)+jdim) += norm_dist_exp2 * N2_i[gp2](i) * dist_dist_T(idim,jdim) * N2_i[gp2](j) *  q1q2_JacFac_GaussWeights;
            }
          }
        }

      }

    } // end gauss quadrature loop (element 2)
  } // end gauss quadrature loop (element 1)


  // apply constant prefactor
  for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i) fpot1_(i)*=prefactor;
  for (unsigned int i=0; i<3*numnodes*numnodalvalues; ++i) fpot2_(i)*=prefactor;

  stiffpot1_.Scale(prefactor);
  stiffpot2_.Scale(prefactor);

  return;
}

/*----------------------------------------------------------------------*
 |  assemble fpot1/2                                         grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tobeampotential<numnodes,numnodalvalues>::AssembleFpot(Epetra_Vector& fint)
{
  // get dimensions for vectors fpot1 and fpot2
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  // temporary vectors for DOF-GIDs and owning procs
  std::vector<int> lm1(dim1);
  std::vector<int> lm2(dim2);
  std::vector<int> lmowner1(dim1);
  std::vector<int> lmowner2(dim2);

  // node ids of both elements
  const int* node_ids1 = element1_->NodeIds();
  const int* node_ids2 = element2_->NodeIds();

  // prepare assembly for element1 (fill lm1 and lmowner1)
  for (int i=0;i<numnodes;++i)
  {
    // get node pointer and dof ids
    DRT::Node* node = cdiscret_.gNode(node_ids1[i]);
    std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

    // prepare assembly
    for (int j=0;j<3*numnodalvalues;++j)
    {
      lm1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
      lmowner1[3*numnodalvalues*i+j] = node->Owner();
    }
  }

  // prepare assembly for element2
  for(int i=0;i<numnodes;++i)
  {
    // get node pointer and dof ids
    DRT::Node* node = cdiscret_.gNode(node_ids2[i]);
    std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

    // compute force vector Fc2 and prepare assembly
    for (int j=0;j<3*numnodalvalues;++j)
    {
      lm2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
      lmowner2[3*numnodalvalues*i+j] = node->Owner();
    }
  }

  // create copy (Epetra_Vector) of fpot1/2 for LINALG::Assemble() fcn
  Epetra_SerialDenseVector fpot1_copy(dim1);
  Epetra_SerialDenseVector fpot2_copy(dim2);
  for (int i=0; i<dim1; i++) fpot1_copy[i]= - FADUTILS::CastToDouble(fpot1_(i));           // ATTENTION: negative sign because forces are treated as EXTERNAL forces in beam3contact_manager
  for (int i=0; i<dim2; i++) fpot2_copy[i]= - FADUTILS::CastToDouble(fpot2_(i));

  // assemble into global force vector
  LINALG::Assemble(fint,fpot1_copy,lm1,lmowner1);
  LINALG::Assemble(fint,fpot2_copy,lm2,lmowner2);

  return;
}

/*----------------------------------------------------------------------*
 |  assemble stiffpot1/2                                     grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tobeampotential<numnodes,numnodalvalues>::AssembleStiffpot(LINALG::SparseMatrix& stiffmatrix)
{
  // get dimensions for vectors fc1 and fc2
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3*numnodes*numnodalvalues;

  std::vector<int> lmrow1(dim1);
  std::vector<int> lmrow2(dim2);
  std::vector<int> lmrowowner1(dim1);
  std::vector<int> lmrowowner2(dim2);
  std::vector<int> lmcol1(dim1+dim2);
  std::vector<int> lmcol2(dim1+dim2);

  // prepare assembly

  // node ids of both elements
  const int* node_ids1 = element1_->NodeIds();
  const int* node_ids2 = element2_->NodeIds();

  // fill lmrow1 and lmrowowner1
  for (int i=0;i<numnodes;++i)
  {
    // get pointer and dof ids
    DRT::Node* node = cdiscret_.gNode(node_ids1[i]);
    std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

    for (int j=0;j<3*numnodalvalues;++j)
    {
      lmrow1[3*numnodalvalues*i+j]=NodeDofGIDs[j];
      lmrowowner1[3*numnodalvalues*i+j]=node->Owner();
    }
  }

  // fill lmrow2 and lmrowowner2
  for (int i=0;i<numnodes;++i)
  {
    // get pointer and node ids
    DRT::Node* node = cdiscret_.gNode(node_ids2[i]);
    std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

    for (int j=0;j<3*numnodalvalues;++j)
    {
      lmrow2[3*numnodalvalues*i+j]=NodeDofGIDs[j];
      lmrowowner2[3*numnodalvalues*i+j]=node->Owner();
    }
  }

  // fill lmcol1 and lmcol2
  for (int i=0;i<numnodes;++i)
  {
    // get pointer and node ids
    DRT::Node* node = cdiscret_.gNode(node_ids1[i]);
    std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

    for (int j=0;j<3*numnodalvalues;++j)
    {
      lmcol1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
      lmcol2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
    }
  }

  // fill lmcol1 and lmcol2
  for (int i=0;i<numnodes;++i)
  {
    // get pointer and node ids
    DRT::Node* node = cdiscret_.gNode(node_ids2[i]);
    std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

    for (int j=0;j<3*numnodalvalues;++j)
    {
      lmcol1[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
      lmcol2[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
    }
  }


  Epetra_SerialDenseMatrix stiffpot1_copy(dim1,dim1+dim2);
  Epetra_SerialDenseMatrix stiffpot2_copy(dim2,dim1+dim2);

#ifndef AUTOMATICDIFF
  for (int j=0;j<dim1+dim2;j++)
  {
    for (int i=0;i<dim1;i++)
      stiffpot1_copy(i,j) = FADUTILS::CastToDouble(stiffpot1_(i,j));
    for (int i=0;i<dim2;i++)
      stiffpot2_copy(i,j) = FADUTILS::CastToDouble(stiffpot2_(i,j));
  }
#else    // automatic differentiation via FAD for debugging
  for (int j=0;j<dim1+dim2;j++)
  {
    for (int i=0;i<dim1;i++)
      stiffpot1_copy(i,j) = FADUTILS::CastToDouble(fpot1_(i).dx(j));
    for (int i=0;i<dim2;i++)
      stiffpot2_copy(i,j) = FADUTILS::CastToDouble(fpot2_(i).dx(j));
  }

#ifdef FADCHECKS
  std::cout << "BTB Potential Pair: " << element1_->Id() << " / " << element2_->Id() << std::endl;

  std::cout << "stiffpot1: " << std::endl;
  for (int i=0;i<dim1;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      std::cout << stiffpot1_(i,j).val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "stiffpot1_FAD: " << std::endl;
  for (int i=0;i<dim1;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      std::cout << stiffpot1_copy(i,j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "stiffpot2: " << std::endl;
  for (int i=0;i<dim2;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      std::cout << stiffpot2_(i,j).val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout << "stiffpot2_FAD: " << std::endl;
  for (int i=0;i<dim2;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      std::cout << stiffpot2_copy(i,j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
#endif // FADCHECKS

#endif //AUTOMATICDIFF

  // assemble into global tangential stiffness matrix
  stiffmatrix.Assemble(0,stiffpot1_copy,lmrow1,lmrowowner1,lmcol1);
  stiffmatrix.Assemble(0,stiffpot2_copy,lmrow2,lmrowowner2,lmcol2);

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                 meier 01/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::GetShapeFunctions( std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i,
                                                                                 std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N2_i,
                                                                                 std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i_xi,
                                                                                 std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N2_i_xi,
                                                                                 DRT::UTILS::IntegrationPoints1D& gausspoints)
{
  // get both discretization types
  const DRT::Element::DiscretizationType distype1 = element1_->Shape();
  const DRT::Element::DiscretizationType distype2 = element2_->Shape();

  if (numnodalvalues==1)
  {
    for (int gp=0; gp<gausspoints.nquad; ++gp)
    {
      // get values and derivatives of shape functions
      DRT::UTILS::shape_function_1D(N1_i[gp], gausspoints.qxg[gp][0], distype1);
      DRT::UTILS::shape_function_1D(N2_i[gp], gausspoints.qxg[gp][0], distype2);
      DRT::UTILS::shape_function_1D_deriv1(N1_i_xi[gp], gausspoints.qxg[gp][0], distype1);
      DRT::UTILS::shape_function_1D_deriv1(N2_i_xi[gp], gausspoints.qxg[gp][0], distype2);
    }
  }
  else if (numnodalvalues==2)
  {

    if ( element1_->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance() )
      dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");

    if ( element2_->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance() )
      dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");

    double length1 = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(element1_))->jacobi();
    double length2 = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(element2_))->jacobi();

    for (int gp=0; gp<gausspoints.nquad; ++gp)
    {
      // get values and derivatives of shape functions
      DRT::UTILS::shape_function_hermite_1D(N1_i[gp], gausspoints.qxg[gp][0], length1, distype1);
      DRT::UTILS::shape_function_hermite_1D(N2_i[gp], gausspoints.qxg[gp][0], length2, distype2);
      DRT::UTILS::shape_function_hermite_1D_deriv1(N1_i_xi[gp], gausspoints.qxg[gp][0], length1, distype1);
      DRT::UTILS::shape_function_hermite_1D_deriv1(N2_i_xi[gp], gausspoints.qxg[gp][0], length2, distype2);
    }
  }
  else
    dserror("Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) values are valid!");

  return;
}

/*----------------------------------------------------------------------*
 |  Compute Jacobi factors at Gauss points                    grill 11/13|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::ComputeJacobiFacs(std::vector<double>& jacobifac1,
                                                          std::vector<double>& jacobifac2,
                                                          std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N1_i_xi,
                                                          std::vector<LINALG::Matrix<1, numnodes*numnodalvalues> >& N2_i_xi,
                                                          DRT::UTILS::IntegrationPoints1D& gausspoints)
{
  // node ids of both elements
  const int* node_ids1 = element1_->NodeIds();
  const int* node_ids2 = element2_->NodeIds();

  //derivative of curve in physical space with respect to curve parameter xi \in [-1;1] on element level
  LINALG::Matrix<3,1> drdxi;

  // loop over gps
  for (int ngp=0; ngp<gausspoints.nquad; ++ngp)
  {

    drdxi.Clear();

    // ele1
    // calculate vector drdxi
    for(int nnode=0; nnode<numnodes; nnode++)
    {
      // get pointer and dof ids
      DRT::Node* node1 = cdiscret_.gNode(node_ids1[nnode]);

      for(int dof=0; dof<3 ; dof++) drdxi(dof) += N1_i_xi[ngp](nnode) * node1->X()[dof];

      //Store Jacobi determinant with respect to reference configuration
      jacobifac1[ngp]= drdxi.Norm2();
    }

    drdxi.Clear();

    // ele2
    // calculate vector drdxi
    for(int nnode=0; nnode<numnodes; nnode++)
    {
      // get pointer and dof ids
      DRT::Node* node2 = cdiscret_.gNode(node_ids2[nnode]);

      for(int dof=0; dof<3 ; dof++) drdxi(dof) += N2_i_xi[ngp](nnode) * node2->X()[dof];

      //Store Jacobi determinant with respect to reference configuration
      jacobifac2[ngp]= drdxi.Norm2();
    }

  } //end: loop over gps

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute Jacobi factor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | compute gauss point coordinates                           grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::ComputeCoords(LINALG::TMatrix<TYPE,3,1>& r,
                                                                            const LINALG::Matrix<1,numnodes*numnodalvalues>& N_i,
                                                                            const LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues,1> elepos)
{
  r.Clear();

  // compute output variable
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<numnodes*numnodalvalues;j++)
    {
      r(i)+=N_i(j)*elepos(3*j+i);
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Get global dofs of a node                                 meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
std::vector<int> CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::GetGlobalDofs(const DRT::Node* node)
{
  // get dofs in beam contact discretization
  std::vector<int> cdofs = ContactDiscret().Dof(node);

  // get dofs in problem discretization via offset
  std::vector<int> pdofs((int)(cdofs.size()));
  for (int k=0;k<(int)(cdofs.size());++k)
  {
    pdofs[k]=(dofoffsetmap_.find(cdofs[k]))->second;
  }

  return pdofs;
}

/*----------------------------------------------------------------------*
 |  Update nodal coordinates (public)                        meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tobeampotential<numnodes, numnodalvalues>::UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                                                                      Epetra_SerialDenseMatrix& newele2pos)
{
  for (int i=0;i<3*numnodalvalues;i++)
  {
    for (int j=0;j<numnodes;j++)
    {
      ele1pos_(3*numnodalvalues*j+i)=newele1pos(i,j);
      ele2pos_(3*numnodalvalues*j+i)=newele2pos(i,j);
    }
  }

  return;
}


Teuchos::RCP<CONTACT::Beam3tobeampotentialinterface> CONTACT::Beam3tobeampotentialinterface::Impl( const int numnodes,
                                                                      const int numnodalvalues,
                                                                      const DRT::Discretization& pdiscret,
                                                                      const DRT::Discretization& cdiscret,
                                                                      const std::map<int,int>& dofoffsetmap,
                                                                      DRT::Element* element1,
                                                                      DRT::Element* element2,
                                                                      Teuchos::ParameterList beampotparams,
                                                                      std::vector<DRT::Condition*> linechargeconds)
{

  switch (numnodalvalues)
  {
    case 1:
    {
      switch (numnodes)
      {
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tobeampotential<2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beampotparams,linechargeconds));
        }
        case 3:
        {
          return Teuchos::rcp (new CONTACT::Beam3tobeampotential<3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beampotparams,linechargeconds));
        }
        case 4:
        {
          return Teuchos::rcp (new CONTACT::Beam3tobeampotential<4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beampotparams,linechargeconds));
        }
        case 5:
        {
          return Teuchos::rcp (new CONTACT::Beam3tobeampotential<5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beampotparams,linechargeconds));
        }
        default:
          dserror("No valid template parameter for the number of nodes (numnodes = 2,3,4,5 for Reissner beams) available!");
          break;
      }
      break;
    }
    case 2:
    {
      switch (numnodes)
      {
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tobeampotential<2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beampotparams,linechargeconds));
        }
        default:
          dserror("No valid template parameter for the number of nodes (only numnodes = 2 for Kirchhoff beams valid so far) available!");
          break;
      }
      break;
    }
    default:
      dserror("No valid template parameter for the Number of nodal values (numnodalvalues = 1 for Reissner beams, numnodalvalues = 2 for Kirchhoff beams) available!");
      break;
  }
  return Teuchos::null;

}


//Possible template cases: this is necessary for the compiler
template class CONTACT::Beam3tobeampotential<2,1>;
template class CONTACT::Beam3tobeampotential<3,1>;
template class CONTACT::Beam3tobeampotential<4,1>;
template class CONTACT::Beam3tobeampotential<5,1>;
template class CONTACT::Beam3tobeampotential<2,2>;
