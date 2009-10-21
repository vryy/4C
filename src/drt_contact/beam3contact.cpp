/*!-----------------------------------------------------------------------------------------------------------
\file beam3contact.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*-----------------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

//compile only if beam3 element is complied, too, as beam3 element required for member variables of this class
#ifdef D_BEAM3

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "beam3contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                      cyron 10/09|
 *----------------------------------------------------------------------*/
Beam3contact::Beam3contact(const DRT::ELEMENTS::Beam3& element1, const DRT::ELEMENTS::Beam3& element2):
element1_(element1),
element2_(element2)
{
  return;
}
/*----------------------------------------------------------------------*
 |  copy-constructor (public)                                 cyron 10/09|
 *----------------------------------------------------------------------*/
Beam3contact::Beam3contact(const Beam3contact& old):
element1_(old.element1_),
element2_(old.element2_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  destructor (public)                                     cyron 10/09 |
 *----------------------------------------------------------------------*/
Beam3contact::~Beam3contact()
{
  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                                 cyron 01/08|
 *----------------------------------------------------------------------------------------------------------*/
int Beam3contact::Evaluate(ParameterList&            params,
                           DRT::Discretization&      discretization,
                           vector<int>&              lm,
                           Epetra_SerialDenseMatrix& stiffmatrix,
                           Epetra_SerialDenseVector& fint)
{

      //extract displacement components of all nodes involved and arrange them in vector mydisp
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) dserror("Cannot get state vectors 'displacement'");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);


      //evaluate internal forces and nonlinear contact stiffness
      nlnstiff(params,mydisp,stiffmatrix,fint);

  return 0;

}//Beam3contact::Evaluate


/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness and internal forces (private)                                               cyron 10/09|
 *-----------------------------------------------------------------------------------------------------------*/
void Beam3contact::nlnstiff(ParameterList& params,
              vector<double>&           disp,  //!< element displacement vector
              Epetra_SerialDenseMatrix& stiffmatrix,  //!< element stiffness matrix
              Epetra_SerialDenseVector& fint)  //!< element internal force vector
{
  
  //the following code is just an example how an assembly of fint and stiffmatrix might look like


  //compute internal forces (dummy)
  for(int k = 0; k<6; k++)
    fint(k) += 0;


  //compute stiffness matrix (dummy)
  for(int i = 0; i < 6; i++)
    for(int j = 0; j < 6; j++)
      stiffmatrix(i,j) += 0;




  return;
} //Beam3contact::nlnstiff


#endif  // #ifdef D_BEAM3
#endif  // #ifdef CCADISCRET 
