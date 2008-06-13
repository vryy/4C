 /*!-------------------------------------------------------------------
\file drt_potential_manager.cpp

\brief Class controlling surface stresses due to potential forces
between interfaces
and containing all necessary history data

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <blitz/array.h>
#include "drt_potential_manager.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_xfem/xfsi_searchtree.H"
#include <cstdlib>


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 06/08|
 *-------------------------------------------------------------------*/
DRT::PotentialManager::PotentialManager(DRT::Discretization& discret):
discret_(discret)
{
  surfrowmap_ = DRT::UTILS::GeometryElementMap(discret, "Potential", false);
  RCP<Epetra_Map> surfcolmap = DRT::UTILS::GeometryElementMap(discret, "Potential", true);

  // we start with interfacial area (A_old_) and concentration
  // (con_quot_) = 0. this is wrong but does not make a difference
  // since we apply the equilibrium concentration gradually, thus we
  // do not need these history variables needed for the dynamic model
  // in the beginning.
  A_old_temp_   = rcp(new Epetra_Vector(*surfcolmap,true));
  A_old_        = rcp(new Epetra_Vector(*surfcolmap,true));
  xTree_        = rcp(new XFEM::XSearchTree());
}


/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/

void DRT::PotentialManager::EvaluatePotential(  ParameterList& p,
                                                RefCountPtr<Epetra_Vector> disp,
                                                RefCountPtr<Epetra_Vector> fint,
                                                RefCountPtr<LINALG::SparseMatrix> stiff)
{
    // action for elements
    p.set("action","calc_potential_stiff");

    discret_.ClearState();
    discret_.SetState("displacement",disp);
    discret_.EvaluateCondition(p,stiff,null,fint,null,null,"Potential");

    return;
}

/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| update surface area                                                |
*--------------------------------------------------------------------*/

void DRT::PotentialManager::Update()
{
  A_old_->Update(1.0, *A_old_temp_, 0.0);
}

/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| write restart                                                      |
*--------------------------------------------------------------------*/

void DRT::PotentialManager::GetHistory(RCP<Epetra_Vector> A_old_temp_row)
{
  // Note that the temporal vectors need to be written since in the
  // final ones we still have the data of the former step. The column
  // map based vector used for calculations is exported to a row map
  // based one needed for writing.

  LINALG::Export(*A_old_temp_, *A_old_temp_row);
}



/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level for Lennard-Jones potential interaction forces    |
*--------------------------------------------------------------------*/

void DRT::PotentialManager::StiffnessAndInternalForces(const int curvenum,
                                                        const double& A,
                                                        Epetra_SerialDenseVector& fint,
                                                        Epetra_SerialDenseMatrix& K_surf,
                                                        const int ID,
                                                        const double time,
                                                        const double dt,
                                                        const double label,
                                                        const double depth,
                                                        const double rootDist,
                                                        const double cutOff)
{

  double traction;
  int LID = A_old_->Map().LID(ID);
  double t_end = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).end();

  dserror("not yet implemented");
  /*------------------------------------------------- initialization */
  (*A_old_temp_)[LID] = A;

  /*-----------calculation of current surface stress and its partial
   *-----------------derivative with respect to the interfacial area */

  if (time <= t_end)         /* gradual application of surface stress */
  {
    traction = 0;
  }
  else
  {
    traction = 1;
  }

  double curvefac = 1.;

  /*------------gradual application of surface stresses via time curve*/
  if (time <= t_end)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  double ndof = 5; //Adiff.Length();

  for (int i=0;i<ndof;++i)
    for (int j=0;j<ndof;++j)
      K_surf(i,j) = curvefac;
  
  /*------calculation of current internal force due to surface energy*/
  for (int i=0;i<ndof;++i)
    fint[i] = curvefac;


  return;
}

#endif

