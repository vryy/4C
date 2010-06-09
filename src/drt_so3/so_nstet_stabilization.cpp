/*!----------------------------------------------------------------------*
\file so_nstet_stabilization.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_nstet.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "Epetra_SerialDenseSolver.h"

#if 0
#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/mooneyrivlin.H"
#endif

using namespace std;

/*----------------------------------------------------------------------*
 |  vol stabilization (protected)                              gee 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::VolDevStabLinear(
                        LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& defgrd,
                        LINALG::Matrix<NUMDOF_NSTET,1>* force)
{
  if (!force) return;

  //------------------------------------ get nodal pressure, J and defgrd
  const int numnode = NumNode();
  vector<double> npressure(numnode);
  vector<double> nJ(numnode);
  LINALG::Matrix<9,1> Fvec(false);
  vector<LINALG::Matrix<3,3> > nF(numnode);
  for (int i=0; i<numnode; ++i)
  {
    const int lid = ElementType().pnodecol_->Map().LID(Nodes()[i]->Id());
    npressure[i] = (*ElementType().pnodecol_)[lid];
    nJ[i]        = (*ElementType().Jnodecol_)[lid];
    //for (int j=0; j<9; ++j) Fvec(j) = (*(*myregister_->Fnodecol_)(j))[lid];
    //myregister_->VectortoMatrix(Fvec,nF[i]);
    //nF[i].Clear(); nF[i](0,0) = 1.0; nF[i](1,1) = 1.0; nF[i](2,2) = 1.0;
    //nF[i].Scale(pow(nJ[i],1./3.));
    //cout << nF[i];
  }

  //------------------------------------- compute average pressure
  //double averJ = 0.0;
  double averpressure = 0.0;
  for (int i=0; i<numnode; ++i)
  {
    //averJ += nJ[i];
    averpressure += npressure[i];
  }
  //averJ /= numnode;
  averpressure /= numnode;

  //------------------------------------------------------- get (p - P p)
  vector<double> nPpressure(npressure);
  for (int i=0; i<numnode; ++i)
    nPpressure[i] -= averpressure;
  //cout << averpressure << " linear " << nPpressure[0] << " " << nPpressure[1] << " " << nPpressure[2] << " " << nPpressure[3] << endl;

  //-------------------------------------------- prepare integration loop
  const double alpha  = (5.0 + 3.0*sqrt(5.0))/20.0;
  const double beta   = (5.0 - sqrt(5.0))/20.0;
  const double weight = 0.25;
  const double V = Volume();
  double xsi[4][4];
  xsi[0][0] = alpha;   xsi[0][1] = beta ;   xsi[0][2] = beta ;   xsi[0][3] = beta ;
  xsi[1][0] = beta ;   xsi[1][1] = alpha;   xsi[1][2] = beta ;   xsi[1][3] = beta ;
  xsi[2][0] = beta ;   xsi[2][1] = beta ;   xsi[2][2] = alpha;   xsi[2][3] = beta ;
  xsi[3][0] = beta ;   xsi[3][1] = beta ;   xsi[3][2] = beta ;   xsi[3][3] = alpha;
  const double fac = V * weight;
  LINALG::Matrix<NUMNOD_NSTET,1> funct;
  for (int gp=0; gp<4; ++gp)
  {
    ShapeFunction(funct,xsi[gp][0],xsi[gp][1],xsi[gp][2],xsi[gp][3]);
    // value of linear pressure distribution at this gaussian point
    double p = 0.0;
    double J = 0.0;
    LINALG::Matrix<3,3> F;
    F.Clear();
    for (int i=0; i<numnode; ++i)
    {
      p += nPpressure[i] * funct(i);
      J += nJ[i]         * funct(i);
      //F.Update(funct(i),nF[i],1.0);
    }
    F.Clear(); F(0,0) = 1.0; F(1,1) = 1.0; F(2,2) = 1.0;
    F.Scale(pow(J,1./3.));

    // righ cauchy green and its inverse
    LINALG::Matrix<3,3> cg;
    cg.MultiplyTN(F,F);
    cg.Invert();
    // stress = -p J C^-1
    cg.Scale(-p*J);
    LINALG::Matrix<6,1> stress(true);
    stress(0) = cg(0,0);
    stress(1) = cg(1,1);
    stress(2) = cg(2,2);
    stress(3) = cg(0,1);
    stress(4) = cg(1,2);
    stress(5) = cg(0,2);

    LINALG::Matrix<NUMSTR_NSTET,NUMDOF_NSTET> bop;
    for (int i=0; i<NUMNOD_NSTET; i++)
    {
      bop(0,NODDOF_NSTET*i+0) = F(0,0)*nxyz_(i,0);
      bop(0,NODDOF_NSTET*i+1) = F(1,0)*nxyz_(i,0);
      bop(0,NODDOF_NSTET*i+2) = F(2,0)*nxyz_(i,0);
      bop(1,NODDOF_NSTET*i+0) = F(0,1)*nxyz_(i,1);
      bop(1,NODDOF_NSTET*i+1) = F(1,1)*nxyz_(i,1);
      bop(1,NODDOF_NSTET*i+2) = F(2,1)*nxyz_(i,1);
      bop(2,NODDOF_NSTET*i+0) = F(0,2)*nxyz_(i,2);
      bop(2,NODDOF_NSTET*i+1) = F(1,2)*nxyz_(i,2);
      bop(2,NODDOF_NSTET*i+2) = F(2,2)*nxyz_(i,2);
      /* ~~~ */
      bop(3,NODDOF_NSTET*i+0) = F(0,0)*nxyz_(i,1) + F(0,1)*nxyz_(i,0);
      bop(3,NODDOF_NSTET*i+1) = F(1,0)*nxyz_(i,1) + F(1,1)*nxyz_(i,0);
      bop(3,NODDOF_NSTET*i+2) = F(2,0)*nxyz_(i,1) + F(2,1)*nxyz_(i,0);
      bop(4,NODDOF_NSTET*i+0) = F(0,1)*nxyz_(i,2) + F(0,2)*nxyz_(i,1);
      bop(4,NODDOF_NSTET*i+1) = F(1,1)*nxyz_(i,2) + F(1,2)*nxyz_(i,1);
      bop(4,NODDOF_NSTET*i+2) = F(2,1)*nxyz_(i,2) + F(2,2)*nxyz_(i,1);
      bop(5,NODDOF_NSTET*i+0) = F(0,2)*nxyz_(i,0) + F(0,0)*nxyz_(i,2);
      bop(5,NODDOF_NSTET*i+1) = F(1,2)*nxyz_(i,0) + F(1,0)*nxyz_(i,2);
      bop(5,NODDOF_NSTET*i+2) = F(2,2)*nxyz_(i,0) + F(2,0)*nxyz_(i,2);
    }

    force->MultiplyTN(fac,bop,stress,1.0);


  } // for (int gp=0; gp<4; ++gp)
  cout << *force;



  return;
}

/*----------------------------------------------------------------------*
 | voldev stabilization (protected)                           gee 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::VolDevStab(
                  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& cauchygreen,
                  LINALG::Matrix<NUMSTR_NSTET,1>& stress,
                  LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& cmat,
                  double& pressure)
{
  LINALG::Matrix<6,6> cmatdev;
  LINALG::Matrix<6,1> stressdev;

  // compute deviatoric stress and tangent
  DevStressTangent(stressdev,cmatdev,cmat,stress,cauchygreen,pressure);

  // reduce deviatoric stresses
  stress.Update(-ALPHA_NSTET,stressdev,1.0);
  // reduce deviatoric tangent
  cmat.Update(-ALPHA_NSTET,cmatdev,1.0);
  return;
}


/*----------------------------------------------------------------------*
 |  voldev stabilization (protected)                           gee 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::VolDevStab(
                  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& defgrd,
                  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& cauchygreen,
                  LINALG::Matrix<NUMSTR_NSTET,1>& stress,
                  LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& cmat,
                  double& density)
{
#if 1
  // do deviatoric F, C, E
  const double J = defgrd.Determinant();
  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> Cbar(cauchygreen);
  Cbar.Scale(pow(J,-2./3.));
  LINALG::Matrix<6,1> glstrainbar(false);
  glstrainbar(0) = 0.5 * (Cbar(0,0) - 1.0);
  glstrainbar(1) = 0.5 * (Cbar(1,1) - 1.0);
  glstrainbar(2) = 0.5 * (Cbar(2,2) - 1.0);
  glstrainbar(3) = Cbar(0,1);
  glstrainbar(4) = Cbar(1,2);
  glstrainbar(5) = Cbar(2,0);
  LINALG::Matrix<3,3> Fbar(false);
  Fbar.SetCopy(defgrd.A());
  Fbar.Scale(pow(J,-1./3.));

  //SelectMaterial(stress,cmat,density,glstrainbar,Fbar,0);
  SelectMaterial(stress,cmat,density,glstrainbar,Fbar,0);

  // define stuff we need to do the split
  LINALG::Matrix<6,6> cmatdev;
  LINALG::Matrix<6,1> stressdev;

  // do separation of deviatoric/volumetric components
  double elepressure = 0.0;
  NStetType::DevStressTangent(stressdev,cmatdev,cmat,stress,
                              cauchygreen,elepressure);
  stress.Update(ALPHA_NSTET,stressdev,0.0);
  cmat.Update(ALPHA_NSTET,cmatdev,0.0);

#else
  double averpressure = 0.0;
  double averJ = 0.0;
  const int numnode = NumNode();
  for (int i=0; i<numnode; ++i)
  {
    int lid = myregister_->pnodecol_->Map().LID(Nodes()[i]->Id());
    averpressure += (*myregister_->pnodecol_)[lid];
    averJ        += (*myregister_->Jnodecol_)[lid];
  }
  //printf("\n"); fflush(stdout);
  averpressure /= numnode;
  averJ        /= numnode;

  // do deviatoric F, C, E and replace the volumetric component accordingly
  const double J = defgrd.Determinant();
  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> Cbar(cauchygreen);
  Cbar.Scale(pow(J,-2./3.));
  Cbar.Scale(pow(averJ,2./3.));
  //printf("Owner %d Id %d Nodal/Ele/Delta J %10.5e %10.5e       %10.5e\n",Owner(),Id(),averJ,J,averJ-J);

  LINALG::Matrix<3,3> Fbar(false);
  Fbar.SetCopy(defgrd.A());
  Fbar.Scale(pow(J,-1./3.));
  Fbar.Scale(pow(averJ,1./3.));

  LINALG::Matrix<6,1> glstrainbar(false);
  glstrainbar(0) = 0.5 * (Cbar(0,0) - 1.0);
  glstrainbar(1) = 0.5 * (Cbar(1,1) - 1.0);
  glstrainbar(2) = 0.5 * (Cbar(2,2) - 1.0);
  glstrainbar(3) = Cbar(0,1);
  glstrainbar(4) = Cbar(1,2);
  glstrainbar(5) = Cbar(2,0);

  //SelectMaterial(stress,cmat,density,glstrainbar,Fbar,0);
  SelectMaterial(stress,cmat,density,glstrainbar,Fbar,0);
  //SelectMaterial(stress,cmat,density,glstrain,defgrd,0);

  stress.Scale(ALPHA_NSTET);
  cmat.Scale(ALPHA_NSTET);

  // define stuff we need to do the split
  LINALG::Matrix<6,6> cmatdev;
  LINALG::Matrix<6,6> cmatvol;
  LINALG::Matrix<6,1> stressdev;
  LINALG::Matrix<6,1> stressvol;

  // do separation of deviatoric/volumetric components
  double elepressure = 0.0;
  NStetType::DevVolStressTangent(stressdev,stressvol,
                                 cmatdev,cmatvol,cmat,stress,
                                 cauchygreen,elepressure);

  stress.Update(ALPHA_NSTET,stressdev,0.0);
  cmat.Update(ALPHA_NSTET,cmatdev,0.0);

  stress.Update(ALPHA_NSTET,stressvol,1.0);
  cmat.Update(ALPHA_NSTET,cmatvol,1.0);

  //printf("Nodal/Ele/Delta %10.5e %10.5e       %10.5e\n",averpressure,elepressure,averpressure-elepressure);
  //printf("Nodal/Ele/Delta J %10.5e %10.5e       %10.5e\n",averJ,J,averJ-J);
#endif
  return;
}

/*----------------------------------------------------------------------*
 | voldev stabilization (protected)                           gee 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStetType::DevStab(
                  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& cauchygreen,
                  LINALG::Matrix<NUMSTR_NSTET,1>& stress,
                  LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& cmat)
{
  LINALG::Matrix<6,6> cmatdev;
  LINALG::Matrix<6,1> stressdev;

  // compute deviatoric stress and tangent
  double pressure = 0.0;
  DevStressTangent(stressdev,cmatdev,cmat,stress,cauchygreen,pressure);

  // reduce deviatoric stresses
  stress.Update(-ALPHA_NSTET,stressdev,1.0);
  // reduce deviatoric tangent
  cmat.Update(-ALPHA_NSTET,cmatdev,1.0);
  return;
}

/*----------------------------------------------------------------------*
 |  vol stabilization (protected)                              gee 05/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NStet::DevStab(
                  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& defgrd,
                  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET>& cauchygreen,
                  LINALG::Matrix<NUMSTR_NSTET,1>& stress,
                  LINALG::Matrix<NUMSTR_NSTET,NUMSTR_NSTET>& cmat,
                  double& density
                                  )
{
  // do deviatoric F, C, E
  const double J = defgrd.Determinant();
  LINALG::Matrix<NUMDIM_NSTET,NUMDIM_NSTET> Cbar(cauchygreen);
  Cbar.Scale(pow(J,-2./3.));
  LINALG::Matrix<6,1> glstrainbar(false);
  glstrainbar(0) = 0.5 * (Cbar(0,0) - 1.0);
  glstrainbar(1) = 0.5 * (Cbar(1,1) - 1.0);
  glstrainbar(2) = 0.5 * (Cbar(2,2) - 1.0);
  glstrainbar(3) = Cbar(0,1);
  glstrainbar(4) = Cbar(1,2);
  glstrainbar(5) = Cbar(2,0);
  LINALG::Matrix<3,3> Fbar(false);
  Fbar.SetCopy(defgrd.A());
  Fbar.Scale(pow(J,-1./3.));

  //SelectMaterial(stress,cmat,density,glstrainbar,Fbar,0);
  SelectMaterial(stress,cmat,density,glstrainbar,Fbar,0);

  // define stuff we need to do the split
  LINALG::Matrix<6,6> cmatdev;
  LINALG::Matrix<6,1> stressdev;

  // do separation of deviatoric/volumetric components
  double elepressure = 0.0;
  NStetType::DevStressTangent(stressdev,cmatdev,cmat,stress,
                              cauchygreen,elepressure);
  stress.Update(ALPHA_NSTET,stressdev,0.0);
  cmat.Update(ALPHA_NSTET,cmatdev,0.0);


  return;
}






#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
