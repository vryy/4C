/*!----------------------------------------------------------------------
\file inversedesign_evaluate.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)

#include "inversedesign.H"
#include "so_hex8.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  hex8 integration method                                    (public) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::soh8_nlnstiffmass(
      DRT::ELEMENTS::So_hex8*   ele,            ///< this element
      vector<int>&              lm,             ///< location matrix
      vector<double>&           disp,           ///< current displacements
      vector<double>&           residual,       ///< current residuum
      Epetra_SerialDenseMatrix* stiffmatrix,    ///< element stiffness matrix
      Epetra_SerialDenseMatrix* massmatrix,     ///< element mass matrix
      Epetra_SerialDenseVector* force,          ///< element internal force vector
      Epetra_SerialDenseMatrix* elestress,      ///< stresses at GP
      Epetra_SerialDenseMatrix* elestrain,      ///< strains at GP
      ParameterList&            params,         ///< algorithmic parameters e.g. time
      const bool                cauchy,         ///< stress output option
      const bool                euler_almansi)  ///< strain output option
{
   const static DRT::ELEMENTS::So_hex8::Integrator_So_hex8 int_hex8;
  
  //---------------------------------------------------------------------
  // element geometry (note that this is inverse!)
  LINALG::SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8);  // material coord. of element
  LINALG::SerialDenseMatrix xcurr(NUMNOD_SOH8,NUMDIM_SOH8);  // current  coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    xcurr(i,0) = ele->Nodes()[i]->X()[0];
    xcurr(i,1) = ele->Nodes()[i]->X()[1];
    xcurr(i,2) = ele->Nodes()[i]->X()[2];

    xrefe(i,0) = xcurr(i,0) + disp[i*NODDOF_SOH8+0];
    xrefe(i,1) = xcurr(i,1) + disp[i*NODDOF_SOH8+1];
    xrefe(i,2) = xcurr(i,2) + disp[i*NODDOF_SOH8+2];
  }

#if 0
  if (ele->Id()==63)
  {
    cout << "\nxcurr\n" << xcurr;
    
    cout << "\ndisp\n";
    for (int i=0; i<24; ++i) printf("%10.5e ",disp[i]);
    cout << "\n";
    
    cout << "\nxrefe\n" << xrefe;
  }
#endif  
  
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  // loop gaussian points
  for (int gp=0; gp<NUMGPT_SOH8; ++gp) 
  {
    //----------------------------------- get inverse of Jacobian mapping
    const double              detj = ele->detJ_[gp];
    Epetra_SerialDenseMatrix& invj = ele->invJ_[gp];
    
    //------------------------------ compute derivs wrt to spatial coords
    LINALG::SerialDenseMatrix n_xyz(NUMDIM_SOH8,NUMNOD_SOH8);
    n_xyz.Multiply('N','N',1.0,invj,int_hex8.deriv_gp[gp],0.0);
  
    //--------------------------- build defgrd of inverse mapping dX / dx
    LINALG::SerialDenseMatrix f(NUMDIM_SOH8,NUMDIM_SOH8);
    f.Multiply('T','T',1.0,xrefe,n_xyz,0.0);
    
    //--------------------------- build defgrd of forward mapping dx / dX
    LINALG::SerialDenseMatrix F(f);
    const double j = LINALG::NonsymInverse3x3(F);
    
    //-------------------------------------------- build determinant of F
    const double J = 1.0/j;

    //------------------------------------ build F and F^T as vectors 9x1
    LINALG::SerialDenseVector Fvec(9);
    Fvec(0) = F(0,0);
    Fvec(1) = F(1,0);
    Fvec(2) = F(2,0);
    Fvec(3) = F(0,1);
    Fvec(4) = F(1,1);
    Fvec(5) = F(2,1);
    Fvec(6) = F(0,2);
    Fvec(7) = F(1,2);
    Fvec(8) = F(2,2);
    
    LINALG::SerialDenseVector FTvec(9);
    FTvec(0) = F(0,0);
    FTvec(1) = F(0,1);
    FTvec(2) = F(0,2);
    FTvec(3) = F(1,0);
    FTvec(4) = F(1,1);
    FTvec(5) = F(1,2);
    FTvec(6) = F(2,0);
    FTvec(7) = F(2,1);
    FTvec(8) = F(2,2);
    
    //--------------------------------------------- build operator Lambda
    LINALG::SerialDenseMatrix Lambda(9,9);
    double Lambda4[3][3][3][3];
    BuildLambda(Lambda4,Lambda,F);
#if 0
    if (ele->Id()==63 && !gp)
    {
      LINALG::SerialDenseMatrix FD_Lambda(9,9);
      FDLambda(FD_Lambda,f); // FD did include the minus
      //FD_Lambda += Lambda;
      //const double norm = FD_Lambda.NormInf();
      //printf("Delta Lambda     norm %15.10e\n",norm);
      //cout << "\n Delta FD_Lambda \n" << FD_Lambda;
      //cout << "\n Lambda \n" << Lambda;
      //FD_Lambda.Scale(-1.0);
      //Lambda = FD_Lambda;
    }
#endif    
    
    LINALG::SerialDenseMatrix LambdaT(9,9);
    BuildLambdaT(LambdaT,F);
#if 0
    if (ele->Id()==63 && !gp)
    {
      LINALG::SerialDenseMatrix FD_LambdaT(9,9);
      FDLambdaT(FD_LambdaT,f);
      FD_LambdaT += LambdaT;
      const double norm = FD_LambdaT.NormInf();
      printf("Delta LambdaT norm %15.10e\n",norm);
      //cout << "\n Delta FD_LambdaT \n" << FD_LambdaT;
      //cout << "\n Lambda \n" << Lambda;
    }
#endif    
    
    //------------------------------------------------- build operator IF
    // this has been analytically tested: IF*S == F S F^T and ok
    LINALG::SerialDenseMatrix IF(6,6);
    BuildIF(IF,F);
    
    //---------------------------------------------- build operator Theta
    double Theta4[3][3][3][3];
    LINALG::SerialDenseMatrix Theta(6,9);
    BuildTheta(Theta4,Theta,F);
#if 0
    LINALG::SerialDenseMatrix FD_Theta(6,9);
    double Theta4FD[3][3][3][3];
    FDTheta(Theta4FD,FD_Theta,F);
    if (ele->Id()==63 && !gp)
    {
      //FD_Theta.Scale(-1.0);
      //FD_Theta += Theta;
      //const double norm = FD_Theta.NormInf();
      //printf("Delta FD_Theta norm %15.10e\n\n",norm);
      cout << "\n FD_Theta \n" << FD_Theta;
      cout << "\n Theta \n" << Theta;
    }
    Theta = FD_Theta;
#endif    
    
    //--------------- build right Cauchy-Green and Green-Lagrange strains
    LINALG::SerialDenseMatrix cauchygreen(NUMDIM_SOH8,NUMDIM_SOH8);
    cauchygreen.Multiply('T','N',1.0,F,F,0.0);
    
    //-- Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    LINALG::SerialDenseVector glstrain(NUMSTR_SOH8);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

    //---- build B-operator (wrt to spatial, that is known configuration)
    LINALG::SerialDenseMatrix B(NUMSTR_SOH8,NUMDOF_SOH8);
    for (int i=0; i<NUMNOD_SOH8; ++i) 
    {
      B(0,NODDOF_SOH8*i+0) = n_xyz(0,i);
      B(0,NODDOF_SOH8*i+1) = 0.0;
      B(0,NODDOF_SOH8*i+2) = 0.0;

      B(1,NODDOF_SOH8*i+0) = 0.0;
      B(1,NODDOF_SOH8*i+1) = n_xyz(1,i);
      B(1,NODDOF_SOH8*i+2) = 0.0;

      B(2,NODDOF_SOH8*i+0) = 0.0;
      B(2,NODDOF_SOH8*i+1) = 0.0;
      B(2,NODDOF_SOH8*i+2) = n_xyz(2,i);

      B(3,NODDOF_SOH8*i+0) = n_xyz(1,i);
      B(3,NODDOF_SOH8*i+1) = n_xyz(0,i);
      B(3,NODDOF_SOH8*i+2) = 0.0;

      B(4,NODDOF_SOH8*i+0) = 0.0;
      B(4,NODDOF_SOH8*i+1) = n_xyz(2,i);
      B(4,NODDOF_SOH8*i+2) = n_xyz(1,i);

      B(5,NODDOF_SOH8*i+0) = n_xyz(2,i);
      B(5,NODDOF_SOH8*i+1) = 0.0;
      B(5,NODDOF_SOH8*i+2) = n_xyz(0,i);
    }
    
    //--------------------------- build N_x operator (wrt spatial config)
    Epetra_SerialDenseMatrix N_x(9,NUMDOF_SOH8);
    for (int i=0; i<NUMNOD_SOH8; ++i) 
    {
      N_x(0,3*i+0) = n_xyz(0,i);
      N_x(1,3*i+1) = n_xyz(0,i);
      N_x(2,3*i+2) = n_xyz(0,i);
      
      N_x(3,3*i+0) = n_xyz(1,i);
      N_x(4,3*i+1) = n_xyz(1,i);
      N_x(5,3*i+2) = n_xyz(1,i);
      
      N_x(6,3*i+0) = n_xyz(2,i);
      N_x(7,3*i+1) = n_xyz(2,i);
      N_x(8,3*i+2) = n_xyz(2,i);
    }
    
    //------------------------------------------------- call material law
    Epetra_SerialDenseMatrix cmat(NUMSTR_SOH8,NUMSTR_SOH8);
    Epetra_SerialDenseVector stress(NUMSTR_SOH8);
    double density;
    ele->soh8_mat_sel(&stress,&cmat,&density,&glstrain,&F,gp,params);
        
    //------------------------------------------- compute cauchy stresses
    // (IF has been tested ok against jFSF^T)
    Epetra_SerialDenseVector cstress(NUMSTR_SOH8);
    cstress.Multiply('N','N',j,IF,stress,0.0);
    
    //-------------------------------------------- build operator Ypsilon
    LINALG::SerialDenseMatrix S(3,3);
    S(0,0) = stress(0);
    S(0,1) = stress(3);
    S(0,2) = stress(5);
    S(1,0) = stress(3);
    S(1,1) = stress(1);
    S(1,2) = stress(4);
    S(2,0) = stress(5);
    S(2,1) = stress(4);
    S(2,2) = stress(2);
    LINALG::SerialDenseMatrix Ypsilon(6,9);
    double Ypsilon4[3][3][3][3];
    BuildYpsilon(Ypsilon4,Ypsilon,F,S);
#if 0
    LINALG::SerialDenseMatrix FD_Ypsilon(6,9);
    double Ypsilon4FD[3][3][3][3];
    if (ele->Id()==63 && !gp)
    {
      FDYpsilon(Ypsilon4FD,FD_Ypsilon,F,S);
      //FD_Ypsilon.Scale(-1.0);
      //FD_Ypsilon += Ypsilon;
      //const double norm = FD_Ypsilon.NormInf();
      //printf("Delta FD_Ypsilon norm %15.10e\n",norm);
      cout << "\n FD_Ypsilon \n" << FD_Ypsilon;
      cout << "\n Ypsilon \n" << Ypsilon;
    }
    //Ypsilon = FD_Ypsilon;
#endif    
    
    //----------------- integration factor dV (spatial) * Gaussian weight
    const double intfac = detj * int_hex8.weights(gp);

    // assemble internal forces
    if (force)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      (*force).Multiply('T','N',intfac,B,cstress,1.0);
    }
    
    //*******************************************************************
    if (stiffmatrix)
    {
      Epetra_SerialDenseMatrix sum(6,9);
      Epetra_SerialDenseMatrix tmp2(1,9);

      //================================== 1.) B^T * cstress * F^T * N_x
#if 1
      {
        sum.Multiply('N','T',1.0,cstress,FTvec,0.0);
        LINALG::SerialDenseMatrix tmp3(6,NUMDOF_SOH8);
        tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
#if 0 // FD of dj / dX
        if (ele->Id()==63 && !gp)
          cout << "\n analytical \n" << tmp3;
        if (ele->Id()==63 && !gp)
        {
          Epetra_SerialDenseMatrix djdX(1,NUMDOF_SOH8);
          FD_djdX(djdX,disp,gp,ele);
          Epetra_SerialDenseMatrix tmp(6,NUMDOF_SOH8);
          tmp.Multiply('N','N',1./j,cstress,djdX,0.0);
          cout << "\n finite diff \n" << tmp;
          //(*stiffmatrix).Multiply('T','N',intfac,B,tmp,1.0);
        }
#endif
        (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);
      }
#endif
      
      //============================== 2.) -j * B^T * S * F^T * Lambda * N_x
#if 0
      {
        tmp2.Multiply('T','N',1.0,Fvec,Lambda,0.0);
        sum.Multiply('N','N',-j,stress,tmp2,0.0);
        LINALG::SerialDenseMatrix tmp3(6,NUMDOF_SOH8);
        tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
        (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);
      }
#endif
      
      //=========================== 4.) -j * B^T * S * F^T * Lambda^T * N_x     
#if 0
      {
        tmp2.Multiply('T','N',1.0,FTvec,LambdaT,0.0);
        sum.Multiply('N','N',-j,stress,tmp2,0.0);
        LINALG::SerialDenseMatrix tmp3(6,NUMDOF_SOH8);
        tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
        (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);
      }
#endif

      //=================================== -j * B^T * Ypsilon * Lambda * N_x
      // Anstatt 2.) und 4.):
#if 1
      {
        sum.Multiply('N','N',-j,Ypsilon,Lambda,0.0);
        //if (ele->Id()==63 && !gp) cout << "\n analytical Voigt \n"  << sum;
        //TensorMultiply(sum,-j,Ypsilon4,Lambda4);
        //if (ele->Id()==63 && !gp) cout << "\n tensor Voigt \n"  << sum;
        LINALG::SerialDenseMatrix tmp3(6,NUMDOF_SOH8);
#if 0 // FD of dIF/df * S
        //if (ele->Id()==63 && !gp)
        {
          Epetra_SerialDenseMatrix dISdf(6,9);
          FD_dISdf(dISdf,stress,f);
          dISdf.Scale(j);
          if (ele->Id()==63 && !gp) cout << "\n analytical \n"  << sum;
          if (ele->Id()==63 && !gp) cout << "\n finite diff \n" << dISdf;
          tmp3.Multiply('N','N',1.0,dISdf,N_x,0.0);
        }
#endif
        tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
#if 0 // FD of dIF/dX * S
        if (ele->Id()==63 && !gp) cout << "\n analytical \n" << tmp3;
        //if (ele->Id()==63 && !gp)
        {
          Epetra_SerialDenseMatrix dISdX(6,NUMDOF_SOH8);
          FD_dISdX(dISdX,disp,gp,ele,params);
          if (ele->Id()==63 && !gp) cout << "\n fd2 \n" << dISdX;
          (*stiffmatrix).Multiply('T','N',intfac,B,dISdX,1.0);
        }
#endif
        (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);
      }
#endif
      
      //======================= 3.) -j * B^T * IF * cmat * Theta * Lambda * N_x
#if 1
      {
        LINALG::SerialDenseMatrix sixnine(6,9);
        LINALG::SerialDenseMatrix sixsix(6,6);
        sixnine.Multiply('N','N',1.0,Theta,Lambda,0.0);
        sixsix.Multiply('N','N',1.0,IF,cmat,0.0);
        sum.Multiply('N','N',-j,sixsix,sixnine,0.0);
        LINALG::SerialDenseMatrix tmp3(6,NUMDOF_SOH8);
        tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
        (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);
      }
#endif



#if 0
      //================ put everything together: K = intfac * B*T * sum * N_x
      {
        LINALG::SerialDenseMatrix tmp3(6,NUMDOF_SOH8);
        tmp3.Multiply('N','N',1.0,sum,N_x,0.0);
#if 0 // FD of whole matrix
        Epetra_SerialDenseMatrix fdstiff(6,NUMDOF_SOH8);
        FDstiffmatrix(fdstiff,disp,gp,ele,params);
        if (ele->Id()==63 && !gp)
        {
          //cout << "\n analytical \n" << tmp3;
          //cout << "\n finite diff \n" << fdstiff;
          //fdstiff.Scale(-1.0);
          //fdstiff += tmp3;
          //double norm = fdstiff.NormInf();
          //printf("Delta Norm %15.10e\n",norm);
        }
        (*stiffmatrix).Multiply('T','N',intfac,B,fdstiff,1.0);
#endif
        (*stiffmatrix).Multiply('T','N',intfac,B,tmp3,1.0);
      }
#endif


    } // if (stiffmatrix)
    //*******************************************************************

    
    //*******************************************************************
    // Strictly, inverse design analysis is stationary and should not have
    // a mass term. Loosely, if no Dirichlet-BCs are present a small
    // mass term might be used. Note that we use density by unit deformed
    // volume here!
    if (massmatrix)
    {
      const double fac = density * detj * int_hex8.weights(gp);
      for (int inod=0; inod<NUMNOD_SOH8; ++inod) 
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod) 
        {
          const double massfactor = (int_hex8.shapefct_gp[gp])(inod) 
                                  * (int_hex8.shapefct_gp[gp])(jnod)
                                  * fac;     
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
    } // if (massmatrix)
    //*******************************************************************

  } // for (int gp=0; gp<NUMGPT_SOH8; ++gp) 
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------
  
  
} // DRT::ELEMENTS::InvDesign::soh8_nlnstiffmass


/*----------------------------------------------------------------------*
 |  Lambda tensor                                             (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildLambda(double Lambda4[][3][3][3],
                                           LINALG::SerialDenseMatrix& L, 
                                           const LINALG::SerialDenseMatrix& F) const

{
  //double Lambda4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int p=0; p<3; ++p)
      for (int q=0; q<3; ++q)
        for (int m=0; m<3; ++m)
          Lambda4[k][m][p][q] = F(k,p)*F(q,m);
  
  L(0,0) = Lambda4[0][0][0][0];
  L(0,1) = Lambda4[0][0][1][0];
  L(0,2) = Lambda4[0][0][2][0];
  L(0,3) = Lambda4[0][0][0][1];
  L(0,4) = Lambda4[0][0][1][1];
  L(0,5) = Lambda4[0][0][2][1];
  L(0,6) = Lambda4[0][0][0][2];
  L(0,7) = Lambda4[0][0][1][2];
  L(0,8) = Lambda4[0][0][2][2];
  
  L(1,0) = Lambda4[1][0][0][0];
  L(1,1) = Lambda4[1][0][1][0];
  L(1,2) = Lambda4[1][0][2][0];
  L(1,3) = Lambda4[1][0][0][1];
  L(1,4) = Lambda4[1][0][1][1];
  L(1,5) = Lambda4[1][0][2][1];
  L(1,6) = Lambda4[1][0][0][2];
  L(1,7) = Lambda4[1][0][1][2];
  L(1,8) = Lambda4[1][0][2][2];
  
  L(2,0) = Lambda4[2][0][0][0];
  L(2,1) = Lambda4[2][0][1][0];
  L(2,2) = Lambda4[2][0][2][0];
  L(2,3) = Lambda4[2][0][0][1];
  L(2,4) = Lambda4[2][0][1][1];
  L(2,5) = Lambda4[2][0][2][1];
  L(2,6) = Lambda4[2][0][0][2];
  L(2,7) = Lambda4[2][0][1][2];
  L(2,8) = Lambda4[2][0][2][2];
  
  L(3,0) = Lambda4[0][1][0][0];
  L(3,1) = Lambda4[0][1][1][0];
  L(3,2) = Lambda4[0][1][2][0];
  L(3,3) = Lambda4[0][1][0][1];
  L(3,4) = Lambda4[0][1][1][1];
  L(3,5) = Lambda4[0][1][2][1];
  L(3,6) = Lambda4[0][1][0][2];
  L(3,7) = Lambda4[0][1][1][2];
  L(3,8) = Lambda4[0][1][2][2];
  
  L(4,0) = Lambda4[1][1][0][0];
  L(4,1) = Lambda4[1][1][1][0];
  L(4,2) = Lambda4[1][1][2][0];
  L(4,3) = Lambda4[1][1][0][1];
  L(4,4) = Lambda4[1][1][1][1];
  L(4,5) = Lambda4[1][1][2][1];
  L(4,6) = Lambda4[1][1][0][2];
  L(4,7) = Lambda4[1][1][1][2];
  L(4,8) = Lambda4[1][1][2][2];
  
  L(5,0) = Lambda4[2][1][0][0];
  L(5,1) = Lambda4[2][1][1][0];
  L(5,2) = Lambda4[2][1][2][0];
  L(5,3) = Lambda4[2][1][0][1];
  L(5,4) = Lambda4[2][1][1][1];
  L(5,5) = Lambda4[2][1][2][1];
  L(5,6) = Lambda4[2][1][0][2];
  L(5,7) = Lambda4[2][1][1][2];
  L(5,8) = Lambda4[2][1][2][2];
  
  L(6,0) = Lambda4[0][2][0][0];
  L(6,1) = Lambda4[0][2][1][0];
  L(6,2) = Lambda4[0][2][2][0];
  L(6,3) = Lambda4[0][2][0][1];
  L(6,4) = Lambda4[0][2][1][1];
  L(6,5) = Lambda4[0][2][2][1];
  L(6,6) = Lambda4[0][2][0][2];
  L(6,7) = Lambda4[0][2][1][2];
  L(6,8) = Lambda4[0][2][2][2];
  
  L(7,0) = Lambda4[1][2][0][0];
  L(7,1) = Lambda4[1][2][1][0];
  L(7,2) = Lambda4[1][2][2][0];
  L(7,3) = Lambda4[1][2][0][1];
  L(7,4) = Lambda4[1][2][1][1];
  L(7,5) = Lambda4[1][2][2][1];
  L(7,6) = Lambda4[1][2][0][2];
  L(7,7) = Lambda4[1][2][1][2];
  L(7,8) = Lambda4[1][2][2][2];
  
  L(8,0) = Lambda4[2][2][0][0];
  L(8,1) = Lambda4[2][2][1][0];
  L(8,2) = Lambda4[2][2][2][0];
  L(8,3) = Lambda4[2][2][0][1];
  L(8,4) = Lambda4[2][2][1][1];
  L(8,5) = Lambda4[2][2][2][1];
  L(8,6) = Lambda4[2][2][0][2];
  L(8,7) = Lambda4[2][2][1][2];
  L(8,8) = Lambda4[2][2][2][2];
  
  return; 
} // DRT::ELEMENTS::InvDesign::BuildLambda

/*----------------------------------------------------------------------*
 |  Lambda tensor                                             (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildLambdaT(LINALG::SerialDenseMatrix& L, 
                                           const LINALG::SerialDenseMatrix& F) const

{
  double Lambda4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int p=0; p<3; ++p)
      for (int q=0; q<3; ++q)
        for (int m=0; m<3; ++m)
          Lambda4[k][m][p][q] = F(m,p)*F(q,k);
  
  L(0,0) = Lambda4[0][0][0][0];
  L(0,1) = Lambda4[0][0][1][0];
  L(0,2) = Lambda4[0][0][2][0];
  L(0,3) = Lambda4[0][0][0][1];
  L(0,4) = Lambda4[0][0][1][1];
  L(0,5) = Lambda4[0][0][2][1];
  L(0,6) = Lambda4[0][0][0][2];
  L(0,7) = Lambda4[0][0][1][2];
  L(0,8) = Lambda4[0][0][2][2];
  
  L(1,0) = Lambda4[1][0][0][0];
  L(1,1) = Lambda4[1][0][1][0];
  L(1,2) = Lambda4[1][0][2][0];
  L(1,3) = Lambda4[1][0][0][1];
  L(1,4) = Lambda4[1][0][1][1];
  L(1,5) = Lambda4[1][0][2][1];
  L(1,6) = Lambda4[1][0][0][2];
  L(1,7) = Lambda4[1][0][1][2];
  L(1,8) = Lambda4[1][0][2][2];
  
  L(2,0) = Lambda4[2][0][0][0];
  L(2,1) = Lambda4[2][0][1][0];
  L(2,2) = Lambda4[2][0][2][0];
  L(2,3) = Lambda4[2][0][0][1];
  L(2,4) = Lambda4[2][0][1][1];
  L(2,5) = Lambda4[2][0][2][1];
  L(2,6) = Lambda4[2][0][0][2];
  L(2,7) = Lambda4[2][0][1][2];
  L(2,8) = Lambda4[2][0][2][2];
  
  L(3,0) = Lambda4[0][1][0][0];
  L(3,1) = Lambda4[0][1][1][0];
  L(3,2) = Lambda4[0][1][2][0];
  L(3,3) = Lambda4[0][1][0][1];
  L(3,4) = Lambda4[0][1][1][1];
  L(3,5) = Lambda4[0][1][2][1];
  L(3,6) = Lambda4[0][1][0][2];
  L(3,7) = Lambda4[0][1][1][2];
  L(3,8) = Lambda4[0][1][2][2];
  
  L(4,0) = Lambda4[1][1][0][0];
  L(4,1) = Lambda4[1][1][1][0];
  L(4,2) = Lambda4[1][1][2][0];
  L(4,3) = Lambda4[1][1][0][1];
  L(4,4) = Lambda4[1][1][1][1];
  L(4,5) = Lambda4[1][1][2][1];
  L(4,6) = Lambda4[1][1][0][2];
  L(4,7) = Lambda4[1][1][1][2];
  L(4,8) = Lambda4[1][1][2][2];
  
  L(5,0) = Lambda4[2][1][0][0];
  L(5,1) = Lambda4[2][1][1][0];
  L(5,2) = Lambda4[2][1][2][0];
  L(5,3) = Lambda4[2][1][0][1];
  L(5,4) = Lambda4[2][1][1][1];
  L(5,5) = Lambda4[2][1][2][1];
  L(5,6) = Lambda4[2][1][0][2];
  L(5,7) = Lambda4[2][1][1][2];
  L(5,8) = Lambda4[2][1][2][2];
  
  L(6,0) = Lambda4[0][2][0][0];
  L(6,1) = Lambda4[0][2][1][0];
  L(6,2) = Lambda4[0][2][2][0];
  L(6,3) = Lambda4[0][2][0][1];
  L(6,4) = Lambda4[0][2][1][1];
  L(6,5) = Lambda4[0][2][2][1];
  L(6,6) = Lambda4[0][2][0][2];
  L(6,7) = Lambda4[0][2][1][2];
  L(6,8) = Lambda4[0][2][2][2];
  
  L(7,0) = Lambda4[1][2][0][0];
  L(7,1) = Lambda4[1][2][1][0];
  L(7,2) = Lambda4[1][2][2][0];
  L(7,3) = Lambda4[1][2][0][1];
  L(7,4) = Lambda4[1][2][1][1];
  L(7,5) = Lambda4[1][2][2][1];
  L(7,6) = Lambda4[1][2][0][2];
  L(7,7) = Lambda4[1][2][1][2];
  L(7,8) = Lambda4[1][2][2][2];
  
  L(8,0) = Lambda4[2][2][0][0];
  L(8,1) = Lambda4[2][2][1][0];
  L(8,2) = Lambda4[2][2][2][0];
  L(8,3) = Lambda4[2][2][0][1];
  L(8,4) = Lambda4[2][2][1][1];
  L(8,5) = Lambda4[2][2][2][1];
  L(8,6) = Lambda4[2][2][0][2];
  L(8,7) = Lambda4[2][2][1][2];
  L(8,8) = Lambda4[2][2][2][2];
  
  return; 
} // DRT::ELEMENTS::InvDesign::BuildLambdaT


/*----------------------------------------------------------------------*
 |  IF tensor                                                 (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildIF(LINALG::SerialDenseMatrix& IF, 
                                       const LINALG::SerialDenseMatrix& F) const

{
  double IF4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int l=0; l<3; ++l)
      for (int m=0; m<3; ++m)
        for (int n=0; n<3; ++n)
          IF4[k][l][m][n] = 0.5 * ( F(k,m)*F(l,n) + F(k,n)*F(l,m) );
  
  IF(0,0) = IF4[0][0][0][0];
  IF(0,1) = IF4[0][0][1][1];
  IF(0,2) = IF4[0][0][2][2];
  IF(0,3) = IF4[0][0][0][1] * 2.0;
  IF(0,4) = IF4[0][0][1][2] * 2.0;
  IF(0,5) = IF4[0][0][2][0] * 2.0;

  IF(1,0) = IF4[1][1][0][0];
  IF(1,1) = IF4[1][1][1][1];
  IF(1,2) = IF4[1][1][2][2];
  IF(1,3) = IF4[1][1][0][1] * 2.0;
  IF(1,4) = IF4[1][1][1][2] * 2.0;
  IF(1,5) = IF4[1][1][2][0] * 2.0;

  IF(2,0) = IF4[2][2][0][0];
  IF(2,1) = IF4[2][2][1][1];
  IF(2,2) = IF4[2][2][2][2];
  IF(2,3) = IF4[2][2][0][1] * 2.0;
  IF(2,4) = IF4[2][2][1][2] * 2.0;
  IF(2,5) = IF4[2][2][2][0] * 2.0;

  IF(3,0) = IF4[0][1][0][0];
  IF(3,1) = IF4[0][1][1][1];
  IF(3,2) = IF4[0][1][2][2];
  IF(3,3) = IF4[0][1][0][1] * 2.0;
  IF(3,4) = IF4[0][1][1][2] * 2.0;
  IF(3,5) = IF4[0][1][2][0] * 2.0;

  IF(4,0) = IF4[1][2][0][0];
  IF(4,1) = IF4[1][2][1][1];
  IF(4,2) = IF4[1][2][2][2];
  IF(4,3) = IF4[1][2][0][1] * 2.0;
  IF(4,4) = IF4[1][2][1][2] * 2.0;
  IF(4,5) = IF4[1][2][2][0] * 2.0;

  IF(5,0) = IF4[2][0][0][0];
  IF(5,1) = IF4[2][0][1][1];
  IF(5,2) = IF4[2][0][2][2];
  IF(5,3) = IF4[2][0][0][1] * 2.0;
  IF(5,4) = IF4[2][0][1][2] * 2.0;
  IF(5,5) = IF4[2][0][2][0] * 2.0;

  return; 
} // DRT::ELEMENTS::InvDesign::BuildIF

/*----------------------------------------------------------------------*
 |  Theta tensor                                              (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildTheta(double Theta4[][3][3][3],
                                          LINALG::SerialDenseMatrix& Theta, 
                                       const LINALG::SerialDenseMatrix& F) const

{
  //double Theta4[3][3][3][3];
  Epetra_SerialDenseMatrix K(3,3);
  K(0,0) = 1.0;
  K(1,1) = 1.0;
  K(2,2) = 1.0;
  
  for (int i=0; i<3; ++i)
    for (int j=0; j<3; ++j)
      for (int k=0; k<3; ++k)
        for (int l=0; l<3; ++l)
          // aus Micha/Christiane
          Theta4[i][j][k][l] = 0.5 * ( K(i,l)*F(k,j) + K(j,l)*F(k,i) );
  
  Theta(0,0) = Theta4[0][0][0][0]; //ok
  Theta(0,1) = Theta4[0][0][1][0];
  Theta(0,2) = Theta4[0][0][2][0];
  Theta(0,3) = Theta4[0][0][0][1];
  Theta(0,4) = Theta4[0][0][1][1];
  Theta(0,5) = Theta4[0][0][2][1];
  Theta(0,6) = Theta4[0][0][0][2];
  Theta(0,7) = Theta4[0][0][1][2];
  Theta(0,8) = Theta4[0][0][2][2];
  
  Theta(1,0) = Theta4[1][1][0][0];
  Theta(1,1) = Theta4[1][1][1][0];
  Theta(1,2) = Theta4[1][1][2][0];
  Theta(1,3) = Theta4[1][1][0][1];
  Theta(1,4) = Theta4[1][1][1][1]; //ok
  Theta(1,5) = Theta4[1][1][2][1];
  Theta(1,6) = Theta4[1][1][0][2];
  Theta(1,7) = Theta4[1][1][1][2];
  Theta(1,8) = Theta4[1][1][2][2];

  Theta(2,0) = Theta4[2][2][0][0];
  Theta(2,1) = Theta4[2][2][1][0];
  Theta(2,2) = Theta4[2][2][2][0];
  Theta(2,3) = Theta4[2][2][0][1];
  Theta(2,4) = Theta4[2][2][1][1];
  Theta(2,5) = Theta4[2][2][2][1];
  Theta(2,6) = Theta4[2][2][0][2];
  Theta(2,7) = Theta4[2][2][1][2];
  Theta(2,8) = Theta4[2][2][2][2]; //ok

  Theta(3,0) = Theta4[0][1][0][0] * 2.0; 
  Theta(3,1) = Theta4[0][1][1][0] * 2.0; //ok
  Theta(3,2) = Theta4[0][1][2][0] * 2.0; //ok
  Theta(3,3) = Theta4[0][1][0][1] * 2.0;
  Theta(3,4) = Theta4[0][1][1][1] * 2.0;
  Theta(3,5) = Theta4[0][1][2][1] * 2.0;
  Theta(3,6) = Theta4[0][1][0][2] * 2.0;
  Theta(3,7) = Theta4[0][1][1][2] * 2.0;
  Theta(3,8) = Theta4[0][1][2][2] * 2.0;

  Theta(4,0) = Theta4[1][2][0][0] * 2.0; 
  Theta(4,1) = Theta4[1][2][1][0] * 2.0;
  Theta(4,2) = Theta4[1][2][2][0] * 2.0;
  Theta(4,3) = Theta4[1][2][0][1] * 2.0; //ok
  Theta(4,4) = Theta4[1][2][1][1] * 2.0;
  Theta(4,5) = Theta4[1][2][2][1] * 2.0; //ok
  Theta(4,6) = Theta4[1][2][0][2] * 2.0;
  Theta(4,7) = Theta4[1][2][1][2] * 2.0;
  Theta(4,8) = Theta4[1][2][2][2] * 2.0;

  Theta(5,0) = Theta4[2][0][0][0] * 2.0;
  Theta(5,1) = Theta4[2][0][1][0] * 2.0;
  Theta(5,2) = Theta4[2][0][2][0] * 2.0;
  Theta(5,3) = Theta4[2][0][0][1] * 2.0;
  Theta(5,4) = Theta4[2][0][1][1] * 2.0;
  Theta(5,5) = Theta4[2][0][2][1] * 2.0;
  Theta(5,6) = Theta4[2][0][0][2] * 2.0; //ok
  Theta(5,7) = Theta4[2][0][1][2] * 2.0; //ok
  Theta(5,8) = Theta4[2][0][2][2] * 2.0;

  return; 
} // DRT::ELEMENTS::InvDesign::BuildTheta


/*----------------------------------------------------------------------*
 |  Ypsilon tensor                                            (private) |
 |                                                             gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::BuildYpsilon(double Y4[][3][3][3],
                                            LINALG::SerialDenseMatrix& Y, 
                                      const LINALG::SerialDenseMatrix& F,
                                      const LINALG::SerialDenseMatrix& S) const

{
  Epetra_SerialDenseMatrix K(3,3);
  K(0,0) = 1.0;
  K(1,1) = 1.0;
  K(2,2) = 1.0;
  
  //double Y4[3][3][3][3];
  for (int k=0; k<3; ++k)
    for (int l=0; l<3; ++l)
      for (int p=0; p<3; ++p)
        for (int q=0; q<3; ++q)
        {
          Y4[k][l][p][q] = 0.0;
          for (int m=0; m<3; ++m)
            for (int n=0; n<3; ++n)
                Y4[k][l][p][q] += 
                  0.5 * (
                          K(k,p) * ( K(m,q)*F(l,n) + K(n,q)*F(l,m) ) +
                          K(l,p) * ( K(n,q)*F(k,m) + K(m,q)*F(k,n) )
                        ) * S(m,n);
         }
                
                
                
                
  Y(0,0) = Y4[0][0][0][0];
  Y(0,1) = Y4[0][0][1][0];
  Y(0,2) = Y4[0][0][2][0];
  Y(0,3) = Y4[0][0][0][1];
  Y(0,4) = Y4[0][0][1][1];
  Y(0,5) = Y4[0][0][2][1];
  Y(0,6) = Y4[0][0][0][2];
  Y(0,7) = Y4[0][0][1][2];
  Y(0,8) = Y4[0][0][2][2];
  
  Y(1,0) = Y4[1][1][0][0];
  Y(1,1) = Y4[1][1][1][0];
  Y(1,2) = Y4[1][1][2][0];
  Y(1,3) = Y4[1][1][0][1];
  Y(1,4) = Y4[1][1][1][1];
  Y(1,5) = Y4[1][1][2][1];
  Y(1,6) = Y4[1][1][0][2];
  Y(1,7) = Y4[1][1][1][2];
  Y(1,8) = Y4[1][1][2][2];

  Y(2,0) = Y4[2][2][0][0];
  Y(2,1) = Y4[2][2][1][0];
  Y(2,2) = Y4[2][2][2][0];
  Y(2,3) = Y4[2][2][0][1];
  Y(2,4) = Y4[2][2][1][1];
  Y(2,5) = Y4[2][2][2][1];
  Y(2,6) = Y4[2][2][0][2];
  Y(2,7) = Y4[2][2][1][2];
  Y(2,8) = Y4[2][2][2][2];

  Y(3,0) = Y4[0][1][0][0];
  Y(3,1) = Y4[0][1][1][0];
  Y(3,2) = Y4[0][1][2][0];
  Y(3,3) = Y4[0][1][0][1];
  Y(3,4) = Y4[0][1][1][1];
  Y(3,5) = Y4[0][1][2][1];
  Y(3,6) = Y4[0][1][0][2];
  Y(3,7) = Y4[0][1][1][2];
  Y(3,8) = Y4[0][1][2][2];

  Y(4,0) = Y4[1][2][0][0];
  Y(4,1) = Y4[1][2][1][0];
  Y(4,2) = Y4[1][2][2][0];
  Y(4,3) = Y4[1][2][0][1];
  Y(4,4) = Y4[1][2][1][1];
  Y(4,5) = Y4[1][2][2][1];
  Y(4,6) = Y4[1][2][0][2];
  Y(4,7) = Y4[1][2][1][2];
  Y(4,8) = Y4[1][2][2][2];

  Y(5,0) = Y4[2][0][0][0];
  Y(5,1) = Y4[2][0][1][0];
  Y(5,2) = Y4[2][0][2][0];
  Y(5,3) = Y4[2][0][0][1];
  Y(5,4) = Y4[2][0][1][1];
  Y(5,5) = Y4[2][0][2][1];
  Y(5,6) = Y4[2][0][0][2];
  Y(5,7) = Y4[2][0][1][2];
  Y(5,8) = Y4[2][0][2][2];

  return; 
} // DRT::ELEMENTS::InvDesign::BuildYpsilon




#endif  // #if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#endif  // #ifdef CCADISCRET
