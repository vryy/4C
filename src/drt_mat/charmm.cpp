/*!----------------------------------------------------------------------
\file charmm.cpp
\brief

<pre>
Maintainer: Robert Metzke
            metzke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "sys/types.h"
#include "sys/stat.h"
#include <math.h>
#include "charmm.H"
#include "../drt_so3/so_hex8.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::CHARMM::CHARMM( Teuchos::RCP<MAT::PAR::Material> matdata )
: Parameter(matdata),
  density_(matdata->GetDouble("DENS"))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::CHARMM::CHARMM()
  : params_(NULL)
{
}

MAT::CHARMM::CHARMM(MAT::PAR::CHARMM* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::CHARMM::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = params_->Id();
  AddtoPack(data,matid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::CHARMM::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::CHARMM*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


/*----------------------------------------------------------------------*/
//! Setup CHARmm history variables
// Actual and history variables, which need to be stored
// The updated version will be written in every iteration step.
// his_charmm = (    updated time, lasttime,
//                updated lambda(1), lambda(1)(t-dt),
//                updated lambda(2), lambda(2)(t-dt),
//                updated lambda(3), lambda(3)(t-dt),
//                I1,  I1(t-dt))
// his_mat[0] = c1 Neohookean from CHARMM for complete element
/*---------------------------------------------------------------------*/
void MAT::CHARMM::Setup(DRT::Container& data_)
{

    vector<double> his_charmm(10);
    his_charmm[0] = 0.0; // actual time
    his_charmm[1] = 0.0; // time at last timestep
    his_charmm[2] = 1.0;  // updated lambda(1)(t)
    his_charmm[3] = 1.0;  // lambda(1)(t-dt)
    his_charmm[4] = 1.0;  // updated lambda(2)(t)
    his_charmm[5] = 1.0;  // lambda(2)(t-dt)
    his_charmm[6] = 1.0;  // updated lambda(3)(t)
    his_charmm[7] = 1.0;  // lambda(3)(t-dt)
    his_charmm[8] = 3.0;  // updated I1(t)
    his_charmm[9] = 3.0;  // I1(t-dt)
    data_.Add("his_charmm",his_charmm);

    vector<double> his_mat(1);
    data_.Add("his_mat",his_mat); // material property from CHARmm

    return;
}

/*----------------------------------------------------------------------*/
//! Evaluate second PK and constitutive tensor
/*----------------------------------------------------------------------*/
void MAT::CHARMM::Evaluate( const LINALG::Matrix<NUM_STRESS_3D,1>* glstrain,
                            LINALG::Matrix<NUM_STRESS_3D,NUM_STRESS_3D>* cmat,
                            LINALG::Matrix<NUM_STRESS_3D,1>* stress,
                            const int ele_ID,
                            const int gp,
                            DRT::Container& data_,
                            const double time,
                            const LINALG::SerialDenseMatrix& xrefe,
                            const LINALG::SerialDenseMatrix& xcurr)
{

    // Parameter collection
    // evaluate lamda at origin or at gp
    bool origin = false;  // change only of xref and xcurr really working!!!!
    // length of the protein in the main pulling direction [A]
    double characteristic_length = 40.625; //50; //originally 44


    // Identity Matrix
    LINALG::Matrix<3,3> I(true);
    for (int i = 0; i < 3; ++i) I(i,i) = 1.0;

    // Green-Lagrange Strain Tensor
    LINALG::Matrix<3,3> E(false);
    E(0,0) = (*glstrain)(0);
    E(1,1) = (*glstrain)(1);
    E(2,2) = (*glstrain)(2);
    E(0,1) = 0.5 * (*glstrain)(3);  E(1,0) = 0.5 * (*glstrain)(3);
    E(1,2) = 0.5 * (*glstrain)(4);  E(2,1) = 0.5 * (*glstrain)(4);
    E(0,2) = 0.5 * (*glstrain)(5);  E(2,0) = 0.5 * (*glstrain)(5);

    // Right Cauchy-Green Tensor  C = 2 * E + I
    LINALG::Matrix<3,3> C(E);
    C.Scale(2.0);
    C += I;

    // Principal Invariants I1 = tr(C) and I3 = det(C)
    double I1 = C(0,0)+C(1,1)+C(2,2);
    double I3 = C(0,0)*C(1,1)*C(2,2) + C(0,1)*C(1,2)*C(2,0)
              + C(0,2)*C(1,0)*C(2,1) - (C(0,2)*C(1,1)*C(2,0)
              + C(0,1)*C(1,0)*C(2,2) + C(0,0)*C(1,2)*C(2,1));

    // Calculation of C^-1 (Cinv)
    LINALG::Matrix<3,3> Cinv(C);
    Cinv.Invert();

    /////////////////////////////////////////////////////////////////////CHARmm
    // CHARmm things come here
    if (gp == 0) {
        // Get the strains in the characteristic directions
        LINALG::SerialDenseVector lambda(3);
        EvalStrain(origin,C,xrefe,xcurr,lambda);
        //cout << C << lambda(0) << " : " << lambda(1) << " : " << lambda(2) << endl;

        // Update and reconfigure history
        vector<double>* his;
        his = data_.GetMutable<vector<double> >("his_charmm");
        if ( (*his)[0] < time ) {
            (*his)[1] = (*his)[0]; // time
            (*his)[3] = (*his)[2]; // lambda(0)
            (*his)[5] = (*his)[4]; // lambda(1)
            (*his)[7] = (*his)[6]; // lambda(2)
            (*his)[9] = (*his)[8]; // I1
        }
        (*his)[0] = time;
        (*his)[2] = lambda(0);
        (*his)[4] = lambda(1);
        (*his)[6] = lambda(2);
        (*his)[8] = I1;
        //cout  << (*his)[0] << " : " << (*his)[1] << " : "
        //      << (*his)[2] << " : " << (*his)[3] << " : "
        //      << (*his)[4] << " : " << (*his)[5] << " : "
        //      << (*his)[6] << " : " << (*his)[7] << " : "
        //      << (*his)[8] << " : " << (*his)[9] << " : "
        //      <<  endl;

        // Prepare and call CHARmm in its beauty itself
        // get lambda t-dt information
        double lambda_his = (*his)[7];

        // calculate STARTD and ENDD for CHARMM
        double STARTD = characteristic_length * (1 - lambda_his);
        double ENDD = characteristic_length * (1 - lambda(2)); // Check for better way to choose!!!!

        // Call API to CHARMM
        // Results vector: charmm_result
        // (Energy STARTD, Energy ENDD, #Atoms STARTD, #Atoms ENDD, Volume STARTD, Volume ENDD)
        LINALG::SerialDenseVector direction(3);
        LINALG::SerialDenseVector charmm_result(6);
        //if (STARTD != ENDD) charmmfileapi(STARTD,ENDD,direction,charmm_result);


    }

    //
    ///////////////////////////////////////////////////////////////////////////

    // Material Constants c1 and beta
    double ym = 1000; // intermediate for testing purpose only
    double nu = 0.3; // intermediate for testing purpose only
    double c1 = 0.5 * ym/(2*(1+nu)); // intermediate for testing purpose only
    double beta = nu/(1-2*nu);
    //if (time > 0.0) {
    //    vector<double>* c_c1;
    //    c_c1 = data_.GetMutable<vector<double> >("c_c1");
    //    c1 = (*c_c1)[0];
    //} else {
    //    c1 = 0.0;
    //}

    // Energy
    //double W = c1/beta * (pow(I3,-beta) - 1) + c1 * (I1-3);

    // PK2 Stresses
    LINALG::Matrix<3,3> PK2(false);
    int i,j;
    for (i=0; i<3; i++)
    for (j=0; j<3; j++)
    {
        PK2(i,j)= 2 * c1 * ( I(i,j) - pow(I3,-beta) * Cinv(i,j) );
    }

    // Transfer PK2 tensor to stress vector
    (*stress)(0) = PK2(0,0);
    (*stress)(1) = PK2(1,1);
    (*stress)(2) = PK2(2,2);
    (*stress)(3) = PK2(0,1);
    (*stress)(4) = PK2(1,2);
    (*stress)(5) = PK2(0,2);

    // Elasticity Tensor
    double delta6 = 4. * c1 * beta * pow(I3,-beta);
    double delta7 = 4. * c1 * pow(I3,-beta);

    int k,l;
    LINALG::Matrix<9,9> ET(false);


    for (k=0; k<3; k++)
    for (l=0; l<3; l++)
    {
        ET(k,l) =       delta6 * (Cinv(0,0) * Cinv(k,l)) +
                        delta7 * 0.5 *(Cinv(0,k) * Cinv(0,l) + Cinv(0,l) * Cinv(0,k));
        ET(k+3,l) =     delta6 * (Cinv(1,0) * Cinv(k,l)) +
                        delta7 * 0.5 *(Cinv(1,k) * Cinv(0,l) + Cinv(1,l) * Cinv(0,k));
        ET(k+3,l+3) =   delta6 * (Cinv(1,1) * Cinv(k,l)) +
			delta7 * 0.5 *(Cinv(1,k) * Cinv(1,l) + Cinv(1,l) * Cinv(1,k));
        ET(k+6,l) =     delta6 * (Cinv(2,0) * Cinv(k,l)) +
			delta7 * 0.5 *(Cinv(2,k) * Cinv(0,l) + Cinv(2,l) * Cinv(0,k));
        ET(k+6,l+3) =   delta6 * (Cinv(2,1) * Cinv(k,l)) +
			delta7 * 0.5 *(Cinv(2,k) * Cinv(1,l) + Cinv(2,l) * Cinv(1,k));
        ET(k+6,l+6) =   delta6 * (Cinv(2,2) * Cinv(k,l)) +
			delta7 * 0.5 *(Cinv(2,k) * Cinv(2,l) + Cinv(2,l) * Cinv(2,k));
    }
    
    (*cmat)(0,0)=ET(0,0);
    (*cmat)(0,1)=ET(1,1);
    (*cmat)(0,2)=ET(2,2);
    (*cmat)(0,3)=ET(1,0);
    (*cmat)(0,4)=ET(2,1);
    (*cmat)(0,5)=ET(2,0);

    (*cmat)(1,0)=ET(3,3);
    (*cmat)(1,1)=ET(4,4);
    (*cmat)(1,2)=ET(5,5);
    (*cmat)(1,3)=ET(4,3);
    (*cmat)(1,4)=ET(5,4);
    (*cmat)(1,5)=ET(5,3);

    (*cmat)(2,0)=ET(6,6);
    (*cmat)(2,1)=ET(7,7);
    (*cmat)(2,2)=ET(8,8);
    (*cmat)(2,3)=ET(7,6);
    (*cmat)(2,4)=ET(8,7);
    (*cmat)(2,5)=ET(8,6);

    (*cmat)(3,0)=ET(3,0);
    (*cmat)(3,1)=ET(4,1);
    (*cmat)(3,2)=ET(5,2);
    (*cmat)(3,3)=ET(4,0);
    (*cmat)(3,4)=ET(5,1);
    (*cmat)(3,5)=ET(5,0);

    (*cmat)(4,0)=ET(6,3);
    (*cmat)(4,1)=ET(7,4);
    (*cmat)(4,2)=ET(8,5);
    (*cmat)(4,3)=ET(7,3);
    (*cmat)(4,4)=ET(8,4);
    (*cmat)(4,5)=ET(8,3);

    (*cmat)(5,0)=ET(6,0);
    (*cmat)(5,1)=ET(7,1);
    (*cmat)(5,2)=ET(8,2);
    (*cmat)(5,3)=ET(7,0);
    (*cmat)(5,4)=ET(8,1);
    (*cmat)(5,5)=ET(8,0);

    return;
}


/*----------------------------------------------------------------------*/
//! Evaluate strains in the charateristic directions
/*----------------------------------------------------------------------*/
void MAT::CHARMM::EvalStrain( const bool& origin,
                        const LINALG::Matrix<3,3>& C,
                        const LINALG::SerialDenseMatrix& xrefe,
                        const LINALG::SerialDenseMatrix& xcurr,
                        LINALG::SerialDenseVector& lambda )
{
    LINALG::SerialDenseVector lambda2(3);
    if (origin) {
          // vector of dN/dxsi |r=s=t=0.0
          double dN0_vector[24] =
               {-0.125,-0.125,-0.125,
                +0.125,-0.125,-0.125,
                +0.125,+0.125,-0.125,
                -0.125,+0.125,-0.125,
                -0.125,-0.125,+0.125,
                +0.125,-0.125,+0.125,
                +0.125,+0.125,+0.125,
                -0.125,+0.125,+0.125};

          // shape function derivatives, evaluated at origin (r=s=t=0.0)
          Epetra_DataAccess CV = Copy;
          Epetra_SerialDenseMatrix dN0(CV,dN0_vector,3,3,8);

          // compute Jacobian, evaluated at element origin (r=s=t=0.0)
          LINALG::SerialDenseMatrix invJacobian0(3,3);
          invJacobian0.Multiply('N','N',1.0,dN0,xrefe,0.0);
          const double detJacobian0 = LINALG::NonsymInverse3x3(invJacobian0);
          if (detJacobian0 < 0.0) dserror("Jacobian at origin negativ (CHARMMAPI)");

          //cout << invJacobian0 << endl;
          LINALG::SerialDenseMatrix N_XYZ(3,8);
          //compute derivatives N_XYZ at gp w.r.t. material coordinates
          // by N_XYZ = J^-1 * N_rst
          N_XYZ.Multiply('N','N',1.0,invJacobian0,dN0,0.0);
          // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
          LINALG::SerialDenseMatrix defgrd0(3,3);
          defgrd0.Multiply('T','T',1.0,xcurr,N_XYZ,0.0);
          // Right Cauchy-Green tensor = F^T * F
          LINALG::SerialDenseMatrix C0(3,3);
          C0.Multiply('T','N',1.0,defgrd0,defgrd0,0.0);

          // compute current eigenvalues of gaussian point C
          LINALG::SerialDenseMatrix Ctmp0(C0);
          LINALG::SymmetricEigen(Ctmp0,lambda2,'V',false);
    } else {
        // compute current eigenvalues of gaussian point C
        LINALG::SerialDenseMatrix Ctmp(3,3);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) Ctmp(i,j) = C(i,j);
        LINALG::SymmetricEigen(Ctmp,lambda2,'V',false);
    }
    for (int i=0;i<3;i++) lambda(i) = sqrt(lambda2(i));

}


/*----------------------------------------------------------------------*
 |  File based API to CHARMM                                    rm 03/08|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector MAT::CHARMM::charmmfileapi ( const double STARTD,
                                    const double ENDD,
                                    const LINALG::SerialDenseVector direction,
                                    LINALG::SerialDenseVector& charmm_result)
{
	cout << "Too be implemented....." << endl;
	exit(1);
	return(0);
}


#endif

