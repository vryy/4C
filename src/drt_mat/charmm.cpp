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
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
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

    // The following needs to come from the parameter in the final version
    vector<string> strain_type;
    strain_type.push_back ("principal");
    strain_type.push_back ("vector");

    vector<double> his_charmm;
    his_charmm.push_back(0.0); // actual time
    his_charmm.push_back(0.0); // time at last timestep
    for(int i=0;i<(int)strain_type.size();i++) {
        his_charmm.push_back(1.0);  // updated lambda(1)(t)
        his_charmm.push_back(1.0);  // lambda(1)(t-dt)
        his_charmm.push_back(1.0);  // updated lambda(2)(t)
        his_charmm.push_back(1.0);  // lambda(2)(t-dt)
        his_charmm.push_back(1.0);  // updated lambda(3)(t)
        his_charmm.push_back(1.0);  // lambda(3)(t-dt)
        his_charmm.push_back(3.0);  // updated I1(t)
        his_charmm.push_back(3.0);  // I1(t-dt)
    }
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
    #ifdef DEBUG
    if (!glstrain || !cmat || !stress)
    dserror("Data missing upon input in material neo hooke");
    #endif

    // Parameter collection
    // evaluate lamda at origin or at gp
    bool origin = false;  // change only of xref and xcurr really working!!!!
    // length of the protein in the main pulling direction [A]
    double characteristic_length = 40.625; //50; //originally 44
    // characteristic direction of the protein
    // Possible selcetions:
    // principal = main strain direction (biggest eigenvalue)
    // vector = using the given vector
    // none = don't use the direction
    vector<string> strain_type;
    strain_type.push_back ("principal");
    strain_type.push_back ("vector");
    vector<LINALG::SerialDenseVector> d;
    LINALG::SerialDenseVector d_1(3);
    LINALG::SerialDenseVector d_2(3);
    d_1(0) = 0;
    d_1(1) = 1;
    d_1(2) = 0;
    d_2(0) = 1;
    d_2(1) = 0;
    d_2(2) = 0;
    d.push_back(d_1);
    d.push_back(d_2);
    // Use the hard coded charmm results (charmmfakeapi == true) or call charmm really (charmmfakeapi == false)
    bool charmmhard = true;

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
        LINALG::Matrix<3,3> V(C);
        LINALG::SerialDenseVector lambda(3);
        vector<LINALG::SerialDenseVector> dir_lambdas;
        vector<LINALG::Matrix<3,3> > dir_eigenv;
        // go through number of directions
        for(int i=0;i<(int)strain_type.size();i++) {
            if ( strain_type[i].compare("principal") == 0 ) {
                V.SetCopy(C);
                EvalStrain(origin,xrefe,xcurr,V,lambda);
                dir_lambdas.push_back(lambda);
                dir_eigenv.push_back(V);
            } else if ( strain_type[i].compare("vector") == 0 ) {
                V.SetCopy(C);
                for (int k=0;k<3;k++) for (int l=0;l<3;l++) V(k,l) = d[i](k) * V(k,l) * d[i](l);
                EvalStrain(origin,xrefe,xcurr,V,lambda);
                dir_lambdas.push_back(lambda);
                dir_eigenv.push_back(V);
            } else if ( strain_type[i].compare("none") == 0 ) {
                V.Clear();
                lambda.Zero();
                dir_lambdas.push_back(lambda);
                dir_eigenv.push_back(V);
            } else {
                dserror("No valid strain type given for CHARmm!");
            }
        }
        //cout << dir_lambdas[1](0) << " : " << dir_lambdas[1](1) << " : " << dir_lambdas[1](2) << endl;
        //cout << dir_eigenv[1] << endl;


        // Update and reconfigure history
        vector<double>* his;
        his = data_.GetMutable<vector<double> >("his_charmm");
        if ( (*his)[0] < time ) {
            (*his)[1] = (*his)[0]; // time
            for(int i=0;i<(int)strain_type.size();i++) {
                (*his)[3+(i*8)] = (*his)[2+(i*8)]; // lambda(0)
                (*his)[5+(i*8)] = (*his)[4+(i*8)]; // lambda(1)
                (*his)[7+(i*8)] = (*his)[6+(i*8)]; // lambda(2)
                (*his)[9+(i*8)] = (*his)[8+(i*8)]; // I1
            }
        }
        (*his)[0] = time;
        for(int i=0;i<(int)strain_type.size();i++) {
            (*his)[2+(i*8)] = dir_lambdas[i](0);
            (*his)[4+(i*8)] = dir_lambdas[i](1);
            (*his)[6+(i*8)] = dir_lambdas[i](2);
            (*his)[8+(i*8)] = I1;
        }
        //cout  << (*his)[0] << " : " << (*his)[1] << " : ";
        //for(int i=0;i<(int)strain_type.size();i++) {
        //    cout    << (*his)[2+(i*8)] << " : " << (*his)[3+(i*8)] << " : "
        //            << (*his)[4+(i*8)] << " : " << (*his)[5+(i*8)] << " : "
        //            << (*his)[6+(i*8)] << " : " << (*his)[7+(i*8)] << " : "
        //            << (*his)[8+(i*8)] << " : " << (*his)[9+(i*8)] << " : ";
        //}
        //cout <<  endl;

        // Prepare and call CHARmm in its beauty itself
        // get lambda t-dt information
        vector<double> lambda_his;
        for(int i=0;i<(int)strain_type.size();i++) {
            lambda_his.push_back((*his)[7+(i*8)]);
        }

        // calculate STARTD and ENDD for CHARmm
        double STARTD = characteristic_length * (1 - lambda_his[0]);
        double ENDD = characteristic_length * (1 - dir_lambdas[0](2)); // Check for better way to choose!!!!
        //cout << "STARTD: " << STARTD << " ENDD: " << ENDD << endl;
        
        // Check if results actually can be computed by CHARmm
        //if (STARTD != ENDD) dserror("STARTD and ENDD identical! CHARmm will not produce any results.");

        // Call API to CHARMM
        // Results vector: charmm_result
        // (Energy STARTD, Energy ENDD, #Atoms STARTD, #Atoms ENDD, Volume STARTD, Volume ENDD)
        LINALG::SerialDenseVector direction(3);
        LINALG::SerialDenseVector charmm_result(6);
        if (charmmhard) {
            charmmfakeapi(STARTD,ENDD,direction,charmm_result);
        } else {
            //if (STARTD != ENDD) charmmfileapi(STARTD,ENDD,direction,charmm_result);
        }


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
                        const LINALG::SerialDenseMatrix& xrefe,
                        const LINALG::SerialDenseMatrix& xcurr,
                        LINALG::Matrix<3,3>& C,
                        LINALG::SerialDenseVector& lambda )
{
    LINALG::SerialDenseVector lambda2(3);
    LINALG::SerialDenseMatrix Ctmp(3,3);
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
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) Ctmp(i,j) = C0(i,j);
        LINALG::SymmetricEigen(Ctmp,lambda2,'V',false);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) C(i,j) = Ctmp(i,j);
    } else {
        // compute current eigenvalues of gaussian point C
        LINALG::SerialDenseMatrix Ctmp(3,3);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) Ctmp(i,j) = C(i,j);
        LINALG::SymmetricEigen(Ctmp,lambda2,'V',false);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) C(i,j) = Ctmp(i,j);
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


/*----------------------------------------------------------------------*/
//! Hard coupling without calling CHARMm
/*----------------------------------------------------------------------*/
void MAT::CHARMM::charmmfakeapi ( const double STARTD,
                                    const double ENDD,
                                    const LINALG::SerialDenseVector direction,
                                    LINALG::SerialDenseVector& charmm_result)
{
    // Define the number n of steps / results from CHARmm (or any MD simluation)
    // If n=2 then it is assumes that always the same values should be used for all steps.
    const int n = 2;
    // Define roundoff for choosing in which step we are
    const double roundoff = 0.005;
    // Hard coded data structure (second variable):
    // (STARTD, Energy, # of Atoms, Volume)
    LINALG::SerialDenseMatrix MD(n,4);

    // Hard coded results from MD
    MD(0,0) = 0.0;
    MD(0,1) = -330.912;
    MD(0,2) = 1202;
    MD(0,3) = 9954.29;

    MD(1,0) = -0.8125;
    MD(1,1) = -321.671;
    MD(1,2) = 1141;
    MD(1,3) = 9441.08;

    // Compute the charmm_result vector
    // (Energy STARTD, Energy ENDD, #Atoms STARTD, #Atoms ENDD, Volume STARTD, Volume ENDD)
    ios_base::fmtflags flags = cout.flags( ); // Save original flags

    for (int i=n-1; i>=0; i--) {
        //cout << ENDD << " " << MD(i,0);
        if (abs(ENDD) == 0.0) { // start call at the beginning; just to give some information
            i=0;
            cout << setw(4) << left << "MD (" << showpoint << STARTD << setw(2)  << "->" << ENDD << setw(3) << "): " << flush;
            charmm_result[0] = NAN;
            charmm_result[1] = MD(i,1);
            charmm_result[2] = NAN;
            charmm_result[3] = MD(i,2);
            charmm_result[4] = NAN;
            charmm_result[5] = MD(i,3);
            cout << setw(4) << "V(0):" << setw(15) << left << scientific << setprecision(6) << (charmm_result[1]);
            cout << setw(8) << "#Atoms:" << setw(10) << left << fixed << setprecision(0) << charmm_result[3] << setw(8) << "Volume:" << setw(12) << left << setprecision(2) << charmm_result[5] << endl;
            i=-1; //break loop
        } else if (abs(ENDD) < (abs(MD(i,0)) + roundoff) && abs(ENDD) > (abs(MD(i,0)) - roundoff)) {
            // main loop where basically at every step the data is given
            cout << setw(4) << left << "MD (" << showpoint << STARTD << setw(2)  << "->" << ENDD << setw(3) << "): " << flush;
            charmm_result[0] = MD(i-1,1);
            charmm_result[1] = MD(i,1);
            charmm_result[2] = MD(i-1,2);
            charmm_result[3] = MD(i,2);
            charmm_result[4] = MD(i-1,3);
            charmm_result[5] = MD(i,3);
            cout << setw(4) << "dV:" << setw(15) << left << scientific << setprecision(6) << (charmm_result[1] - charmm_result[0]);
            cout << setw(8) << "#Atoms:" << setw(10) << left << fixed << setprecision(0) << charmm_result[3] << setw(8) << "Volume:" << setw(12) << left << setprecision(2) << charmm_result[5] << endl;
            i=-1; //break loop
        } else {
            // in case that only one dV is given, use it for all. If more then break.
            if (n==2) {
                cout << setw(4) << left << "MD (" << showpoint << STARTD << setw(2)  << "->" << ENDD << setw(3) << "): " << flush;
                charmm_result[0] = MD(0,1);
                charmm_result[1] = MD(1,1);
                charmm_result[2] = MD(0,2);
                charmm_result[3] = MD(1,2);
                charmm_result[4] = MD(0,3);
                charmm_result[5] = MD(1,3);
                cout << setw(4) << "dV:" << setw(15) << left << scientific << setprecision(6) << (charmm_result[1] - charmm_result[0]);
                cout << setw(8) << "#Atoms:" << setw(10) << left << fixed << setprecision(0) << charmm_result[3] << setw(8) << "Volume:" << setw(12) << left << setprecision(2) << charmm_result[5] << endl;
                i=-1; //break loop
            } else {
                dserror("No appropriate MD result found for ENDD");
            }
        }
    }
    cout.flags(flags);
}


#endif

