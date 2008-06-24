/*!----------------------------------------------------------------------
\file contchainnetw.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "contchainnetw.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_utils.H"

extern struct _MATERIAL *mat;


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         06/08|
 *----------------------------------------------------------------------*/
MAT::ContChainNetw::ContChainNetw()
  : matdata_(NULL)
{
  isinit_=false;
  mytime_=0.0;
  li_ = rcp(new vector<vector<double> >);
  li0_ = rcp(new vector<vector<double> >);
  lambda_ = rcp(new vector<vector<double> >);
  ni0_ = rcp(new vector<Epetra_SerialDenseMatrix>);
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          06/08|
 *----------------------------------------------------------------------*/
MAT::ContChainNetw::ContChainNetw(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
  int histsize;
  if (!Initialized())
  {
    histsize=0;
  }
  else 
  {
    histsize = li_->size();
  }
  AddtoPack(data,histsize);  // lenght of history vector(s)
  for (int var = 0; var < histsize; ++var) 
  {
    AddtoPack(data,li_->at(var));
    AddtoPack(data,li0_->at(var));
    AddtoPack(data,ni0_->at(var));
  }
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         06/08|
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];     // unpack pointer to my specific matdata_

  // history data
  isinit_ = true;
  int histsize;
  ExtractfromPack(position,data,histsize);
  
  if (histsize == 0) isinit_=false;
  li_ = rcp(new vector<vector<double> >);
  li0_ = rcp(new vector<vector<double> >);
  ni0_ = rcp(new vector<Epetra_SerialDenseMatrix>);
  for (int var = 0; var < histsize; ++var) {
    vector<double> li;
    vector<double> li0;
    Epetra_SerialDenseMatrix tmp(3,3);
    ExtractfromPack(position,data,li);
    ExtractfromPack(position,data,li0);
    ExtractfromPack(position,data,tmp);
    li_->push_back(li);
    li0_->push_back(li);
    ni0_->push_back(tmp);
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ContChainNetw::Initialize(const int numgp) 
{
  const double isotropy  = 1/sqrt(3.0) * matdata_->m.contchainnetw->r0;
  srand ( time(NULL) );

  li0_ = rcp(new vector<vector<double> > (numgp));
  li_ = rcp(new vector<vector<double> > (numgp));
  lambda_ = rcp(new vector<vector<double> > (numgp));
  ni0_ = rcp(new vector<Epetra_SerialDenseMatrix>);
  // initial basis is identity
  Epetra_SerialDenseMatrix id(3,3);
  for (int i=0; i<3; ++i) id(i,i) = 1.0;
  
  // initialize cell dimensions
  for(int gp=0; gp<numgp; ++gp){
    li0_->at(gp).resize(3);
    li_->at(gp).resize(3);
    lambda_->at(gp).resize(3);
//    for (int i = 0; i < 3; ++i) {
//      li0_->at(gp)[i] = ( (double)rand() / ((double)(RAND_MAX)+(double)(1)) );
//    }
    li0_->at(gp)[0] = isotropy;
    li0_->at(gp)[1] = isotropy;
    li0_->at(gp)[2] = isotropy;
    for (int i = 0; i < 3; ++i){
      li_->at(gp)[i] = li0_->at(gp)[i];
      lambda_->at(gp)[i] = 0.0;
    }
    ni0_->push_back(id);
  }

  mytime_ = 0.0;
  isinit_ = true;
  
  return ;
  
}


/*----------------------------------------------------------------------*
 |  Evaluate Material                             (public)         06/08|
 *----------------------------------------------------------------------*

*/

void MAT::ContChainNetw::Evaluate(const Epetra_SerialDenseVector* glstrain,
                                  const int gp,
                                  Teuchos::ParameterList& params,
                                  Epetra_SerialDenseMatrix* cmat,
                                  Epetra_SerialDenseVector* stress,
                                  int eleId)

{
  // bulk (isotropic) NeoHooke material parameters (Lame constants)
  const double lambda = matdata_->m.contchainnetw->lambda;
  const double mue = matdata_->m.contchainnetw->mue;
  // chain network unit cell material parameters
  const double nchain = matdata_->m.contchainnetw->nchain; // chain density ~= cell stiffness
  const double abstemp = matdata_->m.contchainnetw->abstemp; // absolute temperature (K)
  const double L = matdata_->m.contchainnetw->contl_l;  // chain contour length
  const double A = matdata_->m.contchainnetw->persl_a;  // chain persistence length
  double r0 = matdata_->m.contchainnetw->r0;      // initial chain length
  const double boltzmann = 1.3806503E-23;
  const int dim = 3;
  
  // right Cauchy-Green Tensor  C = 2 * E + I
  // build identity tensor I
  Epetra_SerialDenseVector Id(6);
  for (int i = 0; i < dim; i++) Id(i) = 1.0;
  Epetra_SerialDenseVector C(*glstrain);
  C.Scale(2.0);
  C += Id;
  
  // we need the 3 by 3 matrix as well later on -> needs improvement
  Epetra_SerialDenseMatrix CG(3,3);
  CG(0,0) = C(0); CG(1,1) = C(1); CG(2,2) = C(2);
  CG(0,1) = C(3); CG(1,0) = C(3);
  CG(1,2) = C(4); CG(2,1) = C(4);
  CG(0,2) = C(5); CG(2,0) = C(5);

  // invariants
  const double IC3 = C(0)*C(1)*C(2)
        + 0.25 * C(3)*C(4)*C(5)
        - 0.25 * C(1)*C(5)*C(5)
        - 0.25 * C(2)*C(3)*C(3)
        - 0.25 * C(0)*C(4)*C(4);    // 3rd invariant, determinant
  const double J = sqrt(IC3);
  const double lJ = log(J);
  
  // invert C
  Epetra_SerialDenseVector Cinv(6);

  Cinv(0) = C(1)*C(2) - 0.25*C(4)*C(4);
  Cinv(1) = C(0)*C(2) - 0.25*C(5)*C(5);
  Cinv(2) = C(0)*C(1) - 0.25*C(3)*C(3);
  Cinv(3) = 0.25*C(5)*C(4) - 0.5*C(3)*C(2);
  Cinv(4) = 0.25*C(3)*C(5) - 0.5*C(0)*C(4);
  Cinv(5) = 0.25*C(3)*C(4) - 0.5*C(5)*C(1);
  Cinv.Scale(1.0/IC3);

  // isotropic part: NeoHooke  ************************************************
  // W = 1/2 lambda ln^2(J) + 1/2 mue (I1-3) - mue ln(J)
  // S = (lambda ln(J) - mue) Cinv + mue Id
  // Elasticity = lambda (Cinv x Cinv) + 2(mue - lambda ln(J))(Cinv o Cinv) 
  Epetra_SerialDenseVector Siso1(Cinv);
  Siso1.Scale(lambda*lJ-mue);
//  if (gp==0){
//    cout << "before isotrop" << endl <<(*stress);
//  }
  *stress += Siso1;
  Siso1 = Id;
  Siso1.Scale(mue);
  *stress += Siso1;
//  if (gp==0){
//    cout << "after isotrop" << endl <<(*stress);
//  }
  
  AddtoCmatHolzapfelProduct((*cmat),Cinv,2*(mue-lambda*lJ));
  for (int i = 0; i < 6; ++i) {
    for (int j = 0; j < 6; ++j) {
      (*cmat)(i,j) += lambda * Cinv(i) * Cinv(j); // add lambda Cinv x Cinv
    }
  }
  
  // anisotropic part *********************************************************
  
  // the chain stiffness factor
  const double chn_stiffact = boltzmann * abstemp * nchain / (4.0*A);
  
  const double time = params.get("total time",-1.0);
  double rem_time = time;
  const double kappa = matdata_->m.contchainnetw->relax; // relaxation time for remodeling
  const double decay = exp(-kappa*rem_time);

  // initial cell dimensions (isotropy)
  vector<double> li0sq(dim);
  for (int i=0; i<dim; ++i) li0sq[i] = li0_->at(gp)[i]*li0_->at(gp)[i];
  r0 = sqrt(li0sq[0] + li0sq[1] + li0sq[2]);

  // scalar to arrive at stressfree reference conf
  const double stressfree = - chn_stiffact * ( 1.0/L + 1.0/(4.0*r0*(1.0-r0/L)*(1.0-r0/L)) - 1.0/(4.0*r0) );

  // structural tensors Ni0
  vector<LINALG::SerialDenseMatrix> Ni0(dim);
  for (int k=0; k<dim; ++k){
    LINALG::SerialDenseMatrix N0(3,3);
    for (int i=0; i<dim; ++i)
      for (int j=0; j<dim; ++j)
        N0(i,j) = (ni0_->at(gp)(i,k)) * (ni0_->at(gp)(j,k));
    Ni0.at(k) = N0;
  }
  // 'non-standard' invariants representing stretch^2 in n0_i direction
  vector<double> I(dim,0);
  for (int i=0; i<dim; ++i){
    LINALG::SerialDenseMatrix CNi0(3,3);
    CNi0.Multiply('N','N',1.0,CG,Ni0.at(i),0.0);
    I[i] = CNi0(0,0) + CNi0(1,1) + CNi0(2,2); // trace(C:Ni0)
  }
//  for(int j=0; j<8; ++j){
//    for (int i=0; i<3; ++i){
//      cout << li_->at(j)[i] << ",";
//    }
//    cout << "; gp="<<j<<endl;
//  }

  // trial S to compute eigenvalues *******************************************
  Epetra_SerialDenseMatrix trial(dim,dim);
  vector<double> lisq(dim);
  for (int i = 0; i < dim; ++i) lisq[i] = li_->at(gp)[i] * li_->at(gp)[i]; 
  
  // current chain length
  double r = sqrt(I[0]*lisq[0] + I[1]*lisq[1] + I[2]*lisq[2]);
  
  double s_chn = chn_stiffact*(4.0/L + 1.0/(r*(1.0-r/L)*(1.0-r/L)) - 1.0/r);
  
  for (int i = 0; i < dim; ++i) trial(i,i) = (*stress)(i);// + lisq[i]*s_chn + lisq[i]/I[i]*4*stressfree;
  trial(0,1) = (*stress)(3); trial(1,0) = (*stress)(3);
  trial(1,2) = (*stress)(4); trial(2,1) = (*stress)(4);
  trial(0,2) = (*stress)(5); trial(2,0) = (*stress)(5);
  
  for (int i=0; i<dim; ++i){
    Epetra_SerialDenseMatrix temp(Ni0.at(i));
    temp.Scale(lisq[i]*s_chn);
    trial += temp;
    temp = Ni0.at(i);
    temp.Scale(lisq[i]/I[i] * 4*stressfree);
    trial += temp;
  }
  
  
  Epetra_SerialDenseVector eig_sp(dim);  // lambda^(sigma+)
//  if (trial.NormOne() > 1.0E-13)
////    if (gp==0) cout << Strial;
//    //LINALG::SymmetricEigenValues(trial,eig_sp);
//    LINALG::SymmetricEigenProblem(trial,eig_sp);
//    if (gp==0){
////      cout << eig_sp;
////      cout << trial;
//    }
////  if ((params.get("total time",-1.0)>=0) && (eleId==0) && (gp==0)){
////    cout << "rem_time: " << rem_time; //params.get("total time",-1.0);
////    cout << ",eigvs: " << eig_sp(0) << "," << eig_sp(1) << "," << eig_sp(2) << endl;
////  }
//  
//  for (int i = 0; i < dim; ++i) {
////    if (eig_sort[i] > 1.0E-13) eig_sp(i) = 1.0;
////    else eig_sp(i) = 0.0;
//    lambda_->at(gp)[i] = eig_sp(i);
//    if (eig_sp(i) > 1.0E-12) eig_sp(i) = 1.0;
//    else eig_sp(i) = 0.0;
//  }
  // end of trial S to compute eigenvalues ************************************
  
  // update cell dimensions (remodeling!)
  if (rem_time > mytime_){
    mytime_ = rem_time;
    if (trial.NormOne() > 1.0E-13) LINALG::SymmetricEigenProblem(trial,eig_sp);
    for (int i = 0; i < dim; ++i) {
      lambda_->at(gp)[i] = eig_sp(i);
      if (eig_sp(i) > 1.0E-12) eig_sp(i) = 1.0;
      else eig_sp(i) = 0.0;
    }
    for (int i = 0; i < dim; ++i){
//      if (gp==0) cout << li_->at(gp)[i] << "->";
      li_->at(gp)[i] = (eig_sp(i) - li0_->at(gp)[i]/r0)*(1-decay)*r0 + li0_->at(gp)[i];
      ni0_->at(gp) = trial;
//      if (gp==0) cout << li_->at(gp)[i] << endl;
      lisq[i] = li_->at(gp)[i] * li_->at(gp)[i]; 
    }
    // structural tensors Ni0
    for (int k=0; k<dim; ++k){
      LINALG::SerialDenseMatrix N0(3,3);
      for (int i=0; i<dim; ++i)
        for (int j=0; j<dim; ++j)
          N0(i,j) = (ni0_->at(gp)(i,k)) * (ni0_->at(gp)(j,k));
      Ni0.at(k) = N0;
    }
  }
  
//  if (gp==0){
//    cout << "isotrop" << endl <<(*stress);
//    for (int i = 0; i < dim; ++i) cout << lisq[i] << "," << I[i] << ";";
//    cout << endl;
//  }
  
  // evaluate chain stress and add to isotropic stress
  r = sqrt(I[0]*lisq[0] + I[1]*lisq[1] + I[2]*lisq[2]);
  s_chn = chn_stiffact*(4.0/L + 1.0/(r*(1.0-r/L)*(1.0-r/L)) - 1.0/r);
//  for (int i = 0; i < dim; ++i){
//    (*stress)(i) += lisq[i]*s_chn + lisq[i]/I[i]*4*stressfree;
////    if (gp==0){
////      cout << "lisq=" << lisq[i] << " I=" << I[i] << "; " << lisq[i]*s_chn << " zu " << lisq[i]/I[i]*4*stressfree << endl;
////    }
//  }
  
  // stress based on remodeled parameters
  Epetra_SerialDenseMatrix rem_stress(dim,dim);
  for (int i=0; i<dim; ++i){
    Epetra_SerialDenseMatrix temp(Ni0.at(i));
    temp.Scale(lisq[i]*s_chn);
    rem_stress += temp;
    temp = Ni0.at(i);
    temp.Scale(lisq[i]/I[i] * 4*stressfree);
    rem_stress += temp;
  }
  for (int i = 0; i < dim; ++i) (*stress)(i) += rem_stress(i,i);
  (*stress)(3) += rem_stress(0,1);
  (*stress)(4) += rem_stress(1,2);
  (*stress)(5) += rem_stress(0,2);
  
  if (gp==0){
//    cout << "anisotrop" << endl << (*glstrain) <<(*stress) << endl;
  }
  
  //evaluate chain tangent and add to isotropic cmat
//  if (gp==0){
//    cout << "isotrop" << endl <<(*cmat);
//  }
  double c_chn = chn_stiffact/(r*r*r) * (1.0 - 1.0/((1.0-r/L)*(1.0-r/L)) + 2.0*r/(L*(1.0-r/L)*(1.0-r/L)*(1.0-r/L)) );
//  // chain part only affects upper left block of cmat (no shear)
//  for (int i = 0; i < dim; ++i) {
//    for (int j = 0; j < dim; ++j) (*cmat)(i,j) += c_chn * lisq[i] * lisq[j];
//    (*cmat)(i,i) +=  - 8.0*stressfree*lisq[i]/(I[i]*I[i]);
//  }
  
  Epetra_SerialDenseMatrix sumN0i(3,3);
  for (int k=0; k<dim; ++k)
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < dim; ++j) sumN0i(i,j) += lisq[i] * (Ni0.at(k))(i,j);
  
  ElastSymTensorMultiply((*cmat),c_chn,sumN0i,sumN0i,1.0);  // fac sumN0i x sumN0i
  
  for (int i=0; i<dim; ++i)
    ElastSymTensorMultiply((*cmat),-8.0*stressfree*lisq[i]/(I[i]*I[i]),Ni0[i],Ni0[i],1.0);
  
  
//  if (gp==0){
//    cout << "anisotrop" << endl <<(*cmat) << endl;
//  }
  
  return;
}


void MAT::ChainOutputToGmsh(const Teuchos::RCP<DRT::Discretization> dis,
                                      const double time,
                                      const int iter)
{
  std::stringstream filename;
  filename << allfiles.outputfile_kenner << "_ContChainMat" << std::setw(3) << setfill('0') << time << std::setw(2) << setfill('0') << iter << ".pos";
  std::ofstream f_system(filename.str().c_str());

  stringstream gmshfilecontent;
  gmshfilecontent << "View \" Time: " << time << " Iter: " << iter << " \" {" << endl;
  for (int iele=0; iele<dis->NumMyColElements(); ++iele)
  {
    const DRT::Element* actele = dis->lColElement(iele);
    
    // build current configuration
    vector<int> lm;
    vector<int> lmowner;
    actele->LocationVector(*dis,lm,lmowner);
    RCP<const Epetra_Vector> disp = dis->GetState("displacement");
    vector<double> mydisp(lm.size(),0);
    //DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
    const int numnode = actele->NumNode();
    const int numdof = 3;
    blitz::Array<double,2> xyze(3, numnode);
    for (int inode = 0; inode < numnode; ++inode)
    {
      xyze(0, inode) = actele->Nodes()[inode]->X()[0]+ mydisp[inode*numdof+0];
      xyze(1, inode) = actele->Nodes()[inode]->X()[1]+ mydisp[inode*numdof+1];
      xyze(2, inode) = actele->Nodes()[inode]->X()[2]+ mydisp[inode*numdof+2];
    }
    gmshfilecontent << IO::GMSH::cellWithScalarToString(actele->Shape(),
        1.0, xyze) << endl;
    
    vector<double> elecenter = MAT::MatPointCoords(actele,mydisp);
    RefCountPtr<MAT::Material> mat = actele->Material();
    MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
    Epetra_SerialDenseMatrix ni0 = chain->Getni0()->at(0);
    vector<double> lamb0 = chain->Getlambdas()->at(0);
    
//    // material plot at element center
//    const int dim=3;
//    for (int k=0; k<dim; ++k){
//      gmshfilecontent << "VP(" << scientific << elecenter[0] << ",";
//      gmshfilecontent << scientific << elecenter[1] << ",";
//      gmshfilecontent << scientific << elecenter[2] << ")";
//      gmshfilecontent << "{" << scientific << 
//      ni0(0,k) * lamb0[k]
//      << "," << ni0(1,k) * lamb0[k] << "," << ni0(2,k) * lamb0[k] << "};" << endl;
//    }
    
    // material plot at gauss points
    int ngp = chain->Getni0()->size();
    for (int gp = 0; gp < ngp; ++gp){
      vector<double> point = MAT::MatPointCoords(actele,mydisp,gp);
      double scalar = 1;
      vector<double> length(3);
      vector<double> gpli =  chain->Getli()->at(gp); 
      vector<double> gpli0 = chain->Getli0()->at(gp);
//      length[0] = scalar * gpli[0]/gpli0[0];
//      length[1] = scalar * gpli[1]/gpli0[1];
//      length[2] = scalar * gpli[2]/gpli0[2];
      length[0] = scalar * chain->Getlambdas()->at(gp)[0];
      length[1] = scalar * chain->Getlambdas()->at(gp)[1];
      length[2] = scalar * chain->Getlambdas()->at(gp)[2];
//      for (int i=0; i<3; ++i){
//        cout << gpli[i] << ":" << gpli0[i] << endl;
//      }
      
      for (int k=0; k<3; ++k){
        // draw eigenvectors
        gmshfilecontent << "VP(" << scientific << point[0] << ",";
        gmshfilecontent << scientific << point[1] << ",";
        gmshfilecontent << scientific << point[2] << ")";
        gmshfilecontent << "{" << scientific
        << ((chain->Getni0())->at(gp))(0,k)
        << "," << ((chain->Getni0())->at(gp))(1,k)
        << "," << ((chain->Getni0())->at(gp))(2,k) << "};" << endl;
        
//        // draw fiber cell vectors
//        gmshfilecontent << "VP(" << scientific << point[0] << ",";
//        gmshfilecontent << scientific << point[1] << ",";
//        gmshfilecontent << scientific << point[2] << ")";
//        gmshfilecontent << "{" << scientific
//        <<        ((chain->Getni0())->at(gp))(0,k) * length[k] 
//        << "," << ((chain->Getni0())->at(gp))(1,k) * length[k] 
//        << "," << ((chain->Getni0())->at(gp))(2,k) * length[k] 
//        << "};" << endl;
      }
    }
  }
  gmshfilecontent << "};" << endl;

  f_system << gmshfilecontent.str();
  f_system.close();
  
  return;
}

const vector<double> MAT::MatPointCoords(const DRT::Element* actele,const vector<double>& mydisp, int gp)
{
  // update element geometry
  const int numnode = actele->NumNode();
  const int numdof = 3;
  Epetra_SerialDenseMatrix xrefe(numnode,3);  // material coord. of element
  for (int i=0; i<numnode; ++i){
    xrefe(i,0) = actele->Nodes()[i]->X()[0]+ mydisp[i*numdof+0];
    xrefe(i,1) = actele->Nodes()[i]->X()[1]+ mydisp[i*numdof+1];
    xrefe(i,2) = actele->Nodes()[i]->X()[2]+ mydisp[i*numdof+2];
  }
  const DRT::Element::DiscretizationType distype = actele->Shape();
  Epetra_SerialDenseVector funct(numnode);
  // Element midpoint at r=s=t=0.0
  if (gp==-1) DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  else{
    const DRT::UTILS::GaussRule3D gaussrule_ = DRT::UTILS::intrule_hex_8point;
    const DRT::UTILS::IntegrationPoints3D intpoints = getIntegrationPoints3D(gaussrule_);
    DRT::UTILS::shape_function_3D(funct, intpoints.qxg[gp][0], intpoints.qxg[gp][1], intpoints.qxg[gp][2], distype);
  }
  Epetra_SerialDenseMatrix point(1,3);
  point.Multiply('T','N',1.0,funct,xrefe,0.0);
  vector<double> coords(3);
  coords[0] = point(0,0);
  coords[1] = point(0,1);
  coords[2] = point(0,2);
  return coords;
}


#endif

