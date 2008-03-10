/*!----------------------------------------------------------------------
\file so_sh8_evaluate.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "so_sh8.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../io/gmsh.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/material.H"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra


/*----------------------------------------------------------------------*
 |  find shell-thickness direction via Jacobian                maf 07/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh8::ThicknessDirection DRT::ELEMENTS::So_sh8::sosh8_findthickdir() 
{
  // update element geometry
  Epetra_SerialDenseMatrix xrefe(NUMNOD_SOH8,NUMDIM_SOH8); // material coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i) {
    xrefe(i, 0) = this->Nodes()[i]->X()[0];
    xrefe(i, 1) = this->Nodes()[i]->X()[1];
    xrefe(i, 2) = this->Nodes()[i]->X()[2];
  }
  // vector of df(origin)
  double df0_vector[NUMDOF_SOH8*NUMNOD_SOH8] =
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
  Epetra_SerialDenseMatrix df0(CV, df0_vector, NUMDIM_SOH8, NUMDIM_SOH8, 
  NUMNOD_SOH8);

  // compute Jacobian, evaluated at element origin (r=s=t=0.0)
  Epetra_SerialDenseMatrix jac0(NUMDIM_SOH8,NUMDIM_SOH8);
  jac0.Multiply('N', 'N', 1.0, df0, xrefe, 0.0);
  // compute inverse of Jacobian at element origin
  Epetra_SerialDenseSolver solve_for_inverseJ0;
  Epetra_SerialDenseMatrix iJ0(jac0);
  solve_for_inverseJ0.SetMatrix(iJ0);
  int err = solve_for_inverseJ0.Invert();
  if (err != 0) dserror("Inversion of Jacobian0 failed");
  
  // separate "stretch"-part of J-mapping between parameter and global space
  Epetra_SerialDenseMatrix jac0stretch(3,3); 
  jac0stretch.Multiply('T','N',1.0,iJ0,iJ0,0.0); // jac0stretch = J^{-T}J
  double r_stretch = sqrt(jac0stretch(0,0));
  double s_stretch = sqrt(jac0stretch(1,1));
  double t_stretch = sqrt(jac0stretch(2,2));
  
  // minimal stretch equivalents with "thinnest" direction
  double min_stretch = min(r_stretch, min(s_stretch, t_stretch));
  
  ThicknessDirection thickdir; // of actual element
  int thick_index;
  
  if (min_stretch == r_stretch) {
    if ((min_stretch / s_stretch <= 0.5) || (min_stretch / t_stretch <=0.5)) {
      dserror("Solid-Shell element geometry has not a shell aspect ratio");
    }
    thickdir = autor;
    thick_index = 0;
  }
  else if (min_stretch == s_stretch) {
    if ((min_stretch / r_stretch <= 0.5) || (min_stretch / t_stretch <=0.5)) {
      dserror("Solid-Shell element geometry has not a shell aspect ratio");
    }
    thickdir = autos;
    thick_index = 1;
  }
  else if (min_stretch == t_stretch) {
    if ((min_stretch / r_stretch <= 0.5) || (min_stretch / s_stretch <=0.5)) {
      dserror("Solid-Shell element geometry has not a shell aspect ratio");
    }
    thickdir = autot;
    thick_index = 2;
  }

  // thickness-vector in parameter-space, has 1.0 in thickness-coord
  Epetra_SerialDenseVector loc_thickvec(3);
  Epetra_SerialDenseVector glo_thickvec(3);
  loc_thickvec(thick_index) = 1.0;
  // thickness-vector in global coord is J times local thickness-vector
  glo_thickvec.Multiply('T','N',1.0,jac0,loc_thickvec,0.0);
  // return doubles of thickness-vector
  thickvec_.resize(3);
  thickvec_[0] = glo_thickvec(0); thickvec_[1] = glo_thickvec(1); thickvec_[2] = glo_thickvec(2);
  

  // set fiber direction for anisotropic material law
  sosh8_setcylinderfiberdirection(xrefe);
  //cout << thickdir;
  return thickdir;
}


void DRT::ELEMENTS::So_sh8::sosh8_setcylinderfiberdirection(const Epetra_SerialDenseMatrix xrefe)
{
//  // calculate hex8 midpoint
//  Epetra_SerialDenseVector midpoint(3);
//  Epetra_SerialDenseVector diagonal(3);
//  for (int i = 0; i < 3; ++i){
//    diagonal(i) = xrefe(0,i) - xrefe(4,i);
//    midpoint(i) = xrefe(0,i);
//  }
//  diagonal.Scale(0.5);
//  midpoint += diagonal;
//  // cylinder coordinates
//  //double radius = sqrt(midpoint(0)*midpoint(0) + midpoint(1)*midpoint(1));
  Epetra_SerialDenseVector midpoint(3);
  midpoint(0) = soh8_ElementCenterRefeCoords()[0];
  midpoint(1) = soh8_ElementCenterRefeCoords()[1];
  midpoint(2) = soh8_ElementCenterRefeCoords()[2];
  midpoint.Scale(midpoint.Norm2());
  //double theta = atan(midpoint(1)/midpoint(0));
  // tangent circumferential direction = - r_unit
  Epetra_SerialDenseVector circ(3);
  //circ(0) = - cos(theta); circ(1) = - sin(theta);
  circ(0) = - midpoint(1); circ(1) = midpoint(0);
  // rotate by 45°
  fiberdirection_.resize(3);
  fiberdirection_[0] = circ(0);
  fiberdirection_[1] = circ(1);
  fiberdirection_[2] = sqrt(circ(0)*circ(0) + circ(1)*circ(1));
  
//  fiberdirection_[0] = circ(0);
//  fiberdirection_[1] = circ(1);
//  fiberdirection_[2] = circ(2);
  
  return;
  
}

void DRT::ELEMENTS::So_sh8::sosh8_gmshplotlabeledelement(const int LabelIds[NUMNOD_SOH8])
{
  stringstream filename;
  filename << "solidelement" << this->Id() << ".gmsh";
  ofstream f_system("solidelement.gmsh");
  stringstream gmshfilecontent;
  gmshfilecontent << "View \" One Solid Element \" {" << endl;
  gmshfilecontent << IO::GMSH::elementToString(this->thickdir_, this) << endl;
  // plot vector from 1st node to 5th node which is parametric t-dir
  vector<double> X15(3);
  X15[0] = this->Nodes()[4]->X()[0] - this->Nodes()[0]->X()[0];
  X15[1] = this->Nodes()[4]->X()[1] - this->Nodes()[0]->X()[1];
  X15[2] = this->Nodes()[4]->X()[2] - this->Nodes()[0]->X()[2];
  gmshfilecontent << "VP(" << scientific << this->Nodes()[0]->X()[0] << ",";
  gmshfilecontent << scientific << this->Nodes()[0]->X()[1] << ",";
  gmshfilecontent << scientific << this->Nodes()[0]->X()[2] << ")";
  gmshfilecontent << "{" << scientific << X15[0] << "," << X15[1] << "," << X15[2] << "};" << endl;
  gmshfilecontent << "};" << endl;
  gmshfilecontent << "View \" LabelIds \" {" << endl;
  for (int i=0; i<NUMNOD_SOH8; ++i) {
    gmshfilecontent << "SP(" << scientific << this->Nodes()[i]->X()[0] << ",";
    gmshfilecontent << scientific << this->Nodes()[i]->X()[1] << ",";
    gmshfilecontent << scientific << this->Nodes()[i]->X()[2] << ")";
    gmshfilecontent << "{" << LabelIds[i] << "};" << endl;
  }
  gmshfilecontent << "};" << endl;
  gmshfilecontent << "View \" I order \" {" << endl;
  for (int i=0; i<NUMNOD_SOH8; ++i) {
    gmshfilecontent << "SP(" << scientific << this->Nodes()[i]->X()[0] << ",";
    gmshfilecontent << scientific << this->Nodes()[i]->X()[1] << ",";
    gmshfilecontent << scientific << this->Nodes()[i]->X()[2] << ")";
    gmshfilecontent << "{" << i << "};" << endl;
  }
  gmshfilecontent << "};" << endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
}

void DRT::ELEMENTS::Sosh8Register::sosh8_gmshplotdis(const DRT::Discretization& dis)
{
  ofstream f_system("solidelements.gmsh");
  stringstream gmshfilecontent;
  gmshfilecontent << "View \" Solid Elements thickdir_ \" {" << endl;
  // plot elements
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_sosh8) continue;
    DRT::ELEMENTS::So_sh8* actele = dynamic_cast<DRT::ELEMENTS::So_sh8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8* failed");
    // plot elements
    gmshfilecontent << IO::GMSH::elementToString(actele->thickdir_, actele) << endl;
    // plot vector from 1st node to 5th node which is parametric t-dir
    vector<double> X15(3);
    X15[0] = actele->Nodes()[4]->X()[0] - actele->Nodes()[0]->X()[0];
    X15[1] = actele->Nodes()[4]->X()[1] - actele->Nodes()[0]->X()[1];
    X15[2] = actele->Nodes()[4]->X()[2] - actele->Nodes()[0]->X()[2];
    gmshfilecontent << "VP(" << scientific << actele->Nodes()[0]->X()[0] << ",";
    gmshfilecontent << scientific << actele->Nodes()[0]->X()[1] << ",";
    gmshfilecontent << scientific << actele->Nodes()[0]->X()[2] << ")";
    gmshfilecontent << "{" << scientific << X15[0] << "," << X15[1] << "," << X15[2] << "};" << endl;
  }
  gmshfilecontent << "};" << endl;
  gmshfilecontent << "View \" Solid Elements Id \" {" << endl;
  // plot elements
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_sosh8) continue;
    DRT::ELEMENTS::So_sh8* actele = dynamic_cast<DRT::ELEMENTS::So_sh8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_sh8* failed");
    // plot elements
    gmshfilecontent << IO::GMSH::elementToString(actele->LID(), actele) << endl;
  }
  gmshfilecontent << "};" << endl;
  // plot vectors
  gmshfilecontent << "View \" Thickness Vectors \" {" << endl;
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_sosh8) continue;
    DRT::ELEMENTS::So_hex8* actele = dynamic_cast<DRT::ELEMENTS::So_hex8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8* failed");
    // plot vector in center of elements
    const vector<double> pv = actele->GetThickvec();
    vector<double> ec = actele->soh8_ElementCenterRefeCoords();
    //gmshfilecontent << "VP(0,0,0){1,1.3,1.7};" << endl;
    gmshfilecontent << "VP(" << scientific << ec[0] << "," << ec[1] << "," << ec[2] << ")";
    gmshfilecontent << "{" << scientific << pv[0] << "," << pv[1] << "," << pv[2] << "};" << endl; 
  }
  gmshfilecontent << "};" << endl;
  // plot vectors
  gmshfilecontent << "View \" Fiber Direction Vectors \" {" << endl;
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->Type() != DRT::Element::element_sosh8) continue;
    DRT::ELEMENTS::So_hex8* actele = dynamic_cast<DRT::ELEMENTS::So_hex8*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex8* failed");
    // plot vector in center of elements
    vector<double> ec = actele->soh8_ElementCenterRefeCoords();
    RefCountPtr<MAT::Material> mat = actele->Material();
    MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());
    vector<double> pv(3);
    if (anba->GlobalFiberDirection()){
      pv = anba->ReturnGlobalFiberDirection();
    }
    else{
      pv = actele->GetFibervec();
    }
      
    //gmshfilecontent << "VP(0,0,0){1,1.3,1.7};" << endl;
    gmshfilecontent << "VP(" << scientific << ec[0] << "," << ec[1] << "," << ec[2] << ")";
    gmshfilecontent << "{" << scientific << pv[0] << "," << pv[1] << "," << pv[2] << "};" << endl; 
  }
  gmshfilecontent << "};" << endl;
  f_system << gmshfilecontent.str();
  f_system.close();
  return;
}

#endif
#endif
