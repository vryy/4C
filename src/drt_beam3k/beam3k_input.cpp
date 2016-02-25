/*!----------------------------------------------------------------------
\file beam3k_input.cpp

\brief three dimensional nonlinear Kirchhoff beam element based on a C1 curve

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*-----------------------------------------------------------------------------------------------------------*/

#include "../drt_beam3k/beam3k.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/largerotations.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam3k::ReadElement(const std::string& eletype,
                                       const std::string& distype,
                                       DRT::INPUT::LineDefinition* linedef)
{

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  int rotvec = 0;
  linedef->ExtractInt("ROTVEC",rotvec);

  int wk = 0;
  linedef->ExtractInt("WK",wk);

  if (rotvec==0)
    rotvec_=false;
  else if (rotvec==1)
    rotvec_=true;
  else
    dserror("The variable ROTVEC can only take on the values 0 (tangent vectors as nodal DoFs) and "
            "1 (rotation vectors as nodal DoFs)!");

  if (wk==0)
    weakkirchhoff_=false;
  else if (wk==1)
  {
    weakkirchhoff_=true;
    #ifdef CONSISTENTSPINSK
      dserror("The flag CONSISTENTSPINSK is only possible for strong Kirchhoff constraint enforcement (weakkirchhoff_=false)");
    #endif
  }
  else
    dserror("The variable WK can only take on the values 0 (Kirchhoff constraint enforced in a strong manner) and "
            "1 (Kirchhoff constraint enforced in a weak manner)!");

  linedef->ExtractDouble("CROSS",crosssec_);

  /*read beam moments of inertia of area; currently the Beam3k element works only with rotationally symmetric
   * crosssection so that the moment of inertia of area around both principal axes can be expressed by one input
   * number I_; however, the implementation itself is a general one and works also for other cases;*/

  //currently only rotationally symmetric profiles for beam --> Iyy = Izz
  linedef->ExtractDouble("MOMINY",Iyy_);
  linedef->ExtractDouble("MOMINZ",Izz_);

  //torsional moment of inertia
  linedef->ExtractDouble("MOMINPOL",Irr_);

  linedef->ExtractDouble("IT",inertscaletrans_);
  linedef->ExtractDouble("IR1",inertscalerot1_);
  linedef->ExtractDouble("IR2",inertscalerot2_);

  //extract triads at element nodes in reference configuration as rotation vectors and save them as quaternions at each node, respectively
  std::vector<double> nodal_thetas;
  linedef->ExtractDoubleVector("TRIADS",nodal_thetas);
  theta0_.resize(COLLOCATION_POINTS);
  for(int i=0; i<COLLOCATION_POINTS; i++)
  {
    for(int j=0; j<3; j++)
      (theta0_[i])(j) = nodal_thetas[3*i+j];
  }

  return true;
}
