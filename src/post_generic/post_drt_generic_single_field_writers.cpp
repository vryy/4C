/*!
  \file post_drt_generic_single_field_writers.cpp

  \brief main routine of the Ensight filter

  <pre>
 Maintainer: Ulrich Kuettler
 kuettler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/kuettler
 089 - 289-15238
 </pre>

*/



#include "post_drt_generic_single_field_writers.H"
#include <string>



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("displacement", "displacement", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("prolongated_gauss_2PK_stresses_xyz", "prolongated_gauss_2PK_stresses_xyz", nodebased,6);
  //EnsightWriter::WriteResult("velocity", "velocity", dofbased, field->problem()->num_dim());
  //EnsightWriter::WriteResult("acceleration", "acceleration", dofbased, field->problem()->num_dim());

  // contact and meshtying results
  EnsightWriter::WriteResult("norcontactstress", "norcontactstress", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("tancontactstress", "tancontactstress", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("interfacetraction", "interfacetraction", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("slaveforces", "slaveforces", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("masterforces", "masterforces", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("norslaveforce", "norslaveforce", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("tanslaveforce", "tanslaveforce", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("normasterforce", "normasterforce", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("tanmasterforce", "tanmasterforce", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("wear", "wear", dofbased, field->problem()->num_dim());
 
  // thermo results
  EnsightWriter::WriteResult("temperature", "temperature", nodebased, 1);
  
  // one-dimensional artery
  EnsightWriter::WriteResult("one_d_artery_pressure", "pressure", nodebased, 1);
  EnsightWriter::WriteResult("one_d_artery_flow", "flow", nodebased, 1);

  // reduced dimensional airway
  EnsightWriter::WriteResult("pnp", "pressure", dofbased, 1);
  EnsightWriter::WriteResult("NodeIDs", "NodeIDs", dofbased, 1);
  EnsightWriter::WriteResult("radii", "radii", dofbased, 1);
  EnsightWriter::WriteResult("acini_volume", "acini_volume", dofbased, 1);
  EnsightWriter::WriteResult("acin_bc", "acini_bc", elementbased, 1);
  EnsightWriter::WriteResult("qin_np", "flow_in", elementbased, 1);
  EnsightWriter::WriteResult("qout_np", "flow_out", elementbased, 1);
  EnsightWriter::WriteResult("generations", "generations", elementbased, 1);

  // additional forces due to lung fsi (volume constraint)
  EnsightWriter::WriteResult("Add_Forces", "Add_Forces", dofbased, field->problem()->num_dim());

  EnsightWriter::WriteElementResults(field); //To comment
  if (stresstype_!="none")
  {
    // although appearing here twice, only one function call to PostStress
    // is really postprocessing Gauss point stresses, since only _either_
    // Cauchy _or_ 2nd Piola-Kirchhoff stresses are written during simulation!
    PostStress("gauss_cauchy_stresses_xyz", stresstype_);
    PostStress("gauss_2PK_stresses_xyz", stresstype_);
    // in case of coupled problem (e.g. TSI) write the stresses arising due to
    // coupling to another field
    PostStress("gauss_cauchy_coupling_stresses_xyz", stresstype_);
    PostStress("gauss_2PK_coupling_stresses_xyz", stresstype_);
  }
  if (straintype_!="none")
  {
    // although appearing here twice, only one function call to PostStress
    // is really postprocessing Gauss point strains, since only _either_
    // Green-Lagrange _or_ Euler-Almansi strains are written during simulation!
    PostStress("gauss_GL_strains_xyz", straintype_);
    PostStress("gauss_EA_strains_xyz", straintype_);
  }
}



