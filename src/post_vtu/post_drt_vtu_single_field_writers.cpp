/*!
  \file post_drt_vtu_single_field_writers.cpp

  \brief main routine of the VTU filter

  <pre>
  Maintainer: Martin Kronbichler
  kronbichler@lnm.mw.tum.de
  http://www.lnm.mw.tum.de
  089 - 289-15235
  </pre>

*/



#include "post_drt_vtu_single_field_writers.H"
#include "../post_drt_common/post_drt_common.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret.H"
#include <string>



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureVtuWriter::WriteAllResults(PostField* field)
{
  VtuWriter::WriteResult("displacement", "displacement", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("prolongated_gauss_2PK_stresses_xyz", "prolongated_gauss_2PK_stresses_xyz", nodebased,6);
  VtuWriter::WriteResult("prolongated_gauss_GL_strains_xyz", "prolongated_gauss_GL_strains_xyz", nodebased,6);
  if(field->problem()->struct_vel_acc() == "yes")
  {
    VtuWriter::WriteResult("velocity", "velocity", dofbased, field->problem()->num_dim());
    VtuWriter::WriteResult("acceleration", "acceleration", dofbased, field->problem()->num_dim());
  }
  if (field->problem()->struct_mat_disp() == "yes")
    VtuWriter::WriteResult("material_displacement", "material_displacement", dofbased, field->problem()->num_dim());

  // Statistical Output from MLMC
  VtuWriter::WriteResult("mean_displacements", "mean_displacement", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("variance_displacements", "variance_displacement", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("mean_gauss_2PK_stresses_xyz", "mean_gauss_2PK_stresses_xyz", nodebased,6);
  VtuWriter::WriteResult("variance_gauss_2PK_stresses_xyz", "variance_gauss_2PK_stresses_xyz", nodebased,6);
  VtuWriter::WriteResult("mean_gauss_GL_strain_xyz", "mean_gauss_GL_strain_xyz", nodebased,6);
  VtuWriter::WriteResult("variance_gauss_GL_strain_xyz", "variance_gauss_GL_strain_xyz", nodebased,6);

  VtuWriter::WriteResult("diff_to_ll_displacement", "diff_to_ll_displacement", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("diff_to_ll_prolongated_gauss_2PK_stresses_xyz", "diff_to_ll_prolongated_gauss_2PK_stresses_xyz", nodebased,6);
  VtuWriter::WriteResult("diff_to_ll_prolongated_gauss_GL_strains_xyz", "diff_to_ll_prolongated_gauss_GL_strains_xyz", nodebased,6);

  VtuWriter::WriteResult("diff_mean_displacements", "diff_mean_displacement", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("diff_variance_displacements", "diff_variance_displacement", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("diff_mean_gauss_2PK_stresses_xyz", "diff_mean_gauss_2PK_stresses_xyz", nodebased,6);
  VtuWriter::WriteResult("diff_variance_gauss_2PK_stresses_xyz", "diff_variance_gauss_2PK_stresses_xyz", nodebased,6);
  VtuWriter::WriteResult("diff_mean_gauss_GL_strain_xyz", "diff_mean_gauss_GL_strain_xyz", nodebased,6);
  VtuWriter::WriteResult("diff_variance_gauss_GL_strain_xyz", "diff_variance_gauss_GL_strain_xyz", nodebased,6);


   // contact and meshtying results
  VtuWriter::WriteResult("activeset", "activeset", nodebased,1);
  VtuWriter::WriteResult("norcontactstress", "norcontactstress", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("tancontactstress", "tancontactstress", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("interfacetraction", "interfacetraction", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("slaveforces", "slaveforces", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("masterforces", "masterforces", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("norslaveforce", "norslaveforce", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("tanslaveforce", "tanslaveforce", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("normasterforce", "normasterforce", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("tanmasterforce", "tanmasterforce", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("wear", "wear", dofbased, field->problem()->num_dim());

  // one-dimensional artery
  VtuWriter::WriteResult("one_d_artery_pressure", "pressure", nodebased, 1);
  VtuWriter::WriteResult("one_d_artery_flow", "flow", nodebased, 1);
  VtuWriter::WriteResult("one_d_artery_area", "area", nodebased, 1);

  // reduced dimensional airway
  VtuWriter::WriteResult("pnp", "pressure", dofbased, 1);
  VtuWriter::WriteResult("NodeIDs", "NodeIDs", dofbased, 1);
  VtuWriter::WriteResult("radii", "radii", dofbased, 1);
//  VtuWriter::WriteResult("juncVolMix", "juncVolMix", dofbased, 1);
//  VtuWriter::WriteResult("forwardVolume","forwardVolume", dofbased, 1);
  VtuWriter::WriteResult("scatraO2np", "scatraO2", dofbased, 1);
  VtuWriter::WriteResult("PO2", "PO2", dofbased, 1);
  VtuWriter::WriteResult("dVO2","dVO2", dofbased, 1);
    VtuWriter::WriteResult("AcinarPO2","AcinarPO2", dofbased, 1);
//  VtuWriter::WriteResult("e1scatranp","e1scatranp",elementbased, 1);
//  VtuWriter::WriteResult("e2scatranp","e2scatranp",elementbased, 1);
  VtuWriter::WriteResult("juncVolMix","juncVolMix", dofbased, 1);

  VtuWriter::WriteResult("acini_vnp", "acini_volume",elementbased, 1);
  VtuWriter::WriteResult("acini_volumetric_strain", "acini_volume_strain",elementbased, 1);
  VtuWriter::WriteResult("acin_bc", "acini_bc",elementbased, 1);
  VtuWriter::WriteResult("acini_v0","acini_volume0",elementbased, 1);

  VtuWriter::WriteResult("qin_np", "flow_in", elementbased, 1);
  VtuWriter::WriteResult("qout_np", "flow_out", elementbased, 1);
  VtuWriter::WriteResult("generations", "generations", elementbased, 1);
  VtuWriter::WriteResult("elemVolumenp", "elemVolumenp", elementbased, 1);


  // additional forces due to lung fsi (volume constraint)
  VtuWriter::WriteResult("Add_Forces", "Add_Forces", dofbased, field->problem()->num_dim());

  // Lagrange multiplier at the interface in monolithic fsi
  VtuWriter::WriteResult("fsilambda", "fsilambda", dofbased, field->problem()->num_dim());

  //additional output for biofilm problems
  VtuWriter::WriteResult("str_growth_displ", "str_growth_displ", dofbased, field->problem()->num_dim());

  //additional output for poro problems
  VtuWriter::WriteResult("porosity_p1", "porosity_p1", dofbased, 1);

  VtuWriter::WriteElementResults(field); //To comment
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
    PostStress("gauss_LOG_strains_xyz", straintype_);
    // the same for plastic strains
    PostStress("gauss_pl_GL_strains_xyz", straintype_);
    PostStress("gauss_pl_EA_strains_xyz", straintype_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidVtuWriter::WriteAllResults(PostField* field)
{
  VtuWriter::WriteResult("averaged_pressure", "averaged_pressure", dofbased, 1);
  VtuWriter::WriteResult("averaged_velnp", "averaged_velocity", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("averaged_scanp", "averaged_scalar_field", dofbased, 1,field->problem()->num_dim());
  VtuWriter::WriteResult("filteredvel", "filteredvel", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("fsvelaf", "fsvel" , dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("velnp", "velocity", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("pressure", "pressure", dofbased, 1);
  VtuWriter::WriteResult("scalar_field", "scalar_field", dofbased, 1);
  VtuWriter::WriteResult("residual", "residual", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("dispnp", "ale_displacement", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("idispnfull", "ale_idisp", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("traction", "traction", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("wss", "wss", dofbased, field->problem()->num_dim());
  //  VtuWriter::WriteResult("radii", "radii", nodebased, 1);
  VtuWriter::WriteResult("par_vel", "par_vel", dofbased, field->problem()->num_dim());

  // additional output for turbulent flows (subfilter/-gridstress)
  VtuWriter::WriteResult("sfs11", "sfs11", nodebased, 1);
  VtuWriter::WriteResult("sfs12", "sfs12", nodebased, 1);
  VtuWriter::WriteResult("sfs13", "sfs13", nodebased, 1);
  VtuWriter::WriteResult("sfs22", "sfs22", nodebased, 1);
  VtuWriter::WriteResult("sfs23", "sfs23", nodebased, 1);
  VtuWriter::WriteResult("sfs33", "sfs33", nodebased, 1);

  // additional forces due to lung fsi (volume constraint)
  VtuWriter::WriteResult("Add_Forces", "Add_Forces", dofbased, field->problem()->num_dim());

  // additional output for poro problems
  VtuWriter::WriteResult("convel", "convective_velocity", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("gridv", "grid_velocity", dofbased, field->problem()->num_dim());

  // additional fields due to adjoint equations
  VtuWriter::WriteResult("adjoint_velnp", "adjoint_velocity", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("adjoint_pressure", "adjoint_pressure", dofbased, 1);

  // additional fields for meshfree problems
  VtuWriter::WriteResult("velatmeshfreenodes", "velocity_atnodes", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("pressureatmeshfreenodes", "pressure_atnodes", dofbased, 1);

  // Lagrange multiplier at the interface in monolithic fsi
  VtuWriter::WriteResult("fsilambda", "fsilambda", dofbased, field->problem()->num_dim());

  //additional output for biofilm problems
  VtuWriter::WriteResult("fld_growth_displ", "fld_growth_displ", dofbased, field->problem()->num_dim());

  WriteElementResults(field);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidVtuWriter::WriteAllResults(PostField* field)
{
  // XFEM has changing number of degrees of freedoms
  // - restart vectors are of changing size
  // - output vectors (e.g. Paraview) have fixed size with 4 DOF per node)
  //   and are named "*_smoothed" for this reason (no integration cells)
  // calling both at the same time will crash, since restart vectors do not fit
  // the 4 DOF per node pattern.

  VtuWriter::WriteResult("velocity_smoothed", "velocity_smoothed", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("pressure_smoothed", "pressure_smoothed", dofbased, 1);
  VtuWriter::WriteResult("averaged_velnp", "averaged_velocity", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("averaged_pressure", "averaged_pressure", dofbased, 1);
  VtuWriter::WriteResult("fsvelocity", "fsvelocity", dofbased, field->problem()->num_dim());

  WriteElementResults(field);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void InterfaceVtuWriter::WriteAllResults(PostField* field)
{
  VtuWriter::WriteResult("idispnp", "idispnp", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("idispn", "idispn", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("ivelnp", "ivelnp", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("iveln", "iveln", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("ivelnm", "ivelnm", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("iaccn", "iaccn", dofbased, field->problem()->num_dim());
  VtuWriter::WriteResult("itrueresnp", "itrueresnp", dofbased, field->problem()->num_dim());
  WriteElementResults(field);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleVtuWriter::WriteAllResults(PostField* field)
{
  VtuWriter::WriteResult("dispnp", "displacement", dofbased, field->problem()->num_dim());
  WriteElementResults(field);
}


/*----------------------------------------------------------------------*
|                                                           gjb 12/07   |
\*----------------------------------------------------------------------*/
void ScaTraVtuWriter::WriteAllResults(PostField* field)
{
  //compute number of dofs per node (ask the first node)
  int numdofpernode = field->discretization()->NumDof(field->discretization()->lRowNode(0));

  // write results for each transported scalar
  for(int k = 1; k <= numdofpernode; k++)
  {
    std::ostringstream temp;
    temp << k;
    std::string name = "phi_"+temp.str();
    VtuWriter::WriteResult("phinp", name, dofbased, 1,k-1);
    VtuWriter::WriteResult("averaged_phinp", "averaged_"+name, dofbased, 1,k-1);
    // additional fields for meshfree problems
    VtuWriter::WriteResult("phiatmeshfreenodes", name+"_atnodes", dofbased, 1,k-1);
    // intermediate work-around for nurbs discretizations (no normal vectors applied)
    VtuWriter::WriteResult("normalflux","normalflux"+name,dofbased,1,k-1);
    // write flux vectors (always 3D)
    VtuWriter::WriteResult("flux_"+name, "flux_"+name, nodebased, 3);
    VtuWriter::WriteResult("activation_time_np","act_time",dofbased,1);
  }

  // write velocity field (always 3D)
  VtuWriter::WriteResult("convec_velocity", "velocity", nodebased, 3);

  // write displacement field (always 3D)
  VtuWriter::WriteResult("dispnp", "ale-displacement", nodebased, 3);

  // write optimization field
  VtuWriter::WriteResult("x_mma","optimization_variable", dofbased,1);

  //additional output for biofilm problems
  VtuWriter::WriteResult("scfld_growth_displ", "scfld_growth_displ", nodebased, 3);
  VtuWriter::WriteResult("scstr_growth_displ", "scstr_growth_displ", nodebased, 3);

  // write element results (e.g. element owner)
  WriteElementResults(field);
}


/*----------------------------------------------------------------------*
|                                                             gjb 09/08 |
\*----------------------------------------------------------------------*/
void ElchVtuWriter::WriteAllResults(PostField* field)
{
  //compute number of dofs per node (ask the first node)
  int numdofpernode = field->discretization()->NumDof(field->discretization()->lRowNode(0));

  // write results for each transported scalar
  if (numdofpernode == 1)
  {
    // do the single ion concentration
    std::string name = "c_1";
    VtuWriter::WriteResult("phinp", name, dofbased, 1, 0);
    // write flux vectors (always 3D)
    VtuWriter::WriteResult("flux", "flux", nodebased, 3);

    // there is no electric potential in this special case

    // temporal mean field from turbulent statistics (if present)
    VtuWriter::WriteResult("averaged_phinp", "averaged_"+name, dofbased, 1, 0);
  }
  else
  {
    // do the ion concentrations first
    for(int k = 1; k < numdofpernode; k++)
    {
      std::ostringstream temp;
      temp << k;
      std::string name = "c_"+temp.str();
      VtuWriter::WriteResult("phinp", name, dofbased, 1,k-1);
      // write flux vectors (always 3D)
      VtuWriter::WriteResult("flux_phi_"+temp.str(), "flux_"+name, nodebased, 3);

      // temporal mean field from turbulent statistics (if present)
      VtuWriter::WriteResult("averaged_phinp", "averaged_"+name, dofbased, 1, k-1);
    }
    // finally, handle the electric potential
    VtuWriter::WriteResult("phinp", "phi", dofbased, 1,numdofpernode-1);
    // temporal mean field from turbulent statistics (if present)
    VtuWriter::WriteResult("averaged_phinp", "averaged_phi", dofbased, 1,numdofpernode-1);
  }

  // write velocity field (always 3D)
  VtuWriter::WriteResult("convec_velocity", "velocity", nodebased, 3);

  // write displacement field (always 3D)
  VtuWriter::WriteResult("dispnp", "ale-displacement", nodebased, 3);

  // write magnetic field (always 3D)
  VtuWriter::WriteResult("magnetic_field", "B", nodebased, 3);

  // write element results (e.g. element owner)
  WriteElementResults(field);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ThermoVtuWriter::WriteAllResults(PostField* field)
{
  // number of dofs per node in thermal problems is always 1
  const int numdofpernode = 1;

  // write temperature
  VtuWriter::WriteResult("temperature", "temperature", dofbased, numdofpernode);

  // write temperature rate
  //VtuWriter::WriteResult("rate", "rate", dofbased, numdofpernode);

  // write element results (e.g. element owner)
  VtuWriter::WriteElementResults(field);

  if (heatfluxtype_ != "none")
  {
    // although appearing here twice, only one function call to PostHeatflux
    // is really postprocessing Gauss point heatfluxes, since only _either_
    // Current _or_ Initial heatfluxes are written during simulation!
    PostHeatflux("gauss_current_heatfluxes_xyz", heatfluxtype_);
    PostHeatflux("gauss_initial_heatfluxes_xyz", heatfluxtype_);
    VtuWriter::WriteResult("heatflux", "heatflux", nodebased, field->problem()->num_dim());
  }
  if (tempgradtype_ != "none")
  {
    // although appearing here twice, only one function call to PostHeatflux
    // is really postprocessing Gauss point temperature gradients, since only _either_
    // Initial _or_ Current temperature gradients are written during simulation!
    PostHeatflux("gauss_current_tempgrad_xyz", tempgradtype_);
    PostHeatflux("gauss_initial_tempgrad_xyz", tempgradtype_);
    VtuWriter::WriteResult("tempgrad", "tempgrad", nodebased, field->problem()->num_dim());
  }

  // write displacement field
  VtuWriter::WriteResult("displacement", "displacement", nodebased, field->problem()->num_dim());

} // ThermoVtuWriter::WriteAllResults

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ParticleVtuWriter::WriteAllResults(PostField* field)
{
  bool fillzeros = true;
  VtuWriter::WriteResult("displacement", "displacement", dofbased, field->problem()->num_dim(), 0, fillzeros);
  VtuWriter::WriteResult("velocity", "velocity", dofbased, field->problem()->num_dim(), 0, fillzeros);
  VtuWriter::WriteResult("acceleration", "acceleration", dofbased, field->problem()->num_dim(), 0, fillzeros);
  VtuWriter::WriteResult("radius", "radius", nodebased, 1);
  VtuWriter::WriteResult("sign", "sign", nodebased, 1);
  VtuWriter::WriteResult("orientation", "orientation", dofbased, field->problem()->num_dim(), 0, fillzeros);
}

/*----------------------------------------------------------------------*
  | write particles as points                               ghamm 02/13 |
  *----------------------------------------------------------------------*/
void ParticleVtuWriter::WriteCells(
  std::ofstream& geofile,
  const Teuchos::RCP<DRT::Discretization> dis,
  const Teuchos::RCP<Epetra_Map>& proc0map
  ) const
{
  CurrentlyNotImplemented();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AcouVtuWriter::WriteAllResults(PostField* field)
{
  VtuWriter::WriteResult("velnp", "velocity", nodebased, 3);
  VtuWriter::WriteResult("pressure", "pressure", nodebased, 1);
  VtuWriter::WriteResult("par_vel", "par_vel", nodebased, 1);
  VtuWriter::WriteResult("pressure_avg", "pressure_avg", elementbased, 1);
  VtuWriter::WriteResult("error", "error", elementbased, 1);
  WriteElementResults(field);
}



