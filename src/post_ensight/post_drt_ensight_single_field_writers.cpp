/*!
  \file post_drt_ensight_single_field_writers.cpp

  \brief main routine of the Ensight filter

  <pre>
  Maintainer: Axel Gerstenberger
  gerstenberger@lnm.mw.tum.de
  http://www.lnm.mw.tum.de/Members/gerstenberger
  089 - 289-15236
  </pre>

*/

#ifdef CCADISCRET

#include "post_drt_ensight_single_field_writers.H"
#include <string>



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("displacement", "displacement", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("norcontactstress", "norcontactstress", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("tancontactstress", "tancontactstress", dofbased, field->problem()->num_dim());
  //EnsightWriter::WriteResult("velocity", "velocity", dofbased, field->problem()->num_dim());
  //EnsightWriter::WriteResult("acceleration", "acceleration", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteElementResults(field); //To comment
  if (stresstype_!="none")
  {
    // although appearing here twice, only one function call to PostStress
    // is really postprocessing Gauss point stresses, since only _either_
    // Cauchy _or_ 2nd Piola-Kirchhoff stresses are written during simulation!
    PostStress("gauss_cauchy_stresses_xyz", stresstype_);
    PostStress("gauss_2PK_stresses_xyz", stresstype_);
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("averaged_pressure", "averaged_pressure", dofbased, 1);
  EnsightWriter::WriteResult("averaged_velnp", "averaged_velocity", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("velnp", "velocity", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("pressure", "pressure", dofbased, 1);
  EnsightWriter::WriteResult("density", "density", dofbased, 1);
  EnsightWriter::WriteResult("residual", "residual", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("dispnp", "ale_displacement", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("traction", "traction", dofbased, field->problem()->num_dim());
  WriteElementResults(field);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteAllResults(PostField* field)
{
  // XFEM has changing number of degrees of freedoms
  // - restart vectors are of changing size
  // - output vectors (e.g. Paraview) have fixed size with 4 DOF per node)
  //   and are named "*_smoothed" for this reason (no integration cells)
  // calling both at the same time will crash, since restart vectors do not fit
  // the 4 DOF per node pattern. BACI will produce consistent naming now, but for old
  // data you can switch the naming convention here (old naming will be removed soon)

  const bool consistent_naming = true;

  if (consistent_naming)
  {
    EnsightWriter::WriteResult("velocity_smoothed", "velocity_smoothed", dofbased, field->problem()->num_dim());
    EnsightWriter::WriteResult("pressure_smoothed", "pressure_smoothed", dofbased, 1);

    // for diffusion problem, if ever needed
    EnsightWriter::WriteResult("temperature_smoothed", "temperature_smoothed", dofbased, 1);

  }
  else
  {
    cout << "Depreciated naming convention!!!" << endl;
    // note old output files might still use the name names velnp and pressure
    // just turn the following lines on
    EnsightWriter::WriteResult("velnp", "velocity", dofbased, field->problem()->num_dim());
    EnsightWriter::WriteResult("pressure", "pressure", dofbased, 1);
  }

  WriteElementResults(field);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void InterfaceEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("idispnp", "idispnp", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("idispn", "idispn", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("ivelnp", "ivelnp", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("iveln", "iveln", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("ivelnm", "ivelnm", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("iaccn", "iaccn", dofbased, field->problem()->num_dim());
  EnsightWriter::WriteResult("itrueresnp", "itrueresnp", dofbased, field->problem()->num_dim());
  WriteElementResults(field);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DGFEMFluidEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("velnp", "elemean_velocity", elementdof, field->problem()->num_dim());
  EnsightWriter::WriteResult("velnp", "elemean_pressure", elementdof, 1, field->problem()->num_dim());
  WriteElementResults(field);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleEnsightWriter::WriteAllResults(PostField* field)
{
  EnsightWriter::WriteResult("dispnp", "displacement", dofbased, field->problem()->num_dim());
  WriteElementResults(field);
}


/*----------------------------------------------------------------------*
|                                                           gjb 12/07   |
\*----------------------------------------------------------------------*/
void ScaTraEnsightWriter::WriteAllResults(PostField* field)
{
  //compute number of dofs per node (ask the first node)
  int numdofpernode = field->discretization()->NumDof(field->discretization()->lRowNode(0));

  // write results for each transported scalar
  if (numdofpernode == 1)
  {
    EnsightWriter::WriteResult("phinp","phi",dofbased,1);
    // write flux vectors (always 3D)
    EnsightWriter::WriteResult("flux", "flux", nodebased, 3);
  }
  else
  {
    for(int k = 1; k <= numdofpernode; k++)
    {
      ostringstream temp;
      temp << k;
      string name = "phi_"+temp.str();
      EnsightWriter::WriteResult("phinp", name, dofbased, 1,k-1);
      // write flux vectors (always 3D)
      EnsightWriter::WriteResult("flux_"+name, "flux_"+name, nodebased, 3);
    }
  }

  // write velocity field (always 3D)
  EnsightWriter::WriteResult("convec_velocity", "velocity", nodebased, 3);

  // write displacement field (always 3D)
  EnsightWriter::WriteResult("dispnp", "ale-displacement", nodebased, 3);

  // write element results (e.g. element owner)
  WriteElementResults(field);
}


/*----------------------------------------------------------------------*
|                                                             gjb 09/08 |
\*----------------------------------------------------------------------*/
void ElchEnsightWriter::WriteAllResults(PostField* field)
{
  //compute number of dofs per node (ask the first node)
  int numdofpernode = field->discretization()->NumDof(field->discretization()->lRowNode(0));

  // write results for each transported scalar
  if (numdofpernode < 3)
    dserror("Problemtype ELCH has at least 3 DOF per node, but got: %d",numdofpernode);
  else
  {
    // do the ion concentrations first
    for(int k = 1; k < numdofpernode; k++)
    {
      ostringstream temp;
      temp << k;
      string name = "c_"+temp.str();
      EnsightWriter::WriteResult("phinp", name, dofbased, 1,k-1);
      // write flux vectors (always 3D)
      EnsightWriter::WriteResult("flux_phi_"+temp.str(), "flux_"+name, nodebased, 3);
    }
    // finally, handle the electric potential
    EnsightWriter::WriteResult("phinp", "phi", dofbased, 1,numdofpernode-1);
  }

  // write velocity field (always 3D)
  EnsightWriter::WriteResult("convec_velocity", "velocity", nodebased, 3);

  // write displacement field (always 3D)
  EnsightWriter::WriteResult("dispnp", "ale-displacement", nodebased, 3);

  // write density distribution
  EnsightWriter::WriteResult("elchdensnp", "density", nodebased, 1);

  // write element results (e.g. element owner)
  WriteElementResults(field);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ThermoEnsightWriter::WriteAllResults(PostField* field)
{
  // number of dofs per node in thermal problems is always 1
  const int numdofpernode = 1;

  // write temperature
  EnsightWriter::WriteResult("temperature", "temperature", dofbased, numdofpernode);

  // write temperature rate
  //EnsightWriter::WriteResult("rate", "rate", dofbased, numdofpernode);

  // write element results (e.g. element owner)
  EnsightWriter::WriteElementResults(field);

  // Start new section 06/11/09
  if (heatfluxtype_ != "none")
  {
    // although appearing here twice, only one function call to PostHeatflux
    // is really postprocessing Gauss point heatfluxes, since only _either_
    // Current _or_ Initial heatfluxes are written during simulation!
    PostHeatflux("gauss_current_heatfluxes_xyz", heatfluxtype_);
    PostHeatflux("gauss_initial_heatfluxes_xyz", heatfluxtype_);
    //EnsightWriter::WriteResult("heatflux", "heatflux", nodebased, 3);
  }
  if (tempgradtype_ != "none")
  {
    // although appearing here twice, only one function call to PostHeatflux
    // is really postprocessing Gauss point temperature gradients, since only _either_
    // Initial _or_ Current temperature gradients are written during simulation!
    PostHeatflux("gauss_current_tempgrad_xyz", tempgradtype_);
    PostHeatflux("gauss_initial_tempgrad_xyz", tempgradtype_);
    //EnsightWriter::WriteResult("tempgrad", "tempgrad", nodebased, 3);
  }
  // end new section 06/11/09

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AnyEnsightWriter::WriteAllResults(PostField* field)
{
  WriteDofResults(field);
  WriteNodeResults(field);
  WriteElementResults(field);
}


#endif
