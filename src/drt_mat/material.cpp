#ifdef CCADISCRET

#include "material.H"
#include "newtonianfluid.H"
#include "stvenantkirchhoff.H"
#include "micromaterial.H"
#include "hyperpolyconvex.H"
#include "neohooke.H"
#include "convecdiffus.H"

extern struct _MATERIAL *mat;



Teuchos::RefCountPtr<MAT::Material> MAT::Material::Factory(int matnum)
{
  // We still query the global variable here because we have no idea
  // which Problem instance we are working at. There might be more
  // that one...
  MATERIAL* actmat = &(mat[matnum-1]);

  switch (actmat->mattyp)
  {
  case m_fluid:
  {
#if 0
    // Newtonian fluid does not store any data. So we could cache the
    // objects and just return the copy we already have. But we have
    // to keep in mind that there could be many Problems and therefore
    // many sets of materials. So the material number is not a
    // sufficient criteria.
    std::map<MATERIAL*,Teuchos::RefCountPtr<Material> > mats;
    if (mats.find(actmat)==mats.end())
    {
      mats[actmat] = Teuchos::rcp(new NewtonianFluid(actmat));
    }
    return mats[actmat];
#else
    return Teuchos::rcp(new NewtonianFluid(actmat));
#endif
  }
  case m_stvenant:
    return Teuchos::rcp(new StVenantKirchhoff(actmat));
  case m_struct_multiscale:
    return Teuchos::rcp(new MicroMaterial(actmat));
  case m_hyper_polyconvex:
    return Teuchos::rcp(new HyperPolyconvex(actmat));
  case m_neohooke:
    return Teuchos::rcp(new NeoHooke(actmat));
  case m_condif:
    return Teuchos::rcp(new ConvecDiffus(actmat));
  case m_pl_mises_3D:
  case m_pl_mises:
  case m_pl_hoff:
  case m_damage:
  case m_pl_foam:
  case m_pl_mises_ls:
  case m_pl_dp:
  case m_pl_epc:
  case m_pl_epc3D:
  case m_stvenpor:
  case m_pl_por_mises:
  case m_compogden:
  case m_viscohyper:
  case m_pl_hash:
  case m_el_orth:
  case m_mfoc:
  case m_mfcc:
  case m_nhmfcc:
  case m_multi_layer:
  case m_ifmat:
  case m_interf_therm:
  case m_dam_mp:
  case m_damage_ge:
  case m_th_fourier_iso:
  case m_th_fourier_gen:
  case m_vp_robinson:
  default:
    dserror("unknown material type %d", actmat->mattyp);
  }

  return Teuchos::null;
}


#endif
