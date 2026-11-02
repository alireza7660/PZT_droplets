"""Dimensionless tables and concise scientific text for the co-directed package."""
import json
import numpy as np
from .combined_flow_prediction import ROOT,write_csv,write_json,sha
from .droplet_dynamics import relaxation_time,reynolds
from .plotting import plt


def dimensionless(out,cfg,p,D0,alignment):
    config_path=(ROOT/cfg['znmf_configuration']).resolve();gcfg=json.loads(config_path.read_text())
    source=config_path.parent.parent/'tables/FOV_1_2_3_centreline_and_width.csv'
    z=np.genfromtxt(source,names=True,delimiter=',',dtype=None,encoding='utf-8')
    valid=(z['status']=='measured')&np.isfinite(z['Uc_m_per_s'])&(z['Uc_m_per_s']>0)
    uc=np.where(valid,z['Uc_m_per_s'],np.nan);f=float(gcfg['forcing_Hz'])
    d=cfg['droplet_diameters_um'][0]*1e-6;tau=float(relaxation_time(d,p));nu=p['mu_air_Pa_s']/p['rho_air_kg_m3']
    y=z['x_over_D0']*D0*1e3+alignment['translation_m'][1]*1e3
    values=dict(x_ZNMF_over_D0=z['x_over_D0'],aligned_y_mm=y,Uc_measured_m_s=z['Uc_m_per_s'],supported=valid,
                Re_D_c=uc*D0/nu,Sr_D_c=f*D0/uc,Stk_D_c=tau*uc/D0,status=z['status'])
    write_csv(out/'data/jet_local_centreline_numbers.csv',values)
    scalars=dict(nominal_diameter_um=[d*1e6],tau_Stokes_s=[tau],tau_Stokes_ms=[tau*1e3],
        jet_actuation_frequency_Hz=[f],acquisition_frequency_Hz_not_used=[gcfg['acquisition_Hz']],D0_m=[D0],
        Stk_f=[tau*f],rho_air_kg_m3=[p['rho_air_kg_m3']],rho_water_kg_m3=[p['rho_water_kg_m3']],
        mu_air_Pa_s=[p['mu_air_Pa_s']],nu_air_m2_s=[nu])
    write_csv(out/'data/scalar_quantities.csv',scalars)
    # Retain the PZT-only diagnostic as explicit supporting reference data, not gas measurement.
    with np.load(ROOT/cfg['pzt_statistics']) as a:
        re=reynolds(d,np.stack([a['mean_u'],a['mean_v']],axis=-1),[0,0],p)
        write_csv(out/'data/supporting_PZT_zero_background_Re.csv',dict(x_mm=a['x_mm'].ravel(),y_mm=a['y_mm'].ravel(),
            mean_u_m_s=a['mean_u'].ravel(),mean_v_m_s=a['mean_v'].ravel(),Re_p_zero_background_reference=re.ravel()))
    fig,axes=plt.subplots(3,1,figsize=(6.6,7.4),sharex=True,layout='constrained')
    for ax,key,label in zip(axes,['Re_D_c','Sr_D_c','Stk_D_c'],[r'$Re_{D,c}\ (-)$',r'$Sr_{D,c}\ (-)$',r'$Stk_{D,c}\ (-)$']):
        ax.plot(y,values[key],color='tab:blue');ax.set_ylabel(label)
    axes[-1].set_xlabel(r'$y\ (\mathrm{mm})$');stem='jet_local_dimensionless_numbers'
    for ext in ('pdf','png'):fig.savefig(out/'figures'/(stem+'.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)
    return dict(tau_Stokes_s=tau,Stk_f=tau*f,actuation_frequency_Hz=f,D0_m=D0,nu_air_m2_s=nu,
        positive_measured_centreline_rows=int(valid.sum()),unsupported_rows=int((~valid).sum()),
        source=str(source),source_sha256=sha(source),definitions='Re_D,c=Uc D0/nu; Sr_D,c=f D0/Uc; Stk_D,c=tau_Stokes Uc/D0; trajectory Stk_D=tau_Stokes norm(gas vector)/D0'),stem


def write_report(out,cfg,result,stems):
    from .interaction_report import write_report as current_report
    return current_report(out,cfg,result,stems)
