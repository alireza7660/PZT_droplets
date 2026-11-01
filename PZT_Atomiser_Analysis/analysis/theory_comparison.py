"""Spherical-droplet references from final SI data; no experimental files rewritten.

Run from workspace root: python -m PZT_Atomiser_Analysis.analysis.theory_comparison
All integrations use SI units; plotted heights alone are converted to mm.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .plotting import plt
from .vertical_profile import extract_vertical_profile

ROOT=Path(__file__).resolve().parents[2]
PROJECT=ROOT/'PZT_Atomiser_Analysis'
FINAL=PROJECT/'Final_Results/Batch4_Cam1'
SOURCES={
 'water_density': 'https://www.nist.gov/system/files/documents/2017/05/09/NISTIR_7383_2013-04-24_2015-07-13Rev.pdf',
 'air_density': 'https://physics.nist.gov/cgi-bin/Star/compos.pl?matno=104&refer=ap',
 'air_viscosity': 'https://ntrs.nasa.gov/api/citations/20250005707/downloads/main.pdf?attachment=true',
 'drag_and_response': 'https://doc.comsol.com/6.4/doc/com.comsol.help.particle/particle_ug_fluid_flow.08.36.html'}


def fluid_properties(config):
    """20 C water density (NIST GLP10 polynomial), dry-air reference and Sutherland mu."""
    t=float(config['temperature_C'])
    if t!=20: raise ValueError('Re-source density properties before using another temperature')
    # NIST GLP10 air-free water, Tanaka equation at atmospheric pressure.
    water=999.974950*(1-((t-3.983035)**2*(t+301.797))/(522528.9*(t+69.34881)))
    # Rounded dry-air reference used as an explicitly configured property approximation.
    air=1.20479*(config['pressure_Pa']/101325.)
    if config['relative_humidity'] not in (None,0): raise ValueError('Humid-air properties require a supplied property model')
    mu=1.716e-5*((t+273.15)/273.15)**1.5*(273.15+110.4)/(t+273.15+110.4)
    return dict(rho_water_kg_m3=water,rho_air_kg_m3=air,mu_air_Pa_s=mu,g_m_s2=9.80665,
                air_density_note='NIST dry-air reference 1.20479 kg/m3; approximate 20 C, pressure scaled; not measured density',
                mu_air_note='Sutherland: mu0=1.716e-5 Pa s, T0=273.15 K, S=110.4 K',sources=SOURCES)


def relaxation_time(d_m,properties):
    d=np.asarray(d_m,float)
    if np.any(d<=0): raise ValueError('Diameter must be positive, in metres')
    return properties['rho_water_kg_m3']*d*d/(18*properties['mu_air_Pa_s'])


def reynolds(d_m,droplet_velocity,gas_velocity,properties):
    """Two-component in-plane slip Re; no unavailable out-of-plane component invented."""
    return properties['rho_air_kg_m3']*np.asarray(d_m)*np.linalg.norm(
        np.asarray(gas_velocity)-np.asarray(droplet_velocity),axis=-1)/properties['mu_air_Pa_s']


def drag_factor(re,config):
    """C_D Re/24, with exact Stokes branch and finite zero-slip limit."""
    re=np.asarray(re,float)
    if np.any(re<0) or np.any(re>=config['schiller_naumann_re_limit_exclusive']):
        raise ValueError('Outside configured spherical drag range 0<=Re<800')
    return np.where(re<=config['stokes_drag_re_threshold'],1.,1.+.15*re**.687)


def acceleration(velocity,d_m,gas_velocity,properties,config,drag=True,buoyancy=True):
    slip=np.asarray(gas_velocity,float)-np.asarray(velocity,float)
    a=np.zeros_like(slip)
    a[...,1]=-properties['g_m_s2']*(1-properties['rho_air_kg_m3']/properties['rho_water_kg_m3'] if buoyancy else 1)
    if drag:
        re=reynolds(d_m,velocity,gas_velocity,properties)
        a=a+slip*drag_factor(re,config)/relaxation_time(d_m,properties)
    return a


def upward_trajectory(d_m,initial,top_m,gas,properties,config,dt=None):
    """RK4 in time; bracket/bisect first turning or top-domain event.

    initial=(x_m,y_m,u_m_s,v_m_s). gas(position_m) returns (u_air,v_air).
    Returns ascending branch only. No empirical-profile interpolation is used.
    """
    state=np.asarray(initial,float).copy()
    if state[3]<=0: raise ValueError('Upward reference needs a positive measured initial v')
    step=min(dt or config['trajectory_time_step_s'],float(relaxation_time(d_m,properties))/30)
    def rhs(s): return np.r_[s[2:],acceleration(s[2:],d_m,gas(s[:2]),properties,config)]
    def advance(s,h):
        a=rhs(s);b=rhs(s+h*a/2);c=rhs(s+h*b/2);d=rhs(s+h*c)
        return s+h*(a+2*b+2*c+d)/6
    rows=[np.r_[0.,state]];time=0.;status='time_limit'
    while time<config['trajectory_max_time_s']:
        h=min(step,config['trajectory_max_time_s']-time);new=advance(state,h)
        events=[]
        for name,index,target in [('turning_point',3,0.),('domain_top',1,top_m)]:
            crossed=(new[index]<=target if index==3 else new[index]>=target)
            if crossed:
                lo,hi=0.,h
                for _ in range(42):
                    mid=(lo+hi)/2;s=advance(state,mid)
                    if (s[index]<=target if index==3 else s[index]>=target): hi=mid
                    else:lo=mid
                events.append((hi,name,index,target))
        if events:
            h,status,index,target=min(events);new=advance(state,h);new[index]=target
        time+=h;rows.append(np.r_[time,new]);state=new
        if events:break
    result=np.asarray(rows)
    if np.any(np.diff(result[:,2])<-1e-12) or np.any(result[:,4]<0):raise AssertionError('Descending branch')
    return result,status


def write_csv(path,columns):
    keys=list(columns);arrays=[np.asarray(columns[k]) for k in keys]
    if len({len(a) for a in arrays})!=1:raise ValueError('Unaligned CSV columns')
    with path.open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(keys);writer.writerows(zip(*arrays))


def figure(out,name,xlabel,ylabel):
    fig,ax=plt.subplots(figsize=(6.4,4.5),layout='constrained');ax.set(xlabel=xlabel,ylabel=ylabel)
    return fig,ax


def export(fig,ax,out,name):
    if name=='velocity_theory_comparison':
        ax.legend(frameon=False,fontsize=8,loc='upper center',bbox_to_anchor=(.5,-.17),ncol=2)
    else:ax.legend(frameon=False,fontsize=8)
    for ext in ('pdf','svg','png'):fig.savefig(out/(name+'.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)


def load_gas(config,profile):
    """Optional simultaneous gas data must be native-y aligned; never auto-register."""
    path=config.get('measured_gas_profile_csv')
    if not path:return None
    if not config.get('measured_gas_profile_registration'):raise ValueError('Document simultaneous gas registration')
    z=np.genfromtxt(ROOT/path,names=True,delimiter=',')
    np.testing.assert_allclose(z['y_mm'],profile['y_mm'],atol=1e-10,rtol=0)
    a=np.c_[z['u_air_m_s'],z['v_air_m_s']]
    if not np.isfinite(a).all():raise ValueError('Registered gas profile contains gaps')
    return a


def validate_physics(p,c):
    d=20e-6;zero=np.zeros(2);g=p['g_m_s2']*(1-p['rho_air_kg_m3']/p['rho_water_kg_m3'])
    np.testing.assert_allclose(acceleration(zero,d,zero,p,c),[0,-g])
    up=acceleration([0,1],d,zero,p,c);down=acceleration([0,-1],d,zero,p,c)
    assert up[1]<-g and down[1]>-g
    tiny=np.array([1e-9,2e-9]);a=acceleration(tiny,d,zero,p,c)-np.array([0,-g])
    np.testing.assert_allclose(a,-tiny/relaxation_time(d,p),rtol=1e-6,atol=1e-14)
    np.testing.assert_allclose(relaxation_time(2*d,p),4*relaxation_time(d,p))
    np.testing.assert_allclose(reynolds(d,[3,4],[0,0],p),p['rho_air_kg_m3']*d*5/p['mu_air_Pa_s'])
    try:drag_factor(800,c)
    except ValueError:pass
    else:raise AssertionError('Drag validity not enforced')
    tiny_d=5e-6;v0=.001;tau=float(relaxation_time(tiny_d,p))
    trajectory,status=upward_trajectory(tiny_d,[0,0,0,v0],1.,lambda q:zero,p,c)
    analytic_time=tau*np.log1p(v0/(g*tau))
    analytic_height=tau*v0-g*tau**2*np.log1p(v0/(g*tau))
    assert status=='turning_point'
    np.testing.assert_allclose(trajectory[-1,[0,2]],[analytic_time,analytic_height],rtol=1e-6,atol=1e-14)
    return dict(zero_slip_finite=True,gravity_downward=True,drag_opposes_slip=True,
                low_Re_Stokes_limit=True,analytic_Stokes_turning_time_height=True,tau_scales_as_d_squared=True,vector_norm_3_4_5=True,Re_validity_enforced=True)


def run(config_path=None):
    from .theory_package import run as build_package
    return build_package(config_path)


def notebook_cells():
    from .theory_package import notebook_cells as make_cells
    return make_cells()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path)
    print(json.dumps(run(parser.parse_args().config),indent=2))
