"""Spherical-droplet forces in SI units; gas velocity must be supplied separately."""
import numpy as np

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
    gravity=np.asarray(config.get('gravity_vector_m_s2',[0.,-properties['g_m_s2']]),float)
    if gravity.shape!=(2,):raise ValueError('Gravity vector must have two components in the integration frame')
    a=np.broadcast_to(gravity,slip.shape).copy()*(1-properties['rho_air_kg_m3']/properties['rho_water_kg_m3'] if buoyancy else 1)
    if drag:
        re=reynolds(d_m,velocity,gas_velocity,properties)
        coefficient=drag_factor(re,config)/relaxation_time(d_m,properties)
        a=a+slip*np.asarray(coefficient)[...,None]
    return a


