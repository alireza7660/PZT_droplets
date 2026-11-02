"""Two-component measured gas field and one-way droplet transport in SI units."""
import numpy as np
from .droplet_dynamics import acceleration,relaxation_time,reynolds,drag_factor


class UnsupportedPosition(Exception):
    """Outside measured gas support or in a cell touching an invalid vector."""


class DragValidityLimit(Exception):
    pass


def rotation(degrees):
    angle=np.deg2rad(float(degrees));c,s=np.cos(angle),np.sin(angle)
    matrix=np.array([[c,-s],[s,c]])
    # Exact cardinal rotations avoid false boundary rejections from cos(pi/2).
    matrix[np.abs(matrix)<1e-15]=0.
    return matrix


class MeasuredGasField:
    """Bilinear gas evaluation only inside four-finite-vector native cells.

    Local coordinates are metres relative to the ZNMF nozzle centre, not the
    saved velocity-peak reference. Velocities enter already in m/s. The same
    orthogonal rotation acts on positions and vectors, translation on positions.
    """
    def __init__(self,x_m,y_m,u,v,nozzle_position_m,rotation_degrees):
        self.x=np.asarray(x_m,float);self.y=np.asarray(y_m,float)
        self.u=np.asarray(u,float);self.v=np.asarray(v,float)
        if self.u.shape!=self.v.shape or self.u.shape!=(len(self.y),len(self.x)):
            raise ValueError('Gas field shape/axes disagree')
        if np.any(np.diff(self.x)<=0) or np.any(np.diff(self.y)<=0):
            raise ValueError('Gas coordinates must be strictly increasing')
        self.R=rotation(rotation_degrees);self.origin=np.asarray(nozzle_position_m,float)
        if self.origin.shape!=(2,) or not np.isfinite(self.origin).all():raise ValueError('Explicit finite nozzle position required')
        self.valid=np.isfinite(self.u)&np.isfinite(self.v)
        self.cell_size=float(min(np.min(np.diff(self.x)),np.min(np.diff(self.y))))

    def to_common(self,local_positions):return np.asarray(local_positions)@self.R.T+self.origin
    def vectors_to_common(self,local_vectors):return np.asarray(local_vectors)@self.R.T
    def to_local(self,common_positions):return (np.asarray(common_positions)-self.origin)@self.R

    def sample(self,position):
        x,y=self.to_local(position)
        if not (self.x[0]<=x<=self.x[-1] and self.y[0]<=y<=self.y[-1]):
            raise UnsupportedPosition('outside_measured_jet_domain')
        j=min(max(int(np.searchsorted(self.x,x,side='right'))-1,0),len(self.x)-2)
        i=min(max(int(np.searchsorted(self.y,y,side='right'))-1,0),len(self.y)-2)
        if not self.valid[i:i+2,j:j+2].all():raise UnsupportedPosition('cell_touches_invalid_gas_vector')
        a=(x-self.x[j])/(self.x[j+1]-self.x[j]);b=(y-self.y[i])/(self.y[i+1]-self.y[i])
        weight=np.array([[(1-a)*(1-b),a*(1-b)],[(1-a)*b,a*b]])
        uv=np.array([np.sum(weight*self.u[i:i+2,j:j+2]),np.sum(weight*self.v[i:i+2,j:j+2])])
        return self.vectors_to_common(uv)

    def common_mesh(self):
        xx,yy=np.meshgrid(self.x,self.y);q=self.to_common(np.stack([xx,yy],axis=-1))
        w=self.vectors_to_common(np.stack([self.u,self.v],axis=-1))
        return q[...,0],q[...,1],w[...,0],w[...,1]


def integrate_droplet(field,initial_state,diameter_m,properties,config,step_scale=1.):
    """RK4 with evolving vector slip and strict measured-domain termination.

    No upward-only stop is imposed: side-jet deflection and later settling may
    occur. The PZT field supplies only initial state, never a subsequent force
    or prescribed particle velocity. If a stage leaves support, the last
    supported substep is located conservatively; no gas extrapolation is used.
    """
    state=np.asarray(initial_state,float).copy();field.sample(state[:2])
    momentum=bool(config.get('save_momentum_integrals',False))
    if momentum:state=np.r_[state,np.zeros(4)]
    gravity=np.asarray(config.get('gravity_vector_m_s2',[0.,-properties['g_m_s2']]))*(1-properties['rho_air_kg_m3']/properties['rho_water_kg_m3'])
    controls=config['solver'];tau=float(relaxation_time(diameter_m,properties))
    time=0.;rows=[np.r_[time,state]];termination='maximum_time'
    def rhs(s):
        gas=field.sample(s[:2]);re=float(reynolds(diameter_m,s[2:4],gas,properties))
        if re>=config['schiller_naumann_re_limit_exclusive']:raise DragValidityLimit('Re_exceeds_drag_validity')
        net=acceleration(s[2:4],diameter_m,gas,properties,config)
        if momentum:return np.r_[s[2:4],net,net-gravity,gravity]
        return np.r_[s[2:4],net]
    def step(s,h):
        a=rhs(s);b=rhs(s+h*a/2);c=rhs(s+h*b/2);d=rhs(s+h*c)
        result=s+h*(a+2*b+2*c+d)/6
        rhs(result)
        # Check the traversed segment too: do not jump across a masked cell.
        for fraction in (.25,.5,.75):field.sample(s[:2]+fraction*(result[:2]-s[:2]))
        return result
    while time<controls['maximum_time_s']:
        gas=field.sample(state[:2])
        speed=max(np.linalg.norm(state[2:4]),np.linalg.norm(gas),1e-3)
        h=min(controls['time_step_s']*step_scale,tau/30*step_scale,
              controls['cell_fraction_per_step']*field.cell_size/speed,
              controls['maximum_time_s']-time)
        try:new=step(state,h)
        except UnsupportedPosition as problem:
            low,high=0.,h;best=state
            for _ in range(40):
                mid=(low+high)/2
                if mid<controls['minimum_time_step_s']:break
                try:trial=step(state,mid)
                except UnsupportedPosition:high=mid
                else:low=mid;best=trial
            if low>=controls['minimum_time_step_s']:
                time+=low;state=best;rows.append(np.r_[time,state])
            termination=str(problem);break
        except DragValidityLimit:
            termination='Re_exceeds_drag_validity';break
        time+=h;state=new;rows.append(np.r_[time,state])
    return np.asarray(rows),termination


def trajectory_quantities(field,trajectory,diameter_m,properties,D0_m,config):
    gas=np.array([field.sample(q) for q in trajectory[:,1:3]])
    drop=trajectory[:,3:5];slip=gas-drop;re=reynolds(diameter_m,drop,gas,properties)
    tau=float(relaxation_time(diameter_m,properties));gas_norm=np.linalg.norm(gas,axis=1)
    flow_time=np.divide(D0_m,gas_norm,out=np.full(len(gas_norm),np.inf),where=gas_norm>0)
    return dict(t_s=trajectory[:,0],x_mm=trajectory[:,1]*1e3,y_mm=trajectory[:,2]*1e3,
        droplet_u_m_s=drop[:,0],droplet_v_m_s=drop[:,1],gas_u_m_s=gas[:,0],gas_v_m_s=gas[:,1],
        slip_u_m_s=slip[:,0],slip_v_m_s=slip[:,1],slip_magnitude_m_s=np.linalg.norm(slip,axis=1),Re_in_plane=re,
        tau_stokes_s=np.full(len(re),tau),tau_drag_coefficient_s=tau/drag_factor(re,config),
        tau_flow_local_D0_over_gas_magnitude_s=flow_time,St_D_local=tau*gas_norm/D0_m)
