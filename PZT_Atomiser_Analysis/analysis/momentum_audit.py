"""Single-droplet force, impulse and initial-condition coverage audits in SI units."""
import copy
import numpy as np
from .combined_flow_prediction import write_csv,write_json
from .droplet_dynamics import acceleration,relaxation_time,reynolds,drag_factor
from .prescribed_gas import integrate_droplet


def cumulative_trapezoid(t,values):
    values=np.asarray(values);return np.vstack([np.zeros(values.shape[1]),np.cumsum(.5*(values[:-1]+values[1:])*np.diff(t)[:,None],axis=0)])


def add_force_and_impulse_fields(q,tr,d,props,cfg,net):
    mass=props['rho_water_kg_m3']*np.pi*d**3/6;area=np.pi*d**2/4
    gravity=np.array(cfg['gravity_vector_m_s2'])*(1-props['rho_air_kg_m3']/props['rho_water_kg_m3']);drag=net-gravity
    slip=np.c_[q['slip_u_p_m_s'],q['slip_v_p_m_s']];norm=np.linalg.norm(slip,axis=1);re=q['Re_p'];factor=drag_factor(re,cfg)
    cd=np.divide(24*factor,re,out=np.full_like(re,np.nan),where=re>0)
    classic=np.zeros_like(drag);nz=re>0
    classic[nz]=.5*props['rho_air_kg_m3']*area/mass*cd[nz,None]*norm[nz,None]*slip[nz]
    np.testing.assert_allclose(classic,drag,rtol=1e-12,atol=1e-10)
    assert tr.shape[1]==9,'RK4 passive drag/gravity impulse states required'
    force_drag=mass*drag;force_gravity=np.broadcast_to(mass*gravity,drag.shape);force_total=mass*net
    impulse_drag=mass*tr[:,5:7];impulse_gravity=mass*tr[:,7:9]
    trap=cumulative_trapezoid(q['t_s'],force_drag);trapg=cumulative_trapezoid(q['t_s'],force_gravity)
    dv=np.c_[q['droplet_u_p_m_s']-q['droplet_u_p_m_s'][0],q['droplet_v_p_m_s']-q['droplet_v_p_m_s'][0]];momentum=mass*dv
    for k,component in enumerate(('x_p','y_p')):
        q['drag_a_'+component+'_m_s2']=drag[:,k];q['effective_gravity_a_'+component+'_m_s2']=np.full(len(tr),gravity[k]);q['net_a_'+component+'_m_s2']=net[:,k]
        q['force_drag_'+component+'_N']=force_drag[:,k];q['force_effective_gravity_'+component+'_N']=force_gravity[:,k];q['force_total_'+component+'_N']=force_total[:,k]
        q['drag_impulse_'+component+'_N_s']=impulse_drag[:,k];q['effective_gravity_impulse_'+component+'_N_s']=impulse_gravity[:,k];q['momentum_change_'+component+'_kg_m_s']=momentum[:,k]
        q['momentum_residual_RK4_'+component+'_N_s']=momentum[:,k]-impulse_drag[:,k]-impulse_gravity[:,k]
        q['drag_impulse_trapezoid_'+component+'_N_s']=trap[:,k];q['gravity_impulse_trapezoid_'+component+'_N_s']=trapg[:,k]
        q['momentum_residual_trapezoid_'+component+'_N_s']=momentum[:,k]-trap[:,k]-trapg[:,k]
        q['force_balance_residual_'+component+'_m_s2']=net[:,k]-drag[:,k]-gravity[k]
    q['droplet_mass_kg']=np.full(len(tr),mass);q['projected_area_m2']=np.full(len(tr),area);q['drag_response_factor_F']=factor;q['C_D_undefined_at_zero_slip']=cd
    np.testing.assert_allclose(q['effective_gravity_impulse_y_p_N_s'],mass*gravity[1]*q['t_s'],atol=1e-25,rtol=1e-12)


def budget_record(h,upper=None):
    q=h['q'];n=h['peak']+1
    def value(key):return float(q[key][-1]) if upper is None else float(np.interp(upper,q['h_mm'][:n],q[key][:n]))
    mass=float(q['droplet_mass_kg'][0]);record=dict(release_number=h['number'],is_reference=h['reference'],scope='full_trajectory' if upper is None else 'common_supported_interval',entry_x_p_mm=float(q['x_p_mm'][0]),entry_h_mm=float(q['h_mm'][0]),end_h_mm=value('h_mm'),elapsed_time_s=value('t_s'),droplet_mass_kg=mass)
    for component in ('x_p','y_p'):
        for key in ('drag_impulse_'+component+'_N_s','effective_gravity_impulse_'+component+'_N_s','momentum_change_'+component+'_kg_m_s','momentum_residual_RK4_'+component+'_N_s','drag_impulse_trapezoid_'+component+'_N_s','momentum_residual_trapezoid_'+component+'_N_s'):record[key]=value(key)
        record['drag_equivalent_delta_v_'+component+'_m_s']=value('drag_impulse_'+component+'_N_s')/mass;record['gravity_equivalent_delta_v_'+component+'_m_s']=value('effective_gravity_impulse_'+component+'_N_s')/mass
        record['residual_RK4_equivalent_delta_v_'+component+'_m_s']=value('momentum_residual_RK4_'+component+'_N_s')/mass;record['residual_trapezoid_equivalent_delta_v_'+component+'_m_s']=value('momentum_residual_trapezoid_'+component+'_N_s')/mass
    return record


def write_momentum_tables(ctx,histories,lower,upper):
    out=ctx['out'];rows=[budget_record(h,limit) for limit in (None,upper) for h in histories]
    for r in rows:
        assert max(abs(r['residual_RK4_equivalent_delta_v_'+c+'_m_s']) for c in ('x_p','y_p'))<1e-10
        assert max(abs(r['residual_trapezoid_equivalent_delta_v_'+c+'_m_s']) for c in ('x_p','y_p'))<ctx['cfg']['momentum_audit']['independent_trapezoid_velocity_residual_tolerance_m_s']
    write_csv(out/'data/momentum_budget.csv',{k:[r[k] for r in rows] for k in rows[0]})
    endpoints=[]
    for h in histories:
        q=h['q']
        for index,stage in [(0,'entry'),(-1,'exit')]:
            endpoints.append(dict(release_number=h['number'],is_reference=h['reference'],stage=stage,t_s=float(q['t_s'][index]),x_p_mm=float(q['x_p_mm'][index]),h_mm=float(q['h_mm'][index]),y_p_mm=float(q['y_p_mm'][index]),droplet_u_p_m_s=float(q['droplet_u_p_m_s'][index]),droplet_v_p_m_s=float(q['droplet_v_p_m_s'][index]),gas_u_p_m_s=float(q['gas_u_p_m_s'][index]),gas_v_p_m_s=float(q['gas_v_p_m_s'][index]),slip_u_p_m_s=float(q['slip_u_p_m_s'][index]),slip_v_p_m_s=float(q['slip_v_p_m_s'][index]),Re_p=float(q['Re_p'][index]),tau_stokes_s=float(q['tau_stokes_s'][index]),tau_drag_coefficient_s=float(q['tau_drag_coefficient_s'][index]),termination_reason=h['termination'] if stage=='exit' else 'measured_entry_initial_condition'))
    write_csv(out/'data/entry_exit_conditions.csv',{k:[r[k] for r in endpoints] for k in endpoints[0]})
    reference=next(h for h in histories if h['reference']);q=reference['q'];n=reference['peak']+1;ys=q['h_mm'][:n];target=np.r_[ys[(ys>=lower)&(ys<upper)],upper]
    columns={k:np.interp(target,ys,v[:n]) for k,v in q.items()};write_csv(out/'data/reference_force_momentum_common_interval.csv',columns)
    validity=[]
    for h in histories:
        r=h['q']['Re_p'];validity.append(dict(release_number=h['number'],minimum_Re_p=float(r.min()),maximum_Re_p=float(r.max()),Stokes_branch_samples=int(np.sum(r<=ctx['cfg']['stokes_drag_re_threshold'])),SN_branch_samples=int(np.sum(r>ctx['cfg']['stokes_drag_re_threshold'])),every_RK_stage_checks_upper_Re_limit=True,Re_limit_exclusive=ctx['cfg']['schiller_naumann_re_limit_exclusive']))
    write_json(out/'validation/drag_validity.json',validity)
    return dict(reference_common=next(r for r in rows if r['scope']=='common_supported_interval' and r['is_reference']),maximum_full_path_RK4_equivalent_velocity_residual_m_s=max(float(np.nanmax(abs(h['q']['momentum_residual_RK4_'+c+'_N_s']/h['q']['droplet_mass_kg']))) for h in histories for c in ('x_p','y_p')),maximum_full_path_trapezoid_equivalent_velocity_residual_m_s=max(float(np.nanmax(abs(h['q']['momentum_residual_trapezoid_'+c+'_N_s']/h['q']['droplet_mass_kg']))) for h in histories for c in ('x_p','y_p')),quadrature_note='RK4-stage impulse closure is expected to roundoff because the same stages update velocity; independent trapezoidal quadrature of saved forces is also reported and checked. No residual is forced to zero.',validity_max_Re_p=max(r['maximum_Re_p'] for r in validity))


def release_coverage(ctx,states):
    from .aligned_codirected import pzt_sample
    a=ctx['a'];height=ctx['base']['y_p_mm'];corners=np.array([[a['x_p_mm'][i,j],a['y_p_mm'][i,j]] for i,j in [(0,0),(0,-1),(-1,-1),(-1,0)]])
    intersections=[]
    for p,r in zip(corners,np.roll(corners,-1,axis=0)):
        if min(p[1],r[1])-1e-12<=height<=max(p[1],r[1])+1e-12 and r[1]!=p[1]:intersections.append(p[0]+(height-p[1])*(r[0]-p[0])/(r[1]-p[1]))
    lo,hi=min(intersections),max(intersections);spacing=float(np.median(np.diff(a['x_p_mm'],axis=1)))/4
    xs=np.unique(np.r_[0,np.arange(np.ceil((lo+1e-8)/spacing),np.floor((hi-1e-8)/spacing)+1)*spacing]);values=[]
    for x in xs:
        try:v,stencil=pzt_sample(ctx,[x,height])
        except ValueError:v=[np.nan,np.nan]
        values.append(v)
    values=np.asarray(values);seeds=[float(state['initial'][0]*1e3) for state in states]
    write_csv(ctx['out']/'data/release_coverage_profile.csv',dict(x_p_mm=xs,y_p_mm=np.full(len(xs),height),h_mm=np.zeros(len(xs)),mean_u_p_m_s=values[:,0],mean_v_p_m_s=values[:,1]))
    return dict(entry_h_mm=0.,entry_original_y_p_mm=height,selection='Audit original seven projected native-grid offsets; retain only four-finite-corner PZT and gas support on the same measured base plane. No extrapolation or plane shift to retain candidates.',initial_x_p_mm=seeds,accepted_count=len(states),candidate_count=len(ctx['cfg']['co_directed']['release_column_offsets']),supported_native_plane_x_p_bounds_mm=[float(lo),float(hi)],whole_plume_distribution=False,measured_population_weights_available=False,entry_FWHM_mm=None,seed_span_fraction_of_FWHM=None,interpretation='The tilted native support intersects h=0 asymmetrically. Only supported starts are used. No whole-plume width or right-hand field is invented at this boundary plane.')


def synthetic_force_checks(cfg,props,d):
    config=copy.deepcopy(cfg);tau=float(relaxation_time(d,props));g=np.array(cfg['gravity_vector_m_s2'])*(1-props['rho_air_kg_m3']/props['rho_water_kg_m3'])
    velocity=np.array([.2,.4]);np.testing.assert_allclose(acceleration(velocity,d,velocity,props,cfg)-g,[0,0],atol=1e-15)
    assert float(reynolds(d,velocity,velocity,props))==0.
    for slip in ([0,.02],[0,-.02],[.01,-.03]):
        drag=acceleration(velocity,d,velocity+slip,props,cfg)-g;assert np.dot(drag,slip)>0
        np.testing.assert_allclose(drag,np.array(slip)*float(drag_factor(reynolds(d,velocity,velocity+slip,props),cfg))/tau,rtol=1e-13,atol=1e-13)
    weak_slip=np.array([0,.5*abs(g[1])*tau]);weak_net=acceleration([0,0],d,weak_slip,props,cfg);assert weak_net[1]<0 and weak_net[1]-g[1]>0
    class ConstantGas:
        cell_size=1.
        def sample(self,point):return np.array([.001,.002])
    gas=np.array([.001,.002]);initial=np.array([0,0,-.001,.001]);config['solver']['maximum_time_s']=5*tau;config['solver']['time_step_s']=tau/50;config['save_momentum_integrals']=True
    errors=[]
    for scale in (1.,.5):
        trajectory,reason=integrate_droplet(ConstantGas(),initial,d,props,config,scale);t=trajectory[:,0];terminal=gas+tau*g;deviation=initial[2:]-terminal;decay=np.exp(-t/tau)
        analytic_v=terminal+decay[:,None]*deviation;analytic_r=initial[:2]+t[:,None]*terminal+tau*(1-decay[:,None])*deviation
        verr=float(np.max(abs(trajectory[:,3:5]-analytic_v)));perr=float(np.max(abs(trajectory[:,1:3]-analytic_r)))
        re=reynolds(d,trajectory[:,3:5],gas,props);assert re.max()<cfg['stokes_drag_re_threshold'];assert verr<1e-9 and perr<1e-12
        errors.append(dict(step_scale=scale,max_position_error_m=perr,max_velocity_error_m_s=verr,max_Re_p=float(re.max()),termination=reason))
    return dict(zero_slip_drag_zero=True,Re_uses_vector_slip_not_particle_velocity=True,drag_parallel_to_slip=True,faster_particle_drag_decelerates_axially=True,faster_gas_drag_accelerates_axially=True,upward_drag_below_effective_weight_gives_downward_net_acceleration=True,analytical_Stokes_constant_gas_2D_position_and_velocity_pass=True,analytical_solution='v=ug+tau*g_eff+(v0-ug-tau*g_eff)*exp(-t/tau); r=r0+(ug+tau*g_eff)*t+tau*(v0-ug-tau*g_eff)*(1-exp(-t/tau))',analytical_checks=errors,existing_piecewise_drag_model_retained=True)


def export_placed_gas(ctx,field):
    """Export the exact native prescribed gas used by trajectory and contour evaluation."""
    x,y,u,v=field.common_mesh();xx,yy=np.meshgrid(field.x,field.y);rows,cols=np.indices(u.shape)
    data=dict(native_row=rows,native_column=cols,ZNMF_Xj_mm=xx*1e3,ZNMF_Yj_mm=yy*1e3,x_p_mm=x*1e3,h_mm=y*1e3,y_p_mm=y*1e3+ctx['base']['y_p_mm'],mean_gas_u_p_m_s=u,mean_gas_v_p_m_s=v,finite_pair=field.valid)
    write_csv(ctx['out']/'data/prescribed_gas_plume_frame.csv',{key:value.ravel() for key,value in data.items()})
    np.savez_compressed(ctx['out']/'data/prescribed_gas_plume_frame.npz',**data,rotation_matrix=field.R,translation_m=field.origin,scope='Prescribed analysis copy of independent measured ZNMF mean, not combined-flow measurements')
