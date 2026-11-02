"""Supported dense-state event search: axial slip and net acceleration are distinct."""
import numpy as np
from .droplet_dynamics import acceleration
from .combined_flow_prediction import write_csv,write_json


def events_from_q(ctx,field,q,upper):
    d=ctx['cfg']['droplet_diameters_um'][0]*1e-6;gravity=np.array(ctx['cfg']['gravity_vector_m_s2'])*(1-ctx['props']['rho_air_kg_m3']/ctx['props']['rho_water_kg_m3']);events=[]
    def describe(state,time,kind,i=-1):
        gas=field.sample(state[:2]);net=acceleration(state[2:],d,gas,ctx['props'],ctx['cfg']);slip=gas-state[2:];lab=np.array([state[0]*1e3,state[1]*1e3+ctx['base']['y_p_mm']])@ctx['Q']+ctx['origin']
        return dict(event=kind,lower_trajectory_index=i,t_s=float(time),x_p_mm=float(state[0]*1e3),h_mm=float(state[1]*1e3),y_p_mm=float(state[1]*1e3+ctx['base']['y_p_mm']),PZT_lab_frame_x_mm=float(lab[0]),PZT_lab_frame_y_mm=float(lab[1]),droplet_u_p_m_s=float(state[2]),droplet_v_p_m_s=float(state[3]),gas_u_p_m_s=float(gas[0]),gas_v_p_m_s=float(gas[1]),slip_u_p_m_s=float(slip[0]),slip_v_p_m_s=float(slip[1]),slip_magnitude_m_s=float(np.linalg.norm(slip)),drag_a_y_p_m_s2=float(net[1]-gravity[1]),net_a_y_p_m_s2=float(net[1]),within_common_interval=bool(0<=state[1]*1e3<=upper))
    state=np.c_[q['x_p_mm']*1e-3,q['h_mm']*1e-3,q['droplet_u_p_m_s'],q['droplet_v_p_m_s']];derivative=np.c_[state[:,2:],q['net_a_x_p_m_s2'],q['net_a_y_p_m_s2']]
    for key,kind in [('net_a_y_p_m_s2','net_acceleration_zero'),('slip_v_p_m_s','axial_slip_zero')]:
        series=q[key]
        for i in np.flatnonzero(series[:-1]*series[1:]<0):
            dt=q['t_s'][i+1]-q['t_s'][i]
            def dense(z):return (2*z**3-3*z*z+1)*state[i]+(z**3-2*z*z+z)*dt*derivative[i]+(-2*z**3+3*z*z)*state[i+1]+(z**3-z*z)*dt*derivative[i+1]
            left,right=0.,1.;fleft=series[i]
            for _ in range(45):
                mid=(left+right)/2;point=dense(mid);gas=field.sample(point[:2]);value=acceleration(point[2:],d,gas,ctx['props'],ctx['cfg'])[1] if key=='net_a_y_p_m_s2' else gas[1]-point[3]
                if fleft*value<=0:right=mid
                else:left=mid;fleft=value
            z=(left+right)/2;direction='negative_to_positive' if series[i]<0 else 'positive_to_negative';event=describe(dense(z),q['t_s'][i]+z*dt,kind+'_'+direction,int(i));events.append(event)
    minima=[e for e in events if e['event']=='net_acceleration_zero_negative_to_positive'];slips=[e for e in events if e['event']=='axial_slip_zero_negative_to_positive']
    first_min=min(minima,key=lambda e:e['t_s']) if minima else None;first_slip=min(slips,key=lambda e:e['t_s']) if slips else None
    common_nodes=np.flatnonzero(q['h_mm']<=upper);j=int(common_nodes[np.argmin(q['droplet_v_p_m_s'][common_nodes])]);node=describe(state[j],q['t_s'][j],'minimum_saved_node',j)
    options=[node]+[e for e in minima if e['within_common_interval']];global_min=min(options,key=lambda e:e['droplet_v_p_m_s'])
    return dict(first_velocity_minimum=first_min,minimum_over_common_interval=global_min,first_axial_zero_slip=first_slip,zero_vector_slip_identified=any(e['slip_magnitude_m_s']<1e-9 for e in events if e['event'].startswith('axial_slip_zero')),minimum_saved_vector_slip_m_s=float(np.min(np.hypot(q['slip_u_p_m_s'],q['slip_v_p_m_s']))),events=sorted(events,key=lambda e:e['t_s']))


def analyse_reference_events(ctx,field,histories,upper):
    h=next(h for h in histories if h['reference']);q=h['q'];result=events_from_q(ctx,field,q,upper)
    fine=h['fine_trajectory'];gas=np.array([field.sample(p) for p in fine[:,1:3]]);net=acceleration(fine[:,3:5],ctx['cfg']['droplet_diameters_um'][0]*1e-6,gas,ctx['props'],ctx['cfg'])
    fq=dict(t_s=fine[:,0],x_p_mm=fine[:,1]*1e3,h_mm=fine[:,2]*1e3,droplet_u_p_m_s=fine[:,3],droplet_v_p_m_s=fine[:,4],net_a_x_p_m_s2=net[:,0],net_a_y_p_m_s2=net[:,1],slip_u_p_m_s=gas[:,0]-fine[:,3],slip_v_p_m_s=gas[:,1]-fine[:,4])
    refined=events_from_q(ctx,field,fq,upper);checks={}
    for key in ('first_velocity_minimum','first_axial_zero_slip'):
        a,b=result[key],refined[key]
        if a is not None and b is not None:
            delta=dict(time_difference_s=abs(a['t_s']-b['t_s']),height_difference_mm=abs(a['h_mm']-b['h_mm']),velocity_difference_m_s=abs(a['droplet_v_p_m_s']-b['droplet_v_p_m_s']))
            assert delta['height_difference_mm']*1e-3<ctx['cfg']['solver']['position_convergence_tolerance_m'];assert delta['velocity_difference_m_s']<ctx['cfg']['solver']['velocity_convergence_tolerance_m_s'];checks[key]=delta
    result['half_step_event_checks']=checks;result['event_method']='Bracketed sign changes using cubic Hermite interpolation of supported RK4 states and endpoint derivatives; reevaluate native gas and drag, then bisect. Full vector slip and axial slip are reported separately.'
    if result['events']:write_csv(ctx['out']/'data/reference_events.csv',{k:[e[k] for e in result['events']] for k in result['events'][0]})
    write_json(ctx['out']/'validation/event_refinement.json',checks);write_json(ctx['out']/'data/reference_event_summary.json',result)
    np.savez_compressed(ctx['out']/'data/reference_half_step_state.npz',trajectory=fine)
    return result
