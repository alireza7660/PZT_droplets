# Figure captions

Coordinates are Cartesian: $x$ is horizontal/rightward and $y$ vertical/upward; $u$ is transverse and $v$ axial. The in-plane vector is $\mathbf{w}=(u,v)$; primes denote fluctuations about the respective means.

**mean_v**: Mean axial velocity $\overline{v}$. All finite signed axial measurements are displayed, including negative values. Only missing or invalid data are masked. Blue denotes negative velocity, white zero and red positive velocity; limits follow the actual range without symmetric endpoints.

**mean_u**: Mean transverse velocity $\overline{u}$; signed measurements are shown.

**rms_v**: Axial fluctuation RMS $v_{\mathrm{rms}}$, using denominator $N$.

**rms_u**: Transverse fluctuation RMS $u_{\mathrm{rms}}$, using denominator $N$.

**magnitude_mean_vector**: Magnitude of the mean in-plane velocity vector, $|\overline{\mathbf{w}}|=\sqrt{\overline{u}^{2}+\overline{v}^{2}}$.

**mean_instantaneous_velocity_magnitude**: Mean instantaneous in-plane velocity magnitude, $\overline{|\mathbf{w}|}=\overline{\sqrt{u^2+v^2}}$.

**covariance**: Signed cross-covariance $\overline{u\prime v\prime}$, using denominator $N$.

**correlation**: Signed correlation $\rho_{uv}$; zero-variance locations are masked.

**frame_1_u**: Numeric field key 1: instantaneous transverse velocity $u$. Shared limits are calculated from the three selected frames.

**frame_1_v**: Numeric field key 1: instantaneous axial velocity $v$. Shared limits are calculated from the three selected frames. All finite signed axial measurements are displayed, including negative values. Only missing or invalid data are masked. Blue denotes negative velocity, white zero and red positive velocity; limits follow the actual range without symmetric endpoints.

**frame_1_velocity_magnitude**: Numeric field key 1: instantaneous in-plane velocity magnitude. Shared limits are calculated from the three selected frames.

**frame_1_quiver**: Instantaneous signed velocity vectors for numeric field key 1; every fourth native node is displayed, with a common arrow scale.

**frame_879_u**: Numeric field key 879: instantaneous transverse velocity $u$. Shared limits are calculated from the three selected frames.

**frame_879_v**: Numeric field key 879: instantaneous axial velocity $v$. Shared limits are calculated from the three selected frames. All finite signed axial measurements are displayed, including negative values. Only missing or invalid data are masked. Blue denotes negative velocity, white zero and red positive velocity; limits follow the actual range without symmetric endpoints.

**frame_879_velocity_magnitude**: Numeric field key 879: instantaneous in-plane velocity magnitude. Shared limits are calculated from the three selected frames.

**frame_879_quiver**: Instantaneous signed velocity vectors for numeric field key 879; every fourth native node is displayed, with a common arrow scale.

**frame_1757_u**: Numeric field key 1757: instantaneous transverse velocity $u$. Shared limits are calculated from the three selected frames.

**frame_1757_v**: Numeric field key 1757: instantaneous axial velocity $v$. Shared limits are calculated from the three selected frames. All finite signed axial measurements are displayed, including negative values. Only missing or invalid data are masked. Blue denotes negative velocity, white zero and red positive velocity; limits follow the actual range without symmetric endpoints.

**frame_1757_velocity_magnitude**: Numeric field key 1757: instantaneous in-plane velocity magnitude. Shared limits are calculated from the three selected frames.

**frame_1757_quiver**: Instantaneous signed velocity vectors for numeric field key 1757; every fourth native node is displayed, with a common arrow scale.

**sampling_stations**: Configured PDF/CDF nodes on the mean axial-velocity field; profiles use rows 287, 159 and 32. All finite signed axial measurements are displayed, including negative values. Only missing or invalid data are masked. Blue denotes negative velocity, white zero and red positive velocity; limits follow the actual range without symmetric endpoints.

**finite_pair_count**: Paired finite sample count, without flag selection.

**finite_pair_fraction**: Finite-pair fraction; not a physical validation measure.

**flag_0_count**: Code 0 count; flag meanings are not assigned.

**flag_0_fraction**: Code 0 frequency.

**flag_0_mean_u**: Flag 0 population sensitivity: Mean transverse velocity $\overline{u}$; signed measurements are shown. Diagnostic limits are independent of ensemble limits.

**flag_0_mean_v**: Flag 0 population sensitivity: Mean axial velocity $\overline{v}$. Diagnostic limits are independent of ensemble limits. All finite signed axial measurements are displayed, including negative values. Only missing or invalid data are masked. Blue denotes negative velocity, white zero and red positive velocity; limits follow the actual range without symmetric endpoints.

**flag_0_rms_u**: Flag 0 population sensitivity: Transverse fluctuation RMS $u_{\mathrm{rms}}$, using denominator $N$. Diagnostic limits are independent of ensemble limits.

**flag_0_rms_v**: Flag 0 population sensitivity: Axial fluctuation RMS $v_{\mathrm{rms}}$, using denominator $N$. Diagnostic limits are independent of ensemble limits.

**flag_1_count**: Code 1 count; flag meanings are not assigned.

**flag_1_fraction**: Code 1 frequency.

**flag_1_mean_u**: Flag 1 population sensitivity: Mean transverse velocity $\overline{u}$; signed measurements are shown. Diagnostic limits are independent of ensemble limits.

**flag_1_mean_v**: Flag 1 population sensitivity: Mean axial velocity $\overline{v}$. Diagnostic limits are independent of ensemble limits. All finite signed axial measurements are displayed, including negative values. Only missing or invalid data are masked. Blue denotes negative velocity, white zero and red positive velocity; limits follow the actual range without symmetric endpoints.

**flag_1_rms_u**: Flag 1 population sensitivity: Transverse fluctuation RMS $u_{\mathrm{rms}}$, using denominator $N$. Diagnostic limits are independent of ensemble limits.

**flag_1_rms_v**: Flag 1 population sensitivity: Axial fluctuation RMS $v_{\mathrm{rms}}$, using denominator $N$. Diagnostic limits are independent of ensemble limits.

**profile_mean_v**: Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.

**profile_mean_u**: Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.

**profile_magnitude_mean_vector**: Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.

**profile_mean_instantaneous_velocity_magnitude**: Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.

**profile_rms_u**: Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.

**profile_rms_v**: Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.

**profile_covariance**: Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.

**profile_correlation**: Native transverse profiles at the existing three axial rows; signed measurements are retained, without interpolation or fitting.

**u_pdf**: Separate native-node full-support empirical histogram PDFs; no spatial pooling. All signed component samples are retained, including negative axial velocities.

**u_cdf**: Separate native-node empirical CDFs; no spatial pooling. All signed component samples are retained, including negative axial velocities.

**v_pdf**: Separate native-node full-support empirical histogram PDFs; no spatial pooling. All signed component samples are retained, including negative axial velocities.

**v_cdf**: Separate native-node empirical CDFs; no spatial pooling. All signed component samples are retained, including negative axial velocities.

**velocity_magnitude_pdf**: Separate native-node full-support empirical histogram PDFs; no spatial pooling. All signed component samples are retained, including negative axial velocities.

**velocity_magnitude_cdf**: Separate native-node empirical CDFs; no spatial pooling. All signed component samples are retained, including negative axial velocities.

**u_normalized_pdf**: Separate native-node full-support empirical histogram PDFs; no spatial pooling. All signed component samples are retained, including negative axial velocities.

**v_normalized_pdf**: Separate native-node full-support empirical histogram PDFs; no spatial pooling. All signed component samples are retained, including negative axial velocities.
