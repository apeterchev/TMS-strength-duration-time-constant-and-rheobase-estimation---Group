# TMS-strength-duration-time-constant-and-rheobase-estimation---Group

FILE: membconstestall.m

OBJECTIVE:
The objective of this code is to estimate the strength-duration (SD) time
constant and rheobase from cTMS neural strength-duration curve data
across a population of subjects when assuming a single, common SD time 
constant across subjects and individual rheobase for each subject. 

BACKGROUND:
Peterchev et al. 2013 used two methods for estimation of the
strength-duration time constant with cTMS. The second approach jointly
estimated a single time constant across all subjects, with individual
rheobases for each subject. The rationale is that the rheobase will vary
because of all kinds if effects like the coil-to-cortex distance, but the
time constant is intrinsic to neurons.

INSTRUCTIONS:
Run this file (membconstestall.m) for an example joint time constant
analysis. 

The recorded E-field waveforms are in data array ctms1_wvfrm.wvfrm in
example file ctms1_wvfrm_11_26_2008.mat.

One needs to use recordings of the electric field (E-field) pulse
waveform used in the specific device and study. The waveform should then
be stored in MATLAB in the same way as the ctms1_wvfrm structure. The
E-field waveforms may have to be preprocessed to ensure two things:
1) the baseline before the pulse is zero mean, and
2) the peak amplitude of the E-field waveforms (excluding small transient
switching spikes) is the same across all waveforms and is normalized to
unity (1). 

The example is for 3 waveform (pulse width) conditions, but can be
extended to more than three conditions.

REFERENCES:

Peterchev AV, Goetz SM, Westin GG, Luber B, Lisanby SH. Pulse width
dependence of motor threshold and input-output curve characterized with
controllable pulse parameter transcranial magnetic stimulation. Clin
Neurophysiol. 2013 Jul;124(7):1364-72.
doi: https://doi.org/10.1016/j.clinph.2013.01.011

Menon, P., Pavey, N., Aberra, A. S., van den Bos, M. A. J., Wang, R.,
Kiernan, M. C., Peterchev, A. V.*, and Vucic, S.* Dependence of cortical
neuronal strength-duration properties on TMS pulse shape. Clin
Neurophysiol 2023; 150: 106-118 * equal contribution.
doi: https://doi.org/10.1016/j.clinph.2023.03.012

AUTHOR:  Angel V. Peterchev
