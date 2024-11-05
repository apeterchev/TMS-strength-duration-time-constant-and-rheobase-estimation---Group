%% MEMBCONSTECTALL: Estimate strength-duration (neural membrane) time
%% constant (tau_m) and rheobase (mt_b) from MT data and cTMS pulse
%% waveforms--assuming distinct individual values of mt_b but a common
%% (joint) value of tau_m
% 
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

AUTHOR:  Angel V. Peterchev (c) 2005-2024
% VERSION:  08/23/2024
%

%% global variables used in subroutines

global wvfrm_all
global subjid_mt
   

%% load cTMS electric field waveform data

load ctms1_wvfrm_11_26_2008.mat;    % variables: ctms1_wvfrm
wvfrm_all = ctms1_wvfrm;
    

%% load MT related data

% subjid_mt - [subject id mt_pw1 mt_pw2 mt_pw3 ... ]
% pw - [pw1 pw2 pw3 ...]
% V_MSO - (V/%MSO) scaling factor of pulse amplitude from % maximum stimuator output to device capacitor voltage in volts; this changes the units of the MT to initial coil voltage (optional)  
load subjid_mt30_mt60_mt120;        % variables: subjid_mt30_mt60_mt120, pw, V 
subjid_mt = subjid_mt30_mt60_mt120;
 

%% MT descriptive statistics
n_pw = length(pw);             % number of PWs
n_subj = size(subjid_mt,1);    % number of subjects

mt_mean = mean(subjid_mt(:,2:n_pw+1));
mt_std = std(subjid_mt(:,2:n_pw+1));
mt_max = max(subjid_mt(:,2:n_pw+1));
mt_min = min(subjid_mt(:,2:n_pw+1));
%mt_skew = skewness(subjid_mt(:,2:n_pw+1));
%mt_kurt = kurtosis(subjid_mt(:,2:n_pw+1))-3;


%% Specify parameters ranges constraining the optimization space
% This can help with global convergence, but it has to be ensured the range is not
% overly constraining

% time constant constraints:
tau_m0 = 200e-6; % optimization starting point (initial guess)
tau_m_lb = tau_m0/100; % optimization lower bound
tau_m_ub = tau_m0*100; % optimization upper bound

% rheobase constraints:
mt_b0 = min(mt_mean)/2; % optimization starting point (initial guess)
mt_b_lb = mt_b0/100; % optimization lower bound
mt_b_ub = mt_b0*100; % optimization upper bound


%% Parameter estimation

[param,resnorm,residual,exitflag,output,lambda,jacobian] = lsqcurvefit(@mt_mtcalcall,[tau_m0 mt_b0.*ones(1,n_subj)],pw,ones(n_subj*n_pw,1),...
    [tau_m_lb mt_b_lb.*ones(1,n_subj)],[tau_m_ub mt_b_ub.*ones(1,n_subj)]);
ci = nlparci(param,residual,'jacobian',jacobian);  % tau_m confidence interval
tau_m = param(1);
mt_b = param(2:end);


%% Calculate depolarization (r) factors for optimal tau_m (optional)

r = zeros(1,n_pw);           % membrane depolarization factor
fs = wvfrm_all.fs;
b = 1/(1+2*tau_m*fs).*[1 1];
a = [1 (1-2*tau_m*fs)/(1+2*tau_m*fs)];
for i = 1:n_pw
    i_wvfrm = find(wvfrm_all.pw == pw(i));
    r(i) = max(filter(b,a,wvfrm_all.wvfrm(:,i_wvfrm)));
end


%% Rheobase descriptive statistics (optional)

mt_b_mean = mean(mt_b);
mt_b_std = std(mt_b);
mt_b_min = min(mt_b);
mt_b_max = max(mt_b); 
%mt_b_skew = skewness(mt_b);
%mt_b_kurt = kurtosis(mt_b)-3;


%% Printing & plotting of results (optional)

% results in %MSO
disp(' ')
disp(['   ' sprintf('  %s\t','subj#',string(cellstr(strcat('MT',num2str(pw')))'),'mt_b (%)','tau_m(us)')])
disp([subjid_mt(:,:) mt_b(:) ones(n_subj,1).*tau_m*1e6])
disp(' ')
%disp('parameter       mean    std     min     max ')
disp(sprintf('%s \t\t%3.1f','tau_m (us)',tau_m*1e6));
disp(sprintf('%s \t\t%3.1f \t%3.1f \t%3.1f \t%3.1f \t%3.2f \t%3.2f','mt_b (%)',mt_b_mean,mt_b_std,mt_b_min,mt_b_max));
for i = 1:n_pw
   disp(sprintf('mt[%3dus](%%) \t%3.1f \t%3.1f \t%3.1f \t%3.1f \t%3.2f \t%3.2f',pw(i),mt_mean(i),mt_std(i),mt_min(i),mt_max(i))); 
end    
disp(' ')
disp('----------------------------------------------')
disp(' ')

% results in volts (V)
disp(' ')
disp(sprintf('\t\t%s','subj#',string(cellstr(strcat('MT',num2str(pw')))'),'mt_b(V)','tau_m(us)'))
disp([subjid_mt(:,1) round(V_MSO.*subjid_mt(:,(2:(n_pw+1)))) round(V_MSO.*mt_b(:)),ones(n_subj,1).*round(tau_m*1e6)])
disp(' ')
disp('parameter       mean    std     min     max')
disp(sprintf('%s \t\t%4d','tau_m (us)',round(tau_m*1e6)));
disp(sprintf('%s \t\t%4d \t%4d \t%4d \t%4d \t%3.2f \t%3.2f','mt_b (%)',...
    round(V_MSO.*mt_b_mean),round(V_MSO.*mt_b_std),round(V_MSO.*mt_b_min),round(V_MSO.*mt_b_max)));
for i = 1:n_pw
   disp(sprintf('mt[%3dus](%%) \t%4d \t%4d \t%4d \t%4d \t%3.2f \t%3.2f',pw(i),round(...
       V_MSO.*mt_mean(i)),round(V_MSO.*mt_std(i)),round(V_MSO.*mt_min(i)),round(V_MSO.*mt_max(i)))); 
end  
disp(' ')
disp('----------------------------------------------')
disp(' ')
disp(['N = ' num2str(n_subj)]);
disp(['PW = ' num2str(pw)]);
disp(['r = ' num2str(r)]);
disp(['tau_m 95% CI (us) = ' num2str(ci(1,:)*1e6)]);
disp(' ')
disp('----------------------------------------------')
disp(' ')

% Plot MTs
figure(1); clf
plot(pw,subjid_mt(:,2:n_pw+1))
title('MTs');
xlabel('PW (\mus)');
ylabel('Pulse amplitude (%MSO)');

% Plot E-field pulses scaled by MT for all PWs
figure(2); clf
hold on
clr = ['b' 'r' 'g' 'm' 'c' 'y'];
for i = 1:n_pw
    i_wvfrm = find(wvfrm_all.pw == pw(i));
    plot(wvfrm_all.t.*1e6,mt_mean(i).*wvfrm_all.wvfrm(:,i_wvfrm),clr(i),'LineWidth',2);
end
%title('Pulse Shapes for Average MTs');
xlabel('Time (\mus)');
ylabel('Electric field (normalized to MT)');
legend(cellstr(strcat(num2str(pw'),' \mus'))')
hold off

%  Rheobase histogram
figure(3); clf
hist(mt_b)
title('Human Motor Cortex Motor Threshold Rheobase')
xlabel('MT rheobase (%MSO)')
ylabel('# of occurrences')

%% Neuron membrane leaky integrate-and-fire model function

function mt_mt = mt_mtcalcall(param,pw)
    
    % Set to global and load wvfrm_all and subjid_mt before using this function.
    % param(1) = membrane time constant (tau_m)
    % param(2:n) = depolarization threshold voltage relative to pulse amplitude (%)
    %           (rheobase)
    
    global wvfrm_all
    global subjid_mt
    
    tau_m = param(1);                   % membrane time const. (same for all subjects)
    mt_b = param(2:length(param));      % rheobases (%MSO)
    
    n_pw = length(pw);                  % number of pulse widths
    n_subj = size(subjid_mt,1);    % number of subjects
    
    r = zeros(1,n_pw);           % membrane depolarization factor
    mt_mt = zeros(n_subj,n_pw);         % mt(estimated)/mt(measured)
    
    % calculate depolarization factor r
    fs = wvfrm_all.fs;
    b = 1/(1+2*tau_m*fs).*[1 1];
    a = [1 (1-2*tau_m*fs)/(1+2*tau_m*fs)];
    for i = 1:n_pw
        i_wvfrm = find(wvfrm_all.pw == pw(i));
        r(i) = max(filter(b,a,wvfrm_all.wvfrm(:,i_wvfrm)));
    end
    
    mt_mt = mt_b'*ones(1,n_pw);
    mt_mt = mt_mt./r;
    mt_mt = mt_mt./subjid_mt(:,2:(1+n_pw));
    mt_mt = mt_mt';
    mt_mt = mt_mt(:);

end
