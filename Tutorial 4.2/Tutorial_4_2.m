% Tutorial 4.2
% Email: jlindbloom@smu.edu

clear all
close all

% Define parameters.
global g_leak g_na_max g_k_max g_t_max e_na e_k e_ca e_leak c_membrane

g_leak = 10e-9;
g_na_max = 3.6e-6;
g_k_max = 1.6e-6;
g_t_max = 0.22e-6;
e_na = 55e-3;
e_k = -90e-3;
e_ca = 120e-3;
e_leak = -70e-3;
c_membrane = 100e-12;


% Setup time vector.
global dt
dt = 0.01e-3;     
t = 0:dt:0.75;


% Create matrices for storing data.
numspikes = zeros(21,21);
interspikes = zeros(21,21);
for j=1:21
    for n=1:21
        [I_app] = applied_current(-200e-12 + (n-1)*20e-12, (j-1)*5e-12, t);
        [v_sim, h_sim, n_sim, h_t_sim, spikes] = PIR(t, I_app);
        nspikes = sum(spikes);
        results = find(spikes);
        min_spike_time = 10.0;
        for m=1:length(results)-1
            if (results(m+1)-results(m))*dt < min_spike_time
                min_spike_time = (results(m+1)-results(m))*dt;
            end
        end
        if min_spike_time == 10.0
            min_spike_time = 0;
        end

        numspikes(22-j,n) = nspikes;
        interspikes(22-j,n) = min_spike_time;  
    end
end

numspikes
interspikes

f1 = figure;
figure(f1);
imagesc(numspikes);
c = colorbar;
c.Label.String = 'Total Number of Spikes';
xlabel("Baseline Current (pA)");
ylabel("Current Step (pA)");
xticks([1:2:21]);
xticklabels({'-200','-160','-120','-80','-40','0','40','80','120','160','200'});
yticks([1:2:21]);
yticklabels({'100','90','80','70','60','50','40','30','20','10','0'});
title("Total Number of Spikes");
saveas(f1,"numspikes.png");

f2 = figure;
figure(f2);
imagesc(interspikes);
c = colorbar;
c.Label.String = 'Minimum ISI (seconds)';
xlabel("Baseline Current (pA)");
ylabel("Current Step (pA)");
xticks([1:2:21]);
xticklabels({'-200','-160','-120','-80','-40','0','40','80','120','160','200'});
yticks([1:2:21]);
yticklabels({'100','90','80','70','60','50','40','30','20','10','0'});
title("Minimum ISIs");
saveas(f2,"interspikes.png");


% Plot qualitatively distinct behaviors.

% Plot A: baseline=-150, step=15.
[I_app] = applied_current(-150e-12, 15e-12, t);
[v_sim, h_sim, n_sim, h_t_sim, spikes] = PIR(t, I_app);
nspikes = sum(spikes);
results = find(spikes);
min_spike_time = 10.0;
for m=1:length(results)-1
    if (results(m+1)-results(m))*dt < min_spike_time
        min_spike_time = (results(m+1)-results(m))*dt;
    end
end
if min_spike_time == 10.0
    min_spike_time = 0;
end

fA = figure;
figure(fA);
subplot(2,1,1);
plot(t, I_app);
xlabel("Time (seconds)");
ylabel("Applied Current");
title("Applied Current vs. Time");
subplot(2,1,2);
plot(t, v_sim);
xlabel("Time (seconds)")
ylabel("Membrane Potential");
title(sprintf("Plot A. Spikes = %d, Min ISI = %f",nspikes,min_spike_time));
saveas(fA,"PlotA.png");

% Plot B: baseline=-80, step=90.
[I_app] = applied_current(-80e-12, 90e-12, t);
[v_sim, h_sim, n_sim, h_t_sim, spikes] = PIR(t, I_app);
nspikes = sum(spikes);
results = find(spikes);
min_spike_time = 10.0;
for m=1:length(results)-1
    if (results(m+1)-results(m))*dt < min_spike_time
        min_spike_time = (results(m+1)-results(m))*dt;
    end
end
if min_spike_time == 10.0
    min_spike_time = 0;
end

fB = figure;
figure(fB);
subplot(2,1,1);
plot(t, I_app);
xlabel("Time (seconds)");
ylabel("Applied Current");
title("Applied Current vs. Time");
subplot(2,1,2);
plot(t, v_sim);
xlabel("Time (seconds)")
ylabel("Membrane Potential");
title(sprintf("Plot B. Spikes = %d, Min ISI = %f",nspikes,min_spike_time));
saveas(fB,"PlotB.png");

% Plot C: baseline=-20, step=45.
[I_app] = applied_current(-20e-12, 45e-12, t);
[v_sim, h_sim, n_sim, h_t_sim, spikes] = PIR(t, I_app);
nspikes = sum(spikes);
results = find(spikes);
min_spike_time = 10.0;
for m=1:length(results)-1
    if (results(m+1)-results(m))*dt < min_spike_time
        min_spike_time = (results(m+1)-results(m))*dt;
    end
end
if min_spike_time == 10.0
    min_spike_time = 0;
end

fC = figure;
figure(fC);
subplot(2,1,1);
plot(t, I_app);
xlabel("Time (seconds)");
ylabel("Applied Current");
title("Applied Current vs. Time");
subplot(2,1,2);
plot(t, v_sim);
xlabel("Time (seconds)")
ylabel("Membrane Potential");
title(sprintf("Plot C. Spikes = %d, Min ISI = %f",nspikes,min_spike_time));
saveas(fC,"PlotC.png");

% Plot D: baseline=100, step=45.
[I_app] = applied_current(100e-12, 45e-12, t);
[v_sim, h_sim, n_sim, h_t_sim, spikes] = PIR(t, I_app);
nspikes = sum(spikes);
results = find(spikes);
min_spike_time = 10.0;
for m=1:length(results)-1
    if (results(m+1)-results(m))*dt < min_spike_time
        min_spike_time = (results(m+1)-results(m))*dt;
    end
end
if min_spike_time == 10.0
    min_spike_time = 0;
end

fD = figure;
figure(fD);
subplot(2,1,1);
plot(t, I_app);
xlabel("Time (seconds)");
ylabel("Applied Current");
title("Applied Current vs. Time");
subplot(2,1,2);
plot(t, v_sim);
xlabel("Time (seconds)")
ylabel("Membrane Potential");
title(sprintf("Plot D. Spikes = %d, Min ISI = %f",nspikes,min_spike_time));
saveas(fD,"PlotD.png");


% Plot E: baseline=180, step=45.
[I_app] = applied_current(180e-12, 45e-12, t);
[v_sim, h_sim, n_sim, h_t_sim, spikes] = PIR(t, I_app);
nspikes = sum(spikes);
results = find(spikes);
min_spike_time = 10.0;
for m=1:length(results)-1
    if (results(m+1)-results(m))*dt < min_spike_time
        min_spike_time = (results(m+1)-results(m))*dt;
    end
end
if min_spike_time == 10.0
    min_spike_time = 0;
end

fE = figure;
figure(fE);
subplot(2,1,1);
plot(t, I_app);
xlabel("Time (seconds)");
ylabel("Applied Current");
title("Applied Current vs. Time");
subplot(2,1,2);
plot(t, v_sim);
xlabel("Time (seconds)")
ylabel("Membrane Potential");
title(sprintf("Plot E. Spikes = %d, Min ISI = %f",nspikes,min_spike_time));
saveas(fE,"PlotE.png");




%%%%%%%%%%%%%%%%%%%%%%%
% Function Definitions:

function [v_sim, h_sim, n_sim, h_t_sim, spikes] = PIR(t, i_applied, v_init, h_init, n_init, h_t_init)
% Simulates the thalamocortical neuron model with a T-type calcium current given the input time vector and
% applied current.

global dt g_leak g_na_max g_k_max g_t_max e_na e_k e_ca e_leak c_membrane

% Default parameters if not inputted.
if (~exist('v_init'))
    v_init = e_leak;
end
if (~exist('h_init'))
    h_init = 0;
end
if (~exist('n_init'))
    n_init = 0;
end
if (~exist('h_t_init'))
    h_t_init = 0;
end

% Setup vectors.
v_sim = zeros(1, length(t));
h_sim = zeros(1, length(t));
n_sim = zeros(1, length(t));
h_t_sim = zeros(1, length(t));
spikes = zeros(1, length(t));

v_sim(1) = v_init;
h_sim(1) = h_init;
n_sim(1) = n_init;
h_t_sim(1) = h_t_init;


% To count the number of spikes, a spike will be recorded when the membrane
% potential exceeds v_exceeds. The variable blocking will be set to 1,
% which will prevent more spikes from being recorded until the membrane
% potential falls below v_unblock, upon which the variable blocking will
% be set back to 0.

blocking = 0;
v_exceeds = 0.0;
v_unblock = -0.06;

% March forward in time.
for n = 1:(length(t)-1)
    if v_sim(n) > v_exceeds
        if blocking == 0
            spikes(n) = 1;
            blocking = 1;
        end
    end
    
    if v_sim(n) < v_unblock
        blocking = 0;
    end
    
    % Update v_sim.
    alpha = ((10^5)*(v_sim(n) + 0.035))/(1 - exp(-100*(v_sim(n)+0.035)));
    beta = 4000*exp((-(v_sim(n)+0.06))/(0.018));
    m = alpha/(alpha+beta);
    m_t = 1/(1 + exp((-(v_sim(n)+0.052))/(0.0074)));
    term1 = g_leak*(e_leak-v_sim(n));
    term2 = g_na_max*(m^3)*h_sim(n)*(e_na-v_sim(n));
    term3 = g_k_max*(n_sim(n)^4)*(e_k-v_sim(n));
    term4 = g_t_max*(m_t^2)*h_t_sim(n)*(e_ca-v_sim(n));
    v_sim(n+1) = v_sim(n) + (dt/c_membrane)*(term1+term2+term3+term4+i_applied(n));
    
    % Update h_sim.
    alpha = 350*exp(-50*(v_sim(n)+0.058));
    beta = 5000/(1 + exp(-100*(v_sim(n)+0.028)));
    term1 = alpha*(1-h_sim(n));
    term2 = -beta*h_sim(n);
    h_sim(n+1) = h_sim(n) + dt*(term1+term2);
    
    % Update n_sim.
    alpha = (( 5*(10^4) )*(v_sim(n)+0.034))/(1 - exp(-100*(v_sim(n)+0.034)));
    beta = 625*exp(-12.5*(v_sim(n)+0.044));
    term1 = alpha*(1-n_sim(n));
    term2 = -beta*n_sim(n);
    n_sim(n+1) = n_sim(n) + dt*(term1+term2);
    
    % Update h_t_sim.
    h_t_inf = 1/(1 + exp(500*(v_sim(n)+0.076)));
    if v_sim(n) < -0.080
        tau_h_t = 0.001*exp(15*(v_sim(n)+0.467));
    else
        tau_h_t = 0.028 + 0.001*exp((-(v_sim(n)+0.022))/(0.0105));
    end
    h_t_sim(n+1) = h_t_sim(n) + dt*((h_t_inf-h_t_sim(n))/tau_h_t);
end
end

function [I_applied] = applied_current(baseline, step, t)
% Returns a vector for applied current given input baseline and step
% currents.
global dt
    I_applied = zeros(1, length(t));
    third = floor(length(t)/3);
    twothird = 2*third;
    I_applied(1:third) = baseline;
    I_applied(third+1:twothird) = baseline+step;
    I_applied(twothird+1:length(I_applied)) = baseline;
end

