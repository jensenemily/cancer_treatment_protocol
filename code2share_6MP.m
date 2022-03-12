%% Simulation of 6-MP PK model from Jayachandran PlosOne 2014 paper

% PK parameters
kab = 4.2; % per day
kel = 3.8; % per day
kcm = 39.4; % pmol 6MP converted/ day
kme = 0.08; % per day
K = 15.11; % pmol
vcm = 1; % pmol 6TGN produced/ pmol 6MP/ 8x10^8 RBCs
b = 10^6/170.20; % pmol 6MP/mg 6MP

% PK dynamics
dose_vals = 100; % daily dose of 6-MP administered (mg/ m^2)
BSA = 1.78; % m^2 (body surface area)
cycle_length = 21; % days

% simulate full PK dynamics and compare to lower order approximation
for i = 1:numel(dose_vals)
    num_days = 14;
    d = @(t) dose_vals(i)*BSA*(t<=num_days) + 0*(t>num_days); % dose in mg/day 
    dpk = @(t,xpk) [-kab*xpk(1) + b*d(t);...
                    kab*xpk(1) - kel*xpk(2) - kcm*xpk(2)./(K+xpk(2));...
                    vcm*kcm*xpk(2)./(K+xpk(2)) - kme*xpk(3)]; % full dynamics for PK with states xpk = [x_g;x_c;x_t]

    [tpk, xpk] = ode45(dpk, [0,200], [0;0;0]); % simulate full dynamics from zero initial condition (zero drug to start)
    
    % reduced dynamics for PK 
    ac = .20*b; ag = .24*b; % approx scale factor between dose and xc, xg
    dxm_approx = @(t, xm) vcm*kcm*ac*d(t)./(K+ac*d(t)) - kme*xm; % state = approx of x_t
    [txm_approx, xm_approx] = ode45(dxm_approx, [0,cycle_length], 0);
    
    % plot:
    figure;
    subplot(3,1,1)
    plot(tpk, d(tpk), 'LineWidth', 2);
    set(gca,'fontsize',14)
    ylabel('dose (mg)', 'FontSize', 16);
    title('6-MP & 6-TGN dynamics for 50mg/m^2 dose', 'FontSize', 18)
    xlim([0 cycle_length])
    subplot(3,1,2)
    plot(tpk, xpk(:,3), 'LineWidth', 2); hold on; plot(txm_approx, xm_approx, '--', 'LineWidth', 2);
    set(gca,'fontsize',14)
    legend( 'x_t (pmol/10^8 RBC)', 'approx x_t (pmol/10^8 RBC)', 'FontSize', 12)
    ylabel('6-TGN (pmol/10^8 RBC)', 'FontSize', 16)
    xlim([0 cycle_length])
    subplot(3,1,3)
    plot(tpk, xpk(:,1), 'LineWidth', 2); hold on;
    set(gca,'fontsize',14)
    plot(tpk, xpk(:,2), 'LineWidth', 2);
    legend( 'x_g (pmol)', 'x_c (pmol)', 'FontSize', 12)
    ylabel('6-MP (pmol)', 'FontSize', 16)
    xlabel('Time (days)', 'FontSize', 16)
    xlim([0 cycle_length])
end

%% Simulation of Drug Effect & Approximations 

% PD parameters
mu = .3287; % per day
p = 8.2*10^9; % cells/ L blood
yl = 0.4368; % dimensionless
ELmax = 0.0782; % per day
ECL50 = 84; % pmol/ 8x10^8 RBCs
ktl = 0.1207; % per day
kdl = 0.5346; % per day

ELdrug = @(xt) ELmax*xt./(ECL50 + xt); % full drug effect model from Jayachandran PlosOne 2014 paper

% plot full drug effect model 
figure
plot(tpk, ELdrug(xpk(:,3)), 'LineWidth', 1.5)
set(gca,'fontsize',14)
hold on

% log linear/affine approximation of drug effect
plot(tpk, -.02 + (1/70)*log(1+xpk(:,3)), 'LineWidth', 1.5)
plot(tpk, (1/100)*log(1+xpk(:,3)), 'LineWidth', 1.5)

% linear approximation of drug effect
plot(tpk, 0 + .00019*xpk(:,3), 'LineWidth', 1.5)

legend('full model','log linear approx','log affine approx','linear approx')
xlabel('Time (days)', 'FontSize', 18)
ylabel('Drug effect (1/days)', 'FontSize', 18)
xlim([0 cycle_length])

%% Simulation of 6-MP PD model from Jayachandran PlosOne 2014 paper
% using full drug effect model 

k_fb = @(xl) mu*(p^yl./(p^yl + xl.^yl)); % feedback term
E_drug = @(xT) ELmax*xT./ (ECL50 + xT); % drug effect 
dpkpd = @(t, x) [-kab*x(1) + b*d(t);...
                    kab*x(1) - kel*x(2) - kcm*x(2)./(K+x(2));...
                    vcm*kcm*x(2)./(K+x(2)) - kme*x(3);...
                    (k_fb(x(5+3)) - E_drug(x(3)) -ktl).*x(1+3);...
                    ktl*(x(1+3)-x(2+3));...
                    ktl*(x(2+3)-x(3+3));...
                    ktl*(x(3+3) -x(4+3));...
                    ktl*x(4+3)- kdl*x(5+3)]; % full combined PK-PD dynamics
L0 = ((mu*p^yl)/ktl - p^yl)^(1/yl); % (baseline)
C0 = kdl*L0/ktl;          
pk_init = [0;0;0]; pd_init = [C0;C0;C0;C0;L0]; % initial conditions - zero drug, equilibrium blood meas.
[tpd, xpd] = ode45(dpkpd, [0,200], [pk_init; pd_init]);

% plot PK-PD dynamics
figure;
subplot(4,1,1)
plot(tpd, xpd(:,3), 'LineWidth', 2); set(gca,'fontsize',14)
legend( 'x_T', 'FontSize', 12)
ylabel({'TGN Concentration', '(pmol/10^8 RBC)'}, 'FontSize', 14)
subplot(4,1,2)
plot(tpd, xpd(:,1), 'LineWidth', 2); hold on;
plot(tpd, xpd(:,2), 'LineWidth', 2); set(gca,'fontsize',14)
legend( 'x_g (gut)', 'x_c (plasma)', 'FontSize', 12)
ylabel('6-MP (pmol)', 'FontSize', 14)
subplot(4,1,3)
plot(tpd, xpd(:,4), 'LineWidth', 2); hold on;
plot(tpd, xpd(:,5), 'LineWidth', 2);
plot(tpd, xpd(:,6), 'LineWidth', 2);
plot(tpd, xpd(:,7), 'LineWidth', 2); set(gca,'fontsize',14)
legend( 'x_s','x_1','x_2','x_3', 'FontSize', 12)
ylabel({'Transit Compartments', '(cells/L)'}, 'FontSize', 14)
subplot(4,1,4)
plot(tpd, xpd(:,8), 'LineWidth', 2); set(gca,'fontsize',14)
legend( 'x_l', 'FontSize', 12)
ylabel({'Mature WBCs', '(cells/L)'}, 'FontSize', 14)
xlabel('Time (days)', 'FontSize', 18)

%% Compare to reduced model dynamics 

% parameters for reduced model
slope = .00019; % 8x10^8 RBCs/(day x pmol 6TGN)
ac = .20*b; ag = .24*b; % approx proportionality constant to dose input (pmol x day/ mg)

dapprox_pkpd = @(t, x) [slope*vcm*kcm*ac*b*d(t)./(K+ ac*b*d(t)) - kme*x(1);...
                         mu*p^yl./( p^yl + (L0*exp(x(2))).^yl ) - x(1) - kdl];
[tapprox1, xapprox1] = ode45(dapprox_pkpd, [0,50], [0;0]);

% compare reduced drug effect model to full drug effect
figure; hold on
plot(tpk, ELdrug(xpk(:,3)), 'LineWidth', 1.5)
plot(tpk, 0 + .00019*xpk(:,3), 'LineWidth', 1.5)
plot(tapprox1, xapprox1(:,1), 'LineWidth', 1.5);
xlim([0 50])
set(gca,'fontsize',14)
xlabel('Time (days)', 'FontSize', 18)
ylabel('Drug Effect (1/days)', 'FontSize', 18)
legend( 'full drug effect model E_{drug}(x_T)', 'linear approx of E_{drug}(x_T)','x_M(t)')

% compare reduced WBC model to full (first using full drug effect model)
dpk_approxpd = @(t, x) [-kab*x(1) + b*d(t);...
                    kab*x(1) - kel*x(2) - kcm*x(2)./(K+x(2));...
                    vcm*kcm*x(2)./(K+x(2)) - kme*x(3);...
                    mu*p^yl./( p^yl + (L0*exp(x(2))).^yl ) - E_drug(x(3)) - kdl];
[tapprox2, xapprox2] = ode45(dpk_approxpd, [0,50], [0;0;0;0]);
figure
plot(tapprox2, exp(xapprox2(:,4))*L0);

% 2 state reduced model
approx_2state = @(t,x) [(vcm*kcm*ac*b*d(t))./(K+ac*b*d(t)) - kme*x(1);...
                        (k_fb(x(2)) - E_drug(x(1)) - ktl).* x(2)];
[tapprox3,xapprox3] = ode45(approx_2state, [0,200], [0;L0]);

% plot
figure;
subplot(2,1,1); plot(tapprox3,xapprox3(:,1), 'LineWidth', 1.5); hold on
plot(tpk, xpk(:,3), 'LineWidth', 1.5); set(gca,'fontsize',14)
ylabel({'6-TGN Concentration',' (pmol/ 10^8 RBCs)'}, 'FontSize', 18)
legend('Approximate ($\tilde{x}_T$)','Full Model ($x_T$)', 'Interpreter', 'latex')
subplot(2,1,2); plot(tapprox3,xapprox3(:,2), 'LineWidth', 1.5); hold on
plot(tpd, xpd(:,8), 'LineWidth', 2); set(gca,'fontsize',14)
ylabel({'Mature WBCs','(cells/L)'}, 'FontSize', 18)
xlabel('time (days)', 'FontSize', 18)
legend('Approximate ($\tilde{x}_{\ell}$)','Full Model ($x_{\ell}$)', 'Interpreter', 'latex')
%%






