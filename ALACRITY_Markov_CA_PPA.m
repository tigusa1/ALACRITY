function ALACRITY_Markov_CA_PPA
%----------------------------------------------------------------------------------------------
% Simple Markov model to produce the CA and PPA results of the EAGLES trials
%----------------------------------------------------------------------------------------------
% Parameters
%   j = person (1~N)
%   i = week (1~24)
%   Prob[person j will abstain in week i] = probs(j) is constant (for each person j)
%   probs(j) is shown in the right of figure 100
%   modelNumber = 1 or 2 depending on the model for probs(j)
%               = 1 for a polynomial model (probs(j) is a polynomial function of j)
%               = 2 for a step-wise model  (probs(j) is a step function of j)
%   S(i,j) = state in week i, person j
%          = 1 if abstain or 0 if smoke
%----------------------------------------------------------------------------------------------
nt   = 24;                             % number of weeks
N    = 1026;                           % number of patients

S0   = zeros(nt,N);                    % states, steady-state probabilities
S    = S0;                             % non-steady-state

modelNumber = 2;                       % choose the type of model for the probabilities
models = {'polynomial','step'};        % model names
model  = models{modelNumber};          % model

switch model
    case 'polynomial'
        p0   = [0.00 0.99];            % parameters for polynomial model
    case 'step'
        p0   = [0.022 1.00             % Prob[abstain] = 100% for  4% of population
                0.12  0.945            % Prob[abstain] =  97% for  8% of population
                1.00  0.045];          % Prob[abstain] =   5% for 88% of population
end

[ probs,rhos ] = probA((1:N)/N,p0,model);  % probs = piA, rhos = pS/pA for each patient

rate = 0.4;                            % rate towards equilibrium, pS = rho*pA, pA = piA*rate

for j=1:N                              % for each patient
    for i=1:nt                         % for each week
        p  = probs(j);                 % get the probability of A
        pA = p*rate;
        pS = rhos(j)*pA;
        S0(i,j) = rand(1,1)<p;         % indicator for A (1=abstain, 0=smoke during week i)
        if i==1
            S(i,j) = 0;                % initial condition
        else
            S1 = S(i-1,j);               % previous condition
            if S1                        % abstain
                S(i,j) = rand(1,1)<1-pS; %   Prob[abstain|abstain] = 1-pS
            else                         % smoking
                S(i,j) = rand(1,1)<pA;   %   Prob[abstain|smoking] = pA
            end
        end
    end
end

CA12 = all(S0(9:12,:));                % indicator for CA (weeks 9~12)
CA24 = all(S0(9:24,:));                % indicator for CA (weeks 9~24)
PPA  = mean(S0,2);                     % weekly PPA
ca12 = all(S( 9:12,:));                % indicator for CA (weeks 9~12)
ca24 = all(S( 9:24,:));                % indicator for CA (weeks 9~24)
ppa  = mean(S, 2);                     % weekly PPA

t    = 1:nt;                           % time (weeks)

%----------------------------------------------------------------------------------------------
% Compare with placebo result (PPA = 15%, CA12 = 11%, CA24 = 8%)
%----------------------------------------------------------------------------------------------
fig = figure(100); fig.Name = 'PPA and prob'; fig.Color = 'w';
subplot(1,3,[1 2])
plot(t,100*PPA,'bo-',t,100*ppa,'ro-'), xlabel('week'), ylabel('7-day PPA')
title(sprintf('PPA versus time, CA (weeks 9~12) = %.1f%%, %.1f%%, CA (weeks 9~24) = %.1f%%, %.1f%%',...
    mean(CA12)*100,mean(ca12)*100,mean(CA24)*100,mean(ca24)*100))
ax = gca; ax.YLim = [0 50];

subplot(1,3,3)
area((1:N)/N*100,probs*100)
xlabel('proportion of patients (%)'), ylabel('probability'), title('Probability of 1-week abstinance')
ax = gca; ax.YLim = [0 100];

fig = figure(110); fig.Name = 'PPA'; fig.Color = 'w';
plot(t,100*ppa,'ro-')
ax = gca; ax.YLim = [0 50]; ax.XLim = [0 24]; ax.YTickLabel = {}; ax.XTickLabel = {};

figure(99), imagesc(S');

fig = figure(120); fig.Name = 'individuals'; fig.Color = 'w';
iks = [90 92 94 100 120 123 124];           % indices of individuals to plot
for i=1:length(iks)                      % for each individual
    ik  = iks(i);                        % get the index
    Sik = S(:,ik);                       % get the states over time
    for j=1:nt                           % for each time interval
        if Sik(j)==1
            plot([j-1 j],[1 1]*i,'b-','Linewidth',6), hold on
        end
    end
end
ax = gca; ax.XLim = [0 nt]; ax.XTick = [0:4:24]; hold off
ax.FontSize = 16; xlabel('Week'); ylabel('simulated individual ID'), grid on
title(sprintf('Sample of microsimulation results\n'))

keyboard


function [ p,rho ] = probA( x, param, type )
%----------------------------------------------------------------------------------------------
% Probability of abstaining in one week
%----------------------------------------------------------------------------------------------

switch type
    case 'polynomial'
        p = param(2)*x.^4 + param(1)*(1-x.^4);
    case 'step'
        for i=1:length(x)
            j    = find(x(i)<=param(:,1),1,'first'); % last column of param(1,:) less than x(i)
            p(i) = param(j,2);                       % get Prob[abstain]
        end
end

rho = (1-p)./p;                                      % rho = pS / pA

