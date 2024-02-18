%Simulates growth of a filament in a limited pool of monomers.
% clear all;
% close all;

r1=0.5; %rate of growth
g1=10; %rate of decay
Ntot=1000;%total number of available monomers
MaxS=200000; %max steps
t1=2500; %get out at this time
MaxTraj=1000; %number of trajectories
p1= zeros(1, MaxS); %probability
p2 = zeros(1, MaxS);
figure;
for j=1:MaxTraj
    
    m1=nan(1,MaxS);
    m2 = nan(1, MaxS);
    m1(1)=1;
    m2(1) = 1;
    monomers=Ntot;
    T=zeros(1,MaxS);
    
    for i=1:MaxS
        
        k1=r1*(monomers-m1(1));
        k2=g1;
        k3 = r1*(monomers-m2(1));
        k4=g1;
        
        k0=k1+k2+k3+k4;
        
        % Determine time spent
        CoinFlip1=rand;
        tau(i)=(1/k0)*log(1/CoinFlip1); %also, tau(i)=exprnd(1/k0);
        
        T(i+1)= T(i)+tau(i);
        % Determine reaction
        CoinFlip2=rand;
        if CoinFlip2 <= k1/k0
            m1(i+1) = m1(i) + 1;
            m2(i+1) = m2(i);
            monomers = monomers-1;
        elseif CoinFlip2 <= (k1+k2) / k0
            m1(i+1) = m1(i) - 1;
            m2(i+1) = m2(i);
            monomers = monomers+1;
        elseif CoinFlip2 <= (k1+k2+k3) / k0
            m2(i+1) = m2(i) + 1;
            m1(i+1) = m1(i);
            monomers = monomers-1;
        else
            m2(i+1) = m2(i) - 1;
            m1(i+1) = m1(i);
            monomers = monomers+1;
        end
     
        if T(i+1)>=t1
            break;
        end
    end
    plot(T,m1, 'g.-', T, m2, 'b.-', 'MarkerSize',10)
    title("Filament lengths")
    xlabel('time')
    ylabel('filament length')
    p1(m1(i+1))= p1(m1(i+1))+1; %calculating probability of a specific length
    p2(m2(i+1))= p2(m2(i+1))+1;
   hold on;
end
legend('Polymer 1', 'Polymer 2')

p1 = p1/sum(p1);
p2 = p2/sum(p2);
x = 1:1:MaxS;
Avg= sum(x.*p1);
variance= sum((x.^2).*p1)-(sum(x.*p1))^2;

% time_steady = 40;
alm = m1(400:i); % the actual length measurements in sampling intervals of lag_t(i) minutes
mean_y= mean(alm);
var_y= std(alm)^2;
% 
%analytical solution
T1= 0:0.1:2000;
Lss = (Ntot-g1/r1);
L = Lss*(1- exp(-r1.*T1));

% plot(T1,L,'k','LineWidth', 2)
% hold off;
% xlim([0, 100])
% title('Length trajectory')

%probability
figure;
plot(x,p1);
xlabel('Filament length')
ylabel('Probability')
title('Probability vs length 1')
xlim([0 Ntot])
hold on

figure;
plot(x,p2);
xlabel('Filament length')
ylabel('Probability')
title('Probability vs length 2')
xlim([0 Ntot])
hold on
