%Simulates growth of a filament in a limited pool of monomers.
% clear all;
% close all;

r1=0.5; %rate of growth
g1=10; %rate of decay
Ntot=1000;%total number of available monomers
MaxS=200000; %max steps
t1=100; %get out at this time
MaxTraj=1000; %number of trajectories
p1= zeros(1, Ntot); %probability
t20 = zeros(1, MaxS);
t60 = zeros(1, MaxS);
t100 = zeros(1, MaxS);
p2 = zeros(1, MaxS);
p5 = zeros(1, MaxS);
p10 = zeros(1, MaxS);
p20 = zeros(1, MaxS);
p60 = zeros(1, MaxS);
p100 = zeros(1, MaxS);
c20 = 1;
c60 = 1;
c100 = 1;
figure;
for j=1:MaxTraj
    
    m1=nan(1,MaxS);
    m1(1)=1;
    monomers=Ntot;
    T=zeros(1,MaxS);
    T(1)=0;
    
    for i=1:MaxS
        
        k1=r1*(monomers-m1(1));
        k2=g1;
        
        k0=k1+k2;
        
        % Determine time spent
        CoinFlip1=rand;
        tau(i)=(1/k0)*log(1/CoinFlip1); %also, tau(i)=exprnd(1/k0);
        
        T(i+1)= T(i)+tau(i);
        % Determine reaction
        CoinFlip2=rand;
        if CoinFlip2<=k1/k0
            m1(i+1)=m1(i)+1;
            monomers=monomers-1;
        else
            if m1(i)==1
               m1(i+1)= m1(i);
            else
               m1(i+1)=m1(i)-1;
               monomers=monomers+1;
            end
            
                
        end

        if round(T(i))==2
            p2(m1(i+1)) = p2(m1(i+1))+1;
        elseif round(T(i))==5
                p5(m1(i+1)) = p5(m1(i+1))+1;
        elseif round(T(i))==10
                p10(m1(i+1)) = p10(m1(i+1))+1;
        elseif round(T(i))==20
                t20(c20) = m1(i+1);
                p20(m1(i+1)) = p20(m1(i+1))+1;
                c20 = c20+1;
        elseif round(T(i))==60
                t60(c60) = m1(i);
                p60(m1(i+1)) = p60(m1(i+1))+1;
                c60 = c60+1;
        elseif round(T(i))==100
                t100(c100) = m1(i);
                p100(m1(i+1)) = p100(m1(i+1))+1;
                c100 = c100+1;
        end
     
        if T(i+1)>=t1
            break;
        end
    end
    plot(T,m1,'.-', 'MarkerSize',10)
    xlabel('time')
    ylabel('filament length')
   hold on;
   % xlim([0 100])
    ylim([0 Ntot+500])
    p1(m1(i+1))= p1(m1(i+1))+1; %calculating probability of a specific length
end

p1 = p1/sum(p1);
p2 = p2/sum(p2);
p5 = p5/sum(p5);
p10 = p10/sum(p10);
p20 = p20/sum(p20);
p60 = p60/sum(p60);
x = 1:1:Ntot;
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
newcolors = [1 0 0
             0 1 0
             0 0 1
             1 1 0
             1 0 1
             0 1 1];
colororder(newcolors)
plot(x,p1);
xlabel('Filament length')
ylabel('Probability')
title('Probability vs length')
xlim([0 Ntot])
hold on
plot(p2);
plot(p5);
plot(p10);
plot(p20);
plot(p60);
legend("Prob at t=100", "Prob at t=2", "Prob at t=5", "Prob at t=10", "Prob at t=20", "Prob at t=60")

figure;
t20 = t20(t20~=0);
histogram(t20, EdgeColor="r", FaceAlpha=0);
xlabel("Filament length at steady-state times");
ylabel("Occurrences");
title("# of filament lengths at t=20");
xlim([950 Ntot]);
hold on
t60 = t60(t60~=0);
histogram(t60, EdgeColor="b", FaceAlpha=0);
t100 = t100(t100~=0);
histogram(t100, EdgeColor="g", FaceAlpha=0);
legend("Filaments at t=20", "Filaments at t=60", "Filaments at t=100")

