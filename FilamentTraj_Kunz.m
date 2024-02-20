%Simulates growth of a filament in a limited pool of monomers.
% clear all;
% close all;

r1=0.5; %rate of growth
g1=10; %rate of decay
Ntot=1000;%total number of available monomers
MaxS=200000; %max steps
times = [2, 5, 10, 20, 60, 100];
MaxTraj=100; %number of trajectories
p1= zeros(1, Ntot); %probability
probfig = figure;
for k = 1:length(times)
    figure;
    for j=1:MaxTraj
    
        m1=nan(1,MaxS);
        m1(1)=1;
        monomers=Ntot;
        T=zeros(1,MaxS);
        T(1)=0;
    
        for i=1:MaxS
        
            k1=r1*(monomers-m1(1));
            if m1(i) == 1
                k2=0;
            else
                k2=g1;
            end
        
            k0=k1+k2;
        
            % Determine time spent
            CoinFlip1=rand;
            tau(i)=(1/k0)*log(1/CoinFlip1); %also, tau(i)=exprnd(1/k0);
        
            T(i+1)= T(i)+tau(i);
            % Determine reaction
            CoinFlip2=rand;
            if CoinFlip2<=k1/k0 || m1(i) == 1
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

            if T(i+1) >= times(k)
                break;
            end
        end
        plot(T,m1,'.-', 'MarkerSize',10)
        title("Filament lengths")
        xlabel('time')
        ylabel('filament length')
        hold on;
   % xlim([0 100])
        ylim([0 Ntot+500])
        p1(m1(i+1))= p1(m1(i+1))+1; %calculating probability of a specific length
    end

    p1 = p1/sum(p1);
    x = 1:1:Ntot;
    Avg= sum(x.*p1);
    variance= sum((x.^2).*p1)-(sum(x.*p1))^2;

    figure(probfig);
    dispstring = sprintf('Probability at time %d', times(k));
    plot(x,p1, "DisplayName", dispstring);
    xlabel('Filament length')
    ylabel('Probability')
    title('Probability vs length')
    xlim([0 Ntot])
    legend('Location', 'northwest');
    hold on
    p1= zeros(1, Ntot);
end
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
