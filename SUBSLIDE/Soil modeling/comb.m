%Simple DP Mob
function [] = comb()
K = 10000;
G = 10000;
N0 = -0.2;
M = 1.5;
Mrho0 = 0.5;
lambda_c = 20;
lambda_rho = 200;
sig0 = [150; 0];
D = [K 0;0 3*G];
Sigma(:,1) = sig0;
Epsilon(:,1) = [0;0];
Nk(1,1) = N0;
Mrhok(1,1) = Mrho0;
%strain controlled
deps = [0.00; 0.004];
%mean stress controlled
p = 100;
%**********
N=N0;
sig = sig0;Mrho = Mrho0;
for k = 2:100
    x = deps(1);
    options = optimoptions('fsolve','Display','none');
    deps(1) = fsolve(@mob_one,x,options);
    sig0 = sig;
    N0 = N;
    Mrho0 = Mrho;
    Epsilon(:,k) = Epsilon(:,k-1)+deps;
    Sigma(:,k) = sig;
    Nk(1,k) = N;
    Mrhok(1,k) = Mrho;
end
figure (1)    
subplot(2,2,1); hold on;
plot(Sigma(1,:),Sigma(2,:)); axis equal;
xlabel('$p$','interpreter','latex')
ylabel('$q$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,2); hold on;
plot(Epsilon(2,:),Sigma(2,:))
xlabel('${\varepsilon_q}$','interpreter','latex')
ylabel('$q$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,3); hold on; axis ij;
plot(Sigma(1,:),Epsilon(1,:))
xlabel('$p$','interpreter','latex')
ylabel('${\varepsilon_p}$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,4); hold on; axis ij
plot(Epsilon(2,:),Epsilon(1,:))
xlabel('${\varepsilon_q}$','interpreter','latex')
ylabel('${\varepsilon_p}$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')


figure (2)
plot(Epsilon(2,:),Mrhok(1,:))
xlabel('${\varepsilon_q}$','interpreter','latex')
ylabel('$M_{\rho}$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
figure (3)
plot(Epsilon(2,:),Nk(1,:))
xlabel('${\varepsilon_q}$','interpreter','latex')
ylabel('$N$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
function [r] = mob_one(x)
    deps(1) = x;
%Elastic trial
    lambda = 0;
    sig = sig0 + D * deps;
    N = N0;
    Mrho = Mrho0;
    y = abs(sig(2)) - (Mrho+N)*sig(1);
    if y > 0
        i = 0;
        while abs(y) > 1d-6
            i = i + 1
            dlambda = y/(sig(1)*(lambda_rho*(M-Mrho)-lambda_c*(M+N)*N)+K*N*(Mrho+N)+3*G);
            lambda = lambda + dlambda;
            depsp(1,1) = -lambda*N;
            depsp(2,1) = lambda*sign(sig(2));
            depse = deps - depsp;
            sig = sig0 + D * depse;
            N = N0 + lambda_c*(M+N)*depsp(1,1);
            Mrho = Mrho0 + lambda_rho*(M-Mrho)*abs(depsp(2,1));
            if sig(1) < 0 %Liquefaction cut-off on p
                sig(1) = 0d0;
                sig(2) = 0d0;
            end
            y = abs(sig(2)) - (Mrho+N)*sig(1);
         end
    end
    % strain controlled
    r = 0;
    %mean stress controlled
    %r = p-sig(1)+1/3*sig(2);
end
end

