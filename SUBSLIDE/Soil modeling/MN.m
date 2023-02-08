% Simple DP Mob
function [] = MN()
E = 10000;                                      % Young's modulus
nuy = 0.3;                                      % Poisson's ratio
lambda_c = 20;                                  % Dilatancy parameters
M = 2 * sqrt(2) * tan (45/180*pi());            % Friction parameters
Epsilon_v_cs = 0.2;                             % Critical state volume strain
N0 = M*(1-exp(-lambda_c*(Epsilon_v_cs-0)));     % Dilatanct state variables

% Intial stress
Sigma = [150 0 0; 0 150 0;0 0 150];        
% 3 x 3 stress tensor

% Elastic matrix
D = E/(1+nuy)/(1-2*nuy)*[ 1-nuy nuy nuy 0 0 0
                        ; nuy 1-nuy nuy 0 0 0
                        ; nuy nuy 1-nuy 0 0 0
                        ; 0 0 0 (1-2*nuy)/2 0 0
                        ; 0 0 0 0 (1-2*nuy)/2 0
                        ; 0 0 0 0 0 (1-2*nuy)/2];
% Strain controled
Epsilon = zeros(6,1);
deps = [0.001; -0.0005; -0.0005; 0; 0; 0];


% Stress variant
I1 = trace(Sigma);
I2 = 0.5*(trace(Sigma)^2 - trace(Sigma^2);
I3 = det (Sigma);
X = sqrt (I1*I2/I3 - 9);

%Elastic trial
Nk(1,1) = N0;
N=N0;
sig0 =  [Sigma(1,1); Sigma(2,2); Sigma(3,3); Sigma(1,2); Sigma(1,3); Sigma(2,3)];
sig = sig0 + D * deps;      % Stress vector
y = X - (M+N);              % Yield function

if y > 0
        i = 0;
        while abs(y) > 1d-6
            i = i + 1
            dlambda = 0;
            lambda = lambda + dlambda;
            
            % Flow rule
            depsp(1) = -lambda*N;
            depsp(2) = lambda*sign(sig(2));
            depsp(3)
            depsp(4)
            depsp(5)
            depsp(6)
            
            % Update stress
            depse = deps - depsp;
            sig = sig0 + D * depse;
            
            % Update state varialbes
            N = N0 + lambda_c*(M+N)*depsp(1,1);
            
            
            if sig(1) < 0 %Liquefaction cut-off on p
                sig(1) = 0d0;
                sig(2) = 0d0;
            end
            
            y = X - (M+N);
        end
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
plot(Epsilon(2,:),Nk(1,:))
xlabel('${\varepsilon_q}$','interpreter','latex')
ylabel('$N$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
