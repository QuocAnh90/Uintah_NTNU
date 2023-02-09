% Simple DP Mob
function [] = MN()
E = 10000;                                      % Young's modulus
nuy = 0.3;                                      % Poisson's ratio
lambda_c = 100;                                  % Dilatancy parameters
M = 2 * sqrt(2) * tan (45/180*pi());            % Friction parameters
Epsilon_v_cs = 0.05;                             % Critical state volume strain
% N0 = M*(1-exp(-lambda_c*(Epsilon_v_cs-0)));     % Dilatancy state variables
N0 = 0.0;

% Intial stress
Sigma = [150 0 0; 0 150 0;0 0 150];        
% 3 x 3 Stress tensor
Sigma1(1,1) = 150; Sigma1(1,2) = 0;
% Elastic matrix
D = E/(1+nuy)/(1-2*nuy)*[ 1-nuy nuy nuy 0 0 0
                        ; nuy 1-nuy nuy 0 0 0
                        ; nuy nuy 1-nuy 0 0 0
                        ; 0 0 0 (1-2*nuy)/2 0 0
                        ; 0 0 0 0 (1-2*nuy)/2 0
                        ; 0 0 0 0 0 (1-2*nuy)/2];
% Strain controled
Epsilon = zeros(6,1);
deps = [0.0001; -0.00005; -0.00005; 0; 0; 0];
N=N0;
Nk(1,1) = N0;

for k = 1:1000
    k
% Plastic multiplier 
lambda = 0;    

% Stress invariants
I1 = trace(Sigma);
I2 = 0.5*(trace(Sigma)^2 - trace(Sigma^2));
I3 = det (Sigma);
X = sqrt (I1*I2/I3 - 9);

if X == 0 % No shearing
    % Elastic stress update
    sig0 =  [Sigma(1,1); Sigma(2,2); Sigma(3,3); Sigma(1,2); Sigma(1,3); Sigma(2,3)];
    % Stress vector
    sig = sig0 + D * deps;   
    % Stress tensor
    Sigma = [sig(1) sig(4) sig(5); sig(4) sig(2) sig(6);sig(5) sig(6) sig(3)];     
    
    continue 
end    % Skip the isotropic pressure

% Elastic trial
sig0 =  [Sigma(1,1); Sigma(2,2); Sigma(3,3); Sigma(1,2); Sigma(1,3); Sigma(2,3)];
sig = sig0 + D * deps;        % Stress vector
y = X - (M+N);                % Yield function

if y > 0
        i = 0;
        while abs(y) > 0.001
            i = i + 1
            
            % Vector B
            Sigma_Inverse = inv(Sigma);
            Sigma_Inverse_Dev = Sigma_Inverse - 1/3*eye(3,3)*trace(Sigma_Inverse);
            Sigma_Inverse_Dev_Vector = [Sigma_Inverse_Dev(1,1); Sigma_Inverse_Dev(2,2); Sigma_Inverse_Dev(3,3); Sigma_Inverse_Dev(1,2); Sigma_Inverse_Dev(1,3); Sigma_Inverse_Dev(2,3)];
            N_vector = [N/3;N/3;N/3;0;0;0];
            B = N_vector + I1 / X * Sigma_Inverse_Dev_Vector;
            
            % A
            A = lambda_c * (M - N) * N;
            
            % Derivative
            dYdI1       = I2/(I3*2*X);
            dYdI2       = I1/(I3*2*X);
            dYdI3       = -I1*I2/(I3^2*2*X);
            dI1dSigma   = eye(3,3);
            dI2dSigma   = [ Sigma(2,2)+Sigma(3,3)   -2*Sigma(1,2)           -2*Sigma(1,3)
                          ; -2*Sigma(1,2)           Sigma(1,1)+Sigma(3,3)   -2*Sigma(2,3)
                          ; -2*Sigma(1,3)           -2*Sigma(2,3)           Sigma(2,2)+Sigma(1,1)];
            dI3dSigma   = [ Sigma(2,2)*Sigma(3,3)-Sigma(2,3)^2              2*(Sigma(2,3)*Sigma(1,3)-Sigma(3,3)*Sigma(1,2)) 2*(Sigma(2,3)*Sigma(1,2)-Sigma(2,2)*Sigma(1,3))
                          ; 2*(Sigma(2,3)*Sigma(1,3)-Sigma(3,3)*Sigma(1,2)) Sigma(1,1)*Sigma(3,3)-Sigma(1,3)^2              2*(Sigma(1,3)*Sigma(1,2)-Sigma(1,1)*Sigma(2,3))
                          ; 2*(Sigma(2,3)*Sigma(1,2)-Sigma(2,2)*Sigma(1,3)) 2*(Sigma(1,3)*Sigma(1,2)-Sigma(1,1)*Sigma(2,3)) Sigma(2,2)*Sigma(1,1)-Sigma(1,3)^2];
            
            % Tensor 3x3
            dYdSigma = dYdI1 * dI1dSigma + dYdI2 * dI2dSigma + dYdI3 * dI3dSigma;
            
            % Vector
            dYdSigma = [dYdSigma(1,1) dYdSigma(2,2) dYdSigma(3,3) dYdSigma(1,2) dYdSigma(1,3) dYdSigma(2,3)];
            
%             dlambda = y / (A - dYdSigma * D * B);
            dlambda = y / (-dYdSigma * D * B);
            
            lambda = lambda + dlambda;
            
            % Flow rule           
            depsp = - lambda *B;
            
            % Update stress
            depse = deps - depsp;
            sig = sig0 + D * depse;
            
            % Update state variables
%             Epsilon_v_plastic = depsp(1) + depsp(2) + depsp(3);
%             N =  N0 + lambda_c*(M+N)*Epsilon_v_plastic;
            
            % Liquefaction cut-off on p
            p = sig(1)+sig(2)+sig(3);
            if p < 0 
                sig = zeros(6,1);
            end
            
            % Recompute yield function
            % Stress invariants
            Sigma = [sig(1) sig(4) sig(5); sig(4) sig(2) sig(6);sig(5) sig(6) sig(3)];     
            I1 = trace(Sigma);
            I2 = 0.5*(trace(Sigma)^2 - trace(Sigma^2));
            I3 = det (Sigma);
            X = sqrt (I1*I2/I3 - 9);
            
            y = X - (M+N);
            
        end
        
end
    sig0 = sig;
    N0 = N;
    Sigma = [sig(1) sig(4) sig(5); sig(4) sig(2) sig(6);sig(5) sig(6) sig(3)];     
    Epsilon(:,k) = Epsilon(:,k-1)+deps;
    % Mean stress
    Sigma1(1,k) = (sig(1)+sig(2)+sig(3))/3;   
    % Deviatoric stress
    Sigma1(2,k) = sqrt(((sig(2)-sig(1))^2 + (sig(3)-sig(1))^2 + (sig(2)-sig(3))^2) / 2 + 3*(sig(4)^2+sig(5)^2+sig(6)^2));
    % Volume strain 
    Epsilon1(1,k) = (Epsilon(1,k)+Epsilon(2,k)+Epsilon(3,k))/3; 
    % Deviatoric stress
    Epsilon1(2,k) = 1/3*sqrt(((Epsilon(2,k)-Epsilon(1,k))^2 + (Epsilon(3,k)-Epsilon(1,k))^2 + (Epsilon(2,k)-Epsilon(3,k))^2) * 2 + 3*(Epsilon(4,k)^2+Epsilon(5,k)^2+Epsilon(6,k)^2));
    Nk(1,k) = N;
    yk(1,k) = y;
end

figure (1)    
subplot(2,2,1); hold on;
plot(Sigma1(1,:),Sigma1(2,:)); axis equal;
xlabel('$p$','interpreter','latex')
ylabel('$q$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,2); hold on;
plot(Epsilon1(2,:),Sigma1(2,:))
xlabel('${\varepsilon_q}$','interpreter','latex')
ylabel('$q$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,3); hold on; axis ij;
plot(Sigma1(1,:),Epsilon1(1,:))
xlabel('$p$','interpreter','latex')
ylabel('${\varepsilon_p}$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
subplot(2,2,4); hold on; axis ij
plot(Epsilon1(2,:),Epsilon1(1,:))
xlabel('${\varepsilon_q}$','interpreter','latex')
ylabel('${\varepsilon_p}$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

% figure (2)
% plot(Epsilon1(2,:),yk(1,:))
% xlabel('${\varepsilon_q}$','interpreter','latex')
% ylabel('$y$','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')

% figure (3)
% plot(Epsilon1(2,:),Nk(1,:))
% xlabel('${\varepsilon_q}$','interpreter','latex')
% ylabel('$N$','interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
