clear all;
clc;
close all;
%System Matrices
 
A = [0,1;-0.89,1.8];
 
B = [0;1];
 
H = [1,0];
 
Q = 0.1*eye(2);
 
R = 0.1;
 
G = eye(2);
 
wk = sqrt(0.1)*randn(2,100);            %Process noise
 
vk = sqrt(0.1)*randn(100,1);            %Measurement Noise

P = 37* eye(2);                       % Initalizing Error covariance

I = eye(2);

k = 100;

x0 = [0 0];                           %Initial value of x
 
uk = ones(k,1);                         %Unit Step
 
x = zeros(100,2);                       %Initalizing states
 
xhat = zeros(2,100);

xhatn = zeros(2,100);

xhat(1,:) = 0;

x(1,:) = x0';

%State Simulation
for k = 1:100
    
    x(k+1,:) = (A*x(k,:)'+ B*uk(k,1)'+ G*wk(:,k))';
    zk(k,:)= (H*x(k,:)'+vk(k,:)')';
    
end


% Steady State DT Kalman Filter


%Riccati Equation
for k1 = 1:99
    
     P = A*(P-P*H'*inv(H*P*H'+R)*H*P)*A'+G*Q*G';
     
end
Ks1 = P*H'*inv(H*P*H'+R)   %steady state kalman gain from Riccati equation


%DLQE for the ARE
[M,P2,Z,E] = dlqe(A,G,H,Q,R); 
 
    Ks2 = P2*H'*inv(H*P2*H'+R)  %Steady state kalman gain 
    
    

for k2=1:99

    xhatm(:,k2+1) =  A*(I-Ks2*H)*x(k2,:)' + B*uk(k2,:)' + A*Ks2*zk(k2,:);
end


%optimal DT Kalman Filter

P = 37* eye(2);                       % Initalizing Error covariance
for k = 1:99
    
    Pm= A*P*A'+ G*Q*G';
    
    xhatn(:,k+1) =  (A*xhat(:,k) + B*uk(k,:)');
    
    Ko = Pm*H'*inv(H*Pm*H'+R);
    
    P = (I-Ko*H)*Pm;
            
    xhat(:,k+1) = xhatn(:,k+1)+Ko*(zk(k+1,:)-H*xhatn(:,k+1));
    
end
 
%Comparision of first state 
figure(1)
O = plot(1:100,x(1:100,1),'-r',1:100,xhat(1,:)','-b',1:100,xhatm(1,:)','-g');
title('Kalman Filters comparision for state 1');
set(O(1), 'LineWidth', 1);
set(O(2), 'LineWidth', 1.7);
set(O(3), 'Linewidth', 1.3);
legend('x(1)','x(1) Optimal','x(2) Steady-State')
 
%Comparision of second state
figure(2)
U = plot(1:100,x(1:100,2),'-r',1:100,xhat(2,:)','-b',1:100,xhatm(2,:)','-g');
title('Kalman Filters comparision for state 2');
set(U(1), 'LineWidth', 1);
set(U(2), 'LineWidth', 1.7);
set(U(3), 'LineWidth', 1.3);
legend('x(2)','x(2) Optimal','x(2) Steady-State')
