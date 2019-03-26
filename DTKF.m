clear all
 clc
 close all
%System Matrices
 
A = [0,1;-0.89,1.8];
 
B = [0;1];
 
H = [1,0];
 
Q = 0.1*eye(2);
 
R = 0.1;
 
G = eye(2);
 
wk = sqrt(0.1)*randn(2,100);            %Process noise
 
vk = sqrt(0.1)*randn(100,1);            %Measurement Noise
 
P = 35*eye(2);                        %Initial Error Covariance P
 
I = eye(2);
 
k = 100;
 
x0 = [0 0];                           %Initial value of x
 
uk = ones(k,1);                         %Unit Step
 
x = zeros(100,2);                       %Initalizing states
 
xhat(2,1) = 0;                          %Initial x estimate values
 
x(1,:) = x0';
 
% State Calculations
for k = 1:100
    
    x(k+1,:) = (A*x(k,:)'+ B*uk(k,:)'+G*wk(:,k))';   
    
    zk(k,:)= (H*x(k,:)'+vk(k,:))';
 
end
 
% Optimal Time Varying DT Kalman Filter
 
for k = 1:99
    
    Pm= A*P*A'+ G*Q*G';
    
    xhatn(:,k+1) =  (A*xhat(:,k) + B*uk(k,:)');
 
    K = Pm*H'*inv(H*Pm*H'+R);
    
    P = (I-K*H)*Pm;
            
    xhat(:,k+1) = xhatn(:,k+1)+K*(zk(k+1,:)-H*xhatn(:,k+1));
    
end
 
%Plot for first state 
figure(1)
O= plot(1:100,x(1:100,1),'-r',1:100,xhat(1,:)','-b');
title('Optimal Time Varying DT Kalman Filter for state 1');
set(O(1), 'LineWidth', 1);
set(O(2), 'LineWidth', 1.7);
legend('x(1)','x(1)Hat(Estimate)')
 
%Plot for second state
figure(2)
U = plot(1:100,x(1:100,2),'-r',1:100,xhat(2,:)','-b');
title('Optimal Time Varying DT Kalman Filter for state 2');
set(U(1), 'LineWidth', 1);
set(U(2), 'LineWidth', 1.7);
legend('x(2)','x(2)Hat(Estimate)')

