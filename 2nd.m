
%%%%%%%%                   Date: Fall 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                      Curve fitting
%%%%%%%%              **** Validation ****
%%%%%%%%               Method 1: Least square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
N=20;                   %   New Number of Observation
 %     x=[0----1];      New Observation with non-Uniform Distribution
    x =[0    0.0526    0.1053    0.1579    0.2105    0.2632    0.3158    0.3684    0.4211 0.4737    0.5263    0.5789    0.6316    0.6842    0.7368    0.7895    0.8421    0.8947 0.9474    1.0000];
% x=[0:1/(N-1):1]; % Observation with Uniform Distribution

% Noise=0.1*randn(1,N);
% Add Noise to Measurments
  Noise=[    -0.0272    0.1098   -0.0278    0.0702   -0.2052   -0.0354   -0.0824   -0.1577    0.0508 0.0282    0.0033   -0.1334    0.1127    0.0350   -0.0299    0.0023   -0.0262   -0.1750 -0.0286   -0.0831];
T=x.*x.*sin(10*x)+Noise ;% Targets with Noise
T=T+0.7*Noise;
% Illustration
plot(x,T,'O','LineWidth',1.5);title('Observation');grid on;hold on;% Input to Model
pause
close all
% %%%%% Polynominal order no.: M=3
% %%%%%% TEST
for i=1:N
    a3(:,i)=[1 x(i) x(i)^2 x(i)^3];
    A3(i,:)=a3(:,i)';
end
W3=inv(A3'*A3)*A3'*T';
for i=1:N
    y3(i)=a3(:,i)'*W3;
end
E3=(y3-T);
E3=0.5*E3*E3';
Erms3=sqrt(2*E3/N)
plot(x,y3,'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);title('Curve Fitting: Validation Method 1: Least Square, Polynominal Order no.: M=3');legend('M=3','NewTargets','Location','northwest');grid on
pause
close all
% %%%% Polynominal order no.: M=5
% %%%%%% TEST 
for i=1:N
    a5(:,i)=[1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5];
    A5(i,:)=a5(:,i)';
end
W5=inv(A5'*A5)*A5'*T';
for i=1:N
    y5(i)=a5(:,i)'*W5;
end
E5=(y5-T);
E5=0.5*E5*E5';
Erms5=sqrt(2*E5/N)
plot(x,y5,'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);title('Curve Fitting: Validation Method 1: Least Square, Polynominal Order no.: M=5');legend('M=5','New Targets','Location','northwest');grid on
pause
close all
% %%%% Polynominal order no.: M=9
% %%%%% TEST 
for i=1:N
    a9(:,i)=[1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5  x(i)^6 x(i)^7  x(i)^8 x(i)^9];
    A9(i,:)=a9(:,i)';
end
W9=inv(A9'*A9)*A9'*T';
for i=1:N
    y9(i)=a9(:,i)'*W9;
end
E9=(y9-T);
E9=0.5*E9*E9';
Erms9=sqrt(2*E9/N)
plot(x,y9,'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);title('Curve Fitting: Validation Method 1: Least Square, Polynominal Order no.: M=9');legend('M=9','New Targets','Location','northwest');grid on
pause
close all
%%%% Polynominal order no.: M=11
%%%% TEST 
% for i=1:N
%     a11(:,i)=[1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5  x(i)^6 x(i)^7  x(i)^8 x(i)^9 x(i)^10 x(i)^11];
%     A11(i,:)=a11(:,i)';
% end
% W11=inv(A11'*A11)*A11'*T';
% for i=1:N
%     y11(i)=a11(:,i)'*W11;
% end
% E11=(y11-T);
% E11=0.5*E11*E11';
% Erms11=sqrt(2*E11/N)
%  plot(x,y9,'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);title('Curve Fitting: Validiation Method 1: Least Square, Polynominal Order no.: M=11');legend('M=11','New Targets','Location','northwest');grid on
