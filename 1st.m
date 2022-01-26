                                    

%%%%%%%%                   Date: Fall 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                      Curve fitting
%%%%%%%%              Method 1: Least square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
N=20;% Number of Observation
M=9;% Degree of Polynomials M= 0,3,5,9
x0=[0:0.001:1];
x=[0:1/(N-1):1]; % Observation with Uniform Distribution
% Add Noise to Measurments
Noise=0.1*randn(1,N);
T=x.*x.*sin(10*x)+Noise ;% Targets with Noise

%%%%%%%%%%%%%%%%%%%%%     fixed data=====>>>>Training :  fixed data method 1 =  fixed data method 2
  Noise = [     -0.1089    0.0033    0.0553    0.1101    0.1544    0.0086   -0.1492   -0.0742   -0.1062 0.2350   -0.0616    0.0748   -0.0192    0.0889   -0.0765   -0.1402   -0.1422    0.0488 -0.0177   -0.0196];
  T = [    -0.0762    0.0037    0.0483    0.1020    0.1462    0.0398   -0.1060   -0.1221   -0.2297 -0.0598   -0.2791   -0.1065   -0.0005    0.3104    0.4266    0.5246    0.4986    0.4020   -0.0563   -0.5577];
T=T+0.7*Noise;
%%%%%%%%%%%%%%%%%%%%%
% save Noise;
% load Noise;

TT=x0.*x0.*sin(10*x0);% Targets without Noise
% Illustration
plot(x,T,'O','LineWidth',1.5);title('Observation');grid on;hold on;% Input to Model
pause
plot(x0,x0.*x0.*sin(10*x0),'g','LineWidth',1.5);grid on;hold on;% Underlying Function
pause
close all
% Training Phase: w Computation By minimizing MSE
%%%%% Polynominal order no.: M=0
%%%%%% Training 
E0=0;M0=0;
for i=1:N
    a0(:,i)=[1];
    A0(i,:)=a0(:,i)';
end
W0=inv(A0'*A0)*A0'*T';
for i=1:N
    y0(i)=a0(:,i)'*W0;
end
E0=(y0-T);
E0=0.5*E0*E0';
Erms0=sqrt(2*E0/N)
%%%%%% Test 
for i=1:length(x0)
     b0(:,i)=[1]';
     yyb(i)=b0(:,i)'*W0;
     Eb0(i,1)=yyb(i)-TT(i);
end
Eb0=0.5*Eb0'*Eb0;
Ermsb0=sqrt(2*Eb0/(length(x0)))
plot(x,y0,'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);title('Curve Fitting Method 1: Least Square, Polynominal Order no.: M=0');legend('M=0','Targets','Location','northwest');grid on
pause
close all
%%%%% Polynominal order no.: M=3
%%%%%% Training 
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
%%%%%% Test 
for i=1:length(x0)
     b3(:,i)=[1 x0(i) x0(i)^2 x0(i)^3]';
     yyb(i)=b3(:,i)'*W3;
     Eb3(i,1)=yyb(i)-TT(i);
end
Eb3=0.5*Eb3'*Eb3;
Ermsb3=sqrt(2*Eb3/(length(x0)))
plot(x,y3,'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);title('Curve Fitting Method 1: Least Square, Polynominal Order no.: M=3');legend('M=3','Targets','Location','northwest');grid on
pause
close all
%%%% Polynominal order no.: M=5
%%%%%% Training 
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
%%%%%% Test 
for i=1:length(x0)
     b5(:,i)=[1 x0(i) x0(i)^2 x0(i)^3 x0(i)^4 x0(i)^5]';
     yyb(i)=b5(:,i)'*W5;
     Eb5(i,1)=yyb(i)-TT(i);
end
Eb5=0.5*Eb5'*Eb5;
Ermsb5=sqrt(2*Eb5/(length(x0)))
plot(x,y5,'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);title('Curve Fitting Method 1: Least Square, Polynominal Order no.: M=5');legend('M=5','Targets','Location','northwest');grid on
pause
close all
%%%% Polynominal order no.: M=9
%%%%% Training 
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
%%%%%% Test 
for i=1:length(x0)
    b9(:,i)=[1 x0(i) x0(i)^2 x0(i)^3 x0(i)^4 x0(i)^5 x0(i)^6 x0(i)^7 x0(i)^8 x0(i)^9]';
     yyb(i)=b9(:,i)'*W9;
     Eb9(i,1)=yyb(i)-TT(i);
end
Eb9=0.5*Eb9'*Eb9;
Ermsb9=sqrt(2*Eb9/(length(x0)))
plot(x,y9,'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);title('Curve Fitting Method 1: Least Square, Polynominal Order no.: M=9');legend('M=9','Targets','Location','northwest');grid on
pause
close all

%%%%%%%% Training and Test Error Curve Fitting Method 1: Least Square
M=[0 3 5 9 ];
Erms=[Erms0 Erms3 Erms5 Erms9 ];
Ermsb=[Ermsb0 Ermsb3 Ermsb5 Ermsb9 ];
plot(M,Erms,'-rs',M,Ermsb,'-gs','LineWidth',2);legend('Training','Test','Location','northwest');grid on;hold on;
xlabel('Polynominal order no.: M');ylabel('Erms');title('Training and Test Error');
W0'
W3'
W5'
W9'



