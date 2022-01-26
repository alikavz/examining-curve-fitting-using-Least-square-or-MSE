
%%%%%%%%                   Date: Fall 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                      Curve fitting
%%%%%%%%              **** Validation ****
%%%%%%%%          Method 2: Shrinkage Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
N=20;                   %  fixed data=====>>>>Test :  fixed data method 1 (validation) =  fixed data method 2 (validation)
 %     x=[0----1];      New Observation with non-Uniform Distribution
  x =[0    0.0526    0.1053    0.1579    0.2105    0.2632    0.3158    0.3684    0.4211 0.4737    0.5263    0.5789    0.6316    0.6842    0.7368    0.7895    0.8421    0.8947 0.9474    1.0000];
 x0=x;
% Add Noise to Measurments
  Noise=[    -0.0272    0.1098   -0.0278    0.0702   -0.2052   -0.0354   -0.0824   -0.1577    0.0508 0.0282    0.0033   -0.1334    0.1127    0.0350   -0.0299    0.0023   -0.0262   -0.1750 -0.0286   -0.0831]; 
  TT=x.*x.*sin(10*x)+Noise;% Targets with Noise
  T=TT;
% Illustration
plot(x,TT,'O','LineWidth',1.5);title('Observation');grid on;hold on;% Input to Model
pause
close all
% Test  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%% Polynominal order no.: M=0
% %%%%%% Test
j=0
for landa=0:0.25:0.75;
    j=j+1;
    for i=1:N;
        a0(:,i)=[1];
        A0(i,:)=a0(:,i)';
    end
    Wl0=inv(A0'*A0+landa*eye(1))*A0'*T';
    WL0(j,:)=Wl0;
    for i=1:N
        y0(i)=a0(:,i)'*Wl0;
    end
    Y0(j,:)=y0;
    E0(j,:)=(y0-T);
    E0(j,:)=0.5*E0(j,:)*E0(j,:)'+(landa/2)*Wl0'*Wl0;
    ErmsL0(j)=sqrt(2*E0(j)/N);
    %%%%%% Test
    for i=1:length(x0);
        b0(:,i)=[1]';
        B(i)=b0(:,i)'*Wl0;
        Eb0(i,1)=B(i)-TT(i);
        Eb0(j)=Eb0(i);
    end
end
subplot(2,2,1);plot(x,Y0(1,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=0')
text(.4,6,' Landa= 0','FontSize',16)
subplot(2,2,2);plot(x,Y0(2,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=0')
text(.4,6,' Landa= 0.25', 'FontSize',16)
subplot(2,2,3);plot(x,Y0(3,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=0')
text(.4,6,' Landa= 0.5','FontSize',16)
subplot(2,2,4);plot(x,Y0(4,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=0')
text(.4,6,' Landa= 0.75','FontSize',16)
pause
close all
%%%% Polynominal order no.: M=3
%%%%%% Test
j=0
for landa=0:0.25:0.75; % landa just 4 numbers
    j=j+1;
    for i=1:N;
        a3(:,i)=[1 x(i) x(i)^2 x(i)^3];
        A3(i,:)=a3(:,i)';
    end
    Wl3=inv(A3'*A3+landa*eye(4))*A3'*T';
    WL3(j,:)=Wl3;
    for i=1:N
        y3(i)=a3(:,i)'*Wl3;
    end
    Y3(j,:)=y3;
    E3(j,:)=(Y3(j,:)-T);
    E3(j,:)=0.5*E3(j,:)*E3(j,:)'+(landa/2)*Wl3'*Wl3;
    ErmsL3(j)=sqrt(2* E3(j)/N);
    %%%%%% Test
    for i=1:length(x0);
        b3(:,i)=[1 x0(i) x0(i)^2 x0(i)^3]';
        B(i)=b3(:,i)'*Wl3;
        Eb3(i,1)=B(i)-TT(i);
        Eb(j)=Eb3(i);
    end
  end
subplot(2,2,1);plot(x,Y3(1,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=3')
text(.4,6,' Landa= 0','FontSize',16)
subplot(2,2,2);plot(x,Y3(2,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=3')
text(.4,6,' Landa= 0.25', 'FontSize',16)
subplot(2,2,3);plot(x,Y3(3,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=3')
text(.4,6,' Landa= 0.5','FontSize',16)
subplot(2,2,4);plot(x,Y3(4,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=3')
text(.4,6,' Landa= 0.75','FontSize',16)
pause
close all
%%%% Polynominal order no.: M=5
%%%%%% Test
j=0
for landa=0:0.25:0.75;
    j=j+1;
    for i=1:N;
        a5(:,i)=[1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5];
        A5(i,:)=a5(:,i)';
    end
    Wl5=inv(A5'*A5+landa*eye(6))*A5'*T';
    WL5(j,:)=Wl5;
    for i=1:N
        y5(i)=a5(:,i)'*Wl5;
    end
    Y5(j,:)=y5;
    E5(j,:)=(Y5(j,:)-T);
    E5(j,:)=0.5*E5(j,:)*E5(j,:)'+(landa/2)*Wl5'*Wl5;
    ErmsL5(j)=sqrt(2*E5(j)/N);
    %%%%%% Test
    for i=1:length(x0);
        b5(:,i)=[1 x0(i) x0(i)^2 x0(i)^3 x0(i)^4 x0(i)^5]';
        B(i)=b5(:,i)'*Wl5;
        Eb5(i,1)=B(i)-TT(i);
        Eb5(j)=Eb5(i);
    end
   end
subplot(2,2,1);plot(x,Y5(1,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=5')
text(.4,6,' Landa=0','FontSize',16)
subplot(2,2,2);plot(x,Y5(2,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=5')
text(.4,6,' Landa=0.25','FontSize',16)
subplot(2,2,3);plot(x,Y5(3,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=5')
text(.4,6,' Landa=0.5','FontSize',16)
subplot(2,2,4);plot(x,Y5(4,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=5')
text(.4,6,' Landa=0.75','FontSize',16)
pause
close all
%%%% Polynominal order no.: M=9
%%%%%% Test
j=0
for landa=0:0.25:0.75;
    j=j+1;
    for i=1:N;
        a9(:,i)=[1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5  x(i)^6 x(i)^7  x(i)^8 x(i)^9];
        A9(i,:)=a9(:,i)';
    end
    Wl9=inv(A9'*A9+landa*eye(10))*A9'*T';
    WL9(j,:)=Wl9;
    for i=1:N
        y9(i)=a9(:,i)'*Wl9;
    end
    Y9(j,:)=y9;
    E9(j,:)=(Y9(j,:)-T);
    E9(j,:)=0.5*E9(j,:)*E9(j,:)'+(landa/2)*Wl9'*Wl9;
    ErmsL9(j)=sqrt(2*E9(j)/N);
    %%%%%% Test
    for i=1:length(x0);
        b9(:,i)=[1 x0(i) x0(i)^2 x0(i)^3 x0(i)^4 x0(i)^5 x0(i)^6 x0(i)^7 x0(i)^8 x0(i)^9]';
        B(i)=b9(:,i)'*Wl9;
        Eb9(i,1)=B(i)-TT(i);
        Eb9(j)=Eb9(i);
    end
    end
subplot(2,2,1);plot(x,Y9(1,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=9')
text(.4,6,' Landa=0','FontSize',16)
subplot(2,2,2);plot(x,Y9(2,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=9')
text(.4,6,' Landa=0.25','FontSize',16)
subplot(2,2,3);plot(x,Y9(3,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=9')
text(.4,6,' Landa=0.5','FontSize',16)
subplot(2,2,4);plot(x,Y9(4,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=9')
text(.4,6,' Landa=0.75','FontSize',16)
pause
close all
%%%% Polynominal order no.: M=11
%%%%%% Test
% j=0;
% for landa=0:0.25:0.75;
%     j=j+1;
%     for i=1:N;
%         a11(:,i)=[1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5  x(i)^6 x(i)^7  x(i)^8 x(i)^9 x(i)^10 x(i)^11];
%         A11(i,:)=a11(:,i)';
%     end
%     Wl11=inv(A11'*A11+landa*eye(12))*A11'*T';
%     WL11(j,:)=Wl11;
%     for i=1:N
%         y11(i)=a11(:,i)'*Wl11;
%     end
%     Y11(j,:)=y11;
%     E11(j,:)=(Y11(j,:)-T);
%     E11(j,:)=0.5*E11(j,:)*E11(j,:)'+(landa/2)*Wl11'*Wl11;
%     ErmsL11(j)=sqrt(2*E11(j)/N);
%     %%%%%% Test
%     for i=1:length(x0);
%         b11(:,i)=[1 x0(i) x0(i)^2 x0(i)^3 x0(i)^4 x0(i)^5 x0(i)^6 x0(i)^7 x0(i)^8 x0(i)^9 x0(i)^10 x0(i)^11]';
%         B(i)=b11(:,i)'*Wl11;
%         Eb11(i,1)=B(i)-TT(i);
%         Eb11(j)=Eb11(i);
%     end
%     end
% subplot(2,2,1);plot(x,Y11(1,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=11')
% text(.4,6,' Landa=0','FontSize',16)
% subplot(2,2,2);plot(x,Y11(2,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=11')
% text(.4,6,' Landa=0.25','FontSize',16)
% subplot(2,2,3);plot(x,Y11(3,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=11')
% text(.4,6,' Landa=0.5','FontSize',16)
% subplot(2,2,4);plot(x,Y11(4,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=11')
% text(.4,6,' Landa=0.75','FontSize',16)
% pause
% close all
% %%%% Polynominal order no.: M=13
% %%%%%% Test
% j=0;
% for landa=0:0.25:0.75;
%     j=j+1;
%     for i=1:N;
%         a13(:,i)=[1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5  x(i)^6 x(i)^7  x(i)^8 x(i)^9 x(i)^10 x(i)^11 x(i)^12 x(i)^13];
%         A13(i,:)=a13(:,i)';
%     end
%     Wl13=inv(A13'*A13+landa*eye(14))*A13'*T';
%     WL13(j,:)=Wl13;
%     for i=1:N
%         y13(i)=a13(:,i)'*Wl13;
%     end
%     Y13(j,:)=y13;
%     E13(j,:)=(Y13(j,:)-T);
%     E13(j,:)=0.5*E13(j,:)*E13(j,:)'+(landa/2)*Wl13'*Wl13;
%     ErmsL13(j)=sqrt(2*E13(j)/N);
%     %%%%%% Test
%     for i=1:length(x0);
%         b13(:,i)=[1 x0(i) x0(i)^2 x0(i)^3 x0(i)^4 x0(i)^5 x0(i)^6 x0(i)^7 x0(i)^8 x0(i)^9 x0(i)^10 x0(i)^11 x0(i)^12 x0(i)^13]';
%         B(i)=b13(:,i)'*Wl13;
%         Eb13(i,1)=B(i)-TT(i);
%         Eb13(j)= Eb13(i);
%     end
%   end
% subplot(2,2,1);plot(x,Y13(1,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=13')
% text(.4,6,' Landa=0','FontSize',16)
% subplot(2,2,2);plot(x,Y13(2,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=13')
% text(.4,6,' Landa=0.25','FontSize',16)
% subplot(2,2,3);plot(x,Y13(3,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=13')
% text(.4,6,' Landa=0.5','FontSize',16)
% subplot(2,2,4);plot(x,Y13(4,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=13')
% text(.4,6,' Landa=0.75','FontSize',16)
% pause
% close all
% %%%% Polynominal order no.: M=15
% %%%%%% Test
% j=0;
% for landa=0:0.25:0.75;
%     j=j+1;
%     for i=1:N;
%         a15(:,i)=[1 x(i) x(i)^2 x(i)^3 x(i)^4 x(i)^5  x(i)^6 x(i)^7  x(i)^8 x(i)^9 x(i)^10 x(i)^11 x(i)^12 x(i)^13 x(i)^14 x(i)^15];
%         A15(i,:)=a15(:,i)';
%     end
%     Wl15=inv(A15'*A15+landa*eye(16))*A15'*T';
%     WL15(j,:)=Wl15;
%     for i=1:N
%         y15(i)=a15(:,i)'*Wl15;
%     end
%     Y15(j,:)=y15;
%     E15(j,:)=(Y15(j,:)-T);
%     E15(j,:)=0.5*E15(j,:)*E15(j,:)'+(landa/2)*Wl15'*Wl15;
%     ErmsL15(j)=sqrt(2*E15(j)/N);
%     %%%%%% Test
%     for i=1:length(x0);
%         b15(:,i)=[1 x0(i) x0(i)^2 x0(i)^3 x0(i)^4 x0(i)^5 x0(i)^6 x0(i)^7 x0(i)^8 x0(i)^9 x0(i)^10 x0(i)^11 x0(i)^12 x0(i)^13  x0(i)^14 x0(i)^15]';
%         B(i)=b15(:,i)'*Wl15;
%         Eb15(i,1)=B(i)-TT(i);
%      Eb15(j)= Eb15(i);
%     end
%    end
% subplot(2,2,1);plot(x,Y15(1,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=15')
% text(.4,6,' Landa=0','FontSize',16)
% subplot(2,2,2);plot(x,Y15(2,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=15')
% text(.4,6,' Landa=0.25','FontSize',16)
% subplot(2,2,3);plot(x,Y15(3,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=15')
% text(.4,6,' Landa=0.5','FontSize',16)
% subplot(2,2,4);plot(x,Y15(4,:),'LineWidth',1.5);hold on; plot(x,T,'or','LineWidth',1.5);grid on;title('Curve Fitting Method 2: Shrinkage, Polynominal Order no.: M=15')
% text(.4,6,' Landa=0.75','FontSize',16)
% 
