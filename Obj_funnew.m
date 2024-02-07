function y = Obj_funnew(z)   
% function y = Obj_funnew(z)   
syms x x1  y1  y2  
% alpha = 0.4;  beta = 4;  LC=20; LP=15; Cins=3; CPM=20; CCM=40;  a1=0.5; b1=1;a2=2; b2=5;  q1=20; q2=30; Cunit=3; rou=2.5;Ch=0.05;
% h=3; S=100;TE=60;
alpha = 0.4;  beta = 4;  LC=20; LP=15; Cins=0.2; CPM=1; CCM=4;  a1=0.5; b1=1;a2=2; b2=5;  q1=20; q2=30; Cunit=0.3; rou=0.5;Ch=0.005;TE=40;
% h=3; S=100;
h=z(1); S=z(2);
u1(y1)=unifpdf(y1, a1, b1);
u2(y2)=unifpdf(y2, a2, b2);
ETPM=(a1+b1)/2;
ETCM=(a2+b2)/2;
% PA=int((int((gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),x1,LP-x,LC-x)),x,0,LP);
% PB=int((int((gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),x1,LC-x,inf)),x,0,LP);

PN0 = @(x) integral(@(x1) (gampdf(x, alpha*max(1,(floor((TE)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S/q1-(floor((TE)/h))*h)), beta)), 0, LC-x);
x_range = 0:0.01:LP;  % 定义 x 的范围和步长
result3 = zeros(size(x_range));
for i = 1:length(x_range)
      result3(i) = PN0(x_range(i));
end
 % 数值积分结果
PN0 = trapz(x_range, result3);
% PN0
CIS0= floor((TE)/h)*Cins;
CPS0=(q1*(TE-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
CHS0=(S*(TE)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
EC=PN0*(CIS0+CPS0+CHS0);
if TE<h
    n=0;
else
    T1=0; C1I=0; C1M=0;C1H=0;C1S=0;
    for j =2:max(2,floor((TE)/h))
      PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
      PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
      % 外部积分
      x_range = 0:0.01:LP;  % 定义 x 的范围和步长
      result1 = zeros(size(x_range));
      for i = 1:length(x_range)
        result1(i) = PA(x_range(i));
      end
      result2 = zeros(size(x_range));
      for i = 1:length(x_range)
        result2(i) = PB(x_range(i));
      end
      % 数值积分结果
      PA = trapz(x_range, result1);
      PB = trapz(x_range, result2);
      T1=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T1;
      C1I=(PA*j+PB*(j+1))*Cins+C1I;
      C1M=PA*CPM+PB*CCM+C1M;
      C1H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C1H;
      C1S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C1S;
    end
    T1=double(T1);
%     T1
    C1P=0;
    for j=2:max(2,floor((TE)/h))
      PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
      PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
      % 外部积分
      x_range = 0:0.01:LP;  % 定义 x 的范围和步长
      result1 = zeros(size(x_range));
      for i = 1:length(x_range)
        result1(i) = PA(x_range(i));
      end
      result2 = zeros(size(x_range));
      for i = 1:length(x_range)
        result2(i) = PB(x_range(i));
      end
      % 数值积分结果
      PA = trapz(x_range, result1);
      PB = trapz(x_range, result2);
      C1P=(PA*((T1-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T1-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C1P;
    end
    S1=T1;
    PN1 = @(x) integral(@(x1) (gampdf(x, alpha*max(1,(floor((TE-S1)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-T1-S/q1-(floor((TE-T1)/h))*h)), beta)), 0, LC-x);
    x_range = 0:0.01:LP;  % 定义 x 的范围和步长
    result3 = zeros(size(x_range));
    for i = 1:length(x_range)
      result3(i) = PN1(x_range(i));
    end
 % 数值积分结果
    PN1 = trapz(x_range, result3);
%     PN1
    CIS1= floor((TE-S1)/h)*Cins;
    CPS1=(q1*(TE-S1-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
    CHS1=(S*(TE-S1)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
    EC=PN1*(CIS1+CPS1+CHS1+C1I+C1H+C1S+C1M+C1P)+EC;
    if isAlways(TE-S1<h) 
        n=1;
    elseif T1==0
        n=0;
    else
        T2=0;  C2I=0; C2M=0;C2H=0;C2S=0;
        for j =2:max(2,floor((TE-S1)/h))
           PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
           PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
           % 外部积分
           x_range = 0:0.01:LP;  % 定义 x 的范围和步长
           result1 = zeros(size(x_range));
           for i = 1:length(x_range)
             result1(i) = PA(x_range(i));
           end
           result2 = zeros(size(x_range));
           for i = 1:length(x_range)
             result2(i) = PB(x_range(i));
           end
           % 数值积分结果
           PA = trapz(x_range, result1);
           PB = trapz(x_range, result2);
           T2=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T2;
           C2I=(PA*j+PB*(j+1))*Cins+C2I;
           C2M=PA*CPM+PB*CCM+C2M;
           C2H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C2H;
           C2S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C2S;
        end
           T2=double(T2);
           C2P=0;
           for j=2:max(2,floor((TE-S1)/h))
                PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                % 外部积分
                x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                result1 = zeros(size(x_range));
                for i = 1:length(x_range)
                  result1(i) = PA(x_range(i));
                end
                result2 = zeros(size(x_range));
                for i = 1:length(x_range)
                  result2(i) = PB(x_range(i));
                end
                % 数值积分结果
                PA = trapz(x_range, result1);
                PB = trapz(x_range, result2);
                C2P=(PA*((T2-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T2-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C2P;
           end
           S2=T1+T2;
           PN2 = @(x) integral(@(x1) (gampdf(x, alpha*max(1,(floor((TE-S2)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S2-S/q1-(floor((TE-S2)/h))*h)), beta)), 0, LC-x);
           x_range = 0:0.01:LP;  % 定义 x 的范围和步长
           result3 = zeros(size(x_range));
           for i = 1:length(x_range)
              result3(i) = PN2(x_range(i));
           end
           PN2 = trapz(x_range, result3);
%            PN2
           CIS2= floor((TE-S2)/h)*Cins;
           CPS2=(q1*(TE-S2-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
           CHS2=(S*(TE-S2)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
           EC=EC+PN2*(CIS2+CPS2+CHS2+C2I+C2H+C2S+C2M+C2P);
           if isAlways(TE-S2<h) 
               n=2;
           elseif T2==0
               n=1;
           else
               T3=0;C3I=0; C3M=0;C3H=0;C3S=0;
               for j =2:max(2,floor((TE-S2)/h))
                  PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                  PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                  % 外部积分
                  x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                  result1 = zeros(size(x_range));
                  for i = 1:length(x_range)
                    result1(i) = PA(x_range(i));
                  end
                  result2 = zeros(size(x_range));
                  for i = 1:length(x_range)
                     result2(i) = PB(x_range(i));
                  end
                  % 数值积分结果
                  PA = trapz(x_range, result1);
                  PB = trapz(x_range, result2);
                  T3=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T3;
                  C3I=(PA*j+PB*(j+1))*Cins+C3I;
                  C3M=PA*CPM+PB*CCM+C3M;
                  C3H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C3H;
                  C3S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C3S;
               end
               T3=double(T3);
               C3P=0;
               for j=2:max(2,floor((TE-S2)/h))
                  PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                  PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                  % 外部积分
                  x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                  result1 = zeros(size(x_range));
                  for i = 1:length(x_range)
                    result1(i) = PA(x_range(i));
                  end
                  result2 = zeros(size(x_range));
                  for i = 1:length(x_range)
                     result2(i) = PB(x_range(i));
                  end
                  % 数值积分结果
                  PA = trapz(x_range, result1);
                  PB = trapz(x_range, result2);
                  C3P=(PA*((T3-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T3-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C3P;
               end
               S3=T1+T2+T3;
%                T3
%                S3
               PN3 = @(x) integral(@(x1) (gampdf(x, alpha*max(1,(floor((TE-S3)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S3-S/q1-(floor((TE-S3)/h))*h)), beta)), 0, LC-x);
               x_range = 0:0.01:LP;  % 定义 x 的范围和步长
               result3 = zeros(size(x_range));
               for i = 1:length(x_range)
                  result3(i) = PN3(x_range(i));
               end
               PN3 = trapz(x_range, result3);
%                PN3
               CIS3= floor((TE-S3)/h)*Cins;
               CPS3=(q1*(TE-S3-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
               CHS3=(S*(TE-S3)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
               EC=EC+PN3*(CIS3+CPS3+CHS3+C3I+C3H+C3S+C3M+C3P);
               if isAlways(TE-S3<h) 
                   n=3;
               elseif T3==0
                   n=2;
               else
                   T4=0;C4I=0; C4M=0;C4H=0;C4S=0;
                   for j =2:max(2,floor((TE-S3)/h))
                     PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                     PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                     % 外部积分
                     x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                     result1 = zeros(size(x_range));
                     for i = 1:length(x_range)
                        result1(i) = PA(x_range(i));
                     end
                     result2 = zeros(size(x_range));
                     for i = 1:length(x_range)
                        result2(i) = PB(x_range(i));
                     end
                     % 数值积分结果
                     PA = trapz(x_range, result1);
                     PB = trapz(x_range, result2);
                     T4=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T4;
                     C4I=(PA*j+PB*(j+1))*Cins+C4I;
                     C4M=PA*CPM+PB*CCM+C4M;
                     C4H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C4H;
                     C4S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C4S;
                   end
                   T4=double(T4);
                   C4P=0;
                   for j=2:max(2,floor((TE-S3)/h))
                        PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                        PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                         % 外部积分
                        x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                        result1 = zeros(size(x_range));
                        for i = 1:length(x_range)
                          result1(i) = PA(x_range(i));
                        end
                        result2 = zeros(size(x_range));
                        for i = 1:length(x_range)
                            result2(i) = PB(x_range(i));
                        end
                        % 数值积分结果
                        PA = trapz(x_range, result1);
                        PB = trapz(x_range, result2);
                        C4P=(PA*((T4-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T4-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C4P;
                   end
                   S4=T1+T2+T3+T4;
%                    T4
%                    S4
                   PN4 = @(x) integral(@(x1) (gampdf(x, alpha*max(1,(floor((TE-S4)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S4-S/q1-(floor((TE-S4)/h))*h)), beta)), 0, LC-x);
                   x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                   result3 = zeros(size(x_range));
                   for i = 1:length(x_range)
                      result3(i) = PN4(x_range(i));
                   end
                   PN4 = trapz(x_range, result3);
%                    PN4
                   CIS4= floor((TE-S4)/h)*Cins;
                   CPS4=(q1*(TE-S4-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
                   CHS4=(S*(TE-S4)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
                   EC=EC+PN4*(CIS4+CPS4+CHS4+C4I+C4H+C4S+C4M+C4P);
                   if isAlways(TE-S4<h)
                       n=4;
                   elseif T4==0
                       n=3;
                   else
                       T5=0;C5I=0; C5M=0;C5H=0;C5S=0;
                       for j =2:max(2,floor((TE-S4)/h))
                           PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                           PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                            % 外部积分
                           x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                           result1 = zeros(size(x_range));
                           for i = 1:length(x_range)
                              result1(i) = PA(x_range(i));
                           end
                           result2 = zeros(size(x_range));
                           for i = 1:length(x_range)
                             result2(i) = PB(x_range(i));
                           end
                           % 数值积分结果
                           PA = trapz(x_range, result1);
                           PB = trapz(x_range, result2);
                           T5=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T5;
                           C5I=(PA*j+PB*(j+1))*Cins+C5I;
                           C5M=PA*CPM+PB*CCM+C5M;
                           C5H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C5H;
                           C5S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C5S;
                       end
                       T5=double(T5);
                       C5P=0;
                       for j=2:max(2,floor((TE-S4)/h))
                           PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                           PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                            % 外部积分
                           x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                           result1 = zeros(size(x_range));
                           for i = 1:length(x_range)
                              result1(i) = PA(x_range(i));
                           end
                           result2 = zeros(size(x_range));
                           for i = 1:length(x_range)
                             result2(i) = PB(x_range(i));
                           end
                           % 数值积分结果
                           PA = trapz(x_range, result1);
                           PB = trapz(x_range, result2);
                           C5P=(PA*((T5-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T5-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C5P;
                       end
                       S5=T1+T2+T3+T4+T5;
%                        T5
%                        S5
                       PN5 = @(x) integral(@(x1) (gampdf(x, alpha*max(1,(floor((TE-S5)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S5-S/q1-(floor((TE-S5)/h))*h)), beta)), 0, LC-x);
                       x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                       result3 = zeros(size(x_range));
                       for i = 1:length(x_range)
                          result3(i) = PN5(x_range(i));
                       end
                       PN5 = trapz(x_range, result3);
%                        PN5
                       CIS5= floor((TE-S5)/h)*Cins;
                       CPS5=(q1*(TE-S5-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
                       CHS5=(S*(TE-S5)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
                       EC=EC+PN5*(CIS5+CPS5+CHS5+C5I+C5H+C5S+C5M+C5P);
                       if isAlways(TE-S5<h)
                           n=5;
                       elseif T5==0
                           n=4;
                       else
                           T6=0;C6I=0; C6M=0;C6H=0;C6S=0;
                           for j=2:max(2,floor((TE-S5)/h))
                                PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                                PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                                % 外部积分
                                x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                result1 = zeros(size(x_range));
                                for i = 1:length(x_range)
                                  result1(i) = PA(x_range(i));
                                end
                                result2 = zeros(size(x_range));
                                for i = 1:length(x_range)
                                  result2(i) = PB(x_range(i));
                                end
                                % 数值积分结果
                                PA = trapz(x_range, result1);
                                PB = trapz(x_range, result2);
                                T6=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T6;
                                C6I=(PA*j+PB*(j+1))*Cins+C6I;
                                C6M=PA*CPM+PB*CCM+C6M;
                                C6H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C6H;
                                C6S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C6S;
                           end
                           T6=double(T6);
                           C6P=0;
                           for j=2:max(2,floor((TE-S5)/h))
                                PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                                PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                                % 外部积分
                                x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                result1 = zeros(size(x_range));
                                for i = 1:length(x_range)
                                  result1(i) = PA(x_range(i));
                                end
                                result2 = zeros(size(x_range));
                                for i = 1:length(x_range)
                                  result2(i) = PB(x_range(i));
                                end
                                % 数值积分结果
                                PA = trapz(x_range, result1);
                                PB = trapz(x_range, result2);
                                C6P=(PA*((T6-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T6-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C6P;
                           end
                           S6=T1+T2+T3+T4+T5+T6;
                           PN6 = @(x) integral(@(x1) (gampdf(x, alpha*max(1,(floor((TE-S6)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S6-S/q1-(floor((TE-S6)/h))*h)), beta)), 0, LC-x);
                           x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                           result3 = zeros(size(x_range));
                           for i = 1:length(x_range)
                              result3(i) = PN6(x_range(i));
                           end
                           PN6 = trapz(x_range, result3);
                           CIS6= floor((TE-S6)/h)*Cins;
                           CPS6=(q1*(TE-S6-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
                           CHS6=(S*(TE-S6)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
                           EC=EC+PN6*(CIS6+CPS6+CHS6+C6I+C6H+C6S+C6M+C6P);
                           if isAlways(TE-S6<h)
                               n=6;
                           elseif T6==0
                               n=5;
                           else
                               T7=0;C7I=0; C7M=0;C7H=0;C7S=0;
                               for j=2:max(2,floor((TE-S6)/h))
                                   PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                                   PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                                   % 外部积分
                                   x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                   result1 = zeros(size(x_range));
                                   for i = 1:length(x_range)
                                      result1(i) = PA(x_range(i));
                                   end
                                   result2 = zeros(size(x_range));
                                   for i = 1:length(x_range)
                                     result2(i) = PB(x_range(i));
                                   end
                                   % 数值积分结果
                                   PA = trapz(x_range, result1);
                                   PB = trapz(x_range, result2);
                                   T7=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T7;
                                   C7I=(PA*j+PB*(j+1))*Cins+C7I;
                                   C7M=PA*CPM+PB*CCM+C7M;
                                   C7H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C7H;
                                   C7S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C7S;
                               end
                               T7=double(T7);
                               C7P=0;
                               for j=2:max(2,floor((TE-S6)/h))
                                   PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                                   PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                                   % 外部积分
                                   x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                   result1 = zeros(size(x_range));
                                   for i = 1:length(x_range)
                                      result1(i) = PA(x_range(i));
                                   end
                                   result2 = zeros(size(x_range));
                                   for i = 1:length(x_range)
                                     result2(i) = PB(x_range(i));
                                   end
                                   % 数值积分结果
                                   PA = trapz(x_range, result1);
                                   PB = trapz(x_range, result2);
                                   C7P=(PA*((T7-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T7-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C7P;
                               end
                               S7=T1+T2+T3+T4+T5+T6+T7;
                               n=7;
                               PN7 = @(x) integral(@(x1) (gampdf(x, alpha*(max(1,floor((TE-S7)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S7-S/q1-(floor((TE-S7)/h))*h)), beta)), 0, LC-x);
                               x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                               result3 = zeros(size(x_range));
                               for i = 1:length(x_range)
                                 result3(i) = PN7(x_range(i));
                               end
                               PN7 = trapz(x_range, result3);
                               CIS7= floor((TE-S7)/h)*Cins;
                               CPS7=(q1*(TE-S7-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
                               CHS7=(S*(TE-S7)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
                               EC=EC+PN7*(CIS7+CPS7+CHS7+C7I+C7H+C7S+C7M+C7P);
                               if isAlways(TE-S7<h)
                                 n=7;
                               elseif T7==0
                                 n=6;
                               else
                                 T8=0;C8I=0; C8M=0;C8H=0;C8S=0;
                                 for j=2:max(2,floor((TE-S7)/h))
                                   PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                                   PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                                   % 外部积分
                                   x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                   result1 = zeros(size(x_range));
                                   for i = 1:length(x_range)
                                      result1(i) = PA(x_range(i));
                                   end
                                   result2 = zeros(size(x_range));
                                   for i = 1:length(x_range)
                                     result2(i) = PB(x_range(i));
                                   end
                                   % 数值积分结果
                                   PA = trapz(x_range, result1);
                                   PB = trapz(x_range, result2);
                                   T8=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T8;
                                   C8I=(PA*j+PB*(j+1))*Cins+C8I;
                                   C8M=PA*CPM+PB*CCM+C8M;
                                   C8H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C8H;
                                   C8S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C8S;
                                 end
                                 T8=double(T8);
                                 C8P=0;
                                 for j=2:max(2,floor((TE-S7)/h))
                                   PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                                   PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                                   % 外部积分
                                   x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                   result1 = zeros(size(x_range));
                                   for i = 1:length(x_range)
                                      result1(i) = PA(x_range(i));
                                   end
                                   result2 = zeros(size(x_range));
                                   for i = 1:length(x_range)
                                     result2(i) = PB(x_range(i));
                                   end
                                   % 数值积分结果
                                   PA = trapz(x_range, result1);
                                   PB = trapz(x_range, result2);
                                   C8P=(PA*((T8-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T8-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C8P;
                                 end
                                 S8=T1+T2+T3+T4+T5+T6+T7+T8;
                                 n=8;
                                 PN8 = @(x) integral(@(x1) (gampdf(x, alpha*(max(1,floor((TE-S8)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S8-S/q1-(floor((TE-S8)/h))*h)), beta)), 0, LC-x);
                                 x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                 result3 = zeros(size(x_range));
                                 for i = 1:length(x_range)
                                   result3(i) = PN8(x_range(i));
                                 end
                                 PN8 = trapz(x_range, result3);
%                                PN8
                                 CIS8= floor((TE-S8)/h)*Cins;
                                 CPS8=(q1*(TE-S8-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
                                 CHS8=(S*(TE-S8)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
                                 EC=EC+PN8*(CIS8+CPS8+CHS8+C8I+C8H+C8S+C8M+C8P);
                                 if isAlways(TE-S8<h)
                                     n=8;
                                 elseif T8==0
                                     n=7;
                                 else
                                     T9=0;C9I=0; C9M=0;C9H=0;C9S=0;
                                     for j=2:max(2,floor((TE-S8)/h))
                                       PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                                       PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                                       % 外部积分
                                       x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                       result1 = zeros(size(x_range));
                                       for i = 1:length(x_range)
                                         result1(i) = PA(x_range(i));
                                       end
                                       result2 = zeros(size(x_range));
                                       for i = 1:length(x_range)
                                         result2(i) = PB(x_range(i));
                                       end
                                        % 数值积分结果
                                       PA = trapz(x_range, result1);
                                       PB = trapz(x_range, result2);
                                       T9=PA*(j*h+max(ETPM,S/q1))+PB*((j-1)*h+PB*h+max(ETCM,S/q1))+T9;
                                       C9I=(PA*j+PB*(j+1))*Cins+C9I;
                                       C9M=PA*CPM+PB*CCM+C9M;
                                       C9H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C9H;
                                       C9S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C9S;
                                     end
                                     T9=double(T9);
                                     C9P=0;
                                     for j=2:max(2,floor((TE-S8)/h))
                                       PA = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)), LP-x, LC-x);
                                       PB = @(x) integral(@(x1) (gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),LC-x,inf);
                                       % 外部积分
                                       x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                       result1 = zeros(size(x_range));
                                       for i = 1:length(x_range)
                                         result1(i) = PA(x_range(i));
                                       end
                                       result2 = zeros(size(x_range));
                                       for i = 1:length(x_range)
                                         result2(i) = PB(x_range(i));
                                       end
                                        % 数值积分结果
                                       PA = trapz(x_range, result1);
                                       PB = trapz(x_range, result2);
                                       C9P=(PA*((T9-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T9-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C9P;
                                     end
                                     S9=T1+T2+T3+T4+T5+T6+T7+T8+T9;
                                     n=9;
                                     PN9 = @(x) integral(@(x1) (gampdf(x, alpha*(max(1,floor((TE-S9)/h)))*h, beta)*gampdf(x1, alpha*max(1,(TE-S9-S/q1-(floor((TE-S9)/h))*h)), beta)), 0, LC-x);
                                     x_range = 0:0.01:LP;  % 定义 x 的范围和步长
                                     result3 = zeros(size(x_range));
                                     for i = 1:length(x_range)
                                         result3(i) = PN9(x_range(i));
                                     end
                                     PN9 = trapz(x_range, result3);
%                                    PN9
                                     CIS9= floor((TE-S9)/h)*Cins;
                                     CPS9=(q1*(TE-S9-S/(q2-q1)-S/q1)+q2*(S/(q2-q1)))*Cunit;
                                     CHS9=(S*(TE-S9)+S*S/2/(q2-q1)-S*S/2/q1)*Ch;
                                     EC=EC+PN9*(CIS9+CPS9+CHS9+C9I+C9H+C9S+C9M+C9P);
                                  end
                               end
                           end 
                       end
                   end
               end
           end
    end
end
n
EC=double(EC);
y=EC;
y
end