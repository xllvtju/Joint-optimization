function y = Obj_funnew1(z)   

syms x x1  y1  y2  
% alpha = 0.4;  beta = 4;  LC=20; LP=15; Cins=3; CPM=20; CCM=40;  a1=0.5; b1=1;a2=2; b2=5;  q1=20; q2=30; Cunit=3; rou=2.5;Ch=0.05;
% h=3; S=100;TE=60;
alpha = 0.4;  beta = 4;  LC=20; LP=15; Cins=0.2; CPM=1; CCM=4;  a1=0.5; b1=1;a2=2; b2=5;  q1=20; q2=30; Cunit=0.3; rou=0.25;Ch=0.005;TE=60;
% h=3; S=100;
h=z(1); S=z(2);
u1(y1)=unifpdf(y1, a1, b1);
u2(y2)=unifpdf(y2, a2, b2);
ETPM=(a1+b1)/2;
ETCM=(a2+b2)/2;
% PA=int((int((gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),x1,LP-x,LC-x)),x,0,LP);
% PB=int((int((gampdf(x, alpha*(j-1)*h, beta)*gampdf(x1, alpha*h, beta)),x1,LC-x,inf)),x,0,LP);

T1=0;
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
  C1I=0;
  C1I=(PA*j+PB*(j+1))*Cins+C1I;
  C1M=0;
  C1M=PA*CPM+PB*CCM+C1M;
  C1H=0;
  C1H=(PA*((j*h*S)-S*S/2/(q2-q1)-S*S/2/q1)+PB*((((j-1)*h+PB*h)*S)-S*S/2/(q2-q1)-S*S/2/q1))*Ch+C1H;
  C1S=0;
  C1S=(PA*max(0,int((q1*y1-S)*u1(y1),y1,S/q1,b1))*rou+PB*max(0,int((q1*y2-S)*u2(y2),y2,S/q1,b2)))*rou+C1S;
end
T1=double(T1);
C1P=0;
for j=2:max(2,floor((TE)/h))
   C1P=(PA*((T1-max(ETPM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1))+PB*((T1-max(ETCM,S/q1)-S/(q2-q1))*q1+S*q2/(q2-q1)))*Cunit+C1P;
end
EC=(C1I+C1H+C1S+C1M+C1P);
EC=EC/T1;
EC=double(EC);
y=EC;
T1
y
