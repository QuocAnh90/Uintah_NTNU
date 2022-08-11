x = 0:0.1:1;

% Unit m s Pa N

E = 100000000;
Ev = E*(1-0.3)/(1+0.3)/(1-2*0.3);

n = 0.3;
d = 0.0001;
mu = 0.001;

% Hydraulic conductivity
K = n*n*n*d*d/180/(1-n)/(1-n)/mu;

Cv = K*Ev;
t = 7.5;
Tv = Cv * t/1/1;

P=zeros(11,1);

for pos = 1:11
  for i=1:100
     M_i = (i-0.5)*pi;
     P(pos) = P(pos) + 2*10000/M_i*sin(M_i*x(pos)/1)*exp(-M_i*M_i*Tv);
  end
end

  figure 
  plot (x,P)