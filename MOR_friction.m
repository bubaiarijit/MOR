clc
clear all
close all

%% System Description 
load('mechanical_example_data.mat');
clear n
n = 40; %dimension of original system
save('mechanical_example_data.mat')
global n
tspan=0:.01:100;
Nr1=40; %dimension of reduced order model
Nr2=30;
K1=K(1:n,1:n);
M1=M(1:n,1:n);
M_inv=inv(M1);
alpha=10^-20;
epsilon_o=1000   ;
epsilon_c=5000;
global beta
beta = [];
for i=1:n
    beta = [beta 0.2*10^(-10)*M1(i,i)*9.81]; 
end
G=[0;1;zeros(n-2,1)];
global M1
global K1
global G
B = [zeros(n,1); G];
H = [K1 zeros(n,n);zeros(n,n) M_inv];
%% Generalized differential controllability Gramian
P=sdpvar(2*n);
A=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n,n)]*H;
JR=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n,n)];
F=[];
F=[F, 
    A*P+P*A'+B*B'<=-10^(-3)*P
    ];
for i=1:n
    r=0.1*ones(n,1);
    r(i)= r(i)+(beta(i)/(sqrt(alpha))* M1(i,i));
    A=[zeros(n,n) eye(n,n);-eye(n,n) -diag(r)]*H;
    F=[F,
        A*P+P*A'+B*B'<=-10^(-3)*P
        ];
end
F=[F, 
    P>=10^(-10)*eye(2*n)];
options=sdpsettings('solver','sedumi');
sol=optimize(F,trace(P),options);
if sol.problem == 0
disp('Feasible')
elseif sol.problem == 1
disp('Infeasible')
else
disp('Something else happened')    
end

P=value(P);
%% Linear controllability Gramian around (q,p) = (0,0)
P2=sdpvar(2*n);
r=0.1*ones(n,1);
for i=1:n
    r(i)= r(i)+(beta(i)/(sqrt(alpha))* M1(i,i));
end    
A=[zeros(n,n) eye(n,n);-eye(n,n) -diag(r)]*H;
F=[];
F=[F, 
    A*P2+P2*A'+B*B'<=-10^(-3)*P2
    ];
F=[F, 
    P2>=10^(-10)*eye(2*n)];
options=sdpsettings('solver','sedumi');
sol=optimize(F,trace(P2),options);
if sol.problem == 0
disp('Feasible')
elseif sol.problem == 1
disp('Infeasible')
else
disp('Something else happened')    
end

P2=value(P2);

%% Generalized differential observability Gramian
Q=sdpvar(2*n);
A=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n)]*H;
JR=[zeros(n,n) eye(n,n);-eye(n,n) -0.1*eye(n)];

F=[];
F=[F,
    Q*A+A'*Q+H*B*B'*H<=-10^(-3)*Q
    ];
for i=1:n
    r=0.1*ones(n,1);
    r(i)= r(i)+(beta(i)/(sqrt(alpha))* M1(i,i));
    A=[zeros(n,n) eye(n,n);-eye(n,n) -diag(r)]*H;
    F=[F,
         Q*A+A'*Q+H*B*B'*H<=-10^(-3)*Q
         ];
end
F=[F, 
    Q>=10^(-10)*eye(2*n)];
options=sdpsettings('solver','sedumi');
sol=optimize(F,trace(Q),options);
if sol.problem == 0
disp('Feasible')
elseif sol.problem == 1
disp('Infeasible')
else
disp('Something else happened')    
end

Q=value(Q);
%% Linear observability Gramian around (q,p) = (0,0)
Q2=sdpvar(2*n);
r=0.1*ones(n,1);
for i=1:n
    r(i)= r(i)+(beta(i)/(sqrt(alpha))* M1(i,i));
end    
A=[zeros(n,n) eye(n,n);-eye(n,n) -diag(r)]*H;
F=[];
F=[F, 
    A*Q2+Q2*A'+B*B'<=-10^(-3)*Q2
    ];
F=[F, 
    Q2>=10^(-10)*eye(2*n)];
options=sdpsettings('solver','sedumi');
sol=optimize(F,trace(Q2),options);
if sol.problem == 0
disp('Feasible')
elseif sol.problem == 1
disp('Infeasible')
else
disp('Something else happened')    
end

Q2=value(Q2);

%% Extended differential Gramians
Pchk=inv(P);
Gamma_o=epsilon_o*Q;
Gamma_c=-epsilon_c*P;
S=Q*inv(6000*Q+Gamma_o)*Q;
Ti=(6000*P+Gamma_c);
%% Extended differential balancing transformation

PhTi=chol(Ti);
[UTS,S2TS]=svd(PhTi*S*PhTi');
STS=sqrt(S2TS);
global We
We=PhTi'*UTS*sqrt(inv(STS));

%% Generalized differential balancing transformation
PhP=chol(P);
[UPQ,Q2PQ]=svd(PhP*Q*PhP');
QPQ=sqrt(Q2PQ);
global W
W=PhP'*UPQ*sqrt(inv(QPQ));
%% Linear balancing transformation
P2hP2=chol(P2);
[U2P2Q,Q2P2Q]=svd(P2hP2*Q2*P2hP2');
Q2P2Q2=sqrt(Q2P2Q);
global W2
W2=P2hP2'*U2P2Q*sqrt(inv(Q2P2Q2));

%% Figure 1 (Comparison of outputs of original systems and reduced models)
z0=zeros(2*n,1);

Cbal2=B'*H*We;
[t,zeg]=ode45(@(t,z) func(t,z,We), tspan, z0);
yeg=[];
for i=1:length(zeg(:,1))
    yeg=[yeg;Cbal2*zeg(i,:)'];
end
figure(1)
plot(t,yeg,'Linewidth',3)
set(gca,'FontSize',40)
xlabel('time', 'fontsize',50)
ylabel('output', 'fontsize',50)
hold on
legend('Original System','FontSize',30)
hold on
[t,zreg]=ode45(@(t,z) funcred(t,z,We,Nr1), tspan, z0);
yreg=[];
for i=1:length(zreg(:,1))
    yreg=[yreg;Cbal2*zreg(i,:)'];
end    
plot(t,yreg,'--','DisplayName',strcat('Reduced system(EDBT), dim=',num2str(Nr1)),'Linewidth',3)
hold on
[t,zreg]=ode45(@(t,z) funcred(t,z,We,Nr2), tspan, z0);
yreg=[];
for i=1:length(zreg(:,1))
    yreg=[yreg;Cbal2*zreg(i,:)'];
end    
plot(t,yreg,'-.','DisplayName',strcat('Reduced system(EDBT), dim=',num2str(Nr2)),'Linewidth',3)
hold on
Cbal3=B'*H*W2;
[t,zrlg]=ode45(@(t,z) funcred(t,z,W2,Nr1), tspan, z0);
yrlg=[];
for i=1:length(zrlg(:,1))
    yrlg=[yrlg;Cbal3*zrlg(i,:)'];
end    
plot(t,yrlg,'--','DisplayName',strcat('Reduced system(LBT), dim=',num2str(Nr1)),'Linewidth',3)
hold on
Cbal3=B'*H*W2;
[t,zrlg]=ode45(@(t,z) funcred(t,z,W2,Nr2), tspan, z0);
yrlg=[];
for i=1:length(zrlg(:,1))
    yrlg=[yrlg;Cbal3*zrlg(i,:)'];
end    
plot(t,yrlg,'--','DisplayName',strcat('Reduced system(LBT), dim=',num2str(Nr2)),'Linewidth',3)
hold off
%% Figure 2 (Error Bounds)
figure(3)
p=[];
bound=[];
for i=1:2*n
    p=[p i];
    SumSTS=0;
    SumSTS=trace(STS(i+1:2*n,i+1:2*n));
    bound=[bound 2*SumSTS*sqrt(0.5)];
end
plot(p,bound,'Linewidth',3.5)
set(gca,'FontSize',40)
xlabel('Dimension of reduced model','fontsize',50)
ylabel('Upper bound for error','fontsize',50)
legend('EDBT', 'FontSize',50)
hold on
q=[];
bound=[];
for i=1:2*n
    q=[q i];
    SumQPQ=0;
    SumQPQ=trace(QPQ(i+1:2*n,i+1:2*n));
    bound=[bound 2*SumQPQ*sqrt(0.5)];
end
plot(q,bound,'DisplayName','GDBT','Linewidth',3.5)
hold off

%% Figure 3(Difference between actual errors incurred by EDBT and GDBT)
error1=[];
error2=[];
error3=[];
for i=1:length(t(:,1))
    error1=[error1;yeg(i)-yreg(i)];
    error2=[error2;yeg(i)-yrg(i)];
    error3=[error3;error2(i)-error1(i)];
end
figure(4)
subplot(2,1,1)
plot(t,error1,'Linewidth',3)
title('Errors')
set(gca,'FontSize',40)
xlabel('time','fontsize',50)
ylabel('Output error','fontsize',50)
legend(strcat('EDBT, dim=',num2str(Nr1)), 'FontSize',30)
hold on
plot(t,error2,'--','DisplayName',strcat('GDBT, dim=',num2str(Nr1)),'Linewidth',3)
hold off
subplot(2,1,2)
plot(t,error3,'Linewidth',2.5)
set(gca,'FontSize',40)
xlabel('time','fontsize',50)
ylabel('error diference','fontsize',50)
title('Difference between errors')

%% Function to simulate original system
function zdot = func(t,z,W)
alpha=10^-20;
global n
global K1
global M1
global beta
M_inv=inv(M1);
global G
    u = sin(t);
    q = z(1:n);
    p = z(n:2*n);
    x = W*z;
    r=[];
    for k = n+1:2*n
        r = [r beta(k-n)/sqrt(alpha+x(k)^2/M_inv(k-n,k-n)^2)];
    end
    R = diag(r);
    zdot = inv(W)*[zeros(n,n) eye(n,n);-eye(n,n) -(0.1*eye(n)+R)]*[K1 zeros(n,n);zeros(n,n) M_inv]*W*z + inv(W)*[zeros(n,1); G]*u;
end  
%% Function to simulate reduced order model
function zrdot = funcred(t,z,W,Nr)
alpha=10^-20;
global n
global K1
global M1
global beta
M_inv=inv(M1);
global G
    u = sin(t);
    z(Nr+1:2*n)=zeros(2*n-Nr,1);
    x = W*z;
    r=[];
    for k = n+1:2*n
        r = [r beta(k-n)/sqrt(alpha+x(k)^2/M_inv(k-n,k-n)^2)];
    end   
    R = diag(r);
    Abal=inv(W)*[zeros(n,n) eye(n,n);-eye(n,n) -(0.1*eye(n)+R)]*[K1 zeros(n,n);zeros(n,n) M_inv]*W;
    Abal(Nr+1:2*n,:)=zeros(2*n-Nr,2*n);
    Ar=Abal;
    Bbal=inv(W)*[zeros(n,1); G];
    Bbal(Nr+1:2*n)=zeros(2*n-Nr,1);
    Br=Bbal;
    zrdot = Ar*z + Br*u;
end  
