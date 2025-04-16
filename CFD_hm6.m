%%  参数设置
clear;
gamma = 1.4;       % 绝热指数
x_left = -0.5;     % 计算域左边界
x_right = 0.5;     % 计算域右边界
t_end = 0.25;      % 计算终止时间
CFL = 0.8;         % CFL数
Mx = 200;           % 网格数
delta_x = (x_right - x_left)/Mx; % 空间步长
Nt = 200;
delta_t = t_end/Nt;
U = zeros(3, Mx + 5); % 3行（ρ, ρu, E），列对应单元
F = zeros(3, Mx + 5);

% 填充实际单元初始值
for i = 1:Mx+5
    [rho,u,p] = InitialConditions2(x_left+(i-3)*delta_x);
    U(1,i) = rho;
    U(2,i) = rho*u;
    U(3,i) = rho*(p/(rho*(gamma-1))+0.5*u^2);
    A = [0 1 0;
        (gamma-3)*0.5*u^2 (3-gamma)*u gamma-1;
        (gamma-1)*u^3-gamma*u*(p/(rho*(gamma-1))+0.5*u^2) -3/2*(gamma-1)*u^2+gamma*(p/(rho*(gamma-1))+0.5*u^2) gamma*u];
    F(:,i) = A*U(:,i);
end

%% Rusanov格式
M = zeros(3, Mx + 5); 
M(:,1) = U(:,1);
M(:,2) = U(:,2);
M(:,3) = U(:,3);
M(:,Mx+3) = U(:,Mx+3);
M(:,Mx+4) = U(:,Mx+4);
M(:,Mx+5) = U(:,Mx+5);
for t=1:Nt
    for x=4:Mx+2
        rho = U(1,x);
        u = U(2,x)/U(1,x);
        p = rho*(gamma-1)*(U(3,x)/U(1,x)-0.5*u^2);
        F12 = [rho*u;rho*u^2+p;u*(U(3,x)+p)];

        rho_l = U(1,x-1);
        u_l = U(2,x-1)/U(1,x-1);
        p_l = rho_l*(gamma-1)*(U(3,x-1)/U(1,x-1)-0.5*u_l^2);
        F11 = [rho_l*u_l;rho_l*u_l^2+p_l;u_l*(U(3,x-1)+p_l)];
        F21 = F12;

        rho_m = U(1,x+1);
        u_m = U(2,x+1)/U(1,x+1);
        p_m = rho_m*(gamma-1)*(U(3,x+1)/U(1,x+1)-0.5*u_m^2);
        F22 = [rho_m*u_m;rho_m*u_m^2+p_m;u_m*(U(3,x+1)+p_m)];

        a_l = sqrt(gamma*p_l/rho_l);
        a_m = sqrt(gamma*p_m/rho_m);
        a = sqrt(gamma*p/rho);
        lamda1 = 0.5*(abs(u_l)+abs(u)+a_l+a);
        lamda2 = 0.5*(abs(u)+abs(u_m)+a+a_m);

        F1 = 0.5*(F11+F12) - lamda1*0.5*(U(:,x)-U(:,x-1));
        F2 = 0.5*(F21+F22) - lamda2*0.5*(U(:,x+1)-U(:,x));%两个数值通量
        M(:,x) = delta_t/delta_x*(F1-F2)+U(:,x);
    end
    U = M;
end

%% Jameson格式
M = zeros(3, Mx + 5); 
M(:,1) = U(:,1);
M(:,2) = U(:,2);
M(:,3) = U(:,3);
M(:,Mx+3) = U(:,Mx+3);
M(:,Mx+4) = U(:,Mx+4);
M(:,Mx+5) = U(:,Mx+5);
k2 = 0.6;
k4 = 1/64;
for t=1:Nt
    for x=4:Mx+2
        for j = x-2:x+3
            rho(j) = U(1,j);
            u(j) = U(2,j)/U(1,j);
            p(j) = rho(j)*(gamma-1)*(U(3,j)/U(1,j)-0.5*u(j)^2);
        end
        for j = x-2:x+1
            v1(j) = abs(p(j+1)-2*p(j)+p(j-1))/abs(p(j+1)+2*p(j)+p(j-1));
        end
        for j = x-1:x+2
            v2(j) = abs(p(j+1)-2*p(j)+p(j-1))/abs(p(j+1)+2*p(j)+p(j-1));
        end
        epsilon21 = k2*max(v1);
        epsilon41 = max(0,k4-epsilon21);
        epsilon22 = k2*max(v2);
        epsilon42 = max(0,k4-epsilon22);


        F11 = [rho(x-1)*u(x-1);rho(x-1)*u(x-1)^2+p(x-1);u(x-1)*(U(3,x-1)+p(x-1))];
        F12 = [rho(x)*u(x);rho(x)*u(x)^2+p(x);u(x)*(U(3,x)+p(x))];
        F21 = F12;
        F22 = [rho(x+1)*u(x+1);rho(x+1)*u(x+1)^2+p(x+1);u(x+1)*(U(3,x+1)+p(x+1))];

        a_l = sqrt(gamma*p(x-1)/rho(x-1));
        a_m = sqrt(gamma*p(x+1)/rho(x+1));
        a = sqrt(gamma*p(x)/rho(x));
        lamda1 = 0.5*(abs(u(x-1))+abs(u(x))+a_l+a);
        lamda2 = 0.5*(abs(u(x))+abs(u(x+1))+a+a_m);

        F1 = 0.5*(F11+F12) - lamda1*epsilon21*(U(:,x)-U(:,x-1)) + lamda1*epsilon41*(U(:,x+1)-3*U(:,x)+3*U(:,x-1)-U(:,x-2));
        F2 = 0.5*(F21+F22) - lamda2*epsilon22*(U(:,x+1)-U(:,x)) + lamda2*epsilon42*(U(:,x+2)-3*U(:,x+1)+3*U(:,x)-U(:,x-1));%两个数值通量
        M(:,x) = delta_t/delta_x*(F1-F2)+U(:,x);
    end
    U = M;
end


%% 作图
x = x_left:delta_x:x_right;
UU = zeros(3,Mx+1);
for i = 3:Mx+3
    UU(:,i-2) = U(:,i);
end
for i = 1:3
    plot(x,UU(i,:),'LineWidth',1.5);
    hold on;
end
legend('\rho','\rhou','\rhoE');
xlabel('x');
ylabel('U的各分量');
title(['Jameson格式(t = ',num2str(t_end),'s)']);
