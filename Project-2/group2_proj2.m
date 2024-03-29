close all;
foo=[4.4, 8.0, 2, 120, -84, -60, 0.02, -1.2, 18, 2, 30, 2, 30, 20];
foocell=num2cell(foo);
[ gca, gk, gl, vca, vk, vl, phi, V1, V2, V3, V4, V5, V6, C]=foocell{:};


syms V W 
dvdt= -gca * m_inf(V,V1,V2) * (V-vca) ./ C - gk * W * (V-vk) ./C - gl * (V-vl)./C ; %+ I./C; (I=0) 
dwdt= phi * (w_inf(V,V3,V4)- W)/ tau_w(V,V3,V4) ;

sol = solve([dvdt, dwdt], [V, W]); %Solve for equilibrium point
VSol = sol.V;
WSol = sol.W;
fprintf("Equilibrium point of the system is at V = %f mV and w = %f",VSol,WSol);

figure %Plotting nullclines
fimplicit(dvdt, [-100 100 -0.2 1])
hold on
fimplicit(dwdt, [-100 100 -0.2 1], '-r')
hold on
title("Nullclines plotted for Morris-Lecar equations");
xlabel("V in milivolts");
ylabel("w");
hold on

my_phase(0, 0.02, 0, WSol) 
hold off
legend('V-nullcline','W-nullcline','Phase plane plot')
jacob_matrix=jacobian([dvdt,dwdt],[V,W]);
jacob_eqbm=subs(jacob_matrix, {V, W}, {VSol, WSol});
eigenvalue=eig(jacob_eqbm)


iter=20;
[V0, W0]=meshgrid(-70:((100+1)/iter ):30 , 0:(1/(iter-1) ):1);
U=zeros(iter,iter);
Vq=zeros(iter,iter);
for i=1:iter*iter
U(i)=phi*(0.5*(1+tanh((V0(i)-V3)/V4))-W0(i))*cosh((V0(i)-V3)/V4);
Vq(i)=(0- gca*0.5*(1 + tanh((V0(i)-V1)/V2))*(V0(i)-vca)-gk*W0(i)*(V0(i)-vk)-gl*(V0(i)-V1))/C; %dV/dt
end
figure
quiver(V0,W0*100,Vq,U*100);
title('Quiver plot for V vs w');
xlabel('V');
ylabel('w*100');

%%%PART 5

figure %Plotting for different phi values
fimplicit(dvdt, [-100 100 -0.2 1])
hold on
fimplicit(dwdt, [-100 100 -0.2 1], '-r')
hold on
my_phase(0,0.02,-10,WSol); 
hold on
my_phase(0,0.04,-10,WSol);
hold on
my_phase(0,0.01,-10,WSol);
title('Phase plane plots for different phi values')
xlabel('V');
ylabel('W');
legend('V-nullcline','W-nullcline','phi=0.02','phi=0.04','phi=0.01')
hold off


%PART6
size=50;
peak_volt=zeros(size,1);
init_volt=zeros(size,1);
i=0;
start=-20; %%%Major jump is found to be in this region by trial and error
ending=-10;
for V_iter=start:(ending-start)/size:ending
    i=i+1;
    peak_volt(i)=my_phase_2(0, 0.02, V_iter, WSol);
    init_volt(i)=V_iter;
end
figure
plot(init_volt,peak_volt)
title("Peak voltage vs initial voltage");
xlabel("Initial voltage");
ylabel("Peak voltage of action potential");


%%%%%%PART 7 AND 8
syms V_new2 W_new2 
dvdt_new2= -gca * m_inf(V_new2,V1,V2) * (V_new2-vca) ./ C - gk * W_new2 * (V_new2-vk) ./C - gl * (V_new2-vl)./C + 86/C; %Added 86uA/cm2 Iext 
dwdt_new2= phi * (w_inf(V_new2,V3,V4)- W_new2)/ tau_w(V_new2,V3,V4) ;

sol = solve([dvdt_new2, dwdt_new2], [V_new2, W_new2]); %Solve for equilibrium point
VSol_new2 = sol.V_new2;
WSol_new2 = sol.W_new2;
fprintf("Equilibrium point of the system is at V = %f mV and w = %f",VSol_new2,WSol_new2);

figure %Plotting nullclines
fimplicit(dvdt_new2, [-100 100 -0.2 0.5])
hold on
fimplicit(dwdt_new2, [-100 100 -0.2 0.5], '-r')
hold on
title("Phase plane plots with different initial conditions with Iext=86mA");
xlabel("V in milivolts");
ylabel("w");
hold on
my_phase(86,0.02,VSol,WSol);
hold on
my_phase(86,0.02,VSol_new2,WSol_new2);
hold on
my_phase(86,0.02,-27.9,0.17);
hold on
my_phase_negtime(86,0.02,-27.9,0.17);
legend('V-nullcline','w-nullcline','Case 1','Case 2','Case 3','Part 8 (-ve time)');
hold off



jacob_matrix_new2=jacobian([dvdt_new2,dwdt_new2],[V_new2,W_new2]);
jacob_eqbm_new2=subs(jacob_matrix_new2, {V_new2, W_new2}, {VSol_new2, WSol_new2});
eigenvalue_new2=eig(jacob_eqbm_new2)



%%%PART9

Iext_pt9_1=[80 86 90];
for i=1:3
Iext_new3=Iext_pt9_1(i);
syms V_new3 W_new3 
dvdt_new3= -gca * m_inf(V_new3,V1,V2) * (V_new3-vca) ./ C - gk * W_new3 * (V_new3-vk) ./C - gl * (V_new3-vl)./C + Iext_new3/C;  
dwdt_new3= phi * (w_inf(V_new3,V3,V4)- W_new3)/ tau_w(V_new3,V3,V4) ;

sol = solve([dvdt_new3, dwdt_new3], [V_new3, W_new3]); %Solve for equilibrium point
VSol_new3 = sol.V_new3;
WSol_new3 = sol.W_new3;
fprintf("Equilibrium point of the system for Iext=%f is at V = %f mV and w = %f",Iext_new3,VSol_new3,WSol_new3);
jacob_matrix_new3=jacobian([dvdt_new3,dwdt_new3],[V_new3,W_new3]);
jacob_eqbm_new3=subs(jacob_matrix_new3, {V_new3, W_new3}, {VSol_new3, WSol_new3});
eigenvalue_new3=eig(jacob_eqbm_new3) %%On printing this, eigenvalue changes sign after sometime (although complex)
end



firing=zeros(31);

Iext_pt9=80:1:110;
for i=1:31
Iext_new3=Iext_pt9(i);
syms V_new3 W_new3 
dvdt_new3= -gca * m_inf(V_new3,V1,V2) * (V_new3-vca) ./ C - gk * W_new3 * (V_new3-vk) ./C - gl * (V_new3-vl)./C + Iext_new3/C;  
dwdt_new3= phi * (w_inf(V_new3,V3,V4)- W_new3)/ tau_w(V_new3,V3,V4) ;

sol = solve([dvdt_new3, dwdt_new3], [V_new3, W_new3]); %Solve for equilibrium point
VSol_new3 = sol.V_new3;
WSol_new3 = sol.W_new3;
firing(i)=my_phase_firing(Iext_new3,0.02,VSol_new3,WSol_new3); %THIS returns number of spikes based on sign changes in V-t plot
end
figure
plot(Iext_pt9,firing);
title("Firing rate vs Iext");
xlabel("Iext");
ylabel("Firing rate / second");



%%%%%%%%%%PART 10

disp("START OF PART 10");
foo=[4, 8.0, 2, 120, -84, -60, 0.0667, -1.2, 18, 12, 17.4, 12, 17.4, 20];
foocell=num2cell(foo);
[ gca, gk, gl, vca, vk, vl, phi, V1, V2, V3, V4, V5, V6, C]=foocell{:};

syms V W
Iext_part10=30;
dvdt=-gca * m_inf(V,V1,V2) * (V-vca) ./ C - gk * W * (V-vk) ./C - gl * (V-vl)./C +Iext_part10./C;
dwdt=phi * (w_inf(V,V3,V4)- W)/ tau_w(V,V3,V4);
eqns=[-gca * m_inf(V,V1,V2) * (V-vca) ./ C - gk * W * (V-vk) ./C - gl * (V-vl)./C +Iext_part10./C==0,phi * (w_inf(V,V3,V4)- W)/ tau_w(V,V3,V4)==0 ];
vars=[V W];
[VSol_part10,WSol_part10]=solve(eqns,vars); %Solve for equilibrium point%%%ONLY 1 found although there are 3
%fprintf("Equilibrium point of the system is at V = %f mV and w = %f",VSol_part10,WSol_part10);

figure %Plotting nullclines
nullclineV=@(V,W)dvdt;
fimplicit(dvdt, [-100 100 -0.2 1])
hold on
nullclineW=@(V,W)dwdt;
fimplicit(dwdt, [-100 100 -0.2 1], '-r')
hold on

hold on

jacob_matrix=jacobian([dvdt,dwdt],[V,W]);
V_part10_1= -41.85;  %Values manually found from graph for 3 eqbm points as only one solution obtained from solve function
V_part10_2= -19.56;
V_part10_3= 3.874;
W_part10_1= 0.002047;
W_part10_2= 0.02588;
W_part10_3= 0.2821;


jacob_eqbm=subs(jacob_matrix, {V, W}, {V_part10_1, W_part10_1});
eigenvalue=double(eig(jacob_eqbm)) %BOTH NEGATIVE AND REAL

jacob_eqbm=subs(jacob_matrix, {V, W}, {V_part10_2, W_part10_2});
eigenvalue=double(eig(jacob_eqbm)) %BOTH REAL, 1 +ve, 1 -ve

jacob_eqbm=subs(jacob_matrix, {V, W}, {V_part10_3, W_part10_3});
eigenvalue=double(eig(jacob_eqbm))%BOTH COMPLEX CONJ with +ve real parts

%%%Manifolds plotted about saddle point
manifolds(30,0.02,V_part10_2-0.001,W_part10_2-0.001,0,1000,'y');
hold off

title("Nullclines and manifolds plotted for Morris-Lecar equations");
xlabel("V in milivolts");
ylabel("w");
legend('V-nullcline','W-nullcline','Manifold about saddle point');



%%%%PART 11
I_p11=[30,39,39.5,40,45,50];
figure
for i=1:6
Iext_part10=I_p11(i);
foo=[4, 8.0, 2, 120, -84, -60, 0.0667, -1.2, 18, 12, 17.4, 12, 17.4, 20];
foocell=num2cell(foo);
[ gca, gk, gl, vca, vk, vl, phi, V1, V2, V3, V4, V5, V6, C]=foocell{:};
syms V W
dvdt=-gca * m_inf(V,V1,V2) * (V-vca) ./ C - gk * W * (V-vk) ./C - gl * (V-vl)./C +Iext_part10/C;
dwdt=phi * (w_inf(V,V3,V4)- W)/ tau_w(V,V3,V4);
 
eqns=[-gca * m_inf(V,V1,V2) * (V-vca) ./ C - gk * W * (V-vk) ./C - gl * (V-vl)./C +Iext_part10/C==0,phi * (w_inf(V,V3,V4)- W)/ tau_w(V,V3,V4)==0 ];
vars=[V W];
[VSol_part10,WSol_part10]=solve(eqns,vars); 
nullclineV=@(V,W)dvdt;
fimplicit(dvdt, [-100 100 -0.2 0.5])
hold on
if(i==1)
    nullclineW=@(V,W)dwdt;
    fimplicit(dwdt, [-100 100 -0.2 0.5], '-r')
    hold on
end
title("Nullclines plotted for Morris-Lecar equations (PART 11)");
xlabel("V in milivolts");
ylabel("w");
legend('V-nullcline(Iext=30mA)','W-nullcline','V-nullcline(Iext=39mA)','V-nullcline(Iext=39.5mA)','V-nullcline(Iext=40mA)','V-nullcline(Iext=45mA)','V-nullcline(Iext=50mA)');

hold on
jacob_matrix=jacobian([dvdt,dwdt],[V,W]);
jacob_eqbm=subs(jacob_matrix, {V, W}, {VSol_part10,WSol_part10});
eigenvalue=double(eig(jacob_eqbm)); %BOTH NEGATIVE AND REAL
end



function manifolds(Iext,phi, VSol, WSol, t_start,t_end,clr)
[~,X] = ode15s(@EOM,[t_start t_end],[double(VSol); double(WSol)]);
Vnew = X(:,1);
Wnew = X(:,2);
plot(Vnew,Wnew,clr)
    function dX = EOM(t,y) %Nested function to ensure value of Iext is passed
    Vnew  = y(1);
    Wnew  = y(2);
    foo=[4.4, 8.0, 2, 120, -84, -60, -1.2, 18, 2, 30, 2, 30, 20]; %phi taken as parameter of function separately
    foocell=num2cell(foo);
    [ gca, gk, gl, vca, vk, vl, V1, V2, V3, V4, V5, V6, C]=foocell{:};
    dvdt_I=-1.*gca .* m_inf(Vnew,V1,V2) .* (Vnew-vca) ./ C - gk .* Wnew .* (Vnew-vk) ./C - gl .* (Vnew-vl)./C +Iext ./C;
    dwdt_I=phi .* (w_inf(Vnew,V3,V4)- Wnew)./ tau_w(Vnew,V3,V4);
    dX = [dvdt_I; dwdt_I];
    end
end

function my_phase(Iext,phi, VSol, WSol)
%[~,X] = ode45(@EOM,[-500 500],[-0.2 1]);
[~,X] = ode15s(@EOM,[0 300],[double(VSol); double(WSol)]);
Vnew = X(:,1);
Wnew = X(:,2);
plot(Vnew,Wnew)
    function dX = EOM(t,y) %Nested function to ensure value of Iext is passed
    Vnew  = y(1);
    Wnew  = y(2);
    foo=[4.4, 8.0, 2, 120, -84, -60, -1.2, 18, 2, 30, 2, 30, 20]; %phi taken as parameter of function separately
    foocell=num2cell(foo);
    [ gca, gk, gl, vca, vk, vl, V1, V2, V3, V4, V5, V6, C]=foocell{:};
    dvdt_I=-1.*gca .* m_inf(Vnew,V1,V2) .* (Vnew-vca) ./ C - gk .* Wnew .* (Vnew-vk) ./C - gl .* (Vnew-vl)./C +Iext ./C;
    dwdt_I=phi .* (w_inf(Vnew,V3,V4)- Wnew)./ tau_w(Vnew,V3,V4);
    dX = [dvdt_I; dwdt_I];
    end
end

function count=my_phase_firing(Iext,phi, V1, W1)
[~,X] = ode45(@EOM,[0 1000],[double(V1); double(W1)]);
Vnew = X(:,1);
Wnew = X(:,2);
count=0;
for i=1:size(Vnew)-1
    if(Vnew(i)<0 && Vnew(i+1)>=0)
        count=count+1;
    end
end
    function dX = EOM(t,y) %Nested function to ensure value of Iext is passed
    Vnew  = y(1);
    Wnew  = y(2);
    foo=[4.4, 8.0, 2, 120, -84, -60, -1.2, 18, 2, 30, 2, 30, 20]; %phi taken as parameter of function separately
    foocell=num2cell(foo);
    [ gca, gk, gl, vca, vk, vl, V1, V2, V3, V4, V5, V6, C]=foocell{:};
    dvdt_I=-1.*gca .* m_inf(Vnew,V1,V2) .* (Vnew-vca) ./ C - gk .* Wnew .* (Vnew-vk) ./C - gl .* (Vnew-vl)./C +Iext ./C;
    dwdt_I=phi .* (w_inf(Vnew,V3,V4)- Wnew)./ tau_w(Vnew,V3,V4);
    dX = [dvdt_I; dwdt_I];
    end
end

function my_phase_negtime(Iext,phi, VSol, WSol)
%[~,X] = ode45(@EOM,[-500 500],[-0.2 1]);
[~,X] = ode15s(@EOM,[0 -300],[double(VSol); double(WSol)]);
Vnew = X(:,1);
Wnew = X(:,2);
plot(Vnew,Wnew)
    function dX = EOM(t,y) %Nested function to ensure value of Iext is passed
    Vnew  = y(1);
    Wnew  = y(2);
    foo=[4.4, 8.0, 2, 120, -84, -60, -1.2, 18, 2, 30, 2, 30, 20]; %phi taken as parameter of function separately
    foocell=num2cell(foo);
    [ gca, gk, gl, vca, vk, vl, V1, V2, V3, V4, V5, V6, C]=foocell{:};
    dvdt_I=-1.*gca .* m_inf(Vnew,V1,V2) .* (Vnew-vca) ./ C - gk .* Wnew .* (Vnew-vk) ./C - gl .* (Vnew-vl)./C +Iext ./C;
    dwdt_I=phi .* (w_inf(Vnew,V3,V4)- Wnew)./ tau_w(Vnew,V3,V4);
    dX = [dvdt_I; dwdt_I];
    end
end

function peak_V=my_phase_2(Iext,phi, V1, W1)
%[~,X] = ode45(@EOM,[-500 500],[-0.2 1]);
[~,X] = ode45(@EOM,[0 1000],[double(V1); double(W1)]);
Vnew = X(:,1);
Wnew = X(:,2);
peak_V=max(Vnew);
% figure
% plot(Vnew)
% title("V vs t plot");
% xlabel("t");
% ylabel("V");
% figure
% plot(Vnew,Wnew,'-g')
% xlabel('V')
% ylabel('W')
% grid

    function dX = EOM(t,y) %Nested function to ensure value of Iext is passed
    Vnew  = y(1);
    Wnew  = y(2);
    foo=[4.4, 8.0, 2, 120, -84, -60, -1.2, 18, 2, 30, 2, 30, 20]; %phi taken as parameter of function separately
    foocell=num2cell(foo);
    [ gca, gk, gl, vca, vk, vl, V1, V2, V3, V4, V5, V6, C]=foocell{:};
    dvdt_I=-1.*gca .* m_inf(Vnew,V1,V2) .* (Vnew-vca) ./ C - gk .* Wnew .* (Vnew-vk) ./C - gl .* (Vnew-vl)./C +Iext ./C;
    dwdt_I=phi .* (w_inf(Vnew,V3,V4)- Wnew)./ tau_w(Vnew,V3,V4);
    dX = [dvdt_I; dwdt_I];
    end
end

function peak_V=my_phase_2_withplot(Iext,phi, V1, W1)
%[~,X] = ode45(@EOM,[-500 500],[-0.2 1]);
[~,X] = ode45(@EOM,[0 1000],[double(V1); double(W1)]);
Vnew = X(:,1);
Wnew = X(:,2);
peak_V=max(Vnew);
t=0:1000/size(Vnew,1):1000-1000/size(Vnew,1);
figure
plot(t,Vnew)
title("V vs t plot");
xlabel("t");
ylabel("V");
figure
plot(Vnew,Wnew,'-g')
xlabel('V')
ylabel('W')
grid

    function dX = EOM(t,y) %Nested function to ensure value of Iext is passed
    Vnew  = y(1);
    Wnew  = y(2);
    foo=[4.4, 8.0, 2, 120, -84, -60, -1.2, 18, 2, 30, 2, 30, 20]; %phi taken as parameter of function separately
    foocell=num2cell(foo);
    [ gca, gk, gl, vca, vk, vl, V1, V2, V3, V4, V5, V6, C]=foocell{:};
    dvdt_I=-1.*gca .* m_inf(Vnew,V1,V2) .* (Vnew-vca) ./ C - gk .* Wnew .* (Vnew-vk) ./C - gl .* (Vnew-vl)./C +Iext ./C;
    dwdt_I=phi .* (w_inf(Vnew,V3,V4)- Wnew)./ tau_w(Vnew,V3,V4);
    dX = [dvdt_I; dwdt_I];
    end
end

function peak_V=my_phase_part10(Iext,phi, V1, W1)
%[~,X] = ode45(@EOM,[-500 500],[-0.2 1]);
[~,X] = ode45(@EOM,[0 2000],[double(V1); double(W1)]);
Vnew = X(:,1);
Wnew = X(:,2);
peak_V=max(Vnew);
t=0:2000/size(Vnew,1):2000-2000/size(Vnew,1);
figure
plot(t,Vnew)
title("V vs t plot");
xlabel("t");
ylabel("V");
figure
plot(Vnew,Wnew,'-g')
xlabel('V')
ylabel('W')
grid

    function dX = EOM(t,y) %Nested function to ensure value of Iext is passed
    Vnew  = y(1);
    Wnew  = y(2);
    foo=[4, 8.0, 2, 120, -84, -60,  -1.2, 18, 12,17.4, 12, 17.4, 20]; %phi taken as parameter of function separately
    foocell=num2cell(foo);
    [ gca, gk, gl, vca, vk, vl, V1, V2, V3, V4, V5, V6, C]=foocell{:};
    dvdt_I=-1.*gca .* m_inf(Vnew,V1,V2) .* (Vnew-vca) ./ C - gk .* Wnew .* (Vnew-vk) ./C - gl .* (Vnew-vl)./C +Iext ./C;
    dwdt_I=phi .* (w_inf(Vnew,V3,V4)- Wnew)./ tau_w(Vnew,V3,V4);
    dX = [dvdt_I; dwdt_I];
    end
end

function res=m_inf(V,V1,V2)
res=0.5*(1+tanh((V-V1)/V2));
end
function res=w_inf(V,V3,V4)
res=0.5*(1+tanh((V-V3)/V4));
end
function res=tau_w(V,V3,V4)
res=1./(cosh((V-V3)/V4));
end