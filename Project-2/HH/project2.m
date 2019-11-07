%%
global Iext
global Phi
global Gca
global Gk Gna
global Gl
global Vca Vna
global Vk
global Vl
global V1
global V2
global V3
global V4
global C
global fni
global hinit

%%



%%% Question 12
Iext = 0;
Gk = 36;
Gl = 0.3;
Gna = 120;
Vna = 55;
Vk = -72;
C = 1;



%%
%%% Question 13

V = -60;
alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12);
betan = 0.125*exp(-(V+60)/80);
alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12);
betam = 4*exp(-(V+60)/18);
alphah = 0.07*(exp(-(V+60)/20));
betah = 1/(exp(-(V+30)/10)+1);

n = alphan/(alphan + betan);
m = alpham/(alpham + betam);
h = alphah/(alphah + betah);

Vl = V - (Iext -Gna*m^3*h*(V-Vna) - Gk*n^4*(V-Vk))/Gl;
display(Vl)
%%
%%% Question 14
disp('Question 14')
figure
Iext = 0;
y_init = [-60; 0;0;0];
equiPoint3 = fsolve(@hh,y_init,optionsfsolve);
plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)

nc1 = @(V)((Iext - Gk*equiPoint3(2)^4*(V-Vk) - Gl*(V-Vl))./(Gna*equiPoint3(4)*(V-Vna))).^1/3;
nc2 = @(V) (-0.1*(V+35)./(exp(-(V+35)/10)-1 + 10^-12))./((-0.1*(V+35)./(exp(-(V+35)/10)-1 + 10^-12))+ (4*exp(-(V+60)/18)));

V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));
equiPoint1 = equiPoint3;

for i = 0:1:20
    equiPoint1(1) = equiPoint3(1) + i ;
    deltaT = [0,1000];
    [t, y] = ode45(@hh4ode,deltaT,equiPoint1);
    plot(y(:,1),100*y(:,3))
end

figure
for i = 0:1:20
    equiPoint1(1) = equiPoint3(1) + i ;
    deltaT = [0,1000];
    [t, y] = ode45(@hh4ode,deltaT,equiPoint1);
    plot(t,y(:,1)); hold on;
end
xlim([0 10])
syms V m
J = jacobian([ (Iext -(Gna*m^3*(equiPoint3(4))*(V-Vna)) - Gk*(equiPoint3(2))^4*(V-Vk) - Gl*(V-Vl))/C, (-0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12))*(1-m) - (4*exp(-(V+60)/18))*m],[V,m]);
equi = zeros(1,2);
equi(1) = equiPoint3(1);
equi(2) = equiPoint3(3);
A = double(subs(J,[V,m],equi));
[V,D] = eig(A)

%% Question 15
disp('Question 15')
figure
for i = 8:1:12
    Iext = i;
    y_init = [-60; 0;0;0];
    equiPoint3 = fsolve(@hh,y_init,optionsfsolve);
    disp(equiPoint3)

    syms V m
    J = jacobian([ (Iext -(Gna*m^3*((0.07*(exp(-(V+60)/20)))/((0.07*(exp(-(V+60)/20)))+ (1/(exp(-(V+30)/10)+1))))*(V-Vna)) - Gk*((-0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12))/((-0.01*(V+50)/(exp(-(V+50)/10)-1+10^-12))+ (0.125*exp(-(V+60)/80))))^4*(V-Vk) - Gl*(V-Vl))/C, (-0.1*(V+35)/(exp(-(V+35)/10)-1 + 10^-12))*(1-m) - (4*exp(-(V+60)/18))*m],[V,m]);
    equi = zeros(1,2);
    equi(1) = equiPoint3(1);
    equi(2) = equiPoint3(3);
    A = double(subs(J,[V,m],equi));
    [V,D] = eig(A)
    if D(1,1) < 0 && D(2,2) < 0
        display(i)
        display('Stable')
    else
        display(i)
        display('Unstable')
    end
    deltaT = [0,1000];
    [t, y] = ode15s(@hh4ode,deltaT,equiPoint1);
    plot(t,y(:,1)); hold on;
end    

xlim([0 100])
%%
%%% Question 16
figure(11)
fni = 0
Iext = 0;
y_init = [-60; 0;0;0];
equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
%plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)
equiPoint1 = equiPoint3;
equiPoint1(1) = equiPoint3(1) + 10 ;
deltaT = [0,1000];
[t, y] = ode45(@hhPara4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,3)); hold on;

fni = 0.1

equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
%plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)
equiPoint1 = equiPoint3;
equiPoint1(1) = equiPoint3(1) + 10 ;
deltaT = [0,1000];
[t, y] = ode45(@hhPara4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,3))

fni = 0.17

equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
%plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)
equiPoint1 = equiPoint3;
equiPoint1(1) = equiPoint3(1) + 10 ;
deltaT = [0,1000];
[t, y] = ode45(@hhPara4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,3))

fni = 0.2
equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
%plot(equiPoint3(1),100*equiPoint3(3),'b*')
disp(equiPoint3)
equiPoint1 = equiPoint3;
equiPoint1(1) = equiPoint3(1) + 10 ;
deltaT = [0,1000];
[t, y] = ode45(@hhPara4ode,deltaT,equiPoint1);
plot(y(:,1),100*y(:,3))


%% Question 17 

figure

Iext = 0;

y_init = [-60; 0;0;0];
equiPoint3 = fsolve(@hh,y_init,optionsfsolve);
plot(equiPoint3(1),100*equiPoint3(2),'b*')
disp(equiPoint3)
hinit = equiPoint3(4);
equiPoint1 = zeros(2,1);
equiPoint1(1) =  equiPoint3(1);
equiPoint1(2) =  equiPoint3(2);
for i = 0:1:20
    equiPoint1(1) = equiPoint3(1) + i ;
    deltaT = [0,1000];
    [t, y] = ode45(@hhReduced4ode,deltaT,equiPoint1);
    plot(t, y(:,1));hold on;
end
ylim([-100 100])
xlim([0 10])

%% Question 18
figure

am = @(V) -0.1 * (35+V) ./ (exp(-0.1*(35+V)) - 1);
bm = @(V) 4 * exp(-(60+V)/18);

an = @(V) 0.01 * (-(50+V)) ./ (exp(-0.1*(50+V)) - 1);
bn = @(V) 0.125 * exp(-(60+V)/80);


nc1 = @(V)((Iext - Gk*equiPoint3(4)*((am(V))./((am(V))+(bm(V)))).^3.*(V-Vna) - Gl*(V-Vl))./(Gk*(V-Vk))).^1/4;
nc2 = @(V) an(V)./(an(V)+bn(V));

V = -100:1:100;
plot(V,100*nc1(V)); hold on;

ylim([0 100])
xlim([-100 100])

plot(V,100*nc2(V));



for i = 0.2:0.1:0.4
    fni = i;
    equiPoint3 = fsolve(@hhPara,y_init,optionsfsolve);
    %plot(equiPoint3(1),100*equiPoint3(3),'b*')
    disp(equiPoint3)
    hinit = equiPoint3(4);
    equiPoint1 = zeros(2,1);
    equiPoint1(1) =  equiPoint3(1);
    equiPoint1(2) =  equiPoint3(2);
    equiPoint1(1) = equiPoint3(1) + 10 ;
    deltaT = [0,1000];
    [t, y] = ode45(@hhReducedPara4ode,deltaT,equiPoint1);
    plot(y(:,1),100*y(:,2)); hold on;
end


%% Question 19
Iext = 0;
Gk = 36;
Gl = 0.3;
Gna = 120;
Vna = 55;
Vk = -72;
C = 1;

am = @(V) -0.1 * (35+V) ./ (exp(-0.1*(35+V)) - 1);
bm = @(V) 4 * exp(-(60+V)/18);

ah = @(V) 0.07 * exp(-(60+V)/20);
bh = @(V) 1 ./ (exp(-(30+V)/10) + 1);
 

an = @(V) 0.01 * (-(50+V)) ./ (exp(-0.1*(50+V)) - 1);
bn = @(V) 0.125 * exp(-(60+V)/80);


Iext = 0;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               am(y(1)).*(1-y(2))-bm(y(1)).*y(2); ...
               an(y(1)).*(1-y(3))-bn(y(1)).*y(3); ...
               ah(y(1)).*(1-y(4))-bh(y(1)).*y(4) ];
[T1,Y1]=ode15s(f, [0 100], [-60 0.052 0.317 0.596]);
Iext = -3;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               am(y(1)).*(1-y(2))-bm(y(1)).*y(2); ...
               an(y(1)).*(1-y(3))-bn(y(1)).*y(3); ...
               ah(y(1)).*(1-y(4))-bh(y(1)).*y(4) ];
[T2,Y2]=ode15s(f, [T1(length(T1)) 20 + T1(length(T1))], [Y1(length(Y1),1),Y1(length(Y1),2),Y1(length(Y1),3),Y1(length(Y1),4)]);
Iext = 0;
f = @(t,y) [ (-Gk*y(3).^4.*(y(1) - Vk) - Gna*y(2).^3.*y(4).*(y(1)-Vna)-Gl.*(y(1)-Vl) + Iext)/C; ...
               am(y(1)).*(1-y(2))-bm(y(1)).*y(2); ...
               an(y(1)).*(1-y(3))-bn(y(1)).*y(3); ...
               ah(y(1)).*(1-y(4))-bh(y(1)).*y(4) ];
[T3,Y3]=ode15s(f, [T2(length(T2)) 150 + T2(length(T2))], [Y2(length(Y2),1),Y2(length(Y2),2),Y2(length(Y2),3),Y2(length(Y2),4)]);
figure
plot([T1;T2;T3], [Y1(:,1); Y2(:,1); Y3(:,1)])
xlabel('t'); ylabel('V(t)');
title('Membrane potential vs time');

%% Question 20

n = equiPoint3(2);
h = equiPoint3(4);

myfun1 = @(V,m) (-Gk*n.^4.*(V - Vk) - Gna*m.^3.*h.*(V-Vna)-Gl.*(V-Vl) + Iext)/C;
myfun2 = @(V,m) am(V).*(1-m)-bm(V).*m;
figure
a1=ezplot(@(V,m) myfun1(V,m), [-72 120 0 1])
set(a1,'Color','red', 'LineStyle', '-', 'LineWidth', 1);
hold on;
set(gca, 'fontsize', 10)
a2=ezplot(@(V,m) myfun2(V,m), [-72 120 0 1])
title('Phase plane for n and h fixed at rest values');

n = Y2(length(Y2),3);
h = Y2(length(Y2),4);

myfun1 = @(V,m) (-Gk*n.^4.*(V - Vk) - Gna*m.^3.*h.*(V-Vna)-Gl.*(V-Vl) + Iext)/C;
myfun2 = @(V,m) am(V).*(1-m)-bm(V).*m;
figure
a1=ezplot(@(V,m) myfun1(V,m), [-72 120 0 1])
set(a1,'Color','red', 'LineStyle', '-', 'LineWidth', 1);
hold on;
set(gca, 'fontsize', 10)
a2=ezplot(@(V,m) myfun2(V,m), [-72 120 0 1])
title('Phase plane for n and h fixed at end of stimulus values');

