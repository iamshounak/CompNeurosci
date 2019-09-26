%%
function dYdt = mle4ode(t,y)
    global Iext
    global Phi
    global Gca
    global Gk
    global Gl
    global Vca
    global Vk
    global Vl
    global V1
    global V2
    global V3
    global V4
    global C
    V = y(1);
    w = y(2);
    dYdt1 = (Iext -Gca*0.5*(1+tanh((V-V1)/V2)).*(V-Vca) - Gk*w.*(V-Vk) - Gl*(V-Vl))/C;
    dYdt2 = (0.5*(1+tanh((V-V3)/V4)) - w)*Phi.*cosh((V-V3)/(V4*2));
    dYdt = [dYdt1;dYdt2]; 
end
