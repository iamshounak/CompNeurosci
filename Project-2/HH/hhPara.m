function dYdt = hhPara(y)
    global Iext
    global Gk
    global Gl Gna
    global Vna
    global Vk
    global Vl
    global C
    global fni
    V = y(1);
    n = y(2);
    m = y(3);
    h = y(4);
    alphan = -0.01*(V+50)/(exp(-(V+50)/10)-1);
    if isnan(alphan)
        alphan = 1;
    end
    betan = 0.125*exp(-(V+60)/80);
    alpham = -0.1*(V+35)/(exp(-(V+35)/10)-1);
    if isnan(alpham)
        alpham = 1;
    end
    betam = 4*exp(-(V+60)/18);
    alphah = 0.07*(exp(-(V+60)/20));
    betah = 1/(exp(-(V+30)/10)+1);
    
    dVdt = (Iext -Gna*m^3*(V-Vna)*fni -Gna*m^3*h*(V-Vna)*(1-fni) - Gk*n^4*(V-Vk) - Gl*(V-Vl))/C;
    dmdt = alpham*(1-m) - betam*m;
    dndt = alphan*(1-n) - betan*n;
    dhdt = alphah*(1-h) - betah*h;
    dYdt = [dVdt;dndt;dmdt;dhdt]; 
end



