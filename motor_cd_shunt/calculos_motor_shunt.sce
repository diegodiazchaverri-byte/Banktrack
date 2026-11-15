// Motor CD shunt - control de velocidad por campo
// Diego Diaz - 2025105226
// Maquinas electricas para mecatronica

clear; clc;
try
    carpeta = get_absolute_file_path("calculos_motor_shunt.sce");
catch
    carpeta = pwd() + "/";
end
chdir(carpeta);
if ~isdir("figuras") then mkdir("figuras"); end

// datos
Va = 120; Ra = 0.5; Rf = 120; IL = 20;
n1 = 1600; n2 = 1800;
La = 0.5; J = 0.5; B = 0.075;
TLB = 5;   // carga caso B

rad2rpm = 60/(2*%pi);
w1 = n1*2*%pi/60;
w2 = n2*2*%pi/60;

mprintf("\n--- PARTE A ---\n");

// A.1 punto inicial 1600 rpm
If1 = Va/Rf;
Ia1 = IL - If1;
Ea1 = Va - Ia1*Ra;
Ke1 = Ea1/w1;
Kt1 = Ke1;   // Kt = Ke en SI
T1 = Kt1*Ia1;
kf = Ke1/If1;

mprintf("\nA.1) n1 = %d rpm, w1 = %.4f rad/s\n", n1, w1);
mprintf("   If = %.4f A\n   Ia = %.4f A\n   Ea = %.4f V\n", If1, Ia1, Ea1);
mprintf("   Ke = %.6f V.s/rad\n   Kt = %.6f N.m/A\n", Ke1, Kt1);
mprintf("   T = %.4f N.m   (Ea*Ia = %.2f W, T*w1 = %.2f W)\n", T1, Ea1*Ia1, T1*w1);

// A.2 Rf para 1800 rpm, par constante
// Ke2*Ia2 = T1  y  Va = Ke2*w2 + Ia2*Ra  ->  w2*Ke2^2 - Va*Ke2 + T1*Ra = 0
a = w2; b = -Va; c = T1*Ra;
disc = b^2 - 4*a*c;
r1 = (-b + sqrt(disc))/(2*a);
r2 = (-b - sqrt(disc))/(2*a);
Ke2 = r1;   // la otra raiz da Ia muy grande
Kt2 = Ke2;
Ia2 = T1/Ke2;
Ea2 = Va - Ia2*Ra;
If2 = Ke2/kf;
Rf2 = Va/If2;
IL2 = Ia2 + If2;

mprintf("\nA.2) n2 = %d rpm, w2 = %.4f rad/s\n", n2, w2);
mprintf("   %.4f*Ke2^2 - %d*Ke2 + %.4f = 0\n", a, Va, c);
mprintf("   raices: %.6f (Ia = %.3f A)  y  %.6f (Ia = %.3f A, se descarta)\n", r1, T1/r1, r2, T1/r2);
mprintf("   Ke2 = %.6f\n   Ia2 = %.4f A\n   Ea2 = %.4f V\n", Ke2, Ia2, Ea2);
mprintf("   If2 = %.5f A\n   Rf'' = %.4f Ohm\n   IL2 = %.4f A\n", If2, Rf2, IL2);

// A.3 estado estable del modelo (derivadas = 0)
function [w, ia] = estado_estable(Ke, Kt, TL, Va, Ra, B)
    w  = (Kt*Va - Ra*TL)/(Kt*Ke + Ra*B);
    ia = (B*w + TL)/Kt;
endfunction

nombres = ["Rf=120 (1600 rpm), Caso A TL=0"; "Rf=120 (1600 rpm), Caso B TL=5"; ..
   "Rf''=" + msprintf("%.2f", Rf2) + " (1800 rpm), Caso A TL=0"; ..
   "Rf''=" + msprintf("%.2f", Rf2) + " (1800 rpm), Caso B TL=5"];
// [Rf Ke TL]
casos = [Rf Ke1 0; Rf Ke1 TLB; Rf2 Ke2 0; Rf2 Ke2 TLB];
nc = size(casos, 1);
wss = zeros(nc,1); iss = zeros(nc,1);

mprintf("\nA.3) Estado estable\n");
mprintf("   %-40s %10s %10s %10s %10s\n", "Caso", "w [rad/s]", "n [rpm]", "ia [A]", "Tm [N.m]");
for k = 1:nc
    [wss(k), iss(k)] = estado_estable(casos(k,2), casos(k,2), casos(k,3), Va, Ra, B);
    mprintf("   %-40s %10.4f %10.2f %10.4f %10.4f\n", nombres(k), wss(k), wss(k)*rad2rpm, iss(k), casos(k,2)*iss(k));
end
mprintf("   nota: con TL=0 da %.1f rpm y no 1800 porque B*w crece con la velocidad\n", wss(3)*rad2rpm);


mprintf("\n--- PARTE B ---\n");
tf = 20; dt = 0.01;
t = (0:dt:tf)';

// B.1 ecuaciones con ode (para comparar con xcos)
function dx = motor_cd(t, x, Ke, Kt, TL)
    dx = [(Va - Ra*x(1) - Ke*x(2))/La;
          (Kt*x(1) - B*x(2) - TL)/J];
endfunction

IA = zeros(length(t), nc); TM = IA; NR = IA;
for k = 1:nc
    K = casos(k,2); TL = casos(k,3);
    x = ode([0; 0], 0, t', 1e-9, 1e-11, list(motor_cd, K, K, TL));
    IA(:,k) = x(1,:)';
    TM(:,k) = K*x(1,:)';
    NR(:,k) = x(2,:)'*rad2rpm;
end

mprintf("\nB.1) ode()\n");
mprintf("   %-40s %8s %7s %9s %8s %8s %9s %9s\n", "Caso", "ia_max", "t_pico", "n_max", "Mp_n", "ts(2%)", "ia(tf)", "n(tf)");
for k = 1:nc
    [imax, kmax] = max(IA(:,k));
    nss = wss(k)*rad2rpm;
    fuera = find(abs(NR(:,k) - nss) > 0.02*nss);
    ts = t(fuera($) + 1);
    mprintf("   %-40s %8.2f %7.2f %9.1f %8.1f %8.2f %9.4f %9.2f\n", nombres(k), imax, t(kmax), ..
        max(NR(:,k)), 100*(max(NR(:,k)) - nss)/nss, ts, IA($,k), NR($,k));
end

// B.2 simulacion del diagrama xcos
// si no existe o no abre en esta version se crea con generar_xcos.sce
archivo_xcos = "Diaz_Diego.zcos";
hay_xcos = %f;
hay_gen = isfile("generar_xcos.sce");
if ~isfile(archivo_xcos) & hay_gen then
    exec("generar_xcos.sce", -1);
end
intento = 0;
while isfile(archivo_xcos) & ~hay_xcos & intento < 2
    intento = intento + 1;
    try
        mprintf("\nB.2) Xcos: %s\n", archivo_xcos);
        loadXcosLibs(); loadScicos();
        importXcosDiagram(archivo_xcos);
        ctx0 = scs_m.props.context;
        XIA = zeros(length(t), nc); XTM = XIA; XNR = XIA;
        for k = 1:nc
            ctx = ctx0;
            i = grep(ctx, "/^Rf\s*=/", "r"); ctx(i) = msprintf("Rf = %.10g;", casos(k,1));
            i = grep(ctx, "/^TL\s*=/", "r"); ctx(i) = msprintf("TL = %.10g;", casos(k,3));
            scs_m.props.context = ctx;
            scs_m.props.tf = tf;
            ventanas = winsid();
            xcos_simulate(scs_m, 4);    // guarda "res"
            nuevas = setdiff(winsid(), ventanas);
            if nuevas <> [] then close(nuevas); end
            XIA(:,k) = interp1(res.time, res.values(:,1), t, "linear", "extrap");
            XTM(:,k) = interp1(res.time, res.values(:,2), t, "linear", "extrap");
            XNR(:,k) = interp1(res.time, res.values(:,3), t, "linear", "extrap");
        end
        hay_xcos = %t;
    catch
        mprintf("   no se pudo simular, se genera de nuevo el diagrama\n");
        if intento == 1 & hay_gen then exec("generar_xcos.sce", -1); end
    end
end

if hay_xcos then
    mprintf("   diferencia xcos vs ode:\n");
    mprintf("   %-40s %14s %14s %14s\n", "Caso", "max|dia| [A]", "max|dTm| [Nm]", "max|dn| [rpm]");
    for k = 1:nc
        mprintf("   %-40s %14.2e %14.2e %14.2e\n", nombres(k), max(abs(XIA(:,k) - IA(:,k))), ..
            max(abs(XTM(:,k) - TM(:,k))), max(abs(XNR(:,k) - NR(:,k))));
    end
    GIA = XIA; GTM = XTM; GNR = XNR; fuente = "Xcos";
else
    mprintf("\nB.2) no se encontro %s, se grafica con ode()\n", archivo_xcos);
    GIA = IA; GTM = TM; GNR = NR; fuente = "ode()";
end

// B.3 graficas ia, Tm, n
function grafica_caso(nfig, cols, titulo, archivo, t, GIA, GTM, GNR, iss, Tss, nss, etiquetas, fuente)
    f = scf(nfig); clf(f);
    f.figure_size = [1000, 900];
    f.background = color("white");
    colores = [color("blue"), color("red")];
    datos = list(GIA, GTM, GNR);
    finales = list(iss, Tss, nss);
    ylab = ["i_a(t)  [A]", "T_m(t)  [N.m]", "n(t)  [rpm]"];
    tit = ["Corriente de armadura", "Torque electromagnetico", "Velocidad angular"];
    for s = 1:3
        subplot(3, 1, s);
        D = datos(s); F = finales(s);
        for j = 1:2
            plot(t, D(:, cols(j)));
            e = gce(); e.children.foreground = colores(j); e.children.thickness = 2;
        end
        for j = 1:2   // valor teorico
            plot([t(1) t($)], [F(cols(j)) F(cols(j))]);
            e = gce(); e.children.foreground = colores(j); e.children.line_style = 2;
        end
        a = gca(); a.grid = [color("gray") color("gray")];
        a.data_bounds(1,1) = 0; a.data_bounds(2,1) = t($);
        xlabel("t  [s]"); ylabel(ylab(s));
        title(titulo + " - " + tit(s) + "  (" + fuente + ")");
        legend([etiquetas(1); etiquetas(2); "valor teorico estable"; "valor teorico estable"], "in_upper_right");
    end
    xs2png(f, archivo);
endfunction

Tss = casos(:,2).*iss;
nssr = wss*rad2rpm;
et = ["Rf = 120 Ohm (1600 rpm)", "Rf'' = " + msprintf("%.2f", Rf2) + " Ohm (1800 rpm)"];
grafica_caso(1, [1 3], "Caso A (TL = 0)", "figuras/caso_A_TL0.png", t, GIA, GTM, GNR, iss, Tss, nssr, et, fuente);
grafica_caso(2, [2 4], "Caso B (TL = 5 N.m)", "figuras/caso_B_TL5.png", t, GIA, GTM, GNR, iss, Tss, nssr, et, fuente);

mprintf("\nfiguras en %sfiguras/\n", carpeta);
