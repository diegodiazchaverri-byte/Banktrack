// =====================================================================
//  calculos_motor_shunt.sce
//  Control de velocidad de un motor CD Shunt por ajuste del campo (Rf)
//
//  PARTE A : cálculos analíticos (memoria de cálculo impresa en consola)
//  PARTE B : simulación dinámica
//            B.1  integración numérica de las EDO con ode()  (referencia)
//            B.2  simulación del diagrama Xcos Diaz_Diego.zcos (se genera
//                 automáticamente con generar_xcos.sce si no existe)
//                 (si existe en la misma carpeta) y comparación con B.1
//            B.3  gráficas ia(t), Tm(t), n(t)  ->  carpeta figuras/
//            B.4  (extra) transición 1600 -> 1800 rpm en caliente
//
//  Ejecutar desde la carpeta del archivo:  exec("calculos_motor_shunt.sce", -1)
//  Probado en Scilab 2024.0.0
// =====================================================================
clear; clc;
// Carpeta del script: get_absolute_file_path funciona en exec(); si no, pwd()
try
    carpeta = get_absolute_file_path("calculos_motor_shunt.sce");
catch
    carpeta = pwd() + "/";
end
chdir(carpeta);
if ~isdir("figuras") then mkdir("figuras"); end

// ---------------------------------------------------------------------
// DATOS DEL PROBLEMA
// ---------------------------------------------------------------------
Va  = 120;     // tensión terminal [V]
Ra  = 0.5;     // resistencia de armadura [Ohm]
Rf  = 120;     // resistencia de campo inicial [Ohm]
IL  = 20;      // corriente de línea a 1600 rpm [A]
n1  = 1600;    // velocidad inicial [rpm]
n2  = 1800;    // velocidad deseada [rpm]
La  = 0.5;     // inductancia de armadura [H]
J   = 0.5;     // momento de inercia [kg.m^2]
B   = 0.075;   // fricción viscosa [N.m.s/rad]
TLB = 5;       // par de carga adicional del Caso B [N.m]

rpm2rad = 2*%pi/60;   // 1 rpm = 2*pi/60 rad/s
rad2rpm = 60/(2*%pi);
w1 = n1*rpm2rad;
w2 = n2*rpm2rad;

linea = "------------------------------------------------------------------";
mprintf("\n%s\n PARTE A - CALCULOS ANALITICOS\n%s\n", linea, linea);

// ---------------------------------------------------------------------
// A.1  Punto de operación inicial (1600 rpm)
// ---------------------------------------------------------------------
// Conexión shunt: el campo está en paralelo con la armadura, ambos a Va.
If1 = Va/Rf;                  // corriente de campo
Ia1 = IL - If1;               // LCK en el nodo: IL = Ia + If
Ea1 = Va - Ia1*Ra;            // LVK en la armadura: Va = Ea + Ia*Ra
Ke1 = Ea1/w1;                 // Ea = Ke*w   ->  Ke = K*phi  [V.s/rad]
Kt1 = Ke1;                    // en SI, Kt [N.m/A] = Ke [V.s/rad] (balance de potencia)
T1  = Kt1*Ia1;                // par electromagnético desarrollado
Pc1 = Ea1*Ia1;                // potencia convertida (eléctrica -> mecánica)
kf  = Ke1/If1;                // Ke por amperio de campo (zona lineal: phi ~ If)

mprintf("\nA.1) Punto de operacion inicial (n1 = %d rpm, w1 = %.4f rad/s)\n", n1, w1);
mprintf("   If  = Va/Rf          = %d/%d           = %.4f A\n", Va, Rf, If1);
mprintf("   Ia  = IL - If        = %d - %.1f        = %.4f A\n", IL, If1, Ia1);
mprintf("   Ea  = Va - Ia*Ra     = %d - %.1f*%.1f   = %.4f V\n", Va, Ia1, Ra, Ea1);
mprintf("   Ke  = Ea/w1          = %.1f/%.4f   = %.6f V.s/rad\n", Ea1, w1, Ke1);
mprintf("   Kt  = Ke (SI)                          = %.6f N.m/A\n", Kt1);
mprintf("   Tm  = Kt*Ia                            = %.4f N.m\n", T1);
mprintf("   Pconv = Ea*Ia = %.2f W  ;  Tm*w1 = %.2f W  (deben coincidir)\n", Pc1, T1*w1);
mprintf("   Ke por amperio de campo: kf = Ke/If    = %.6f V.s/(rad.A)\n", kf);
mprintf("   Observacion: B*w1 = %.4f N.m ~ Tm = %.4f N.m  -> con TL = 0 la\n", B*w1, T1);
mprintf("   friccion viscosa B ES la carga del punto inicial (ver guia).\n");

// ---------------------------------------------------------------------
// A.2  Nueva resistencia de campo para 1800 rpm a par constante
// ---------------------------------------------------------------------
// Incógnitas: Ke2 (= K*phi2) e Ia2. Dos ecuaciones:
//   (1) par constante : Ke2*Ia2 = Ke1*Ia1 = T1
//   (2) armadura      : Va = Ke2*w2 + Ia2*Ra
// Sustituyendo Ia2 = T1/Ke2 en (2) y multiplicando por Ke2:
//       w2*Ke2^2 - Va*Ke2 + T1*Ra = 0     (ecuación cuadrática en Ke2)
a = w2; b = -Va; c = T1*Ra;
disc = b^2 - 4*a*c;
Ke2_alta = (-b + sqrt(disc))/(2*a);
Ke2_baja = (-b - sqrt(disc))/(2*a);
Ke2 = Ke2_alta;               // raíz físicamente aceptable (ver abajo)
Kt2 = Ke2;
Ia2 = T1/Ke2;
Ea2 = Va - Ia2*Ra;
If2 = Ke2/kf;                 // phi proporcional a If (circuito magnético lineal)
Rf2 = Va/If2;                 // = Rf*Ke1/Ke2
IL2 = Ia2 + If2;

mprintf("\nA.2) Nueva Rf para n2 = %d rpm (w2 = %.4f rad/s) a par constante T = %.4f N.m\n", n2, w2, T1);
mprintf("   Ecuacion: w2*Ke2^2 - Va*Ke2 + T*Ra = 0\n");
mprintf("             %.4f*Ke2^2 - %d*Ke2 + %.4f = 0\n", a, Va, c);
mprintf("   Discriminante = %.4f\n", disc);
mprintf("   Raices: Ke2 = %.6f  (Ia2 = %.3f A)   |   Ke2 = %.6f  (Ia2 = %.3f A)\n", ..
        Ke2_alta, T1/Ke2_alta, Ke2_baja, T1/Ke2_baja);
mprintf("   -> Se descarta la raiz pequena: exige Ia = %.0f A (~%.0f veces la nominal),\n", T1/Ke2_baja, (T1/Ke2_baja)/Ia1);
mprintf("      Ea casi nula y un flujo de apenas %.1f%% del inicial: punto inestable/irreal.\n", 100*Ke2_baja/Ke1);
mprintf("   Ke2 = Kt2 = %.6f   (phi2/phi1 = %.5f)\n", Ke2, Ke2/Ke1);
mprintf("   Ia2 = T/Ke2 = %.4f A ;  Ea2 = Va - Ia2*Ra = %.4f V  (check Ke2*w2 = %.4f V)\n", Ia2, Ea2, Ke2*w2);
mprintf("   If2 = If1*(Ke2/Ke1) = %.5f A\n", If2);
mprintf("   Rf'' = Va/If2 = %.4f Ohm   (aumento de %.2f Ohm)\n", Rf2, Rf2 - Rf);
mprintf("   IL2 = Ia2 + If2 = %.4f A ;  Pconv2 = Ea2*Ia2 = %.2f W = T*w2 = %.2f W\n", IL2, Ea2*Ia2, T1*w2);
mprintf("   Aproximacion rapida (Ea ~ cte):  Rf'' ~ Rf*n2/n1 = %.2f Ohm  (error %.2f%%)\n", Rf*n2/n1, 100*abs(Rf*n2/n1 - Rf2)/Rf2);

// ---------------------------------------------------------------------
// A.3  Valores teóricos en estado estable del MODELO DINÁMICO
// ---------------------------------------------------------------------
// En estado estable dia/dt = 0 y dw/dt = 0:
//   0 = Va - Ra*ia - Ke*w      ->  ia = (Va - Ke*w)/Ra
//   0 = Kt*ia - B*w - TL
//   => w_ss  = (Kt*Va - Ra*TL)/(Kt*Ke + Ra*B)
//      ia_ss = (B*w_ss + TL)/Kt
function [w, ia] = estado_estable(Ke, Kt, TL, Va, Ra, B)
    w  = (Kt*Va - Ra*TL)/(Kt*Ke + Ra*B);
    ia = (B*w + TL)/Kt;
endfunction

// Casos: [Rf  Ke  TL]
nombres = ["Rf=120 (1600 rpm), Caso A TL=0"; "Rf=120 (1600 rpm), Caso B TL=5"; ..
           "Rf''=" + msprintf("%.2f", Rf2) + " (1800 rpm), Caso A TL=0"; ..
           "Rf''=" + msprintf("%.2f", Rf2) + " (1800 rpm), Caso B TL=5"];
casos = [Rf  Ke1 0;  Rf  Ke1 TLB;  Rf2 Ke2 0;  Rf2 Ke2 TLB];
nc = size(casos, 1);
wss = zeros(nc,1); iss = zeros(nc,1);

mprintf("\nA.3) Estado estable del modelo dinamico (incluye friccion B)\n");
mprintf("   w_ss = (Kt*Va - Ra*TL)/(Kt*Ke + Ra*B) ;  ia_ss = (B*w_ss + TL)/Kt\n");
mprintf("   %-40s %10s %10s %10s %10s\n", "Caso", "w_ss", "n_ss", "ia_ss", "Tm_ss");
mprintf("   %-40s %10s %10s %10s %10s\n", "", "[rad/s]", "[rpm]", "[A]", "[N.m]");
for k = 1:nc
    [wss(k), iss(k)] = estado_estable(casos(k,2), casos(k,2), casos(k,3), Va, Ra, B);
    mprintf("   %-40s %10.4f %10.2f %10.4f %10.4f\n", nombres(k), wss(k), wss(k)*rad2rpm, iss(k), casos(k,2)*iss(k));
end

// Dinámica: La*J*s^2 + (Ra*J + La*B)*s + (Ra*B + Kt*Ke) = 0
mprintf("\n   Dinamica (polinomio caracteristico La*J*s^2 + (Ra*J+La*B)*s + (Ra*B+Kt*Ke)):\n");
mprintf("   tau_e = La/Ra = %.3f s ;  tau_m = J/B = %.3f s\n", La/Ra, J/B);
for k = [1 3]
    K  = casos(k,2);
    p2 = La*J; p1 = Ra*J + La*B; p0 = Ra*B + K*K;
    wn = sqrt(p0/p2); zeta = p1/(2*sqrt(p0*p2));
    polos = roots(poly([p0 p1 p2], "s", "coeff"));
    Mp = exp(-zeta*%pi/sqrt(1 - zeta^2));
    mprintf("   Ke=%.4f: polos = %.4f +/- j%.4f ; wn = %.4f rad/s ; zeta = %.4f\n", K, real(polos(1)), abs(imag(polos(1))), wn, zeta);
    mprintf("             sobreimpulso teorico de w (TL=0) = %.1f%% ;  ts(2%%) ~ 4/(zeta*wn) = %.2f s\n", 100*Mp, 4/(zeta*wn));
end

// Nota sobre la diferencia "par constante" vs "carga viscosa"
// Si la única carga es B*w (TL = 0), ¿qué Rf da exactamente 1800 rpm en el modelo?
Ke3 = (Va + sqrt(Va^2 - 4*w2*Ra*B*w2))/(2*w2);
Rf3 = Rf*Ke1/Ke3;
mprintf("\n   NOTA: con Rf'' = %.2f Ohm el modelo con TL=0 llega a %.1f rpm (no 1800)\n", Rf2, wss(3)*rad2rpm);
mprintf("   porque la carga viscosa B*w CRECE con la velocidad (%.2f -> %.2f N.m) y\n", B*w1, B*w2);
mprintf("   no es de par constante. Para 1800 rpm exactos en el modelo viscoso: Rf = %.2f Ohm.\n", Rf3);

// =====================================================================
//  PARTE B - SIMULACIÓN DINÁMICA
// =====================================================================
mprintf("\n%s\n PARTE B - SIMULACION DINAMICA\n%s\n", linea, linea);
tf = 20; dt = 0.01;
t  = (0:dt:tf)';

// ---------------- B.1  Integración de las EDO con ode() ----------------
// Estado x = [ia; w]
//   dia/dt = (Va - Ra*ia - Ke*w)/La
//   dw/dt  = (Kt*ia - B*w - TL)/J
function dx = motor_cd(t, x, Ke, Kt, TL)
    ia = x(1); w = x(2);
    dx = [(Va - Ra*ia - Ke*w)/La;
          (Kt*ia - B*w - TL)/J];
endfunction

IA = zeros(length(t), nc); TM = IA; NR = IA;
for k = 1:nc
    K = casos(k,2); TL = casos(k,3);
    x = ode([0; 0], 0, t', 1e-9, 1e-11, list(motor_cd, K, K, TL));  // rtol, atol
    IA(:,k) = x(1,:)';
    TM(:,k) = K*x(1,:)';
    NR(:,k) = x(2,:)'*rad2rpm;
end

// Indicadores de la respuesta
mprintf("\nB.1) Resultados de ode() (t = 0 .. %d s, condiciones iniciales nulas)\n", tf);
mprintf("   %-40s %8s %7s %9s %8s %8s %9s %9s\n", "Caso", "ia_max", "t_pico", "n_max", "Mp_n", "ts(2%)", "ia(tf)", "n(tf)");
mprintf("   %-40s %8s %7s %9s %8s %8s %9s %9s\n", "", "[A]", "[s]", "[rpm]", "[%]", "[s]", "[A]", "[rpm]");
for k = 1:nc
    [imax, kmax] = max(IA(:,k));
    nss = wss(k)*rad2rpm;
    fuera = find(abs(NR(:,k) - nss) > 0.02*nss);
    ts = t(fuera($) + 1);
    mprintf("   %-40s %8.2f %7.2f %9.1f %8.1f %8.2f %9.4f %9.2f\n", nombres(k), imax, t(kmax), ..
            max(NR(:,k)), 100*(max(NR(:,k)) - nss)/nss, ts, IA($,k), NR($,k));
end

// ---------------- B.2  Simulación del diagrama Xcos ----------------
// El diagrama lo construye generar_xcos.sce CON TU PROPIA VERSIÓN de Scilab
// (los .zcos guardados por otra versión pueden no abrir bien). Si el archivo
// no existe, o si existe pero falla al simularse, se vuelve a generar.
archivo_xcos = "Diaz_Diego.zcos";
hay_xcos = %f;
hay_generador = isfile("generar_xcos.sce");
if ~isfile(archivo_xcos) & hay_generador then
    mprintf("\nB.2) No existe %s: generandolo con tu version de Scilab...\n", archivo_xcos);
    exec("generar_xcos.sce", -1);
end
intento = 0;
while isfile(archivo_xcos) & ~hay_xcos & intento < 2
    intento = intento + 1;
    try
        mprintf("\nB.2) Simulando el diagrama %s con xcos_simulate...\n", archivo_xcos);
        loadXcosLibs(); loadScicos();
        importXcosDiagram(archivo_xcos);
        ctx0 = scs_m.props.context;
        XIA = zeros(length(t), nc); XTM = XIA; XNR = XIA;
        for k = 1:nc
            ctx = ctx0;
            i = grep(ctx, "/^Rf\s*=/", "r");  ctx(i) = msprintf("Rf = %.10g;", casos(k,1));
            i = grep(ctx, "/^TL\s*=/", "r");  ctx(i) = msprintf("TL = %.10g;", casos(k,3));
            scs_m.props.context = ctx;
            scs_m.props.tf = tf;
            ventanas = winsid();
            xcos_simulate(scs_m, 4);               // deja la variable "res" (bloque To workspace)
            nuevas = setdiff(winsid(), ventanas);  // cerrar el osciloscopio del diagrama
            if nuevas <> [] then close(nuevas); end
            // res.time = 0, dt, ..., tf-dt  ->  interpolar sobre t
            XIA(:,k) = interp1(res.time, res.values(:,1), t, "linear", "extrap");
            XTM(:,k) = interp1(res.time, res.values(:,2), t, "linear", "extrap");
            XNR(:,k) = interp1(res.time, res.values(:,3), t, "linear", "extrap");
        end
        hay_xcos = %t;
    catch
        mprintf("   Aviso: no se pudo simular %s:\n   %s\n", archivo_xcos, strcat(lasterror(), " "));
        if intento == 1 & hay_generador then
            mprintf("   -> Se vuelve a generar el diagrama con tu version de Scilab y se reintenta.\n");
            exec("generar_xcos.sce", -1);
        end
    end
end

if hay_xcos then
    mprintf("   %-40s %14s %14s %14s\n", "Caso", "max|dia| [A]", "max|dTm| [Nm]", "max|dn| [rpm]");
    for k = 1:nc
        mprintf("   %-40s %14.2e %14.2e %14.2e\n", nombres(k), max(abs(XIA(:,k) - IA(:,k))), ..
                max(abs(XTM(:,k) - TM(:,k))), max(abs(XNR(:,k) - NR(:,k))));
    end
    mprintf("   (diferencias pequenas = el diagrama Xcos implementa correctamente las EDO)\n");
    // Para las gráficas se usan los datos del diagrama Xcos
    GIA = XIA; GTM = XTM; GNR = XNR; fuente = "Xcos";
else
    mprintf("\nB.2) No se pudo usar %s: se grafican los resultados de ode().\n", archivo_xcos);
    mprintf("     (Pon generar_xcos.sce en esta misma carpeta y vuelve a ejecutar)\n");
    GIA = IA; GTM = TM; GNR = NR; fuente = "ode()";
end

// ---------------- B.3  Gráficas ----------------
function grafica_caso(nfig, cols, titulo, archivo, t, GIA, GTM, GNR, iss, Tss, nss, etiquetas, fuente)
    f = scf(nfig); clf(f);
    f.figure_size = [1000, 900];
    f.background = color("white");
    colores = [color("blue"), color("red")];
    datos = list(GIA, GTM, GNR);
    finales = list(iss, Tss, nss);
    ylab = ["i_a(t)  [A]", "T_m(t)  [N.m]", "n(t)  [rpm]"];
    tit  = ["Corriente de armadura", "Torque electromagnetico", "Velocidad angular"];
    for s = 1:3
        subplot(3, 1, s);
        D = datos(s); F = finales(s);
        for j = 1:2
            plot(t, D(:, cols(j)));
            e = gce(); e.children.foreground = colores(j); e.children.thickness = 2;
        end
        for j = 1:2
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

// Figura resumen: velocidad de los 4 casos
f = scf(3); clf(f); f.figure_size = [1000, 500]; f.background = color("white");
cl = [color("blue"), color("cyan"), color("red"), color("magenta")];
for k = 1:nc
    plot(t, GNR(:,k)); e = gce(); e.children.foreground = cl(k); e.children.thickness = 2;
end
a = gca(); a.grid = [color("gray") color("gray")];
xlabel("t  [s]"); ylabel("n(t)  [rpm]");
title("Velocidad - comparacion de los cuatro casos (" + fuente + ")");
legend(nombres, "in_lower_right");
xs2png(f, "figuras/velocidad_4_casos.png");

// ---------------- B.4  (Extra) Transición en caliente 1600 -> 1800 rpm ----------------
// El motor arranca con Rf = 120 y, ya en régimen (t = 15 s), se cambia a Rf''.
// Así se ve el "control por campo" tal como ocurre en la práctica.
t_sw = 15; tf4 = 40; t4 = (0:dt:tf4)';
// Se integra en dos tramos (Ke constante en cada uno) para respetar la
// discontinuidad exacta en t_sw.
ta = t4(t4 < t_sw)'; tb = t4(t4 >= t_sw)';
xa = ode([0;0], 0, ta, 1e-9, 1e-11, list(motor_cd, Ke1, Ke1, 0));
x_sw = ode([0;0], 0, t_sw, 1e-9, 1e-11, list(motor_cd, Ke1, Ke1, 0));   // estado en t_sw
xb = ode(x_sw, t_sw, tb, 1e-9, 1e-11, list(motor_cd, Ke2, Ke2, 0));
x4 = [xa, xb];
Kt4 = Ke1*ones(t4); Kt4(t4 >= t_sw) = Ke2;
ia4 = x4(1,:)'; n4 = x4(2,:)'*rad2rpm; Tm4 = Kt4.*ia4;
[ipk4, kpk4] = max(ia4(t4 >= t_sw));
mprintf("\nB.4) Extra: cambio de Rf = 120 -> %.2f Ohm en t = %d s (TL = 0)\n", Rf2, t_sw);
mprintf("   Antes: n = %.1f rpm, ia = %.2f A.  Despues: n = %.1f rpm, ia = %.2f A\n", ..
        n4(find(t4 < t_sw)($)), ia4(find(t4 < t_sw)($)), n4($), ia4($));
mprintf("   Pico transitorio de ia tras el cambio = %.2f A (en t = %.2f s)\n", ipk4, t_sw + (kpk4 - 1)*dt);

f = scf(4); clf(f); f.figure_size = [1000, 900]; f.background = color("white");
datos4 = list(ia4, Tm4, n4); yl4 = ["i_a(t)  [A]", "T_m(t)  [N.m]", "n(t)  [rpm]"];
for s = 1:3
    subplot(3,1,s); plot(t4, datos4(s)); e = gce(); e.children.foreground = color("blue"); e.children.thickness = 2;
    a = gca(); a.grid = [color("gray") color("gray")];
    xlabel("t  [s]"); ylabel(yl4(s));
    title("Extra: aumento de Rf en t = 15 s (debilitamiento de campo) - " + yl4(s));
end
xs2png(f, "figuras/extra_transicion_1600_1800.png");

mprintf("\nFiguras guardadas en %sfiguras/\n", carpeta);
