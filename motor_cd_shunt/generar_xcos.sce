// =====================================================================
//  generar_xcos.sce
//  Construye por código el diagrama de bloques del motor CD Shunt y lo
//  guarda como Diaz_Diego.zcos.
//
//  Ecuaciones implementadas (dominio del tiempo):
//     La * dia/dt = Va - Ra*ia - Ke*w          (malla de armadura)
//     J  * dw/dt  = Kt*ia - B*w - TL            (2a ley de Newton rotacional)
//     Tm = Kt*ia ,   Ea = Ke*w ,   n[rpm] = w*30/pi
//
//  Ke y Kt NO se escriben a mano: el contexto del diagrama los calcula a
//  partir de Rf (control por flujo): Ke = kf*If = kf*Va/Rf.
//  Para cambiar de caso basta con editar Rf y TL en el CONTEXTO
//  (Simulación > Establecer contexto) y volver a simular.
//
//  Uso:  exec("generar_xcos.sce", -1)   (desde la carpeta del archivo)
// =====================================================================
funcprot(0);   // permite re-ejecutar el script sin avisos de funciones redefinidas
loadXcosLibs(); loadScicos();

archivo = "Diaz_Diego.zcos";

// ---------------- Contexto del diagrama (parámetros) -----------------
ctx = [
"// ===== Datos del motor CD Shunt ====="
"Va = 120;      // Tension terminal [V]"
"Ra = 0.5;      // Resistencia de armadura [Ohm]"
"La = 0.5;      // Inductancia de armadura [H]"
"J  = 0.5;      // Momento de inercia [kg.m^2]"
"B  = 0.075;    // Friccion viscosa [N.m.s/rad]"
"// ===== Punto de operacion inicial (calculo de Ke, Kt) ====="
"Rf0 = 120; IL0 = 20; n0 = 1600;"
"If0 = Va/Rf0;                 // 1 A"
"Ia0 = IL0 - If0;              // 19 A"
"Ea0 = Va - Ia0*Ra;            // 110.5 V"
"Ke0 = Ea0/(n0*2*%pi/60);      // 0.6595 V.s/rad"
"kf  = Ke0/If0;                // Ke por amperio de campo (zona lineal)"
"// ===== CASO A SIMULAR (editar aqui) ====="
"Rf = 120;      // 120 -> 1600 rpm  |  136.63 -> 1800 rpm (Rf prima)"
"TL = 0;        // 0 -> Caso A      |  5 -> Caso B  [N.m]"
"// ===== Constantes que dependen del flujo ====="
"If = Va/Rf;    // corriente de campo [A]"
"Ke = kf*If;    // constante de FCEM [V.s/rad]"
"Kt = Ke;       // constante de par [N.m/A] (SI: Kt = Ke)"
"// ===== Simulacion ====="
"tf = 20; dt = 0.01; Npts = round(tf/dt) + 1;"
];

// ---------------- Utilidades ----------------
// Nota: en Xcos el eje Y crece hacia ABAJO. "mirror=true" invierte el
// bloque horizontalmente (entrada a la derecha, salida a la izquierda),
// que es lo habitual para los lazos de realimentación.
function b = colocar(b, x, y, w, h, ex, sty, etiqueta)
    b.graphics.orig = [x, y];
    b.graphics.sz   = [w, h];
    if ex <> [] then b.graphics.exprs = ex; end
    if sty <> "" then b.graphics.style = sty; end
    if etiqueta <> "" then b.graphics.id = etiqueta; end
endfunction

function b = gain(x, y, ex, izq, etiqueta)
    b = GAINBLK("define");
    if izq then
        b = colocar(b, x, y, 40, 40, ex, "GAINBLK;flip=false;mirror=true", etiqueta);
    else
        b = colocar(b, x, y, 40, 40, ex, "GAINBLK;flip=false;mirror=false", etiqueta);
    end
endfunction

function b = suma3(x, y, sgn, etiqueta)
    b = SUMMATION("define");
    b = colocar(b, x, y, 40, 60, sci2exp(sgn), "", etiqueta);
    b.model.in = [-1;-1;-1]; b.model.in2 = [-2;-2;-2]; b.model.intyp = [1;1;1];
    b.model.ipar = sgn;
    b.graphics.pin = [0;0;0]; b.graphics.in_implicit = ["E";"E";"E"];
endfunction

function b = split(x, y, n)
    b = SPLIT_f("define");
    b = colocar(b, x, y, 7, 7, [], "", "");
    b.model.out = -ones(n,1); b.model.out2 = -2*ones(n,1); b.model.outtyp = ones(n,1);
    b.graphics.pout = zeros(n,1); b.graphics.out_implicit = emptystr(n,1) + "E";
endfunction

function b = integrador(x, y, etiqueta)
    b = colocar(INTEGRAL_m("define"), x, y, 40, 40, ["0";"0";"0";"1";"-1"], "", etiqueta);
endfunction

// ---------------- Bloques ----------------
scs_m = scicos_diagram();
scs_m.props.title   = "Diaz_Diego";
scs_m.props.tf      = 20;
scs_m.props.context = ctx;
// [atol rtol ttol deltat scale solver hmax]; solver 0 = LSodar
scs_m.props.tol     = [1e-6; 1e-8; 1e-10; 100001; 0; 0; 0];

o = list();
// ----- Circuito de armadura: La dia/dt = Va - Ke w - Ra ia -----
o(1)  = colocar(CONST_m("define"), 20, 290, 40, 40, "Va", "", "Va");
o(2)  = suma3(110, 290, [1;-1;-1], "Va - Ea - Ra ia");      // in1=Va, in2=Ke w, in3=Ra ia
o(3)  = gain(190, 300, "1/La", %f, "1/La");
o(4)  = integrador(270, 300, "ia");
o(5)  = split(350, 317, 3);                                 // nodo ia
o(6)  = gain(270, 400, "Ra", %t, "Ra");                     // caída resistiva (realim.)
// ----- Conversión electromecánica -----
o(7)  = gain(390, 300, "Kt", %f, "Tm = Kt ia");
o(8)  = split(470, 317, 2);                                 // nodo Tm
// ----- Parte mecánica: J dw/dt = Tm - TL - B w -----
o(9)  = suma3(520, 290, [-1;1;-1], "Tm - TL - B w");        // in1=TL, in2=Tm, in3=B w
o(10) = colocar(CONST_m("define"), 540, 200, 40, 40, "TL", "CONST_m;flip=false;mirror=true", "TL");
o(11) = gain(600, 300, "1/J", %f, "1/J");
o(12) = integrador(680, 300, "w [rad/s]");
o(13) = split(760, 317, 3);                                 // nodo w
o(14) = gain(600, 400, "B", %t, "B");                       // fricción viscosa (realim.)
o(15) = gain(400, 480, "Ke", %t, "Ea = Ke w");              // FCEM (realim.)
o(16) = gain(790, 300, "30/%pi", %f, "rad/s a rpm");
// ----- Medición: MUX -> (To workspace, DEMUX -> osciloscopio) -----
mx = MUX("define"); mx.model.in = [-1;-2;-3]; mx.model.ipar = 3;
mx.graphics.pin = [0;0;0]; mx.graphics.in_implicit = ["E";"E";"E"];
o(17) = colocar(mx, 920, 120, 10, 60, "3", "", "");
o(18) = split(957, 147, 2);
dm = DEMUX("define"); dm.model.out = [-1;-2;-3]; dm.model.ipar = 3;
dm.graphics.pout = [0;0;0]; dm.graphics.out_implicit = ["E";"E";"E"];
o(19) = colocar(dm, 1000, 120, 10, 60, "3", "", "");
sc = CMSCOPE("define");
sc.model.in = [1;1;1]; sc.model.in2 = [1;1;1]; sc.model.intyp = [1;1;1];
sc.graphics.pin = [0;0;0]; sc.graphics.in_implicit = ["E";"E";"E"];
o(20) = colocar(sc, 1060, 120, 60, 60, ..
        ["1 1 1"; "1 3 5"; "-1"; "[]"; "[]"; "-20 -10 0"; "140 90 2400"; ..
         "tf tf tf"; "Npts"; "0"; "ia [A];Tm [N.m];n [rpm]"], "", "ia, Tm, n");
o(21) = colocar(TOWS_c("define"), 1000, 230, 60, 40, ["Npts"; "res"; "0"], "", "res = [ia Tm n]");
ck = CLOCK_c("define");
ck.model.rpar.objs(2).graphics.exprs = ["dt"; "0"];
ck.model.rpar.objs(2).model.rpar = [0.01; 0];
ck.model.rpar.objs(2).model.firing = 0;
o(22) = colocar(ck, 1070, 10, 40, 40, [], "", "dt");
o(23) = colocar(CLKSPLIT_f("define"), 1087, 77, 7, 7, [], "", "");

for k = 1:size(o)
    scs_m.objs(k) = o(k);
end
nb = size(o);

// ---------------- Enlaces ----------------
// (bloque origen, puerto) -> (bloque destino, puerto), con puntos de paso
function [scs_m, nb] = enlazar(scs_m, nb, f, fp, t, tp, xs, ys, etiqueta)
    nb = nb + 1;
    l = scicos_link(xx = xs(:), yy = ys(:), ct = [1, 1], from = [f, fp, 0], to = [t, tp, 1]);
    if etiqueta <> "" then l.id = etiqueta; end
    scs_m.objs(nb) = l;
    g = scs_m.objs(f).graphics; g.pout(fp) = nb; scs_m.objs(f).graphics = g;
    g = scs_m.objs(t).graphics; g.pin(tp)  = nb; scs_m.objs(t).graphics = g;
endfunction
function [scs_m, nb] = enlazar_ev(scs_m, nb, f, fp, t, tp, xs, ys)
    nb = nb + 1;
    l = scicos_link(xx = xs(:), yy = ys(:), ct = [5, -1], from = [f, fp, 0], to = [t, tp, 1]);
    scs_m.objs(nb) = l;
    g = scs_m.objs(f).graphics; g.peout(fp) = nb; scs_m.objs(f).graphics = g;
    g = scs_m.objs(t).graphics; g.pein(tp)  = nb; scs_m.objs(t).graphics = g;
endfunction

// Malla de armadura
[scs_m,nb] = enlazar(scs_m,nb, 1,1, 2,1, [60 110], [310 310], "");
[scs_m,nb] = enlazar(scs_m,nb, 2,1, 3,1, [150 190], [320 320], "");
[scs_m,nb] = enlazar(scs_m,nb, 3,1, 4,1, [230 270], [320 320], "");
[scs_m,nb] = enlazar(scs_m,nb, 4,1, 5,1, [310 350], [320 320], "ia");
[scs_m,nb] = enlazar(scs_m,nb, 5,1, 7,1, [353 390], [320 320], "");
[scs_m,nb] = enlazar(scs_m,nb, 5,2, 6,1, [353 353 310], [320 420 420], "");
[scs_m,nb] = enlazar(scs_m,nb, 6,1, 2,3, [270 90 90 110], [420 420 335 335], "");
[scs_m,nb] = enlazar(scs_m,nb, 5,3, 17,1, [353 353 920], [320 135 135], "");
// Par electromagnético
[scs_m,nb] = enlazar(scs_m,nb, 7,1, 8,1, [430 470], [320 320], "Tm");
[scs_m,nb] = enlazar(scs_m,nb, 8,1, 9,2, [473 520], [320 320], "");
[scs_m,nb] = enlazar(scs_m,nb, 8,2, 17,2, [473 473 920], [320 150 150], "");
[scs_m,nb] = enlazar(scs_m,nb, 10,1, 9,1, [540 505 505 520], [220 220 305 305], "");
// Parte mecánica
[scs_m,nb] = enlazar(scs_m,nb, 9,1, 11,1, [560 600], [320 320], "");
[scs_m,nb] = enlazar(scs_m,nb, 11,1, 12,1, [640 680], [320 320], "");
[scs_m,nb] = enlazar(scs_m,nb, 12,1, 13,1, [720 760], [320 320], "w");
[scs_m,nb] = enlazar(scs_m,nb, 13,1, 16,1, [763 790], [320 320], "");
[scs_m,nb] = enlazar(scs_m,nb, 13,2, 14,1, [763 763 640], [320 420 420], "");
[scs_m,nb] = enlazar(scs_m,nb, 14,1, 9,3, [600 505 505 520], [420 420 335 335], "");
[scs_m,nb] = enlazar(scs_m,nb, 13,3, 15,1, [763 763 440], [320 500 500], "");
[scs_m,nb] = enlazar(scs_m,nb, 15,1, 2,2, [400 80 80 110], [500 500 320 320], "");
// Salidas
[scs_m,nb] = enlazar(scs_m,nb, 16,1, 17,3, [830 873 873 920], [320 320 165 165], "n [rpm]");
[scs_m,nb] = enlazar(scs_m,nb, 17,1, 18,1, [930 957], [150 150], "");
[scs_m,nb] = enlazar(scs_m,nb, 18,1, 19,1, [960 1000], [150 150], "");
[scs_m,nb] = enlazar(scs_m,nb, 18,2, 21,1, [960 960 1000], [150 250 250], "");
[scs_m,nb] = enlazar(scs_m,nb, 19,1, 20,1, [1010 1060], [135 135], "");
[scs_m,nb] = enlazar(scs_m,nb, 19,2, 20,2, [1010 1060], [150 150], "");
[scs_m,nb] = enlazar(scs_m,nb, 19,3, 20,3, [1010 1060], [165 165], "");
// Reloj de muestreo (dt) para osciloscopio y exportación
[scs_m,nb] = enlazar_ev(scs_m,nb, 22,1, 23,1, [1090 1090], [50 77]);
[scs_m,nb] = enlazar_ev(scs_m,nb, 23,1, 20,1, [1090 1090], [84 120]);
[scs_m,nb] = enlazar_ev(scs_m,nb, 23,2, 21,1, [1090 1140 1140 1030 1030], [84 84 210 210 230]);

// ---------------- Guardar ----------------
xcosDiagramToScilab(archivo, scs_m);
mprintf("\nDiagrama guardado en: %s  (%d objetos)\n", archivo, nb);
