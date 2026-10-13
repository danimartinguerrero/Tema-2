#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Electromagnetismo I - Tema 2: Desarrollo multipolar
====================================================

Generador de los datos numéricos de las figuras «tipo Mathematica» del tema
(dipolo y cuadrupolo, puntuales y físicos). El script NO dibuja nada: solo
calcula curvas y las escribe como llamadas a macros TikZ en
``figuras/datos/*.tex``; el estilo (colores, grosores, puntas de flecha) se
fija en LaTeX, en cada fichero ``figuras/fig_*.tex``.

Se generan tres tipos de objetos, que imitan a ContourPlot, StreamPlot y
VectorPlot de Mathematica:

* ``\\EQ{(x,y) (x,y) ...}``    equipotenciales (curvas de nivel de phi);
* ``\\LC{(x,y) (x,y) ...}``    líneas de campo con separación uniforme
                              (algoritmo de Jobard-Lefer simplificado);
* ``\\FL{x}{y}{ang}``          punta de flecha sobre una línea de campo;
* ``\\VF{x}{y}{dx}{dy}{t}``    flecha del campo vectorial (t en [0,1] = color).

Unidades: 1/(4 pi eps0) = 1, q = 1 y longitudes en unidades de l (separación
de las cargas), como en las figuras del tema. El plano representado es el
meridiano y = 0 (coordenadas x, z); en el ejemplo del cuadrupolo plano es el
plano z = 0 (coordenadas x, y).

Uso:  python3 scripts/generar_datos.py      (desde la carpeta latex/)
Requiere numpy y matplotlib (solo se usa su rutina de curvas de nivel).
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

AQUI = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(AQUI, "..", "figuras", "datos")
os.makedirs(SALIDA, exist_ok=True)


# ---------------------------------------------------------------------------
#  Campos (plano meridiano: X = coordenada transversal, Z = eje de simetría)
# ---------------------------------------------------------------------------
def cargas(lista):
    """Campo y potencial de un sistema de cargas puntuales [(q, x, z), ...]
    evaluados en el plano que las contiene (ley de Coulomb en 3D)."""
    def E(X, Z):
        ex = np.zeros_like(X, dtype=float)
        ez = np.zeros_like(X, dtype=float)
        for q, x0, z0 in lista:
            dx, dz = X - x0, Z - z0
            r3 = (dx * dx + dz * dz) ** 1.5
            ex += q * dx / r3
            ez += q * dz / r3
        return ex, ez

    def phi(X, Z):
        v = np.zeros_like(X, dtype=float)
        for q, x0, z0 in lista:
            v += q / np.hypot(X - x0, Z - z0)
        return v
    return E, phi


def dipolo_puntual():
    """p = 1 u_z en el origen:  phi = z/r^3,  E = [3 (p.u_r) u_r - p]/r^3."""
    def E(X, Z):
        r2 = X * X + Z * Z
        r5 = r2 ** 2.5
        return 3 * X * Z / r5, (3 * Z * Z - r2) / r5

    def phi(X, Z):
        return Z / (X * X + Z * Z) ** 1.5
    return E, phi


def cuadrupolo_puntual():
    """Término cuadrupolar con simetría de revolución (eje z):
    phi = P2(cos t)/r^3 = (2 z^2 - x^2)/(2 r^5)."""
    def phi(X, Z):
        r2 = X * X + Z * Z
        return (2 * Z * Z - X * X) / (2 * r2 ** 2.5)

    def E(X, Z):
        r2 = X * X + Z * Z
        r = np.sqrt(r2)
        c, s = Z / r, X / r
        er = 3 * (1.5 * c * c - 0.5) / r ** 4      # E_r  = 3 P2 / r^4
        et = 3 * s * c / r ** 4                    # E_th = 3 sen cos / r^4
        # u_r = (s, c), u_th = (c, -s) en el plano (x, z)
        return er * s + et * c, er * c - et * s
    return E, phi


# ---------------------------------------------------------------------------
#  Líneas de campo con separación uniforme (Jobard & Lefer, 1997)
# ---------------------------------------------------------------------------
class Rejilla:
    """Tabla espacial para consultar rápidamente la distancia a las líneas."""

    def __init__(self, L, celda):
        self.L, self.c = L, celda
        self.n = int(np.ceil(2 * L / celda)) + 1
        self.cel = {}

    def _k(self, x, y):
        return (int((x + self.L) / self.c), int((y + self.L) / self.c))

    def libre(self, x, y, d, ignorar=None):
        i, j = self._k(x, y)
        for a in (i - 1, i, i + 1):
            for b in (j - 1, j, j + 1):
                for (px, py, idl) in self.cel.get((a, b), ()):
                    if idl == ignorar:
                        continue
                    if (px - x) ** 2 + (py - y) ** 2 < d * d:
                        return False
        return True

    def anade(self, pts, idl):
        for x, y in pts:
            self.cel.setdefault(self._k(x, y), []).append((x, y, idl))


def lineas_uniformes(E, L=10.0, r0=1.2, dsep=0.95, dtest=0.45, h=0.025,
                     n_circ=24, semillas_extra=(), polos=(), malla=True,
                     dsep_malla=None, dtest_malla=None):
    """Traza líneas de campo con separación aproximadamente uniforme en
    el cuadrado [-L, L]^2, excluyendo el disco r < r0 (fuente)."""
    def dir_(x, y):
        ex, ey = E(np.array(x), np.array(y))
        m = np.hypot(ex, ey)
        if not np.isfinite(m) or m < 1e-14:
            return None
        return ex / m, ey / m

    def dentro(x, y):
        if any((x - px) ** 2 + (y - py) ** 2 < rp * rp for px, py, rp in polos):
            return False
        return abs(x) <= L and abs(y) <= L and x * x + y * y >= r0 * r0

    rej = Rejilla(L, max(dsep, dsep_malla or 0.0))
    lineas = []

    def integra(x0, y0, sgn, idl, dtest=dtest):
        pts = []
        x, y = x0, y0
        for _ in range(int(6 * L / h)):
            k1 = dir_(x, y)
            if k1 is None:
                break
            k2 = dir_(x + 0.5 * h * sgn * k1[0], y + 0.5 * h * sgn * k1[1])
            if k2 is None:
                break
            k3 = dir_(x + 0.5 * h * sgn * k2[0], y + 0.5 * h * sgn * k2[1])
            if k3 is None:
                break
            k4 = dir_(x + h * sgn * k3[0], y + h * sgn * k3[1])
            if k4 is None:
                break
            x += sgn * h * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
            y += sgn * h * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
            if not dentro(x, y):
                # se recorta en el borde para que la línea llegue al marco
                pts.append((max(-L, min(L, x)), max(-L, min(L, y))))
                break
            if not rej.libre(x, y, dtest, ignorar=idl):
                break
            # evita lazos: no volver cerca de un tramo propio ya lejano
            if len(pts) > 60 and len(pts) % 10 == 0:
                A = np.array(pts[:-40])
                if np.min((A[:, 0] - x) ** 2 + (A[:, 1] - y) ** 2) < dtest * dtest:
                    break
            pts.append((x, y))
        return pts

    def lanza(x0, y0, ds=dsep, dt=dtest):
        if not dentro(x0, y0) or not rej.libre(x0, y0, ds):
            return
        idl = len(lineas)
        adelante = integra(x0, y0, +1, idl, dt)
        atras = integra(x0, y0, -1, idl, dt)
        pts = atras[::-1] + [(x0, y0)] + adelante
        if len(pts) * h < 0.8:
            return
        rej.anade(pts, idl)
        lineas.append(pts)

    # 1) semillas en la frontera del disco central (las líneas «nacen» allí)
    for k in range(n_circ):
        a = 2 * np.pi * (k + 0.5) / n_circ
        lanza((r0 + 0.02) * np.cos(a), (r0 + 0.02) * np.sin(a))
    for (x0, y0) in semillas_extra:
        lanza(x0, y0)
    # 2) semillas en una malla fina, recorrida en orden fijo (reproducible)
    if malla:
        dm = dsep_malla or dsep
        tm = dtest_malla or dtest
        g = np.arange(-L + 0.35 * dm, L, 0.5 * dm)
        orden = np.random.default_rng(7).permutation(len(g) * len(g))
        for idx in orden:
            lanza(g[idx // len(g)], g[idx % len(g)], dm, tm)
    return lineas


def flechas_sobre(lineas, E, paso=2.6, minimo=0.9):
    """Puntas de flecha repartidas por longitud de arco a lo largo de cada línea."""
    out = []
    for pts in lineas:
        P = np.array(pts)
        if len(P) < 3:
            continue
        s = np.concatenate([[0], np.cumsum(np.hypot(*np.diff(P, axis=0).T))])
        Ltot = s[-1]
        if Ltot < minimo:
            continue
        n = max(1, int(Ltot / paso))
        for k in range(n):
            sk = Ltot * (k + 0.5) / n
            i = int(np.searchsorted(s, sk))
            i = min(max(i, 1), len(P) - 1)
            x, y = P[i]
            ex, ey = E(np.array(x), np.array(y))
            out.append((x, y, np.degrees(np.arctan2(ey, ex))))
    return out


# ---------------------------------------------------------------------------
#  Equipotenciales (curvas de nivel) y campo vectorial
# ---------------------------------------------------------------------------
def equipotenciales(phi, niveles, L=10.0, r0=1.2, n=801):
    x = np.linspace(-L, L, n)
    X, Z = np.meshgrid(x, x)
    with np.errstate(divide="ignore", invalid="ignore"):
        V = phi(X, Z)
    V = np.ma.masked_where(X * X + Z * Z < r0 * r0, V)
    fig = plt.figure()
    cs = plt.contour(X, Z, V, levels=sorted(niveles))
    curvas = [seg for segs in cs.allsegs for seg in segs if len(seg) > 2]
    plt.close(fig)
    return curvas


def campo_vectorial(E, L=10.0, n=21, r0=1.4, escala=1.3):
    """Malla de flechas: longitud ~ |E|^(1/3) (saturada) y color ~ log|E|."""
    g = np.linspace(-L, L, n)
    X, Z = np.meshgrid(g, g)
    ex, ez = E(X, Z)
    m = np.hypot(ex, ez)
    mask = X * X + Z * Z >= r0 * r0
    lm = np.log10(m[mask])
    lo, hi = np.percentile(lm, 3), np.percentile(lm, 97)
    out = []
    for x, z, a, b, mm, ok in zip(X.ravel(), Z.ravel(), ex.ravel(), ez.ravel(),
                                  m.ravel(), mask.ravel()):
        if not ok or mm == 0:
            continue
        t = float(np.clip((np.log10(mm) - lo) / (hi - lo), 0, 1))
        long_ = escala * (0.30 + 0.70 * t)            # longitud de la flecha
        out.append((x, z, long_ * a / mm, long_ * b / mm, t))
    return out


# ---------------------------------------------------------------------------
#  Escritura
# ---------------------------------------------------------------------------
def simplifica(pts, d=0.05):
    out = [pts[0]]
    for p in pts[1:-1]:
        if (p[0] - out[-1][0]) ** 2 + (p[1] - out[-1][1]) ** 2 >= d * d:
            out.append(p)
    out.append(pts[-1])
    return out


def coords(pts):
    return " ".join("(%.3f,%.3f)" % (x, y) for x, y in simplifica(list(pts)))


def escribe(nombre, cabecera, cuerpo):
    with open(os.path.join(SALIDA, nombre), "w", encoding="utf-8") as f:
        f.write("%% %s\n%% Generado por scripts/generar_datos.py (no editar a mano).\n"
                % cabecera)
        f.write("\n".join(cuerpo) + "\n")
    print("  ->", nombre, "(%d objetos)" % len(cuerpo))


def exporta(prefijo, titulo, E, phi, niveles, L=10.0, r0_eq=1.2, r0_lc=1.2,
            dsep=0.95, vect=True, eq=True, **kw):
    if eq:
        cur = equipotenciales(phi, niveles, L=L, r0=r0_eq)
        escribe(prefijo + "_eq.tex", titulo + " - equipotenciales",
                ["\\EQ{%s}" % coords(c) for c in cur])
    lin = lineas_uniformes(E, L=L, r0=r0_lc, dsep=dsep, **kw)
    escribe(prefijo + "_lc.tex", titulo + " - líneas de campo",
            ["\\LC{%s}" % coords(c) for c in lin])
    fl = flechas_sobre(lin, E)
    escribe(prefijo + "_fl.tex", titulo + " - puntas de flecha",
            ["\\FL{%.3f}{%.3f}{%.1f}" % f for f in fl])
    if vect:
        vf = campo_vectorial(E, L=L, r0=r0_lc + 0.2)
        escribe(prefijo + "_vf.tex", titulo + " - campo vectorial",
                ["\\VF{%.3f}{%.3f}{%.3f}{%.3f}{%.3f}" % v for v in vf])


def exporta_psi(nombre, titulo, Psi, niveles, L=3.0, r0=1.0, n=601):
    """Líneas de campo exactas de un sistema con simetría de revolución como
    curvas de nivel de la función de corriente de Stokes Psi (flujo constante
    entre dos líneas consecutivas)."""
    cur = equipotenciales(Psi, niveles, L=L, r0=r0, n=n)
    escribe(nombre, titulo, ["\\LC{%s}" % coords(c) for c in cur])


def tabla_discrepancia(nombre, rmin=1.2, rmax=15.0, n=280):
    """Discrepancia relativa |E_puntual - E_fisico| / |E_fisico| (en %) a lo
    largo de las semirrectas theta = 0, pi/4 y pi/2 (Proposición de la
    discrepancia entre el dipolo físico y el puntual)."""
    Ef, _ = cargas([(1, 0, 0.5), (-1, 0, -0.5)])
    Ep, _ = dipolo_puntual()
    r = np.linspace(rmin, rmax, n)
    filas = []
    for th in (0.0, np.pi / 4, np.pi / 2):
        X, Z = r * np.sin(th), r * np.cos(th)
        fx, fz = Ef(X, Z)
        px, pz = Ep(X, Z)
        filas.append(100 * np.hypot(px - fx, pz - fz) / np.hypot(fx, fz))
    with open(os.path.join(SALIDA, nombre), "w", encoding="utf-8") as f:
        f.write("r e0 e45 e90\n")
        for k in range(n):
            f.write("%.4f %.5f %.5f %.5f\n" % (r[k], filas[0][k], filas[1][k], filas[2][k]))
    print("  ->", nombre)


def tabla_error_carga(nombre, n=300):
    """Error relativo (%) al aproximar el potencial de una carga q en a u_z por
    sus términos monopolar y dipolar, en función de r/a (escala logarítmica),
    para theta = 0, pi/4 y pi/2 (sección de estimación del error)."""
    x = np.logspace(np.log10(1.5), 3, n)          # r/a
    cols = []
    for th in (0.0, np.pi / 4, np.pi / 2):
        c = np.cos(th)
        exacto = 1 / np.sqrt(x * x + 1 - 2 * x * c)
        aprox = 1 / x + c / x ** 2
        cols.append(100 * np.abs(exacto - aprox) / exacto)
    with open(os.path.join(SALIDA, nombre), "w", encoding="utf-8") as f:
        f.write("x e0 e45 e90\n")
        for k in range(n):
            f.write("%.5f %.6e %.6e %.6e\n" % (x[k], cols[0][k], cols[1][k], cols[2][k]))
    print("  ->", nombre)


def tabla_error_barra(nombre, n=300):
    """Problema 2.1: error relativo (%) del potencial de una barra de longitud
    l = 1 (lambda/(4 pi eps0) = 1) aproximado por el término monopolar (e0) y
    por monopolar + cuadrupolar (e2), sobre el eje (theta = 0) y en el plano
    ecuatorial (theta = pi/2), frente al potencial exacto."""
    x = np.logspace(0, 2, n)                      # r / l
    exe = np.log((x + 0.5) / (x - 0.5))           # theta = 0 (eje)
    exq = 2 * np.arcsinh(0.5 / x)                 # theta = pi/2 (ecuador)
    phi0 = 1 / x
    e0e = 100 * np.abs(exe - phi0) / exe
    e0q = 100 * np.abs(exq - phi0) / exq
    e2e = 100 * np.abs(exe - phi0 - 2 / (24 * x ** 3)) / exe
    e2q = 100 * np.abs(exq - phi0 + 1 / (24 * x ** 3)) / exq
    with open(os.path.join(SALIDA, nombre), "w", encoding="utf-8") as f:
        f.write("x e0eje e0ecu e2eje e2ecu\n")
        for k in range(1, n):                     # x = 1 excluido (singular en el eje)
            f.write("%.5f %.6e %.6e %.6e %.6e\n" % (x[k], e0e[k], e0q[k], e2e[k], e2q[k]))
    print("  ->", nombre)


def niveles_simetricos(paso, n):
    return [k * paso for k in range(1, n + 1)] + [-k * paso for k in range(1, n + 1)]


if __name__ == "__main__":
    print("Generando datos en", os.path.normpath(SALIDA))

    # Dipolo puntual (p = q l u_z con q = l = 1)
    E, phi = dipolo_puntual()
    exporta("dip_puntual", "Dipolo puntual", E, phi,
            niveles_simetricos(0.0105, 20), r0_eq=1.3, r0_lc=1.2, dsep=0.88)

    # Dipolo físico: +q en z = +l/2, -q en z = -l/2
    E, phi = cargas([(1, 0, 0.5), (-1, 0, -0.5)])
    exporta("dip_fisico", "Dipolo físico", E, phi,
            niveles_simetricos(0.0105, 20), r0_eq=1.3, r0_lc=1.2, dsep=0.88, vect=False)

    # Término cuadrupolar (cuadrupolo puntual con simetría de revolución)
    E, phi = cuadrupolo_puntual()
    exporta("cua_puntual", "Cuadrupolo puntual", E, phi,
            niveles_simetricos(0.00125, 17) + [0.0], r0_eq=1.3, r0_lc=1.2, dsep=0.88)

    # Cuadrupolo físico lineal: +2 en z = +-1/2, -4 en el origen (Q_zz = 4 q d^2)
    E, phi = cargas([(2, 0, 0.5), (2, 0, -0.5), (-4, 0, 0)])
    exporta("cua_fisico", "Cuadrupolo físico (lineal)", E, phi,
            niveles_simetricos(0.00125, 17) + [0.0], r0_eq=1.3, r0_lc=1.2, dsep=0.88, vect=False)

    # Esquema «dipolo físico frente a dipolo puntual» (fuera de una esfera r < 1)
    #   puntual: Psi = sen^2(t)/r = x^2/r^3 ;  físico: Psi = cos(a-) - cos(a+)
    niv = [0.04, 0.10, 0.17, 0.25, 0.35, 0.48, 0.62]
    exporta_psi("esq_dip_puntual_lc.tex", "Dipolo puntual - función de corriente",
                lambda X, Z: X * X / (X * X + Z * Z) ** 1.5, niv)
    exporta_psi("esq_dip_fisico_lc.tex", "Dipolo físico - función de corriente",
                lambda X, Z: (Z + 0.5) / np.hypot(X, Z + 0.5)
                - (Z - 0.5) / np.hypot(X, Z - 0.5), niv)

    # Error de truncamiento para una carga puntual fuera del origen
    tabla_error_carga("error_carga.dat")

    # Problema 2.1: error de la aproximación multipolar de la barra
    tabla_error_barra("error_barra.dat")

    # Discrepancia entre el dipolo físico y el puntual
    tabla_discrepancia("discrepancia.dat")

    # Cuadrupolo plano del ejemplo: +q en (+-a, 0), -q en (0, +-a), plano z = 0.
    # Líneas sembradas solo alrededor de las cargas positivas, a ángulos iguales
    # (flujo igual entre dos líneas consecutivas cerca de la carga).
    E, phi = cargas([(1, 1, 0), (1, -1, 0), (-1, 0, 1), (-1, 0, -1)])
    ang = np.linspace(0, 2 * np.pi, 18, endpoint=False) + np.pi / 18
    semillas = [(sx + 0.1 * np.cos(a), 0.1 * np.sin(a)) for sx in (1, -1) for a in ang]
    exporta("cua_plano", "Cuadrupolo plano (cuatro cargas)", E, phi, [0.0],
            L=4.0, r0_eq=0.0, r0_lc=0.0, dsep=0.02, dtest=0.005, h=0.01, n_circ=0,
            vect=False, eq=False, semillas_extra=semillas, dsep_malla=0.62, dtest_malla=0.3,
            polos=[(1, 0, 0.08), (-1, 0, 0.08), (0, 1, 0.08), (0, -1, 0.08)])
