"""Converte o estado cartesiano do Horizons (Vector Table) em elementos
osculadores: a e i w om M. Script independente (so' precisa do numpy).

Entrada: o trecho do Horizons com X, Y, Z, VX, VY, VZ, por ex.
     X = 1.234567890123456E+09 Y =-2.345678901234567E+09 Z = 3.456789012345678E+08
     VX= 1.234567890123456E+00 VY= 2.345678901234567E+00 VZ=-3.456789012345678E-01
Pode estar num arquivo ou ser colado no terminal. Se tiver mais de um
registro, so' o primeiro e' usado. As unidades (km e km/s, km e km/dia ou
au e au/dia) sao reconhecidas sozinhas e o script diz qual usou.

Uso:
    python3 horizons_para_in.py sol.txt
    python3 horizons_para_in.py            # cola o texto e aperta Ctrl+D

Saida (uma linha): a(km) e i(deg) w(deg) om(deg) M(deg)
"""
import re
import sys

import numpy as np

# ======================= CONFIGURACAO (edite aqui) =======================
GM_PLANETA = 6835099.5            # km^3/s^2  (Netuno; Jupiter = 126686531.9)
GM_CORPO = 132712440041.279419    # km^3/s^2  (Sol; para uma lua, o GM dela)
# =========================================================================

UA_KM = 149597870.7
DIA = 86400.0


def le_estado(texto):
    """X, Y, Z (km) e VX, VY, VZ (km/s) do primeiro registro do texto."""
    valores = {}
    for chave in ("X", "Y", "Z", "VX", "VY", "VZ"):
        m = re.search(rf"(?<![A-Z]){chave}\s*=\s*([-+]?[\d.]+(?:[Ee][-+]?\d+)?)", texto)
        if m is None:
            sys.exit(f"nao achei '{chave} =' no texto")
        valores[chave] = float(m.group(1))
    r = np.array([valores["X"], valores["Y"], valores["Z"]])
    v = np.array([valores["VX"], valores["VY"], valores["VZ"]])
    return r, v


def cart_para_osc(mu, r, v):
    """Estado -> (a, e, i, w, om, M), angulos em graus [0, 360).
    Mesmas formulas do cart_para_osc do nbody.f90."""
    rr = np.linalg.norm(r)
    v2 = v @ v
    h = np.cross(r, v)
    hh = np.linalg.norm(h)

    a = 1.0/(2.0/rr - v2/mu)
    if a <= 0.0:
        v_esc = np.sqrt(2.0*mu/rr)
        sys.exit(f"orbita nao eliptica: |r| = {rr:.6g} km, |v| = {np.sqrt(v2):.6g} km/s, "
                 f"mas a velocidade de escape e' {v_esc:.6g} km/s.\n"
                 "Confira GM_PLANETA e GM_CORPO no topo do script e se X..VZ sao do corpo certo.")
    esinu = (r @ v)/np.sqrt(mu*a)
    ecosu = rr*v2/mu - 1.0
    e = np.hypot(esinu, ecosu)
    u = np.arctan2(esinu, ecosu)
    m = u - e*np.sin(u)
    f = 2.0*np.arctan2(np.sqrt(1.0 + e)*np.sin(0.5*u), np.sqrt(1.0 - e)*np.cos(0.5*u))

    inc = np.arccos(np.clip(h[2]/hh, -1.0, 1.0))
    if np.hypot(h[0], h[1]) < 1e-14*hh:
        om = 0.0                                  # equatorial: nodo indefinido
        ulat = np.arctan2(r[1], r[0])
    else:
        om = np.arctan2(h[0], -h[1])
        cosu = r[0]*np.cos(om) + r[1]*np.sin(om)
        sinu = (r[1]*np.cos(om) - r[0]*np.sin(om))*np.cos(inc) + r[2]*np.sin(inc)
        ulat = np.arctan2(sinu, cosu)
    w = ulat - f
    graus = lambda x: np.degrees(x) % 360.0
    return a, e, np.degrees(inc), graus(w), graus(om), graus(m)


if len(sys.argv) > 1:
    texto = open(sys.argv[1]).read()
else:
    print("Cole o texto do Horizons (X = ... VZ = ...) e aperte Ctrl+D:", file=sys.stderr)
    texto = sys.stdin.read()

r, v = le_estado(texto)
mu = GM_PLANETA + GM_CORPO

# unidades: numa orbita real |v| fica perto da velocidade circular
# sqrt(mu/r); escolhe a unidade em que v/v_circ fica mais perto de 1
# (as alternativas erram por fatores de 86400 ou mais)
opcoes = (("km e km/s", 1.0, 1.0),
          ("km e km/dia", 1.0, 1.0/DIA),
          ("au e au/dia", UA_KM, UA_KM/DIA))
razao = lambda fr, fv: np.linalg.norm(v*fv)/np.sqrt(mu/np.linalg.norm(r*fr))
nome, fr, fv = min(opcoes, key=lambda o: abs(np.log(razao(o[1], o[2]))))
r, v = r*fr, v*fv
print(f"# unidades lidas como {nome}", file=sys.stderr)
a, e, i, w, om, M = cart_para_osc(mu, r, v)
print("# a(km)                 e                      i(deg)               w(deg)               om(deg)              M(deg)")
print(f"{a:22.15f} {e:22.15e} {i:20.15f} {w:20.15f} {om:20.15f} {M:20.15f}")
