"""Converte o estados_iniciais.csv do elementos_spice.py (estados cartesianos)
em um arquivo .in do nbody (elementos osculadores).

Le:  nome, naif_id, GM_km3_s2, x, y, z (km), vx, vy, vz (km/s)
     + o GM_JUPITER e o frame da primeira linha (comentario)
Gera um .in completo, com:
     - gm_planeta = GM do CSV  (o mesmo mu que o SPICE usou)
     - unidade_massa = gm      (massas = GM de cada lua, em km^3/s^2)
     - um satelite por linha: GM a e i w Om M (km, graus)

A conversao estado -> elementos e' a mesma do cart_para_osc do nbody.f90,
com mu = GM_planeta + GM_lua, entao o nbody reconstroi exatamente o estado
do SPICE em t = 0. O frame do CSV tem que ser o EQUADOR_JUPITER (o nbody
poe o J2 no eixo z).

Uso:
    python3 spice_para_in.py saida_elementos/estados_iniciais.csv
    python3 spice_para_in.py estados_iniciais.csv -o jupiter_spice.in --tempo-final 365 --passo 0.01
"""
import argparse

import numpy as np


# ------------------------------------------------------------- conversao
def cart_para_osc(mu, r, v):
    """Estado -> (a, e, i, w, Om, M), angulos em graus [0, 360).
    Mesmas formulas do cart_para_osc do nbody.f90."""
    rr = np.linalg.norm(r)
    v2 = v @ v
    h = np.cross(r, v)
    hh = np.linalg.norm(h)

    a = 1.0/(2.0/rr - v2/mu)
    if a <= 0.0:
        raise ValueError("orbita nao eliptica")
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


def osc_para_cart(mu, a, e, inc, w, om, m):
    """Inverso (como o osc_para_cart do nbody.f90), so' para conferir."""
    inc, w, om, m = np.radians([inc, w, om, m])
    u = m
    for _ in range(100):
        du = (u - e*np.sin(u) - m)/(1.0 - e*np.cos(u))
        u -= du
        if abs(du) < 1e-15:
            break
    n = np.sqrt(mu/a**3)
    fac = np.sqrt(1.0 - e*e)
    xo, yo = a*(np.cos(u) - e), a*fac*np.sin(u)
    vxo = -n*a*np.sin(u)/(1.0 - e*np.cos(u))
    vyo = n*a*fac*np.cos(u)/(1.0 - e*np.cos(u))
    cw, sw, co, so, ci, si = np.cos(w), np.sin(w), np.cos(om), np.sin(om), np.cos(inc), np.sin(inc)
    p = np.array([co*cw - so*sw*ci, so*cw + co*sw*ci, sw*si])
    q = np.array([-co*sw - so*cw*ci, -so*sw + co*cw*ci, cw*si])
    return xo*p + yo*q, vxo*p + vyo*q


# ---------------------------------------------------------------- leitura
def le_csv(arq):
    with open(arq) as f:
        cab = f.readline().lstrip("#").split()
    info = dict(c.split("=", 1) for c in cab if "=" in c)
    d = np.genfromtxt(arq, delimiter=",", comments="#", skip_header=2, dtype=None, encoding=None)
    luas = [(str(l[0]), float(l[2]), np.array([l[3], l[4], l[5]], float),
             np.array([l[6], l[7], l[8]], float)) for l in np.atleast_1d(d)]
    return info, luas


# ------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", nargs="?", default="saida_elementos/estados_iniciais.csv",
                   help="estados_iniciais.csv do elementos_spice.py")
    p.add_argument("-o", "--saida", default="jupiter_spice.in", help="arquivo .in gerado")
    p.add_argument("--planeta", default="Jupiter", help="nome em planetas.dat (raio, J2, J4)")
    p.add_argument("--tempo-final", type=float, default=365.0, help="dias (padrao: 365)")
    p.add_argument("--passo", type=float, default=0.01, help="passo de saida em dias (padrao: 0.01)")
    a = p.parse_args()

    info, luas = le_csv(a.csv)
    gm_p = float(info["GM_JUPITER"])
    if info.get("frame") != "EQUADOR_JUPITER":
        print(f"AVISO: frame do CSV = {info.get('frame')}; o nbody precisa do EQUADOR_JUPITER")

    linhas = []
    for nome, gm, r, v in luas:
        mu = gm_p + gm
        el = cart_para_osc(mu, r, v)
        r2, v2 = osc_para_cart(mu, *el)
        erro = max(np.abs(r2 - r).max(), 86400*np.abs(v2 - v).max())
        linhas.append(f"{nome.capitalize():10s} {gm:22.15e} {el[0]:22.15f} {el[1]:22.15e} "
                      f"{el[2]:20.15f} {el[3]:20.15f} {el[4]:20.15f} {el[5]:20.15f}")
        print(f"{nome:9s} a = {el[0]:12.3f} km  e = {el[1]:.6f}  i = {el[2]:.6f}°   "
              f"(ida e volta: erro {erro:.1e} km)")

    with open(a.saida, "w") as f:
        f.write(f"""! Gerado por spice_para_in.py a partir de {a.csv}
! epoca (t = 0): {info.get('epoca_utc', '?')} UTC   et = {info.get('et_tdb_s', '?')} s TDB
! frame: {info.get('frame', '?')}   centro: {info.get('centro', '?')}

planeta    = {a.planeta}
gm_planeta = {gm_p:.15g}     ! km^3/s^2, o mesmo GM do SPICE

tempo_final = {a.tempo_final}
passo_saida = {a.passo}
precisao    = 12

saida_osculadores = sim
saida_geometricos = sim
pasta_saida       = .

unidade_massa = gm          ! massas abaixo = G*m em km^3/s^2

satelites:
! nome     GM(km^3/s^2)           a(km)                  e                      i(deg)               w(deg)               Om(deg)              M(deg)
""")
        f.write("\n".join(linhas) + "\nfim\n")
    print(f"\nSalvo: {a.saida}")


if __name__ == "__main__":
    main()
