"""Testa se o polo usado no frame EQUADOR_JUPITER (polo IAU do pck) e' o
mesmo eixo em torno do qual as luas do SPICE precessionam.

Ideia: o vetor inclinacao (p, q) = (i sin Om, i cos Om) de uma lua descreve
um circulo (precessao nodal) em torno do vetor da inclinacao forcada. Para
as luas pequenas, dominadas pelo J2, o plano forcado e' o equador de
Jupiter, entao o circulo tem que ser centrado na origem se o eixo z do
frame for o polo verdadeiro. Se o polo IAU estiver deslocado de delta do
polo dinamico da efemeride, todos os circulos ficam deslocados do mesmo
vetor (e a inclinacao do SPICE oscila entre |i - delta| e i + delta, com
o periodo nodal, enquanto a do nbody fica constante).

Para cada lua: ajusta um circulo aos pontos (p, q) do SPICE (e do nbody, se
houver NOME_osc.dat) e imprime o centro. Centros iguais nas luas pequenas e
diferentes de zero = polo deslocado. Salva polo_spice.png.

Precisa de pelo menos uma volta do nodo no periodo do CSV (Metis/Adrastea
~45 d, Amalthea ~150 d, Thebe ~300 d; Io bem mais).

Uso:
    python3 polo_spice.py                              # Metis Adrastea Amalthea Thebe Io
    python3 polo_spice.py Metis Thebe --pasta-spice saida_elementos --pasta-nbody .
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")                # so' salva a imagem, sem abrir janela
import matplotlib.pyplot as plt
import numpy as np

from compara_spice import CORES, arquivo_spice, le_nbody, le_spice

LUAS = ["Metis", "Adrastea", "Amalthea", "Thebe", "Io"]


def vetor_inclinacao(i, om):
    """(p, q) = (i sin Om, i cos Om), em graus."""
    om = np.radians(om)
    return i*np.sin(om), i*np.cos(om)


def ajusta_circulo(x, y):
    """Circulo de minimos quadrados (Kasa): -> (xc, yc, raio)."""
    A = np.column_stack([x, y, np.ones_like(x)])
    D, E, F = np.linalg.lstsq(A, -(x*x + y*y), rcond=None)[0]
    xc, yc = -D/2, -E/2
    return xc, yc, np.sqrt(max(xc*xc + yc*yc - F, 0.0))


def voltas_do_nodo(om):
    return abs(np.ptp(np.unwrap(np.radians(om))))/(2*np.pi)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("nomes", nargs="*", default=LUAS, help="luas (padrao: as pequenas e Io)")
    p.add_argument("--pasta-spice", default="saida_elementos", help="pasta dos elementos_NOME.csv")
    p.add_argument("--pasta-nbody", default=".", help="pasta dos NOME_osc.dat")
    a = p.parse_args()

    resultados = []
    fig, axes = plt.subplots(1, len(a.nomes), figsize=(4.2*len(a.nomes), 4.6), squeeze=False)
    for ax, nome in zip(axes[0], a.nomes):
        arq_sp = arquivo_spice(nome, a.pasta_spice)
        if arq_sp is None:
            print(f"{nome}: sem CSV do SPICE em {a.pasta_spice}, pulando")
            ax.set_visible(False)
            continue
        _, el = le_spice(arq_sp)
        series = {"SPICE": (el[:, 2], el[:, 3])}                 # i, Omega
        arq_nb = os.path.join(a.pasta_nbody, f"{nome}_osc.dat")
        if os.path.exists(arq_nb):
            nb = le_nbody(arq_nb)
            series["nbody"] = (nb[:, 3], nb[:, 5])               # i, om

        for metodo, (i, om) in series.items():
            x, y = vetor_inclinacao(i, om)
            xc, yc, raio = ajusta_circulo(x, y)
            voltas = voltas_do_nodo(om)
            ax.plot(x, y, ".", ms=1, color=CORES[metodo], alpha=0.5, label=metodo)
            ax.plot(xc, yc, "+", ms=12, mew=2, color=CORES[metodo])
            aviso = "" if voltas >= 1 else "  (menos de 1 volta do nodo: centro incerto)"
            print(f"{nome:9s} {metodo:5s}  centro = ({xc:+.5f}, {yc:+.5f})°  "
                  f"|centro| = {np.hypot(xc, yc):.5f}°  raio = {raio:.5f}°  "
                  f"voltas do nodo = {voltas:.1f}{aviso}")
            if metodo == "SPICE" and voltas >= 1 and nome.lower() != "io":
                resultados.append((xc, yc))

        ax.plot(0, 0, "k+", ms=8)
        ax.set_aspect("equal", "datalim")
        ax.set_title(nome)
        ax.set_xlabel("$i \\sin\\Omega$ (°)")
        ax.grid(True, linestyle="--", alpha=0.4)
    axes[0, 0].set_ylabel("$i \\cos\\Omega$ (°)")
    axes[0, 0].legend(markerscale=8, frameon=False)
    fig.suptitle("Vetor inclinação — o '+' colorido é o centro do círculo ajustado", fontsize=12)
    fig.tight_layout()
    fig.savefig("polo_spice.png", dpi=200, bbox_inches="tight")
    print("Salvo: polo_spice.png")

    if len(resultados) >= 2:
        c = np.array(resultados)
        xc, yc = c.mean(axis=0)
        dispersao = c.std(axis=0).max()
        print(f"\nCentro medio das luas pequenas (SPICE): ({xc:+.5f}, {yc:+.5f})°, "
              f"|delta| = {np.hypot(xc, yc):.5f}°, dispersao = {dispersao:.5f}°")
        if np.hypot(xc, yc) > 3*dispersao:
            # normal do plano forcado no frame: (sin i sin Om, -sin i cos Om, cos i)
            print("=> Os circulos estao deslocados juntos: o polo do frame (IAU) nao e' o eixo\n"
                  "   em torno do qual as luas precessionam na efemeride. Polo dinamico no\n"
                  f"   frame EQUADOR_JUPITER ~ (x, y, z) = ({np.radians(xc):+.3e}, "
                  f"{-np.radians(yc):+.3e}, 1).")
        else:
            print("=> Centros compativeis com zero (ou espalhados): o polo nao explica a diferenca.")


if __name__ == "__main__":
    main()
