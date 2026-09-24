"""Plota os elementos orbitais gerados pelo nbody.

Uso:
    python3 plot_elementos.py                 # todos os satelites (*_osc.dat / *_geo.dat)
    python3 plot_elementos.py Thebe Io        # so' os satelites escolhidos
    python3 plot_elementos.py Thebe --pasta resultados --tmax 200 --nao-mostrar

Para cada satelite gera tres imagens:
    NOME_osc.png         elementos osculadores (azul)
    NOME_geo.png         elementos geometricos (laranja)
    NOME_comparacao.png  os dois sobrepostos
"""
import argparse
import glob
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

COR_OSC = "#2a78d6"   # azul
COR_GEO = "#eb6834"   # laranja

ROTULOS = [
    ("Semieixo $a$ (km)", False),
    ("Excentricidade $e$", False),
    ("Inclinação $i$ (°)", False),
    ("Pericentro $\\omega$ (°)", True),
    ("Nodo $\\Omega$ (°)", True),
    ("Anomalia média $M$ (°)", True),
]


def carrega(arquivo, tmax):
    if not os.path.exists(arquivo):
        return None
    d = np.loadtxt(arquivo, comments="#", ndmin=2)
    if tmax is not None:
        d = d[d[:, 0] <= tmax]
    return d


def figura(nome, series, titulo, arquivo, mostrar):
    """series: lista de (dados, cor, legenda)."""
    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    # coluna 1: a, e, i ; coluna 2: w, Om, M
    posicoes = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]

    for k, ((rotulo, angulo), (lin, col)) in enumerate(zip(ROTULOS, posicoes)):
        ax = axes[lin, col]
        for dados, cor, legenda in series:
            t, y = dados[:, 0], dados[:, k + 1]
            if angulo:
                # pontos, para nao ligar 359° com 0° por uma linha vertical
                ax.plot(t, np.mod(y, 360.0), ".", ms=1.5, color=cor, label=legenda)
            else:
                ax.plot(t, y, lw=0.8, color=cor, label=legenda)
        ax.set_ylabel(rotulo)
        if angulo:
            ax.set_ylim(0, 360)
            ax.set_yticks([0, 90, 180, 270, 360])
        else:
            ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.margins(x=0)

    axes[2, 0].set_xlabel("Tempo (dias)")
    axes[2, 1].set_xlabel("Tempo (dias)")
    if len(series) > 1:
        axes[0, 0].legend(loc="best", markerscale=6, frameon=False)
    fig.suptitle(f"{titulo} — {nome}", fontsize=14, fontweight="bold")
    fig.tight_layout()
    fig.savefig(arquivo, dpi=200, bbox_inches="tight")
    print(f"Salvo: {arquivo}")
    if mostrar:
        plt.show()
    plt.close(fig)


def plota(nome, pasta, tmax, mostrar):
    """Gera ate' 3 imagens: NOME_osc.png, NOME_geo.png e NOME_comparacao.png."""
    osc = carrega(os.path.join(pasta, f"{nome}_osc.dat"), tmax)
    geo = carrega(os.path.join(pasta, f"{nome}_geo.dat"), tmax)
    if osc is None and geo is None:
        print(f"Nada encontrado para {nome} em {pasta}")
        return

    base = os.path.join(pasta, nome)
    if osc is not None:
        figura(nome, [(osc, COR_OSC, "osculador")],
               "Elementos orbitais osculadores", f"{base}_osc.png", mostrar)
    if geo is not None:
        figura(nome, [(geo, COR_GEO, "geométrico")],
               "Elementos orbitais geométricos", f"{base}_geo.png", mostrar)
    if osc is not None and geo is not None:
        figura(nome, [(osc, COR_OSC, "osculador"), (geo, COR_GEO, "geométrico")],
               "Osculadores × geométricos", f"{base}_comparacao.png", mostrar)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("satelites", nargs="*", help="nomes dos satelites (padrao: todos)")
    p.add_argument("--pasta", default=".", help="pasta com os .dat (padrao: .)")
    p.add_argument("--tmax", type=float, default=None, help="tempo maximo no grafico (dias)")
    p.add_argument("--nao-mostrar", action="store_true", help="so' salva o PNG, sem abrir janela")
    a = p.parse_args()

    nomes = a.satelites
    if not nomes:
        arquivos = glob.glob(os.path.join(a.pasta, "*_osc.dat")) + glob.glob(os.path.join(a.pasta, "*_geo.dat"))
        nomes = sorted({os.path.basename(f).rsplit("_", 1)[0] for f in arquivos})
    for nome in nomes:
        plota(nome, a.pasta, a.tmax, not a.nao_mostrar)


if __name__ == "__main__":
    main()
