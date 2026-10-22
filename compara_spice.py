"""Compara os elementos osculadores do nbody com os do SPICE (efemerides).

Le, para cada satelite:
    nbody : NOME_osc.dat              -> t(dias) a(km) e i w om M  (graus)
    SPICE : elementos_NOME.csv        -> data_utc, et_tdb_s, a, e, i, Omega, omega, M
            (gerado pelo elementos_spice.py)

e faz duas figuras:
    NOME_spice.png      os 6 elementos dos dois metodos sobrepostos
                        (nbody em azul, SPICE em vermelho)
    NOME_spice_dif.png  diferenca nbody - SPICE de cada elemento
                        (angulos levados para [-180, 180) graus)

O t = 0 do nbody e' a primeira linha do CSV do SPICE (a epoca INICIO do
elementos_spice.py, a mesma do estados_iniciais.csv). Se a condicao inicial
do nbody veio de outra epoca, use --epoca-et.

Para a comparacao fazer sentido, os dois precisam estar no mesmo referencial
(FRAME = "EQUADOR_JUPITER" no elementos_spice.py, que e' o plano do equador
usado pelo nbody) e com o mesmo mu (SOMAR_GM_LUA = True <-> G(M+m) do nbody).
Lembre que w e Om ficam mal definidos quando e ou i sao muito pequenos.

Uso:
    python3 compara_spice.py                       # todos os NOME_osc.dat com CSV do SPICE
    python3 compara_spice.py Thebe Io --tmax 100
    python3 compara_spice.py Io --pasta-nbody . --pasta-spice saida_elementos
"""
import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")                # so' salva a imagem, sem abrir janela
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

DIA = 86400.0

# cores de cada metodo (mude aqui se quiser)
CORES = {"nbody": "#1f4fd6", "SPICE": "#d62728"}

# transparencia das curvas (1 = opaco). Onde as curvas se sobrepoem, as cores
# se misturam (azul + vermelho fica roxo). Tambem da' para mudar com --alpha.
ALPHA = 0.5

ROTULOS = [
    ("Semieixo $a$ (km)", False),
    ("Excentricidade $e$", False),
    ("Inclinação $i$ (°)", False),
    ("Pericentro $\\omega$ (°)", True),
    ("Nodo $\\Omega$ (°)", True),
    ("Anomalia média $M$ (°)", True),
]
POSICOES = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]
ROTULOS_DIF = ["$\\Delta a$ (km)", "$\\Delta e$", "$\\Delta i$ (°)",
               "$\\Delta\\omega$ (°)", "$\\Delta\\Omega$ (°)", "$\\Delta M$ (°)"]

# nome no .in do nbody (minusculo) -> nome SPICE (minusculo, como no CSV)
APELIDOS = {
    "amalteia": "amalthea",
    "tebe": "thebe",
    "adrasteia": "adrastea",
    "ganimedes": "ganymede",
    "calisto": "callisto",
}


# ------------------------------------------------------------------ leitura
def le_nbody(arq):
    return np.loadtxt(arq, comments="#", ndmin=2)


def le_spice(arq):
    """CSV do elementos_spice.py -> (et, colunas a e i Omega omega M)."""
    d = np.genfromtxt(arq, delimiter=",", skip_header=1, usecols=range(1, 8), ndmin=2)
    return d[:, 0], d[:, 1:]


def arquivo_spice(nome, pasta):
    alvo = APELIDOS.get(nome.lower(), nome.lower())
    arq = os.path.join(pasta, f"elementos_{alvo}.csv")
    return arq if os.path.exists(arq) else None


def spice_no_formato_nbody(et, el, epoca):
    """-> t(dias) a e i w om M, na mesma ordem das colunas do nbody."""
    t = (et - epoca)/DIA
    a, e, i, om, w, m = el.T
    return np.column_stack([t, a, e, i, w, om, m])


# ------------------------------------------------------------- diferencas
def diferencas(nb, sp):
    """nbody - SPICE nos tempos do nbody (SPICE interpolado). Os angulos sao
    desenrolados antes de interpolar; o passo do SPICE precisa ser bem menor
    que o periodo orbital (para M)."""
    dentro = (nb[:, 0] >= sp[0, 0]) & (nb[:, 0] <= sp[-1, 0])
    nb = nb[dentro]
    dif = np.empty_like(nb)
    dif[:, 0] = nb[:, 0]
    for k, (_, angulo) in enumerate(ROTULOS):
        y = sp[:, k + 1]
        if angulo:
            y = np.unwrap(y, period=360.0)
            d = nb[:, k + 1] - np.interp(nb[:, 0], sp[:, 0], y)
            dif[:, k + 1] = np.mod(d + 180.0, 360.0) - 180.0
        else:
            dif[:, k + 1] = nb[:, k + 1] - np.interp(nb[:, 0], sp[:, 0], y)
    return dif


# ------------------------------------------------------------------ figuras
def eixos_base(titulo):
    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=True)
    fig.suptitle(titulo, fontsize=14, fontweight="bold")
    for ax in axes.flat:
        ax.grid(True, linestyle="--", alpha=0.4)
        ax.margins(x=0)
    axes[2, 0].set_xlabel("Tempo (dias)")
    axes[2, 1].set_xlabel("Tempo (dias)")
    return fig, axes


def salva(fig, arquivo):
    fig.tight_layout()
    fig.savefig(arquivo, dpi=200, bbox_inches="tight")
    print(f"Salvo: {arquivo}")
    plt.close(fig)


def figura_elementos(nome, series, arquivo, alpha=ALPHA):
    fig, axes = eixos_base(f"Elementos osculadores — {nome} (nbody × SPICE)")
    for k, ((rotulo, angulo), (lin, col)) in enumerate(zip(ROTULOS, POSICOES)):
        ax = axes[lin, col]
        for metodo, d in series.items():
            if angulo:
                # pontos, para nao ligar 359° com 0° por uma linha vertical
                ax.plot(d[:, 0], np.mod(d[:, k + 1], 360.0), ".", ms=1.5,
                        color=CORES[metodo], alpha=alpha, label=metodo)
            else:
                ax.plot(d[:, 0], d[:, k + 1], lw=0.8, color=CORES[metodo], alpha=alpha, label=metodo)
        ax.set_ylabel(rotulo)
        if angulo:
            ax.set_ylim(0, 360)
            ax.set_yticks([0, 90, 180, 270, 360])
        else:
            ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    leg = axes[0, 0].legend(loc="best", markerscale=6, frameon=False)
    for h in getattr(leg, "legend_handles", None) or leg.legendHandles:
        h.set_alpha(1.0)                 # legenda sem transparencia
    salva(fig, arquivo)


def figura_diferencas(nome, dif, arquivo):
    fig, axes = eixos_base(f"Diferença nbody − SPICE — {nome}")
    for k, (rotulo, (lin, col)) in enumerate(zip(ROTULOS_DIF, POSICOES)):
        ax = axes[lin, col]
        ax.plot(dif[:, 0], dif[:, k + 1], lw=0.8, color="k")
        ax.axhline(0.0, color="0.5", lw=0.6)
        ax.set_ylabel(rotulo)
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    salva(fig, arquivo)


# -------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("nomes", nargs="*", help="satelites (padrao: todos os NOME_osc.dat da pasta do nbody)")
    p.add_argument("--pasta-nbody", default=".", help="pasta dos NOME_osc.dat (padrao: .)")
    p.add_argument("--pasta-spice", default="saida_elementos",
                   help="pasta dos elementos_NOME.csv (padrao: saida_elementos)")
    p.add_argument("--epoca-et", type=float, default=None,
                   help="et TDB (s) do t = 0 do nbody (padrao: primeira linha do CSV)")
    p.add_argument("--tmax", type=float, default=None, help="tempo maximo no grafico (dias)")
    p.add_argument("--alpha", type=float, default=ALPHA, help="transparencia das curvas, 0 a 1 (padrao: 0.5)")
    a = p.parse_args()

    nomes = a.nomes or sorted(os.path.basename(f)[:-len("_osc.dat")]
                              for f in glob.glob(os.path.join(a.pasta_nbody, "*_osc.dat")))
    if not nomes:
        p.error(f"nenhum *_osc.dat em {a.pasta_nbody}")

    for nome in nomes:
        arq_nb = os.path.join(a.pasta_nbody, f"{nome}_osc.dat")
        arq_sp = arquivo_spice(nome, a.pasta_spice)
        if not os.path.exists(arq_nb):
            print(f"{nome}: sem {arq_nb}, pulando")
            continue
        if arq_sp is None:
            print(f"{nome}: sem CSV do SPICE em {a.pasta_spice}, pulando")
            continue

        et, el = le_spice(arq_sp)
        epoca = et[0] if a.epoca_et is None else a.epoca_et
        series = {"nbody": le_nbody(arq_nb), "SPICE": spice_no_formato_nbody(et, el, epoca)}
        if a.tmax is not None:
            series = {m: d[(d[:, 0] >= 0.0) & (d[:, 0] <= a.tmax)] for m, d in series.items()}

        figura_elementos(nome, series, f"{nome}_spice.png", a.alpha)

        dif = diferencas(series["nbody"], series["SPICE"])
        if len(dif) == 0:
            print(f"{nome}: os tempos do nbody e do SPICE nao se sobrepoem (confira --epoca-et)")
            continue
        figura_diferencas(nome, dif, f"{nome}_spice_dif.png")
        print(f"   max |Δa| = {np.abs(dif[:, 1]).max():.4g} km   "
              f"max |Δe| = {np.abs(dif[:, 2]).max():.3g}   max |Δi| = {np.abs(dif[:, 3]).max():.3g}°")


if __name__ == "__main__":
    main()
