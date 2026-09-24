"""Compara os elementos orbitais osculadores de um satelite calculados por
tres codigos: NLSN (codigo antigo), nbody (este) e Mercury (.aei).

Formatos lidos:
    NLSN    : dois arquivos, "t a(km) e i(deg)" e "t w om M" (graus)
    nbody   : NOME_osc.dat  -> t a(km) e i w om M
    Mercury : NOME.aei      -> 4 linhas de cabecalho, depois
              t a(UA) e i peri node M   (peri = argumento do pericentro)

Uso (todos os arquivos sao opcionais; use os que tiver):
    python3 compara_codigos.py --nome Thebe \
        --nlsn saida3.dat saida4.dat \
        --nbody Thebe_osc.dat \
        --mercury Thebe.aei \
        --tmax 500 --alpha 0.4

Gera:
    NOME_codigos.png      os elementos dos tres codigos sobrepostos
    NOME_diferencas.png   diferencas de cada codigo em relacao a' referencia
                          (--ref nbody | nlsn | mercury; padrao: nbody)
e imprime na tela a maior diferenca de cada elemento.
"""
import argparse

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

UA_KM = 1.495978707e8          # UA em km (IAU 2012)

# cores de cada codigo (mude aqui se quiser)
CORES = {"nbody": "#1f4fd6", "NLSN": "#d62728", "Mercury": "#1a9e3a"}

# transparencia das curvas: 1 = opaco, 0 = invisivel.
# Com ~0.5, onde as curvas se sobrepoem as cores se misturam
# (azul + vermelho fica roxo). Tambem da' para mudar com --alpha.
ALPHA = 0.5

ROTULOS = [
    ("Semieixo $a$ (km)", "a (km)", False),
    ("Excentricidade $e$", "e", False),
    ("Inclinação $i$ (°)", "i (°)", False),
    ("Pericentro $\\omega$ (°)", "ω (°)", True),
    ("Nodo $\\Omega$ (°)", "Ω (°)", True),
    ("Anomalia média $M$ (°)", "M (°)", True),
]
POSICOES = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)]


# ------------------------------------------------------------------ leitura
def le_nlsn(arq_aei, arq_ang):
    """saida3.dat (t a e i) + saida4.dat (t w om M) -> t a e i w om M."""
    d1 = np.loadtxt(arq_aei, ndmin=2)
    d2 = np.loadtxt(arq_ang, ndmin=2)
    n = min(len(d1), len(d2))
    return np.column_stack([d1[:n, 0:4], d2[:n, 1:4]])


def le_nbody(arq):
    return np.loadtxt(arq, comments="#", ndmin=2)


def le_mercury(arq):
    """Le o .aei do Mercury (pula o cabecalho ate' a linha com 'Time')."""
    with open(arq) as f:
        linhas = f.readlines()
    inicio, em_anos = 0, False
    for k, lin in enumerate(linhas):
        if "Time" in lin:
            inicio = k + 1
            em_anos = "year" in lin.lower()
            break
    # genfromtxt aguenta campos estourados ("*****") como NaN
    d = np.genfromtxt(linhas[inicio:], invalid_raise=False)
    d = d[:, :7]
    d[:, 1] *= UA_KM                 # a: UA -> km
    if em_anos:
        d[:, 0] *= 365.25            # tempo em dias
    return d


# ------------------------------------------------------------- comparacao
def diferenca(dados, ref):
    """Diferenca dados - ref nos tempos de 'dados'.

    Usa os instantes em comum quando existem; senao interpola a referencia
    (angulos sao desenrolados antes de interpolar)."""
    t, tr = dados[:, 0], ref[:, 0]
    tol = 1e-6 * max(1.0, np.abs(tr).max())
    idx = np.searchsorted(tr, t)
    idx = np.clip(idx, 1, len(tr) - 1)
    viz = np.where(np.abs(tr[idx - 1] - t) < np.abs(tr[idx] - t), idx - 1, idx)
    comum = np.abs(tr[viz] - t) < tol

    if comum.sum() >= 10:
        tt, a, b = t[comum], dados[comum], ref[viz[comum]]
        interpolado = False
    else:
        dentro = (t >= tr[0]) & (t <= tr[-1])
        tt, a = t[dentro], dados[dentro]
        b = np.empty_like(a)
        b[:, 0] = tt
        for k in range(1, 7):
            y = ref[:, k]
            if ROTULOS[k - 1][2]:
                y = np.degrees(np.unwrap(np.radians(y)))
            b[:, k] = np.interp(tt, tr, y)
        interpolado = True

    dif = a[:, 1:] - b[:, 1:]
    for k in range(6):
        if ROTULOS[k][2]:
            dif[:, k] = (dif[:, k] + 180.0) % 360.0 - 180.0
    return tt, dif, interpolado


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


def figura_elementos(nome, series, arquivo, mostrar, alpha=ALPHA):
    fig, axes = eixos_base(f"Elementos osculadores — {nome}")
    for k, ((rotulo, _, angulo), (lin, col)) in enumerate(zip(ROTULOS, POSICOES)):
        ax = axes[lin, col]
        for cod, d in series.items():
            if angulo:
                ax.plot(d[:, 0], np.mod(d[:, k + 1], 360.0), ".", ms=1.5,
                        color=CORES[cod], alpha=alpha, label=cod)
            else:
                ax.plot(d[:, 0], d[:, k + 1], lw=0.8, color=CORES[cod], alpha=alpha, label=cod)
        ax.set_ylabel(rotulo)
        if angulo:
            ax.set_ylim(0, 360)
            ax.set_yticks([0, 90, 180, 270, 360])
        else:
            ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    leg = axes[0, 0].legend(loc="best", markerscale=6, frameon=False)
    for h in getattr(leg, "legend_handles", None) or leg.legendHandles:
        h.set_alpha(1.0)                 # legenda sem transparencia
    fig.tight_layout()
    fig.savefig(arquivo, dpi=200, bbox_inches="tight")
    print(f"Salvo: {arquivo}")
    if mostrar:
        plt.show()
    plt.close(fig)


def figura_diferencas(nome, series, ref, arquivo, mostrar, alpha=ALPHA):
    outros = [c for c in series if c != ref]
    if not outros:
        return
    fig, axes = eixos_base(f"Diferenças em relação ao {ref} — {nome}")
    print(f"\nMaior |diferença| em relação ao {ref}:")
    print(f"{'':10s}" + "".join(f"{r[1]:>14s}" for r in ROTULOS))
    algum_interp = False
    for cod in outros:
        tt, dif, interp = diferenca(series[cod], series[ref])
        algum_interp = algum_interp or interp
        aviso = "  <- INTERPOLADO (passos de saida diferentes)" if interp else ""
        print(f"{cod:10s}" + "".join(f"{np.nanmax(np.abs(dif[:, k])):14.4e}" for k in range(6)) + aviso)
        for k, (lin, col) in enumerate(POSICOES):
            axes[lin, col].plot(tt, dif[:, k], lw=0.8, color=CORES[cod], alpha=alpha,
                                label=f"{cod} − {ref}")
    for k, (lin, col) in enumerate(POSICOES):
        axes[lin, col].set_ylabel("Δ" + ROTULOS[k][1])
        axes[lin, col].axhline(0.0, color="0.5", lw=0.6)
    leg = axes[0, 0].legend(loc="best", frameon=False)
    for h in getattr(leg, "legend_handles", None) or leg.legendHandles:
        h.set_alpha(1.0)
    fig.tight_layout()
    fig.savefig(arquivo, dpi=200, bbox_inches="tight")
    if algum_interp:
        print("AVISO: interpolar elementos que oscilam rapido cria diferencas falsas.\n"
              "       Para comparar de verdade, rode com o mesmo passo de saida nos dois codigos.")
    print(f"Salvo: {arquivo}")
    if mostrar:
        plt.show()
    plt.close(fig)


# -------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--nome", default="satelite", help="nome do satelite (titulo e arquivos)")
    p.add_argument("--nlsn", nargs=2, metavar=("AEI", "ANG"), help="saida3.dat saida4.dat do NLSN")
    p.add_argument("--nbody", metavar="OSC", help="NOME_osc.dat do nbody")
    p.add_argument("--mercury", metavar="AEI", help="NOME.aei do Mercury")
    p.add_argument("--ref", default="nbody", choices=["nbody", "nlsn", "mercury"],
                   help="codigo de referencia para as diferencas (padrao: nbody)")
    p.add_argument("--tmax", type=float, default=None, help="tempo maximo no grafico (dias)")
    p.add_argument("--alpha", type=float, default=ALPHA, help="transparencia das curvas, 0 a 1 (padrao: 0.5)")
    p.add_argument("--nao-mostrar", action="store_true", help="so' salva os PNG")
    a = p.parse_args()

    series = {}
    if a.nbody:
        series["nbody"] = le_nbody(a.nbody)
    if a.nlsn:
        series["NLSN"] = le_nlsn(*a.nlsn)
    if a.mercury:
        series["Mercury"] = le_mercury(a.mercury)
    if not series:
        p.error("de' pelo menos um de --nbody, --nlsn, --mercury")

    if a.tmax is not None:
        series = {c: d[d[:, 0] <= a.tmax] for c, d in series.items()}

    ref = {"nbody": "nbody", "nlsn": "NLSN", "mercury": "Mercury"}[a.ref]
    if ref not in series:
        ref = next(iter(series))

    mostrar = not a.nao_mostrar
    figura_elementos(a.nome, series, f"{a.nome}_codigos.png", mostrar, a.alpha)
    figura_diferencas(a.nome, series, ref, f"{a.nome}_diferencas.png", mostrar, a.alpha)


if __name__ == "__main__":
    main()
