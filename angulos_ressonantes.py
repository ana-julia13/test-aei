"""Plota os angulos criticos de uma ressonancia de movimento medio p:q
entre dois satelites, usando os arquivos gerados pelo nbody
(NOME_geo.dat e/ou NOME_osc.dat).

Os argumentos sao os da parte direta da expansao da funcao perturbadora
de Murray & Dermott (tabelas B.4, B.8, B.12 e B.16), para ressonancias de
ordem 1 a 4:

    phi = j*lambda' + (k - j)*lambda + (combinacao de varpi, varpi', Omega, Omega')

com  j = p,  k = p - q  (ordem).  Sem linha = satelite INTERNO,
com linha (') = satelite EXTERNO.

Uso:
    python3 angulos_ressonantes.py 73:69 Thebe Io
    python3 angulos_ressonantes.py 2:1 Metis Io --tipo osc --tmax 2000
    python3 angulos_ressonantes.py 73:69 Thebe Io --ids 4D4.1 4D4.9
    python3 angulos_ressonantes.py                 (pergunta a ressonancia e os satelites)

    - a ordem dos satelites nao importa: o interno e' escolhido pelo semieixo;
    - --tipo geo | osc | ambos  (padrao: ambos, se os arquivos existirem);
    - --ids escolhe so' alguns angulos (padrao: todos da ordem).

Salva ress_P_Q_INTERNO_EXTERNO_TIPO.png (nao abre janela).
"""
import argparse
import math
import os
import re

import matplotlib
matplotlib.use("Agg")                # so' salva a imagem, sem abrir janela
import matplotlib.pyplot as plt
import numpy as np

COR = "#2a78d6"

# ---------------------------------------------------------------------------
# Tabelas de Murray & Dermott (parte direta).
# Cada entrada: ID -> (c_varpi, c_varpi', c_Omega, c_Omega')
# ou seja, os coeficientes que somam ao j*lambda' + (k-j)*lambda.
#   varpi  = pericentro do INTERNO     varpi' = pericentro do EXTERNO
#   Omega  = nodo do INTERNO           Omega' = nodo do EXTERNO
# ---------------------------------------------------------------------------
TABELAS = {
    # Tabela B.4 - primeira ordem
    1: {
        "4D1.1":  (-1,  0,  0,  0),
        "4D1.2":  ( 0, -1,  0,  0),
        "4D1.3":  (-2,  1,  0,  0),
        "4D1.4":  ( 1, -2,  0,  0),
        "4D1.5":  ( 1,  0, -2,  0),
        "4D1.6":  ( 0,  1, -2,  0),
        "4D1.7":  (-1,  0,  1, -1),
        "4D1.8":  (-1,  0, -1,  1),
        "4D1.9":  ( 1,  0, -1, -1),
        "4D1.10": ( 0, -1,  1, -1),
        "4D1.11": ( 0, -1, -1,  1),
        "4D1.12": ( 0,  1, -1, -1),
        "4D1.13": ( 1,  0,  0, -2),
        "4D1.14": ( 0,  1,  0, -2),
    },
    # Tabela B.8 - segunda ordem
    2: {
        "4D2.1":  (-2,  0,  0,  0),
        "4D2.2":  (-1, -1,  0,  0),
        "4D2.3":  ( 0, -2,  0,  0),
        "4D2.4":  ( 0,  0, -2,  0),
        "4D2.5":  ( 0,  0, -1, -1),
        "4D2.6":  ( 0,  0,  0, -2),
        "4D2.7":  (-3,  1,  0,  0),
        "4D2.8":  ( 1, -3,  0,  0),
        "4D2.9":  ( 1, -1, -2,  0),
        "4D2.10": (-1,  1, -2,  0),
        "4D2.11": (-2,  0,  1, -1),
        "4D2.12": (-2,  0, -1,  1),
        "4D2.13": (-1, -1,  1, -1),
        "4D2.14": (-1, -1, -1,  1),
        "4D2.15": ( 1, -1, -1, -1),
        "4D2.16": (-1,  1, -1, -1),
        "4D2.17": ( 0, -2,  1, -1),
        "4D2.18": ( 0, -2, -1,  1),
        "4D2.19": ( 0,  0, -3,  1),
        "4D2.20": ( 1, -1,  0, -2),
        "4D2.21": (-1,  1,  0, -2),
        "4D2.22": ( 0,  0,  1, -3),
    },
    # Tabela B.12 - terceira ordem
    3: {
        "4D3.1":  (-3,  0,  0,  0),
        "4D3.2":  (-2, -1,  0,  0),
        "4D3.3":  (-1, -2,  0,  0),
        "4D3.4":  ( 0, -3,  0,  0),
        "4D3.5":  (-1,  0, -2,  0),
        "4D3.6":  ( 0, -1, -2,  0),
        "4D3.7":  (-1,  0, -1, -1),
        "4D3.8":  ( 0, -1, -1, -1),
        "4D3.9":  (-1,  0,  0, -2),
        "4D3.10": ( 0, -1,  0, -2),
    },
    # Tabela B.16 - quarta ordem
    4: {
        "4D4.1":  (-4,  0,  0,  0),
        "4D4.2":  (-3, -1,  0,  0),
        "4D4.3":  (-2, -2,  0,  0),
        "4D4.4":  (-1, -3,  0,  0),
        "4D4.5":  ( 0, -4,  0,  0),
        "4D4.6":  (-2,  0, -2,  0),
        "4D4.7":  (-1, -1, -2,  0),
        "4D4.8":  ( 0, -2, -2,  0),
        "4D4.9":  ( 0,  0, -4,  0),
        "4D4.10": (-2,  0, -1, -1),
        "4D4.11": (-1, -1, -1, -1),
        "4D4.12": ( 0, -2, -1, -1),
        "4D4.13": ( 0,  0, -3, -1),
        "4D4.14": (-2,  0,  0, -2),
        "4D4.15": (-1, -1,  0, -2),
        "4D4.16": ( 0, -2,  0, -2),
        "4D4.17": ( 0,  0, -2, -2),
        "4D4.18": ( 0,  0, -1, -3),
        "4D4.19": ( 0,  0,  0, -4),
    },
}

# Conferencia das tabelas (regras de d'Alembert): a soma de todos os
# coeficientes e' zero e a soma dos coeficientes dos nodos e' par.
for _k, _tab in TABELAS.items():
    for _id, (_w, _wl, _o, _ol) in _tab.items():
        assert _k + _w + _wl + _o + _ol == 0, _id
        assert (_o + _ol) % 2 == 0, _id

TITULO_TIPO = {"geo": "elementos geométricos", "osc": "elementos osculadores"}


# ---------------------------------------------------------------------------
def le_ressonancia(texto):
    m = re.fullmatch(r"\s*(\d+)\s*[:/]\s*(\d+)\s*", texto)
    if not m:
        raise SystemExit(f"Ressonancia invalida: '{texto}'. Use o formato p:q, ex.: 73:69")
    p, q = int(m.group(1)), int(m.group(2))
    if p < q:
        p, q = q, p
    k = p - q
    if not 1 <= k <= 4:
        raise SystemExit(f"{p}:{q} e' de ordem {k}; so' tenho as tabelas de ordem 1 a 4.")
    if math.gcd(p, q) > 1:
        print(f"Obs.: {p}:{q} nao esta' na forma reduzida "
              f"({p // math.gcd(p, q)}:{q // math.gcd(p, q)}); uso j = {p}, ordem {k} como pedido.")
    return p, q, k


def carrega(pasta, nome, tipo):
    arq = os.path.join(pasta, f"{nome}_{tipo}.dat")
    if not os.path.exists(arq):
        return None
    return np.loadtxt(arq, comments="#", ndmin=2)


def angulos(d_in, d_ex, p, q, k, ids):
    """Calcula os angulos criticos (graus, em [0, 360))."""
    n = min(len(d_in), len(d_ex))
    d_in, d_ex = d_in[:n], d_ex[:n]
    if not np.allclose(d_in[:, 0], d_ex[:, 0]):
        print("AVISO: os tempos dos dois arquivos nao coincidem (sao da mesma simulacao?)")
    t = d_in[:, 0]

    # colunas: 0 t | 1 a | 2 e | 3 i | 4 w | 5 Omega | 6 M
    w, Om, M = d_in[:, 4], d_in[:, 5], d_in[:, 6]
    wl, Oml, Ml = d_ex[:, 4], d_ex[:, 5], d_ex[:, 6]
    varpi, varpil = w + Om, wl + Oml
    lam, laml = M + varpi, Ml + varpil

    base = p * laml - q * lam              # = j*lambda' + (k - j)*lambda
    res = {}
    for id_ in ids:
        cw, cwl, co, col = TABELAS[k][id_]
        phi = base + cw * varpi + cwl * varpil + co * Om + col * Oml
        res[id_] = np.mod(phi, 360.0)
    return t, res


def argumento_texto(id_, k):
    """Texto do argumento, ex.: 73λ' − 69λ − 4ϖ (so' para o titulo do painel)."""
    nomes = ["ϖ", "ϖ'", "Ω", "Ω'"]
    s = ""
    for c, n in zip(TABELAS[k][id_], nomes):
        if c == 0:
            continue
        s += (" + " if c > 0 else " − ") + (f"{abs(c)}" if abs(c) != 1 else "") + n
    return s


def figura(t, res, p, q, k, interno, externo, tipo, tmax, arquivo):
    n = len(res)
    ncol = 4 if n > 3 else n
    nlin = math.ceil(n / ncol)
    fig, axes = plt.subplots(nlin, ncol, figsize=(4.5 * ncol, 2.6 * nlin + 1),
                             sharex=True, sharey=True, squeeze=False)
    axes = axes.flatten()

    for ax, (id_, phi) in zip(axes, res.items()):
        ax.scatter(t, phi, s=0.8, c=COR, alpha=0.7, linewidths=0)
        ax.set_title(f"{id_}:  …{argumento_texto(id_, k)}", fontsize=9)
        ax.set_ylim(0, 360)
        ax.set_yticks([0, 90, 180, 270, 360])
        ax.set_xlim(t.min(), tmax if tmax is not None else t.max())
        ax.grid(True, linestyle="--", alpha=0.5)
    for ax in axes[n:]:
        fig.delaxes(ax)

    fig.supxlabel("Tempo (dias)", fontweight="bold")
    fig.supylabel("Ângulo ressonante (°)", fontweight="bold")
    fig.suptitle(f"Ressonância {p}:{q} (ordem {k}) — {interno} (interno) × {externo} (externo) — "
                 f"{TITULO_TIPO[tipo]}\n"
                 f"φ = {p}λ' − {q}λ + …", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(arquivo, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Salvo: {arquivo}")


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("ressonancia", nargs="?", help="ressonancia p:q, ex.: 73:69")
    ap.add_argument("satelites", nargs="*", help="nomes dos dois satelites (como nos .dat)")
    ap.add_argument("--tipo", default="ambos", choices=["geo", "osc", "ambos"],
                    help="elementos usados (padrao: ambos)")
    ap.add_argument("--ids", nargs="+", help="so' esses angulos, ex.: 4D4.1 4D4.9")
    ap.add_argument("--pasta", default=".", help="pasta com os .dat (padrao: .)")
    ap.add_argument("--tmax", type=float, default=None, help="tempo maximo no grafico (dias)")
    a = ap.parse_args()

    # pergunta o que faltar
    if a.ressonancia is None:
        a.ressonancia = input("Ressonancia (p:q, ex. 73:69): ")
    if len(a.satelites) != 2:
        a.satelites = input("Os dois satelites (ex. Thebe Io): ").split()
        if len(a.satelites) != 2:
            raise SystemExit("Preciso de exatamente dois satelites.")

    p, q, k = le_ressonancia(a.ressonancia)

    ids = list(TABELAS[k])
    if a.ids:
        ids = [i.upper() for i in a.ids]
        ruins = [i for i in ids if i not in TABELAS[k]]
        if ruins:
            raise SystemExit(f"{', '.join(ruins)} nao existe(m) na ordem {k}. "
                             f"Disponiveis: {', '.join(TABELAS[k])}")

    tipos = ["geo", "osc"] if a.tipo == "ambos" else [a.tipo]
    feito = False
    for tipo in tipos:
        s1, s2 = a.satelites
        d1, d2 = carrega(a.pasta, s1, tipo), carrega(a.pasta, s2, tipo)
        if d1 is None or d2 is None:
            falta = [f"{s}_{tipo}.dat" for s, d in ((s1, d1), (s2, d2)) if d is None]
            print(f"Pulando {tipo}: nao achei {', '.join(falta)} em {a.pasta}")
            continue

        # interno = menor semieixo medio
        if np.mean(d1[:, 1]) > np.mean(d2[:, 1]):
            s1, s2, d1, d2 = s2, s1, d2, d1
        if a.tmax is not None:
            d1, d2 = d1[d1[:, 0] <= a.tmax], d2[d2[:, 0] <= a.tmax]

        t, res = angulos(d1, d2, p, q, k, ids)
        arq = os.path.join(a.pasta, f"ress_{p}_{q}_{s1}_{s2}_{tipo}.png")
        figura(t, res, p, q, k, s1, s2, tipo, a.tmax, arq)
        feito = True

    if not feito:
        raise SystemExit("Nenhuma figura gerada.")


if __name__ == "__main__":
    main()
