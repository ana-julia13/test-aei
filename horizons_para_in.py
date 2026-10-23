"""Converte a saida "Vector Table" do Horizons (estado cartesiano) em uma
linha do bloco "satelites:" do .in do nbody (elementos osculadores).

Serve para qualquer corpo, e principalmente para o Sol, que o Horizons
nao deixa sair como "Osculating Orbital Elements" em relacao a um planeta.

No Horizons (https://ssd.jpl.nasa.gov/horizons/app.html):
    1  Ephemeris Type : Vector Table
    2  Target Body    : Sun  (ou o corpo que quiser)
    3  Center         : o centro do planeta, ex. 500@899 (Netuno), 500@599 (Jupiter)
    4  Time           : a MESMA epoca das condicoes iniciais dos satelites (TDB)
    5  Table Settings : Reference plane = o mesmo usado nos satelites
                        (para o nbody: "Body mean equator and node of date")
                        Output units = km e km/s (au e dia tambem funciona)
Depois copie todo o texto da resposta para um arquivo, ex. sol.txt.

A conversao e' a mesma do cart_para_osc do nbody.f90 (spice_para_in.py),
com mu = GM_planeta + GM_corpo; o GM do planeta vem do planetas.dat.

Uso:
    python3 horizons_para_in.py sol.txt --planeta Netuno
    python3 horizons_para_in.py sol.txt --planeta Netuno --unidade kg
    python3 horizons_para_in.py lua.txt --planeta Netuno --massa 2.14e22 --nome Triton
"""
import argparse
import re

import numpy as np

from spice_para_in import cart_para_osc, osc_para_cart

UA_KM = 149597870.7
DIA = 86400.0
GM_SOL = 132712440041.279419        # km^3/s^2 (DE440)
G_KM = 6.67430e-20                  # km^3 kg^-1 s^-2 (o mesmo do nbody.f90)


def gm_do_planeta(nome, arquivo="planetas.dat"):
    """GM (km^3/s^2) do planeta em planetas.dat, como o nbody le."""
    alvo = nome.lower()
    for linha in open(arquivo, encoding="utf-8"):
        linha = linha.split("!")[0].split("#")[0].split()
        if len(linha) >= 2 and alvo in [n.lower() for n in linha[0].split("/")]:
            return float(linha[1])
    raise SystemExit(f"planeta '{nome}' nao encontrado em {arquivo}")


def le_horizons(arquivo):
    """-> (nome do alvo, r em km, v em km/s, linhas do cabecalho que importam)."""
    texto = open(arquivo, encoding="utf-8", errors="replace").read()
    if "$$SOE" not in texto:
        raise SystemExit(f"{arquivo}: nao achei $$SOE; copie a resposta inteira do Horizons")
    cabecalho, resto = texto.split("$$SOE", 1)
    registro = resto.split("$$EOE")[0]

    valores = {}
    for chave in ("X", "Y", "Z", "VX", "VY", "VZ"):
        m = re.search(rf"(?<![A-Z]){chave}\s*=\s*([-+]?[\d.]+(?:[Ee][-+]?\d+)?)", registro)
        if m is None:
            raise SystemExit(f"{arquivo}: nao achei '{chave} =' no primeiro registro; "
                             "use Vector Table com posicao e velocidade")
        valores[chave] = float(m.group(1))
    r = np.array([valores["X"], valores["Y"], valores["Z"]])
    v = np.array([valores["VX"], valores["VY"], valores["VZ"]])

    if re.search(r"Output units\s*:\s*AU-D", cabecalho):
        r, v = r*UA_KM, v*UA_KM/DIA

    m = re.search(r"Target body name:\s*(\S+)", cabecalho)
    nome = m.group(1) if m else "Corpo"
    info = [l.strip() for l in cabecalho.splitlines()
            if re.search(r"Center body name|Reference frame|Ref\. plane|Reference plane|Output units", l)]
    data = registro.strip().splitlines()[0].strip() if registro.strip() else "?"
    return nome, r, v, info, data


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("arquivo", help="texto da Vector Table do Horizons")
    p.add_argument("--planeta", required=True, help="planeta central, como em planetas.dat (ex. Netuno)")
    p.add_argument("--unidade", choices=["gm", "kg", "razao"], default="gm",
                   help="a mesma unidade_massa do seu .in (padrao: gm)")
    p.add_argument("--massa", type=float, default=None,
                   help="massa do corpo na --unidade escolhida (para o Sol nao precisa)")
    p.add_argument("--nome", default=None, help="nome no .in (padrao: o do Horizons)")
    p.add_argument("--planetas", default="planetas.dat", help="banco de planetas do nbody")
    a = p.parse_args()

    nome, r, v, info, data = le_horizons(a.arquivo)
    nome = a.nome or ("Sol" if nome.lower() == "sun" else nome)
    gm_p = gm_do_planeta(a.planeta, a.planetas)

    if a.massa is not None:
        massa = a.massa
        gm = {"gm": massa, "kg": G_KM*massa, "razao": massa*gm_p}[a.unidade]
    elif nome.lower() in ("sun", "sol"):
        gm = GM_SOL
        massa = {"gm": GM_SOL, "kg": GM_SOL/G_KM, "razao": GM_SOL/gm_p}[a.unidade]
    else:
        raise SystemExit(f"de' a massa de {nome} com --massa (na unidade '{a.unidade}')")

    mu = gm_p + gm
    el = cart_para_osc(mu, r, v)
    r2, v2 = osc_para_cart(mu, *el)
    erro = max(np.abs(r2 - r).max(), DIA*np.abs(v2 - v).max())

    print(f"Horizons: {data}")
    for l in info:
        print(f"          {l}")
    print(f"Planeta: {a.planeta}  GM = {gm_p} km^3/s^2   mu = {mu:.10g} km^3/s^2")
    print(f"a = {el[0]:.6g} km ({el[0]/UA_KM:.4f} UA)   e = {el[1]:.6f}   i = {el[2]:.6f}°   "
          f"(ida e volta: erro {erro:.1e} km)")
    print(f"\nLinha para o bloco satelites: (unidade_massa = {a.unidade}); "
          "ordem: nome massa a e i w Om M\n")
    print(f"{nome.capitalize():10s} {massa:22.15e} {el[0]:22.15f} {el[1]:22.15e} "
          f"{el[2]:20.15f} {el[3]:20.15f} {el[4]:20.15f} {el[5]:20.15f}")
    if nome.lower() in ("sun", "sol"):
        print("\nLembre: 'salvar = <so' os satelites>' no .in para o Sol nao gerar .dat.")


if __name__ == "__main__":
    main()
