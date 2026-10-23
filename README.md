# nbody — satélites em torno de um planeta achatado (Fortran)

Integrador de N corpos (RA15, Everhart 1985) para satélites em torno de um
planeta com J2 e J4. Gera, para cada satélite, a evolução temporal dos
**elementos osculadores** e dos **elementos geométricos** (Renner & Sicardy 2006).

O código não muda entre simulações: tudo vem do arquivo `.in` (estilo Mercury), e
os dados dos planetas (GM, raio, J2, J4) ficam em `planetas.dat`. Basta escrever
`planeta = Jupiter`.

## Compilar e rodar

```bash
make                         # ou: gfortran -O2 nbody.f90 -o nbody
./nbody                      # lê simulacao.in
./nbody netuno.in            # ou outro arquivo de entrada
python3 plot_elementos.py    # gráficos de todos os satélites (precisa de numpy e matplotlib)
```

## Arquivos

| arquivo | o que é |
|---|---|
| `nbody.f90` | o código (forças, RA15, conversões de elementos, leitura do `.in`) |
| `planetas.dat` | banco de planetas: `nome  GM(km³/s²)  R_eq(km)  J2  J4`. Para adicionar um planeta, acrescente uma linha. |
| `simulacao.in` | exemplo de entrada (sistema de Júpiter do `ips_mas_pra_aei_jup.f`) |
| `plot_elementos.py` | gera 3 imagens por satélite: `NOME_osc.png`, `NOME_geo.png` e `NOME_comparacao.png` |
| `compara_codigos.py` | compara os osculadores do nbody, do NLSN (`saida3.dat`+`saida4.dat`) e do Mercury (`.aei`) |
| `elementos_spice.py` | elementos osculadores das luas de Júpiter a partir do SPICE (+ `estados_iniciais.csv`); `POLO_DINAMICO` corrige o polo do frame |
| `compara_spice.py` | compara os osculadores do nbody (`NOME_osc.dat`) com os do SPICE (`elementos_NOME.csv` do `elementos_spice.py`) |
| `spice_para_in.py` | converte o `estados_iniciais.csv` (cartesiano, do `elementos_spice.py`) num `.in` com elementos osculadores |
| `horizons_para_in.py` | converte uma "Vector Table" do Horizons (ex. o Sol visto do planeta) na linha do `.in` com elementos osculadores |
| `polo_spice.py` | testa se o polo IAU do frame `EQUADOR_JUPITER` é o eixo em torno do qual as luas do SPICE precessionam (círculos do vetor inclinação) |
| `angulos_ressonantes.py` | ângulos críticos de uma ressonância p:q de ordem 1 a 4 (tabelas B.4, B.8, B.12, B.16 de Murray & Dermott) |

## O arquivo `.in`

```
planeta     = Jupiter        ! procurado em planetas.dat
tempo_final = 1000.0         ! dias
passo_saida = 0.05           ! dias entre linhas da saída
precisao    = 12             ! LL do RA15

saida_osculadores = sim
saida_geometricos = sim
pasta_saida       = .

unidade_massa = kg           ! kg | razao (m/M_planeta) | gm (km^3/s^2)
satelites:
! nome   massa     a(km)      e        i(deg)   w(deg)   Om(deg)  M(deg)
Thebe    4.517e17  222372.0   0.01827  1.0677   38.43    326.80   292.16
Io       8.927e22  351740.64  0.00378  0.0363   281.87   30.90    238.11
fim
```

* Os elementos iniciais são **osculadores, planetocêntricos, no plano do equador do planeta**.
* Massa `0` é permitida (partícula-teste).
* `salvar = Thebe Io` grava `.dat` só desses satélites; os outros continuam na
  integração (perturbando normalmente), só não geram arquivos.
* Dá para sobrescrever os valores do banco no próprio `.in` (`gm_planeta`,
  `raio_planeta`, `j2`, `j4`), ou usar `planeta = personalizado` e dar todos eles.

## Saída

Para cada satélite, `NOME_osc.dat` e `NOME_geo.dat`, com 7 colunas:

```
tempo(dias)  a(km)  e  inc(deg)  w(deg)  om(deg)  M(deg)
```

Esse é o mesmo formato do `sat1_geo.dat` do `osc_geo.f90`. Os elementos
osculadores usam μ = G(M + m); os geométricos usam G·M, J2 e J4 do planeta (J6 = 0),
como no `geo` original. Satélites retrógrados (i > 90°) só geram o `_osc.dat`,
porque o método de Renner & Sicardy não se aplica a eles.

## Comparar com NLSN e Mercury

```bash
python3 compara_codigos.py --nome Thebe --nbody Thebe_osc.dat \
    --nlsn saida3.dat saida4.dat --mercury Thebe.aei --tmax 500
```

Salva `Thebe_codigos.png` com os elementos dos três códigos sobrepostos (curvas
semitransparentes; mude com `--alpha`). Não abre janela.

## Ângulos ressonantes

```bash
python3 angulos_ressonantes.py 73:69 Thebe Io            # geo e osc
python3 angulos_ressonantes.py 2:1 Metis Io --tipo geo --tmax 2000
python3 angulos_ressonantes.py 73:69 Thebe Io --ids 4D4.1 4D4.9
```

A ordem (p − q) escolhe a tabela automaticamente; o satélite interno é o de menor
semieixo. Salva `ress_P_Q_INTERNO_EXTERNO_geo.png` / `_osc.png`.

## Física e unidades

* Equações planetocêntricas (Cowell), com termos indiretos, interação mútua
  entre os satélites, J2/J4 e a reação do planeta ao próprio achatamento.
  Isso é equivalente às variáveis canônicas do código original.
* Unidades internas iguais às do original: raio do planeta, dia e G·m.

## Validação

* Comparado com `ips_mas_pra_aei_jup.f` (mesmas condições iniciais e constantes,
  120 dias): diferenças no nível do arredondamento da saída do original.
* A rotina `geo` dá resultados idênticos ao `osc_geo.f90`.
* A aceleração de J2/J4 confere com o gradiente numérico do potencial.
* A energia se conserva em ~1e-14 (1000 dias, 5 satélites).

## Comparar com o SPICE

Para o nbody partir exatamente do estado do SPICE, gere o `.in` a partir do
`estados_iniciais.csv` (frame `EQUADOR_JUPITER`):

```bash
python3 spice_para_in.py saida_elementos/estados_iniciais.csv -o jupiter_spice.in --tempo-final 365 --passo 0.01
./nbody jupiter_spice.in
```

A conversão estado → elementos é a mesma do `cart_para_osc` do nbody, com
μ = GM_Júpiter + GM_lua e o GM de Júpiter do próprio SPICE (`gm_planeta` no `.in`),
então em t = 0 os dois coincidem.

O nbody só tem as forças de Júpiter (J2, J4) e dos corpos do `.in`; a efeméride
do SPICE tem também o Sol, cuja maré mexe na inclinação de Calisto (período de
~2170 dias, metade do ano de Júpiter). O `spice_para_in.py` inclui o Sol por padrão: edite `KERNELS` no topo dele
(o mesmo caminho do `elementos_spice.py`) e rode `python3 spice_para_in.py`
(ou `--sem-sol` para só as luas). O script pega o estado do Sol relativo a Júpiter na época e no frame do CSV e o põe
na lista de satélites. Visto de Júpiter (o referencial do nbody), o Sol percorre uma
órbita de a ≈ 5,2 UA e período de 11,86 anos; o termo indireto do nbody faz cada lua
sentir só a maré solar. A linha `salvar` deixa o Sol fora dos `.dat`. Precisa de um
SPK planetário (de440s.bsp ou similar) além do `jup365.bsp`.

```bash
python3 compara_spice.py                          # todos os NOME_osc.dat que tem CSV do SPICE
python3 compara_spice.py Thebe Io --tmax 100 --pasta-spice saida_elementos
```

Lê `NOME_osc.dat` (nbody) e `elementos_NOME.csv` (gerado pelo `elementos_spice.py`)
e salva `NOME_spice.png` (os 6 elementos sobrepostos: nbody em azul, SPICE em
vermelho) e `NOME_spice_dif.png` (nbody − SPICE, com o SPICE interpolado nos
tempos do nbody). O t = 0 do nbody é a primeira linha do CSV (a época `INICIO`);
use `--epoca-et` se for outra. Use `FRAME = "EQUADOR_JUPITER"` e
`SOMAR_GM_LUA = True` no `elementos_spice.py`, para ficar no mesmo referencial e
com o mesmo μ do nbody.

### Polo do frame

O polo IAU (pck) não é exatamente o eixo em torno do qual as luas precessionam no
`jup365`. Com o polo IAU, a inclinação das luas pequenas no SPICE oscila com o período
do nodo e a do nbody não. O `polo_spice.py` mede o deslocamento (círculos do vetor
inclinação deslocados juntos); copie o `(x, y)` que ele imprime para `POLO_DINAMICO`
no `elementos_spice.py` e rode de novo `elementos_spice.py` → `spice_para_in.py` →
`nbody` → `compara_spice.py`. O valor vai no cabeçalho do `estados_iniciais.csv` e o
`spice_para_in.py` usa o mesmo frame para o Sol.

Com `--janela 200`, o `polo_spice.py` também ajusta o círculo em janelas de tempo e salva
`polo_spice_tempo.png` (centro × tempo). Se o centro das luas pequenas do SPICE andar junto
ao longo dos anos, é o polo de Júpiter se movendo na efeméride (precessão/nutação), que o
nbody não tem (J2 fixo no eixo z); aí `POLO_DINAMICO` só remove a média.
