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
| `plot_elementos.py` | plota osculadores × geométricos |

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
* Dá para sobrescrever os valores do banco no próprio `.in` (`gm_planeta`,
  `raio_planeta`, `j2`, `j4`), ou usar `planeta = personalizado` e dar todos eles.

## Saída

Para cada satélite, `NOME_osc.dat` e `NOME_geo.dat`, com 7 colunas:

```
tempo(dias)  a(km)  e  inc(deg)  w(deg)  om(deg)  M(deg)
```

Esse é o mesmo formato do `sat1_geo.dat` do `osc_geo.f90`. Os elementos
osculadores usam μ = G(M + m); os geométricos usam G·M, J2 e J4 do planeta (J6 = 0),
como no `geo` original. O programa também mostra na tela o erro relativo da
energia, para controlar a qualidade da integração.

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
