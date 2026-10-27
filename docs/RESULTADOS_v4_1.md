# Resultados da v4.1: rodada principal (2020–2025) e rodadas P1/P2 (30/09/2026)

Todos os números saem dos arquivos das três rodadas:
- `outputs_v4*/manuscript_numbers.json`;
- `outputs_v4*/tables/`;
- `outputs_v4*/gee/PROVENANCE.json`, com o SHA-256 de cada arquivo e as propriedades de cada asset do Earth Engine;
- `outputs_v4*/gee/gee_tasks.csv`, com o custo de cada tarefa.

As conferências de desenvolvimento estão em `validacao_v4_1/`. Este documento serve de guia para reescrever o manuscrito; os valores oficiais são os dos arquivos.

**Histórico de revisão.** O documento passou por três auditorias independentes com verificação adversarial:
- **1ª auditoria:** 34 de 37 apontamentos confirmados. O principal foi um estimador de pareamento com peso errado.
- **2ª auditoria:** 32 de 36 confirmados. O principal foi o bootstrap do pareamento, que não reamostrava os blocos dos controles; com a correção, a maior parte dos resultados pareados perdeu significância. Também foram corrigidas sínteses que contradiziam as tabelas.
- **3ª auditoria:** 15 de 15 confirmados. Oito eram de redação ou arredondamento. Sete eram de código e afetavam intervalos em casos de borda; nenhum muda as diferenças de média.

O código foi corrigido e as análises refeitas. O que segue incorpora as três rodadas.

## 1. Rodadas

| | Principal | P1 | P2 |
|---|---|---|---|
| Período | 2020-01-01 a 2025-12-31 (72 meses) | 2020-01-01 a 2022-07-16 (31 meses) | 2022-07-17 a 2025-12-31 (42 meses) |
| Versões do processador (órbitas) | 10 | 4 (1.3.2, 1.4.0, 2.2.0, 2.3.1) | 7 (2.3.1 [2], 2.4.0 … 02.09.01) |
| Mínimo de recuperações por célula / de meses | 50 / 24 | 22 / 11 | 29 / 14 |
| Células do NE com dado / na tela estrita | 32.582 / 21.537 | 31.158 / 16.412 | 32.484 / 23.997 |
| Retiradas por contagem baixa | 10.399 (32%) | 14.727 (47%) | 7.634 (24%) |
| Recuperações por célula no NE (mediana) | 100 | 29 | 71 |
| Limiar de detecção (realce, ppb) | 8,8 | 10,1 | 8,9 |
| Células a 0–5 km de um duto (tela estrita) | 168 | 80 | 229 |
| Blocos de 2° efetivos dessas células* | 3,8 | 1,7 | 4,3 |
| EECU-h (`gee_tasks.csv`) | 20,3, dos quais 12,1 em tentativas que falharam ou foram canceladas | 3,7 | 8,7 |

\* 1 / Herfindahl das participações das células próximas nos blocos de 2°. Com menos de ~10 grupos tratados, os erros-padrão robustos por grupo e os intervalos de bootstrap ficam otimistas (Cameron, Gelbach & Miller 2008; MacKinnon & Webb 2017). **P1 não sustenta inferência perto × longe dos dutos.**

- **O que P1 e P2 são.** Não são períodos homogêneos de processador: são os períodos **antes e depois da mudança de correção de viés de 17/07/2022**. P2 também é um subconjunto da rodada principal, com cerca de 2/3 das suas recuperações por célula, então a concordância entre a principal e P2 não é uma réplica independente. A réplica é **P1 × P2** (§3.5, Tabela S19).
- **Por que P1 e P2 usam outro mínimo.** `--min-obs` mantém o **limiar** por mês da regra principal (50 em 72 meses). A precisão por célula e a amostra retida, porém, diferem entre as rodadas: veja as linhas de retirada por contagem e de limiar acima.

## 2. O que mudou no código

### 2.1 Desempenho e robustez da exportação (4.1.0)

Cada item foi comparado com dados reais contra a versão anterior (`validacao_v4_1/LEIA-ME.md`).

1. **Estatísticas por célula numa passagem.** As 17 bandas são idênticas (diferença 0) a uma referência independente com passagens separadas.
2. **Um só `reduceResolution` por órbita.** Mesmas células e diferença 0 em 4 órbitas reais.
3. **Banda `n_months`** com `countDistinctNonNull`. O `countDistinct` comum soma +1 em todas as células: 255 de 255 na caixa de teste. `n_months_equals_monthly_stack` dá `true` nas três rodadas.
4. **Pilha mensal numa passagem,** com um asset por ano. Diferença 0 nas 24 bandas.
5. **Covariáveis calculadas com `computePixels`** em 88 blocos de 32 × 32 células. A exportação em lote estourava o limite de 2^31 pixels.
6. **Exportação que pode ser retomada.** Reaproveita assets e tarefas em execução. O `PROVENANCE.json` registra as propriedades `ch4ne_*` de cada asset. Se a versão do código mudou, o programa avisa: na rodada principal, `legacy` e `monthly_background` são da 4.0.0, com as mesmas configurações. `gee_tasks.csv` registra o custo de cada tarefa. Com `--overwrite`, a execução para antes da fase A se houver tarefa da rodada em andamento.

### 2.2 Correções vindas das auditorias

- **Pareamento (Tabela 7):**
  - estimador por conjunto (Abadie & Imbens 2006);
  - correção de viés com β estimado nos **controles pareados** (principal) e em todas as células distantes (sensibilidade), Abadie & Imbens 2011;
  - **bootstrap conjunto** dos blocos das células próximas e dos controles;
  - caliper com o DP combinado dos dois grupos (Austin 2011);
  - sem intervalo com menos de 5 blocos;
  - relatório da concentração de peso nos controles e do suporte comum do escore.
- **Contraste estratificado por faixa de distância à costa** (Tabela S18). Mais simples e mais transparente que o pareamento. Também há a versão por faixa × estado.
- **Atribuição (Tabela 4).** Estratificada pela costa; estratos não informativos ficam de fora; com um só estrato informativo, usa o teste exato de Fisher.
- **Diagnósticos:**
  - S14: espelho por estrato, com o corte de densidade na mediana da rodada e num corte fixo de 3 × `min_obs`;
  - S15: simetria de sinal por quintis de `n_obs` e por faixa de costa, também centrada na mediana de cada estrato;
  - S16: estimativas deixando um bloco de fora;
  - S17: clusters significativos das regras de sensibilidade;
  - S7: as regras escalam com `--min-obs` e trazem a probabilidade de nenhum pico significativo cair a < 90 km de um duto por acaso;
  - S19: contrastes perto × longe nas células comuns a P1 e P2.
- **Tendência:** o erro-padrão do modelo com degraus e uma especificação alternativa com mudança de inclinação, comparadas por AIC.
- **Robustez do código:**
  - estados sem células;
  - pareamento vazio;
  - `legacy.tif` só lido se estiver no `PROVENANCE`;
  - a opção `--jobs` é validada antes da fase A;
  - `wait` com lista exata de tarefas opcionais;
  - a Tabela S17 é gravada com cabeçalho mesmo vazia;
  - `--overwrite --jobs` refaz a fase A só quando pedido;
  - a lista de custos vazia não interrompe o download;
  - os JSON são estritamente válidos (sem `NaN`).
- **Intervalos em casos de borda** (3ª auditoria):
  - o bootstrap exige pelo menos 5 blocos **das células próximas**;
  - nenhum intervalo é dado se mais de 5% das réplicas forem indefinidas;
  - os pesos dos controles são normalizados em cada conjunto sorteado;
  - a correção de viés pelos controles pareados exige ≥ 5 controles distintos por parâmetro e posto completo;
  - a estratificação conta só os blocos dos estratos informativos.

  Com isso, o pareamento euclidiano (Tabela S13) e o estrato costa × estado de P1 ficam sem intervalo.

## 3. Resultados

Nesta seção, **±** indica o **erro-padrão (EP)**, salvo quando se escreve IC 95%.

### 3.1 Background regional (Tabela S5, Figura 3)

| Modelo (rodada principal, 71 meses; 2020-03 descartado) | Tendência | AIC |
|---|---|---|
| Tendência + harmônicos | 10,1 ppb/ano (EP Newey–West 0,6; IC 95% ± 1,2) | 462,5 |
| + degraus impostos em 2021-07 e 2022-07 | 6,7 (EP 0,9); degraus +5,9 (EP 3,1) e +9,6 (EP 2,4) | 443,7 |
| + mudança de inclinação em 2022-07-17 | 15,5 (EP 1,0) antes; variação −8,6 (EP 1,6) depois | **431,5** |

- **Os dados não separam degrau de mudança de ritmo.** O modelo com mudança de inclinação é o preferido pelo AIC, e os ajustes por período concordam com ele: P1 dá 12,7 ppb/ano (EP 1,2) e P2, 5,8 (EP 0,6).
- **O degrau de 2021-07 não é identificável:** dá +5,9 na rodada principal e −7,0 (EP 3,2) dentro de P1.
- **A tendência de ~10 ppb/ano do período inteiro é uma média** de ritmos diferentes, misturada com as mudanças de produto.
- **Amplitude sazonal: 6,7 ppb.** Ela é afetada pela mudança sazonal da área observada: a latitude mediana das células observadas varia de −13,2° a −7,1°.

### 3.2 Estados (Tabela 2, Figura 6)

Anomalia sem tendência, tela estrita, rodada principal:

| RN | PB | CE | MA | AL | PE | PI | BA | SE |
|---|---|---|---|---|---|---|---|---|
| +6,3 | +5,3 | +5,1 | +3,4 | +2,0 | +2,0 | −1,0 | −3,5 | −4,1 |

- **Gradiente norte–sul.** RN > BA se repete em P1 e P2 (RN +6,3/+6,2; BA −3,8/−3,4). Os outros estados mudam de ordem; em P1, AL tem 11 células e SE, 1.
- **XCH4 bruto.** O ranking bruto põe AL em primeiro por causa do horário de amostragem (`bg_sampled`), então o gradiente deve ser descrito com a anomalia.
- **Maior DP espacial:** SE (6,86 ppb, n = 84).

### 3.3 Detecção de clusters (Tabela 3, Figuras 4 e 8, Tabelas S6–S8, S14, S15, S17)

- **Regra principal:** 86 clusters positivos e 49 do espelho, todos com q entre 0,54 e 0,58; **nenhum com q ≤ 0,10**. P1 dá 67 e 20; P2, 87 e 63; também nenhum.
- **Clusters significativos só aparecem nas regras mais densas:**
  - na principal, com ≥ 200 recuperações: 18–21 clusters conforme o z;
  - no estrato de alta densidade de P1 (Tabela S14): 25 clusters.

  O resultado é frágil e depende do corte. No estrato de alta densidade da principal, o corte na mediana (141 recuperações) dá 0 clusters com q ≤ 0,10, e o corte fixo em 150 (3 × `min_obs`) dá 22. Em P1, os cortes 54 e 66 dão 25 e 18 (Tabela S14).
- **O espelho não é calibrado** (Tabela S15). Ele supõe desvios positivos e negativos igualmente prováveis sem emissão:
  - nos dois quintis superiores de densidade, z > 2,5 é 9 a 12 vezes mais frequente que z < −2,5 (principal), e até 28 vezes em P1;
  - a 0–50 km da costa, é 4,6 a 13 vezes menos frequente;
  - centrando cada estrato na sua própria mediana, a assimetria costeira diminui mas não some (razões 0,56 e 0,17). A de alta densidade cai de 9–12 para ~8,6, mas continua muito acima de 1.

  Os q-valores são uma estimativa heurística.
- **Poder perto dos dutos nas regras densas.** Os clusters significativos das regras de sensibilidade (Tabela S17) ficam ≥ 90,8 km dos dutos. Essas regras, porém, retêm poucas células perto de dutos: na regra de ≥ 200, 33 de 5.123 a ≤ 20 km e 3,4% a < 90 km. A chance de nenhum dos 18–21 picos cair a < 90 km por acaso é de 0,48–0,53 (Tabela S7). **O resultado não é informativo sobre os dutos.**
- **Onde estão os clusters da regra principal:**
  - 20 dos 86 picos ficam na borda oeste da BA (a oeste de 45,4°W, 10–15°S);
  - os **6 mais fortes ficam entre 45,0°W e 38,6°W**, fora dessa faixa, com 55–161 recuperações no pico (a mediana da tela é 141);
  - a mediana da distância do pico ao duto mais próximo é 464 km.
- **Três picos ficam a ≤ 5 km de um duto:**
  - cluster 4, no CE: **CEGAS, 0,16 km**;
  - cluster 18, na BA: Bahiagás, 2,8 km;
  - cluster 63, no RN: BEN-A/ET-B, 1,5 km.

  Os três estão a 11–15 km da costa, com q = 0,54.
- **O cluster 4 é o principal candidato a acompanhamento:**
  - é o 4º mais forte (z = 4,97; realce +17,3 ppb);
  - fica na região metropolitana de Fortaleza, sobre a linha da CEGAS;
  - tem a maior anomalia mensal média entre os 6 mais fortes (+18,0 ppb; +19,2 ponderando pelas recuperações);
  - a anomalia é positiva em 96% dos meses com dado;
  - reaparece por inteiro em P2 (7 de 7 células no cluster);
  - em P1 tem z entre 2,1 e 4,1, mas só 1 célula passa na tela, por falta de meses.

  Na mesma área podem existir outras fontes (urbanas, aterros, indústria). A atribuição exige dados de maior resolução.
- **Persistência.** 42 dos 67 clusters de P1 se sobrepõem a um cluster de P2 (Tabela S12). Dos 6 mais fortes da rodada principal, 3 não podem ser avaliados em P1, porque quase nenhuma célula passa na tela.
- **Séries temporais** (Tabela S6, Figura 8). A anomalia regional dos 6 mais fortes é positiva em 77–96% dos meses com dado. As médias ficam entre +7 e +18 ppb (média simples; +9 a +19 ponderada). 87% dos meses com dado caem entre maio e novembro, e o cluster 3 só tem dados a partir de 2022-06. Médias mensais não resolvem intermitência em escala menor que um mês.

### 3.4 Atribuição aos dutos (Tabela 4)

- **Hipótese nula de área** (rodada principal, 5 km): 3 picos observados contra 0,67 esperados (enriquecimento 4,5; p binomial = 0,03). Essa nula ignora a variação espacial da detectabilidade.
- **Contra o espelho sem estratificação:** OR = 0,26 (p = 0,99). O resultado é confundido pela costa: os picos do espelho são costeiros e os positivos, interiores.
- **Estratificado pela costa.** Só o estrato ≤ 50 km é informativo, então vale o teste exato de Fisher:
  - principal: 3 de 4 picos positivos a ≤ 5 km, contra 6 de 21 do espelho (OR 7,5; p = 0,12);
  - P2: 4 de 6 contra 8 de 23 (OR 3,75; p = 0,20);
  - P1: nenhum pico a ≤ 20 km na tela estrita. O cluster 4 existe em P1, mas não passa na tela (§3.3).

**Conclusão: inconclusivo.** São 4 a 6 picos positivos costeiros.

### 3.5 Proximidade (Tabelas 5–8, S4, S9–S11, S13, S16, S18, S19; Figura 7)

- **Confundimento pelo litoral** (Tabela S9). 78% das células a ≤ 5 km de um duto estão a < 50 km da costa, contra 10% no NE. A distância mediana à costa é 28 km, contra 299 km.

**Diferença perto (0–5 km) − longe (≥ 50 km), em ppb, rodada principal.** IC 95% entre parênteses; * = IC exclui zero, ou p < 0,05. Os IC e os p são otimistas, porque há 3,8 blocos efetivos.

| Variável | Anéis, sem ajuste [faixa deixando 1 bloco de fora] | Regressão com covariáveis (coef.; p) | Pareamento com correção de viés | Estratificado por faixa de costa | Por faixa de costa × estado (82 de 168 células; 5 blocos) |
|---|---|---|---|---|---|
| Anomalia média | +1,2 (−10,4 a 7,1) [−3,3; +3,5] | +0,4; 0,76 | +2,2 (−2,5 a 9,2) | −3,6 (−13,4 a 26,6) | +1,6 (−2,4 a 17,5) |
| Realce local | −1,7 (−6,0 a 1,7) [−2,5; −0,9] | +0,3; 0,83 | −3,0 (−5,8 a 2,4) | −1,4 (−5,1 a 8,5) | +0,5 (−3,2 a 9,8) |
| DP temporal | +4,6* (1,0 a 12,1) [+3,2; +8,3] | +2,0; 0,13 | −1,7 (−4,3 a 0,6) | +1,7 (−13,5 a 7,5) | +0,1 (−6,2 a 1,2) |
| Cauda superior | +3,7* (1,5 a 8,2) [+2,8; +6,1] | +1,8*; 0,024 | +0,5 (−2,0 a 3,0) | +1,5 (−12,6 a 5,2) | +0,8 (−5,1 a 1,6) |
| Cauda inferior | +7,1* (1,0 a 20,2) [+4,8; +13,3] | +2,4; 0,33 | −5,5* (−13,5 a −0,7) | +2,0 (−23,9 a 12,2) | −1,7 (−20,7 a 0,1) |
| Assimetria (sup. − inf.) | −3,3 (−12,1 a 0,7)† [−7,2; −1,9] | −0,6; 0,72 | +6,1* (1,6 a 12,8) | −0,5 (−7,6 a 9,8) | +2,5* (1,3 a 16,1) |

† O IC exclui zero só com blocos de 1° (Tabela S10).

**Pareamento** (Tabelas 7, S4, S13):
- 158 de 168 células foram pareadas, mas o suporte comum é mínimo: 99% das células próximas têm escore acima do percentil 99 das distantes.
- 3 controles do litoral sul da BA, a 10,7–16,5 km da costa (coordenadas em `manuscript_numbers.json`), carregam 28% do peso. Os conjuntos têm tamanhos {1: 15, 2: 29, 3: 3, 5: 111}, e o número efetivo de controles é 29,5.
- Equilíbrio residual (DMP): água −0,62, albedo +0,49, log n_obs +0,27, DP do relevo +0,20, AOD −0,17, fração jan–jun −0,14.
- **O pareamento depende do modelo de ajuste.** A anomalia média corrigida vai de +2,2 (β nos pareados) a +4,4 (β em todas as distantes); nenhuma exclui zero.

**P1 e P2** (mesmas colunas):
- **P2** (4,3 blocos):
  - os anéis mostram de novo DP e caudas maiores perto dos dutos, com realce −4,1* (−9,9 a −0,2);
  - o pareamento só é significativo no realce, e negativo: −3,6* (−5,9 a −0,2);
  - as duas estratificações dão tudo nulo.
- **P1** (1,7 blocos, só descritivo):
  - **nos anéis**, a cauda superior é mais larga (+1,4*) e a inferior não muda (−0,6); assimetria +2,1* (1,1 a 8,0);
  - **nos estimadores ajustados**, a assimetria positiva vem sobretudo de uma cauda **inferior mais estreita**:
    - regressão: superior +0,9 (p = 0,21), inferior −2,6 (p < 0,001), assimetria +3,5 (p < 0,001);
    - estratificação por faixa de costa: inferior −1,8* (−4,5 a −0,7);
    - costa × estado: 20 células em 4 blocos, sem intervalo;
  - a regressão também dá anomalia média +4,2 (p = 0,049);
  - no pareamento, assimetria +1,9 e anomalia −2,6, ambas não significativas.

  Em P1 o padrão "só a cauda superior" aparece apenas nos anéis sem ajuste, e com menos de 2 blocos efetivos.

**Réplica P1 × P2 nas mesmas células** (Tabela S19; 78 células próximas e 15.747 distantes que passam na tela nos dois períodos):

| | Anomalia média | Realce local | Assimetria |
|---|---|---|---|
| P1 | +8,8 | −0,1 | +1,7 |
| P2 | +6,7 | −0,7 | +0,0 |

Na amostra comum, a diferença de anomalia é parecida nos dois períodos. Ela reflete o gradiente regional, já que as células próximas estão no RN e no CE (§3.2), e não um realce local. A mudança da diferença nas amostras completas (P1 +8,4; P2 −0,9) é sobretudo efeito de composição: as células que só passam na tela em P2 ficam sobretudo em SE, AL e BA.

**Leitura:**
1. **Realce local.** Nenhum estimador mostra realce local positivo perto dos dutos, em nenhum período; em P2 o realce é até ligeiramente negativo.
2. **Anomalia regional.** Na principal e em P2, nenhum estimador ajustado mostra diferença significativa. Em P1, a regressão dá +4,2 (p = 0,049), com 1,7 blocos efetivos. A diferença das células comuns vem do gradiente norte–sul.
3. **Dispersão e caudas.**
   - Contra todas as células distantes, as duas caudas são mais largas perto dos dutos, e a inferior mais, com assimetria negativa, na principal e em P2.
   - Estratificando por faixa de costa, tudo fica nulo na principal e em P2, compatível com um efeito de costa. A exceção é P1, em que a cauda inferior fica mais estreita (−1,8; IC −4,5 a −0,7).
   - O pareamento e o estrato costa × estado dão assimetria positiva, mas com controles concentrados e pouca sobreposição.
   - Em P1, os anéis mostram o padrão de emissão intermitente (só a cauda superior). Os estimadores ajustados atribuem a assimetria a uma cauda inferior mais estreita. Tudo com menos de 2 blocos efetivos.
   - **Conclusão: a evidência de assinatura intermitente é fraca e inconsistente entre estimadores, e não pode ser afirmada nem descartada.**

### 3.6 Reprodução da v3 (Tabela S2)

- **μ:** 1880,15 ppb, contra 1880,14 publicado.
- **σ:** 9,05 ppb ponderado (−3,7%) e 9,09 pelo centroide, contra 9,39.
- **n:** 32.263 pelo centroide e 32.894 por toque, contra 32.904.

As duas checagens ficam `false` por pouco. Causas prováveis: mudanças na coleção depois da rodada da v3, e a emulação aproximada da ponderação de polígonos. O composto v3 (`legacy`) foi exportado pela 4.0.0. Isso não afeta a v4.

## 4. Síntese para o manuscrito (formulação sugerida)

1. **Nenhum realce local de XCH4 perto dos dutos.** Nenhuma comparação mostra realce local positivo na escala do TROPOMI agregado a ~7 km e em médias de vários anos. Isso vale para a rodada principal e para os períodos antes e depois da mudança de 17/07/2022, com qualquer estimador. Nenhum cluster da regra principal passa na estimativa de falsas descobertas, que é heurística (§3.3).
2. **Ausência de evidência, não evidência de ausência.**
   - Os intervalos para a anomalia média chegam a ~±10 ppb nos anéis e no pareamento. Na estratificação por faixa de costa vão de −13 a +27 ppb, da ordem ou acima do limiar de detecção de 8,8 ppb.
   - A inferência perto dos dutos conta com ~4 blocos espaciais independentes.
   - O pareamento tem suporte comum mínimo.
3. **Dispersão dominada pela costa.** As diferenças de dispersão perto dos dutos seguem a posição costeira deles e desaparecem na estratificação por faixa de costa, na principal e em P2. A assimetria positiva das caudas aparece só em estimadores frágeis: o pareamento, com suporte mínimo, e P1, com menos de 2 blocos efetivos e vindo sobretudo de uma cauda inferior mais estreita. Por isso a assinatura de emissão intermitente é inconclusiva.
4. **O cluster 4 merece acompanhamento.** Fica na região metropolitana de Fortaleza, a 0,16 km da CEGAS, e é persistente em P2. Deve ser citado como candidato para observação de maior resolução, sem atribuição.
5. **Clusters do interior.** A maioria dos clusters fica no interior, longe dos dutos. Os que só são significativos nas regras densas indicam efeitos de superfície ou relevo, ou fontes não relacionadas a O&G, a investigar.
6. **Crescimento regional.** Mudou em torno de jul/2022 (~15 → ~6–7 ppb/ano), e os dados não o separam da mudança de produto.

## 5. Pontos de atenção

- **Espelho não calibrado** (Tabela S15); descrever os q-valores como heurísticos.
- **Background mensal afetado pela cobertura.** Alternativa futura: calculá-lo sobre um conjunto fixo de células bem observadas.
- **Poucos blocos independentes.** Uma inferência mais conservadora seria bootstrap de grupos selvagem (*wild cluster bootstrap*) ou inferência por aleatorização.
- **Limitações que continuam:** sem vento; sem fontes naturais (manguezais, reservatórios, queimadas no MATOPIBA); sem colunas de maior resolução (Sentinel-2, EMIT, PRISMA), a começar pelo cluster 4.

## 6. Pendências

- [ ] Reescrever o manuscrito com estes números (ver [PENDENCIAS_ATUALIZADO.md](PENDENCIAS_ATUALIZADO.md)).
- [ ] Arquivar no Zenodo, em nova versão, `outputs_v4*/gee/` (rasters, `PROVENANCE.json`, `gee_tasks.csv`) e as tabelas e figuras das três rodadas.
