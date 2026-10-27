# Auditoria do código e do processamento da versão v3

Artigo: *Satellite Screening of Methane Anomalies over Onshore Oil and Gas Infrastructure in Northeast Brazil* (Cardoso, Antonino, Reichert; submissão ao *Remote Sensing*, MDPI).

Objeto: `codigo para o artigo 7 v3.ipynb` (idêntico ao depositado no Zenodo 10.5281/zenodo.18717652, v1.0.0), as tabelas e figuras publicadas e o manuscrito `Cardoso_et_al_Methane_TROPOMI_NE_Brazil_RemoteSensing_MDPI.docx`.

Data: 23/09/2026.

**Como os itens foram verificados.** Cada problema traz a evidência. Há três tipos de evidência:
- **[código]**: leitura linha a linha do notebook;
- **[dados reais]**: recálculo a partir das tabelas publicadas, das camadas ANP e dos limites estaduais (`validacao_offline/validacao_offline.py` e `tests/test_real_inputs.py` reproduzem tudo);
- **[documentação GEE]**: documentação oficial e código-fonte do `earthengine-api`.

Não se usou nenhum dado sintético ou inventado.

---

## A. Erros que mudam números publicados

| # | Problema | Evidência | Impacto |
|---|---|---|---|
| A1 | **Estados trocados.** Os polígonos foram enviados ao GEE sem atributos (`gdf_to_ee(ne_gdf[['geometry']])`). Os nomes foram depois reconstruídos por um índice na ordem IBGE, mas o GeoJSON está em ordem alfabética (AL, BA, CE, MA, PB, PE, PI, RN, SE). | [dados reais] Com os rótulos corrigidos, 19 de 20 picos cabem no máximo do seu estado; com os originais, só 11 de 20. | 6 dos 9 estados receberam estatísticas de outro estado (Tabela 2, Figura 6/violino). Já estava identificado no PENDENCIAS; a correção está confirmada. |
| A2 | **Todas as distâncias "0 km" são artefato.** `ee.Join.saveBest` grava a distância no objeto *casado* (`nearest_infra`), mas o código lia `props.get('dist_m', 0)` do objeto primário e recebia sempre o valor-padrão 0. | [documentação GEE] "measureKey: the key used to save the measure of the join condition on the match". [dados reais] Distâncias reais dos 4 centroides ao duto onshore mais próximo: **4,9; 12,4; 14,0; 19,2 km**. | Não houve nenhuma "interseção direta". O join só prova que os 4 centroides estão a ≤ 20 km de um duto. A Tabela 3 e as Seções 3.3, 4.3, 4.7 e 5 precisam ser reescritas. |
| A3 | **O ranking não é de pixels.** A "Tabela 5 – vinte anomalias mais intensas" é o máximo de cada *cluster* (`reduceRegions` sobre os polígonos), com a latitude/longitude do **centroide** do polígono, não do pico. São 20 dos 21 clusters. | [código] célula 10; `Target_ID` = 1 em todas as linhas (é a propriedade `label` do `reduceToVectors`). O rank 12 tem centroide na PB, mas o valor 1904,28 é exatamente o máximo de PE. | A legenda, o texto ("twenty most intense anomalous pixels") e a coluna de estado da Tabela 5 estão errados. O "alvo 1" da série temporal é um centroide, não o pixel de pico. |
| A4 | **μ/σ calculados a 14 km e aplicados a 7 km.** | [código] `reduceRegion(scale=14000)`. | Limiar 1899,63 em vez de 1898,91 ppb (já estava no PENDENCIAS). |
| A5 | **Inventário: 417 segmentos "no domínio de estudo" é incorreto.** A seleção usou os estados do NE com um buffer de 1° (~110 km). | [dados reais] Dos 417: **262 onshore no NE, 127 majoritariamente offshore, 28 em outro estado (Espírito Santo, 18,8–20,2° S)**. | Tabela 1 e resumo. A Figura 2 da v3 mostra dutos até −21°. |
| A6 | **Figura de clusters vazia no Zenodo.** O download do thumbnail falhou ("Computation timed out"), e o arquivo `Figure5_Final_Vectorized_Clusters.png` do Zenodo não tem nenhum cluster. A figura do manuscrito veio de outra execução e mostra dezenas de manchas, incompatíveis com "21 clusters". | [código] saída da célula 6. | Figura 3 do manuscrito não é reprodutível a partir do código depositado. |

## B. Problemas de método que invalidam interpretações

| # | Problema | Evidência | Consequência |
|---|---|---|---|
| B1 | **O "sinal de variância" mede a coisa errada.** O desvio-padrão da Tabela 4 é o desvio **espacial** entre pixels de uma média de 6 anos, dentro de um buffer. O texto o interpreta como variabilidade **temporal** (emissões intermitentes "diluídas na média, mas preservadas na variância"). Uma média temporal não guarda a variância temporal. | [código] `reduceRegion(stdDev)` sobre `mean_methane`. | A principal conclusão do artigo não é sustentada pela estatística calculada. A v4 calcula a variância temporal real de cada pixel (`anom_sd`) e a frequência de desvios fortes positivos **e** negativos (`frac_hi`, `frac_lo`): emissão intermitente só aumentaria a cauda superior. |
| B2 | **Confundimento com o litoral.** | [dados reais] 78% das células a ≤ 5 km de um duto estão a < 50 km do litoral, contra 10% de todas as células do NE (distância mediana ao litoral: 28 km contra 299 km). O TROPOMI tem artefatos documentados em pixels mistos terra/água (PRF, Balasus et al. 2023). | Um desvio-padrão maior perto dos dutos é esperado mesmo sem nenhuma emissão. A v4 compara com controles pareados por distância ao litoral, albedo, água, altitude e densidade de dados, e usa máscara costeira. |
| B3 | **Zonas aninhadas.** O "background regional" inclui as zonas de 5 e 15 km, e a zona de 5 km está contida na de 15 km. O texto diz "remaining pixels". | [código] célula 9. | Comparações não independentes. A v4 usa anéis exclusivos (0–5, 5–15, 15–50, ≥50 km). |
| B4 | **Sem remoção de tendência.** O XCH₄ cresceu cerca de 8–15 ppb/ano no período. Um pixel observado mais em 2025 que em 2020 parece "anômalo". | [dados reais] o próprio manuscrito relata ~7,6 ppb/ano nos alvos 2 e 3. | A v4 subtrai o background regional mensal de cada órbita antes de compor. |
| B5 | **Sem remoção do gradiente espacial de grande escala.** O mapa Z da v3 é alto no norte e baixo no sul da BA: um Z contra uma média regional única marca como "hotspot" o gradiente latitudinal. | [dados reais] Figura 4 publicada. | A v4 remove um background local (mediana móvel de ~220 km) antes do limiar. |
| B6 | **Sem número mínimo de observações.** O alvo 1 (pico máximo, 1937,9 ppb) tem uma única observação mensal válida na série temporal. | [manuscrito] Seção 4.4. | Picos podem vir de poucas medições. A v4 exige ≥ 50 recuperações e ≥ 24 meses com dados. |
| B7 | **A "limpeza morfológica" não removeu ruído.** `z.gt(2).selfMask().focal_min(1).focal_max(1)`: os redutores do GEE ignoram pixels mascarados, então `focal_min` nunca vê um zero e não erode. O resultado é a identidade ou uma dilatação, nunca a "closing" do texto (nem uma abertura). Os kernels ainda rodaram na escala de vetorização (10 km), não na grade de 7 km. | [documentação GEE] | A Seção 2.5 descreve uma operação que não ocorreu. Os "21 clusters" dependem da escala de 10 km. A v4 rotula componentes conexos na grade nativa, com tamanho mínimo explícito. |
| B8 | **Z > 2 sem controle de falsos positivos.** Num campo espacialmente correlacionado, Z > 2 gera muitas manchas por acaso. | [método] | A v4 estima a taxa de falsas descobertas pelo "espelho" (z < −limiar) e dá um q-valor a cada cluster. |
| B9 | **Atribuição sem taxa-base.** "4 de 21 clusters perto de dutos" não diz nada sem a fração da área que fica perto de dutos. | [método] | A v4 compara o número observado com o esperado ao acaso (Tabela 4: enriquecimento e p binomial). |

## C. Descrições do Método que não correspondem ao código

| # | Texto do manuscrito | O que o código faz |
|---|---|---|
| C1 | "uncertainty exceeding 15%" | Filtro de incerteza ≤ **15 ppb** (valor absoluto). Na prática não remove quase nada: o GEE já mantém só qa_value > 0,5 na ingestão, o que elimina ruído > 10 ppb. |
| C2 | "low quality-assurance flags" removidos pelos autores | Não há filtro de QA no código; ele é aplicado pelo GEE na conversão para L3. |
| C3 | "reprojected to a fixed 7 km grid… no artificial spatial detail" | `.reproject(scale=7000)` sem `reduceResolution`: amostragem por vizinho mais próximo de um nível da pirâmide, não média de área. |
| C4 | "morphological closing" | Ver B7. |
| C5 | "Euclidean distance from the cluster centroid" | Ver A2: a distância nunca foi lida. |
| C6 | Figura 7: "Grey points are individual TROPOMI observations" | São pixels do composto médio, amostrados a 1 km de uma imagem de 7 km (até 49 amostras repetidas por pixel). |
| C7 | Figura 8: "monthly medians… 2020–2024" | Mediana mensal da coleção **sem** o filtro de incerteza, extraída no centroide do polígono, com 180 chamadas `getInfo` sequenciais. |
| C8 | "Z > 2 ≈ upper 2.5%" | P(Z > 2) = 2,28%. |
| C9 | Composto "2020–2025" | `filterDate('2020-01-01','2025-12-31')` exclui 31/12/2025 (fim exclusivo). |
| C10 | Coleção homogênea | A coleção OFFL do GEE mistura versões do processador (correção de viés mudou em 01/07/2021; correção espectral de superfície desde 17/07/2022; lacuna de dados 26/07–31/08/2022). A v4 registra as versões e permite rodar só o período homogêneo (`--start 2022-07-17 --label P2`). |

## D. Inconsistências internas do manuscrito

| # | Onde | Problema |
|---|---|---|
| D1 | Seção 4.2 | Diz que RN tem o maior desvio-padrão e a cauda baixa mais longa. Após a correção dos estados, o maior desvio é do MA; RN tem um dos menores (6,56 ppb). Contradiz a Seção 3.2. |
| D2 | Legenda da Fig. 5 | Usa 1880,01 / 1899,63 ppb (14 km); o texto usa 1880,14 / 1898,91 ppb (7 km). |
| D3 | Tabela 2 | Anomalia calculada contra 1880,01; o restante do texto usa 1880,14. |
| D4 | Resumo e 4.5 | 57,9 ppb (contra 1880,01) × 57,79 ppb na Tabela 5 (contra 1880,14). |
| D5 | Seções 3.3, 4.3, 5 | "Direct intersection (0 km)" (ver A2). |
| D6 | Tabela 5 | "most intense anomalous pixels" (ver A3); estado do rank 12 é do centroide, não do pico. |
| D7 | Conclusões | "prioritised target list for field campaigns" apoiada em A2/A3. |

## E. Figuras

- Títulos e números embutidos nas imagens (a MDPI pede sem título).
- Figura 1 da v3: paleta tipo *rainbow* com barra de cores de outra paleta (`jet`), eixos desligados, sem barra de escala.
- Thumbnails do GEE desenhados com interpolação bilinear (mistura cores e cria bordas semitransparentes que não são dados).
- A barra de escala de "1,8° = 200 km" varia com a latitude.
- A v4 gera todas as figuras a partir dos valores exportados (GeoTIFF), com barra de escala geodésica, seta norte, grade em °S/°W, paletas acessíveis e 600 dpi + PDF vetorial.

## F. Reprodutibilidade

- Resultados dependiam de estado do notebook (variáveis de células anteriores, `mu`/`sigma` com *fallback* silencioso 1880/9, `except:` genérico).
- Nenhum registro de versões, datas de execução ou proveniência.
- CSVs com `;` e rótulos vazios (`UF` vazio, `Plume_ID = ??`).
- A v4 corrige tudo isso:
  - configuração única em `ch4ne/config.py`;
  - exportação em lote para assets, sem limite de tempo;
  - registro `PROVENANCE.json` com checksums; a análise se recusa a rodar sem ele;
  - testes automáticos com dados reais;
  - tabelas em CSV padrão;
  - `manuscript_numbers.json` com todos os números citados no texto.

---

## O que NÃO foi possível verificar sem novo processamento

Todas as conclusões quantitativas da v4 (número de clusters, efeito de proximidade, tendência regional) dependem de rodar a exportação no Earth Engine com a sua conta. Até lá, os únicos números novos citáveis são os calculados com dados reais já disponíveis:
- inventário (A5);
- distâncias reais (A2);
- validação dos estados (A1);
- colocalização dutos × litoral (B2).
