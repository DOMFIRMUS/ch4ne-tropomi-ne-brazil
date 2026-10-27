# Triagem de anomalias de metano (TROPOMI) sobre a infraestrutura de óleo e gás do Nordeste do Brasil — pipeline v4

Código (versão 4.1.0) do artigo *Detectability of pipeline-scale methane enhancements from TROPOMI under false-positive control: a screening protocol for onshore oil and gas networks in Northeast Brazil* (Cardoso et al., submetido em 2026; 25 autores, listados com ORCID e afiliação em [CITATION.cff](CITATION.cff)).

Esta é a versão em português. A versão principal, em inglês, é o [README.md](README.md).

Esta versão (v4) substitui os notebooks v2/v3. A auditoria que motivou a reescrita está em [docs/AUDITORIA_v3.md](docs/AUDITORIA_v3.md).

> **Política de dados:** este repositório não contém nem gera dados sintéticos ou inventados. A análise só roda sobre rasters baixados do Google Earth Engine, com registro de proveniência (`PROVENANCE.json`, com checksums). Os testes usam apenas dados reais: camadas da ANP, limites estaduais e as tabelas publicadas da v3.

## O que mudou em relação à v3

| v3 | v4 |
|---|---|
| Nomes dos estados reconstruídos por índice (6/9 trocados) | Estado vem sempre do atributo do polígono, com checagem independente ponto-em-polígono |
| Distâncias cluster–duto sempre 0 (lidas do objeto errado) | Distâncias geodésicas reais (UTM SIRGAS 2000 por zona / azimutal equidistante por cluster) |
| Média de 6 anos sem remover a tendência | Anomalia por órbita = XCH₄ − mediana regional do mês |
| Z contra média e desvio regionais únicos | Background local (mediana móvel ~220 km) + z robusto (mediana/MAD) |
| Sem número mínimo de observações | ≥ 50 recuperações e ≥ 24 meses com dados por célula |
| `reproject` (amostragem) a 7 km | `reduceResolution` (média de área) numa grade fixa de 0,0625° |
| "Limpeza morfológica" que não erodia, a 10 km | Componentes conexos na grade nativa + q-valor de falsa descoberta (teste do espelho, knockoff+) |
| Zonas aninhadas; desvio espacial lido como temporal | Anéis exclusivos; variância temporal real; frequência de desvios + e −; controles pareados (litoral, albedo, água, altitude); bootstrap espacial em blocos; regressões com erro robusto por bloco |
| Dutos com buffer de 1° (inclui ES e offshore) | Inventário classificado (onshore NE / offshore / outro estado); distâncias só às partes em terra |
| Figuras de thumbnails com título embutido | Figuras a partir dos GeoTIFFs, 600 dpi + PDF, escala geodésica, seta norte, grade °S/°W, paletas acessíveis |

## Estrutura

```
ch4ne/                  pacote Python
  config.py             TODOS os parâmetros (o Método pode ser escrito a partir daqui)
  gee_export.py         etapa Earth Engine: órbitas -> anomalia -> estatísticas por célula -> assets -> download
  provenance.py         registro e verificação da proveniência dos rasters
  boundaries.py         limites estaduais (atributos UF/nome) e litoral
  infra.py              inventário ANP harmonizado e classificado
  raster.py             leitura dos GeoTIFFs, rasterização, distâncias
  detection.py          background local, z robusto, clusters, q-valores, atribuição
  zones.py              anéis de distância, bootstrap em blocos, pareamento, regressões
  states.py             tabela por estado + checagem independente
  timeseries.py         tendência regional e séries dos clusters
  figures.py            figuras de publicação
  analysis.py           orquestra a análise local e grava tabelas/figuras/números
scripts/run_pipeline.py interface de linha de comando
notebooks/              notebook para Google Colab
tests/                  testes com dados reais (pytest)
data/anp/               camadas públicas da ANP (dutos) usadas no artigo
data/limites/           limites estaduais (cópia fixa; checksums em data/SHA256SUMS.txt)
validacao_offline/      reprodução dos erros da v3 a partir das tabelas publicadas
validacao_v4_1/         conferências da v4.1 com dados reais (versão nova × antiga, área pequena)
docs/                   auditoria, referências verificadas, recomendações de método, pendências, resultados (RESULTADOS_v4_1.md)
outputs_v4/             saídas (tabelas, figuras, números); rasters em outputs_v4/gee/ (fora do Git)
```

## Como rodar

### 1. Ambiente

```bash
python -m pip install -r requirements.txt
```

### 2. Autenticação no Earth Engine (uma vez por máquina)

```bash
python -c "import ee; ee.Authenticate(auth_mode='localhost', force=True)"
```

O navegador abre. Entre com a conta Google do Earth Engine e autorize. O projeto do Google Cloud padrão está em `ch4ne/config.py` e pode ser trocado de duas formas:
- com `--project <id>` em qualquer etapa;
- com a variável de ambiente `CH4NE_GEE_PROJECT`.

O projeto precisa estar registrado em https://code.earthengine.google.com/register.

### 3. Processamento (na raiz do repositório)

```bash
python scripts/run_pipeline.py check
```

```bash
python scripts/run_pipeline.py all
```

- `check` (cerca de 1 minuto): bandas, órbitas, versões do processador, uma órbita preparada e o background de um mês.
- `all`: exportação em duas fases. A **fase A** grava os backgrounds mensais num asset (tabela). A **fase B** usa esse asset e grava os assets de imagem:
  - `stats`: todas as estatísticas por célula, numa única passagem pelas órbitas;
  - `monthly_2020` … `monthly_2025`: anomalia e contagem mensais, um asset por ano;
  - `legacy`: o composto da v3.

  Depois vêm `wait`, `download` e `analyze`. O `download` calcula as covariáveis estáticas (albedo MODIS, relevo SRTM e água JRC) direto na grade, em blocos, com `ee.data.computePixels` (`covariates.tif`). Em seguida monta `stats.tif`, `monthly_anom.tif` e `monthly_count.tif`, além do CSV e do `PROVENANCE.json`.
- A exportação pode ser retomada. Um asset que já existe com a mesma identidade da rodada é reaproveitado, e uma tarefa ainda em execução é aguardada; rodar de novo nunca refaz trabalho concluído. `--jobs stats monthly` exporta só os assets indicados.
- Assets existentes nunca são apagados sem `--overwrite`. Mesmo com a opção, o asset antigo só é substituído quando a nova tarefa conclui.
- Se só as séries mensais falharem, a análise roda assim mesmo e omite apenas as séries temporais dos clusters. O número de meses com dado já vem de `stats`.
- As etapas podem ser rodadas separadamente. A etapa `export` grava `run_config.json`, e as etapas seguintes param se a linha de comando pedir outros parâmetros.
- **Cota (Community Tier).** O custo de cada tarefa fica registrado em `outputs_v4*/gee/gee_tasks.csv`. Na rodada principal (72 meses), as tarefas que deram certo custaram 8,2 EECU-h: `stats` 4,6 e os seis `monthly_*` juntos 3,5. P1 custou 3,7 EECU-h e P2, 8,7; o custo de `stats` varia bastante entre execuções. Os parâmetros `*_tile_px` em `config.py` controlam só a memória, nunca os resultados. Falhas medidas: como exportação em lote, a covariável de água do JRC (30 m) ultrapassou o limite de 2^31 pixels, mesmo com blocos de 32 px (por isso ela é calculada com `computePixels`); e um único asset mensal de 144 bandas ficou sem memória.

**Rodadas de robustez** por versão do processador (o rótulo é obrigatório para datas diferentes do padrão):

```bash
python scripts/run_pipeline.py all --start 2020-01-01 --end 2022-07-16 --label P1 --min-obs 22
```

```bash
python scripts/run_pipeline.py all --start 2022-07-17 --end 2025-12-31 --label P2 --min-obs 29
```

```bash
python scripts/run_pipeline.py compare --compare P1 P2
```

`--min-obs` mantém nos períodos curtos a mesma densidade da regra principal (50 recuperações em 72 meses): P1 tem 31 meses e P2, 42. Com o mínimo absoluto de 50, só 6 células perto de dutos sobram em P1.

### 4. Testes

```bash
python -m pytest tests -q
```

### Google Colab

Abra `notebooks/CH4_NE_v4_colab.ipynb`. Ele clona ou monta este repositório, instala as dependências, autentica e roda as mesmas etapas.

## Saídas (`outputs_v4/`)

| Arquivo | Conteúdo |
|---|---|
| `tables/Table1_infrastructure_inventory.csv` | segmentos ANP por classe e localização |
| `tables/Table2_states_strict.csv` | estatísticas por estado (triagem estrita) |
| `tables/Table3_clusters.csv` | clusters: pico, centroide, estado, z, q-valor, distâncias, ativo mais próximo |
| `tables/Table4_attribution_enrichment.csv` | observado × esperado de picos a 5/10/20 km de dutos |
| `tables/Table5_proximity_rings_*.csv` | estatísticas por anel de distância |
| `tables/Table6_near_vs_far_bootstrap_*.csv` | contrastes perto × longe, IC 95% por bootstrap em blocos |
| `tables/Table7_near_vs_matched_controls.csv` | contrastes com controles pareados |
| `tables/Table8_regressions_cluster_robust.csv` | regressões (média e dispersão) com erro robusto por bloco |
| `tables/TableS*.csv` | inventário por segmento, reprodução da v3, estados (triagem relaxada), balanço do pareamento, background mensal, séries dos clusters, sensibilidades (limiares, janela, tamanho de bloco), colocalização dutos × litoral, comparação P1 × P2 |
| `figures/Fig0*.png/.pdf` | figuras do artigo (sem título; legendas ficam no texto) |
| `manuscript_numbers.json` | todos os números citados no texto, com a proveniência, os parâmetros e as versões de software da rodada |
| `checks.json` | resultado das checagens automáticas |
| `clusters_v4.geojson`, `derived_fields.tif` | camadas para SIG |

## Dados

- **TROPOMI/Sentinel-5P:** coleção `COPERNICUS/S5P/OFFL/L3_CH4` do Google Earth Engine, banda `CH4_column_volume_mixing_ratio_dry_air_bias_corrected`. O GEE mantém apenas qa_value > 0,5.
- **Dutos:** camadas públicas da ANP (gasodutos de escoamento, transporte e distribuição; oleodutos).
- **Limites estaduais:** `click_that_hood/brazil-states.geojson`, na mesma versão usada na v3 (cópia fixa com checksum). Pode ser trocado pela malha oficial do IBGE em `config.states_file`.
- **Covariáveis** (no GEE): MODIS MCD43A3 (albedo, banda 7), CGIAR SRTM90 (altitude), JRC Global Surface Water (água).

## Como citar

Cite o artigo (referência acima; será atualizada na publicação) e o arquivo no Zenodo, versão 2.0.0, que reúne este código e os dados derivados das três rodadas. DOI-conceito (todas as versões): https://doi.org/10.5281/zenodo.18717651. A versão 1.0.0 (notebooks v3) foi substituída e **não deve ser usada**. O GitHub mostra a citação pronta em "Cite this repository", a partir do [CITATION.cff](CITATION.cff).

## Licença

- **Código:** licença MIT (ver [LICENSE](LICENSE)).
- **Dados derivados e figuras** (`outputs_v4*/` e o arquivo no Zenodo): CC BY 4.0. Contêm dados Copernicus Sentinel modificados (2020–2025). As camadas da ANP em `data/anp/` são dados públicos da ANP e seguem os termos da ANP; cite a fonte ao reutilizá-las.
