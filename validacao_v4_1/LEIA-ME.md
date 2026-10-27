# Conferências da v4.1 com dados reais (30/09/2026)

Antes de gastar cota em lote no Earth Engine, cada mudança de desempenho da v4.1 foi comparada com a versão anterior. A comparação usou dados reais do TROPOMI numa área pequena, em chamadas interativas. Os scripts e as saídas estão aqui. Para rodar, é preciso o Earth Engine autenticado; o projeto é o de `ch4ne/config.py`.

| Script | O que compara | Resultado (log) |
|---|---|---|
| `verifica_prep_orbita.py` | Agregação de cada órbita com um só `reduceResolution` (XCH4, AOD e fração válida) contra a versão com dois, em 4 órbitas de set/2023 com dados (caixa de 3° × 3° em BA/PE) | Mesmas células, diferença 0 (`verifica_prep_orbita.log`) |
| `verifica_passo_unico_v2.py` | As **17 bandas** de `stats_image` (uma passagem, redutor de múltiplas entradas) contra uma referência independente, com redutores de uma entrada em passagens separadas e `n_months` tirado da pilha mensal (caixa de 1° × 1° na BA, 2023, 255 células) | Diferença 0 em todas as bandas. `countDistinct` comum − pilha = **+1 em 255 de 255 células**; `countDistinctNonNull` − pilha = 0 em 255 de 255 (`verifica_passo_unico_v2.log`) |
| `verifica_mensal_passo_unico.py` | Pilha mensal numa passagem contra 12 subcoleções filtradas (2023) | Diferença 0 nas 24 bandas; células válidas idênticas em todos os meses (`verifica_mensal_passo_unico.log`) |
| `verifica_passo_unico.py`, `verifica_n_months2.py` | Primeiras versões das conferências acima | Ver a nota abaixo |
| `diag_covariaveis.py`, `teste_computepixels.py` | Covariáveis numa caixa pequena e tempo de `computePixels` | Valores plausíveis. Na grade inteira, 88 blocos de 32 × 32 células levaram 418 s com 4 requisições em paralelo, ≈ 19 s por bloco (`../outputs_v4/logs/covariaveis_computepixels.log`) |

**Primeiras versões.**
- `verifica_passo_unico.log` comparou só 9 das 17 bandas. O teste de órbitas desse script caiu em órbitas sem dados na caixa, então não testou nada; a conferência de órbitas válida é `verifica_prep_orbita.log`.
- Esse script também usava a primeira versão do redutor, com `countDistinct`, daí `n_months_minus_stack_maxabs = 1`.
- `verifica_n_months2.log` já é da versão corrigida.

A versão v2 substitui as duas, com todas as bandas e o histograma célula a célula. A checagem `n_months_equals_monthly_stack` do `checks.json` confirma o resultado na grade inteira das três rodadas.
