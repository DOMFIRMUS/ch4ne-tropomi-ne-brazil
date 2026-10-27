# Revisão independente do código v4 (23/09/2026)

**Como foi feita.** Quatro revisores independentes olharam o código, cada um por um ângulo:
- API do Earth Engine;
- estatística;
- geometria;
- reprodutibilidade.

Cada achado passou depois por um verificador cético, encarregado de tentar refutá-lo. Das 33 constatações:
- 21 foram confirmadas pelo verificador e 1 foi refutada;
- 11 ficaram sem verificação automática por limite de uso e foram avaliadas manualmente;
- **todas as 32 procedentes foram corrigidas.**

A revisão usou somente dados reais.

| Área | Problema | Correção |
|---|---|---|
| GEE | Mediana e percentis exatos em redutores separados (risco de estourar a memória) | Um único `percentile([10,50,90])` com `maxRaw` e `parallelScale` |
| GEE | Background de 72 meses recalculado dentro de cada exportação | Exportação em duas fases: tabela de background (A), depois estatísticas (B) |
| GEE | Assets apagados antes de a nova exportação terminar | `overwrite` do próprio `toAsset`, só com `--overwrite` explícito |
| GEE | Órbitas sem a banda corrigida derrubariam a tarefa | Filtro `listContains` das três bandas usadas |
| GEE | Mediana do background ponderada, mas contagem pelo centroide | Redutor `unweighted()` |
| GEE | A checagem `check` computava tudo interativamente | Checagem de um único mês, calculado direto |
| GEE | `getDownloadURL` descarta `crsTransform` | `ee.data.getPixels` na grade do asset, com verificação da grade |
| Estatística | q-valor sem o "+1": clusters do topo recebiam q = 0 | Estimador knockoff+ (Barber & Candès 2015) |
| Estatística | Modelo de dispersão confundia deslocamento de média com dispersão | Resíduos do modelo completo, centrados por anel (Brown–Forsythe) |
| Estatística | Enriquecimento só contra a fração de área | Nulo empírico pelos clusters-espelho + teste exato de Fisher |
| Estatística | IC do pareamento ignorava a dependência dos controles | Bootstrap em blocos dos dois lados + caliper |
| Estatística | Reprodução da v3 com regra de inclusão diferente da do GEE | Emulação da regra ponderada do `reduceRegion` |
| Estatística | Background local incluía a faixa costeira excluída da detecção | Mesma máscara de superfície no background |
| Estatística | Sensibilidade "min_obs = 1" ainda aplicava a regra de meses | Pares explícitos (observações, meses), incluindo (1, 0) |
| Estatística | XCH₄ bruto usado em contrastes inferenciais | Só variáveis sem tendência; `bg_sampled` como diagnóstico |
| Estatística | Blocos de 1° menores que a janela de ~2° | Blocos de 2° + sensibilidade de 1°, 2° e 3° |
| Estatística | Amplitude sazonal parcial; HAC tratava meses ausentes como adjacentes | Amplitude do ciclo completo; Newey–West com defasagem de calendário |
| Estatística | Coluna de polígonos sem nulo correspondente | Removida da Tabela 4 |
| Geometria | Comprimentos na Policônica ~2% maiores | Comprimentos geodésicos: 9.880 km em terra no NE |
| Geometria | Áreas de célula na esfera (+0,4%) | Áreas no elipsoide WGS84 (idênticas ao pyproj) |
| Geometria | Fernando de Noronha alargava os mapas em ~15% | Enquadramento pelo continente |
| Reprodutibilidade | Parâmetros da exportação não persistiam entre etapas | `run_config.json`, identidade carimbada nos assets e checada no download e na análise |
| Reprodutibilidade | Versões de software e parâmetros não registrados | Bloco `run` em `manuscript_numbers.json` + `requirements-lock.txt` |
| Reprodutibilidade | Saídas antigas misturadas às novas | Limpeza das saídas da rodada antes de gravar |
| Reprodutibilidade | Arquivos de nomes de bandas fora do checksum | Incluídos no `PROVENANCE.json` |
| Reprodutibilidade | Proveniência e versões do processador fora do Git | `.gitignore` versiona `PROVENANCE.json` e o background mensal |
| Reprodutibilidade | Projeto do GEE fixo no código | `--project` e variável `CH4NE_GEE_PROJECT` |
| Reprodutibilidade | Períodos P1/P2 incompletos | `--start/--end/--label` (rótulo obrigatório) e etapa `compare` |
| Reprodutibilidade | Notebook do Colab dependia da ordem das células; saídas se perdiam na VM | Células independentes; saídas no Drive |
