# Análise de Serviços Ecossistêmicos de Estoque de Carbono da Vegetação

Este repositório contém o código e a descrição metodológica da etapa analítica da pesquisa **“Serviços Ecossistêmicos de Estoque de Carbono da Vegetação sob Diferentes Intensidades de Uso Antrópico”**.

O objetivo desta etapa é comparar **GPP**, **NPP** e **Biomassa** em relação aos diferentes **usos do solo** identificados por meio de uma **classificação supervisionada** de imagens Sentinel-2.  
Os dados foram cruzados com produtos MODIS em formato TIFF e analisados por meio do script `modoA_analise_stats_v2.py`.


## 1. Descrição geral do código

O script realiza o **cruzamento espacial** entre:

- **Raster de uso e cobertura do solo** (classificação Sentinel-2 / QGIS);
- **Rasters métricos MODIS** (GPP, NPP e Biomassa);
- **Shapefile de municípios**.

### Principais etapas executadas

1. **Recorte e reprojeção** dos rasters métricos para o CRS e grid do raster de classes (uso do solo);
2. **Extração de estatísticas** (média, mediana, desvio padrão, soma e contagem) por classe e município;
3. **Cálculo de totais** (soma × área de pixel);
4. **Testes estatísticos inferenciais** (ANOVA ou Kruskal–Wallis, conforme pressupostos);
5. **Testes pós-hoc** (Tukey HSD ou Dunn–Holm);
6. **Cálculo de tamanhos de efeito** (η² e ε²);
7. **Geração de gráficos** com anotações de p-valor e tamanho de efeito;
8. **Exportação** dos resultados e gráficos por município e métrica.



## 2. Funções analíticas principais

| Função | Descrição |
|--------|------------|
| `_summarize_by_classes` | Calcula média, mediana, std, soma e contagem por classe |
| `_sample_per_class` | Amostragem aleatória estratificada por classe (até 5000 pixels) |
| `shapiro` | Teste de normalidade (por grupo) |
| `levene` | Teste de homogeneidade de variâncias |
| `f_oneway` | ANOVA one-way (caso normalidade + homocedasticidade satisfeitas) |
| `kruskal` | Teste de Kruskal–Wallis (caso não atendidos) |
| `pairwise_tukeyhsd` | Pós-hoc Tukey HSD (para ANOVA) |
| `posthoc_dunn` | Pós-hoc Dunn com ajuste Holm (para Kruskal) |
| `_effect_sizes_anova` | Calcula η² (ANOVA) |
| `_effect_size_kruskal` | Calcula ε² (Kruskal) |
| `_barplot_with_annotation` | Gera gráfico de barras (média ± std) com anotações de p e efeito |
| `WarpedVRT` | Reprojeção dos rasters MODIS para o CRS do raster de classes |
| `gdf.to_crs()` | Reprojeção do shapefile de municípios para o CRS do raster de classes |



## 3. Correção de CRS

O **raster de classes (Sentinel)** é a **referência** do sistema de coordenadas.  
O código **não reprojeta** esse raster internamente — ele deve estar **em CRS projetado (metros)** antes da execução.

As demais correções são feitas automaticamente:
- **Rasters MODIS:** reprojetados com `WarpedVRT` para o CRS/transform do raster de classes;
- **Shapefile de municípios:** reprojetado com `gdf.to_crs(src_class.crs)`.

> 🔸 Se o raster de classes não estiver em CRS projetado, o script lançará um erro.



## 4. Gráficos e visualizações

### Gráficos gerados pelo código
- Gráficos de **barras (média ± desvio padrão)** com legenda contendo:
  - tipo de teste (ANOVA ou Kruskal);
  - p-value global;
  - tamanho de efeito (η² ou ε²).

### Limitações
- Ocultam a **distribuição completa** dos dados (assimetria e outliers);
- O uso do desvio padrão pode ser substituído pelo erro padrão (SEM);
- Podem não representar adequadamente métricas não-normais.

### Gráficos recomendados
| Objetivo | Gráfico sugerido |
|-----------|------------------|
| Comparar classes | **Boxplot** ou **Violin plot** (com pontos jitter) |
| Visualizar correlação entre métricas | **Scatter plot** (com regressão e R²) |
| Comparar métricas múltiplas | **Heatmap** ou **Faceted boxplots** |
| Destacar diferenças significativas | **Boxplot com letras de agrupamento** (a/b/c) do pós-hoc |



## 5. Métodos estatísticos utilizados

### Diagnóstico de pressupostos
- **Normalidade:** `shapiro`
- **Homocedasticidade:** `levene (center='median')`

### Testes principais
- **ANOVA one-way:** `f_oneway`
- **Kruskal–Wallis:** `kruskal`

### Pós-hoc
- **Tukey HSD:** `pairwise_tukeyhsd`  
- **Dunn–Holm:** `scikit_posthocs.posthoc_dunn(method='holm')`

### Tamanhos de efeito
- **η² (eta squared):** medida de magnitude do efeito (ANOVA)  
- **ε² (epsilon squared):** medida aproximada (Kruskal)



## 6. Confiabilidade dos resultados

Os testes aplicados são **estatisticamente adequados** para comparar métricas entre classes, **desde que**:

- Os dados amostrados representem observações independentes;
- As unidades e escalas (GPP, NPP, Biomassa) sejam compatíveis;
- O CRS e a resolução espacial estejam corrigidos.

### Fatores que podem afetar a confiabilidade
1. **Autocorrelação espacial:** pixels próximos não são independentes.  
   - Solução: calcular Moran’s I, agregar por patch ou usar modelos espaciais.
2. **Amostragem por pixel:** pode superestimar o poder do teste.
3. **Unidades:** confirmar unidade original (gC/m², Mg/ha, etc.) antes de interpretar totais.
4. **Testes múltiplos entre cidades:** ajustar p-values (Bonferroni/FDR) se houver comparações múltiplas.



## 7. Tipos de confiabilidade e justificativa dos testes

| Teste / Procedimento | Finalidade | Necessidade |
|-----------------------|-------------|--------------|
| **Shapiro** | Verifica normalidade dos grupos | Necessário (para ANOVA) |
| **Levene** | Verifica homogeneidade das variâncias | Necessário (para ANOVA) |
| **ANOVA** | Compara médias (paramétrico) | Adequado se pressupostos válidos |
| **Kruskal–Wallis** | Compara medianas (não paramétrico) | Alternativa robusta |
| **Tukey / Dunn** | Identifica diferenças par a par | Indispensável para interpretação |
| **η² / ε²** | Mede magnitude do efeito | Fortemente recomendado |



## 8. Recomendações práticas

1. **Garantir CRS projetado** para o raster de classificação (em metros);
2. **Verificar unidades** de GPP/NPP/Biomassa antes de interpretar totais;
3. **Complementar gráficos** de barras com boxplots ou violin plots;
4. **Testar autocorrelação espacial** (Moran’s I);
5. **Agregação por patch** para reduzir dependência espacial;
6. **Modelos hierárquicos (GLMM)** se quiser comparar todas as cidades simultaneamente;
7. **Apresentar IC95% e tamanhos de efeito** nos resultados.



## 9. Extensões sugeridas

- [ ] Substituir `_barplot_with_annotation` por uma função que gere **boxplots + pontos + letras de pós-hoc**  
- [ ] Adicionar cálculo de **Moran’s I** por município  
- [ ] Criar **gráficos de dispersão** (GPP × NPP × Biomassa) com correlação de Pearson e Spearman  
- [ ] Implementar **ajuste para múltiplos testes** entre municípios  
- [ ] Permitir **modelos lineares mistos** (cidade como efeito aleatório)



## 10. Saídas geradas
   - `./dados_gerados/<Cidade>_stats_por_classe.csv`  
   - `./dados_gerados/todas_cidades_stats_por_classe.csv`  
   - `./dados_gerados/stats/resumo_inferencial_por_cidade.csv`  
   - Gráficos PNG por cidade e métrica, com anotação de significância.


## Estrutura Recomendada

Para organização e testes futuros, a estrutura modular ideal seria:

```
geostats_mode_a/
├─ config.py
├─ io.py
├─ sampling.py
├─ infer.py
├─ plots.py
├─ pipeline.py
└─ cli.py
```

> O `main.py` atual já contém toda a lógica consolidada. A modularização é opcional.



## Requisitos

### Linguagem
- Python **3.10+**

### Bibliotecas
- `numpy`
- `pandas`
- `rasterio`
- `geopandas`
- `shapely`
- `affine`
- `scipy`
- `statsmodels`
- `scikit-posthocs`
- `matplotlib`

## Instalação

Instale todas as dependências globalmente para o seu usuário:

```bash
pip install --user numpy pandas rasterio geopandas shapely affine scipy statsmodels scikit-posthocs matplotlib
```



## Como Executar

1. **Ajuste os caminhos no início do arquivo `main.py`:**

   ```python
   class_raster_path = "classificacao/no_clouds2.tif"
   metrics_rasters = {
       "GPP": "metricas/GPP_sete_cidades.tif",
       "NPP": "metricas/NPP_sete_cidades.tif",
       "Biomassa": "metricas/Biomass_sete_cidades.tif",
   }
   vector_cities_path = "shapefile/sete_cidades.shp"
   city_field = "NM_MUN"
   class_map = {1: "Vegetação", 2: "Urbano", 3: "Água", 4: "Solo"}
   ```

2. **(Opcional) Configure os parâmetros globais no topo:**

   * `RESAMPLE_METRICS = "nearest"`
   * `MAKE_PLOTS = True`
   * `RUN_INFERENTIAL_TESTS = True`
   * `EXCLUDE_CLASSES = [5]`
   * `SAMPLE_PER_CLASS = 5000`
   * `ALPHA = 0.05`

3. **Execute o script:**

   ```bash
   python main.py
   ```

4. **Verifique as saídas na pasta `dados_gerados/`.**



## Interpretação dos Resultados

| Tipo de saída                                                                         | Descrição                                                     |
| ------------------------------------------------------------------------------------- | ------------------------------------------------------------- |
| `*_stats_por_classe.csv`                                                              | Estatísticas descritivas por classe e métrica.                |
| `todas_cidades_stats_por_classe.csv`                                                  | Todas as cidades combinadas.                                  |
| `resumo_inferencial_por_cidade.csv`                                                   | Síntese: cidade, métrica, teste global, p, tamanho de efeito. |
| `pairwise_<cidade>_<metrica>_tukey.csv` / `pairwise_<cidade>_<metrica>_dunn_holm.csv` | Comparações par-a-par (pós-hoc).                              |
| Gráficos PNG                                                                          | Médias por classe com caixa de anotação (teste/p/efeito).     |

**Anotação no gráfico:**

```
ANOVA
p = 0.034   η² = 0.62   ★
```

> A estrela `★` indica significância (p < 0.05).
> η² é o tamanho de efeito da ANOVA, ε² o do Kruskal–Wallis.

---

## Dicas de Solução de Problemas

| Problema                                | Causa provável                       | Solução                                  |
| --------------------------------------- | ------------------------------------ | ---------------------------------------- |
| `ValueError: geometries do not overlap` | Polígono fora do raster              | Verifique o CRS e extensão dos dados     |
| CSVs vazios                             | Classe sem pixels válidos            | Confirme `nodata` e máscara              |
| Gráficos sem barras                     | `MAKE_PLOTS=False` ou falta de dados | Ative flag e revise os rasters           |
| ImportError (rasterio/geopandas)        | Falta de libs GDAL/Fiona             | Instale via conda-forge                  |
| P-valores = NaN                         | Classes pequenas (<10 amostras)      | Aumente `MIN_N_FOR_TESTS` ou amplie área |

---

## Referências estatísticas

- **Field, A. (2018)**. *Discovering Statistics Using R.* SAGE.  
- **Zar, J. H. (2010)**. *Biostatistical Analysis.* Pearson.  
- **Legendre & Fortin (1989)**. *Spatial pattern and ecological analysis.* Vegetatio.  
- **Goslee & Urban (2007)**. *The ecodist package for dissimilarity-based analysis of ecological data.* Journal of Statistical Software.



## Script Autor/coautor

**Rayssa de Oliveira Dias e Luiz Felipe Sá**
Projeto: *Serviços Ecossistêmicos de Estoque de Carbono da Vegetação sob Diferentes Intensidades de Uso Antrópico*  


