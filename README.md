# Estatística por Cidade com Anotação nos Gráficos

Pipeline para recorte por município, sumarização por classe e testes de hipótese globais (ANOVA ou Kruskal–Wallis) com **anotação automática** (teste, p-valor e tamanho de efeito) diretamente nos gráficos.

---

## 🧩 O que o código faz

1. **Validação e preparo**
   - Confere caminhos e garante que o raster de classes está em **CRS projetado** (metros).
   - Reprojeta o **vetor de cidades** para o CRS do raster de classes.

2. **Recorte por cidade**
   - Recorta o raster de classes pelo polígono da cidade.
   - Para cada **métrica** (ex.: GPP, NPP, Biomassa), utiliza `WarpedVRT` para alinhar ao grid do recorte e extrair os valores válidos.

3. **Sumarização por classe**
   - Calcula `mean`, `median`, `std`, `sum`, `count` e `total_kg` (multiplicando o somatório pela área do pixel em m²).
   - Agrega todas as métricas por classe e gera CSV por cidade.

4. **Estatística inferencial (opcional)**
   - Diagnóstico: normalidade (`Shapiro–Wilk`) e homocedasticidade (`Levene`).
   - **Se normal e homocedástico:** executa **ANOVA** e pós-hoc **Tukey**.
   - **Caso contrário:** executa **Kruskal–Wallis** e pós-hoc **Dunn-Holm**.
   - Salva resultados em CSV e adiciona **caixa de anotação** nos gráficos (teste, p, efeito).

5. **Saídas geradas**
   - `./dados_gerados/<Cidade>_stats_por_classe.csv`  
   - `./dados_gerados/todas_cidades_stats_por_classe.csv`  
   - `./dados_gerados/stats/resumo_inferencial_por_cidade.csv`  
   - Gráficos PNG por cidade e métrica, com anotação de significância.

---

## 📂 Estrutura Recomendada

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

---

## ⚙️ Requisitos

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

### Instalação rápida (via pip)
```bash
pip install numpy pandas rasterio geopandas shapely affine scipy statsmodels scikit-posthocs matplotlib
````

### Instalação recomendada (via conda)

```bash
conda install -c conda-forge python=3.10 rasterio geopandas shapely affine scipy statsmodels scikit-posthocs matplotlib pandas numpy
```

---

## ▶️ Como Executar

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

---

## 📊 Interpretação dos Resultados

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

## 💡 Boas Práticas

* O raster de **classes deve estar em CRS projetado** (área em metros).
* `WarpedVRT` garante alinhamento preciso entre rasters.
* A amostragem é **determinística** (seed fixa) — reproduzível.
* Configure `EXCLUDE_CLASSES` para descartar ruídos ou classes inválidas.
* Pós-hoc são automaticamente salvos em `dados_gerados/stats/<cidade>/`.

---

## 🧠 Dicas de Solução de Problemas

| Problema                                | Causa provável                       | Solução                                  |
| --------------------------------------- | ------------------------------------ | ---------------------------------------- |
| `ValueError: geometries do not overlap` | Polígono fora do raster              | Verifique o CRS e extensão dos dados     |
| CSVs vazios                             | Classe sem pixels válidos            | Confirme `nodata` e máscara              |
| Gráficos sem barras                     | `MAKE_PLOTS=False` ou falta de dados | Ative flag e revise os rasters           |
| ImportError (rasterio/geopandas)        | Falta de libs GDAL/Fiona             | Instale via conda-forge                  |
| P-valores = NaN                         | Classes pequenas (<10 amostras)      | Aumente `MIN_N_FOR_TESTS` ou amplie área |

---
