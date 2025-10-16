#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modo A (v2): análises por cidade + estatística inferencial com ANOTAÇÃO NOS GRÁFICOS
- Mesma lógica do modoA_analise_stats, mas agora cada gráfico recebe uma legenda compacta com:
  • teste global usado (ANOVA ou Kruskal)
  • p-valor global
  • tamanho de efeito (η² para ANOVA, ε² para Kruskal)
- A legenda é inserida no canto superior direito, com caixa clara e discreta.
"""

import os
from typing import Dict, List, Optional
from contextlib import ExitStack

import numpy as np
import pandas as pd

import rasterio
from rasterio.mask import mask as rio_mask
from rasterio.warp import Resampling
from rasterio.vrt import WarpedVRT

import geopandas as gpd
from shapely.geometry import mapping
from affine import Affine

from scipy.stats import shapiro, levene, f_oneway, kruskal
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import scikit_posthocs as sp

import matplotlib.pyplot as plt

# -------------------- Config --------------------
RESAMPLE_METRICS = "nearest"
MAKE_PLOTS = True
OUTDIR = "./dados_gerados"
RUN_INFERENTIAL_TESTS = True
EXCLUDE_CLASSES = [5]
SAMPLE_PER_CLASS = 5000
MIN_N_FOR_TESTS = 10
ALPHA = 0.05

def _pixel_area_from_transform(transform: Affine) -> float:
    return abs(transform.a * transform.e)

def _resampling_mode(name: str):
    name = (name or "nearest").lower()
    return {
        "nearest": Resampling.nearest,
        "bilinear": Resampling.bilinear,
        "cubic": Resampling.cubic,
        "cubicspline": Resampling.cubic_spline,
        "lanczos": Resampling.lanczos,
        "average": Resampling.average,
        "mode": Resampling.mode,
        "max": Resampling.max,
        "min": Resampling.min,
        "med": Resampling.med,
        "q1": Resampling.q1,
        "q3": Resampling.q3,
    }.get(name, Resampling.nearest)

def _summarize_by_classes(values: np.ndarray, classes: np.ndarray) -> pd.DataFrame:
    mask = ~np.isnan(values) & ~np.isnan(classes)
    if not mask.any():
        return pd.DataFrame(columns=["classe", "mean", "median", "std", "sum", "count"])
    v = values[mask]; g = classes[mask].astype(int)
    df = pd.DataFrame({"classe": g, "val": v})
    return df.groupby("classe")["val"].agg(["mean","median","std","sum","count"]).reset_index()

def _barplot_with_annotation(df_city: pd.DataFrame, metric: str, label_col: str, outdir: str, city: str,
                             annot: Optional[dict] = None):
    labels = df_city[label_col].astype(str).values if label_col in df_city.columns else df_city["classe"].astype(str).values
    mean_col = f"{metric}_mean"; std_col = f"{metric}_std"
    if mean_col not in df_city.columns: return
    x = np.arange(len(labels))
    plt.figure()
    plt.bar(x, df_city[mean_col].values, yerr=df_city.get(std_col, None), capsize=4)
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.title(f"{city}: Média de {metric} por classe")
    if annot is not None and "p_global" in annot:
        test_name = annot.get("teste_global","")
        p = annot.get("p_global", np.nan)
        eff = annot.get("efeito", np.nan)
        if isinstance(p, float):
            p_txt = "p < 0.001" if p < 0.001 else f"p = {p:.3f}"
        else:
            p_txt = "p = n/a"
        eff_sym = "η²" if ("ANOVA" in str(test_name)) else "ε²"
        eff_txt = f"{eff_sym} = {eff:.2f}" if isinstance(eff, float) and not np.isnan(eff) else f"{eff_sym} = n/a"
        sig_star = "★" if (isinstance(p, float) and p < 0.05) else ""
        txt = f"{test_name}\n{p_txt}   {eff_txt}  {sig_star}"
        plt.gca().text(0.98, 0.98, txt, transform=plt.gca().transAxes, ha="right", va="top",
                       bbox=dict(boxstyle="round", facecolor="white", alpha=0.7, edgecolor="none"))
    plt.margins(x=0.02); plt.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(os.path.join(outdir, f"{city}_{metric}_means.png"), dpi=200); plt.close()

def _sample_per_class(values: np.ndarray, classes: np.ndarray, exclude: List[int], k: int, min_n: int):
    samples = {}
    mask_all = ~np.isnan(values) & ~np.isnan(classes)
    vals = values[mask_all]; cls = classes[mask_all].astype(int)
    rng = np.random.default_rng(42)
    for c in np.unique(cls):
        if c in exclude: continue
        v = vals[cls == c]
        if v.size >= min_n:
            if v.size > k: v = v[rng.choice(v.size, size=k, replace=False)]
            samples[c] = v
    return samples

def _effect_sizes_anova(groups: List[np.ndarray]) -> float:
    all_vals = np.concatenate(groups); grand = np.mean(all_vals)
    ss_between = sum(len(g)*(np.mean(g)-grand)**2 for g in groups)
    ss_total   = sum(((g - grand)**2).sum() for g in groups)
    return float(ss_between/ss_total) if ss_total>0 else np.nan

def _effect_size_kruskal(groups: List[np.ndarray]) -> float:
    k = len(groups); n = sum(len(g) for g in groups)
    H = kruskal(*groups).statistic
    return float((H - k + 1)/(n - k)) if (n - k) > 0 else np.nan

def _infer_tests_for_city_metric(city: str, metric: str, metr_clip: np.ndarray, class_clip: np.ndarray,
                                 class_map: Optional[dict], outdir: str):
    samples = _sample_per_class(metr_clip, class_clip, EXCLUDE_CLASSES, SAMPLE_PER_CLASS, MIN_N_FOR_TESTS)
    if len(samples) < 2:
        return {"cidade": city, "metrica": metric, "teste_global": "—", "p_global": np.nan, "efeito": np.nan}
    codes = sorted(samples.keys())
    labels = [class_map.get(c, str(c)) if class_map else str(c) for c in codes]
    groups = [samples[c] for c in codes]
    # Diagnósticos
    all_norm = True
    for g in groups:
        if len(g) >= 3:
            try:
                if shapiro(g).pvalue < ALPHA: all_norm = False
            except Exception:
                all_norm = False
        else:
            all_norm = False
    try:
        lev_p = levene(*groups, center='median').pvalue
    except Exception:
        lev_p = np.nan
    homo = (lev_p >= ALPHA) if not np.isnan(lev_p) else False
    use_anova = all_norm and homo
    if use_anova:
        p_global = f_oneway(*groups).pvalue; efeito = _effect_sizes_anova(groups); teste_global = "ANOVA"
        # pós-hoc (salvo mas não necessário para a anotação)
        long_df = pd.DataFrame({"valor": np.concatenate(groups), "classe_nome": np.repeat(labels, [len(g) for g in groups])})
        try:
            tuk = pairwise_tukeyhsd(endog=long_df["valor"].values, groups=long_df["classe_nome"].values, alpha=ALPHA)
            tuk_df = pd.DataFrame(tuk._results_table.data[1:], columns=tuk._results_table.data[0])
            os.makedirs(outdir, exist_ok=True); tuk_df.to_csv(os.path.join(outdir, f"pairwise_{city}_{metric}_tukey.csv"), index=False)
        except Exception:
            pass
    else:
        p_global = kruskal(*groups).pvalue; efeito = _effect_size_kruskal(groups); teste_global = "Kruskal–Wallis"
        long_df = pd.DataFrame({"valor": np.concatenate(groups), "classe_nome": np.repeat(labels, [len(g) for g in groups])})
        try:
            dunn = sp.posthoc_dunn(long_df, val_col="valor", group_col="classe_nome", p_adjust="holm")
            os.makedirs(outdir, exist_ok=True); dunn.to_csv(os.path.join(outdir, f"pairwise_{city}_{metric}_dunn_holm.csv"))
        except Exception:
            pass
    return {"cidade": city, "metrica": metric, "teste_global": teste_global, "p_global": p_global, "efeito": efeito}

def run_mode_A_with_stats_annot(class_raster_path: str, metrics_rasters: Dict[str, str], vector_cities_path: str,
                                city_field: str = "municipio", class_map: Optional[dict] = None,
                                resample_mode: str = RESAMPLE_METRICS, make_plots: bool = MAKE_PLOTS,
                                outdir: str = OUTDIR) -> pd.DataFrame:
    for p in [class_raster_path, *metrics_rasters.values(), vector_cities_path]:
        if not os.path.exists(p): raise FileNotFoundError(f"Caminho não encontrado: {p}")
    os.makedirs(outdir, exist_ok=True); stats_dir = os.path.join(outdir, "stats"); os.makedirs(stats_dir, exist_ok=True)
    metrics_order = list(metrics_rasters.keys())
    with rasterio.open(class_raster_path) as src_class, ExitStack() as stack:
        if not (src_class.crs and src_class.crs.is_projected):
            raise ValueError("O raster de classes precisa estar em CRS projetado (metros).")
        src_metrics = {m: stack.enter_context(rasterio.open(p)) for m, p in metrics_rasters.items()}
        resampling = _resampling_mode(resample_mode)
        gdf = gpd.read_file(vector_cities_path)
        if city_field not in gdf.columns: raise ValueError(f"Campo '{city_field}' não encontrado no vetor.")
        gdf = gdf.to_crs(src_class.crs)
        combined_rows = []; infer_rows = []
        for _, row in gdf.iterrows():
            city = str(row[city_field]).strip()
            if row.geometry is None or row.geometry.is_empty: continue
            geom = [mapping(row.geometry)]
            # class clip
            try:
                class_ma, class_transform = rio_mask(src_class, geom, crop=True, filled=False)
            except ValueError:
                continue
            class_clip = class_ma[0].astype("float32", copy=False)
            if np.ma.isMaskedArray(class_ma): class_clip[class_ma.mask[0]] = np.nan
            if src_class.nodata is not None: class_clip[class_clip == float(src_class.nodata)] = np.nan
            pixel_area_m2 = _pixel_area_from_transform(class_transform)
            metric_stats = []; raw_arrays = {}
            for m, srcm in src_metrics.items():
                SENTINEL = np.float32(-3.4e38)
                with WarpedVRT(srcm, crs=src_class.crs, transform=class_transform,
                               width=class_clip.shape[1], height=class_clip.shape[0],
                               resampling=resampling, src_nodata=srcm.nodata, dst_nodata=float(SENTINEL)) as vrt:
                    metr_ma, _ = rio_mask(vrt, geom, crop=False, filled=False)
                    metr = metr_ma[0].astype("float32", copy=False)
                    if np.ma.isMaskedArray(metr_ma): metr[metr_ma.mask[0]] = np.nan
                    metr[metr == SENTINEL] = np.nan
                    if srcm.nodata is not None: metr[metr == float(srcm.nodata)] = np.nan
                    metr[metr <= -1e10] = np.nan
                stats = _summarize_by_classes(metr, class_clip)
                if not stats.empty:
                    stats["total_kg"] = stats["sum"] * pixel_area_m2
                    stats = stats.rename(columns={
                        "mean": f"{m}_mean", "median": f"{m}_median", "std": f"{m}_std",
                        "sum": f"{m}_sum", "count": f"{m}_count", "total_kg": f"{m}_total_kg"
                    })
                metric_stats.append(stats); raw_arrays[m] = metr
            if not metric_stats: continue
            df_city = metric_stats[0]
            for d in metric_stats[1:]:
                df_city = df_city.merge(d, on="classe", how="outer")
            if df_city.empty: continue
            df_city["cidade"] = city
            if class_map: df_city["classe_nome"] = df_city["classe"].map(class_map).fillna(df_city["classe"].astype(str))
            # salvar por cidade
            city_csv = os.path.join(outdir, f"{city}_stats_por_classe.csv"); df_city.to_csv(city_csv, index=False)
            combined_rows.append(df_city)
            # inferência & anotação
            stats_annots = {}
            if RUN_INFERENTIAL_TESTS:
                city_stats_dir = os.path.join(stats_dir, city.replace(" ", "_")); os.makedirs(city_stats_dir, exist_ok=True)
                for m in metrics_order:
                    res = _infer_tests_for_city_metric(city, m, raw_arrays[m], class_clip, class_map, city_stats_dir)
                    infer_rows.append(res); stats_annots[m] = res
            if MAKE_PLOTS:
                label_col = "classe_nome" if "classe_nome" in df_city.columns else "classe"
                for m in metrics_order:
                    _barplot_with_annotation(df_city, m, label_col, outdir, city, annot=stats_annots.get(m))
        if combined_rows:
            combined = pd.concat(combined_rows, ignore_index=True)
            combined_path = os.path.join(outdir, "todas_cidades_stats_por_classe.csv"); combined.to_csv(combined_path, index=False)
            print(f"[OK] Saída combinada: {combined_path}")
        else:
            combined = pd.DataFrame(); print("[Aviso] Nenhuma estatística gerada.")
        if RUN_INFERENTIAL_TESTS and infer_rows:
            infer_df = pd.DataFrame(infer_rows)
            infer_path = os.path.join(outdir, "stats", "resumo_inferencial_por_cidade.csv")
            infer_df.to_csv(infer_path, index=False); print(f"[OK] Resumo inferencial: {infer_path}")
        return combined

def main():
    # >>>>>>>>>>>> EDITE AQUI SEUS CAMINHOS <<<<<<<<<<<<
    class_raster_path = "classificacao/no_clouds2.tif"
    metrics_rasters = {
        "GPP":      "metricas/GPP_sete_cidades.tif",
        "NPP":      "metricas/NPP_sete_cidades.tif",
        "Biomassa": "metricas/Biomass_sete_cidades.tif",
    }
    vector_cities_path = "shapefile/sete_cidades.shp"
    city_field = "NM_MUN"
    class_map = {1:"Vegetação", 2:"Urbano", 3:"Água", 4:"Solo"}

    run_mode_A_with_stats_annot(class_raster_path, metrics_rasters, vector_cities_path,
                                city_field=city_field, class_map=class_map,
                                resample_mode=RESAMPLE_METRICS, make_plots=MAKE_PLOTS, outdir=OUTDIR)

if __name__ == "__main__":
    main()
