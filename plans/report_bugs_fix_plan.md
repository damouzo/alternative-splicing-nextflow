# Plan de corrección — bugs detectados en los reportes HTML

## Contexto

Se revisó el reporte `MNC_Healthy_vs_MNC_DDX41_Patient` de la ejecución
`batch_2026_04` (CRUK_DDX41_DHX34) y se diagnosticaron 7 problemas, validados
contra los ficheros intermedios y el código fuente. Este plan recoge las
correcciones priorizadas (P0→P3) con el archivo y la localización exacta de
cada cambio, más los pasos de verificación.

Bugs confirmados:

| # | Problema | Gravedad | Origen |
|---|---|---|---|
| 1 | LeafCutter: 0/118.785 clusters "OK" | P0 | Rmd filtra `status == "OK"`, LeafCutter escribe `Success` |
| 2 | Sección QC Metrics vacía | P0 | `nfcore_multiqc_dir` no está en `params.yaml` |
| 3 | ISAR: 0/9.654 switches, q-values todos ~0.99 | P1 | Test satuRn (investigar), probable fallo de modelo |
| 4 | UpSet silencia herramientas sin sets | P1 | `overlap_build` + falta aviso en `overlap_note` |
| 5 | rMATS FDR=0.0 masivo (677 SE, 91 A5SS, 141 A3SS, 250 RI) | P2 | Suelo numérico `--cstat` (display engañoso en volcano/tabla) |
| 6 | PEGASAS min_p=0 en 50/50 vías sin corrección múltiple | P2 | Falta BH sobre ~28.205 tests por vía |
| 7 | HTML 32MB (44 PDFs embebidos en base64) | P3 | Fallback del código cuando `pdftoppm` no está en la imagen |

Datos "limpios" del mismo contraste (no son bugs, solo confusión del lector):
- Los KS scores de PEGASAS SÍ varían por muestra (ADIPOGENESIS: 0.674, 0.704,
  0.661, ...). El heatmap parece idéntico porque el rango real es 0.65-0.70.
- No hay BAMs duplicados ni varianza colapsada en los eventos rMATS FDR=0.
- El contraste MNC tiene 3 sanos + 7 pacientes, BAMs únicos.

---

## P0-1 LeafCutter: `status == "OK"` → `status == "Success"`

**Archivo:** `bin/render_report.Rmd`

LeafCutter (`leafcutter_ds.R`) escribe `status = "Success"` en
`*_cluster_significance.txt`. El Rmd busca `"OK"` en tres sitios, por lo que
`n_ok=0`, `Significant=FALSE` siempre, y la tabla "Top 100" queda vacía.
Se verificó que hay 58.678 clusters `Success` y 6.309 con `p.adjust < 0.05`.

Cambios:

1. Línea ≈955 (`leafcutter_load`):
   ```r
   mutate(Significant = p.adjust < params$fdr_cutoff & status == "OK")
   ```
   → `status == "Success"`.

2. Línea ≈970 (`leafcutter_summary`):
   ```r
   n_ok <- sum(leafcutter_clusters$status == "OK", na.rm = TRUE)
   ```
   → `status == "Success"`.

3. Línea ≈1008 (`leafcutter_top_table`):
   ```r
   top_clusters <- leafcutter_clusters %>% filter(status == "OK")
   ```
   → `filter(status == "Success")`.

4. Añadir comprobación de cordura en `leafcutter_load` (después del filtro de
   `Significant`) para que un futuro mismatch no falle en silencio:
   ```r
   if (nrow(leafcutter_clusters) > 0 &&
       sum(leafcutter_clusters$status %in% c("OK", "Success")) == 0) {
       warning("LeafCutter: ningun cluster con status OK/Success — posible cambio de vocabulario en leafcutter_ds.R")
   }
   ```

**Re-render:** el Rmd es entrada del proceso `RENDER_REPORT`, así que un
`nextflow run ... -resume` tras tocar el Rmd re-renderiza los 10 reportes.

Comando (dir de trabajo: `.../alt_splicing`):
```
bash run_command.sh
```
(pasa `-resume`; el cambio del Rmd invalida solo `RENDER_REPORT`).

**Verificación:** tras el render, comprobar en el HTML:
- `grep -o "Clusters passing noise filter (status = Success)[^<]*"`
- `grep -o "Significant clusters[^<]*"` → debe mostrar >0
- Esperado: ~6.309 clusters significativos; el UpSet debe volver a incluir LeafCutter.

---

## P0-2 Activar QC Metrics (`nfcore_multiqc_dir`)

**Archivo:** `params.yaml` de la ejecución (`.../alt_splicing/params.yaml`).

El Rmd espera que `params$nfcore_multiqc_dir` sea el directorio que CONTIENE
`multiqc_report_data/` (el Rmd hace `file.path(dir, "multiqc_report_data")`).
Ese directorio existe y tiene los ficheros esperados
(`multiqc_general_stats.txt`, `multiqc_star.txt` con `uniquely_mapped_percent`).

Añadir:
```yaml
nfcore_multiqc_dir: /data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/totalRNAseq/batch_2026_04/nfcore_rnaseq/results/multiqc/star_salmon
```

**Nota:** el contraste `MNC_Healthy_vs_MNC_DDX41_Patient` NUNCA pudo mostrar QC
porque `report_params.json` tiene `"nfcore_multiqc_dir": null`. Es global del
run (un solo `params.yaml` para los 10 contrastes), no específico de este
contraste.

**Verificación:**
- Tras re-render: el HTML debe contener la tabla QC (con 10 filas) y los 2
  plots (`Uniquely mapped reads (%)`, `Duplication rate (%)`).
- `grep -c "uniquely_mapped_percent" report.html` > 0.
- Revisar que los sample IDs coincidan con los del samplesheet (la tabla QC
  muestra TODAS las muestras del multiqc, no solo las del contraste).

---

## P1-1 ISAR: q-values pinneados en [0.986, 0.997]

**Estado:** bug confirmado como sistémico (MNC y NB4, potencialmente los 10).
El log del proceso (`ISAR_SWITCH_TEST`) muestra que satuRn corrió dos veces con
0 genes switching, tras un pre-filter que eliminó el 96,37 % de isoformas
(43112 → 2.890 genes / 9.654 isoformas). Hasta isoformas con dIF = 0.63
(U2, ENST00000613956) tienen q = 0.986.

Causa probable: fallo de modelo satuRn (ajuste con muestras desalineadas o
matriz de expresión inconsistente entre el diseño y las cuantificaciones).
Patrón idéntico al bug ya corregido en el commit `64a8231`.

**Pasos de diagnóstico (requiere R + contenedor ISAR):**

1. Inspeccionar el RDS intermedio ya generado
   (`/gpfs/scratch/qp241615/altsplicing-cruk-batch_2026_04/96/baa00ca88388e4878634545b2086de/MNC_Healthy_vs_MNC_DDX41_Patient_tested.rds`)
   con un script R (ver `scripts/diagnose_isar.R` abajo): comprobar
   `isoformSwitchAnalysis`, rango de q-values, y si `dIF` están bien estimados.

2. Verificar el orden muestra↔condición:
   - `isar_import.R` construye `design` y `quant_path` desde el samplesheet CSV.
   - El samplesheet se construye por `ISAR_WRITE_SAMPLESHEET` + `paste` en
     `ISAR_IMPORT`. **Comprobar que `paste -d','` alinea bien las tres
     columnas** (hito del commit 64a8231).
   - Verificar que `importIsoformExpression` devuelve las columnas en el mismo
     orden que `design$sampleID` (mismatch de orden produce q-values uniformes).

3. Re-ejecutar el test a mano con `verbose=TRUE` sobre el RDS raw del contraste
   para ver warnings de convergencia (singular fit / glmmTMB error).

4. Fix probable (según diagnóstico): alinear el orden de `importRdata` con el
   diseño, o usar el modelo correcto en satuRn (verificar en `.opensrc/isar/`
   la versión instalada).

5. Auto-diagnóstico en el reporte (mismo espíritu que LeafCutter): en el chunk
   `isar_summary`, si IQR estrecho de q-values (< 0.05) y hay |dIF| ≥ 0.2,
   emitir aviso visible "posible fallo de convergencia del test".

**Nota:** el fix definitivo de ISAR depende del diagnóstico de saturn/convergencia
y no es un one-liner. Se entrega como tarea de investigación primero, con el
reporte mejorado (aviso automático) como mitigación inmediata.

Script de diagnóstico (`scripts/diagnose_isar.R`):
```r
#!/usr/bin/env Rscript
# Diagnose satuRn q-value pinning (uses the *_tested.rds)
args <- commandArgs(TRUE)
sar <- readRDS(args[1])
feat <- sar$isoformFeatures
q <- feat$isoform_switch_q_value[!is.na(feat$isoform_switch_q_value)]
d <- feat$dIF[!is.na(feat$isoform_switch_q_value)]
cat("n tested:", length(q), "\n")
cat("q range:", range(q), " IQR:", IQR(q), "\n")
cat("max |dIF| non-sig:",
    ifelse(any(q >= 0.05), max(abs(d[q >= 0.05]), na.rm = TRUE), NA), "\n")
# print the design matrix to check ordering
sar$designMatrix
```

---

## P1-2 UpSet: avisar de herramientas activas sin sets

**Archivo:** `bin/render_report.Rmd`

1. Chunk `overlap_build`: ya mantiene `active_tool_sig_sets` por herramienta con
   dir disponible. Añadir, al final del chunk, la lista de sets vacíos:
   ```r
   empty_tool_sets <- setdiff(names(active_tool_sig_sets), names(gene_sets))
   ```

2. Chunk `overlap_note` → convertir a `results='asis'` y añadir:
   ```r
   if (length(empty_tool_sets) > 0) {
     cat(sprintf(
       "<p style='color:#b30000'>ADVERTENCIA: %d herramienta(s) con datos "
       "analizados no aportan genes significativos y quedan fuera del UpSet: %s.</p>\n",
       length(empty_tool_sets), paste(empty_tool_sets, collapse = ", ")
     ))
   }
   ```

3. Chunk `overlap_summary` (línea ≈1201): el texto "Genes significant in all N
   active tools" cuenta `length(active_tool_sig_sets)` pero debe contar
   `length(gene_sets)` para no inflar N cuando una herramienta no aporta genes:
   ```r
   all_tools_genes_keys <- Reduce(intersect, gene_sets)   # gene_sets, no active_tool_sig_sets
   cat(sprintf("Genes significant in all %d tools with signal: %d\n\n",
               length(gene_sets), length(all_tools_genes_keys)))
   ```

**Verificación:** en el reporte re-renderizado del contraste MNC (con LeafCutter
arreglado) debe haber como mucho ISAR fuera del UpSet, y el aviso rojo debe
aparecer solo si hay herramientas sin sets. Con LeafCutter corregido el UpSet
debería volver a tener 3-4 sets.

---

## P2-1 rMATS volcano y tabla: no truncar a 1e-300

**Archivo:** `bin/render_report.Rmd`

Problema: `FDR_plot = pmax(FDR, 1e-300)` crea una banda horizontal en
y≈300 para todos los FDR=0 (que son el 0,2-0,5 % de SE/A5SS/A3SS/RI y son
reales, sólo underflow de `--cstat`), y `formatRound(digits=4)` muestra los
FDR minúsculos como "0.0000" sin contexto.

Cambios:

1. Chunk `rmats_load` (línea ≈276): usar piso más bajo y que la banda se
   pueda identificar. Recomendación: piso `pmax(FDR, .Machine$double.xmin)`
   (~1e-308) y que el volcano recorte la parte superior:
   ```r
   FDR_plot = pmax(FDR, .Machine$double.xmin),
   ```
   Añadir en el chunk `rmats_volcano` (tras el `geom_hline`):
   ```r
   annotate("text", x = -1, y = 300, label = "FDR < 1e-300 (suelo del test, --cstat)", hjust = 0, size = 3, color = "grey40")
   ```

2. Chunk `rmats_top_table`: reemplazar el `formatRound(FDR, digits=4)` por una
   columna texto que distinga 0 exacto:
   ```r
   mutate(FDR_display = ifelse(FDR == 0, "< 1e-300", formatC(FDR, format = "e", digits = 2)))
   ```
   y mostrarla en `datatable` (`esc = FALSE` para el "<").

Opcional pero recomendado (10') — verificación de que `--cstat` no se usa como
clamp de p-value explícito: grep en la fuente de rMATS-turbo
(`.opensrc/rmats/`) por `max(` o `pvalue` junto a `cstat`. En esta instalación
el C no está a la vista, así que basta con confirmar en la doc de `--cstat`
que actúa como margen de equivalencia (H0: |ΔΨ| ≤ c).

**Verificación:** volcanes sin banda de puntos pegada al techo; la tabla top-50
muestra `FDR < 1e-300` para los eventos truncados.

---

## P2-2 PEGASAS: BH por pathway + aclarar los dos tipos de p-value

**Archivo:** `bin/render_report.Rmd` y opcionalmente `bin/cor_matrix_direct_perm.R`

Problema: `pathway_summary` muestra `min_p` (permutación, sin corregir) y con
~28.205 eventos × 50 vías es esperable por azar llegar a p=0 en cada vía.
No obstante los p=0 significan "más extremo que las 1000 permutaciones", que es
correcto a nivel de evento; el problema es presentarlos sin FDR ni contexto.

Cambios en el Rmd (chunk `pegasas_load`):

1. Calcular BH por pathway sobre eventos con p calculado (columna `p_value`):
   ```r
   cor_list <- lapply(cor_list, function(d) {
     d$p_adj <- p.adjust(d$p_value, method = "BH")
     d
   })
   ```

2. `pathway_summary` → añadir columnas `min_p_adj` y `n_sig_adj` (contar
   `p_adj < 0.05`), y reportar esas en vez del `min_p` crudo.

3. `pegasas_summary_table`: mostrar `min_p_adj` y `n_sig_adj`.

4. `pegasas_cor_header`: aclarar que hay DOS p-valores distintos:
   - el **p del enriquecimiento KS** de la vía (en el heatmap, puede llegar a
     -log10(p)~199 legítimamente por el tamaño del gene set), y
   - el **p de correlación evento↔score por permutación** (tablas correlation).
   Añadir un párrafo explícito explicando ambos y que en la tabla de eventos se
   muestra el p de permutación (y ahora su BH).

5. Opcional (mejora de resolución del test en
   `bin/cor_matrix_direct_perm.R`): usar estimador Laplace para evitar p=0
   exacto en la permutación:
   ```r
   p <- (min(mean(obs >= reps, na.rm=TRUE), mean(obs <= reps, na.rm=TRUE)) * N_PERMS + 1) / (N_PERMS + 1)
   ```
   Esto da una resolución inferior de 1/1001 en lugar de 0 y hace la
   corrección BH más significativa. Requiere re-ejecutar `PEGASAS_CORRELATION`
   (solo ese proceso; `-resume` vuelve a calcular las 50×28.205 correlaciones).

**Verificación:** en `pathway_summary` `min_p_adj` > 0 en todas las vías (BH
sobre 28.205 pruebas) y `n_sig_adj` mucho más comedido que `n_significant`.

---

## P3 Contenedor de reporte: asegurar `pdftoppm` (HTML ligero)

**Estado:** el HTML actual (33 MB) contiene 44 PDFs en base64 y 16 PNGs. El
código del Rmd genera PNG cuando `Sys.which("pdftoppm")` encuentra el binario;
la imagen usada para este render NO lo tenía (los PDFs se embebieron como
fallback). El `containers/report/Dockerfile` declara `poppler-utils` (línea 15),
pero hay que confirmar que la imagen que se ejecutó realmente lo trae, o
reconstruirla.

**Pasos:**
1. Verificar la imagen report usada en el render de `batch_2026_04`:
   ```
   apptainer exec /data/BCI-KRP/containers/apptainer_cache/ghcr.io-damouzo-alternative-splicing-nextflow-report-latest.img which pdftoppm
   ```
   (esa imagen en caché es del 1-jun-2026; si no contiene `pdftoppm`, faltan
   `poppler-utils` en la imagen publicada, aunque el Dockerfile lo declare).

2. Si falta: reconstruir el contenedor report
   - Docker: `docker build -t ghcr.io/damouzo/alternative-splicing-nextflow/report:latest containers/report/`
   - Apptainer local: `apptainer build report.sif containers/report/report.def`

3. Re-renderizar con la imagen corregida.

4. Alternativa a largo plazo (más robusta): cambiar el fallback del Rmd para
   que, si `pdftoppm` no está, en vez de incrustar el PDF en base64, publique un
   enlace/thumbnail a los PDFs ya publicados en `outdir/sashimi/...` (así el
   HTML no se vuelve a inflar).

**Verificación objetivo:** tamaño del HTML < 10 MB y
`grep -c "data:application/pdf;base64" report.html` = 0.

---

## Orden de ejecución sugerido

1. **P0-1** LeafCutter fix en Rmd + sanity warning (5').
2. **P0-2** Añadir `nfcore_multiqc_dir` al `params.yaml` (1').
3. Ejecutar `bash run_command.sh` (con `-resume`) → re-render 10 reportes.
   - Comprobar: LeafCutter visible (≈6.309 clusters), QC table presente.
4. **P1-2** UpSet warning + `overlap_summary` count fix (10') → re-render.
5. **P2-1** rMATS volcano/tabla display (20') → re-render.
6. **P2-2** PEGASAS BH + doble p-value note (30') → opcional re-ejecutar
   PEGASAS_CORRELATION si se adopta Laplace.
7. **P0-3/P1-1** ISAR: diagnóstico con `scripts/diagnose_isar.R` sobre los 10
   RDS `*_tested.rds` (cada uno en su work dir) + fix del modelo + aviso
   automático en el reporte.
8. **P3** Verificar/reconstruir imagen report con `pdftoppm` y re-render.

## Paths relevantes

- `params.yaml` ejecución:
  `/data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/totalRNAseq/batch_2026_04/alt_splicing/params.yaml`
- MultiQC:
  `/data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/totalRNAseq/batch_2026_04/nfcore_rnaseq/results/multiqc/star_salmon`
- Salida reportes:
  `/data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/totalRNAseq/batch_2026_04/alt_splicing/results/report/`
- RDS ISAR testeados (work): glob
  `/gpfs/scratch/qp241615/altsplicing-cruk-batch_2026_04/**/*_tested.rds`
- Imagen report en caché local:
  `/data/BCI-KRP/containers/apptainer_cache/ghcr.io-damouzo-alternative-splicing-nextflow-report-latest.img`