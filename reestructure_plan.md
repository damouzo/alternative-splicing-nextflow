# Plan de Reestructura de Outputs

## 0) Objetivo real del plan

Este plan define como reorganizar la publicacion de resultados para que:

1. El equipo wet lab e IP tenga una capa simple y util de consulta.
2. El equipo bioinfo tenga tablas maestras completas para downstream analysis.
3. La pipeline mantenga trazabilidad y reproducibilidad sin ruido innecesario.
4. No se rompa la estructura actual ni los reportes existentes.

Contexto real observado:

- La separacion por herramienta y por contraste ya es buena.
- El problema aparece en el ultimo nivel: demasiados archivos tecnicos, poca capa de consumo.
- Caso ejemplo claro en rMATS:
  /data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/totalRNAseq/batch_2026_04/alt_splicing/results/rmats/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient

Decisiones cerradas para la v1 de la reestructura:

1. Se mantiene la estructura actual por tool y contraste.
2. Se anade una capa oficial de consumo en results/deliverables.
3. Se generan dos tablas nuevas obligatorias en rMATS: master y significant.
4. Para MAJIQ, ISAR, LeafCutter y PEGASAS no se rehace estadistica: se estandariza export y naming.
5. WORK no se usa como almacenamiento de archivo a medio plazo.
6. Lo que no va en CORE pero puede ser util en auditoria se publica en RAW.

---

## 1) Reglas de buenas practicas (primero reglas, luego implementacion)

### Regla 1: separar salida de consumo vs salida tecnica

- Capa de consumo: lo que se consulta de forma habitual.
- Capa tecnica: lo que sirve para auditoria, debug o preguntas puntuales.

### Regla 2: un punto de entrada por contraste

Cada contraste debe tener una carpeta de deliverables con:

1. tablas maestras
2. tablas de significativos
3. plots indexados
4. metadata del contraste

### Regla 3: publicar poco pero con alta utilidad

No publicar por defecto todo lo que produce cada herramienta.
Publicar lo que responde preguntas biologicas y de priorizacion.

### Regla 4: no usar work como archivo historico

work es util para ejecucion y debug inmediato, pero no es una capa de publicacion fiable a medio plazo.
Si un artefacto puede ser necesario en 3-6 meses, debe ir a results/raw y no depender de work.

### Regla 5: trazabilidad minima obligatoria

Toda salida de consumo debe poder responder:

1. que muestras entraron
2. que parametros clave se usaron
3. que umbrales definieron significancia
4. de que archivo original se deriva cada fila

### Regla 6: mantener compatibilidad hacia atras

No romper rutas legacy en esta fase. La reestructura es aditiva.

### Regla 7: nombres canonicos estables

Para todas las herramientas:

- <comparison_id>.<tool>.master.tsv
- <comparison_id>.<tool>.significant.tsv
- <comparison_id>.<tool>.summary.tsv

### Regla 8: preservar metrica nativa de herramienta

La normalizacion de columnas no debe ocultar columnas nativas.

### Regla 9: QA automatica antes de compartir

Ninguna corrida se comparte si falla el validador de estructura y consistencia.

### Regla 10: report HTML como guia, no como contenedor unico

El report debe apuntar a archivos en disco, no sustituirlos.

---

## 2) Modelo de publicacion recomendado

Tres niveles claros, activables por parametro:

1. CORE
Para compartir por defecto con IP/wet lab/bioinfo.

2. RAW
Para conservar evidencia tecnica relevante sin llenar la capa CORE.

3. WORK_ONLY
Artefactos que no se publican. Permanecen solo en work durante la ejecucion.

Decision operativa:

- CORE: obligatorio
- RAW: recomendado
- WORK_ONLY: solo para intermedios de bajo valor

---

## 3) Estructura objetivo en results

La estructura actual por herramienta se mantiene, pero se anade una capa oficial de entrega:

results/
  rmats/
  majiq/
  isoformswitchr/
  leafcutter/
  sashimi/
  pegasas/
  report/
  deliverables/
    metadata/
      run_manifest.yaml
      sample_index.tsv
      tools_matrix.tsv
    contrasts/
      <comparison_id>/
        data_tables/
          <comparison_id>.rmats.master.tsv
          <comparison_id>.rmats.significant.tsv
          <comparison_id>.majiq.master.tsv
          <comparison_id>.majiq.significant.tsv
          <comparison_id>.isar.master.tsv
          <comparison_id>.isar.significant.tsv
          <comparison_id>.leafcutter.master.tsv
          <comparison_id>.leafcutter.significant.tsv
          <comparison_id>.pegasas.master.tsv
          <comparison_id>.pegasas.significant.tsv
          <comparison_id>.cross_tool.master.tsv
        plots/
          sashimi/
            sashimi_index.tsv
        summaries/
          <comparison_id>.summary.tsv
        metadata/
          contrast_manifest.yaml

---

## 4) Ejemplo completo: que hacer con rMATS

### 4.1 Problema observado

En cada contraste rMATS hay una buena cantidad de archivos utiles para computo, pero poco operativos para IP/wet lab.

### 4.2 Salidas nuevas obligatorias en CORE

1. <comparison_id>.rmats.master.tsv
Incluye TODOS los eventos (significativos y no significativos) de SE/A5SS/A3SS/MXE/RI.

2. <comparison_id>.rmats.significant.tsv
Solo eventos que cumplen criterios de significancia.

3. <comparison_id>.rmats.summary.tsv
Resumen por tipo de evento:
- n_total
- n_significant
- n_novel_splice_site
- mediana abs(deltaPSI)
- percentil 90 abs(deltaPSI)

### 4.3 Columnas recomendadas del master rMATS

Minimo:

1. comparison_id
2. tool (rmats)
3. event_type
4. event_id
5. gene_id
6. gene_symbol
7. chr
8. strand
9. inc_level_1
10. inc_level_2
11. inc_level_difference
12. pvalue
13. fdr
14. is_novel_splice_site
15. is_significant
16. significance_rule
17. source_file

### 4.4 Regla de significancia por defecto

is_significant = (FDR <= report_fdr_cutoff) and (abs(IncLevelDifference) >= report_dpsi_cutoff)

### 4.5 Que queda en RAW para rMATS

1. SE.MATS.JC.txt, A5SS.MATS.JC.txt, A3SS.MATS.JC.txt, MXE.MATS.JC.txt, RI.MATS.JC.txt
2. fromGTF.*.txt
3. fromGTF.novelJunction.*.txt
4. fromGTF.novelSpliceSite.*.txt
5. b1_samples.txt, b2_samples.txt

### 4.6 Que puede quedar en WORK_ONLY (no publicar)

1. tmp/
2. JC.raw.input.*
3. JCEC.raw.input.*
4. individualCounts.*
5. *.MATS.JCEC.txt (opcional RAW si se quiere conservar)

Nota bioinfo:
Para cohortes de splicing factors mutados, novelSpliceSite merece conservarse en RAW casi siempre.

---

## 5) Matriz por herramienta: CORE vs RAW vs WORK_ONLY

## 5.1 rMATS

CORE:
1. rmats.master.tsv
2. rmats.significant.tsv
3. rmats.summary.tsv

RAW:
1. MATS.JC por tipo
2. fromGTF.*
3. novelJunction/novelSpliceSite
4. b1/b2 samples

WORK_ONLY:
1. tmp
2. raw input counts
3. JCEC e individualCounts si no se necesita auditoria fina

## 5.2 MAJIQ

CORE:
1. majiq.master.tsv
2. majiq.significant.tsv
3. majiq.summary.tsv

RAW:
1. <comparison>.tsv original completo (si no coincide al 100% con master)
2. built_sg.zarr
3. <comparison>.dpsicov

WORK_ONLY:
1. temporales de ejecucion y conversion

## 5.3 IsoformSwitchAnalyzeR

CORE:
1. isar.master.tsv
2. isar.significant.tsv
3. consequence_summary.csv (si existe)
4. isar.summary.tsv

RAW:
1. <comparison>_final.rds
2. switchplots (si interes para figuras)

WORK_ONLY:
1. intermedios de import y anotacion no consultados

## 5.4 LeafCutter

CORE:
1. leafcutter.master.tsv
2. leafcutter.significant.tsv
3. leafcutter.summary.tsv

RAW:
1. cluster_significance.txt
2. effect_sizes.txt

WORK_ONLY:
1. intermedios de bam2junc y clustering no usados en interpretacion

## 5.5 PEGASAS

CORE:
1. pegasas.master.tsv
2. pegasas.significant.tsv
3. pegasas.summary.tsv

RAW:
1. all_pathways_scores.tsv
2. correlation_out por pathway
3. sig_pathways por contraste

WORK_ONLY:
1. temporales de agregacion

## 5.6 Sashimi

CORE:
1. sashimi_index.tsv
2. PDFs finales de plots

RAW:
1. directorio completo sashimi_out por tipo de evento

WORK_ONLY:
1. archivos auxiliares internos de plotting duplicados

## 5.7 Report

CORE:
1. <comparison>_splicing_report.html
2. seccion Data exports con rutas a deliverables

RAW:
1. assets auxiliares si alguna version futura los separa

WORK_ONLY:
1. objetos temporales del render

## 5.8 Ejemplos reales por herramienta (mismo contraste)

Contraste de referencia:
/data/BCI-KRP/projects/CRUK_DDX41_DHX34/analysis/totalRNAseq/batch_2026_04/alt_splicing/results/*/MNC_Healthy_vs_MNC_DDX41_Patient

rMATS observado:
1. results/rmats/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/SE.MATS.JC.txt
2. results/rmats/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/individualCounts.SE.txt
3. results/rmats/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/tmp/

MAJIQ observado:
1. results/majiq/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient.tsv
2. results/majiq/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/built_sg.zarr/
3. results/majiq/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient.dpsicov/

IsoformSwitchAnalyzeR observado:
1. results/isoformswitchr/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/top_isoform_switches.csv
2. results/isoformswitchr/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/consequence_summary.csv
3. results/isoformswitchr/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient_final.rds

LeafCutter observado:
1. results/leafcutter/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient_cluster_significance.txt
2. results/leafcutter/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient_effect_sizes.txt

PEGASAS observado:
1. results/pegasas/MNC_Healthy_vs_MNC_DDX41_Patient/pegasas_out/all_pathways_scores.tsv
2. results/pegasas/MNC_Healthy_vs_MNC_DDX41_Patient/pegasas_out/correlation_out/
3. results/pegasas/MNC_Healthy_vs_MNC_DDX41_Patient/MNC_Healthy_vs_MNC_DDX41_Patient_sig_pathways.tsv

Sashimi observado:
1. results/sashimi/MNC_Healthy_vs_MNC_DDX41_Patient/sashimi_out/SE/Sashimi_plot/*.pdf
2. results/sashimi/MNC_Healthy_vs_MNC_DDX41_Patient/sashimi_out/SE/Sashimi_index/
3. results/sashimi/MNC_Healthy_vs_MNC_DDX41_Patient/sashimi_out/SE/Sashimi_index_<GENE>_<RANK>/

Report observado:
1. results/report/MNC_Healthy_vs_MNC_DDX41_Patient_splicing_report.html

Interpretacion operativa:
1. Las rutas actuales son tecnicamente correctas.
2. Lo que falta no es estructura por carpetas, sino una capa de consumo canonica y transversal.
3. Por eso la reestructura se centra en deliverables + clasificacion CORE/RAW/WORK_ONLY.

---

## 6) Que tablas extra hay que generar (ademas de rMATS)

No solo rMATS. Para toda la pipeline:

1. Tabla master por herramienta y contraste.
2. Tabla significant por herramienta y contraste.
3. Tabla summary por herramienta y contraste.
4. Tabla cross_tool por contraste (opcional pero muy util).

### 6.1 Estandar minimo de columnas en todos los masters

1. comparison_id
2. tool
3. feature_type
4. feature_id
5. gene_id
6. gene_symbol
7. effect_size
8. pvalue (si aplica)
9. fdr o metrica equivalente
10. is_significant
11. source_file

### 6.2 Cross-tool master recomendado

Unifica filas relevantes de todas las herramientas para consumo rapido del equipo bioinfo.
No sustituye las tablas nativas de cada tool.

---

## 7) Metadata obligatoria de publicacion

## 7.1 run_manifest.yaml (global)

Campos minimos:

1. pipeline_name
2. pipeline_version
3. nextflow_version
4. session_id
5. execution_date
6. outdir
7. params_criticos
8. tools_enabled
9. tools_disabled
10. comparisons
11. results_contract_version

## 7.2 contrast_manifest.yaml (por contraste)

Campos minimos:

1. comparison_id
2. group1_name
3. group2_name
4. group1_n
5. group2_n
6. report_html_path
7. data_tables
8. plots
9. thresholds_used

## 7.3 sample_index.tsv (global)

Debe permitir responder de forma directa:

- que muestras entraron en cada contraste
- a que grupo pertenecen

Campos:

1. sample_id
2. condition
3. replicate
4. comparison_id
5. group
6. bam_path
7. bai_path
8. salmon_dir

---

## 8) Cambios concretos en pipeline y report

## 8.1 Cambios de pipeline (sin romper)

Archivos principales a tocar:

1. workflows/alternative_splicing.nf
2. nextflow.config
3. modules/local/report/main.nf
4. bin/render_report.Rmd

Modulos/scripts nuevos recomendados:

1. modules/local/export_metadata/main.nf
2. modules/local/rmats_master/main.nf
3. modules/local/export_tool_masters/main.nf
4. modules/local/sashimi_index/main.nf
5. bin/build_rmats_master.py
6. bin/export_tool_masters.py
7. bin/build_sashimi_index.py
8. bin/build_results_manifest.py
9. bin/validate_results_contract.py

## 8.2 Cambios de report

Cambios minimos y utiles:

1. anadir seccion Data exports
2. mostrar rutas a:
- master table
- significant table
- sashimi index
- manifest del contraste
3. mantener top-X y visualizacion actual sin cambiar semantica

---

## 9) Parametros nuevos recomendados en nextflow.config

1. publish_level = 'core' | 'core_raw' | 'full'
2. publish_deliverables = true
3. build_cross_tool_master = true
4. publish_rmats_jcec = false
5. publish_rmats_individual_counts = false
6. publish_raw = true

Regla de decision por publish_level:

1. core: publica solo deliverables CORE + report HTML.
2. core_raw: publica CORE + RAW.
3. full: publica CORE + RAW + artefactos ahora marcados como WORK_ONLY.

Recomendacion por defecto para colaboracion de grupo:

1. publish_level = core_raw
2. publish_rmats_jcec = false
3. publish_rmats_individual_counts = false

Interpretacion:

- core: publica solo capa CORE.
- core_raw: publica CORE + RAW.
- full: publica todo (incluyendo artefactos menos utiles).

---

## 10) Plan de implementacion por fases

## Fase A - Base de contrato y metadata (bajo riesgo)

1. Crear deliverables/metadata.
2. Generar run_manifest + sample_index + tools_matrix.
3. Integrar sin afectar ramas opcionales.

Criterio de aceptacion:
- se genera metadata en todas las corridas.

## Fase B - Masters de tools con tabla ya amplia (bajo riesgo)

1. MAJIQ, ISAR, LeafCutter, PEGASAS.
2. Estandarizar naming y significant/summary.

Criterio:
- 1 master + 1 significant + 1 summary por tool activa.

## Fase C - rMATS master y summary (riesgo medio)

1. Consolidar los 5 tipos de evento.
2. Marcar novel splice sites.
3. Generar significant y summary.

Criterio:
- filas del master = suma de filas de los 5 MATS.JC.

## Fase D - Sashimi index y limpieza de publicacion (riesgo medio)

1. indice navegable con rutas reales.
2. evitar duplicados de poco valor en capa CORE.

Criterio:
- cada fila del indice apunta a un PDF existente.

## Fase E - Report y QA final

1. Seccion Data exports en report.
2. Script validate_results_contract.py.
3. Actualizar docs/output.md.

Criterio:
- si falta CORE, QA falla y se bloquea comparticion.

---

## 11) Respuesta concreta a la duda central: "que dejar de publicar"

Si el objetivo es uso real por wet lab/IP y mantener soporte bioinfo:

1. Publicar siempre:
- tablas master
- tablas significant
- summaries
- report HTML
- sashimi index + PDFs
- manifests

2. Publicar como RAW (no en la capa principal):
- outputs nativos extensos por herramienta
- binarios MAJIQ
- fromGTF/novel files de rMATS

3. No publicar (dejar en WORK_ONLY):
- temporales
- archivos de conteo crudo raramente consultados
- duplicados de bajo valor

Mensaje practico:
No hay que publicar todo lo que sale de cada tool. Hay que publicar lo que se usa para decidir y justificar hallazgos.

---

## 12) Riesgos y mitigaciones

1. Riesgo: perder artefactos utiles al reducir publicacion.
Mitigacion: usar nivel RAW para conservar evidencia tecnica importante.

2. Riesgo: confiar en work y luego perder datos.
Mitigacion: todo lo potencialmente relevante a futuro debe ir en results/raw.

3. Riesgo: romper report o joins opcionales.
Mitigacion: cambios aditivos y respeto del patron de sentinelas NO_*.

4. Riesgo: heterogeneidad entre herramientas.
Mitigacion: contrato comun de columnas + manifests por contraste.

---

## 13) Definition of Done

El plan queda correctamente implementado cuando:

1. Existe capa deliverables para todos los contrastes.
2. Cada herramienta activa genera master + significant + summary.
3. rMATS tiene tabla maestra total y tabla de significativos.
4. Sashimi se localiza por indice sin navegar manualmente subcarpetas complejas.
5. Se puede compartir solo CORE y responder preguntas biologicas principales.
6. Se puede acceder a RAW para auditoria sin depender de work.
7. QA automatica valida estructura, columnas y rutas.

---

## 14) Recomendacion final de catedra bioinfo

Una pipeline madura no es la que mas archivos genera, sino la que entrega evidencia util, interpretable y trazable con el menor ruido posible.

En este proyecto, la mejor estrategia es:

1. conservar la buena arquitectura por herramienta/contraste que ya teneis
2. anadir una capa de deliverables de alta utilidad
3. mover lo tecnico a RAW
4. dejar en WORK_ONLY solo lo verdaderamente efimero

Con esto, wet lab e IP leen lo importante, bioinfo reanaliza sin friccion y mantenimiento conserva reproducibilidad.
