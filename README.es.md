# Vina-Ion-Metal

[English](README.md) | **Español**

[![Licencia: Apache 2.0](https://img.shields.io/badge/Licencia-Apache_2.0-blue.svg)](LICENSE)
[![CI](https://github.com/cristhianyaranganana-arch/Vina-Ion-Metal/actions/workflows/build.yml/badge.svg)](https://github.com/cristhianyaranganana-arch/Vina-Ion-Metal/actions/workflows/build.yml)
[![Release](https://img.shields.io/github/v/release/cristhianyaranganana-arch/Vina-Ion-Metal)](https://github.com/cristhianyaranganana-arch/Vina-Ion-Metal/releases)
![Plataforma: Linux x64](https://img.shields.io/badge/plataforma-Linux%20x64-blue.svg)

**Fork de AutoDock Vina con scoring nativo de coordinación de metales de
transición basado en física, detección automática de ocupación del sitio de
unión y un pipeline robusto de preparación de receptores.**

Vina-Ion-Metal 1.0.0 es el sucesor corregido y completado de
[VinaMetal 1.0.0](https://github.com/cristhianyaranganana-arch/VinaMetal).
Está construido sobre AutoDock Vina 1.2.7 e incorpora las correcciones de los
16 problemas identificados por la auditoría técnica de 2026 sobre VinaMetal,
además de defectos adicionales encontrados durante esta revisión (incluido un
trigger de CI corrupto que desactivaba silenciosamente la integración
continua).

> **Para qué sirve:** docking en sitios activos de metaloenzimas donde la
> coordinación ligando-metal es parte del mecanismo de unión (p. ej.
> GlfT2/Mn²⁺ en *M. tuberculosis*, metalo-β-lactamasas, anhidrasas
> carbónicas, HDACs). Vina-Ion-Metal es un **filtro mecanístico que guía las
> poses hacia la coordinación correcta del metal** — no es un predictor de
> afinidad, y sus pesos de scoring aún no están calibrados por benchmark
> (ver [Limitaciones](#limitaciones)).

## El problema

Las funciones de scoring empíricas genéricas tratan los iones metálicos como
átomos pesados ordinarios. No pueden premiar el número de coordinación
correcto, penalizar la sobre-coordinación ni distinguir sitios tetraédricos
de octaédricos. El Vina estándar coloca el tipo genérico de metal a una
distancia óptima de contacto de ~2,9 Å, lejos de las distancias reales de
coordinación (~2,0–2,4 Å), por lo que las campañas de docking contra
metaloenzimas fallan sistemáticamente la geometría catalítica.

## La solución

Vina-Ion-Metal integra la física metálica de AutoDock 4 directamente en la
función de scoring nativa de Vina, con gradientes analíticos consumidos por el
motor ILM+BFGS sin modificar:

- **Van der Waals 12-6 específico del metal** con distancias de coordinación
  de Harding (por ion y elemento donante N/O/S), tomado de las tablas de
  parámetros de AD4.
- **Penalización por número de coordinación (CN)** con dos regímenes:
  - *lenient* (por defecto): solo se penaliza la sobre-coordinación (bisagra
    suave C-infinito) — una esfera incompleta del ligando no se castiga;
  - *strict* (`--metal_cn_strict`): cualquier desviación del objetivo de CN se
    penaliza con la cuadrática pura.
- **Bonus acotado por geometría de coordinación** (tetraédrica, octaédrica,
  plana-cuadrada).
- **Detección automática de ocupación del receptor**: los donantes del
  receptor rígido ya enlazados a cada metal se cuentan, y el objetivo de CN
  del ligando se convierte en
  `max(0, CN_ideal − donantes_enlazados_del_receptor)` — sin configuración por
  sitio. Un conteo geométrico independiente del tipado también se reporta
  como diagnóstico de calidad de preparación cuando ambos conteos difieren.
- **Preparación robusta de receptores** (`scripts/prepare_receptor.py`):
  resolución de altLoc, verificaciones de integridad de cadena, tipado de
  metales **incluido Cu²⁺/Ni²⁺ (que Meeko no puede emitir)** y un reporte
  pre-flight de coordinación.

Los pares metal-donante con iones parametrizados se excluyen de los
potenciales estándar de Vina para evitar doble conteo. Los receptores con
metales **no parametrizados** (Na⁺, K⁺, Co, Cd, Hg, …) conservan el
tratamiento exacto de Vina 1.2.7, preservando la promesa de identidad
byte a byte.

## Validación

| Sistema | Ion | Donantes detectados | Objetivo CN ligando | Resultado (corridas de validación registradas, semilla 42, exhaustividad 16) |
|---|---|---|---|---|
| GlfT2, PDB 4FIY cadena A | Mn²⁺ | 4 | 2 | **redocking cristal v1.0.2**: pose top-1 −8.86 kcal/mol (metal-consciente), RMSD simétrico 1,50 Å contra la pose UDP depositada, coordina el Mn²⁺ con tres oxígenos a 2,12/2,18/2,51 Å (Harding Mn–O = 2,15 Å); los 9 modos coordinan el ion |
| PDB 1s63 | Zn²⁺ | 3 (geométrico 4, con diagnóstico) | 1 | pose top-1 coordina el Zn²⁺ a 2,26 Å (poses sin cambios desde la corrida de la era v1.0.0) |
| PDB 1iep (sin metal) | — | — | — | salida `--local_only` **idéntica byte a byte** a AutoDock Vina 1.2.7 (verificado en CI) |

El ajuste automático se imprime antes del docking, por ejemplo:

    Vina-Ion-Metal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

El RMSD es de átomos pesados con emparejamiento por simetría (RDKit
`GetBestRMS`) contra el ligando cristal. Dado que el régimen de CN por defecto
y los pesos por defecto cambiaron respecto a VinaMetal 1.0.0 (auditorías #2 y
#9), los números de re-docking pueden variar; ejecute `examples/*/config.txt`
para reproducirlos. Desde v1.0.2 el ejemplo GlfT2 incluye el ligando cristal
depositado, por lo que su redocking se mide por RMSD contra la pose cristal
(vea `docs/validation/`), y las afinidades reportadas incluyen el término
metálico (corrección F6).

## Inicio rápido

Requisitos: CMake ≥ 3.12, un compilador C++11 y Boost ≥ 1.65 con las
librerías `system`, `thread`, `serialization`, `filesystem` y
`program_options`. Se prefiere Boost estático por defecto (binarios
portables); use `-DBoost_USE_STATIC_LIBS=OFF` si su sistema solo tiene
librerías dinámicas.

    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build -j
    ctest --test-dir build --output-on-failure

Prepare un receptor (metales tipados automáticamente, con soporte Cu/Ni):

    scripts/prepare_receptor.py --pdb receptor.pdb --out receptor.pdbqt

Ejecute el ejemplo de GlfT2 (el Mn²⁺ se detecta automáticamente):

    cd examples/glft2
    ../../build/vina --config config.txt --metals auto --out glft2_out.pdbqt

El log imprime el objetivo ajustado por ocupación:

    Vina-Ion-Metal: parametrized metal ions in receptor: Mn(x1) - coordination scoring is ACTIVE
    Vina-Ion-Metal: Mn2+ CN target = ideal 6 - receptor-bound 4 = 2

Use `--metals off` para desactivar el scoring de coordinación y recuperar el
comportamiento de Vina puro. Vea [docs/USAGE.md](docs/USAGE.md) para todas las
opciones (incluido el protocolo recomendado multi-semilla),
[docs/INSTALL.md](docs/INSTALL.md) para compilaciones por plataforma y
[docs/LIMITATIONS.md](docs/LIMITATIONS.md) antes de un uso en producción.

## Qué se corrigió respecto a VinaMetal 1.0.0

La trazabilidad completa está en [CHANGELOG.md](CHANGELOG.md). Resumen de los
16 problemas auditados:

| # | Problema | Estado en Vina-Ion-Metal 1.0.0 |
|---|---|---|
| 1 | Sobre-exclusión de grids para cualquier átomo de elemento `Met` | **Corregido** — solo se excluyen iones parametrizados; test de regresión añadido |
| 2 | Opción muerta `--metal_cn_strict` | **Implementada** — penalización CN lenient (por defecto) vs strict, suave C-infinito |
| 3 | Variable muerta `dE_dcos` | **Eliminada** |
| 4 | Branding apuntaba a Vina upstream | **Corregido** — cadena de versión, issue tracker y citas son de Vina-Ion-Metal |
| 5 | Cu²⁺/Ni²⁺ inalcanzables vía Meeko | **Corregido** — `prepare_receptor.py` los tipa nativamente |
| 6 | Fallos silenciosos de `--metal_oxidation` | **Corregido** — validación estricta (elemento, entero, rango) |
| 7 | Fallback silencioso de estado de oxidación | **Corregido** — advertencia con los estados parametrizados al cargar el receptor |
| 8 | `data/AD4Zn.dat` huérfano | **Eliminado** (era una copia mal etiquetada de `AD4.1_bound.dat`) |
| 9 | Pesos sin calibrar (1.0/1.0/1.0) | **Por defecto heurísticos** (1.2/0.5/1.0) con aviso en ejecución; benchmark de calibración pendiente |
| 10 | Detección de ocupación frágil | **Endurecida** — conteo geométrico de donantes + advertencias por divergencia |
| 11 | Preparación de receptor frágil | **Reemplazada** — pipeline propio robusto (altLoc, verificación de integridad, reporte pre-flight) |
| 12 | Binarios no portables por Boost dinámico | **Corregido** — Boost estático por defecto (con fallback elegante) |
| 13 | Test de regresión no ejecutable out-of-the-box | **Corregido** — registrado en CTest; el CI descarga la referencia Vina 1.2.7 |
| 14 | CI solo Linux | **Mantenido** (alcance intencional); bug del trigger de CI encontrado y corregido |
| 15 | RMSD por orden de átomos (benchmark del auditor) | **Documentado** — el protocolo recomendado usa RMSD con simetría |
| 16 | Una semilla / exhaustividad baja | **Documentado** — protocolo recomendado (≥3 semillas, exhaustividad ≥16) |

Correcciones adicionales fuera de la auditoría: YAML corrupto
(`branches: ain]`) en el workflow de CI (el CI nunca se disparaba) y la cita
de AutoDock4Zn en el README con páginas incorrectas (ver abajo).

## Limitaciones

- Los pesos del scoring metálico son **heurísticas provisionales**, no
  calibradas contra un benchmark multiproteína; trate Vina-Ion-Metal como un
  filtro mecanístico de poses, no como un predictor de afinidad.
- Fe²⁺/Fe³⁺ comparten los parámetros genéricos de hierro de AD4; el refinamiento
  por estado de oxidación está pendiente.
- Los parámetros de Cu²⁺/Ni²⁺ siguen siendo estimaciones estilo AD4 (sin
  calibrar en AD4); sobrescríbalos con `--metal_params` para trabajo serio.
- El receptor es rígido; las cadenas laterales coordinantes no se reorganizan.
- Los maps pregenerados deben provenir de Vina-Ion-Metal (los iones metálicos
  se excluyen de los maps de donantes).

Vea la lista completa en [docs/LIMITATIONS.md](docs/LIMITATIONS.md).

## Protocolo de docking recomendado

1. Prepare el receptor con `scripts/prepare_receptor.py`; lea el reporte
   pre-flight de coordinación y resuelva cualquier advertencia.
2. Haga docking con `--exhaustiveness 16` (o más) y **al menos 3 semillas
   distintas**; los rankings por semilla con baja exhaustividad están
   dominados por ruido.
3. Evalúe las poses por RMSD de átomos pesados con simetría (p. ej. RDKit
   `GetBestRMS`), nunca por orden de átomos.
4. Trate los hits coordinantes como **generadores de hipótesis**; re-ranquee
   con MM/GBSA o MD corta, y aplique contraste anti-PAINS/de agregación antes
   de cualquier ensayo bioquímico.

## Citación

    @software{vina_ion_metal2026,
      title        = {Vina-Ion-Metal: AutoDock Vina fork with native
                      transition-metal coordination scoring},
      author       = {Yaranga Nana, Christian},
      year         = {2026},
      version      = {1.0.0},
      url          = {https://github.com/cristhianyaranganana-arch/Vina-Ion-Metal},
      note         = {Corrected successor of VinaMetal 1.0.0; based on
                      AutoDock Vina 1.2.7}
    }

Este trabajo se basa en VinaMetal 1.0.0 (Yaranga Nana, 2026,
DOI 10.5281/zenodo.23499212) — cítelo al comparar resultados entre versiones.
Cite también los fundamentos:

- Trott, O.; Olson, A. J. *J. Comput. Chem.* **2010**, 31, 455–461.
- Eberhardt, J.; Santos-Martins, D.; Tillack, A. F.; Forli, S.
  *J. Chem. Inf. Model.* **2021**, 61, 3891–3898.
- Santos-Martins, D.; Forli, S.; Ramos, M. J.; Olson, A. J. AutoDock4Zn: an
  improved AutoDock force field for small-molecule docking to zinc
  metalloproteins. *J. Chem. Inf. Model.* **2014**, 54(8), 2371–2379.
- Harding, M. M. *Acta Crystallogr. D* **2001**, 57, 401–411.

## Créditos y licencia

Vina-Ion-Metal es desarrollado por **Christian Yaranga Nana** (Universidad
Norbert Wiener, Lima, Perú) como sucesor corregido de VinaMetal 1.0.0. Está
construido sobre AutoDock Vina y AutoDock 4, desarrollados en The Scripps
Research Institute; sus avisos de copyright originales se conservan en los
archivos fuente (ver [NOTICE](NOTICE)). Las entradas de ejemplo derivan de las
entradas PDB **4FIY**, **1s63** y **1iep** (RCSB Protein Data Bank).

El código fuente se distribuye bajo la Licencia Apache, Versión 2.0; ver
[LICENSE](LICENSE). Los archivos de parámetros de AutoDock 4 en
`data/ad4_reference/` se distribuyen bajo la GNU General Public License
v2.0; ver `data/ad4_reference/COPYING`.
