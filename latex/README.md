# Electromagnetismo I — Tema 2: Desarrollo multipolar

Apuntes formalizados (LaTeX) con el mismo estilo que el Tema 0 (*Herramientas matemáticas*) y el Tema 1
(*Leyes básicas de la electrostática*): `amsart`, Palatino, enunciados numerados por sección y la misma
notación y entornos formales. La paleta pasa de azules a **magentas y rojizos claros**; las figuras que
recrean las de Mathematica del tema conservan sus colores originales.

## Compilación

```bash
latexmk -pdf main.tex      # pdfLaTeX (-shell-escape, ya en latexmkrc) + biber, ~35 s
```

**Overleaf**: subir el ZIP como proyecto nuevo (*New Project → Upload Project*), con compilador
*pdfLaTeX* y documento principal `main.tex`. La bibliografía usa *biber*, que Overleaf detecta solo.

* **Figuras externalizadas.** Todas las figuras TikZ/pgfplots se compilan una vez a PDF en
  `figuras/ext/` (librería `external`); en las compilaciones siguientes se insertan directamente si su
  código no ha cambiado (comprobación MD5). Así la compilación baja de ~4 min a ~35 s. Si se modifica una
  figura, se regenera sola en la siguiente compilación (Overleaf permite `-shell-escape`); para forzarlo,
  basta borrar su `main-figureN.pdf`.

* Las figuras «tipo Mathematica» (equipotenciales, líneas de campo y campo vectorial del dipolo y del
  cuadrupolo, puntuales y físicos) se dibujan con TikZ a partir de coordenadas precalculadas en
  `figuras/datos/*.tex`. Esos datos ya están incluidos; para regenerarlos (requiere `numpy` y `matplotlib`):

  ```bash
  python3 scripts/generar_datos.py
  ```

  El script integra las líneas de campo con separación uniforme (algoritmo de Jobard–Lefer, como
  `StreamPlot`), extrae las equipotenciales como curvas de nivel (`ContourPlot`) y muestrea el campo
  (`VectorPlot`). También genera las tablas exactas de error (`*.dat`) que dibuja `pgfplots`.
* El resto de figuras (diagramas, polinomios de Legendre, líneas de campo analíticas de la esfera
  polarizada…) se calculan directamente en TikZ/pgfplots.

## Estructura

| Fichero | Contenido |
|---|---|
| `main.tex`, `portada.tex`, `preambulo.tex`, `latexmkrc` | Documento, portada y preámbulo (estilo de los temas 0 y 1, paleta magenta) |
| `secciones/01_introduccion.tex` | Triángulo ρ–φ–E, planteamiento del problema, momentos multipolares, notación |
| `secciones/02_dipolo.tex` | Dipolo físico y puntual, potencial (Taylor), campo (esféricas y vectorial), líneas de campo, discrepancia físico/puntual |
| `secciones/03_polarizacion.tex` | Polarización, cargas de polarización ρ_P y σ_P, carga total nula |
| `secciones/04_multipolar.tex` | Desarrollo de 1/R, polinomios de Legendre, términos monopolar, dipolar y cuadrupolar, cambio de origen, simetría de revolución, ejemplos |
| `secciones/05_error.tex` | Error de truncamiento: carga fuera del origen y cota (a/r)^n |
| `secciones/06_seminarios.tex` | Problemas 2.1–2.5 resueltos |
| `secciones/07_apendice.tex` | Formulario del Tema 2, polinomios de Legendre e integrales angulares |
| `figuras/` | Todas las figuras (TikZ / pgfplots); `figuras/datos/`, los datos numéricos; `figuras/ext/`, las figuras ya compiladas |
| `scripts/generar_datos.py` | Generador documentado de los datos numéricos |
| `referencias.bib` | Bibliografía |
