# Poznámky v LaTeXu – MFF UK

LaTeX rewrites of handwritten lecture notes from courses at the
**Faculty of Mathematics and Physics, Charles University (MATFYZ / MFF UK)**.

The scanned handwritten notes are transcribed faithfully (in Czech, keeping the original
wording, formulas and diagrams), with all drawings redrawn in TikZ and the most important
parts highlighted, so the notes are easier to read, search and study from.

## Subjects

| Folder | Subject | Status |
|--------|---------|--------|
| [`QM_I/latex`](QM_I/latex) | Kvantová mechanika I (2014/2015) | done – [`QM_I.pdf`](QM_I/latex/QM_I.pdf) |
| [`Termo_I/latex`](Termo_I/latex) | Termodynamika a statistická fyzika I (NTMF043) – přednáška + cvičení | done – [`Termo_I.pdf`](Termo_I/latex/Termo_I.pdf) |

## Structure

```
<Subject>/
└── latex/
    ├── <Subject>.tex   # main document
    ├── preamble.tex    # packages, colours, highlight boxes
    ├── parts/          # one .tex file per original handwritten document
    └── <Subject>.pdf   # compiled output
```

Only the `latex/` directories are tracked; the original scans are kept locally
(see [`.gitignore`](.gitignore)).

## Conventions in the documents

- **Blue boxes** – definitions
- **Red boxes** – theorems, postulates, lemmas, key statements
- **Green boxes** – worked examples
- **Orange frames** – the most important formulas
- *Grey italics in brackets* – faint pencil notes from the original
- Margin tags (e.g. `QM-I-7`) – page label of the original handwritten notes, for cross-checking

## Building

Requires a LaTeX distribution (MiKTeX or TeX Live) with `tcolorbox`, `tikz`, `physics`
and `titlesec`. From a subject's `latex/` directory run `pdflatex` twice (for the
table of contents):

```sh
cd QM_I/latex
pdflatex QM_I.tex
pdflatex QM_I.tex
```

## Disclaimer

These are transcriptions of personal lecture notes, not official course materials.
Some faint pencil annotations could only be read partially; when in doubt, check the
original page referenced in the margin.
