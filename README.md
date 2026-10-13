# The Goats

42057 Introduction to Space Communications and Sensing, Spring 2026.
Assignment 2, Topic 2: Binary Data and Spread-Spectrum Communications
(Ingle & Proakis, *DSP Using MATLAB* 4e, Projects 12.7 and 12.8).

## Layout

| Path | Contents |
| --- | --- |
| `matlab/` | One function per block of Figs 12.16 and 12.18, the Monte-Carlo simulator, theory helpers and the experiment scripts |
| `results/` | Numeric results (`.mat`) and LaTeX table rows written by the scripts |
| `figures/` | Figures saved directly by the MATLAB scripts |
| `tools/make_report_figures.py` | Re-plots the saved results with a consistent style for the report and slides |
| `output/` | Rendered deliverables: `Assignment2_Report.pdf` and `Assignment2_Presentation.pdf` |
| `report/` | LaTeX source of the report (`main.tex`) |
| `slides/` | Beamer source of the presentation (`slides.tex`) |

## Reproducing everything

```matlab
cd matlab
run_all          % Project 12.7 then 12.8, about 4 minutes
```

The code needs no toolboxes and also runs in GNU Octave:

```sh
octave-cli --path matlab --eval run_all
python3 tools/make_report_figures.py          # needs numpy, scipy, matplotlib
cd report && pdflatex main.tex && pdflatex main.tex && cp main.pdf ../output/Assignment2_Report.pdf
cd ../slides && pdflatex slides.tex && pdflatex slides.tex && cp slides.pdf ../output/Assignment2_Presentation.pdf
```

## MATLAB functions

| Function | Block |
| --- | --- |
| `bit_generator(N)` | Binary data generator (equiprobable bits) |
| `modulator(b, M)` | Maps 1 to M samples of +1 and 0 to M samples of -1 |
| `noise_generator(L, a)` | Uniform noise on (-a, a), Pn = a^2/3 |
| `demodulator(r, M)` | Sums each block of M samples (integrate and dump) |
| `detector(y, b)` | Threshold at zero and counts bit errors |
| `pn_generator(L)` | Equiprobable +/-1 PN chips |
| `sine_generator(L, A, w0, n0)` | Jammer A sin(w0 n) with continuous phase |
| `binary_comm_sim(M, a, opts)` | Chains the blocks, with or without PN and jammer |
| `exact_ber_uniform(M, a)` | Exact BER for uniform noise (numerical convolution) |
| `ber_theory_jam(M, A, w0, Pn)` | Analytical BER with a sinusoidal jammer, with and without PN |
