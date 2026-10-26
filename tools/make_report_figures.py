"""Render the report and slide figures from the MATLAB/Octave results.

The simulations live in ../matlab (run_all.m writes results/res127.mat and
results/res128.mat). This script only re-plots those saved numbers with a
consistent style for the LaTeX report and Beamer slides.

    python3 tools/make_report_figures.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.special import erfc
from scipy.stats import beta

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "report" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    # STIX fonts embedded as TrueType keep symbols such as omega and pi as real
    # Unicode glyphs in the PDF (the default cm/Type 3 fonts map them to "!" and
    # "\xbc" when text is extracted or copied).
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "pdf.fonttype": 42,
    "font.size": 10,
    "axes.grid": True,
    "grid.alpha": 0.35,
    "legend.fontsize": 8.5,
    "legend.framealpha": 0.9,
    "savefig.bbox": "tight",
})
C = {"pn": "#1f5fbf", "nopn": "#c8372d", "th": "black",
     "M": ["#1f5fbf", "#c8372d", "#2a9d3c", "#8e44ad"]}


def load(name, key):
    m = loadmat(ROOT / "results" / name, squeeze_me=True, struct_as_record=False)
    return m[key]


def Q(x):
    return 0.5 * erfc(np.asarray(x) / np.sqrt(2))


def ber_ci(k, n, conf=0.95):
    """Clopper-Pearson interval for k errors in n bits (matches matlab/ber_ci.m)."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    a = 1 - conf
    lo = np.where(k > 0, beta.ppf(a / 2, np.maximum(k, 1), n - k + 1), 0.0)
    hi = np.where(k < n, beta.ppf(1 - a / 2, k + 1, np.maximum(n - k, 1)), 1.0)
    return lo, hi


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf")
    fig.savefig(OUT / f"{name}.png", dpi=200)
    plt.close(fig)


# ---------------------------------------------------------------- 12.7
r = load("res127.mat", "res127")

# waveforms through each block
M = int(r.wave_M)
s, rx, y = r.wave_s, r.wave_r, r.wave_y
b, bh = np.atleast_1d(r.wave_b).astype(int), np.atleast_1d(r.wave_bhat).astype(int)
n = np.arange(len(s))
fig, ax = plt.subplots(3, 1, figsize=(6.4, 5.2), sharex=False)
ax[0].step(n, s, where="post", color=C["pn"])
ax[0].set_ylim(-1.6, 1.6); ax[0].set_ylabel("$s(n)$")
ax[0].set_title(f"Modulator output, bits = {' '.join(map(str, b))}, $M={M}$")
ax[1].plot(n, rx, ".-", color="0.35", ms=3, lw=0.7, label="$r(n)=s(n)+w(n)$")
ax[1].step(n, s, where="post", color=C["pn"], lw=1, alpha=0.6, label="$s(n)$")
ax[1].set_ylabel("$r(n)$"); ax[1].legend(loc="upper right", ncol=2)
ax[1].set_title(r"Channel output, $w(n)\sim\mathcal{U}(-\sqrt{3},\sqrt{3})$, $P_n=1$")
k = np.arange(len(y))
ax[2].stem(k, y, basefmt="k-")
ax[2].axhline(0, color="k", lw=0.8, ls="--")
ax[2].set_ylabel("$y(k)$"); ax[2].set_xlabel("bit index $k$  (top two panels: sample index $n$)")
ax[2].set_title(f"Demodulator output, detected bits = {' '.join(map(str, bh))}")
fig.tight_layout()
save(fig, "p127_waveforms")

# BER vs SNR
fig, ax = plt.subplots(figsize=(6.4, 4.2))
snr = np.linspace(-1, 14, 300)
ax.semilogy(snr, Q(np.sqrt(10 ** (snr / 10))), color=C["th"], lw=1.2,
            label=r"Gaussian approx. $Q(\sqrt{\mathrm{SNR}})$")
ax.semilogy(r.snr_a_db, r.exact_a, "--", color=C["pn"], lw=1, label="Exact (Irwin–Hall), sweep (a)")
ax.semilogy(r.snr_b_db, r.exact_b, ":", color=C["nopn"], lw=1.5, label="Exact (Irwin–Hall), sweep (b)")
for snr_db, ber, k, n, mk, col, lab in [
        (r.snr_a_db, r.ber_a, r.nerr_a, r.nbits_a, "o", C["pn"],
         r"Simulated (a): $P_n=1$, $M=1\ldots20$"),
        (r.snr_b_db, r.ber_b, r.nerr_b, r.nbits_b, "s", C["nopn"],
         r"Simulated (b): $M=10$, $P_n$ varied")]:
    lo, hi = ber_ci(k, n)
    ok = ber > 0
    ax.errorbar(snr_db[ok], ber[ok], yerr=[ber[ok] - lo[ok], hi[ok] - ber[ok]], fmt=mk,
                mfc="none", color=col, ms=5, capsize=2, elinewidth=0.8, label=lab)
ax.set_xlim(-1, 14); ax.set_ylim(1e-8, 1)
ax.set_xlabel(r"SNR $= M/P_n$ (dB)"); ax.set_ylabel("Bit error rate")
ax.legend(loc="lower left", title="Error bars: 95% Clopper–Pearson", title_fontsize=8)
save(fig, "p127_ber")

# ---------------------------------------------------------------- 12.8
r = load("res128.mat", "res128")
Ms = np.atleast_1d(r.M_list).astype(int)

# A = M, the brief's experiment
fig, ax = plt.subplots(figsize=(6.4, 3.8))
x = np.arange(len(Ms)); wbar = 0.2
labels = [(0, f"No PN, $\\omega_0={r.w_in:g}$", C["nopn"], None),
          (1, f"PN, $\\omega_0={r.w_in:g}$", C["pn"], None),
          (2, f"No PN, $\\omega_0={r.w_out:g}$", C["nopn"], "//"),
          (3, f"PN, $\\omega_0={r.w_out:g}$", C["pn"], "//")]
for j, (c, lab, col, hatch) in enumerate(labels):
    xs = x + (j - 1.5) * wbar
    ax.bar(xs, r.T[:, c], wbar, color=col if hatch is None else "white",
           edgecolor=col, hatch=hatch, label=lab)
    ax.plot(xs, r.Tth[:, c], "k_", ms=12, mew=1.5)
ax.plot([], [], "k_", ms=12, mew=1.5, label="Theory")
ax.set_xticks(x, [f"$M=A={m}$" for m in Ms])
ax.set_ylabel("Bit error rate"); ax.set_ylim(0, 0.6)
ax.legend(ncol=3, loc="upper center", fontsize=8)
save(fig, "p128_ber_A_eq_M")

# BER vs jammer-to-signal ratio, in-band jammer
fig, ax = plt.subplots(figsize=(6.4, 4.3))
jsr = r.jsr_db
for i, m in enumerate(Ms):
    col = C["M"][i]
    ax.semilogy(jsr, r.Th1[i], "-", color=col, lw=1.2)
    ax.semilogy(jsr, np.where(r.B1[i] > 0, r.B1[i], np.nan), "o", mfc="none", color=col, ms=4)
    ax.semilogy(jsr, r.Th0[i], "--", color=col, lw=0.9)
    ax.semilogy(jsr, np.where(r.B0[i] > 0, r.B0[i], np.nan), "x", color=col, ms=4)
    ax.axvline(10 * np.log10(m ** 2 / 2), color=col, ls=":", lw=0.9)
    ax.plot([], [], "s", color=col, label=f"$M={m}$")
ax.plot([], [], "k-o", mfc="none", ms=4, label="PN (theory / sim)")
ax.plot([], [], "k--x", ms=4, label="No PN (theory / sim)")
ax.plot([], [], "k:", label="$A=M$")
ax.set_xlim(jsr.min(), jsr.max()); ax.set_ylim(1e-5, 1)
ax.set_xlabel(r"Jammer-to-signal power ratio $A^2/2$ (dB)")
ax.set_ylabel("Bit error rate")
ax.set_title(f"In-band jammer, $\\omega_0={r.w_in:g}$ rad/sample, $P_n=1$", fontsize=10)
ax.legend(loc="lower right", ncol=2)
save(fig, "p128_ber_vs_jsr")

# BER vs jammer frequency
fig, axs = plt.subplots(2, 1, figsize=(6.4, 5.4))
Mw = 100
for q, ax in enumerate(axs):
    A = np.atleast_1d(r.A_pair)[q]
    ax.plot(r.w_fine / np.pi, r.Wt0[q], "-", color=C["nopn"], lw=0.6, label="No PN (theory)")
    ax.plot(r.w_list / np.pi, r.W0[q], "x", color=C["nopn"], ms=4, label="No PN (sim)")
    ax.plot(r.w_fine / np.pi, r.Wt1[q], "-", color=C["pn"], lw=1.4, label="PN (theory)")
    ax.plot(r.w_list / np.pi, r.W1[q], "o", mfc="none", color=C["pn"], ms=3.5, label="PN (sim)")
    ax.set_title(f"$M={Mw}$, $A={A:g}$, $P_n=1$", fontsize=10)
    ax.set_ylabel("Bit error rate"); ax.set_xlim(0, 1)
    if q == 0:
        ax.set_ylim(0, 0.72)
        ax.legend(loc="upper right", ncol=4, fontsize=7.5)
    else:
        ax.set_yscale("log"); ax.set_ylim(1e-5, 1)
axs[-1].set_xlabel(r"Jammer frequency $\omega_0/\pi$")
fig.tight_layout()
save(fig, "p128_ber_vs_w0")

# spectra before / after despreading
om = r.psd_om / np.pi
dB = lambda P: 10 * np.log10(np.maximum(P, 1e-6))
fig, axs = plt.subplots(2, 1, figsize=(6.4, 5.0), sharex=True)
axs[0].plot(om, dB(r.Psp), color=C["pn"], lw=0.8, label=r"$s(n)p(n)$ (spread signal)")
axs[0].plot(om, dB(r.Pj), color=C["nopn"], lw=1, label=r"$i(n)$ (jammer)")
axs[0].plot(om, dB(r.Pw), color="0.5", lw=0.6, label=r"$w(n)$ (noise)")
axs[0].set_title("In the channel (before despreading)", fontsize=10)
axs[1].plot(om, dB(r.Ps), color=C["pn"], lw=0.8, label=r"$s(n)$ (despread signal)")
axs[1].plot(om, dB(r.Pjd), color=C["nopn"], lw=0.8, label=r"$i(n)p(n)$ (spread jammer)")
axs[1].plot(om, dB(r.Pwd), color="0.5", lw=0.6, label=r"$w(n)p(n)$")
axs[1].set_title("After despreading (demodulator input)", fontsize=10)
for ax in axs:
    ax.set_ylim(-40, 50); ax.set_ylabel("PSD (dB)"); ax.legend(loc="upper right")
axs[1].set_xlim(0, 1); axs[1].set_xlabel(r"$\omega/\pi$")
fig.suptitle(f"$M={int(r.psd_M)}$, $A={r.psd_A:g}$, $\\omega_0={r.psd_w0/np.pi:g}\\pi$, $P_n=1$",
             fontsize=10)
fig.tight_layout()
save(fig, "p128_spectra")
print("figures written to", OUT)

# jammer-to-signal ratio tolerated at BER = 1e-3 (theory curves, in-band jammer)
rows = []
for i, m in enumerate(Ms):
    vals = []
    for Th in (r.Th0, r.Th1):
        y = np.log10(np.maximum(Th[i], 1e-30))
        j = int(np.argmax(y > -3))
        x = r.jsr_db
        vals.append(x[j - 1] + (x[j] - x[j - 1]) * (-3 - y[j - 1]) / (y[j] - y[j - 1]))
    pg = 10 * np.log10(m)
    rows.append(f"{m} & {pg:.1f} & {vals[0]:.1f} & {vals[1]:.1f} & {vals[1]-vals[0]:.1f} & "
                f"{10*np.log10(m*m/2):.1f} \\\\")
(ROOT / "results" / "table_margin.tex").write_text("\n".join(rows) + "\n")
print("margin table written")
