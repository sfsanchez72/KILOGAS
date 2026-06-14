#!/usr/bin/env python3
"""
Create the WHaD / BPT diagnostic figure from a Pipe3D KILOGAS cube.

This script is a command-line version of the notebook
`WHaD_BPT_Pipe3D_KG.ipynb`.

Required inputs
---------------
--name       Object/cube name, e.g. KG-SAMI-63777
--input-dir  Directory containing <name>.Pipe3D.cube.fits.gz
--output-dir Directory where the diagnostic PDF will be written

Expected files
--------------
<input-dir>/<name>.Pipe3D.cube.fits.gz

Optional background image, searched in this order:
<input-dir>/PanStarss/<name>_g.fits
<input-dir>/<name>_g.fits

If the g-band FITS file is not found, the script uses the Pipe3D pseudo-V map
as the background for the top-left panel.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import rcParams as rc

import scipy.ndimage as ndimage
warnings.simplefilter("ignore")

ALPHA = 0.7


def configure_matplotlib(usetex: bool = False) -> None:
    """Set plotting style close to the original notebook."""
    rc.update({
        "font.size": 20,
        "font.weight": 900,
        "text.usetex": usetex,
        "path.simplify": True,
        "xtick.labelsize": 20,
        "ytick.labelsize": 20,
        "axes.linewidth": 2.0,
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.major.width": 1,
        "ytick.major.width": 1,
        "lines.markeredgewidth": 1,
        "legend.numpoints": 1,
        "legend.frameon": False,
        "legend.handletextpad": 0.3,
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "axes.facecolor": "w",
    })


def find_file(candidates: list[Path], description: str, required: bool = True) -> Path | None:
    """Return the first existing path among candidates."""
    for path in candidates:
        if path.exists():
            return path
    if required:
        msg = f"Could not find {description}. Tried:\n" + "\n".join(f"  - {p}" for p in candidates)
        raise FileNotFoundError(msg)
    return None


def safe_log10(arr: np.ndarray) -> np.ndarray:
    """log10 with invalid values turned into NaN."""
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.log10(arr)
    return out


def my_contour(
    ax,
    x_cont: np.ndarray,
    y_cont: np.ndarray,
    x_min: float,
    x_max: float,
    y_min: float,
    y_max: float,
    c_color: str = "red",
    title: str = "",
    nbins: int = 30,
    zorder: int = 1,
    linewidths: float = 2,
    alpha: float = 0.75,
) -> None:
    """Draw density contours, following the notebook implementation."""
    mask = (
        np.isfinite(x_cont) & np.isfinite(y_cont) &
        (x_cont > x_min) & (x_cont < x_max) &
        (y_cont > y_min) & (y_cont < y_max)
    )
    if np.count_nonzero(mask) < 10:
        return

    x_plt = x_cont[mask]
    y_plt = y_cont[mask]
    counts, xbins, ybins = np.histogram2d(
        x_plt,
        y_plt,
        bins=nbins,
        density=True,
        range=[[x_min, x_max], [y_min, y_max]],
    )
    if not np.isfinite(counts).any() or np.nanmax(counts) <= 0:
        return

    counts = ndimage.gaussian_filter(counts, sigma=1, order=0)
    if np.nanmax(counts) <= 0:
        return
    counts = counts / np.nanmax(counts)
    sum_total = np.nansum(counts)
    if sum_total <= 0:
        return

    vals = []
    levels = []
    for cut in np.arange(0.00, 1.0, 0.01):
        mask_now = counts > cut
        levels.append(cut)
        vals.append(np.nansum(counts[mask_now]) / sum_total)

    vals = np.asarray(vals)
    levels = np.asarray(levels)
    vals_cont = np.array([0.95, 0.65, 0.40])
    levels_cont = np.interp(vals_cont, vals[::-1], levels[::-1])
    levels_cont = np.unique(np.sort(levels_cont))
    if len(levels_cont) == 0:
        return

    counts_rot = np.rot90(counts, 3)
    xbins = xbins + 0.5 * (x_max - x_min) / nbins
    ybins = ybins + 0.5 * (y_max - y_min) / nbins
    flip_counts_rot = np.fliplr(counts_rot)
    p_cont = ax.contour(
        xbins[0:nbins],
        ybins[0:nbins],
        flip_counts_rot,
        levels_cont,
        colors=c_color,
        alpha=alpha,
        linewidths=linewidths,
        zorder=zorder,
    )
    if title and len(p_cont.collections) > 0:
        p_cont.collections[0].set_label(title)


def get_background_image(input_dir: Path, name: str, fallback: np.ndarray) -> tuple[np.ndarray, str]:
    """Read optional Pan-STARRS g image; otherwise use a fallback Pipe3D map."""
    from astropy.io import fits
    bg_path = find_file(
        [
            input_dir / "PanStarss" / f"{name}_g.fits",
            input_dir / "PanSTARRS" / f"{name}_g.fits",
            input_dir / f"{name}_g.fits",
        ],
        "optional g-band FITS image",
        required=False,
    )
    if bg_path is None:
        return np.asarray(fallback, dtype=float), "Pipe3D pseudo-V"
    with fits.open(bg_path) as hdu:
        return np.asarray(hdu[0].data, dtype=float), "g-band"


def read_pipe3d_cube(input_dir: Path, name: str):
    from astropy.io import fits
    pipe3d_path = find_file(
        [input_dir / f"{name}.Pipe3D.cube.fits.gz", input_dir / f"{name}.Pipe3D.cube.fits"],
        "Pipe3D cube",
        required=True,
    )
    return fits.open(pipe3d_path), pipe3d_path


def compute_diagnostics(elines_data: np.ndarray, name: str) -> dict[str, np.ndarray | float]:
    """Compute BPT and WHaD diagnostic maps/masks from Pipe3D emission-line planes."""
    nf = 0.25

    img_Ha = elines_data[45, :, :]
    img_e_Ha = nf * elines_data[273, :, :]
    img_WHa = elines_data[216, :, :]
    img_e_WHa = nf * elines_data[444, :, :]
    img_disp = elines_data[159, :, :]
    img_e_disp = elines_data[387, :, :]
    img_OIII = elines_data[26, :, :]
    img_e_OIII = nf * elines_data[254, :, :]
    img_Hb = elines_data[28, :, :]
    img_e_Hb = nf * elines_data[256, :, :]
    img_NII = elines_data[46, :, :]
    img_e_NII = nf * elines_data[274, :, :]

    inst_disp = 1.0 if "SAMI" in name else 1.4
    with np.errstate(invalid="ignore", divide="ignore"):
        img_disp_kms = np.sqrt((img_disp / 2.354) ** 2 - inst_disp**2) / 6562.0 * 300000.0
        img_disp_kms = img_disp_kms / np.sqrt(2.354)

    O3 = safe_log10(img_OIII / img_Hb)
    N2 = safe_log10(img_NII / img_Ha)

    mask_BPT = (
        (img_Ha > 3 * img_e_Ha) &
        (img_OIII > img_e_OIII) &
        (img_NII > img_e_NII) &
        (img_Hb > img_e_Hb) &
        (np.abs(img_WHa) > img_e_WHa)
    )
    mask_ND = (
        (img_Ha > 3 * img_e_Ha) &
        (np.abs(img_WHa) > img_e_WHa) &
        (img_disp > img_e_disp)
    )

    ny, nx = img_Ha.shape
    img_Ha_NII = img_Ha + img_NII
    img_diag_BPT = np.zeros((ny, nx))
    img_diag_new = np.zeros((ny, nx))

    cut_Kew = 0.61 / (N2 - 0.47) + 1.19
    mask_SF_BPT = ((O3 < cut_Kew) & (N2 < 0.2)) & (np.abs(img_WHa) > 6)
    mask_sAGN_BPT = (O3 >= cut_Kew) & (np.abs(img_WHa) > 10)
    mask_wAGN_BPT = (O3 >= cut_Kew) & (np.abs(img_WHa) > 3) & (np.abs(img_WHa) <= 10)
    mask_RG_BPT = np.abs(img_WHa) < 3

    img_diag_BPT[mask_RG_BPT] = 1
    img_diag_BPT[mask_SF_BPT] = 2
    img_diag_BPT[mask_sAGN_BPT] = 3
    img_diag_BPT[mask_wAGN_BPT] = 4
    img_diag_BPT[~mask_BPT] = 0

    cut_vel = 1.75
    log_disp = safe_log10(img_disp_kms)
    mask_SF_new = (img_disp > 0) & (np.abs(img_WHa) > 6) & (log_disp < cut_vel)
    mask_sAGN_new = (img_disp > 0) & (np.abs(img_WHa) > 10) & (log_disp > cut_vel)
    mask_wAGN_new = (img_disp > 0) & (np.abs(img_WHa) > 3) & (np.abs(img_WHa) < 10) & (log_disp > cut_vel)
    mask_UK_new = (img_disp < 0) | ((np.abs(img_WHa) > 3) & (np.abs(img_WHa) < 6) & (log_disp < cut_vel))
    mask_RG_new = np.abs(img_WHa) < 3

    img_diag_new[mask_RG_new] = 1
    img_diag_new[mask_SF_new] = 2
    img_diag_new[mask_sAGN_new] = 3
    img_diag_new[mask_wAGN_new] = 4
    img_diag_new[mask_UK_new] = 0
    img_diag_new[~mask_ND] = 0

    return {
        "img_Ha": img_Ha,
        "img_WHa": img_WHa,
        "img_disp": img_disp,
        "img_disp_kms": img_disp_kms,
        "img_Ha_NII": img_Ha_NII,
        "O3": O3,
        "N2": N2,
        "mask_BPT": mask_BPT,
        "mask_ND": mask_ND,
        "mask_SF_BPT": mask_SF_BPT,
        "mask_sAGN_BPT": mask_sAGN_BPT,
        "mask_wAGN_BPT": mask_wAGN_BPT,
        "mask_RG_BPT": mask_RG_BPT,
        "mask_SF_new": mask_SF_new,
        "mask_sAGN_new": mask_sAGN_new,
        "mask_wAGN_new": mask_wAGN_new,
        "mask_UK_new": mask_UK_new,
        "mask_RG_new": mask_RG_new,
        "img_diag_BPT": img_diag_BPT,
        "img_diag_new": img_diag_new,
    }


def finite_contour_levels(image: np.ndarray) -> tuple[float, np.ndarray]:
    """Return Halpha+[NII] threshold and contour levels used in the notebook."""
    arr = np.asarray(image, dtype=float)
    valid = np.isfinite(arr) & (arr > 0)
    if np.count_nonzero(valid) == 0:
        return np.nan, np.array([])
    med = np.nanmedian(arr[valid])
    std = np.nanstd(arr[valid])
    ha_min = med + 0.25 * std
    if not np.isfinite(ha_min) or ha_min <= 0 or not np.isfinite(std) or std <= 0:
        return ha_min, np.array([])
    levels = np.log10(ha_min) + abs(np.log10(std)) * np.arange(0, 10, 2)
    return ha_min, levels[np.isfinite(levels)]


def overlay_ha_contours(ax, img_Ha_NII: np.ndarray, levels: np.ndarray, extent, color: str = "black", alpha: float = ALPHA) -> None:
    if levels.size == 0:
        return
    with np.errstate(divide="ignore", invalid="ignore"):
        z = np.log10(img_Ha_NII)
    try:
        ax.contour(z, levels, linestyles="solid", colors=color, linewidths=1.0,
                   alpha=alpha, extent=extent, origin="lower")
    except Exception:
        pass


def plot_diag_map(ax, fig, image, select_reg, img_Ha_NII, ha_min, levels, cmap_t, extent, title):
    image = np.ma.masked_invalid(image)
    image = image * select_reg
    image = np.ma.masked_array(image, img_Ha_NII < ha_min)
    im = ax.imshow(image, interpolation="nearest", cmap=cmap_t, vmin=0, vmax=4,
                   alpha=ALPHA, extent=extent, origin="lower")
    cb = fig.colorbar(im, ax=ax, orientation="vertical", fraction=0.044, pad=0.02)
    cb.set_ticks([0.5, 1.25, 2, 2.75, 3.5])
    cb.ax.set_yticklabels(["Unk.", "Ret.", "SF", "wAGN", "sAGN"])
    overlay_ha_contours(ax, img_Ha_NII, levels, extent, color="black", alpha=ALPHA)
    ax.text(0.04, 0.08, title, fontsize=19, bbox={"facecolor": "white", "pad": 5},
            va="bottom", ha="left", transform=ax.transAxes)
    ax.set_xticks([-15, -7.5, 0, 7.5, 15])
    ax.set_yticks([-15, -7.5, 0, 7.5, 15])
    ax.set_xlabel(r"$\Delta$ RA (arcsec)")
    ax.set_ylabel(r"$\Delta$ DEC (arcsec)")


def plot_bpt_panel(ax, diag: dict, cmap_t, mode: str) -> None:
    O3 = diag["O3"]
    N2 = diag["N2"]
    c_DL = "black"
    x_min, x_max = -1.45, 0.7
    y_min, y_max = -1.65, 1.3
    x = np.linspace(-1.45, 0.3, 100)
    cut_Kau = 0.61 / (x - 0.05) + 1.3
    cut_Kew = 0.61 / (x - 0.47) + 1.19
    ax.plot(x[x < 0], cut_Kau[x < 0], "-", color=c_DL, linewidth=2)
    ax.plot(x, cut_Kew, "-.", color=c_DL, linewidth=2)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    hex_c = [mcolors.rgb2hex(cmap_t(i)[:3]) for i in range(cmap_t.N)]
    spx = 20 if mode == "BPT" else 15

    my_contour(ax, N2, O3, x_min, x_max, y_min, y_max, c_color="black",
               nbins=50, linewidths=1.5, alpha=0.7)

    if mode == "BPT":
        ax.scatter(N2[~diag["mask_BPT"]], O3[~diag["mask_BPT"]], color=hex_c[0], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        ax.scatter(N2[diag["mask_RG_BPT"]], O3[diag["mask_RG_BPT"]], color=hex_c[1], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        ax.scatter(N2[diag["mask_SF_BPT"]], O3[diag["mask_SF_BPT"]], color=hex_c[2], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        ax.scatter(N2[diag["mask_wAGN_BPT"]], O3[diag["mask_wAGN_BPT"]], color=hex_c[3], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        ax.scatter(N2[diag["mask_sAGN_BPT"]], O3[diag["mask_sAGN_BPT"]], color=hex_c[4], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        title = r"BPT+WH$\alpha$"
    else:
        ax.scatter(N2[~diag["mask_ND"]], O3[~diag["mask_ND"]], color=hex_c[0], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        ax.scatter(N2[diag["mask_RG_new"]], O3[diag["mask_RG_new"]], color=hex_c[1], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        ax.scatter(N2[diag["mask_SF_new"]], O3[diag["mask_SF_new"]], color=hex_c[2], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        ax.scatter(N2[diag["mask_wAGN_new"]], O3[diag["mask_wAGN_new"]], color=hex_c[3], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        ax.scatter(N2[diag["mask_sAGN_new"]], O3[diag["mask_sAGN_new"]], color=hex_c[4], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
        title = "WHaD"

    ax.text(x_min + 0.05 * (x_max - x_min), y_min + 0.05 * (y_max - y_min),
            title, fontsize=21, bbox={"facecolor": "white", "pad": 5}, va="bottom", ha="left")
    ax.set_xlabel(r"log([NII]/H$\alpha$)")
    ax.set_ylabel(r"log([OIII]/H$\beta$)")


def plot_whad_panel(ax, diag: dict, cmap_t) -> None:
    hex_c = [mcolors.rgb2hex(cmap_t(i)[:3]) for i in range(cmap_t.N)]
    spx = 15
    x_min, x_max = 0.85, 2.55
    y_min, y_max = -0.49, 2.3

    ax.plot([x_min, 1.75], [0.78, 0.78], "--", c="black")
    ax.plot([1.75, 1.75], [0.47, y_max], "--", c="black")
    ax.plot([x_min, x_max], [0.47, 0.47], "--", c="black")
    ax.plot([2 - 0.243, x_max], [1.0, 1.0], "--", c="black")

    X = safe_log10(diag["img_disp_kms"])
    Y = safe_log10(np.abs(diag["img_WHa"]))

    my_contour(ax, X, Y, x_min, x_max, y_min, y_max, c_color="black", nbins=50, linewidths=1.5, alpha=0.7)
    ax.scatter(X[diag["mask_UK_new"]], Y[diag["mask_UK_new"]], color=hex_c[0], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
    ax.scatter(X[diag["mask_RG_new"]], Y[diag["mask_RG_new"]], color=hex_c[1], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
    ax.scatter(X[diag["mask_SF_new"]], Y[diag["mask_SF_new"]], color=hex_c[2], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
    ax.scatter(X[diag["mask_wAGN_new"]], Y[diag["mask_wAGN_new"]], color=hex_c[3], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)
    ax.scatter(X[diag["mask_sAGN_new"]], Y[diag["mask_sAGN_new"]], color=hex_c[4], alpha=ALPHA, edgecolor="None", s=spx, rasterized=True)

    ax.text(x_max - 0.2 * abs(x_max - x_min), y_min + 0.05 * abs(y_max - y_min), "Ret.")
    ax.text(x_min + 0.05 * abs(x_max - x_min), y_max - 0.15 * abs(y_max - y_min), "SF")
    ax.text(x_min + 0.05 * abs(x_max - x_min), 0.55, "Unk.")
    ax.text(x_max - 0.3 * abs(x_max - x_min), y_max - 0.15 * abs(y_max - y_min), "sAGN")
    ax.text(x_max - 0.3 * abs(x_max - x_min), 0.7, "wAGN")

    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.text(x_min + 0.05 * (x_max - x_min), y_min + 0.05 * (y_max - y_min),
            "WHaD", fontsize=21, bbox={"facecolor": "white", "pad": 5}, va="bottom", ha="left")
    ax.set_xlabel(r"log($\sigma$) [km s$^{-1}$]")
    ax.set_ylabel(r"log$|$EW(H$\alpha$)$|$ [${\rm \AA}$]")


def make_figure(name: str, input_dir: Path, output_dir: Path, usetex: bool = False) -> Path:
    configure_matplotlib(usetex=usetex)
    output_dir.mkdir(parents=True, exist_ok=True)

    hdu, cube_path = read_pipe3d_cube(input_dir, name)
    try:
        select_reg_data = hdu[8].data
        ssp_data = hdu[1].data
        ssp_hdr = hdu[1].header
        elines_data = hdu[5].data

        if elines_data.shape[0] <= 444:
            raise ValueError(
                f"Emission-line extension has only {elines_data.shape[0]} planes; "
                "this script expects the Pipe3D/KILOGAS indexing used in the notebook."
            )

        _, ny, nx = elines_data.shape
        diag = compute_diagnostics(elines_data, name)
        img_Ha_NII = diag["img_Ha_NII"]
        ha_min, levels = finite_contour_levels(img_Ha_NII)
        if not np.isfinite(ha_min):
            ha_min = -np.inf

        bg, bg_label = get_background_image(input_dir, name, fallback=ssp_data[0, :, :])
        bg = np.asarray(bg, dtype=float)
        image_PS = np.fliplr(np.flipud(bg))
        nx_v23 = image_PS.shape[1] * 0.5
        ny_v23 = image_PS.shape[0] * 0.5
        map_extent = [nx_v23 / 4, -nx_v23 / 4, -ny_v23 / 4, ny_v23 / 4]
        bg_extent = [nx_v23, -nx_v23, -ny_v23, ny_v23]

        fig, axes = plt.subplots(2, 3, figsize=(13, 7.2))
        cmap_t = cm.get_cmap("jet_r", 5)

        # Top-left: image background + Halpha+[NII] contours.
        ax = axes[0][0]
        with np.errstate(invalid="ignore"):
            bg_to_plot = np.sqrt(np.clip(image_PS, a_min=0, a_max=None))
        ax.imshow(bg_to_plot, extent=bg_extent, origin="lower", interpolation="none", cmap="Grays")
        ax.plot([15, 15, -15, -15, 15], [-15, 15, 15, -15, -15],
                color="orange", linestyle="-", linewidth=1.2)
        overlay_ha_contours(ax, img_Ha_NII, levels, map_extent, color="white", alpha=0.5 * ALPHA)
        ax.set_xticks([-40, -20, 0, 20, 40])
        ax.set_yticks([-40, -20, 0, 20, 40])
        ax.set_xlabel(r"$\Delta$ RA (arcsec)")
        ax.set_ylabel(r"$\Delta$ DEC (arcsec)")
        ax.text(0.04, 0.08, rf"{bg_label}, H$\alpha$+[NII]", fontsize=19,
                bbox={"facecolor": "white", "pad": 5}, va="bottom", ha="left",
                transform=ax.transAxes)

        plot_bpt_panel(axes[0][1], diag, cmap_t, mode="BPT")
        plot_diag_map(axes[0][2], fig, diag["img_diag_BPT"], select_reg_data,
                      img_Ha_NII, ha_min, levels, cmap_t, map_extent, r"BPT+WH$\alpha$")

        plot_whad_panel(axes[1][0], diag, cmap_t)
        plot_bpt_panel(axes[1][1], diag, cmap_t, mode="WHAD")
        plot_diag_map(axes[1][2], fig, diag["img_diag_new"], select_reg_data,
                      img_Ha_NII, ha_min, levels, cmap_t, map_extent, "WHaD")

        # Red circles at the nominal center, as in the notebook.
        for ax in (axes[0][2], axes[1][2]):
            circle = patches.Circle((0, 0), 1.5, ec="none")
            circle.set(color="None", edgecolor="red", linewidth=3, alpha=0.7)
            ax.add_artist(circle)

        fig.tight_layout(w_pad=0.3, h_pad=0.2)
        out_path = output_dir / f"{name}_diag.pdf"
        fig.savefig(out_path, transparent=False, facecolor="white", edgecolor="white")
        plt.close(fig)
        print(f"Input cube: {cube_path}")
        print(f"Output figure: {out_path}")
        return out_path
    finally:
        hdu.close()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create WHaD/BPT diagnostic figure from a Pipe3D KILOGAS cube."
    )
    parser.add_argument("--name", required=True, help="Object/cube name, e.g. KG-SAMI-63777")
    parser.add_argument("--input-dir", required=True, type=Path, help="Input directory containing the Pipe3D cube")
    parser.add_argument("--output-dir", required=True, type=Path, help="Output directory for the PDF figure")
    parser.add_argument("--usetex", action="store_true", help="Use LaTeX rendering for matplotlib text")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        make_figure(args.name, args.input_dir, args.output_dir, usetex=args.usetex)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
