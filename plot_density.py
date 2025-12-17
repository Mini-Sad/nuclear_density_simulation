#!/usr/bin/env python3
"""plot_density.py

Trace des plots 2D à partir des fichiers générés par main.cpp.

Le main C++ (par défaut) écrit dans le dossier `output/` :
  - grid_r.arma
  - grid_z.arma
  - density_optimized.arma
  - density_naive.arma

Ce script :
  - charge les fichiers Armadillo (arma_ascii),
  - calcule la différence (naive - optimized),
  - produit 3 heatmaps (3 figures séparées) : optimized / naive / diff,
  - affiche (terminal) + annote (figure) une estimation de la densité totale,
  - sauvegarde les images dans un répertoire `plots/` (créé si nécessaire).

Exemples :
  - Affichage interactif :
      python3 plot_density.py

  - Sauvegarde en PNG dans ./plots/ :
      python3 plot_density.py --save

  - Utiliser un autre dossier d'entrée :
      python3 plot_density.py --out-dir output_test --save
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, Tuple, Optional

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Lecture robuste des fichiers Armadillo ASCII (arma_ascii)
# ---------------------------------------------------------------------

def _read_next_nonempty_line(f) -> str:
    for line in f:
        s = line.strip()
        if s:
            return s
    raise ValueError("Fichier Armadillo ASCII incomplet : header/dimensions manquants")


def load_arma_ascii(path: os.PathLike | str) -> np.ndarray:
    """Charge un .arma sauvegardé en `arma::arma_ascii` (mat/vec).

    Format attendu :
      - 1 ligne header (ex: ARMA_MAT_TXT_FN008)
      - 1 ligne dimensions : n_rows n_cols
      - valeurs

    Retour : np.ndarray shape (n_rows, n_cols)
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Fichier introuvable : {path}")

    with path.open("r", encoding="utf-8", errors="replace") as f:
        _ = _read_next_nonempty_line(f)  # header
        dims = _read_next_nonempty_line(f)
        parts = dims.split()
        if len(parts) < 2:
            raise ValueError(f"Ligne dimensions invalide dans {path}: {dims!r}")
        n_rows, n_cols = int(parts[0]), int(parts[1])

        data = np.loadtxt(f)

    if n_rows == 0 or n_cols == 0:
        return np.zeros((n_rows, n_cols), dtype=float)

    arr = np.asarray(data)

    if arr.ndim == 0:
        arr = arr.reshape((1, 1))

    if arr.ndim == 1:
        if arr.size != n_rows * n_cols:
            raise ValueError(
                f"Taille data inattendue dans {path}: got {arr.size}, expected {n_rows*n_cols}"
            )
        arr = arr.reshape((n_rows, n_cols))

    if arr.shape != (n_rows, n_cols):
        flat = np.ravel(arr)
        if flat.size != n_rows * n_cols:
            raise ValueError(
                f"Shape inattendue dans {path}: got {arr.shape}, expected {(n_rows, n_cols)}"
            )
        arr = flat.reshape((n_rows, n_cols))

    return arr


def load_arma_ascii_vec(path: os.PathLike | str) -> np.ndarray:
    """Charge un vecteur Armadillo (arma_ascii) et renvoie un array 1D."""
    mat = load_arma_ascii(path)
    return np.squeeze(mat)


def _clip_small_negative_noise(a: np.ndarray, eps: float = 1e-14) -> Tuple[np.ndarray, bool]:
    """Clip les très petites valeurs négatives (bruit numérique) à 0 pour l'affichage."""
    a2 = np.array(a, copy=True)
    mask = (a2 < 0.0) & (a2 > -eps)
    changed = bool(np.any(mask))
    a2[mask] = 0.0
    return a2, changed


def make_symmetric_density(density: np.ndarray, r_vals: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Crée une vue symétrique de la densité (miroir par rapport à r=0).
    
    La densité nucléaire étant axisymétrique, on peut créer une vue complète
    en miroir pour r allant de -rmax à +rmax.
    
    Retour: (density_sym, r_sym) avec r_sym allant de -rmax à +rmax
    """
    # Créer le miroir de r (partie négative)
    r_neg = -r_vals[::-1]  # Inverser et rendre négatif
    
    # Créer la grille r symétrique: [-rmax, ..., -dr, 0, dr, ..., rmax]
    # On évite de dupliquer r=0 si présent
    if np.abs(r_vals[0]) < 1e-12:  # r commence à 0
        r_sym = np.concatenate([r_neg[:-1], r_vals])
        # Miroir de la densité (sans dupliquer la colonne r=0)
        density_sym = np.concatenate([density[:, ::-1][:, :-1], density], axis=1)
    else:
        r_sym = np.concatenate([r_neg, r_vals])
        density_sym = np.concatenate([density[:, ::-1], density], axis=1)
    
    return density_sym, r_sym


# ---------------------------------------------------------------------
# Métriques (densité totale)
# ---------------------------------------------------------------------

def estimate_total_density_rz(density: np.ndarray, r_vals: np.ndarray, z_vals: np.ndarray) -> float:
    """Estimation de ∬ rho(r,z) dr dz via intégration trapézoïdale."""
    # density.shape = (nbZ, nbR) avec z en lignes, r en colonnes.
    int_over_r = np.trapz(density, x=r_vals, axis=1)
    return float(np.trapz(int_over_r, x=z_vals, axis=0))


def estimate_total_density_axisymmetric(density: np.ndarray, r_vals: np.ndarray, z_vals: np.ndarray) -> float:
    """Estimation de la masse 3D axisymétrique : 2π ∬ rho(r,z) r dr dz."""
    weighted = density * r_vals[None, :]
    int_over_r = np.trapz(weighted, x=r_vals, axis=1)
    return float(2.0 * np.pi * np.trapz(int_over_r, x=z_vals, axis=0))


# ---------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------

def _format_info_lines(lines: Iterable[str]) -> str:
    clean = [str(x).strip() for x in lines if str(x).strip()]
    return "\n".join(clean)


def plot_combined(
    density_opt: np.ndarray,
    density_naive: np.ndarray,
    density_diff: np.ndarray,
    r_vals: np.ndarray,
    z_vals: np.ndarray,
    vmin: float,
    vmax: float,
    diff_abs: float,
    total_opt: float,
    total_naive: float,
    diff_l2: float,
    out_dir: Optional[Path],
    zoom_out_frac: float = 0.15,
    r_label: str = r"$r$ (fm)",
) -> None:
    """Crée une figure combinée avec les 3 heatmaps côte à côte."""
    extent = [float(r_vals.min()), float(r_vals.max()), float(z_vals.min()), float(z_vals.max())]
    x0, x1, y0, y1 = extent
    dx = (x1 - x0) * zoom_out_frac
    dy = (y1 - y0) * zoom_out_frac

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    fig.suptitle("Comparaison des Densités Nucléaires 2D", fontsize=16, fontweight='bold', y=1.02)

    # Plot 1: Optimisé
    im1 = axes[0].imshow(density_opt, extent=extent, origin="lower", aspect="auto",
                          cmap="viridis", vmin=vmin, vmax=vmax, interpolation='bilinear')
    axes[0].set_title("Algorithme Optimisé", fontsize=12, fontweight='bold')
    axes[0].set_xlabel(r_label, fontsize=11)
    axes[0].set_ylabel(r"$z$ (fm)", fontsize=11)
    axes[0].set_xlim(x0 - dx, x1 + dx)
    axes[0].set_ylim(y0 - dy, y1 + dy)
    axes[0].grid(True, linestyle='--', alpha=0.3)
    axes[0].axvline(x=0, color='white', linestyle='-', linewidth=0.5, alpha=0.5)
    cbar1 = plt.colorbar(im1, ax=axes[0], shrink=0.8)
    cbar1.set_label(r"$\rho$", fontsize=10)
    axes[0].text(0.02, 0.98, f"∬ρ = {total_opt:.2e}", transform=axes[0].transAxes,
                 va='top', fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Plot 2: Naïf
    im2 = axes[1].imshow(density_naive, extent=extent, origin="lower", aspect="auto",
                          cmap="viridis", vmin=vmin, vmax=vmax, interpolation='bilinear')
    axes[1].set_title("Algorithme Naïf", fontsize=12, fontweight='bold')
    axes[1].set_xlabel(r_label, fontsize=11)
    axes[1].set_ylabel(r"$z$ (fm)", fontsize=11)
    axes[1].set_xlim(x0 - dx, x1 + dx)
    axes[1].set_ylim(y0 - dy, y1 + dy)
    axes[1].grid(True, linestyle='--', alpha=0.3)
    axes[1].axvline(x=0, color='white', linestyle='-', linewidth=0.5, alpha=0.5)
    cbar2 = plt.colorbar(im2, ax=axes[1], shrink=0.8)
    cbar2.set_label(r"$\rho$", fontsize=10)
    axes[1].text(0.02, 0.98, f"∬ρ = {total_naive:.2e}", transform=axes[1].transAxes,
                 va='top', fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Plot 3: Différence
    diff_vmax = diff_abs if diff_abs > 1e-15 else 1e-15
    im3 = axes[2].imshow(density_diff, extent=extent, origin="lower", aspect="auto",
                          cmap="RdBu_r", vmin=-diff_vmax, vmax=diff_vmax, interpolation='bilinear')
    axes[2].set_title("Différence (Naïf - Optimisé)", fontsize=12, fontweight='bold')
    axes[2].set_xlabel(r_label, fontsize=11)
    axes[2].set_ylabel(r"$z$ (fm)", fontsize=11)
    axes[2].set_xlim(x0 - dx, x1 + dx)
    axes[2].set_ylim(y0 - dy, y1 + dy)
    axes[2].grid(True, linestyle='--', alpha=0.3)
    axes[2].axvline(x=0, color='gray', linestyle='-', linewidth=0.5, alpha=0.5)
    cbar3 = plt.colorbar(im3, ax=axes[2], shrink=0.8)
    cbar3.set_label(r"$\Delta\rho$", fontsize=10)
    status = "✓ OK" if diff_abs < 1e-10 else "⚠"
    axes[2].text(0.02, 0.98, f"||Δρ||₂ = {diff_l2:.2e}\n{status}", transform=axes[2].transAxes,
                 va='top', fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    if out_dir is not None:
        out_path = out_dir / "density_combined.png"
        out_dir.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
        print(f" - {out_path}")


def plot_heatmap(
    data: np.ndarray,
    r_vals: np.ndarray,
    z_vals: np.ndarray,
    title: str,
    source_label: str,
    cbar_label: str,
    out_path: Optional[Path],
    cmap: Optional[str],
    vmin: Optional[float],
    vmax: Optional[float],
    zoom_out_frac: float = 0.15,
    info_lines: Iterable[str] = (),
    show_contours: bool = True,
    show_grid: bool = True,
    r_label: str = r"$r$ (fm)",
) -> None:
    """Affiche une heatmap avec options améliorées pour la visualisation."""
    extent = [float(r_vals.min()), float(r_vals.max()), float(z_vals.min()), float(z_vals.max())]

    # Figure plus grande pour meilleure lisibilité
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Titre principal avec style amélioré
    ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
    ax.set_xlabel(r_label, fontsize=12)
    ax.set_ylabel(r"$z$ (fm)", fontsize=12)

    # Image principale avec interpolation pour lisser l'affichage
    im = ax.imshow(
        data,
        extent=extent,
        origin="lower",
        aspect="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        interpolation='bilinear',
    )
    
    # Colorbar améliorée
    cbar = plt.colorbar(im, ax=ax, label=cbar_label, shrink=0.85, pad=0.02)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(cbar_label, fontsize=11)

    # Contours pour mieux visualiser les niveaux de densité
    if show_contours and np.max(np.abs(data)) > 1e-15:
        try:
            # Créer grilles pour contour
            R, Z = np.meshgrid(r_vals, z_vals)
            # Niveaux de contour automatiques
            levels = np.linspace(np.min(data), np.max(data), 8)
            # Filtrer les niveaux trop proches de 0
            levels = levels[np.abs(levels) > 1e-15]
            if len(levels) > 2:
                contour = ax.contour(R, Z, data, levels=levels, colors='white', 
                                     linewidths=0.5, alpha=0.6)
                ax.clabel(contour, inline=True, fontsize=7, fmt='%.2e')
        except Exception:
            pass  # Ignorer les erreurs de contour

    # Grille pour mieux repérer les positions
    if show_grid:
        ax.grid(True, linestyle='--', alpha=0.3, color='gray')

    # Zoom out : on élargit les limites des axes pour voir clairement les bords
    x0, x1, y0, y1 = extent
    dx = (x1 - x0) * float(zoom_out_frac)
    dy = (y1 - y0) * float(zoom_out_frac)
    ax.set_xlim(x0 - dx, x1 + dx)
    ax.set_ylim(y0 - dy, y1 + dy)

    # Améliorer les ticks
    ax.tick_params(axis='both', which='major', labelsize=10)
    
    # Ajouter des marqueurs pour les limites du domaine
    ax.axhline(y=y0, color='black', linestyle='-', linewidth=0.8, alpha=0.5)
    ax.axhline(y=y1, color='black', linestyle='-', linewidth=0.8, alpha=0.5)
    ax.axvline(x=x0, color='black', linestyle='-', linewidth=0.8, alpha=0.5)
    ax.axvline(x=x1, color='black', linestyle='-', linewidth=0.8, alpha=0.5)
    
    # Axe de symétrie à r=0 si données symétriques
    if x0 < 0 < x1:
        ax.axvline(x=0, color='white', linestyle='--', linewidth=1.0, alpha=0.7)

    # Box d'information avec style amélioré
    info = _format_info_lines(info_lines)
    if info:
        ax.text(
            0.02,
            0.98,
            info,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=10,
            fontfamily='monospace',
            bbox=dict(boxstyle="round,pad=0.5", facecolor="white", 
                      alpha=0.85, edgecolor='gray', linewidth=1),
        )

    # Source en bas de la figure
    fig.text(0.5, 0.01, f"Source: {source_label}", ha='center', fontsize=8, 
             style='italic', color='gray')

    plt.tight_layout(rect=[0, 0.03, 1, 1])

    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')


def main() -> int:
    parser = argparse.ArgumentParser(description="Plots 2D des densités (output Armadillo).")
    parser.add_argument(
        "--out-dir",
        default="output",
        help="Dossier contenant grid_r.arma, grid_z.arma, density_optimized.arma, density_naive.arma (défaut: output)",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Sauvegarde les figures en PNG dans le dossier ./plots/ (créé si nécessaire).",
    )
    parser.add_argument(
        "--plots-dir",
        default="plots",
        help="Répertoire où écrire les PNG si --save (défaut: plots).",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Ne pas afficher la fenêtre matplotlib (utile en batch si --save).",
    )
    parser.add_argument(
        "--no-clip",
        action="store_true",
        help="Ne pas clipper les très petites valeurs négatives (bruit numérique) à 0 pour l'affichage.",
    )
    parser.add_argument(
        "--zoom-out",
        type=float,
        default=0.15,
        help="Fraction (ex: 0.15 = 15%%) pour élargir les limites des axes et voir toute la densité (défaut: 0.15).",
    )
    parser.add_argument(
        "--no-contours",
        action="store_true",
        help="Désactiver les lignes de contour sur les heatmaps.",
    )
    parser.add_argument(
        "--no-grid",
        action="store_true",
        help="Désactiver la grille de fond sur les heatmaps.",
    )
    parser.add_argument(
        "--combined",
        action="store_true",
        help="Générer en plus une figure combinée avec les 3 plots côte à côte.",
    )
    parser.add_argument(
        "--no-symmetric",
        action="store_true",
        help="Ne pas créer la vue symétrique (miroir). Par défaut, on affiche la densité complète de -rmax à +rmax.",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)

    r_path = out_dir / "grid_r.arma"
    z_path = out_dir / "grid_z.arma"
    dens_opt_path = out_dir / "density_optimized.arma"
    dens_naive_path = out_dir / "density_naive.arma"

    print("Chargement des grilles et densités depuis:", out_dir.resolve())

    r_vals = load_arma_ascii_vec(r_path)
    z_vals = load_arma_ascii_vec(z_path)
    density_opt = load_arma_ascii(dens_opt_path)
    density_naive = load_arma_ascii(dens_naive_path)

    if density_opt.ndim != 2 or density_naive.ndim != 2:
        raise ValueError("Les densités doivent être des matrices 2D")

    nbZ, nbR = density_opt.shape
    if density_naive.shape != (nbZ, nbR):
        raise ValueError(
            f"Dimensions incohérentes: optimized={density_opt.shape}, naive={density_naive.shape}"
        )

    if r_vals.shape[0] != nbR:
        raise ValueError(f"Taille de grid_r incompatible: grid_r={r_vals.shape[0]} vs nbR={nbR}")
    if z_vals.shape[0] != nbZ:
        raise ValueError(f"Taille de grid_z incompatible: grid_z={z_vals.shape[0]} vs nbZ={nbZ}")

    density_diff = density_naive - density_opt

    if not args.no_clip:
        density_opt_disp, c1 = _clip_small_negative_noise(density_opt)
        density_naive_disp, c2 = _clip_small_negative_noise(density_naive)
        if c1 or c2:
            print(
                "Note: petites valeurs négatives (bruit numérique) clippées à 0 pour l'affichage "
                "(utilisez --no-clip pour désactiver)."
            )
    else:
        density_opt_disp = density_opt
        density_naive_disp = density_naive

    # Densité totale (estimations) : affichée dans le terminal + annotée sur les figures.
    total_opt_rz = estimate_total_density_rz(density_opt_disp, r_vals, z_vals)
    total_naive_rz = estimate_total_density_rz(density_naive_disp, r_vals, z_vals)
    total_opt_3d = estimate_total_density_axisymmetric(density_opt_disp, r_vals, z_vals)
    total_naive_3d = estimate_total_density_axisymmetric(density_naive_disp, r_vals, z_vals)

    print("--- Totaux (estimations) ---")
    print(f"Optimized  ∬rho dr dz      = {total_opt_rz:.6e}")
    print(f"Naive      ∬rho dr dz      = {total_naive_rz:.6e}")
    print(f"Optimized  2π ∬rho r dr dz = {total_opt_3d:.6e}")
    print(f"Naive      2π ∬rho r dr dz = {total_naive_3d:.6e}")

    # Créer les vues symétriques (miroir) pour affichage complet
    if not args.no_symmetric:
        print("Création des vues symétriques (densité complète de -rmax à +rmax)...")
        density_opt_disp, r_vals_plot = make_symmetric_density(density_opt_disp, r_vals)
        density_naive_disp, _ = make_symmetric_density(density_naive_disp, r_vals)
        density_diff_sym, _ = make_symmetric_density(density_diff, r_vals)
        r_label = r"$r$ (fm)"
    else:
        r_vals_plot = r_vals
        density_diff_sym = density_diff
        r_label = r"$r_\perp$ (fm)"

    # Pour comparer visuellement opt vs naive, on fixe une échelle commune.
    common_vmin = 0.0
    common_vmax = float(max(np.max(density_opt_disp), np.max(density_naive_disp)))

    plots_dir = Path(args.plots_dir)

    # Chemins de sortie (toujours dans plots/ si --save)
    opt_png = (plots_dir / "density_optimized.png") if args.save else None
    naive_png = (plots_dir / "density_naive.png") if args.save else None
    diff_png = (plots_dir / "density_diff.png") if args.save else None

    plot_heatmap(
        density_opt_disp,
        r_vals_plot,
        z_vals,
        title="Densité Nucléaire 2D — Algorithme Optimisé",
        source_label=str(dens_opt_path),
        cbar_label=r"$\rho(r, z)$",
        out_path=opt_png,
        cmap="viridis",
        vmin=common_vmin,
        vmax=common_vmax,
        zoom_out_frac=args.zoom_out,
        info_lines=(
            f"Total ∬ρ dr dz = {total_opt_rz:.3e}",
            f"Total 2π ∬ρ r dr dz = {total_opt_3d:.3e}",
            f"Grille: {nbR} × {nbZ} points",
        ),
        show_contours=not args.no_contours,
        show_grid=not args.no_grid,
        r_label=r_label,
    )

    plot_heatmap(
        density_naive_disp,
        r_vals_plot,
        z_vals,
        title="Densité Nucléaire 2D — Algorithme Naïf",
        source_label=str(dens_naive_path),
        cbar_label=r"$\rho(r, z)$",
        out_path=naive_png,
        cmap="viridis",
        vmin=common_vmin,
        vmax=common_vmax,
        zoom_out_frac=args.zoom_out,
        info_lines=(
            f"Total ∬ρ dr dz = {total_naive_rz:.3e}",
            f"Total 2π ∬ρ r dr dz = {total_naive_3d:.3e}",
            f"Grille: {nbR} × {nbZ} points",
        ),
        show_contours=not args.no_contours,
        show_grid=not args.no_grid,
        r_label=r_label,
    )

    # Différence: diverging pour voir + / -
    diff_abs = float(np.max(np.abs(density_diff)))
    diff_l2 = float(np.linalg.norm(density_diff))
    plot_heatmap(
        density_diff_sym,
        r_vals_plot,
        z_vals,
        title="Différence 2D — (Naïf - Optimisé)",
        source_label=f"{dens_naive_path} - {dens_opt_path}",
        cbar_label=r"$\Delta\rho$",
        out_path=diff_png,
        cmap="RdBu_r",
        vmin=-diff_abs if diff_abs > 1e-15 else -1e-15,
        vmax=diff_abs if diff_abs > 1e-15 else 1e-15,
        zoom_out_frac=args.zoom_out,
        info_lines=(
            f"max |Δρ| = {diff_abs:.3e}",
            f"||Δρ||₂ = {diff_l2:.3e}",
            f"Validation: {'✓ OK' if diff_abs < 1e-10 else '⚠ Écart'}",
        ),
        show_contours=False,  # Contours peu utiles pour diff proche de 0
        show_grid=not args.no_grid,
        r_label=r_label,
    )

    # Figure combinée optionnelle
    if args.combined:
        plot_combined(
            density_opt_disp, density_naive_disp, density_diff_sym,
            r_vals_plot, z_vals,
            common_vmin, common_vmax, diff_abs,
            total_opt_rz, total_naive_rz, diff_l2,
            plots_dir if args.save else None,
            args.zoom_out,
            r_label=r_label,
        )

    if args.save:
        plots_dir.mkdir(parents=True, exist_ok=True)
        print("PNG sauvegardés dans:", plots_dir.resolve())
        print(" -", opt_png)
        print(" -", naive_png)
        print(" -", diff_png)

    if not args.no_show:
        plt.show()
    else:
        plt.close("all")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
