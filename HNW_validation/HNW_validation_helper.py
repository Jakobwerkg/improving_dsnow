import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import pandas as pd
import xarray as xr





def _filter_season(df, full_season=False):
    """Filter DataFrame to snow season (Nov 1 – Apr 30) unless full_season=True."""
    
    df = _set_datime_index(df)

    if full_season:
        return df.copy()
    
    else: 
        mask = (df.index.month >= 11) | (df.index.month <= 4)
        return df.loc[mask].copy()

def _set_datime_index(df):
    """Ensure DataFrame has a DatetimeIndex."""
    
    df.index = pd.to_datetime(df.time)
    return df


def _calculate_metrics(x, y):
    """
    Calculate validation metrics.
    x = observed
    y = modeled
    """

    residuals = y - x

    rmse = np.sqrt(np.mean(residuals**2))
    bias = np.mean(residuals)
    pbias = np.sum(residuals) / np.sum(x) if np.sum(x) != 0 else np.nan

    ss_res = np.sum((x - y)**2)
    ss_tot = np.sum((x - np.mean(x))**2)
    r2 = 1 - ss_res / ss_tot if ss_tot != 0 else np.nan

    return {
        "RMSE": rmse,
        "Bias": bias,
        "Rel_BIAS": pbias,
        "R2": r2
    }


def _plot_validation(x, y, stats, model_name, lim, xlabel, ylabel):
    """
    Generic validation density plot.
    """

    vmax = max(1, len(x) / 10)

    plt.figure(figsize=(8, 7))

    plt.hist2d(
        x, y,
        bins=50,
        range=[lim, lim],
        norm=LogNorm(vmin=1, vmax=vmax),
        cmap="viridis"
    )

    cb = plt.colorbar(label="Number of observations")

    ticks = [t for t in [1, 10, 100, 1000, 10000] if t <= vmax]
    if vmax not in ticks:
        ticks.append(vmax)

    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{int(t)}" for t in ticks[:-1]] + [f"{int(vmax)}"])

    plt.plot(lim, lim, "--", color="gray", linewidth=1.3)

    ticks_xy = np.linspace(lim[0], lim[1], 5)
    plt.xticks(ticks_xy)
    plt.yticks(ticks_xy)

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(model_name, fontsize=15)

    textstr = (
        f"$R^2$: {stats['R2']:.2f}\n"
        f"Bias: {stats['Bias']:.2f}\n"
        f"RMSE: {stats['RMSE']:.1f}\n"
        f"Rel_BIAS: {stats['Rel_BIAS']:.1%}\n"
    )

    plt.text(
        0.03, 0.97, textstr,
        transform=plt.gca().transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    plt.xlim(lim)
    plt.ylim(lim)

    plt.grid(False)
    plt.tight_layout()


def _plot_validation_ax(ax, x, y, stats, title, lim, xlabel, ylabel, fig):
    """Draw a density validation plot into an existing Axes."""

    vmax = max(1, len(x) / 10)

    h, xedges, yedges, img = ax.hist2d(
        x, y,
        bins=50,
        range=[lim, lim],
        norm=LogNorm(vmin=1, vmax=vmax),
        cmap="viridis"
    )

    cb = fig.colorbar(img, ax=ax, label="Number of observations")

    ticks = [t for t in [1, 10, 100, 1000, 10000] if t <= vmax]
    if vmax not in ticks:
        ticks.append(vmax)

    cb.set_ticks(ticks)
    cb.set_ticklabels([f"{int(t)}" for t in ticks[:-1]] + [f"{int(vmax)}"])

    ax.plot(lim, lim, "--", color="gray", linewidth=1.3)

    ticks_xy = np.linspace(lim[0], lim[1], 5)
    ax.set_xticks(ticks_xy)
    ax.set_yticks(ticks_xy)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=13)

    textstr = (
        f"$R^2$: {stats['R2']:.2f}\n"
        f"Bias: {stats['Bias']:.2f}\n"
        f"RMSE: {stats['RMSE']:.1f}\n"
        f"Rel_BIAS: {stats['Rel_BIAS']:.1%}\n"
        f"N: {stats['N']}"
    )

    ax.text(
        0.03, 0.97, textstr,
        transform=ax.transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8)
    )

    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.grid(False)


def validate_hnw_mag25(df,
                       model_name,
                       obs_col="HNW_obs",
                       mod_col="HNW_mod",
                       save_dir=None,
                       filename="hnw_validation.png",
                       full_season=False,
                       drop_weisfluh_joch=True, ax=None):
    if drop_weisfluh_joch:
        df = df[df["station"] != "Weisfluh_Joch"].copy()

    df = _filter_season(df, full_season)

    df_valid = df.dropna(subset=[obs_col, mod_col])

    

    # Caution greate differnce if it is >= 0 or > 0
    df_valid = df_valid[df_valid[obs_col] >= 0]

    df_valid = df_valid[
        np.isfinite(df_valid[obs_col]) &
        np.isfinite(df_valid[mod_col])
    ]

    x = df_valid[obs_col].values
    y = df_valid[mod_col].values

    stats = _calculate_metrics(x, y)
    stats["N"] = len(df_valid)

    print(stats)

    _plot_validation(
        x,
        y,
        stats,
        model_name,
        lim=[0, 100],
        xlabel="Observed HNW (mm)",
        ylabel="Modeled HNW (mm)"
    )

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    plt.show()

    return stats


def validate_swe_mag25(df,
                       model_name,
                       obs_col="SWE_obs",
                       mod_col="SWE_mod",
                       save_dir=None,
                       filename="swe_validation.png",
                       full_season=False,
                       drop_weisfluh_joch=True):


    if drop_weisfluh_joch:
        df = df[df["station"] != "Weisfluh_Joch"].copy()

    df = _filter_season(df, full_season)

    df_valid = df.dropna(subset=[obs_col, mod_col])

    df_valid = df_valid[df_valid[obs_col] >= 0]

    print(f"Number of valid observations after filtering: {len(df_valid)}")

    df_valid = df_valid[
        np.isfinite(df_valid[obs_col]) &
        np.isfinite(df_valid[mod_col])
    ]

    x = df_valid[obs_col].values
    y = df_valid[mod_col].values

    stats = _calculate_metrics(x, y)
    stats["N"] = len(df_valid)

    print(stats)

    _plot_validation(
        x,
        y,
        stats,
        model_name,
        lim=[0, 1000],
        xlabel="Observed SWE (mm)",
        ylabel="Modeled SWE (mm)"
    )

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    plt.show()

    return stats


def validate_hnw_swe_combined(hnw_df, swe_df, model_name,
                               params=None,
                               hnw_obs_col="HNW_obs", hnw_mod_col="HNW_mod",
                               swe_obs_col="SWE_obs", swe_mod_col="SWE_mod",
                               full_season=False, drop_weisfluh_joch=True,
                               save_dir=None, filename="hnw_swe_validation_combined.png"):
    """
    Plot HNW and SWE validation side-by-side in one figure.

    Parameters
    ----------
    hnw_df : pd.DataFrame  — must contain hnw_obs_col, hnw_mod_col, station, time
    swe_df : pd.DataFrame  — must contain swe_obs_col, swe_mod_col, station, time
    model_name : str       — shown as figure suptitle
    params : dict, optional
        SnowToSwe parameters to annotate in the figure, e.g.
        dict(rho_max=381.203, rho_null=106.832, c_ov=0.00055347,
             k_ov=0.403141, k=0.0272148, tau=0.0222341, eta_null=8.65803e6)
    full_season : bool     — if False, restrict to Nov–Apr
    drop_weisfluh_joch : bool
    save_dir : str or Path, optional
    filename : str

    Returns
    -------
    dict with keys "HNW" and "SWE", each a stats dict
    (RMSE, Bias, Rel_BIAS, R2, N).
    """

    # ── HNW preparation ───────────────────────────────────────────────────────
    hnw = hnw_df.copy()
    if drop_weisfluh_joch:
        hnw = hnw[hnw["station"] != "Weisfluh_Joch"]
    hnw = _filter_season(hnw, full_season)
    hnw = hnw.dropna(subset=[hnw_obs_col, hnw_mod_col])
    hnw = hnw[hnw[hnw_obs_col] >= 0]
    hnw = hnw[np.isfinite(hnw[hnw_obs_col]) & np.isfinite(hnw[hnw_mod_col])]

    x_hnw = hnw[hnw_obs_col].values
    y_hnw = hnw[hnw_mod_col].values
    stats_hnw = _calculate_metrics(x_hnw, y_hnw)
    stats_hnw["N"] = len(hnw)

    # ── SWE preparation ───────────────────────────────────────────────────────
    swe = swe_df.copy()
    if drop_weisfluh_joch:
        swe = swe[swe["station"] != "Weisfluh_Joch"]
    swe = _filter_season(swe, full_season=True)   # SWE uses full year
    swe = swe.dropna(subset=[swe_obs_col, swe_mod_col])
    swe = swe[swe[swe_obs_col] >= 0]
    swe = swe[np.isfinite(swe[swe_obs_col]) & np.isfinite(swe[swe_mod_col])]

    x_swe = swe[swe_obs_col].values
    y_swe = swe[swe_mod_col].values
    stats_swe = _calculate_metrics(x_swe, y_swe)
    stats_swe["N"] = len(swe)

    print("HNW stats:", stats_hnw)
    print("SWE stats:", stats_swe)

    # ── Build params annotation string ────────────────────────────────────────
    param_str = ""
    if params is not None:
        parts = [f"{k}={v:.6g}" for k, v in params.items()]
        param_str = "  |  ".join(parts)

    # ── Figure ────────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))

    _plot_validation_ax(
        ax=axes[0], x=x_hnw, y=y_hnw, stats=stats_hnw,
        title="HNW validation",
        lim=[0, 100],
        xlabel="Observed HNW (mm)",
        ylabel="Modeled HNW (mm)",
        fig=fig
    )

    _plot_validation_ax(
        ax=axes[1], x=x_swe, y=y_swe, stats=stats_swe,
        title="SWE validation",
        lim=[0, 1000],
        xlabel="Observed SWE (mm)",
        ylabel="Modeled SWE (mm)",
        fig=fig
    )

    suptitle = model_name
    if param_str:
        suptitle += f"\n{param_str}"

    fig.suptitle(suptitle, fontsize=13, y=1.01)
    plt.tight_layout()

    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    plt.show()

    return {"HNW": stats_hnw, "SWE": stats_swe}

