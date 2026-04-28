from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt


# ============================================================================
# CONFIG
# ============================================================================

BASE_DIR = Path("/Users/jakobwerkgarner/code/mt_dsnow/model_diff/R_comparisson/output")
DSNOW_DIR = BASE_DIR / "HS_SWE_by_station" / "HS_SWE_by_station_dsnow_only"
HS2SWE_DIR = BASE_DIR / "HS_SWE_by_station" / "HS_SWE_by_station_hs2swe_only"

OBS_DIR = Path(
    "/Users/jakobwerkgarner/code/mt_dsnow/calibration/calibration_data/output/HS_SWE_by_station"
)

OUTPUT_DIR = BASE_DIR / "batch_comparison_output"
PLOT_DIR = OUTPUT_DIR / "plots"
NC_OUT_DIR = OUTPUT_DIR / "netcdf_with_obs"

SWE_MEAN_DIR = PLOT_DIR / "swe_mean"
SWE_SEASON_DIR = PLOT_DIR / "swe_by_season"
SWE_SCATTER_DIR = PLOT_DIR / "swe_scatter"
FIRST_APP_SUMMARY_DIR = PLOT_DIR / "first_appearance_summary"
FIRST_APP_SEASON_DIR = PLOT_DIR / "first_appearance_by_season"
FIRST_APP_ALL_STATIONS_DIR = PLOT_DIR / "first_appearance_all_stations"

for d in [
    OUTPUT_DIR,
    PLOT_DIR,
    NC_OUT_DIR,
    SWE_MEAN_DIR,
    SWE_SEASON_DIR,
    SWE_SCATTER_DIR,
    FIRST_APP_SUMMARY_DIR,
    FIRST_APP_SEASON_DIR,
    FIRST_APP_ALL_STATIONS_DIR,
]:
    d.mkdir(parents=True, exist_ok=True)


# ============================================================================
# FILE MATCHING
# ============================================================================

def station_from_filename(path: Path, suffix: str) -> str:
    return path.name.replace(suffix, "").replace(".nc", "")


def build_station_file_table(
    dsnow_dir: Path,
    hs2swe_dir: Path,
    obs_dir: Path,
) -> pd.DataFrame:
    dsnow_files = sorted(dsnow_dir.glob("*_dsnow_allseasons.nc"))
    hs2swe_files = sorted(hs2swe_dir.glob("*_hs2swe_allseasons.nc"))

    ds_map = {
        station_from_filename(f, "_dsnow_allseasons"): f
        for f in dsnow_files
    }
    hs2_map = {
        station_from_filename(f, "_hs2swe_allseasons"): f
        for f in hs2swe_files
    }

    stations = sorted(set(ds_map) & set(hs2_map))

    rows = []
    for station in stations:
        obs_file = obs_dir / f"{station}_hs_swe_obs.csv"
        rows.append(
            {
                "station": station,
                "dsnow_file": ds_map[station],
                "hs2swe_file": hs2_map[station],
                "obs_file": obs_file,
                "obs_exists": obs_file.exists(),
            }
        )

    return pd.DataFrame(rows)


# ============================================================================
# OBSERVATION HANDLING
# ============================================================================

def load_and_prepare_observations(
    obs_file: Path,
    target_seasons: np.ndarray,
    target_dos: np.ndarray,
) -> pd.DataFrame:
    obs = pd.read_csv(obs_file).copy()
    obs["date"] = pd.to_datetime(obs["date"], errors="coerce")
    obs = obs.dropna(subset=["date"])

    obs["season"] = np.where(
        obs["date"].dt.month >= 11,
        obs["date"].dt.year + 1,
        obs["date"].dt.year,
    ).astype(int)

    season_start = pd.to_datetime((obs["season"] - 1).astype(str) + "-11-01")
    obs["dos"] = ((obs["date"] - season_start).dt.days + 1).astype(int)

    obs = obs[
        obs["season"].isin(target_seasons)
        & obs["dos"].between(1, int(np.max(target_dos)))
    ].copy()

    obs_daily = (
        obs.groupby(["season", "dos"], as_index=False)[["hs", "swe_obs"]]
        .mean()
    )

    return obs_daily


def observations_to_xarray(
    obs_daily: pd.DataFrame,
    seasons: np.ndarray,
    dos: np.ndarray,
) -> tuple[xr.DataArray, xr.DataArray]:
    idx = pd.MultiIndex.from_product([seasons, dos], names=["season", "dos"])
    grid = pd.DataFrame(index=idx).join(obs_daily.set_index(["season", "dos"]))

    hs_obs_m = xr.DataArray(
        grid["hs"].to_numpy().reshape(len(seasons), len(dos)),
        coords={"season": seasons, "dos": dos},
        dims=("season", "dos"),
        name="HS_obs_m",
    )

    swe_obs_mm = xr.DataArray(
        grid["swe_obs"].to_numpy().reshape(len(seasons), len(dos)),
        coords={"season": seasons, "dos": dos},
        dims=("season", "dos"),
        name="SWE_obs_mm",
    )

    return hs_obs_m, swe_obs_mm


def add_observations_to_datasets(
    ds: xr.Dataset,
    hs2: xr.Dataset,
    hs_obs_m_hs2: xr.DataArray,
    swe_obs_mm_hs2: xr.DataArray,
) -> tuple[xr.Dataset, xr.Dataset]:
    ds = ds.copy()
    hs2 = hs2.copy()

    hs2["HS_obs_cm"] = hs_obs_m_hs2 * 100.0
    hs2["HS_obs_cm"].attrs["units"] = "cm"

    hs2["SWE_obs_mm"] = swe_obs_mm_hs2
    hs2["SWE_obs_mm"].attrs["units"] = "mm"

    ds["HS_meas_m"] = (
        hs_obs_m_hs2.reindex(season=ds["season"], dos=ds["dos"])
        .transpose("dos", "season")
    )
    ds["HS_meas_m"].attrs["units"] = "m"

    ds["SWE_obs_mm"] = (
        swe_obs_mm_hs2.reindex(season=ds["season"], dos=ds["dos"])
        .transpose("dos", "season")
    )
    ds["SWE_obs_mm"].attrs["units"] = "mm"

    return ds, hs2


# ============================================================================
# MODEL SERIES
# ============================================================================

def prepare_model_series(ds: xr.Dataset, hs2: xr.Dataset) -> dict[str, xr.DataArray]:
    swe_ds = ds["SWE_total"].transpose("season", "dos")
    swe_hs2 = hs2["SWE_mm"].transpose("season", "dos")
    swe_ds, swe_hs2 = xr.align(swe_ds, swe_hs2, join="inner")

    hs_ds = ds["HS"].sum("layer").transpose("season", "dos")
    hs_hs2 = (hs2["HS_layer_cm"].sum("layer") / 100.0).transpose("season", "dos")
    hs_ds, hs_hs2 = xr.align(hs_ds, hs_hs2, join="inner")

    return {
        "swe_ds": swe_ds,
        "swe_hs2": swe_hs2,
        "hs_ds": hs_ds,
        "hs_hs2": hs_hs2,
    }


def running_mean(da: xr.DataArray, window: int = 7) -> xr.DataArray:
    return da.rolling(dos=window, center=True, min_periods=1).mean()


# ============================================================================
# FIRST APPEARANCE
# ============================================================================

def first_appearance_events(
    layer_da: xr.DataArray,
    scale_to_m: float = 1.0,
) -> pd.DataFrame:
    da = (layer_da * scale_to_m).transpose("season", "dos", "layer")

    pos = (da > 0) & np.isfinite(da)
    has_pos = pos.any("dos")
    first_idx = pos.argmax("dos")

    first_dos = da["dos"].isel(dos=first_idx).where(has_pos).rename("first_dos")
    first_inc = da.isel(dos=first_idx).where(has_pos).rename("inc").reset_coords(drop=True)

    ev = xr.Dataset({"first_dos": first_dos, "inc": first_inc}).to_dataframe().reset_index()
    ev = ev.dropna(subset=["first_dos", "inc"]).copy()

    ev["season"] = ev["season"].astype(int)
    ev["first_dos"] = ev["first_dos"].astype(int)
    ev["timestamp"] = (
        pd.to_datetime((ev["season"] - 1).astype(str) + "-11-01")
        + pd.to_timedelta(ev["first_dos"] - 1, unit="D")
    )

    return ev[["season", "timestamp", "inc"]]


def season_curve(events: pd.DataFrame, season: int) -> pd.DataFrame:
    out = (
        events.loc[events["season"] == season]
        .groupby("timestamp", as_index=True)["inc"]
        .sum()
        .sort_index()
        .to_frame("increment")
    )
    out["cumulative"] = out["increment"].cumsum()
    return out


def cumulative_matrix_from_events(
    events: pd.DataFrame,
    seasons: np.ndarray,
    dos: np.ndarray,
) -> pd.DataFrame:
    e = events.copy()
    e["season_start"] = pd.to_datetime((e["season"] - 1).astype(str) + "-11-01")
    e["dos"] = ((e["timestamp"] - e["season_start"]).dt.days + 1).astype(int)

    daily = e.groupby(["season", "dos"], as_index=False)["inc"].sum()

    idx = pd.MultiIndex.from_product([seasons, dos], names=["season", "dos"])
    inc = (
        daily.set_index(["season", "dos"])["inc"]
        .reindex(idx, fill_value=0.0)
        .unstack("season")
    )

    return inc.cumsum(axis=0)


def curve_stats(cum_mat: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=cum_mat.index)
    out["mean"] = cum_mat.mean(axis=1)
    out["median"] = cum_mat.median(axis=1)
    out["std"] = cum_mat.std(axis=1, ddof=1)
    out["q25"] = cum_mat.quantile(0.25, axis=1)
    out["q75"] = cum_mat.quantile(0.75, axis=1)
    out["iqr"] = out["q75"] - out["q25"]
    return out


def final_summary(cum_mat: pd.DataFrame) -> pd.Series:
    vals = cum_mat.iloc[-1].to_numpy(dtype=float)
    q25 = np.quantile(vals, 0.25)
    q75 = np.quantile(vals, 0.75)
    return pd.Series(
        {
            "n": len(vals),
            "mean": np.mean(vals),
            "median": np.median(vals),
            "q25": q25,
            "q75": q75,
            "iqr": q75 - q25,
            "std": np.std(vals, ddof=1),
        }
    )


# ============================================================================
# PLOTTING
# ============================================================================

def save_swe_mean_plot_with_obs(
    station: str,
    swe_ds: xr.DataArray,
    swe_hs2: xr.DataArray,
    swe_obs: xr.DataArray,
    out_file: Path,
    window: int = 7,
) -> None:
    swe_ds_7d = running_mean(swe_ds, window)
    swe_hs2_7d = running_mean(swe_hs2, window)

    swe_ds_mean = swe_ds_7d.mean("season", skipna=True)
    swe_ds_std = swe_ds_7d.std("season", skipna=True)

    swe_hs2_mean = swe_hs2_7d.mean("season", skipna=True)
    swe_hs2_std = swe_hs2_7d.std("season", skipna=True)

    swe_obs_mean = swe_obs.mean("season", skipna=True)

    dos = swe_ds_7d["dos"].values
    obs_vals = swe_obs_mean.values
    mask_obs = np.isfinite(obs_vals)

    fig, ax = plt.subplots(figsize=(12, 7), constrained_layout=True)

    ax.plot(dos, swe_ds_mean.values, label="delta-snow", color="C0", lw=2)
    ax.plot(dos, swe_hs2_mean.values, label="hs2swe", color="C1", lw=2)

    ax.fill_between(
        dos,
        (swe_ds_mean - swe_ds_std).values,
        (swe_ds_mean + swe_ds_std).values,
        color="C0",
        alpha=0.2,
        label="delta-snow ±1 std",
    )
    ax.fill_between(
        dos,
        (swe_hs2_mean - swe_hs2_std).values,
        (swe_hs2_mean + swe_hs2_std).values,
        color="C1",
        alpha=0.2,
        label="hs2swe ±1 std",
    )

    ax.scatter(
        dos[mask_obs],
        obs_vals[mask_obs],
        s=25,
        c="k",
        label="obs (mean scatter)",
        zorder=4,
    )

    ax.set_title(f"{station}: Mean SWE over seasons (7-day running mean)")
    ax.set_xlabel("DOS")
    ax.set_ylabel("SWE [mm]")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.savefig(out_file, dpi=200)
    plt.close(fig)


def save_swe_by_season_plot(
    station: str,
    swe_ds: xr.DataArray,
    swe_hs2: xr.DataArray,
    swe_obs: xr.DataArray,
    out_file: Path,
) -> None:
    seasons = swe_ds["season"].values
    n = len(seasons)

    fig, axes = plt.subplots(
        n, 1, figsize=(13, max(3.0 * n, 4.0)), sharex=True, constrained_layout=True
    )

    if n == 1:
        axes = np.array([axes])

    for i, season in enumerate(seasons):
        ax = axes[i]

        ds_season = swe_ds.sel(season=season)
        hs2_season = swe_hs2.sel(season=season)
        obs_season = swe_obs.sel(season=season)

        mask_obs = np.isfinite(obs_season.values)

        ax.plot(ds_season["dos"], ds_season.values, lw=1.8, label="delta-snow", color="C0")
        ax.plot(hs2_season["dos"], hs2_season.values, lw=1.8, label="hs2swe", color="C1")

        ax.scatter(
            obs_season["dos"].values[mask_obs],
            obs_season.values[mask_obs],
            s=18,
            c="k",
            label="obs",
            zorder=4,
        )

        ax.set_ylabel(f"SWE [mm]\n{int(season)}")
        ax.grid(True, alpha=0.3)

        if i == 0:
            ax.set_title(f"{station}: SWE by season")
            ax.legend()

    axes[-1].set_xlabel("DOS")
    fig.savefig(out_file, dpi=200)
    plt.close(fig)


def save_swe_scatter_plot(
    station: str,
    swe_obs: xr.DataArray,
    swe_ds: xr.DataArray,
    swe_hs2: xr.DataArray,
    out_file: Path,
) -> dict[str, float]:
    swe_obs_al, swe_ds_al, swe_hs2_al = xr.align(swe_obs, swe_ds, swe_hs2, join="inner")

    obs_flat = swe_obs_al.values.ravel()
    ds_flat = swe_ds_al.values.ravel()
    hs2_flat = swe_hs2_al.values.ravel()

    def _metrics(obs_arr, mod_arr):
        mask = np.isfinite(obs_arr) & np.isfinite(mod_arr)
        o = obs_arr[mask]
        m = mod_arr[mask]

        if len(o) == 0:
            return mask, np.nan, np.nan, 0

        rmse = np.sqrt(np.mean((m - o) ** 2))
        denom = np.sum(o)
        rel_bias = np.nan if np.isclose(denom, 0.0) else 100.0 * np.sum(m - o) / denom
        return mask, rmse, rel_bias, len(o)

    mask_ds, rmse_ds, rb_ds, n_ds = _metrics(obs_flat, ds_flat)
    mask_hs2, rmse_hs2, rb_hs2, n_hs2 = _metrics(obs_flat, hs2_flat)

    all_vals = np.concatenate([
        obs_flat[mask_ds], ds_flat[mask_ds],
        obs_flat[mask_hs2], hs2_flat[mask_hs2],
    ])
    vmin, vmax = np.nanmin(all_vals), np.nanmax(all_vals)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharex=True, sharey=True, constrained_layout=True)

    panels = [
        ("delta-snow", ds_flat, mask_ds, rmse_ds, rb_ds, n_ds),
        ("hs2swe", hs2_flat, mask_hs2, rmse_hs2, rb_hs2, n_hs2),
    ]

    for ax, (name, mod_flat, mask, rmse, rb, n) in zip(axes, panels):
        ax.scatter(obs_flat[mask], mod_flat[mask], s=30, alpha=0.75)
        ax.plot([vmin, vmax], [vmin, vmax], "k--", lw=1)
        ax.set_title(f"{station}: {name}")
        ax.set_xlabel("Observed SWE [mm]")
        ax.set_ylabel("Modelled SWE [mm]")
        ax.grid(True, alpha=0.3)
        ax.set_xlim(vmin, vmax)
        ax.set_ylim(vmin, vmax)
        ax.set_aspect("equal", adjustable="box")
        ax.text(
            0.03,
            0.97,
            f"RMSE = {rmse:.2f} mm\nRel. bias = {rb:+.2f}%\nn = {n}",
            transform=ax.transAxes,
            va="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
        )

    fig.savefig(out_file, dpi=200)
    plt.close(fig)

    return {
        "rmse_ds": rmse_ds,
        "rel_bias_ds": rb_ds,
        "n_ds": n_ds,
        "rmse_hs2": rmse_hs2,
        "rel_bias_hs2": rb_hs2,
        "n_hs2": n_hs2,
    }


def save_first_appearance_summary_plot(
    station: str,
    ds: xr.Dataset,
    hs2: xr.Dataset,
    out_file: Path,
    center_stat: str = "mean",
    error_stat: str = "std",
    errorevery: int = 7,
) -> pd.DataFrame:
    hs_ds_layers = ds["HS"]
    hs_hs2_layers = hs2["HS_layer_cm"]

    ev_ds = first_appearance_events(hs_ds_layers, scale_to_m=1.0)
    ev_hs2 = first_appearance_events(hs_hs2_layers, scale_to_m=0.01)

    season_vals = np.intersect1d(ds["season"].values, hs2["season"].values).astype(int)
    dos_vals = np.intersect1d(ds["dos"].values, hs2["dos"].values).astype(int)

    cum_ds = cumulative_matrix_from_events(ev_ds, season_vals, dos_vals)
    cum_hs2 = cumulative_matrix_from_events(ev_hs2, season_vals, dos_vals)

    stats_ds = curve_stats(cum_ds)
    stats_hs2 = curve_stats(cum_hs2)

    summary = pd.DataFrame(
        {
            "delta-snow": final_summary(cum_ds),
            "hs2swe": final_summary(cum_hs2),
        }
    ).T

    fig, ax = plt.subplots(figsize=(12, 5))

    for stats, color, label in [
        (stats_ds, "C0", "delta-snow"),
        (stats_hs2, "C1", "hs2swe"),
    ]:
        center = stats[center_stat]

        if error_stat == "std":
            yerr = stats["std"].to_numpy()
        else:
            yerr = np.vstack([
                (center - stats["q25"]).to_numpy(),
                (stats["q75"] - center).to_numpy(),
            ])

        ax.plot(dos_vals, center.to_numpy(), color=color, lw=2, label=f"{label} {center_stat}")
        ax.errorbar(
            dos_vals,
            center.to_numpy(),
            yerr=yerr,
            fmt="none",
            ecolor=color,
            alpha=0.35,
            elinewidth=1,
            capsize=2,
            errorevery=errorevery,
        )

    ax.set_title(f"{station}: Cumulative first-appearance HS across seasons")
    ax.set_xlabel("DOS")
    ax.set_ylabel("Cumulative HS [m]")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.savefig(out_file, dpi=200)
    plt.close(fig)

    return summary


def save_first_appearance_by_season_plots(
    station: str,
    ds: xr.Dataset,
    hs2: xr.Dataset,
    out_dir: Path,
) -> None:
    hs_ds_layers = ds["HS"]
    hs_hs2_layers = hs2["HS_layer_cm"]

    ev_ds = first_appearance_events(hs_ds_layers, scale_to_m=1.0)
    ev_hs2 = first_appearance_events(hs_hs2_layers, scale_to_m=0.01)

    season_vals = np.intersect1d(ds["season"].values, hs2["season"].values).astype(int)

    station_dir = out_dir / station
    station_dir.mkdir(parents=True, exist_ok=True)

    for s in season_vals:
        ds_df = season_curve(ev_ds, s)
        hs2_df = season_curve(ev_hs2, s)

        fig, ax = plt.subplots(figsize=(12, 4))

        if not ds_df.empty:
            ax.plot(ds_df.index, ds_df["cumulative"], lw=2, label="delta-snow", color="C0")
            ax.scatter(ds_df.index, ds_df["cumulative"], s=18, color="C0", zorder=3)

        if not hs2_df.empty:
            ax.plot(hs2_df.index, hs2_df["cumulative"], lw=2, label="hs2swe", color="C1")
            ax.scatter(hs2_df.index, hs2_df["cumulative"], s=18, color="C1", zorder=3)

        ax.set_title(f"{station}: First-appearance cumulative HS — season {s-1}/{s}")
        ax.set_xlabel("Time")
        ax.set_ylabel("Cumulative HS [m]")
        ax.grid(True, alpha=0.3)
        ax.legend()

        fig.savefig(station_dir / f"{station}_first_appearance_{s-1}_{s}.png", dpi=200)
        plt.close(fig)


def save_first_appearance_all_stations_plot(
    file_table: pd.DataFrame,
    out_file: Path,
    center_stat: str = "mean",
    error_stat: str = "std",
    errorevery: int = 7,
) -> pd.DataFrame:
    ds_cum_list: list[pd.DataFrame] = []
    hs2_cum_list: list[pd.DataFrame] = []

    for row in file_table.itertuples(index=False):
        ds = xr.open_dataset(row.dsnow_file)
        hs2 = xr.open_dataset(row.hs2swe_file)

        try:
            hs_ds_layers = ds["HS"]
            hs_hs2_layers = hs2["HS_layer_cm"]

            ev_ds = first_appearance_events(hs_ds_layers, scale_to_m=1.0)
            ev_hs2 = first_appearance_events(hs_hs2_layers, scale_to_m=0.01)

            season_vals = np.intersect1d(ds["season"].values, hs2["season"].values).astype(int)
            dos_vals = np.intersect1d(ds["dos"].values, hs2["dos"].values).astype(int)

            if len(season_vals) == 0 or len(dos_vals) == 0:
                continue

            cum_ds = cumulative_matrix_from_events(ev_ds, season_vals, dos_vals)
            cum_hs2 = cumulative_matrix_from_events(ev_hs2, season_vals, dos_vals)

            cum_ds.columns = [f"{row.station}_{int(s)}" for s in cum_ds.columns]
            cum_hs2.columns = [f"{row.station}_{int(s)}" for s in cum_hs2.columns]

            ds_cum_list.append(cum_ds)
            hs2_cum_list.append(cum_hs2)

        finally:
            ds.close()
            hs2.close()

    if not ds_cum_list or not hs2_cum_list:
        raise RuntimeError("No valid first-appearance data available for all-station plot.")

    all_cum_ds = pd.concat(ds_cum_list, axis=1, join="outer").sort_index()
    all_cum_hs2 = pd.concat(hs2_cum_list, axis=1, join="outer").sort_index()

    stats_ds = curve_stats(all_cum_ds)
    stats_hs2 = curve_stats(all_cum_hs2)

    summary = pd.DataFrame(
        {
            "delta-snow": final_summary(all_cum_ds.fillna(method="ffill").fillna(0.0)),
            "hs2swe": final_summary(all_cum_hs2.fillna(method="ffill").fillna(0.0)),
        }
    ).T

    dos_vals = stats_ds.index.to_numpy(dtype=int)

    fig, ax = plt.subplots(figsize=(12, 5))

    for stats, color, label in [
        (stats_ds, "C0", "delta-snow"),
        (stats_hs2, "C1", "hs2swe"),
    ]:
        center = stats[center_stat]

        if error_stat == "std":
            yerr = stats["std"].to_numpy()
        else:
            yerr = np.vstack([
                (center - stats["q25"]).to_numpy(),
                (stats["q75"] - center).to_numpy(),
            ])

        ax.plot(dos_vals, center.to_numpy(), color=color, lw=2, label=f"{label} {center_stat}")
        ax.errorbar(
            dos_vals,
            center.to_numpy(),
            yerr=yerr,
            fmt="none",
            ecolor=color,
            alpha=0.35,
            elinewidth=1,
            capsize=2,
            errorevery=errorevery,
        )

    ax.set_title("All stations: Cumulative first-appearance HS across all seasons")
    ax.set_xlabel("DOS")
    ax.set_ylabel("Cumulative HS [m]")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.savefig(out_file, dpi=200)
    plt.close(fig)

    return summary


# ============================================================================
# PER-STATION PROCESSING
# ============================================================================

def process_station(
    station: str,
    dsnow_file: Path,
    hs2swe_file: Path,
    obs_file: Path,
) -> dict[str, float | str]:
    print(f"Processing {station} ...")

    ds = xr.open_dataset(dsnow_file)
    hs2 = xr.open_dataset(hs2swe_file)

    try:
        target_seasons = hs2["season"].values
        target_dos = hs2["dos"].values

        if not obs_file.exists():
            raise FileNotFoundError(f"Observation file not found: {obs_file}")

        obs_daily = load_and_prepare_observations(obs_file, target_seasons, target_dos)
        hs_obs_m_hs2, swe_obs_mm_hs2 = observations_to_xarray(obs_daily, target_seasons, target_dos)

        ds, hs2 = add_observations_to_datasets(ds, hs2, hs_obs_m_hs2, swe_obs_mm_hs2)

        ds_out = NC_OUT_DIR / f"{dsnow_file.stem}_with_obs.nc"
        hs2_out = NC_OUT_DIR / f"{hs2swe_file.stem}_with_obs.nc"
        ds.to_netcdf(ds_out)
        hs2.to_netcdf(hs2_out)

        series = prepare_model_series(ds, hs2)

        swe_obs = ds["SWE_obs_mm"].transpose("season", "dos")
        swe_obs, swe_ds, swe_hs2 = xr.align(
            swe_obs, series["swe_ds"], series["swe_hs2"], join="inner"
        )

        save_swe_mean_plot_with_obs(
            station=station,
            swe_ds=swe_ds,
            swe_hs2=swe_hs2,
            swe_obs=swe_obs,
            out_file=SWE_MEAN_DIR / f"{station}_swe_mean_with_obs.png",
        )

        save_swe_by_season_plot(
            station=station,
            swe_ds=swe_ds,
            swe_hs2=swe_hs2,
            swe_obs=swe_obs,
            out_file=SWE_SEASON_DIR / f"{station}_swe_by_season.png",
        )

        scatter_metrics = save_swe_scatter_plot(
            station=station,
            swe_obs=swe_obs,
            swe_ds=swe_ds,
            swe_hs2=swe_hs2,
            out_file=SWE_SCATTER_DIR / f"{station}_swe_scatter.png",
        )

        first_app_summary = save_first_appearance_summary_plot(
            station=station,
            ds=ds,
            hs2=hs2,
            out_file=FIRST_APP_SUMMARY_DIR / f"{station}_first_appearance_summary.png",
        )

        save_first_appearance_by_season_plots(
            station=station,
            ds=ds,
            hs2=hs2,
            out_dir=FIRST_APP_SEASON_DIR,
        )

        row = {
            "station": station,
            "dsnow_file": str(dsnow_file),
            "hs2swe_file": str(hs2swe_file),
            "obs_file": str(obs_file),
            "ds_out": str(ds_out),
            "hs2_out": str(hs2_out),
            **scatter_metrics,
            "fa_mean_ds": float(first_app_summary.loc["delta-snow", "mean"]),
            "fa_mean_hs2": float(first_app_summary.loc["hs2swe", "mean"]),
            "fa_median_ds": float(first_app_summary.loc["delta-snow", "median"]),
            "fa_median_hs2": float(first_app_summary.loc["hs2swe", "median"]),
            "fa_std_ds": float(first_app_summary.loc["delta-snow", "std"]),
            "fa_std_hs2": float(first_app_summary.loc["hs2swe", "std"]),
            "status": "ok",
        }

        print(f"Finished {station}")
        return row

    except Exception as e:
        print(f"Failed for {station}: {e}")
        return {
            "station": station,
            "dsnow_file": str(dsnow_file),
            "hs2swe_file": str(hs2swe_file),
            "obs_file": str(obs_file),
            "status": f"failed: {e}",
        }

    finally:
        ds.close()
        hs2.close()


# ============================================================================
# MAIN
# ============================================================================

def main() -> None:
    file_table = build_station_file_table(DSNOW_DIR, HS2SWE_DIR, OBS_DIR)

    if file_table.empty:
        raise RuntimeError("No matching station files found.")

    print("Found stations:")
    print(file_table[["station", "obs_exists"]].to_string(index=False))

    results = []
    for row in file_table.itertuples(index=False):
        if not row.obs_exists:
            print(f"Skipping {row.station}: missing observation file")
            results.append(
                {
                    "station": row.station,
                    "dsnow_file": str(row.dsnow_file),
                    "hs2swe_file": str(row.hs2swe_file),
                    "obs_file": str(row.obs_file),
                    "status": "missing observation file",
                }
            )
            continue

        results.append(
            process_station(
                station=row.station,
                dsnow_file=row.dsnow_file,
                hs2swe_file=row.hs2swe_file,
                obs_file=row.obs_file,
            )
        )

    summary_df = pd.DataFrame(results)
    summary_csv = OUTPUT_DIR / "station_comparison_summary.csv"
    summary_df.to_csv(summary_csv, index=False)

    all_station_first_app_summary = save_first_appearance_all_stations_plot(
        file_table=file_table,
        out_file=FIRST_APP_ALL_STATIONS_DIR / "all_stations_first_appearance_summary.png",
    )
    all_station_first_app_summary.to_csv(
        OUTPUT_DIR / "all_stations_first_appearance_summary_stats.csv"
    )

    print("\nBatch processing finished.")
    print(f"Summary written to: {summary_csv}")


if __name__ == "__main__":
    main()