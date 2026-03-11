import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.cm import ScalarMappable



def plot_layer_stack(
    ds,
    var,
    hs_var="HS",
    x="dos",
    norm=False,
    cmap="viridis",
    ax=None,
    centered_cmap=False,
):
    """
    Plot snowpack layers as stacked bars colored by a variable.

    Parameters
    ----------
    ds : xarray.Dataset
    var : str
        Variable used for shading (e.g. "RHO", "AGE", "SWE_layers")
    hs_var : str
        Layer thickness variable
    x : str
        X coordinate
    norm : bool
        If True normalize shading variable to [0,1]
    cmap : str
        Colormap
    centered_cmap : bool
        Center the colormap around zero (useful for anomalies)
    """

    hs = ds[hs_var].transpose("layer", x).values
    shade = ds[var].transpose("layer", x).values
    xvals = ds[x].values

    n_layers, n_times = hs.shape

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))
    else:
        fig = ax.figure

    width = 0.9 * np.min(np.diff(xvals)) if len(xvals) > 1 else 0.8

    # ---- compute limits ----
    valid = np.isfinite(shade)
    vmin = np.nanmin(shade[valid]) if np.any(valid) else 0
    vmax = np.nanmax(shade[valid]) if np.any(valid) else 1

    # ---- normalization ----
    if centered_cmap:
        vmax_abs = np.nanmax(np.abs(shade[valid]))
        color_norm = TwoSlopeNorm(vmin=-vmax_abs, vcenter=0, vmax=vmax_abs)
        shade_plot = shade

    elif norm:
        shade_plot = (shade - vmin) / (vmax - vmin + 1e-20)
        color_norm = Normalize(0, 1)

    else:
        shade_plot = shade
        color_norm = Normalize(vmin=vmin, vmax=vmax)

    cmap_obj = plt.get_cmap(cmap)

    bottom = np.zeros(n_times)

    for i in range(n_layers):
        h = hs[i]
        c = shade_plot[i]

        colors = [
            cmap_obj(color_norm(val)) if np.isfinite(val) and hh > 0 else (0, 0, 0, 0)
            for val, hh in zip(c, h)
        ]

        ax.bar(
            xvals,
            h,
            bottom=bottom,
            width=width,
            color=colors,
            align="center",
            edgecolor="none",
        )

        bottom += np.nan_to_num(h)

    # ---- colorbar ----
    sm = ScalarMappable(norm=color_norm, cmap=cmap_obj)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)

    if norm:
        cbar.set_label(f"{var} (normalized)")
    else:
        cbar.set_label(var)

    ax.set_xlabel(x)
    ax.set_ylabel(hs_var)
    ax.set_title(f"Stacked snow layers colored by {var}")
    ax.set_xlim(xvals.min() - width, xvals.max() + width)

    return fig, ax




def plot_var_seasons(
    var_name,
    mode="layer_sum",
    ds=None,
    spread="std",
    show_layer_std=True,
    running_mean=1,
):
    """
    Plot variables with simple layer/season aggregation options.

    Parameters
    ----------
    var_name : str
        Variable name in dataset.
    mode : str
        "layer_sum"              -> sum over layers, plot all seasons
        "layer_mean"             -> mean over layers, plot all seasons
        "season_mean_layer_sum"  -> sum over layers, then mean over seasons
        "season_mean_layer_mean" -> mean over layers, then mean over seasons
    ds : xarray.Dataset
        Input dataset.
    spread : str
        For seasonal mean modes: "std" or "iqr"
    show_layer_std : bool
        Only used for mode="layer_mean".
        If True, shade ±1 std across layers for each season.
    running_mean : int
        Running mean window along dos for seasonal mean modes.
        Use 1 for no smoothing, 3 for 3-day smoothing, 5 for 5-day smoothing, etc.
    """

    if ds is None:
        raise ValueError("Provide dataset with ds=...")

    if var_name not in ds:
        raise ValueError(f"{var_name!r} not found in dataset.")

    if running_mean < 1:
        raise ValueError("running_mean must be >= 1")

    da = ds[var_name]

    if "dos" not in da.dims or "season" not in da.dims:
        raise ValueError(f"{var_name!r} must have 'dos' and 'season' dimensions.")

    if "layer" not in da.dims:
        raise ValueError(f"{var_name!r} must have a 'layer' dimension for these modes.")

    dos = da["dos"].values
    seasons = da["season"].values

    fig, ax = plt.subplots(figsize=(11, 5))


    if mode not in ["values", "layer_sum", "layer_mean", "season_mean_layer_sum", "season_mean_layer_mean"]:
        raise ValueError(
            "mode must be one of: "
            "'values', 'layer_sum', 'layer_mean', 'season_mean_layer_sum', 'season_mean_layer_mean'"
        )
    


    # 1) Sum over layers, plot every season
    if mode == "layer_sum":
        da_plot = da.sum(dim="layer", skipna=True).transpose("dos", "season")
        colors = plt.cm.tab20(np.linspace(0, 1, len(seasons)))

        for i, s in enumerate(seasons):
            y = da_plot.sel(season=s).values
            ax.plot(dos, y, lw=1.5, alpha=0.9, color=colors[i], label=str(s))

        title = f"{var_name}: sum over layers"

        
    # 1.0 Plot the actual values 
    elif mode == "values":
        da_plot = da.transpose("dos", "season")
        colors = plt.cm.tab20(np.linspace(0, 1, len(seasons)))

        for i, s in enumerate(seasons):
            y = da_plot.sel(season=s).values
            ax.plot(dos, y, lw=1.5, alpha=0.9, color=colors[i], label=str(s))

        title = f"{var_name}: values by layer and season"

    # 2) Mean over layers, plot every season, optionally with layer std
    elif mode == "layer_mean":
        da_mean = da.mean(dim="layer", skipna=True).transpose("dos", "season")
        da_std = da.std(dim="layer", skipna=True).transpose("dos", "season")
        colors = plt.cm.tab20(np.linspace(0, 1, len(seasons)))

        for i, s in enumerate(seasons):
            y = da_mean.sel(season=s).values
            ax.plot(dos, y, lw=1.5, alpha=0.9, color=colors[i], label=str(s))

            if show_layer_std:
                ystd = da_std.sel(season=s).values
                ax.fill_between(dos, y - ystd, y + ystd, color=colors[i], alpha=0.15)

        title = f"{var_name}: mean over layers"
        if show_layer_std:
            title += " with layer std"

    # 3) Sum over layers, then mean over seasons, then smooth
    elif mode == "season_mean_layer_sum":
        da_tmp = da.sum(dim="layer", skipna=True).transpose("dos", "season")

        mean = da_tmp.mean(dim="season", skipna=True)

        if spread == "std":
            spread_da = da_tmp.std(dim="season", skipna=True)
            low = mean - spread_da
            high = mean + spread_da
            spread_label = "±1 std"
        elif spread == "iqr":
            low = da_tmp.quantile(0.25, dim="season")
            high = da_tmp.quantile(0.75, dim="season")
            spread_label = "IQR"
        else:
            raise ValueError("spread must be 'std' or 'iqr'.")

        if running_mean > 1:
            mean = mean.rolling(dos=running_mean, center=True, min_periods=1).mean()
            low = low.rolling(dos=running_mean, center=True, min_periods=1).mean()
            high = high.rolling(dos=running_mean, center=True, min_periods=1).mean()

        ax.plot(dos, mean.values, lw=2, label="mean over seasons")
        ax.fill_between(dos, low.values, high.values, alpha=0.3, label=spread_label)

        title = f"{var_name}: seasonal mean of layer sum"
        if running_mean > 1:
            title += f" ({running_mean}-day running mean)"

    # 4) Mean over layers, then mean over seasons, then smooth
    elif mode == "season_mean_layer_mean":
        da_tmp = da.mean(dim="layer", skipna=True).transpose("dos", "season")

        mean = da_tmp.mean(dim="season", skipna=True)

        if spread == "std":
            spread_da = da_tmp.std(dim="season", skipna=True)
            low = mean - spread_da
            high = mean + spread_da
            spread_label = "±1 std"
        elif spread == "iqr":
            low = da_tmp.quantile(0.25, dim="season")
            high = da_tmp.quantile(0.75, dim="season")
            spread_label = "IQR"
        else:
            raise ValueError("spread must be 'std' or 'iqr'.")

        if running_mean > 1:
            mean = mean.rolling(dos=running_mean, center=True, min_periods=1).mean()
            low = low.rolling(dos=running_mean, center=True, min_periods=1).mean()
            high = high.rolling(dos=running_mean, center=True, min_periods=1).mean()

        ax.plot(dos, mean.values, lw=2, label="mean over seasons")
        ax.fill_between(dos, low.values, high.values, alpha=0.3, label=spread_label)

        title = f"{var_name}: seasonal mean of layer mean"
        if running_mean > 1:
            title += f" ({running_mean}-day running mean)"

    else:
        raise ValueError(
            "mode must be one of: "
            "'layer_sum', 'layer_mean', 'season_mean_layer_sum', 'season_mean_layer_mean'"
        )

    ax.set_title(title)
    ax.set_xlabel("day of season")
    ax.set_ylabel(var_name)
    ax.grid(alpha=0.3)

    if mode in ["layer_sum", "layer_mean"]:
        ax.legend(title="season", bbox_to_anchor=(1.02, 1), loc="upper left")
    else:
        ax.legend()

    plt.tight_layout()
    plt.show()