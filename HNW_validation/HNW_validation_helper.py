import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os


def validate_hnw(df,
                 obs_col="HNW_obs",
                 mod_col="HNW_mod",
                 save_dir=None,
                 filename="hnw_validation.png"):
    """
    Validate modeled vs observed HNW and create a density scatter plot.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing observed and modeled HNW values.

    obs_col : str, default="HNW_obs"
        Column name of observed HNW values.

    mod_col : str, default="HNW_mod"
        Column name of modeled HNW values.

    save_dir : str or None, optional
        Directory where the plot should be saved. If None, the plot is only shown.

    filename : str, default="hnw_validation.png"
        Name of the saved figure file.

    Returns
    -------
    dict
        Dictionary containing validation metrics:
        {
            "N": number of samples,
            "RMSE": root mean square error,
            "Bias": mean bias,
            "PBIAS": relative bias,
            "R2": coefficient of determination
        }

    Notes
    -----
    - Rows with NaN or non-finite values are removed.
    - Observations <= 0 are filtered out.
    - The plot shows a 2D density histogram with a 1:1 reference line.
    """

    # ---------------------------------------------------------
    # DATA FILTERING
    # ---------------------------------------------------------
    df_valid = df.dropna(subset=[obs_col, mod_col])
    df_valid = df_valid[df_valid[obs_col] >= 0]

    df_valid = df_valid[
        np.isfinite(df_valid[obs_col]) &
        np.isfinite(df_valid[mod_col])
    ]

    y = df_valid[obs_col].values
    x = df_valid[mod_col].values

    # ---------------------------------------------------------
    # METRICS
    # ---------------------------------------------------------
    residuals = x - y

    rmse = np.sqrt(np.mean(residuals**2))
    bias = np.mean(residuals)
    pbias = np.sum(residuals) / np.sum(y)

    ss_res = np.sum((y - x)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    r2 = 1 - ss_res / ss_tot

    stats = {
        "N": len(df_valid),
        "RMSE": rmse,
        "Bias": bias,
        "PBIAS": pbias,
        "R2": r2
    }

    print(f"N        = {stats['N']}")
    print(f"RMSE     = {stats['RMSE']:.3f}")
    print(f"Bias     = {stats['Bias']:.3f}")
    print(f"PBIAS    = {stats['PBIAS']:.2f}")
    print(f"R²       = {stats['R2']:.3f}")

    # ---------------------------------------------------------
    # PLOT
    # ---------------------------------------------------------
    plt.figure(figsize=(8, 7))

    plt.hist2d(
        x, y,
        bins=50,
        range=[[0, 100], [0, 100]],
        norm=LogNorm(vmin=1, vmax=1000),
        cmap="jet"
    )

    cb = plt.colorbar(label="Number of observations")
    cb.set_ticks([1, 10, 100, 1000])
    cb.set_ticklabels(["1", "10", "100", ">999"])

    lim = [0, 100]
    plt.plot(lim, lim, "--", color="gray", linewidth=1.3)

    ticks = [0, 25, 50, 75, 100]
    plt.xticks(ticks)
    plt.yticks(ticks)

    plt.xlabel("Modeled HNW (mm)")
    plt.ylabel("Observed HNW (mm)")
    plt.title("ΔSNOW", fontsize=15)

    textstr = (
        f"$R^2$: {r2:.2f}\n"
        f"Rel. bias: {pbias:.2f}\n"
        f"RMSE: {rmse:.1f}"
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

    # ---------------------------------------------------------
    # SAVE (OPTIONAL)
    # ---------------------------------------------------------
    if save_dir is not None:
        os.makedirs(save_dir, exist_ok=True)
        save_path = os.path.join(save_dir, filename)
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Plot saved to: {save_path}")

    plt.show()

    return stats