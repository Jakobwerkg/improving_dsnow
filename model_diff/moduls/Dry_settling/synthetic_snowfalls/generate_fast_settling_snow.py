import pandas as pd
import numpy as np
from pathlib import Path

# settings
start_date = "2024-10-22"
end_date = "2025-04-20"

snowfall_events = {
    "2024-10-27": 0.6,
    "2024-11-21": 0.7,
    "2024-12-25": 0.55,
}

output_file = Path("/Users/jakobwerkgarner/code/mt_dsnow/model_diff/moduls/Dry_settling/synthetic_snowfalls/synthetic_big_snowfalls_100days.csv")
output_file.parent.mkdir(parents=True, exist_ok=True)

dates = pd.date_range(start=start_date, end=end_date, freq="D")

hs = 0.0
data = []

for i, d in enumerate(dates):
    ds = d.strftime("%Y-%m-%d")

    # snowfall event
    if ds in snowfall_events:
        hs += snowfall_events[ds]

    elif hs > 0:
        # --- nonlinear settling ---
        if hs > 0.5:
            decay = 0.75   # fast settling for fresh snow
        elif hs > 0.2:
            decay = 0.88   # medium
        else:
            decay = 0.95   # slow settling for old snow

        hs *= decay

        # --- small variability (plateaus + noise) ---
        noise = np.random.normal(0, 0.01)
        hs += noise

        # prevent negative
        hs = max(hs, 0)

        # create plateaus (like your example)
        if np.random.rand() < 0.2:
            hs = hs  # hold value (no change)

    # rounding → makes it look like real measurements
    hs = round(hs, 2)

    data.append([ds, hs])

df = pd.DataFrame(data, columns=["timestamp", "HS_meas"])
df.to_csv(output_file, index=False)

if hs < 0.05:
    hs = 0.05  # set to zero if very small to mimic measurement threshold
print(f"Saved to: {output_file}")
print(df.head(30))