def load_SLF_data(path):


    import pandas as pd 
    import numpy as np
    """
    Load SLF snow depth observation data from a text file.

    The function searches for the table header (line starting with 'Time,'),
    then reads the data into a pandas DataFrame, skipping any metadata or comments
    before the header.

    Parameters
    ----------
    path : str
        Path to the SLF observation text file.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the loaded data, with columns as specified in the file.

    Example
    -------
    >>> df = load_SLF_data('HN/OBS-HN.txt')
    >>> print(df.head())
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    # Find the line index where the table header starts
    for i, line in enumerate(lines):
        if line.startswith('Time,'):
            header_idx = i
            break

    data = pd.read_csv(path, skiprows=header_idx)
    data = data.replace(-9999, np.nan)
    # data.index = pd.to_datetime(data['Time'], format='%Y-%m-%d %H:%M:%S')
    # data = data.drop(columns=['Time'])
    return data
