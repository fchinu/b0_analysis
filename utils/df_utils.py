"""Module containing utility functions for working with dataframes."""
import pandas as pd
import pyarrow.parquet as pq

def read_parquet_in_batches(file_path, selections=None, batch_size=1000000):
    """
    Read a Parquet file in batches and return a concatenated DataFrame.

    Parameters:
    file_path (str): The path to the Parquet file.
    selections (str, optional): A string representing the selection criteria to apply to each batch. Defaults to None.
    batch_size (int, optional): The number of rows to read per batch. Defaults to 1000000.

    Returns:
    pandas.DataFrame: The concatenated DataFrame.

    """
    parquet_file = pq.ParquetFile(file_path)
    df = []
    for batch in parquet_file.iter_batches(batch_size):
        batch_df = batch.to_pandas()
        if selections is not None:
            batch_df = batch_df.query(selections)
        df.append(batch_df)
    return pd.concat(df)
