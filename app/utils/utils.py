import pandas as pd

def estimate_table_height(df):
    """
    Estimates an appropriate height for displaying a DataFrame in Streamlit,
    capping the height to avoid excessive scrolling.
    """

    if df is None or df.empty:
        return 128  # very little in case it's empty
    
    row_height = 38 # Default row height
    base_padding = 200  # Base padding

    height = len(df) * row_height + base_padding
    height = max(height, 128) # At least 128px always
    height = min (height, 1024) # At most 1024px

    return height