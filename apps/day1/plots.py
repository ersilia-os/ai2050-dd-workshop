import pandas as pd
import altair as alt
import numpy as np

def plot_umap(df_list):
    df = df_list ### WIP
    # Create the Altair scatter plot
    fig = alt.Chart(df).mark_circle(size=60).encode(
        x=alt.X('UMAP1', title='UMAP 1'),
        y=alt.Y('UMAP2', title='UMAP 2'),
        #color=alt.Color('Binary:N', scale=alt.Scale(domain=[0, 1], range=['#0000FF', '#FF0000']), legend=None),
        tooltip=['image'] #'UMAP1', 'UMAP2', 'Binary']
    ).properties(
        title='2D UMAP Projection of Molecules'
    ).configure_title(
        anchor='middle'
    ).interactive()
    return fig