import pandas as pd
import altair as alt
import numpy as np

def plot_umap(df):
    # Create the Altair scatter plot
    fig = alt.Chart(df).mark_circle(size=60).properties(
    width = 700,
    height = 700,
    ).encode(
        x=alt.X('UMAP1', title='UMAP 1'),
        y=alt.Y('UMAP2', title='UMAP 2'),
        color=alt.Color('file_name:N', scale=alt.Scale(scheme="category10"), legend=None),
        tooltip=['image', 'file_name', 'molecule_index']
    ).properties(
        title='2D UMAP Projection of Molecules'
    ).configure_title(
        anchor='middle'
    ).interactive()
    return fig

def plot_pca(df):
    # Create the Altair scatter plot
    fig = alt.Chart(df).mark_circle(size=60).properties(
    width = 700,
    height = 700,
    ).encode(
        x=alt.X('PCA1', title='PCA 1'),
        y=alt.Y('PCA2', title='PCA 2'),
        color=alt.Color('file_name:N', scale=alt.Scale(scheme="category10"), legend=None),
        tooltip=['image', 'file_name', 'molecule_index']
    ).properties(
        title='2D PCA Projection of Molecules'
    ).configure_title(
        anchor='middle'
    ).interactive()
    return fig