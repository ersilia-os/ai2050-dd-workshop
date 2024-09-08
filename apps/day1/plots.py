import pandas as pd
import altair as alt
import numpy as np
import streamlit as st

def plot_umap(df):
    # Create the Altair scatter plot
    fig = alt.Chart(df).mark_circle(size=20).properties(
    width = 500,
    height = 500,
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
    fig = alt.Chart(df).mark_circle(size=20).properties(
    width = 500,
    height = 500,
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

def plot_mol_weights(df):
    fig = alt.Chart(df).transform_density('mol weight', as_=['MOL WEIGHT', 'DENSITY'], groupby=['file_name']).mark_area(opacity=0, line=True).encode(
        x=alt.X('MOL WEIGHT:Q', title="Mol. Weight"),
        y=alt.Y('DENSITY:Q', title="Density"),
        color=alt.Color('file_name:N', scale=alt.Scale(scheme="category10"), legend=None)
    ).interactive()
    return fig

def plot_logp(df):
    fig = alt.Chart(df).transform_density('logp', as_=['LOGP', 'DENSITY'], groupby=['file_name']).mark_area(opacity=0, line=True).encode(
        x=alt.X('LOGP:Q', title="LogP"),
        y=alt.Y('DENSITY:Q', title="Density"),
        color=alt.Color('file_name:N', scale=alt.Scale(scheme="category10"), legend=None)
    ).interactive()
    return fig

def plot_qed(df):
    fig = alt.Chart(df).transform_density('qed', as_=['QED', 'DENSITY'], groupby=['file_name']).mark_area(opacity=0, line=True).encode(
        x=alt.X('QED:Q', title="QED"),
        y=alt.Y('DENSITY:Q', title="Density"),
        color=alt.Color('file_name:N', scale=alt.Scale(scheme="category10"), legend=None)
    ).interactive()
    return fig

def plot_legend(df):
    fig = alt.Chart(df).mark_point(size=0).properties(
    width = 200,
    height = 100,
    ).encode(
        color=alt.Color('file_name:N', scale=alt.Scale(scheme="category10"))
    )
    return fig
    
