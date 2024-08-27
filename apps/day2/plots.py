import pandas as pd
import altair as alt
from sklearn.metrics import confusion_matrix

def plot_act_inact(df): 
    #needs a datagrame with Molecule index, Mean value and Binary
    color_mapping = alt.Color('Binary:N', scale=alt.Scale(domain=[1, 0], range=['#FF0000', '#0000FF']), legend=alt.Legend(title='Molecule Activity', labelExpr="if(datum.label == '1', 'Active', 'Inactive')", orient="bottom"))
    scatter_plot = alt.Chart(df).mark_circle(size=60).encode(
        x=alt.X('Molecule index:Q', title='Molecule Index'),
        y=alt.Y('Mean:Q', title='Mean growth (OD)'),
        color=color_mapping
    ).properties(
        title='Molecule activity for A.baumannii'
    ).configure_title(
            anchor='middle'
        ).interactive()
    return scatter_plot

def plot_roc_curve(df):
    fig = alt.Chart(df).transform_fold(
        ['tpr_cv1', 'tpr_cv2', 'tpr_cv3', 'tpr_cv4', 'tpr_cv5', 'Mean TPR'],
        as_=['Variable', 'Value']
    ).mark_line().encode(
        x=alt.X('FPR:Q', title='False Positive Rate (FPR)'),
        y=alt.Y('Value:Q', title='True Positive Rate (TPR)'),
        color=alt.Color('Variable:N', scale=alt.Scale(range= ['#0000FF']+['#d3d3d3']*5), legend=None),
    ).properties(
        title='ROC Curve'
    ).configure_title(
        anchor='middle'
    ).interactive()
    return fig

def plot_lolp(X,y):
    X_ = [x[:2] for x in X]
    LolP1 = [arr[0] for arr in X_]
    LolP2 = [arr[1] for arr in X_]
    lolp_df = pd.DataFrame({
        'LolP1': LolP1,
        'LolP2': LolP2,
        'Binary': y,
        'Color': ['#0000FF' if x==0 else  '#FF0000'for x in y]
    })  
    lolp_df_sorted = lolp_df.sort_values(by='Binary', ascending=True)
    fig = alt.Chart(lolp_df_sorted).mark_circle(size=60).encode(
    x=alt.X('LolP1', title='LolP1'),
    y=alt.Y('LolP2', title='LolP2'),
    color=alt.Color('Binary:N', scale=alt.Scale(domain=[0, 1], range=['#0000FF', '#FF0000']), legend=None)
    ).properties(
        title='2D Chemical Space'
    ).configure_title(
        anchor='middle'
    ).interactive()
    return fig

def _create_confusion_matrix_df(y_test, y_pred, cutoff):
    y_pred_bin = (y_pred >= cutoff).astype(int)
    cm = confusion_matrix(y_test, y_pred_bin)
    classes = sorted(set(y_test))
    data = []
    for i, actual_class in enumerate(classes):
        for j, predicted_class in enumerate(classes):
            data.append({
                'Real': str(actual_class),
                'Predicted': str(predicted_class),
                'confusion_matrix': cm[i, j]
            })
    
    df_cm = pd.DataFrame(data)
    return df_cm

def plot_contingency_table(y_test, y_pred, cutoff=0.5):
    df = _create_confusion_matrix_df(y_test, y_pred, cutoff)
    color = alt.Color('confusion_matrix:Q', scale=alt.Scale(scheme="blues"), legend=None)
    heatmap = (
        alt.Chart(df).mark_rect()
        .encode(
            alt.X("Predicted:N"),
            alt.Y("Real:N"),
            color=color,
        )
    )

    text = (
        alt.Chart(df)
        .mark_text(baseline="middle", fontSize=16)
        .encode(
            x="Predicted:N",
            y="Real:N",
            text="confusion_matrix:N",
            color=alt.condition(
                alt.datum.confusion_matrix > 0,
                alt.value("black"),
                alt.value("black"),
            ),
        )
    )
    cm_chart = (heatmap + text).properties(width=300, height=400)
    return cm_chart


def plot_cm_chart(df, alt):
    color = alt.Color(scale=alt.Scale(scheme="oranges"))
    heatmap = (
        alt.Chart(df, title="Confusion Matrix")
        .mark_rect()
        .encode(
            alt.X("predicted:N"),
            alt.Y("actual:N"),
            color=color,
        )
    )
    text = (
        alt.Chart(df, title="Confusion Matrix")
        .mark_text(baseline="middle", fontSize=25, fontWeight="bold")
        .encode(
            x="predicted:N",
            y="actual:N",
            text="confusion_matrix:N",
            color=alt.condition(
                alt.datum.confusion_matrix > 0,
                alt.value("black"),
                alt.value("black"),
            ),
        )
    )
    cm_chart = (heatmap + text).properties(width=600, height=480)
    return cm_chart