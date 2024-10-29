import os
from glob import glob
from pathlib import Path
import tifffile
from tqdm import tqdm

import pandas as pd
import numpy as np

import scanpy as sc
import anndata

from skimage.exposure import adjust_gamma
from skimage import filters
import scipy.ndimage as ndi
import scipy

import seaborn as sns
import matplotlib.pyplot as plt

import pathlib
import anndata
import warnings
warnings.simplefilter("ignore", UserWarning)

import yaml
import matplotlib
sc.settings.set_figure_params(dpi=200, dpi_save=300, fontsize=12)
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["axes.grid"] = False
matplotlib.rcParams['figure.dpi'] = 300
#matplotlib.use('Agg')


# measure mean intensity on patient level
metadata_filename = 'metadata/ggo_config.yml'

with open(metadata_filename, "r") as stream:
    try:
        metadata = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        
conditions = metadata['CONDITIONS'] + ['PID', 'Race', 'Smoking', 'Packyrs', 'Age', 'Molecular']

adata_dict = dict()
combined = dict()
for p in ['PANEL_G', 'PANEL_H']:

    n_data = sc.read(metadata[p]['AnnData']['patient_normal_expression_name'])
    t_data = sc.read(metadata[p]['AnnData']['patient_tumor_expression_name'])

    adata_dict[f'{p}_normal'] = n_data
    adata_dict[f'{p}_tumor'] = t_data
    

    n_data.var['Panel'] = p
    t_data.var['Panel'] = p
    # adata_dict[p] = adata_dict[p][:,~adata_dict[p].var.index.str.contains('EMPTY')]
    combined[p] = anndata.concat([n_data, t_data], axis = 0)
    

# adata_dict
intersection = combined['PANEL_G'].obs.index.intersection(combined['PANEL_H'].obs.index)
combined['PANEL_G'] = combined['PANEL_G'][intersection]
combined['PANEL_H'] = combined['PANEL_H'][intersection]

expression = anndata.concat(
    [combined['PANEL_G'], combined['PANEL_H']],
    axis = 1
)

expression.obs = combined['PANEL_G'].obs
expression.var_names_make_unique()
# Running sanity checks
expression.var.index.sort_values()

'''
---------------------
Index(['41BB(Yb172)', 'CCR4(Sm149)', 'CCR7(Gd158)', 'CD103(Nd150)',
       'CD117(Dy164)', 'CD11b(Dy164)', 'CD11c(Sm154)', 'CD14(Nd144)',
       'CD14(Nd144)-1', 'CD15(Gd158)', 'CD16(Eu151)', 'CD16(Nd146)',
       'CD163(Sm147)', 'CD163(Sm147)-1', 'CD20(Dy161)', 'CD20(Dy161)-1',
       'CD206(Er167)', 'CD25(Lu175)', 'CD27(Yb171)', 'CD28(Gd160)',
       'CD3(Er170)', 'CD3(Er170)-1', 'CD31(Eu151)', 'CD31(Nd143)',
       'CD33(Nd145)', 'CD4(Gd156)', 'CD4(Gd156)-1', 'CD45(Sm152)',
       'CD45RA(Er166)', 'CD45RO(Yb173)', 'CD56(Dy163)', 'CD56(Tm169)',
       'CD57(Tm169)', 'CD66b(Sm152)', 'CD68(Tb159)', 'CD68(Tb159)-1',
       'CD8a(Dy162)', 'CD8a(Dy162)-1', 'CTLA4(Yb176)', 'FOLR1(Lu175)',
       'FoxP3(Gd155)', 'GLUT5(Yb171)', 'GranzymeB(Nd142)', 'HLAABC(Er167)',
       'HLADR(Pr141)', 'HLADR(Pr141)-1', 'ICOS(Nd148)', 'IFNg(Gd160)',
       'IL12p40(Ho165)', 'IL17A(Yb174)', 'IL1R1(Eu153)', 'IL1alpha(Yb172)',
       'IL1beta(Er166)', 'IL23R(Sm149)', 'IL23p19(Sm154)', 'Ki67(Er168)',
       'LAG3(Eu153)', 'MMP7(Yb176)', 'NKp44(Nd142)', 'PD1(Ho165)',
       'PDL1(Nd150)', 'PanCytokeratin(Nd148)', 'PanCytokeratin(Pt198)',
       'RAGE(Gd155)', 'SCL5A2SGLT2(Dy163)', 'SFTPC(Er168)', 'TIM3(Yb174)',
       'Tbet(Nd145)', 'VISTA(Yb173)', 'Vimentin(Nd143)', 'aSMA(Nd146)',
       'aSMA(Pt196)'],
      dtype='object')
---------------------
Duplicated features:
[
    ('CD14(Nd144)', 'CD14(Nd144)-1'),
    ('CD16(Eu151)', 'CD16(Nd146)'),
    ('CD163(Sm147)', 'CD163(Sm147)-1'),
    ('CD3(Er170)', 'CD3(Er170)-1',),
    ('CD20(Dy161)', 'CD20(Dy161)-1'),
    ('CD31(Eu151)', 'CD31(Nd143)'),
    ('CD56(Dy163)', 'CD56(Tm169)'),
    ('CD68(Tb159)', 'CD68(Tb159)-1'),
    ('CD8a(Dy162)', 'CD8a(Dy162)-1',),
    ('HLADR(Pr141)', 'HLADR(Pr141)-1',),
    ('PanCytokeratin(Nd148)', 'PanCytokeratin(Pt198)'),
    ('aSMA(Nd146)', 'aSMA(Pt196)'),
]
'''

# Show through scatterplot that these variables are highly correlated
dup_features = [
    ('CD14(Nd144)', 'CD14(Nd144)-1'),
    ('CD16(Eu151)', 'CD16(Nd146)'),
    ('CD163(Sm147)', 'CD163(Sm147)-1'),
    ('CD3(Er170)', 'CD3(Er170)-1',),
    ('CD20(Dy161)', 'CD20(Dy161)-1'),
    ('CD31(Eu151)', 'CD31(Nd143)'),
    ('CD56(Dy163)', 'CD56(Tm169)'),
    ('CD68(Tb159)', 'CD68(Tb159)-1'),
    ('CD8a(Dy162)', 'CD8a(Dy162)-1',),
    ('HLADR(Pr141)', 'HLADR(Pr141)-1',),
    ('PanCytokeratin(Nd148)', 'PanCytokeratin(Pt198)'),
    ('aSMA(Nd146)', 'aSMA(Pt196)'),
]



expression.obs['radio'] = pd.Categorical(expression.obs['radio'], categories = ['N', 'PNS', 'PS', 'S', 'UNK'])
df = expression.to_df()
os.makedirs('figures/colinearity/', exist_ok = True)
df['radio'] = expression.obs['radio']

# sanity check for duplicate markers
for f1, f2 in dup_features:

    # sc.pl.scatter(expression, x = f1, y = f2, size = 20)
    # plt.close()

    g = sns.lmplot(data=df, x=f1, y=f2, scatter_kws={"s": 10}) # hue = marker

    def annotate(data, **kws):
        import scipy as sp
        r, p = sp.stats.pearsonr(data[f1], data[f2])
        ax = plt.gca()
        ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
                transform=ax.transAxes)
    
    g.map_dataframe(annotate)
    plt.tight_layout()
    plt.show()
    plt.savefig(f'figures/colinearity/{f1}_{f2}.pdf', bbox_inches = 'tight')
    plt.close()

    # sc.pl.heatmap(
    #     expression,
    #     var_names = expression.var.index,
    #     groupby = 'radio',
    #     log = True,
    #     standard_scale = 'var'
    # )



# preprocessing
e = expression.copy()
e.X = np.expm1(e.X)
sc.pp.scale(e, max_value = 3)
e.X[e.X < -3] = -3
df = e.to_df()

for f1, f2 in dup_features:

    # sc.pl.scatter(expression, x = f1, y = f2, size = 20)
    # plt.close()

    g = sns.lmplot(data=df, x=f1, y=f2, scatter_kws={"s": 10}) # hue = marker

    def annotate(data, **kws):
        import scipy as sp
        r, p = sp.stats.pearsonr(data[f1], data[f2])
        ax = plt.gca()
        ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p),
                transform=ax.transAxes)
    
    g.map_dataframe(annotate)
    plt.tight_layout()
    plt.show()


# plt.figure(figsize=(20, 15), dpi = 150)
# # define the mask to set the values in the upper triangle to True
# mask = np.triu(np.ones_like(e.to_df().corr(), dtype=np.bool))
# heatmap = sns.heatmap(e.to_df().corr(), mask=mask, vmin=-1, vmax=1, annot=True, cmap='BrBG')
# heatmap.set_title('Triangle Correlation Heatmap', fontdict={'fontsize':18}, pad=16);
# plt.show()

# sc.pl.heatmap(e, e.var.index, groupby='radio', cmap='viridis', dendrogram=True)


# TO DO LIST
# Heatmap stratified per condition
    # Statistical testing per condition (Wil-Cox?)
# Random Forest stratified per condition
# Gene expression comparison


def plot_importance_bar(
    feature_importance: pd.Series,
    std: pd.Series,
    save_dir: pathlib.Path = None,
    show_fig: bool = False,
):
    fig, ax = plt.subplots(1,1,figsize =(6,4), dpi = 300)
    feature_importance.plot.bar(yerr=std.loc[feature_importance.index], ax=ax)
    ax.set_title("Feature importances using MDI")
    ax.set_ylabel("Mean decrease in impurity")
    fig.tight_layout()
    
    if show_fig:
        plt.show()
    if save_dir:
        plt.savefig(save_dir, bbox_inches = 'tight')
    plt.close()

    return fig

def report_permutation_importance(
    importances, features
) -> None:
    for metric in importances:
        print(f"{metric}")
        r = importances[metric]
        for i in r.importances_mean.argsort()[::-1]:
            if r.importances_mean[i] - r.importances_std[i] > 0:
                print(f"    {features[i]:<8}"
                      f"{r.importances_mean[i]:.3f}"
                      f" +/- {r.importances_std[i]:.3f}")


key_maps = {
    'pathology': {
        'Normal': 0,
        'AIS': 1,
        'MIA': 2,
        'IAC': 3,
    },
    'radio': {
        'N': 0,
        'PNS': 1,
        'PS': 2,
        'S': 3
    },
    'Smoking': {
        'No': 0,
        'Yes': 1
    },
    'group': {
        'N': 0,
        'T': 1
    }
}

e.obs.loc[e.obs['Smoking']=='No','Packyrs'] = 0
e.obs['Packyrs'] = e.obs['Packyrs'].fillna(e.obs.loc[e.obs['Smoking']=="Yes", 'Packyrs'].median())

# sample stratification for prediction model and to impute feature importances
# Compare using ridge regression, random forests, and permutation importances

from sklearn.model_selection import train_test_split
X, y = pd.DataFrame(e.X, index = e.obs.index, columns = e.var.index), e.obs
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 42) #stratify = y, 


for f in ['pathology', 'radio']:
    y_train[f'{f}_num'] = y_train[f].replace(key_maps[f])
    y_test[f'{f}_num'] = y_test[f].replace(key_maps[f])

numeric_features = ['radio_num', 'Age', 'Packyrs', 'pathology_num',]
categorical_features = ['group','Smoking']

# Random forest based feature importance
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

forest_score = pd.DataFrame()
forest_imp = dict()
perm_imp = dict()
all_imp = dict()

for feature in numeric_features:
    X_train_ = X_train
    X_test_ = X_test
    y_train_ = y_train[feature]
    y_test_ = y_test[feature]
    if y_train[feature].isna().sum() > 0:
        y_train_ = y_train[feature].fillna(y_train[feature].mean())


    # filter UNK for radiological feature
    if feature == 'radio_num':
        train_idx_drop = y_train_[y_train_=='UNK'].index
        test_idx_drop = y_test_[y_test_=='UNK'].index

        X_train_ = X_train.loc[~X_train.index.isin(train_idx_drop),:]
        X_test_ = X_test.loc[~X_test.index.isin(test_idx_drop),:]
        y_train_ = y_train_.loc[~y_train_.index.isin(train_idx_drop)]
        y_test_ = y_test_.loc[~y_test_.index.isin(test_idx_drop)]
        

    forest = RandomForestRegressor(random_state=0)
    forest.fit(X_train_, y_train_)

    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)

    forest_importances = pd.Series(importances, index=X_train_.columns)
    std = pd.Series(std, index = X_train_.columns)

    forest_importances = forest_importances.sort_values(ascending = False)
    std = std.loc[forest_importances.index]

    forest_imp[feature] = pd.DataFrame({
        'importances_mean':forest_importances,
        'importances_std':std,
    })

    score_ = pd.Series({
        'train':forest.score(X_train_, y_train_),
        'test': forest.score(X_test_, y_test_),
    })

    print(feature)
    print(score_.round(2))
    forest_score[feature] = score_

    # Permutation based feature importances 
    from sklearn.inspection import permutation_importance
    scoring = ['r2']
    r_multi = permutation_importance(
        forest, X_test_, y_test_, n_repeats=100, random_state=0, scoring=scoring)

    report_permutation_importance(r_multi, X_train_.columns)

    del r_multi['r2']['importances']
    perm_imp[feature] = pd.DataFrame(r_multi['r2'], index = X_train_.columns)

for feature in categorical_features:
    y_train_ = y_train[feature].replace(key_maps[feature])
    y_test_ = y_test[feature].replace(key_maps[feature])

    # Random forest based feature importance
    forest = RandomForestClassifier(random_state=0)
    forest.fit(X_train, y_train_)

    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)

    forest_importances = pd.Series(importances, index=X_train.columns)
    std = pd.Series(std, index = X_train.columns)

    forest_importances = forest_importances.sort_values(ascending = False)
    std = std.loc[forest_importances.index]
    fimp = forest_importances.iloc[:10]
    
    # plot_importance_bar(fimp, std, show_fig = True)

    forest_imp[feature] = pd.DataFrame({
        'importances_mean':forest_importances,
        'importances_std':std,
    })

    score_ = pd.Series({
        'train':forest.score(X_train, y_train_),
        'test': forest.score(X_test, y_test_),
    })

    print(feature)
    print(score_.round(2))
    forest_score[feature] = score_

    # Permutation based feature importances
    from sklearn.inspection import permutation_importance
    scoring = ['accuracy']
    r_multi = permutation_importance(
        forest, X_test, y_test_, n_repeats=100, random_state=0, scoring=scoring)

    report_permutation_importance(r_multi, X_train.columns)


    del r_multi['accuracy']['importances']
    perm_imp[feature] = pd.DataFrame(r_multi['accuracy'], index = X_train_.columns)

forest_score = forest_score.T
for feature in forest_imp:
    forest_imp[feature].columns = 'rf_' + forest_imp[feature].columns
    perm_imp[feature].columns = 'perm_' + perm_imp[feature].columns

for feature in forest_imp:
    all_imp[feature] = pd.concat([forest_imp[feature], perm_imp[feature]], axis = 1)


k = 5

fig, axes = plt.subplots(2,2,figsize = (7,4), dpi = 300)

features = ['group', 'radio_num', 'pathology_num', 'Smoking']

for i, ax in enumerate(axes.flatten()):
    df = all_imp[features[i]]

    df['importance_sum'] = df[['perm_importances_mean', 'rf_importances_mean']].mean(axis=1)
    top_k_idx = df['importance_sum'].sort_values(ascending = False).index[:k]
    df_k = df.loc[top_k_idx]
    df_ = df.loc[~df.index.isin(top_k_idx)]

    ax.scatter(
        x = df_k['perm_importances_mean'],
        y = df_k['rf_importances_mean'],
        c = 'k',
        s = 5,
    )

    ax.scatter(
        x = df_['perm_importances_mean'],
        y = df_['rf_importances_mean'],
        c = 'gray',
        s = 5
    )

    ax.set_title(f"{features[i]}, R2: {round(forest_score.loc[features[i], 'test'],4)}")
    # ax.set_xlabel('Permutation score')
    # ax.set_ylabel('Random Forest impurity')

    padding = (df['perm_importances_mean'].max() - df['perm_importances_mean'].min())/100

    for idx in top_k_idx:
        ax.text(
            x = df.loc[idx, 'perm_importances_mean'] + padding,
            y = df.loc[idx, 'rf_importances_mean'],
            s = idx,
            verticalalignment = 'center',
            fontsize = 8)

plt.tight_layout()
#plt.show()
plt.savefig('figures/feature_importance_prediction.pdf')
plt.close()


# Molecular Prediction
e.obs.loc[e.obs['Molecular'].str.lower() == 'not performed', 'Molecular'] = np.nan
e.obs.loc[e.obs['Molecular'].str.lower() == 'nor performed', 'Molecular'] = np.nan
## mutational status
e.obs['EGFR'] = e.obs['Molecular'].str.upper().str.contains('EGFR').astype(float)
e.obs['KRAS'] = e.obs['Molecular'].str.upper().str.contains('KRAS').astype(float)
e.obs['TP53'] = e.obs['Molecular'].str.upper().str.contains('TP53').astype(float)

idx = e.obs[['EGFR', 'KRAS', 'TP53']].isna().any(axis = 1)
molecular = e[idx[~idx].index]

X, y = pd.DataFrame(molecular.X, index = molecular.obs.index, columns = molecular.var.index), molecular.obs
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state = 42) #stratify = y, 

features = ['KRAS', 'TP53', 'EGFR']
forest_score_molecular = pd.DataFrame()
perm_imp_molecular = dict()
forest_imp_molecular = dict()
all_imp_molecular = dict()

for feature in features:
    y_train_ = y_train[feature]
    y_test_ = y_test[feature]

    # Random forest based feature importance
    forest = RandomForestClassifier(random_state=0, max_depth = 5)
    forest.fit(X_train, y_train_)

    importances = forest.feature_importances_
    std = np.std([tree.feature_importances_ for tree in forest.estimators_], axis=0)

    forest_importances = pd.Series(importances, index=X_train.columns)
    std = pd.Series(std, index = X_train.columns)

    forest_importances = forest_importances.sort_values(ascending = False)
    std = std.loc[forest_importances.index]
    fimp = forest_importances.iloc[:10]
    
    # plot_importance_bar(fimp, std, show_fig = True)

    forest_imp_molecular[feature] = pd.DataFrame({
        'importances_mean':forest_importances,
        'importances_std':std,
    })

    score_ = pd.Series({
        'train':forest.score(X_train, y_train_),
        'test': forest.score(X_test, y_test_),
    })

    print(feature)
    print(score_.round(2))
    forest_score_molecular[feature] = score_

    # Permutation based feature importances
    from sklearn.inspection import permutation_importance
    scoring = ['accuracy']
    r_multi = permutation_importance(
        forest, X_test, y_test_, n_repeats=30, random_state=0, scoring=scoring)

    report_permutation_importance(r_multi, X_train.columns)

    del r_multi['accuracy']['importances']
    perm_imp_molecular[feature] = pd.DataFrame(r_multi['accuracy'], index = X_train.columns)
    print(pd.DataFrame(r_multi['accuracy'], index = X_train.columns))


forest_score_molecular = forest_score_molecular.T
for feature in forest_imp_molecular:
    forest_imp_molecular[feature].columns = 'rf_' + forest_imp_molecular[feature].columns
    perm_imp_molecular[feature].columns = 'perm_' + perm_imp_molecular[feature].columns

for feature in forest_imp_molecular:
    all_imp_molecular[feature] = pd.concat([forest_imp_molecular[feature], perm_imp_molecular[feature]], axis = 1)


k = 5

fig, axes = plt.subplots(1,3,figsize = (10,3), dpi = 300)


for i, ax in enumerate(axes.flatten()):
    df = all_imp_molecular[features[i]]

    df['importance_sum'] = df[['perm_importances_mean', 'rf_importances_mean']].mean(axis=1)
    top_k_idx = df['importance_sum'].sort_values(ascending = False).index[:k]
    df_k = df.loc[top_k_idx]
    df_ = df.loc[~df.index.isin(top_k_idx)]

    ax.scatter(
        x = df_k['perm_importances_mean'],
        y = df_k['rf_importances_mean'],
        c = 'k',
        s = 5,
    )

    ax.scatter(
        x = df_['perm_importances_mean'],
        y = df_['rf_importances_mean'],
        c = 'gray',
        s = 5
    )

    ax.set_title(f"{features[i]}, Accuracy: {round(forest_score_molecular.loc[features[i], 'test'],4) * 100}%")
    # ax.set_xlabel('Permutation score')
    # ax.set_ylabel('Random Forest impurity')

    padding = (df['perm_importances_mean'].max() - df['perm_importances_mean'].min())/100

    for idx in top_k_idx:
        ax.text(
            x = df.loc[idx, 'perm_importances_mean'] + padding,
            y = df.loc[idx, 'rf_importances_mean'],
            s = idx,
            verticalalignment = 'center',
            fontsize = 8)

plt.tight_layout()
#plt.show()
plt.savefig('figures/feature_importance_prediction_molecular.pdf', bbox_inches = 'tight')
plt.close()



all_imp, all_imp_molecular
feats = ['radio_num', 'pathology_num', 'group', 'Smoking'] + ['KRAS', 'TP53']

dfs = []
for f in feats:
    if f in all_imp:
        df = all_imp[f]
    elif f in all_imp_molecular:
        df = all_imp_molecular[f]

    df = df[['rf_importances_mean', 'rf_importances_std', 'perm_importances_mean', 'perm_importances_std']]
    df.columns = f + '_' + df.columns
    dfs.append(df)

imp = pd.concat(dfs, axis = 1)
imp['importance_sum'] = imp[imp.columns[imp.columns.str.contains('mean')]].sum(axis=1)
imp['importance_sum'].sort_values(ascending = False)

imp_data = pd.DataFrame([imp[imp.columns[imp.columns.str.contains(f) & imp.columns.str.contains('mean')]].sum(axis = 1) for f in feats])
imp_data.index = ['Radiology', 'Pathology', 'Tumor Status', 'Smoking', 'KRAS', 'TP53']
imp_data = imp_data.T
imp_data['importances_sum'] = imp_data.sum(axis=1)
imp_data = imp_data.sort_values('importances_sum', ascending = False)
del imp_data['importances_sum']
imp_data = imp_data.clip(lower=0)

imp_data.to_csv('metadata/feature_importances.csv')

imp_data = pd.read_csv('metadata/feature_importances.csv', index_col = 0)
from plotly.subplots import make_subplots
import plotly.graph_objects as go

row, col = 4, 5

fig = make_subplots(
    rows=4, cols=5,
    subplot_titles = imp_data.index[:row*col],
    specs=[[{'type': 'polar'}] * col] * row,
    horizontal_spacing = 0,
    vertical_spacing = 0.05)

subplots = []
for i in range(row * col):
    idx = imp_data.index[i]
    
    # px.line_polar(, 
    #     theta='index', line_close=True, 
        
    #     title='Average Pokemon Attributes by Primary Type',
    #     template = 'ggplot2',
    #     ax = ax
    # )
    
    fig.add_trace(
        go.Scatterpolar(
            r = imp_data.loc[idx].tolist() + [imp_data.loc[idx, 'Radiology']],
            theta = imp_data.columns.tolist() + ['Radiology'],
            # range_r=(0,0.4),
            # title = imp_data.loc[idx].name,
            name = idx,
            line_color = 'gray',
            mode = 'lines + markers',
        ),
        row = i//col + 1,
        col = i%col + 1
    )

fig.update_traces(fill='toself')
fig.update_layout(
    polar = dict(radialaxis_range = [0,0.4]),
    polar2 = dict(radialaxis_range = [0,0.4]),
    polar3 = dict(radialaxis_range = [0,0.2]),
    polar4 = dict(radialaxis_range = [0,0.1]),
    polar5 = dict(radialaxis_range = [0,0.1]),
    polar6 = dict(radialaxis_range = [0,0.1]),
    polar7 = dict(radialaxis_range = [0,0.1]),
    polar8 = dict(radialaxis_range = [0,0.1]),
    polar9 = dict(radialaxis_range = [0,0.1]),
    polar10 = dict(radialaxis_range = [0,0.05]),
    polar11 = dict(radialaxis_range = [0,0.05]),
    polar12 = dict(radialaxis_range = [0,0.05]),
    polar13 = dict(radialaxis_range = [0,0.05]),
    polar14 = dict(radialaxis_range = [0,0.05]),
    polar15 = dict(radialaxis_range = [0,0.05]),
    polar16 = dict(radialaxis_range = [0,0.05]),
    polar17 = dict(radialaxis_range = [0,0.05]),
    polar18 = dict(radialaxis_range = [0,0.05]),
    polar19 = dict(radialaxis_range = [0,0.05]),
    polar20 = dict(radialaxis_range = [0,0.05]),
    showlegend=False,
    height=1200, width=2000,
)
fig.write_image("figures/feature_importance_radial_plot.pdf")

fig.show()
