import pandas as pd
import numpy as np
import pickle
from scipy.spatial.distance import pdist, squareform

from scipy import stats
import statsmodels.stats.multitest as multitest

import optuna
from sklearn.manifold import TSNE

import xgboost as xgb
import lightgbm as lgb
import catboost as cb
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_val_predict, cross_validate

from .optuna_objective import *

import pprint
import seaborn as sns
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from plotly.graph_objs.layout import YAxis
import dtale
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from skbio.stats.composition import multiplicative_replacement, ancom
from skbio.diversity.alpha import chao1, shannon
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity


# Utils Functions
def stat_distribution(dd):

    stats = {"min": 0, "max": 0, "mean": 0, "median": 0}
    ft_counts = []
    dd_len = len(dd.index)
    for i in range(0, dd_len):
        ft_counts.append(sum(dd.iloc[i, :].apply(lambda x: 1 if x > 0 else 0)))
    stats["min"] = min(ft_counts)
    stats["max"] = max(ft_counts)
    stats["mean"] = np.mean(ft_counts)
    stats["median"] = np.median(ft_counts)
    stats["ft_count"] = ft_counts
    return stats


def pvals_stats(sorted_pvals, dd):

    p_stats = sorted_pvals

    for k, v in sorted_pvals.items():
        print(dd[k].mean())


def ancom_test(dd, labels):
    dd_adj = multiplicative_replacement(dd)

    dd_adj = pd.DataFrame(
        data=dd_adj,  # values
        index=dd.index,  # 1st column as index
        columns=dd.columns)
    results = ancom(dd_adj,
                    labels,
                    multiple_comparisons_correction='holm-bonferroni')
    results[0]["ft"] = dd.columns
    dtale.show(results[0].set_index(["ft"]),
               ignore_duplicate=True,
               notebook=True)
    return results


def features_stat(dd, labels):

    samples_pos_labels = labels[labels.surpoids > 0].index
    samples_neg_labels = labels[labels.surpoids < 1].index
    dd_pos = dd.loc[samples_pos_labels]
    dd_neg = dd.loc[samples_neg_labels]

    stats_pos = stat_distribution(dd_pos)
    stats_neg = stat_distribution(dd_neg)

    plt.figure(figsize=(12, 6))
    fig = sns.distplot(stats_pos["ft_count"], kde=True, color="red")
    fig = sns.distplot(stats_neg["ft_count"], kde=True, color="blue")
    fig.set(xlabel='Number of features', ylabel='Count')
    plt.show()

    pvals = []
    pvals_stats = {}
    fts = dd_pos.columns
    features_good = []
    for ft in fts:
        if (sorted(set(dd_pos[ft])) == sorted(set(dd_neg[ft]))) == True:
            pvals.append(1)
            continue
        features_good.append(ft)
        stat, p = stats.mannwhitneyu(dd_pos[ft], dd_neg[ft])
        mean_pos = dd_pos[ft].mean()
        std_pos = dd_pos[ft].std()
        mean_neg = dd_neg[ft].mean()
        std_neg = dd_neg[ft].std()
        pvals_stats[ft] = {
            "ft": ft,
            "pval": p,
            "mean_overweight": mean_pos,
            "std_overweight": std_pos,
            "mean_healthy": mean_neg,
            "std_healthy": std_neg,
            "fold_change": (mean_pos / mean_neg)
        }
        pvals.append(p)

    pvals_corrected = multitest.multipletests(pvals,
                                              alpha=0.05,
                                              method="fdr_bh")

    for i, _ft in enumerate(features_good):
        pvals_stats[_ft]["pval_good"] = pvals_corrected[0][i]
        pvals_stats[_ft]["pval_corrected"] = pvals_corrected[1][i]

    pvals_stats_df = pd.DataFrame.from_dict(pvals_stats, orient='index')
    dtale.show(pvals_stats_df.set_index(['ft']), notebook=True)

    # plot significant pvalues boxplot
    max_plot_number = 25
    i = 0
    pval_box_plot = {
        'feature': [],
        'Relative abundances': [],
        'Overweight': []
    }
    for ft in features_good:
        if pvals_stats[ft]["pval_good"] == True and i < max_plot_number:
            i += 1
            pval_box_plot['feature'].extend([ft] * len(dd[ft].values))
            pval_box_plot['Relative abundances'].extend(dd[ft].values)
            pval_box_plot['Overweight'].extend(labels['surpoids'].values)

    pvals_good_df_box = pd.DataFrame(pval_box_plot)
    fig = px.box(
        pvals_good_df_box,
        x='feature',
        y='Relative abundances',
        color='Overweight',
        color_discrete_sequence=['#EF553B', '#636EFA'],
    )
    fig.update_traces(
        quartilemethod="exclusive")  # or "inclusive", or "linear" by default
    #fig.update_layout(yaxis_type="log")
    fig.show()

    #return sorted_pvals


def alpha_diversity_chao1(dd, labels):
    _alpha_diversity = {}
    for index, row in dd.iterrows():
        # print(index)
        # print(chao1(row))
        _alpha_diversity[index] = chao1(row)

    _dd_alpha_diversity = pd.DataFrame.from_dict(_alpha_diversity,
                                                 orient='index',
                                                 columns=['chao1'])
    plot_boxplot(_dd_alpha_diversity, 'chao1', labels)


def alpha_diversity_shannon(dd, labels):
    _alpha_diversity = {}
    for index, row in dd.iterrows():
        _alpha_diversity[index] = shannon(row)

    _dd_alpha_diversity = pd.DataFrame.from_dict(_alpha_diversity,
                                                 orient='index',
                                                 columns=['shannon'])
    plot_boxplot(_dd_alpha_diversity, 'shannon', labels)


def beta_diversity_braycurtis(dd, labels):

    nb_colors = labels.nunique()
    _dm = beta_diversity('braycurtis', dd.to_numpy(), dd.index)
    _dm_pcoa = pcoa(_dm)
    _dm_pcoa.samples["labels"] = labels

    plt.figure(figsize=(8, 6))

    fig = px.scatter(_dm_pcoa.samples,
                     x='PC1',
                     y='PC2',
                     color=labels.sort_values(ascending=True),
                     title=("PCoA %s" % labels.name))

    fig.show()

    fig = px.scatter_3d(_dm_pcoa.samples,
                        x='PC1',
                        y='PC2',
                        z='PC3',
                        color=labels.sort_values(ascending=True),
                        title=("PCoA %s" % labels.name))
    fig.show()


def plot_boxplot(dd,
                 ft,
                 labels,
                 width=500,
                 height=400,
                 title="",
                 title_yaxis="",
                 show_legend=False,
                 y_max=0,
                 y_min=0):

    _title = title
    fig = px.box(dd,
                 y=ft,
                 color=labels.sort_values(ascending=True),
                 points="all",
                 title=_title,
                 width=width,
                 height=height)

    fig.update_layout(dict(boxgroupgap=0.5))
    fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black')
    fig.update_layout(legend_title=('%s') % labels.name)
    fig.update_layout(showlegend=show_legend)
    fig.update_layout(title={
        'text': _title,
        'y': 0.85,
        'x': 0.5,
        'xanchor': 'center',
        'yanchor': 'top'
    })
    if y_max != 0:
        fig.update_layout(yaxis=YAxis(autorange=False,
                                      range=[y_min, y_max],
                                      title=title_yaxis,
                                      showline=True,
                                      linewidth=2,
                                      linecolor='black',
                                      showgrid=True,
                                      gridwidth=1,
                                      gridcolor='lightgrey'))
    else:
        fig.update_layout(yaxis=YAxis(autorange=True,
                                      title=title_yaxis,
                                      showline=True,
                                      linewidth=2,
                                      linecolor='black',
                                      showgrid=True,
                                      gridwidth=1,
                                      gridcolor='lightgrey'))

    fig.show()


def plot_multiboxplot(dd, ft, labels, rows_nb, cols_nb):

    box = px.box(dd,
                 y=ft,
                 color=labels.sort_values(ascending=True),
                 points="all",
                 title=("Boxplot %s - %s" % (labels.name, ft)),
                 width=600,
                 height=400)

    fig = make_subplots(rows=rows_nb, cols=cols_nb)
    fig.add_box(box, row=1, col=1)
    fig.add_box(box, row=1, col=2)
    fig.show()


def build_dist_matrix(dd, dist):
    return pd.DataFrame(squareform(pdist(dd.iloc[:, 1:], dist)),
                        columns=dd.index,
                        index=dd.index)


def print_obj(obj):
    pp = pprint.PrettyPrinter(indent=2)
    pp.pprint(obj)


def tsne_plot_2d(dd, labels):

    nb_colors = labels.nunique()

    # tsne visualization of samples
    tsne = TSNE(n_components=2,
                verbose=0,
                perplexity=50,
                n_iter=5000,
                learning_rate=200)
    tsne_results = tsne.fit_transform(dd)
    _dd = dd.copy()
    _dd["labels"] = labels
    _dd['tsne-2d-one'] = tsne_results[:, 0]
    _dd['tsne-2d-two'] = tsne_results[:, 1]

    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x="tsne-2d-one",
        y="tsne-2d-two",
        hue="labels",
        palette=sns.color_palette("hls", nb_colors),
        data=_dd,
        legend="full",
        alpha=0.3,
    ).set_title('TSNE %s' % labels.name)
    plt.show()


def tsne_plot_3d(dd, labels):

    nb_colors = labels.nunique()

    # tsne visualization of samples
    tsne = TSNE(n_components=3,
                verbose=0,
                perplexity=50,
                n_iter=5000,
                learning_rate=200)
    tsne_results = tsne.fit_transform(dd)
    _dd = dd.copy()
    _dd["labels"] = labels
    _dd['tsne-3d-one'] = tsne_results[:, 0]
    _dd['tsne-3d-two'] = tsne_results[:, 1]
    _dd['tsne-3d-three'] = tsne_results[:, 2]

    plt.figure(figsize=(8, 8))
    fig = px.scatter_3d(_dd,
                        x='tsne-3d-one',
                        y='tsne-3d-two',
                        z='tsne-3d-three',
                        color=labels,
                        title=("TSNE %s" % labels.name))

    fig.show()


def model_scores(data, model, cv):
    print(" # Computing model scores")
    scoring = [
        "accuracy", "balanced_accuracy", "average_precision", "f1",
        "precision", "recall", "roc_auc"
    ]
    scores = {}
    kfold = StratifiedKFold(n_splits=cv, random_state=None, shuffle=True)
    scores_cross_val = cross_validate(
        model,
        data['x'],
        data['y'].ravel(),
        cv=kfold,
        scoring=("accuracy", "balanced_accuracy", "average_precision", "f1",
                 "precision", "recall", "roc_auc"),
        n_jobs=-1)
    return scores_cross_val


def optuna_xgboost_cv(dd, labels, imbalance_ratio, n_trials):
    print(" # Optuna parameters search")
    data = {
        'x': dd.values,
        'y': labels.values,
    }
    optuna.logging.set_verbosity(optuna.logging.CRITICAL)
    pruner = optuna.pruners.MedianPruner(n_warmup_steps=5)
    objective = Objective_xgboost_cv(data)
    study = optuna.create_study(pruner=pruner, direction='maximize')
    study.optimize(objective, n_trials=n_trials, n_jobs=-1)
    print(" # Optuna best trial score")
    print_obj(study.best_trial.value)
    print(" # Optuna best params")
    print_obj(study.best_params)
    model = xgb.XGBClassifier(**study.best_params,
                              scale_pos_weight=imbalance_ratio)
    results = model_scores(data, model, 10)
    print_obj(results)
    for r, a in results.items():
        print("%s : %f" % (r, a.mean()))

    y_pred = cross_val_predict(model, data['x'], data['y'].ravel(), cv=10)
    print(" # Confusion matrix")
    print_obj(confusion_matrix(data['y'].ravel(), y_pred))

    return model


def optuna_xgboost_accuracy(dd, labels, imbalance_ratio, n_trials):
    print(" # Optuna parameters search")
    data = {
        'x': dd.values,
        'y': labels.values,
    }
    optuna.logging.set_verbosity(optuna.logging.CRITICAL)
    pruner = optuna.pruners.MedianPruner(n_warmup_steps=5)
    objective = Objective_xgboost_accuracy(data)
    study = optuna.create_study(pruner=pruner, direction='maximize')
    study.optimize(objective, n_trials=n_trials, n_jobs=-1)
    print(" # Optuna best trial score")
    print_obj(study.best_trial.value)
    print(" # Optuna best params")
    print_obj(study.best_params)
    model = xgb.XGBClassifier(**study.best_params,
                              scale_pos_weight=imbalance_ratio)
    results = model_scores(data, model, 10)
    print_obj(results)
    for r, a in results.items():
        print("%s : %f" % (r, a.mean()))

    y_pred = cross_val_predict(model, data['x'], data['y'].ravel(), cv=10)
    print(" # Confusion matrix")
    print_obj(confusion_matrix(data['y'].ravel(), y_pred))

    return model


def optuna_lightgbm_accuracy(dd, labels, imbalance_ratio, n_trials):
    print(" # Optuna parameters search")
    data = {
        'x': dd.values,
        'y': labels.values,
    }
    optuna.logging.set_verbosity(optuna.logging.CRITICAL)
    pruner = optuna.pruners.MedianPruner(n_warmup_steps=5)
    objective = Objective_lightgbm_accuracy(data)
    study = optuna.create_study(pruner=pruner, direction='maximize')
    study.optimize(objective, n_trials=n_trials, n_jobs=-1)
    print(" # Optuna best trial score")
    print_obj(study.best_trial.value)
    print(" # Optuna best params")
    print_obj(study.best_params)
    model = lgb.LGBMClassifier(**study.best_params,
                               scale_pos_weight=imbalance_ratio)
    results = model_scores(data, model, 10)
    print_obj(results)
    for r, a in results.items():
        print("%s : %f" % (r, a.mean()))

    y_pred = cross_val_predict(model, data['x'], data['y'].ravel(), cv=10)
    print(" # Confusion matrix")
    print_obj(confusion_matrix(data['y'].ravel(), y_pred))

    return model


def optuna_catboost_accuracy(dd, labels, imbalance_ratio, n_trials):
    print(" # Optuna parameters search")
    data = {
        'x': dd.values,
        'y': labels.values,
    }
    optuna.logging.set_verbosity(optuna.logging.CRITICAL)
    pruner = optuna.pruners.MedianPruner(n_warmup_steps=5)
    objective = Objective_catboost_accuracy(data)
    study = optuna.create_study(pruner=pruner, direction='maximize')
    study.optimize(objective, n_trials=n_trials, n_jobs=-1)
    print(" # Optuna best trial score")
    print_obj(study.best_trial.value)
    print(" # Optuna best params")
    print_obj(study.best_params)
    model = cb.CatBoostClassifier(**study.best_params,
                                  scale_pos_weight=imbalance_ratio,
                                  silent=True)
    results = model_scores(data, model, 10)
    print_obj(results)
    for r, a in results.items():
        print("%s : %f" % (r, a.mean()))

    y_pred = cross_val_predict(model, data['x'], data['y'].ravel(), cv=10)
    print(" # Confusion matrix")
    print_obj(confusion_matrix(data['y'].ravel(), y_pred))

    return model


def optuna_RF_accuracy(dd, labels, imbalance_ratio, n_trials):
    print(" # Optuna parameters search")
    data = {
        'x': dd.values,
        'y': labels.values,
    }
    optuna.logging.set_verbosity(optuna.logging.CRITICAL)
    pruner = optuna.pruners.MedianPruner(n_warmup_steps=5)
    objective = Objective_RF_accuracy(data)
    study = optuna.create_study(pruner=pruner, direction='maximize')
    study.optimize(objective, n_trials=n_trials, n_jobs=-1)
    print(" # Optuna best trial score")
    print_obj(study.best_trial.value)
    print(" # Optuna best params")
    print_obj(study.best_params)
    class_weight = None
    if 'imbalance_ratio' != 1:
        class_weight = "balanced"
    model = RandomForestClassifier(**study.best_params,
                                   n_estimators=100,
                                   class_weight=class_weight)
    results = model_scores(data, model, 10)
    print_obj(results)
    for r, a in results.items():
        print("%s : %f" % (r, a.mean()))

    y_pred = cross_val_predict(model, data['x'], data['y'].ravel(), cv=10)
    print(" # Confusion matrix")
    print_obj(confusion_matrix(data['y'].ravel(), y_pred))

    return model


def optuna_SVC_accuracy(dd, labels, imbalance_ratio, n_trials):
    print(" # Optuna parameters search")
    data = {
        'x': dd.values,
        'y': labels.values,
    }
    optuna.logging.set_verbosity(optuna.logging.CRITICAL)
    pruner = optuna.pruners.MedianPruner(n_warmup_steps=5)
    objective = Objective_SVC_accuracy(data)
    study = optuna.create_study(pruner=pruner, direction='maximize')
    study.optimize(objective, n_trials=n_trials, n_jobs=-1)
    print(" # Optuna best trial score")
    print_obj(study.best_trial.value)
    print(" # Optuna best params")
    print_obj(study.best_params)
    class_weight = None
    if 'imbalance_ratio' != 1:
        class_weight = "balanced"
    model = SVC(**study.best_params, class_weight=class_weight)
    results = model_scores(data, model, 10)
    print_obj(results)
    for r, a in results.items():
        print("%s : %f" % (r, a.mean()))

    y_pred = cross_val_predict(model, data['x'], data['y'].ravel(), cv=10)
    print(" # Confusion matrix")
    print_obj(confusion_matrix(data['y'].ravel(), y_pred))

    return model


def optuna_viz(study):
    optuna.visualization.plot_intermediate_values(study)
    optuna.visualization.plot_optimization_history(study)
    #optuna.visualization.plot_contour(study, params=['max_depth', 'alpha', 'gamma'])
    #optuna.visualization.plot_parallel_coordinate(study, params=['max_depth', 'booster'])
    #optuna.visualization.plot_slice(study, params=['max_depth', 'booster'])
    #plt.show()
