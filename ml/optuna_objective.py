import numpy as np
import sklearn.datasets
import sklearn.metrics
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.svm import SVC
import xgboost as xgb
import lightgbm as lgb
import catboost as cb

import optuna


class Objective_xgboost_accuracy(object):
    def __init__(self, data):
        self.data = data

    def __call__(self, trial):
        x, y = self.data['x'], self.data['y']

        imbalance_ratio = 1
        if 'imbalance_ratio' in self.data:
            imbalance_ratio = self.data['imbalance_ratio']

        train_x, test_x, train_y, test_y = train_test_split(x,
                                                            y,
                                                            test_size=0.25)
        dtrain = xgb.DMatrix(train_x, label=train_y)
        dtest = xgb.DMatrix(test_x, label=test_y)

        param = {
            'silent':
            1,
            'objective':
            'binary:logistic',
            'eval_metric':
            'auc',
            'booster':
            trial.suggest_categorical('booster',
                                      ['gbtree', 'gblinear', 'dart']),
            'lambda':
            trial.suggest_loguniform('lambda', 1e-8, 1.0),
            'alpha':
            trial.suggest_loguniform('alpha', 1e-8, 1.0),
            'scale_pose_weight':
            imbalance_ratio
        }

        if param['booster'] == 'gbtree' or param['booster'] == 'dart':
            param['max_depth'] = trial.suggest_int('max_depth', 1, 9)
            param['eta'] = trial.suggest_loguniform('eta', 1e-8, 1.0)
            param['gamma'] = trial.suggest_loguniform('gamma', 1e-8, 1.0)
            param['grow_policy'] = trial.suggest_categorical(
                'grow_policy', ['depthwise', 'lossguide'])
        if param['booster'] == 'dart':
            param['sample_type'] = trial.suggest_categorical(
                'sample_type', ['uniform', 'weighted'])
            param['normalize_type'] = trial.suggest_categorical(
                'normalize_type', ['tree', 'forest'])
            param['rate_drop'] = trial.suggest_loguniform(
                'rate_drop', 1e-8, 1.0)
            param['skip_drop'] = trial.suggest_loguniform(
                'skip_drop', 1e-8, 1.0)

        # Add a callback for pruning.
        pruning_callback = optuna.integration.XGBoostPruningCallback(
            trial, 'validation-auc')

        bst = xgb.train(param,
                        dtrain,
                        evals=[(dtest, 'validation')],
                        callbacks=[pruning_callback],
                        verbose_eval=False)

        preds = bst.predict(dtest)
        pred_labels = np.rint(preds)
        accuracy = sklearn.metrics.accuracy_score(test_y, pred_labels)
        return accuracy


class Objective_xgboost_cv(object):
    def __init__(self, data):
        self.data = data

    def __call__(self, trial):
        x, y = self.data['x'], self.data['y']
        dtrain = xgb.DMatrix(x, label=y)

        imbalance_ratio = 1
        if 'imbalance_ratio' in self.data:
            imbalance_ratio = self.data['imbalance_ratio']

        param = {
            'silent':
            1,
            'objective':
            'binary:logistic',
            'eval_metric':
            'auc',
            'booster':
            trial.suggest_categorical('booster',
                                      ['gbtree', 'gblinear', 'dart']),
            'lambda':
            trial.suggest_loguniform('lambda', 1e-8, 1.0),
            'alpha':
            trial.suggest_loguniform('alpha', 1e-8, 1.0),
            'scale_pose_weight':
            imbalance_ratio
        }

        if param['booster'] == 'gbtree' or param['booster'] == 'dart':
            param['max_depth'] = trial.suggest_int('max_depth', 1, 9)
            param['eta'] = trial.suggest_loguniform('eta', 1e-8, 1.0)
            param['gamma'] = trial.suggest_loguniform('gamma', 1e-8, 1.0)
            param['grow_policy'] = trial.suggest_categorical(
                'grow_policy', ['depthwise', 'lossguide'])
        if param['booster'] == 'dart':
            param['sample_type'] = trial.suggest_categorical(
                'sample_type', ['uniform', 'weighted'])
            param['normalize_type'] = trial.suggest_categorical(
                'normalize_type', ['tree', 'forest'])
            param['rate_drop'] = trial.suggest_loguniform(
                'rate_drop', 1e-8, 1.0)
            param['skip_drop'] = trial.suggest_loguniform(
                'skip_drop', 1e-8, 1.0)

        pruning_callback = optuna.integration.XGBoostPruningCallback(
            trial, 'test-auc')
        history = xgb.cv(param,
                         dtrain,
                         num_boost_round=100,
                         callbacks=[pruning_callback],
                         verbose_eval=False)

        mean_auc = history['test-auc-mean'].values[-1]
        return mean_auc


class Objective_lightgbm_accuracy(object):
    def __init__(self, data):
        self.data = data

    def __call__(self, trial):

        x, y = self.data['x'], self.data['y']
        imbalance_ratio = 1
        if 'imbalance_ratio' in self.data:
            imbalance_ratio = self.data['imbalance_ratio']

        train_x, test_x, train_y, test_y = train_test_split(x,
                                                            y,
                                                            test_size=0.25)

        dtrain = lgb.Dataset(train_x, label=train_y.ravel())
        dtest = lgb.Dataset(test_x, label=test_y.ravel())

        param = {
            'objective': 'binary',
            'metric': 'auc',
            'verbosity': -1,
            'boosting_type': 'gbdt',
            'lambda_l1': trial.suggest_loguniform('lambda_l1', 1e-8, 10.0),
            'lambda_l2': trial.suggest_loguniform('lambda_l2', 1e-8, 10.0),
            'num_leaves': trial.suggest_int('num_leaves', 2, 256),
            'feature_fraction': trial.suggest_uniform('feature_fraction', 0.4,
                                                      1.0),
            'bagging_fraction': trial.suggest_uniform('bagging_fraction', 0.4,
                                                      1.0),
            'bagging_freq': trial.suggest_int('bagging_freq', 1, 7),
            'min_child_samples': trial.suggest_int('min_child_samples', 5,
                                                   100),
            'scale_pos_weight': imbalance_ratio
        }

        # Add a callback for pruning.
        pruning_callback = optuna.integration.LightGBMPruningCallback(
            trial, 'auc')
        gbm = lgb.train(param,
                        dtrain,
                        valid_sets=[dtest],
                        verbose_eval=False,
                        callbacks=[pruning_callback])

        preds = gbm.predict(test_x)
        pred_labels = np.rint(preds)
        accuracy = sklearn.metrics.accuracy_score(test_y, pred_labels)

        return accuracy


class Objective_catboost_accuracy(object):
    def __init__(self, data):
        self.data = data

    def __call__(self, trial):

        x, y = self.data['x'], self.data['y']
        imbalance_ratio = 1
        if 'imbalance_ratio' in self.data:
            imbalance_ratio = self.data['imbalance_ratio']

        train_x, test_x, train_y, test_y = train_test_split(x,
                                                            y,
                                                            test_size=0.25)

        param = {
            'objective':
            trial.suggest_categorical('objective', ['Logloss']),
            'colsample_bylevel':
            trial.suggest_uniform('colsample_bylevel', 0.01, 0.1),
            'depth':
            trial.suggest_int('depth', 1, 12),
            'boosting_type':
            trial.suggest_categorical('boosting_type', ['Ordered', 'Plain']),
            'bootstrap_type':
            trial.suggest_categorical('bootstrap_type',
                                      ['Bayesian', 'Bernoulli', 'MVS']),
            'scale_pos_weight':
            imbalance_ratio,
            # "used_ram_limit":
            # "5gb",
        }

        if param['bootstrap_type'] == 'Bayesian':
            param['bagging_temperature'] = trial.suggest_uniform(
                'bagging_temperature', 0, 10)
        elif param['bootstrap_type'] == 'Bernoulli':
            param['subsample'] = trial.suggest_uniform('subsample', 0.1, 1)

        gbm = cb.CatBoostClassifier(**param, silent=True)

        gbm.fit(train_x,
                train_y,
                eval_set=[(test_x, test_y)],
                silent=True,
                early_stopping_rounds=100)

        preds = gbm.predict(test_x)
        pred_labels = np.rint(preds)
        accuracy = sklearn.metrics.accuracy_score(test_y, pred_labels)

        return accuracy


class Objective_RF_accuracy(object):
    def __init__(self, data):
        self.data = data

    def __call__(self, trial):

        x, y = self.data['x'], self.data['y']

        class_weight = None
        if 'imbalance_ratio' in self.data and self.data['imbalance_ratio'] != 1:
            class_weight = "balanced"

        max_depth = int(trial.suggest_loguniform('max_depth', 2, 32))
        rf_classifier = RandomForestClassifier(max_depth=max_depth,
                                               n_estimators=100,
                                               class_weight=class_weight)

        score = sklearn.model_selection.cross_val_score(rf_classifier,
                                                        x,
                                                        y.ravel(),
                                                        cv=10)
        accuracy = score.mean()

        return accuracy


class Objective_SVC_accuracy(object):
    def __init__(self, data):
        self.data = data

    def __call__(self, trial):

        x, y = self.data['x'], self.data['y']

        class_weight = None
        if 'imbalance_ratio' in self.data and self.data['imbalance_ratio'] != 1:
            class_weight = "balanced"

        svc_c = trial.suggest_loguniform("C", 1e-10, 1e10)
        svc_classifier = sklearn.svm.SVC(C=svc_c,
                                         gamma="auto",
                                         class_weight=class_weight)

        score = sklearn.model_selection.cross_val_score(svc_classifier,
                                                        x,
                                                        y.ravel(),
                                                        cv=10)
        accuracy = score.mean()

        return accuracy
