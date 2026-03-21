from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.model_selection import KFold, LeaveOneOut


@dataclass(slots=True)
class ShufflingTestResult:
    randomized_mse: np.ndarray
    real_mse: float


class MetaboLiteLearner:
    """Partial least squares learner mirroring the MATLAB class."""

    def __init__(
        self,
        x: np.ndarray,
        y: np.ndarray,
        kfold: int = 0,
        max_components: int = 30,
        nrandomized: int = 1000,
        random_state: int = 0,
    ) -> None:
        self.x_full_data = np.asarray(x, dtype=float)
        self.y_full_data = np.asarray(y, dtype=float)
        self.maxn = min(max_components, self.x_full_data.shape[0] - 1, self.x_full_data.shape[1])
        self.nrandomized = nrandomized
        self.random_state = random_state
        self.cv = LeaveOneOut() if kfold == 0 else KFold(n_splits=kfold, shuffle=True, random_state=random_state)
        self.kfold = self.x_full_data.shape[0] if kfold == 0 else kfold

        self.nopt, self.test_sse, self.train_sse = self.optimize_components_and_learn()
        (
            self.model,
            self.beta,
            self.y_pred,
            self.loss,
            self.sse,
            self.x_loadings,
            self.y_loadings,
            self.x_scores,
            self.y_scores,
            self.pctvar,
        ) = self.learn(self.x_full_data, self.y_full_data, self.nopt)

    def map_to_latent_space(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        x = np.asarray(x, dtype=float)
        latent_scores = self.model.transform(x)
        predicted_y = self.model.predict(x)
        return latent_scores, predicted_y

    def shuffling_test(self) -> ShufflingTestResult:
        rng = np.random.default_rng(self.random_state)
        randomized = np.zeros(self.nrandomized, dtype=float)
        for idx in range(self.nrandomized):
            y_rand = self.y_full_data[rng.permutation(self.y_full_data.shape[0]), :]
            _, _, kfold_sse, _ = self.cross_validation_evaluation(self.x_full_data, y_rand, self.nopt)
            randomized[idx] = np.mean(kfold_sse)
        return ShufflingTestResult(randomized_mse=randomized, real_mse=float(np.mean(self.test_sse[:, self.nopt - 1])))

    def learn(
        self,
        x: np.ndarray,
        y: np.ndarray,
        n_components: int,
    ) -> tuple[PLSRegression, np.ndarray, np.ndarray, np.ndarray, float, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        model = PLSRegression(n_components=n_components, scale=False)
        model.fit(x, y)
        y_pred = model.predict(x)
        loss = (y_pred - y) ** 2
        sse = float(loss.sum())

        beta = np.vstack([model.intercept_.reshape(1, -1), model.coef_])
        pctvar = self._calculate_pctvar(x, y, model)
        return (
            model,
            beta,
            y_pred,
            loss,
            sse,
            model.x_loadings_,
            model.y_loadings_,
            model.x_scores_,
            model.y_scores_,
            pctvar,
        )

    def cross_validation_evaluation(
        self,
        x: np.ndarray,
        y: np.ndarray,
        n_components: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        y_pred = np.zeros_like(y, dtype=float)
        train_sse = np.zeros(self.kfold, dtype=float)
        kfold_sse = np.zeros(self.kfold, dtype=float)
        beta_distribution = np.zeros((x.shape[1] + 1, y.shape[1], x.shape[0]), dtype=float)

        for fold_index, (train_idx, test_idx) in enumerate(self.cv.split(x, y)):
            model, beta, _, _, sse, *_ = self.learn(x[train_idx], y[train_idx], n_components)
            predicted = model.predict(x[test_idx])
            y_pred[test_idx] = predicted
            train_sse[fold_index] = sse / train_idx.size
            kfold_sse[fold_index] = np.sum((predicted - y[test_idx]) ** 2) / test_idx.size
            for row_idx in test_idx:
                beta_distribution[:, :, row_idx] = beta

        return y_pred, train_sse, kfold_sse, beta_distribution

    def optimize_components_and_learn(self) -> tuple[int, np.ndarray, np.ndarray]:
        train_sse = np.zeros((self.kfold, self.maxn), dtype=float)
        test_sse = np.zeros((self.kfold, self.maxn), dtype=float)

        for component_index in range(1, self.maxn + 1):
            _, train_sse_component, test_sse_component, _ = self.cross_validation_evaluation(
                self.x_full_data,
                self.y_full_data,
                component_index,
            )
            train_sse[:, component_index - 1] = train_sse_component
            test_sse[:, component_index - 1] = test_sse_component

        mean_test_sse = test_sse.mean(axis=0)
        n_min = int(np.argmin(mean_test_sse)) + 1
        std_error = test_sse.std(axis=0, ddof=0) / np.sqrt(test_sse.shape[0])
        min_error = mean_test_sse[n_min - 1] + std_error[n_min - 1]
        candidate_idx = np.where((mean_test_sse > min_error) & (np.arange(1, self.maxn + 1) < n_min))[0]
        n_opt = int(candidate_idx.max() + 1) if candidate_idx.size else n_min
        return n_opt, test_sse, train_sse

    @staticmethod
    def _calculate_pctvar(x: np.ndarray, y: np.ndarray, model: PLSRegression) -> np.ndarray:
        total_x = np.sum((x - x.mean(axis=0, keepdims=True)) ** 2)
        total_y = np.sum((y - y.mean(axis=0, keepdims=True)) ** 2)

        x_scores = model.x_scores_
        x_loadings = model.x_loadings_
        y_scores = model.y_scores_
        y_loadings = model.y_loadings_

        explained_x = []
        explained_y = []
        for component_index in range(model.n_components):
            x_hat = np.outer(x_scores[:, component_index], x_loadings[:, component_index])
            y_hat = np.outer(y_scores[:, component_index], y_loadings[:, component_index])
            explained_x.append(np.sum(x_hat**2) / total_x if total_x else 0.0)
            explained_y.append(np.sum(y_hat**2) / total_y if total_y else 0.0)
        return np.vstack([explained_x, explained_y])
