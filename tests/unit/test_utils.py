#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (c) 2025 The WfCommons Team.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import numpy as np
import pytest

import wfcommons.utils


class TestUtils:

    @pytest.mark.unit
    def test_best_fit_distribution_uses_normalized_samples_and_density_histogram(
        self,
        monkeypatch: pytest.MonkeyPatch,
    ) -> None:
        # Twenty observations produce two bins under the current
        # ceil(len(data) / 10) rule.
        data = list(range(20))

        expected_normalized = np.asarray(data, dtype=float)
        expected_normalized = (
            (expected_normalized - expected_normalized.min())
            / (expected_normalized.max() - expected_normalized.min())
        )

        histogram_call = {}
        fit_calls = []
        pdf_calls = []

        def fake_histogram(values, bins, density=False):
            histogram_call["values"] = np.asarray(
                values,
                dtype=float,
            ).copy()
            histogram_call["bins"] = bins
            histogram_call["density"] = density

            # A valid two-bin density:
            # 0.5 * 0.5 + 1.5 * 0.5 == 1.0
            return (
                np.array([0.5, 1.5], dtype=float),
                np.array([0.0, 0.5, 1.0], dtype=float),
            )

        class FakeDistribution:

            def __init__(self, index: int) -> None:
                self.index = index

            def fit(self, values):
                fit_calls.append(
                    np.asarray(values, dtype=float).copy()
                )
                return 0.0, 1.0

            def pdf(self, x, *args, loc, scale):
                pdf_calls.append(
                    np.asarray(x, dtype=float).copy()
                )

                if self.index == 0:
                    # Closest to the density histogram [0.5, 1.5].
                    return np.array([0.55, 1.45], dtype=float)

                if self.index == 1:
                    # Closest to the old min-max-scaled histogram
                    # [0.0, 1.0]. This candidate would win if that
                    # obsolete transformation were restored.
                    return np.array([0.05, 0.95], dtype=float)

                return np.array([10.0, 10.0], dtype=float)

        class FakeStats:

            def __init__(self) -> None:
                self.requested_names = []

            def __getattr__(self, name):
                index = len(self.requested_names)
                self.requested_names.append(name)
                return FakeDistribution(index)

        fake_stats = FakeStats()

        monkeypatch.setattr(
            wfcommons.utils.np,
            "histogram",
            fake_histogram,
        )
        monkeypatch.setattr(
            wfcommons.utils.scipy,
            "stats",
            fake_stats,
        )

        distribution_name, params = (
            wfcommons.utils.best_fit_distribution(data)
        )

        assert histogram_call["bins"] == 2
        assert histogram_call["density"] is True
        np.testing.assert_allclose(
            histogram_call["values"],
            expected_normalized,
        )

        # Every candidate must be fitted to the normalized observations,
        # not to the histogram heights.
        assert fit_calls
        for fit_data in fit_calls:
            np.testing.assert_allclose(
                fit_data,
                expected_normalized,
            )

        # The PDF must be evaluated at the centers of [0.0, 0.5]
        # and [0.5, 1.0].
        expected_centers = np.array([0.25, 0.75], dtype=float)

        assert pdf_calls
        for pdf_x in pdf_calls:
            np.testing.assert_allclose(
                pdf_x,
                expected_centers,
            )

        # The first fake distribution is closest to the density
        # histogram. The second would win under the old transformation.
        assert distribution_name == fake_stats.requested_names[0]
        assert params == (0.0, 1.0)

    @pytest.mark.unit
    def test_best_fit_distribution_returns_usable_fit(self) -> None:
        rng = np.random.default_rng(12345)
        data = rng.lognormal(
            mean=0.0,
            sigma=0.5,
            size=30,
        ).tolist()

        distribution_name, params = (
            wfcommons.utils.best_fit_distribution(data)
        )

        assert isinstance(distribution_name, str)
        assert distribution_name
        assert hasattr(
            wfcommons.utils.scipy.stats,
            distribution_name,
        )

        params_array = np.asarray(params, dtype=float)

        # SciPy continuous-distribution fits end with loc and scale.
        assert params_array.size >= 2
        assert np.all(np.isfinite(params_array))
        assert params_array[-1] > 0

    @pytest.mark.unit
    def test_generate_rvs_without_distribution_returns_minimum(
        self,
    ) -> None:
        assert wfcommons.utils.generate_rvs(
            None,
            min_value=10,
            max_value=100,
        ) == 10

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "n,r,combinations",
        [
            pytest.param(2, 2, 1),
            pytest.param(10, 2, 45),
            pytest.param(10, 3, 120),
        ],
    )
    def test_ncr(
        self,
        n: int,
        r: int,
        combinations: int,
    ) -> None:
        assert wfcommons.utils.ncr(n, r) == combinations
