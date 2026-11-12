"""Numerical conditioning of periodic LSO cascades."""

import numpy as np

from .lso_filter import _filter_tf, design_lso_filter, order_periodic_passes


def test_periodic_pass_order_preserves_response_and_limits_intermediate_gain():
    coefficients = design_lso_filter(32 / np.sqrt(2), 19)
    xi = np.linspace(0, np.pi, 4097)
    original = np.array([_filter_tf(np.asarray(c), xi) for c in coefficients])
    ordered = np.array([_filter_tf(np.asarray(c), xi) for c in order_periodic_passes(coefficients)])
    assert np.max(np.abs(np.cumprod(original, axis=0))) > 1e6
    assert np.max(np.abs(np.cumprod(ordered, axis=0))) < 1.01
    np.testing.assert_allclose(np.prod(ordered, axis=0), np.prod(original, axis=0), rtol=1e-12, atol=1e-12)
