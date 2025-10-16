import sys
import numpy as np
import pytest
from ptg_model.core_functions import rate_adj, defaultparameters, release_rate

# --- Mock dependencies from ptg_model submodules ---
from types import SimpleNamespace

mock_utils = SimpleNamespace(
    stim=lambda x, mode=None: np.tanh(x),
    sens=lambda a, b: 1 / (1 + np.exp(-(a + b))),
    smooth_pw=lambda x, endpoints=None: np.clip(x, 0, 1),
)

mock_core = SimpleNamespace(
    rate_adj=lambda x, param=None: 1 + 0.1 * np.tanh(x),
    defaultparameters=lambda mode=None: 1.0,
    release_rate=lambda c, rp: 0.05 * rp * c,
)

sys.modules["ptg_model.utils"] = mock_utils
sys.modules["ptg_model.core_functions"] = mock_core

# --- Import the model after mocks are registered ---
from ptg_model.model import deriv


@pytest.fixture
def base_inputs():
    """Base input fixture for reproducible ODE tests."""
    y = np.ones(23) * 0.5
    t = 1.0
    params = {
    "endpoints_p": [0, 1],
    "endpoints_d": [0, 1],
    "copt": 1.2,
    "dopt": 0.9,
    "popt": 1.1,
    "c_pat": 1.0,
    "p_pat": 1.0,
    "d_pat": 1.0,
    "s0": 1.0,
    "tm": 10.0,
    "gfr_in": 1.0,
    "y_pat": np.ones(23),
    "calcium_clamp": True,
    }

    return t, y, params


def test_deriv_output_shape(base_inputs):
    """Ensure deriv() returns the correct dimensionality."""
    t, y, params = base_inputs
    dydt = deriv(t, y, **params)
    assert isinstance(dydt, np.ndarray)
    assert dydt.shape == (23,)


def test_deriv_finite_values(base_inputs):
    """Ensure no NaNs or infinities appear in outputs."""
    t, y, params = base_inputs
    dydt = deriv(t, y, **params)
    assert np.all(np.isfinite(dydt)), "deriv() produced NaN or inf values."


def test_deriv_zero_threshold_behavior(base_inputs):
    """Verify that very small derivatives are set to zero."""
    t, y, params = base_inputs
    y *= 1e-8  # nearly zero input
    dydt = deriv(t, y, **params)
    assert np.all(
        (np.abs(dydt) == 0) | (np.abs(dydt) > 1e-12)
    ), "Small values not zeroed as expected."


def test_deriv_consistency(base_inputs):
    """Check small time perturbation produces similar derivatives."""
    t, y, params = base_inputs
    dydt1 = deriv(t, y, **params)
    dydt2 = deriv(t + 0.001, y, **params)
    diff = np.linalg.norm(dydt1 - dydt2)
    assert diff < 1.0, f"Derivatives changed too abruptly (diff={diff:.3f})"
