import numpy as np
import pytest
from ptg_model.utils import smooth_pw, stim, sens

# --- Fixtures ----------------------------------------------------------------


@pytest.fixture
def simple_endpoints():
    """Provide simple piecewise-linear endpoints for smooth_pw testing."""
    # 3 segments: (0,0)→(1,1)→(2,0)
    endpoints = np.array([[0, 1, 2], [0, 1, 0]])
    return endpoints


# --- smooth_pw tests ----------------------------------------------------------


def test_smooth_pw_scalar_output(simple_endpoints):
    """smooth_pw should return a scalar for scalar input."""
    x = 0.5
    val = smooth_pw(x, simple_endpoints)
    assert np.isscalar(val) or np.size(val) == 1
    assert np.isfinite(val)


def test_smooth_pw_vectorized_behavior(simple_endpoints):
    """smooth_pw should produce continuous, finite values across range."""
    x_vals = np.linspace(0, 2, 20)
    y_vals = [smooth_pw(x, simple_endpoints) for x in x_vals]
    assert np.all(np.isfinite(y_vals))
    assert np.ptp(y_vals) > 0, "Output range should not be zero"


# --- stim tests ---------------------------------------------------------------


@pytest.mark.parametrize("param", ["c", "p", "d", "unknown"])
def test_stim_behavior(param):
    """stim should produce finite output for all parameter types."""
    x_vals = np.linspace(-5, 5, 100)
    out = stim(x_vals, param)
    assert out.shape == x_vals.shape
    assert np.all(np.isfinite(out))

    # Calcium and phosphate stimuli should be nonzero around midrange
    if param in ["c", "p"]:
        assert np.any(np.abs(out) > 0), f"stim({param}) returned all zeros unexpectedly"


def test_stim_cutoff_behavior():
    """stim should apply cutoff correctly near zero."""
    near_zero = stim(0.0, "c")
    assert np.all(
        np.abs(near_zero) <= 1e-3
    ), "stim did not cut off near zero as expected"


# --- sens tests ---------------------------------------------------------------


def test_sens_scalar_output():
    """sens should return scalar for scalar input."""
    val = sens(1.0, 1.0)
    assert np.isscalar(val)
    assert np.isfinite(val)


def test_sens_array_output():
    """sens should work elementwise and stay normalized near 1."""
    c_vals = np.linspace(0.5, 2.0, 10)
    d_vals = np.linspace(0.5, 2.0, 10)
    out = sens(c_vals, d_vals)
    assert out.shape == c_vals.shape
    assert np.all(np.isfinite(out))
    assert np.all(out > 0)
    assert np.isclose(
        out[5], 1.0, atol=0.2
    ), "sens not approximately normalized at midrange"


def test_sens_monotonic_trend():
    """sens should be non-decreasing overall with calcium/calcitriol."""
    c_vals = np.linspace(0, 10, 50)
    d_vals = np.linspace(0, 10, 50)
    out = sens(c_vals, d_vals)
    diffs = np.diff(out)
    # allow small numerical noise, but no negative trend
    assert np.all(
        diffs >= -1e-6
    ), "sens should not decrease with higher calcium/calcitriol"
