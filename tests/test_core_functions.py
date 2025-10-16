import numpy as np
import pytest
from ptg_model.core_functions import rate_adj, defaultparameters, release_rate


# --- rate_adj tests ---------------------------------------------------------


def test_rate_adj_scalar_behavior():
    """Check that rate_adj behaves correctly for scalar inputs."""
    r, a = 10, 0.2
    # Below threshold
    val_low = rate_adj(0.5, [r, a])
    # Above threshold
    val_high = rate_adj(1.5, [r, a])

    assert np.isclose(val_low, (r - a * r) * 0.5 + a * r)
    assert np.isclose(val_high, r)
    assert val_low < val_high, "rate_adj should increase with c"


def test_rate_adj_vectorized_behavior():
    """Ensure rate_adj works on arrays and matches elementwise behavior."""
    r, a = 8, 0.1
    c_values = np.array([0.2, 0.8, 1.2])
    out = rate_adj(c_values, [r, a])
    expected = np.where(c_values < 1, (r - a * r) * c_values + a * r, r)
    np.testing.assert_allclose(out, expected)


# --- defaultparameters tests ------------------------------------------------


@pytest.mark.parametrize(
    "name,expected_rate_range",
    [
        ("degrad", (0.01, 1.0)),
        ("prolif", (1.0, 5.0)),
        ("prod", (300, 5000)),
        ("clear", (10, 100)),
    ],
)
def test_defaultparameters_valid(name, expected_rate_range):
    """Ensure all defaultparameter presets return valid positive values."""
    params = defaultparameters(name)
    assert isinstance(params, list) and len(params) == 2
    rate, adj = params
    assert expected_rate_range[0] <= rate <= expected_rate_range[1]
    assert adj > 0


def test_defaultparameters_invalid():
    """Ensure an invalid name raises ValueError."""
    with pytest.raises(ValueError):
        defaultparameters("unknown_function")


# --- release_rate tests -----------------------------------------------------


def test_release_rate_monotonicity():
    """Check that release_rate shows a globally decreasing trend with calcium concentration."""
    rp = 1.0
    c_values = np.linspace(0.5, 2.0, 50)
    rates = np.array([release_rate(c, rp) for c in c_values])

    # Allow for small numerical plateaus or jitters (<1e-6)
    diffs = np.diff(rates)
    assert np.all(
        diffs <= 1e-6
    ), f"release_rate should not increase with calcium, got max positive diff {np.max(diffs):.2e}"

    # Ensure the mean of first 5 values is greater than last 5 → global decreasing trend
    assert np.mean(rates[:5]) > np.mean(
        rates[-5:]
    ), "release_rate should show an overall decrease with calcium"


def test_release_rate_scaling_with_rp():
    """Ensure phosphate scaling (rp) influences amplitude as expected."""
    c = 1.0
    val_low = release_rate(c, rp=0.5)
    val_high = release_rate(c, rp=2.0)
    assert val_low < val_high, "Higher rp should increase release rate amplitude"


def test_release_rate_expected_range():
    """Check that release_rate values are finite and within plausible range."""
    c_values = np.linspace(0.5, 2.0, 20)
    rates = release_rate(c_values, rp=1.0)
    assert np.all(np.isfinite(rates)), "release_rate produced NaN or inf"
    assert np.all((rates > 0) & (rates < 20)), "release_rate out of expected range"
