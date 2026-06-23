"""Tests for the AR-emergence -> SFT -> calibration pipeline."""

import numpy as np

from sftmap import calibration as cal
from sftmap.ar_emergence import EmergenceParameters, generate_cycle, mean_tilt_coefficient
from sftmap.sft import SFTModel, cycle_delta_dm


def test_quenching_laws_monotonic():
    p = EmergenceParameters()
    # tilt coefficient decreases with cycle strength (tilt quenching)
    assert mean_tilt_coefficient(80, p) > mean_tilt_coefficient(300, p)
    # mean latitude increases with cycle strength (latitude quenching)
    from sftmap.ar_emergence import mean_latitude
    assert mean_latitude(80, p) < mean_latitude(300, p)


def test_more_regions_in_stronger_cycles():
    rng = np.random.default_rng(0)
    weak = np.mean([len(generate_cycle(80, rng=rng)) for _ in range(50)])
    strong = np.mean([len(generate_cycle(300, rng=rng)) for _ in range(50)])
    assert strong > weak


def test_delta_dm_increases_then_saturates():
    # The surrogate dDM(SN) should rise then flatten (increase-then-saturate).
    rng = np.random.default_rng(0)

    def mean_dm(sn):
        return np.mean([cycle_delta_dm(generate_cycle(sn, rng=rng)) for _ in range(200)])

    low, mid, high = mean_dm(60), mean_dm(180), mean_dm(320)
    assert low < mid                      # rises at first
    assert abs(high - mid) < 0.5 * mid    # saturates (small change at high SN)


def test_calibration_recovers_standard_set():
    # The default emergence model is tuned so the calibration recovers the
    # J20 / Wang standard set (k1~6.94, quench~75.85, stoch~0.17).
    res = cal.calibrate(n_realizations=400, seed=0)
    assert 5.5 < res.k1 < 8.5
    assert 55 < res.quench < 95
    assert 0.12 < res.stoch < 0.24


def test_full_pde_dipole_is_positive_and_coherent():
    # The full PDE solver should give a positive net dDM for a Joy's-law cycle
    # (hemispheres add coherently).
    rng = np.random.default_rng(0)
    ars = generate_cycle(180.0, rng=rng)
    model = SFTModel(n_grid=181, eta=250.0)
    assert model.cycle_delta_dm(ars) > 0.0
