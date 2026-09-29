"""Conservative, descriptive cycle-to-cycle drift screens."""

from __future__ import annotations

import math
import statistics


def cycle_drift(values: list[float], transient_cycles: int,
                relative_tolerance: float = .1) -> dict:
    """Compare early/late cycle means using contiguous two-cycle batches.

    At least 16 eligible cycles provide four batches in each half. Two batch
    standard errors form a descriptive margin, not a confidence interval:
    even the batches can remain correlated. A bounded mean shift does not
    certify a stationary distribution.
    """
    if transient_cycles < 0 or not 0 < relative_tolerance < 1:
        raise ValueError("invalid drift-screen settings")
    eligible=values[transient_cycles:]
    result={"eligible_cycles":len(eligible),"required_cycles":16,
            "relative_tolerance":relative_tolerance,
            "status":"insufficient_cycles","early_mean":None,
            "late_mean":None,"mean_shift":None,
            "relative_mean_shift":None,"screening_standard_error":None,
            "two_se_margin":None,"absolute_tolerance":None,
            "early_batch_means":[],"late_batch_means":[]}
    if len(eligible)<16:
        return result
    midpoint=len(eligible)//2
    early,late=eligible[:midpoint],eligible[midpoint:]

    def batch_means(half):
        return [statistics.mean(half[start:end]) for start,end in
                ((i*len(half)//4,(i+1)*len(half)//4) for i in range(4))]

    early_batches=batch_means(early)
    late_batches=batch_means(late)
    early_mean=statistics.mean(early)
    late_mean=statistics.mean(late)
    shift=late_mean-early_mean
    error=math.sqrt(statistics.variance(early_batches)/4+
                    statistics.variance(late_batches)/4)
    margin=2*error
    scale=max(1.,statistics.mean(abs(value) for value in eligible))
    tolerance=relative_tolerance*scale
    if abs(shift)-margin>tolerance:
        status="drift_signal"
    elif abs(shift)+margin<tolerance:
        status="shift_bounded"
    else:
        status="inconclusive"
    result.update(status=status,early_mean=early_mean,late_mean=late_mean,
                  mean_shift=shift,relative_mean_shift=shift/scale,
                  screening_standard_error=error,two_se_margin=margin,
                  absolute_tolerance=tolerance,
                  early_batch_means=early_batches,
                  late_batch_means=late_batches)
    return result


def gas_energy_drift(cycle_end_energies: list[float],
                     transient_cycles: int) -> dict:
    """Compatibility wrapper for the cycle-end gas-energy drift screen."""
    return cycle_drift(cycle_end_energies,transient_cycles)
