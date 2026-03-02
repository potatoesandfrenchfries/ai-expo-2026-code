"""Metrics and logging for the adversarial feedback loop (Step 5).

Tracks the following signals over a rolling window of the last N molecules:
  • Per-constraint violation rate (fraction of molecules that violated each type)
  • Generator reward curve (mean reward per training batch)
  • Surrogate discriminator accuracy (vs. ground-truth Rocq verdicts)
  • Rocq proof latency per molecule (seconds)
  • Fraction of molecules that pass all constraints over training time

All metrics are stored in Python deques (no external dependencies).
Optional structured output to a JSONL log file is supported.

Usage
─────
    tracker = ConstraintMetricsTracker(window=1000)

    # After each batch:
    tracker.update_violation_rates(violation_rates_dict)
    tracker.update_reward(mean_batch_reward)
    tracker.update_rocq_latency(latency_seconds)
    tracker.update_pass_rate(n_passed, n_total)

    # Optionally, after surrogate retrain:
    tracker.update_surrogate_accuracy(accuracy)

    # Print a summary:
    tracker.print_summary()

    # Get a dict (for external logging):
    d = tracker.summary_dict()
"""

import json
import os
import time
from collections import deque
from typing import Optional

from .constraint_vector import CONSTRAINT_TYPES


# ---------------------------------------------------------------------------
# ConstraintMetricsTracker
# ---------------------------------------------------------------------------

class ConstraintMetricsTracker:
    """Rolling-window metrics tracker for the adversarial feedback loop.

    Parameters
    ----------
    window : int
        Number of most-recent per-molecule observations to keep.
    log_file : str, optional
        If given, each call to log_to_file() appends a JSON record here.
    """

    def __init__(self, window: int = 1000, log_file: Optional[str] = None):
        self._window = window

        # Per-constraint violation history: deque of 0/1 per molecule
        self._violations: dict = {
            c: deque(maxlen=window) for c in CONSTRAINT_TYPES
        }

        # Reward per batch (not per-molecule, so no fixed window)
        self._rewards: deque = deque(maxlen=window)

        # Surrogate discriminator accuracy history
        self._surrogate_acc: deque = deque(maxlen=100)

        # Rocq proof latency (seconds per molecule)
        self._rocq_latency: deque = deque(maxlen=window)

        # All-constraints pass rate: store (n_passed, n_total) tuples
        self._pass_counts: deque = deque(maxlen=window)

        # Training step counter
        self._step = 0

        # Optional JSONL log
        self._log_file = log_file
        if log_file:
            os.makedirs(os.path.dirname(log_file) or ".", exist_ok=True)

    # ------------------------------------------------------------------
    # Update methods
    # ------------------------------------------------------------------

    def update_violation_rates(self, rates: dict) -> None:
        """Record per-constraint violation rates from a batch.

        Parameters
        ----------
        rates : dict[str, float]
            Fraction of molecules in the batch that violated each constraint.
            Keys should match CONSTRAINT_TYPES.
        """
        for ctype in CONSTRAINT_TYPES:
            if ctype in rates:
                # Convert rate to a 0/1 flag per-batch for rolling average
                self._violations[ctype].append(float(rates[ctype]))

    def update_violation_flags(self, flags: dict) -> None:
        """Record per-molecule binary violation flags.

        Parameters
        ----------
        flags : dict[str, bool]
            Single-molecule violation result (from ConstraintViolationVector).
        """
        for ctype in CONSTRAINT_TYPES:
            self._violations[ctype].append(float(bool(flags.get(ctype, False))))

    def update_reward(self, mean_reward: float) -> None:
        """Record the mean reward for the most recent batch."""
        self._rewards.append(mean_reward)
        self._step += 1

    def update_surrogate_accuracy(self, accuracy: float) -> None:
        """Record surrogate discriminator accuracy after a retrain."""
        self._surrogate_acc.append(accuracy)

    def update_rocq_latency(self, latency_seconds: float) -> None:
        """Record per-molecule Rocq proof latency."""
        self._rocq_latency.append(latency_seconds)

    def update_pass_rate(self, n_passed: int, n_total: int) -> None:
        """Record how many molecules in the batch passed all constraints."""
        if n_total > 0:
            self._pass_counts.append(n_passed / n_total)

    # ------------------------------------------------------------------
    # Summary computation
    # ------------------------------------------------------------------

    def violation_rate(self, ctype: str) -> float:
        """Rolling-window violation rate for a single constraint type."""
        buf = self._violations.get(ctype)
        if not buf:
            return float("nan")
        return sum(buf) / len(buf)

    def mean_reward(self) -> float:
        """Mean reward over the rolling window."""
        if not self._rewards:
            return float("nan")
        return sum(self._rewards) / len(self._rewards)

    def mean_surrogate_accuracy(self) -> float:
        """Mean surrogate accuracy over the last 100 retrains."""
        if not self._surrogate_acc:
            return float("nan")
        return sum(self._surrogate_acc) / len(self._surrogate_acc)

    def mean_rocq_latency(self) -> float:
        """Mean Rocq proof latency (seconds) over the rolling window."""
        if not self._rocq_latency:
            return float("nan")
        return sum(self._rocq_latency) / len(self._rocq_latency)

    def pass_rate(self) -> float:
        """Fraction of molecules that passed all constraints (rolling window)."""
        if not self._pass_counts:
            return float("nan")
        return sum(self._pass_counts) / len(self._pass_counts)

    def summary_dict(self) -> dict:
        """Return all current metrics as a flat dict suitable for logging."""
        d = {
            "step":               self._step,
            "mean_reward":        round(self.mean_reward(), 5),
            "pass_rate":          round(self.pass_rate(), 4),
            "surrogate_accuracy": round(self.mean_surrogate_accuracy(), 4),
            "rocq_latency_s":     round(self.mean_rocq_latency(), 4),
        }
        for ctype in CONSTRAINT_TYPES:
            d[f"violation_rate_{ctype}"] = round(self.violation_rate(ctype), 4)
        return d

    def print_summary(self) -> None:
        """Print a formatted metrics summary to stdout."""
        d = self.summary_dict()
        print(
            f"  [Metrics step={d['step']}]"
            f"  reward={d['mean_reward']:.4f}"
            f"  pass_rate={d['pass_rate']:.3f}"
        )
        viol_parts = []
        for ctype in CONSTRAINT_TYPES:
            rate = d[f"violation_rate_{ctype}"]
            if not (rate != rate):  # skip NaN
                viol_parts.append(f"{ctype}={rate:.3f}")
        if viol_parts:
            print(f"  Violation rates: {' '.join(viol_parts)}")
        if not (d["surrogate_accuracy"] != d["surrogate_accuracy"]):
            print(f"  Surrogate acc: {d['surrogate_accuracy']:.4f}")
        if not (d["rocq_latency_s"] != d["rocq_latency_s"]):
            print(f"  Rocq latency:  {d['rocq_latency_s']:.3f}s/mol")

    def log_to_file(self) -> None:
        """Append a JSON record of current metrics to the log file."""
        if self._log_file is None:
            return
        record = {"timestamp": time.time(), **self.summary_dict()}
        with open(self._log_file, "a") as f:
            f.write(json.dumps(record) + "\n")


# ---------------------------------------------------------------------------
# Timing context manager for Rocq latency measurement
# ---------------------------------------------------------------------------

class RocqTimer:
    """Context manager that records elapsed time into a ConstraintMetricsTracker.

    Usage::
        with RocqTimer(tracker):
            result = check_molecule(smiles)
    """

    def __init__(self, tracker: ConstraintMetricsTracker):
        self._tracker = tracker
        self._start: float = 0.0

    def __enter__(self):
        self._start = time.perf_counter()
        return self

    def __exit__(self, *_):
        elapsed = time.perf_counter() - self._start
        self._tracker.update_rocq_latency(elapsed)
