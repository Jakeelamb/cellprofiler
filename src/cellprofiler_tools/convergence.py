"""Per-species convergence tracking via SEM%.

Stops collecting for a species once the standard error of the mean
(as a percentage of the mean) drops below a threshold, meaning
additional measurements won't meaningfully change the estimate.

Used by both the Cellpose cell-size pipeline and the nucleus IOD pipeline.
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np


class ConvergenceTracker:
    """Track per-species measurement convergence via SEM%.

    With typical RBC CV ~15%, convergence at SEM%<3% needs ~25 cells.

    Distinguishes truly converged (SEM% met) from capped (hit max_cells
    without meeting SEM%). Species exempt from the cap get unlimited cells.
    """

    def __init__(self, min_cells: int = 30, max_cells: int = 200, sem_pct: float = 3.0):
        self.min_cells = min_cells
        self.max_cells = max_cells
        self.sem_pct = sem_pct
        self.species_values: dict[str, list[float]] = defaultdict(list)
        self.exempt_species: set[str] = set()

    def add(self, species: str, values: list[float]) -> None:
        self.species_values[species].extend(values)

    def add_cells(self, species: str, measurements: list[dict], key: str = "area_um2") -> None:
        for m in measurements:
            self.species_values[species].append(float(m[key]))

    def exempt(self, species: str) -> None:
        """Remove the cell cap for a species."""
        self.exempt_species.add(species)

    def _sem_pct_value(self, species: str) -> float:
        vals = self.species_values.get(species, [])
        n = len(vals)
        if n < 2:
            return float("inf")
        arr = np.array(vals)
        mean = arr.mean()
        if mean <= 0:
            return float("inf")
        sem = arr.std(ddof=1) / np.sqrt(n)
        return 100 * sem / mean

    def truly_converged(self, species: str) -> bool:
        n = len(self.species_values.get(species, []))
        if n < self.min_cells:
            return False
        return self._sem_pct_value(species) < self.sem_pct

    def is_capped(self, species: str) -> bool:
        if species in self.exempt_species:
            return False
        n = len(self.species_values.get(species, []))
        return n >= self.max_cells and not self.truly_converged(species)

    def is_done(self, species: str) -> bool:
        if self.truly_converged(species):
            return True
        if species in self.exempt_species:
            return False
        n = len(self.species_values.get(species, []))
        return n >= self.max_cells

    def cells_needed(self, species: str) -> int:
        vals = self.species_values.get(species, [])
        n = len(vals)
        cap = None if species in self.exempt_species else self.max_cells
        if cap is not None and n >= cap:
            return 0
        if n < 2:
            return (cap or 10000) - n
        arr = np.array(vals)
        mean = arr.mean()
        if mean <= 0:
            return (cap or 10000) - n
        if self._sem_pct_value(species) < self.sem_pct:
            return 0
        sd = arr.std(ddof=1)
        target_n = int(np.ceil((100 * sd / (mean * self.sem_pct)) ** 2))
        remaining = max(0, target_n - n)
        if cap is not None:
            remaining = min(remaining, cap - n)
        return remaining

    def cell_count(self, species: str) -> int:
        return len(self.species_values.get(species, []))

    def status(self, species: str) -> str:
        vals = self.species_values.get(species, [])
        n = len(vals)
        if n < 2:
            return f"n={n}"
        sem_pct = self._sem_pct_value(species)
        if self.truly_converged(species):
            tag = " CONVERGED"
        elif self.is_capped(species):
            tag = " CAPPED (not converged)"
        elif species in self.exempt_species:
            tag = " EXEMPT (no cap)"
        else:
            tag = ""
        arr = np.array(vals)
        return f"n={n}, mean={arr.mean():.1f}, SEM%={sem_pct:.1f}%{tag}"

    def summary(self) -> str:
        lines = []
        for sp in sorted(self.species_values):
            lines.append(f"  {sp}: {self.status(sp)}")
        return "\n".join(lines)

    def unconverged_species(self) -> list[str]:
        return [sp for sp in self.species_values if self.is_capped(sp)]


def interleave_by_species(jobs: list[dict], species_key: str = "species") -> list[dict]:
    """Round-robin jobs across species for even sampling."""
    species_queues: dict[str, list[dict]] = defaultdict(list)
    for j in jobs:
        species_queues[j[species_key]].append(j)
    interleaved: list[dict] = []
    while any(species_queues.values()):
        for sp in sorted(species_queues):
            if species_queues[sp]:
                interleaved.append(species_queues[sp].pop(0))
        species_queues = {k: v for k, v in species_queues.items() if v}
    return interleaved
