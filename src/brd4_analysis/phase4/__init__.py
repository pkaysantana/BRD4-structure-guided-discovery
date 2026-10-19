"""Phase 4 adversarial validation of the frozen Phase 3 scoring harness.

This package is deliberately isolated: nothing in ``brd4_analysis`` or
``brd4_analysis.phase2`` imports it, and it imports the frozen Phase 2/3
metrics read-only. It therefore cannot alter the analysis-source hashes
recorded in preregistered experiment manifests, and it never changes the
registered acceptance thresholds.
"""
