"""Joint two-marker continuous-time HMM world model.

Pipeline: build_sequences -> build_model -> sample (nutpie) -> export.
See docs/methodology.md sect. 4 for the math.
"""

from .linalg import expm_np, expm_pt, birthdeath_Q  # noqa: F401
from .model import build_sequences, build_model      # noqa: F401
from .simulate import simulate_panel                 # noqa: F401
