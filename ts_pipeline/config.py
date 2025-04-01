import shutil
import warnings

# Tools to check (name: default fallback)
_REQUIRED_TOOLS = {
    "xtb": "xtb",
    "orca": "orca.run",
    "pysis": "pysis",
    "crest": "crest3"
}

# Build resolved paths or fallback
TOOLS = {}

for name, fallback in _REQUIRED_TOOLS.items():
    path = shutil.which(name)
    if path:
        TOOLS[name] = path
    else:
        TOOLS[name] = fallback
        warnings.warn(
            f"[ts_pipeline] Could not find '{name}' in PATH. "
            f"Using fallback: '{fallback}' — check manually if it's correct."
        )

# Convenience access for imports
xtb_path = TOOLS["xtb"]
orca_path = TOOLS["orca"]
pysis_path = TOOLS["pysis"]
crest_path = TOOLS["crest"]