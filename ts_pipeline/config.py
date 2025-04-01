import shutil
import warnings

# Tools to check (name: default fallback)
_REQUIRED_TOOLS = {
    "xtb"  : "",
    "orca" : "",
    "pysis": "",
    "crest": ""
}

# Build resolved paths or fallback
TOOLS = {}

for name, user_defined in _REQUIRED_TOOLS.items():
    path = shutil.which(name)
    if path:
        TOOLS[name] = path
    else:
        TOOLS[name] = user_defined
        warnings.warn(
            f"[ts_pipeline] Could not find '{name}' in PATH. "
            f"Using user-defined value f***REMOVED*** ts_pipeline/config.py: '{user_defined}' — make sure it's correct."
        )

# Convenience access for imports
xtb_path = TOOLS["xtb"]
orca_path = TOOLS["orca"]
pysis_path = TOOLS["pysis"]
crest_path = TOOLS["crest"]