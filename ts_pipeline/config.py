# This is config file for ts_pipeline tool.
# If you have all the tools available in your global PATH, you do not have to change anything here.
# If you use more than one version of any of these tools, we recommend providing valid paths to desired versions.
# If user-defined paths is invalid, the config file will try to find them in your PATH by running "which" with default names (i.e., "orca", "xtb", "pysis" and "crest").
# If user-configured paths are invalid, you will keep seeing user warnings, but you can ignore them.

import shutil, sys, warnings, pprint

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
    path = shutil.which(user_defined)
    if path:
        TOOLS[name] = path
    elif shutil.which(name):
        TOOLS[name] = shutil.which(name)
        warnings.warn(
            f"\n[ts_pipeline] Could not find user defined path for {name} '{user_defined}' in PATH. "
            f"Found \"{name}\" in {shutil.which(name)} in***REMOVED***ad"
        )
    else:
        warnings.warn(
            f"\n[ts_pipeline] Could not find '{name}' in PATH. "
            f"Please add \"{name}\" in your PATH or reconfigure ts_pipeline/config.py manually. "
            f"Aborting now. "
        )
        sys.exit()
        

# Convenience access for imports
xtb_path = TOOLS["xtb"]
orca_path = TOOLS["orca"]
pysis_path = TOOLS["pysis"]
crest_path = TOOLS["crest"]

print('All tools loaded:')
pprint.pprint(TOOLS, indent=4)
