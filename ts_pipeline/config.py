import shutil, sys, warnings, pprint

# Tools to check (name: default fallback)
_REQUIRED_TOOLS = {
    "xtb"  : "xtb",
    "orca" : "orca.run",
    "pysis": "pysis",
    "crest": "crest"
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