#!/usr/bin/env python3
import subprocess
import yaml
import os
import platform
from datetime import datetime, timezone

def get_git_info():
    """Collect commit hash, branch, and remote URL."""
    def run(cmd):
        return subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode().strip()
    try:
        return {
            "commit": run(["git", "rev-parse", "HEAD"]),
            "branch": run(["git", "rev-parse", "--abbrev-ref", "HEAD"]),
            "remote": run(["git", "config", "--get", "remote.origin.url"])
        }
    except subprocess.CalledProcessError:
        return {"error": "Not a Git repository or Git not available"}

def get_env_info():
    """Collect Python version, OS info, and UTC timestamp."""
    return {
        "python_version": platform.python_version(),
        "system": platform.system(),
        "release": platform.release(),
        "timestamp_utc": datetime.now(timezone.utc).isoformat()
    }

def get_dependencies():
    """Capture installed Python packages via pip freeze."""
    try:
        deps = subprocess.check_output(["pip", "freeze"]).decode("utf-8").splitlines()
        return deps
    except Exception:
        return ["Could not retrieve dependencies"]

def save_metadata(data, path="prod/metadata.yaml"):
    """Write metadata to YAML file, creating directories if needed."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        yaml.dump(data, f, sort_keys=False)
    print(f"✅ Metadata saved to {path}")

if __name__ == "__main__":
    metadata = {
        "git_info": get_git_info(),
        "environment": get_env_info(),
        "dependencies": get_dependencies()
    }
    save_metadata(metadata)
