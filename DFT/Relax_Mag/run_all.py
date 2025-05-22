#!/usr/bin/env python3
import os
import subprocess

def run_all_main_py(base_dir=None):
    """Recursively find and run every main.py under base_dir (defaults to cwd)."""
    if base_dir is None:
        base_dir = os.getcwd()

    for root, dirs, files in os.walk(base_dir):
        if "main.py" in files and "slurm_runfile.sh" not in files:
            print(f"▶️  Running main.py in {root}")
            try:
                # Run `python3 main.py` in the directory containing main.py
                subprocess.run(
                    ["python3", "main.py"],
                    cwd=root,
                    check=True
                )
            except subprocess.CalledProcessError as e:
                print(f"❌  main.py in {root} exited with code {e.returncode}")

if __name__ == "__main__":
    run_all_main_py()
