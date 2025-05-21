import sys
from pathlib import Path

# Add project_root/scripts to sys.path
scripts_path = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(scripts_path))
