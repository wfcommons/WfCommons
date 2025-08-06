import subprocess
from pathlib import Path

def main():
    bin_path = Path(__file__).parent / "bin" / "cpu-benchmark"
    subprocess.run([str(bin_path)] + sys.argv[1:])

