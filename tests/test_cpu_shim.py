"""Exercise the compiled preload shim in a fresh Linux process."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


@pytest.mark.skipif(sys.platform != "linux", reason="Linux preload shim")
def test_cpu_shim_maps_restricted_affinity_with_concurrent_first_use(tmp_path):
    compiler = shutil.which("cc")
    if compiler is None:
        pytest.skip("C compiler required")
    cpus = sorted(os.sched_getaffinity(0))
    selected = [cpu for cpu in cpus if cpu > 0][:2]
    if len(selected) < 2:
        pytest.skip("Two nonzero host CPUs required")
    source = Path(__file__).resolve().parents[1] / "containers/hecras-prepare/cpushim.c"
    library = tmp_path / "cpushim.so"
    subprocess.run([compiler, "-shared", "-fPIC", "-O2", "-Wall", "-Wextra",
                    "-Werror", "-o", str(library), str(source), "-ldl", "-pthread"],
                   check=True)
    script = '''
import concurrent.futures, ctypes, os, threading
cpus = sorted(os.sched_getaffinity(0))
getcpu = ctypes.CDLL(None).sched_getcpu
barrier = threading.Barrier(32)
def check(i):
    barrier.wait()
    # Initialize while every thread still has the process CPU set.
    assert 0 <= getcpu() < len(cpus)
    index = i % len(cpus)
    os.sched_setaffinity(0, {cpus[index]})
    for _ in range(1000):
        assert getcpu() == index
with concurrent.futures.ThreadPoolExecutor(max_workers=32) as pool:
    list(pool.map(check, range(32)))
'''
    environment = dict(os.environ, LD_PRELOAD=str(library))
    # Repeat fresh process initialization; setting affinity on Python threads
    # in the child avoids unsafe preexec_fn use in a threaded test runner.
    script = "import os\nos.sched_setaffinity(0, " + repr(set(selected)) + ")\n" + script
    for _ in range(5):
        subprocess.run([sys.executable, "-c", script], env=environment,
                       check=True, timeout=30)
