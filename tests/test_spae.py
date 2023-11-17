import subprocess
from pathlib import Path
import pytest

TEST_ROOTDIR = Path(__file__).parent
EXEC_ROOTDIR = Path(__file__).parent.parent

@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("temp")

@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir) 

def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None

# Add more tests as required
def test_spae_version(tmp_dir):
    exec_command("spae -v")

def test_spae_help(tmp_dir):
    exec_command("spae -h")

def test_spae_config(tmp_dir):
    exec_command("spae config")

def test_spae_citation(tmp_dir):
    exec_command("spae citation")

def test_spae_run_help(tmp_dir):
    exec_command("spae run -h")

def test_spae_config_help(tmp_dir):
    exec_command("spae config -h")

def test_spae_illumina(tmp_dir):
    exec_command(f"spae run --input {TEST_ROOTDIR}/data/illumina-subset --output {tmp_dir} --touch")

def test_spae_longread(tmp_dir):
    exec_command(f"spae run --input {TEST_ROOTDIR}/data/nanopore-subset --output {tmp_dir} --sequencing longread --touch")

