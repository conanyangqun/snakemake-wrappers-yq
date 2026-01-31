import difflib
from pathlib import Path
import subprocess
import os
import tempfile
import shutil
import pytest
import sys
import yaml
import filecmp
from itertools import chain

BRANCH_NAME = 'main'
DIFF_MASTER = os.environ.get("DIFF_MASTER", "false") == "true"
DIFF_LAST_COMMIT = os.environ.get("DIFF_LAST_COMMIT", "false") == "true"

if DIFF_MASTER or DIFF_LAST_COMMIT:
    compare = "HEAD^" if DIFF_LAST_COMMIT else f"origin/{BRANCH_NAME}"

    # check if wrapper is modified compared to master
    DIFF_FILES = set(
        subprocess.check_output(["git", "diff", compare, "--name-only"])
        .decode()
        .split("\n")
    )

CONTAINERIZED = os.environ.get("CONTAINERIZED", "false") == "true"


@pytest.fixture
def tmp_test_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d

        # cleanup environments to save disk space
        subprocess.check_call(
            f"for env in `conda env list | grep -P '{d}' | "
            "cut -f1 | tr -d ' '`; do conda env remove --yes --prefix $env; done",
            shell=True,
        )


@pytest.fixture
def run(tmp_test_dir):
    def _run(wrapper, cmd, check_log=None, compare_results_with_expected=None):
        wrapper_dir = Path(wrapper)

        is_meta_wrapper = wrapper.startswith("meta/")

        tmp_test_subdir = Path(tempfile.mkdtemp(dir=tmp_test_dir))
        origdir = os.getcwd()

        meta_path = os.path.join(wrapper, "meta.yaml")
        try:
            with open(meta_path) as f:
                meta = yaml.load(f, Loader=yaml.BaseLoader)
        except Exception:
            raise ValueError(f"Unable to load or parse {meta_path}.")

        if meta.get("blacklisted"):
            pytest.skip("wrapper blacklisted")

        dst = tmp_test_subdir / BRANCH_NAME

        os.symlink(origdir, dst)

        used_wrappers = []
        wrapper_file = "used_wrappers.yaml"
        if os.path.exists(os.path.join(wrapper, wrapper_file)):
            # is meta wrapper
            with open(os.path.join(wrapper, wrapper_file), "r") as wf:
                wf = yaml.load(wf, Loader=yaml.BaseLoader)
                used_wrappers = wf["wrappers"]
        else:
            used_wrappers.append(wrapper)

        if (DIFF_MASTER or DIFF_LAST_COMMIT) and not any(
            any(f.startswith(w) for f in DIFF_FILES)
            for w in chain(used_wrappers, [wrapper])
        ):
            pytest.skip("wrappers not modified")

        testdir = tmp_test_subdir / "test"

        if is_meta_wrapper:
            # make sure that the meta-wrapper is where we expect it
            for path in wrapper_dir.iterdir():
                if path.is_dir():
                    shutil.copytree(path, tmp_test_subdir / path.name)
                else:
                    shutil.copy(path, tmp_test_subdir)
        else:
            shutil.copytree(wrapper_dir / "test", testdir)

        # switch to test directory
        os.chdir(testdir)
        if os.path.exists(".snakemake"):
            shutil.rmtree(".snakemake")
        cmd += [
            "--conda-cleanup-pkgs",
            "--printshellcmds",
            "--show-failed-logs",
        ]
        if not is_meta_wrapper:
            # meta-wrappers define their specific wrapper versions
            cmd += [
                "--wrapper-prefix",
                f"file://{tmp_test_subdir}/{BRANCH_NAME}/",
            ]


        if CONTAINERIZED:
            # run snakemake in container
            cmd = [
                "sudo",
                "docker",
                "run",
                "-it",
                "-v",
                "{}:{}".format(os.getcwd(), "/workdir"),
                "snakemake/snakemake",
                " ".join(cmd),
            ]

        try:
            subprocess.check_call(cmd)
            if compare_results_with_expected:
                for generated, expected in compare_results_with_expected.items():
                    if not filecmp.cmp(generated, expected, shallow=False):
                        with open(generated) as genf, open(expected) as expf:
                            gen_lines = genf.readlines()
                            exp_lines = expf.readlines()
                        diff = "".join(
                            difflib.Differ().compare(gen_lines, exp_lines)
                        )
                        raise ValueError(
                            f"Unexpected results: {generated} != {expected}."
                            f"Diff:\n{diff}"
                        )
        except Exception as e:
            # go back to original directory
            os.chdir(origdir)
            logfiles = [
                os.path.join(d, f)
                for d, _, files in os.walk(os.path.join(testdir, "logs"))
                for f in files
            ]
            for path in logfiles:
                with open(path) as f:
                    msg = "###### Logfile: " + path + " ######"
                    print(msg, "\n")
                    print(f.read())
                    print("#" * len(msg))
            if check_log is not None:
                for f in logfiles:
                    check_log(open(f).read())
            else:
                raise e
        finally:
            # go back to original directory
            os.chdir(origdir)
        return tmp_test_subdir

    return _run


def test_flye(run):
    run(
        "bio/flye",
        [
            "snakemake",
            "--cores",
            "1",
            "--use-conda",
            "-F",
            "all"
        ]
    )

def test_dnaapler_all(run):
    run(
        "bio/dnaapler/all",
        [
            "snakemake",
            "--cores",
            "1",
            "--use-conda",
            "-F",
            "all"
        ]
    )

def test_bamstats_histogram(run):
    run(
        "bio/fastcat/bamstats",
        [
            "snakemake",
            "--cores",
            "1",
            "--use-conda",
            "-F",
            "all"
        ]
    )
