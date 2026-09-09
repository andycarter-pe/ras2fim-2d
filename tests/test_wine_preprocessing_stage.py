import configparser
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import threading
from unittest.mock import patch

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))
import stage_hecras_for_linux_wine_02b as stage

PREPARE_IMAGE = "rascommander/hec-ras-wine-precompute_6.5@sha256:" + "a" * 64


def load_container_script(name):
    path = REPO_ROOT / "containers" / "hecras-prepare" / (name + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_config(path, shared_root):
    config = configparser.ConfigParser(interpolation=None)
    config["03_run_hec_ras"] = {
        "str_linux_host": "worker@example",
        "str_windows_share": str(shared_root),
        "str_linux_share": "/shared",
        "str_prepare_image": PREPARE_IMAGE,
        "str_prepare_version": "6.5",
        "str_wine_profile": "/runtime/hecras-6.5",
        "int_prepare_cores_per_job": "2",
        "int_prepare_memory_gb": "6",
        "b_use_ntsync": "false",
    }
    with path.open("w", encoding="utf-8") as stream:
        config.write(stream)
    return path


def write_project(models, name):
    project_dir = models / name
    project_dir.mkdir()
    project = project_dir / (name + ".prj")
    project.write_text(
        "Proj Title=" + name + "\nCurrent Plan=p01\nPlan File=p01\n",
        encoding="utf-8",
    )
    project.with_suffix(".p01").write_text(
        "Plan Title=Test\nGeom File=g01\nFlow File=u01\n", encoding="utf-8"
    )
    return project


def write_outputs(project):
    for suffix in stage._REQUIRED_SUFFIXES:
        (project.parent / (project.stem + suffix)).write_bytes(
            (project.stem + suffix).encode()
        )


def write_receipt(project, run_id):
    artifacts = []
    for suffix in stage._REQUIRED_SUFFIXES:
        path = project.parent / (project.stem + suffix)
        artifacts.append({
            "path": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })
    receipt = {
        "schema": "ras-commander-job/v1",
        "run_id": run_id,
        "command": "prepare",
        "status": "succeeded",
        "project": project.name,
        "plan": "01",
        "geometry": "01",
        "runtime": {"kind": "wine", "hec_ras_version": "6.5"},
        "result": {
            "timed_out": False,
            "full_result_copied": False,
            "signal_source": "bco",
        },
        "artifacts": artifacts,
    }
    path = project.parent / ".ras-commander" / "runs" / run_id / "prepare.json"
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(receipt), encoding="utf-8")


@pytest.fixture
def generated_models(tmp_path):
    output = tmp_path / "output"
    models = output / "02_model_copies"
    models.mkdir(parents=True)
    (models / "source_terrain").mkdir()
    (models / "projection").mkdir()
    projects = [write_project(models, "model-one"), write_project(models, "model-two")]
    return write_config(tmp_path / "config.ini", tmp_path), output, models, projects


def test_all_preprocessors_finish_before_any_artifact_is_staged(generated_models):
    config, output, models, projects = generated_models
    completed = set()
    lock = threading.Lock()

    def fake_run(settings, project, timeout, run_id):
        write_outputs(project)
        with lock:
            completed.add(project.stem)

    original_stage = stage._stage_projects

    def assert_barrier(items, destination):
        assert completed == {project.stem for project in projects}
        original_stage(items, destination)

    with patch.object(stage, "_run_preparation", side_effect=fake_run), \
            patch.object(stage, "_stage_projects", side_effect=assert_barrier):
        stage.fn_prepare_hecras_for_linux(
            config, models, output / "02b_prep_for_ras", 4, 900, False
        )

    for project in projects:
        staged = output / "02b_prep_for_ras" / project.stem
        assert {path.name for path in staged.iterdir()} == {
            project.stem + suffix for suffix in stage._REQUIRED_SUFFIXES
        }


def test_failed_preparation_does_not_stage_a_partial_batch(generated_models):
    config, output, models, projects = generated_models

    def fake_run(settings, project, timeout, run_id):
        if project == projects[1]:
            raise RuntimeError("expected failure")
        write_outputs(project)

    destination = output / "02b_prep_for_ras"
    with patch.object(stage, "_run_preparation", side_effect=fake_run), \
            pytest.raises(RuntimeError, match="model-two: expected failure"):
        stage.fn_prepare_hecras_for_linux(config, models, destination, 4, 900, False)
    assert list(destination.iterdir()) == []


def test_receipt_validation_rehashes_all_three_handoff_files(generated_models):
    config, _, _, projects = generated_models
    settings = stage._load_settings(config)
    project = projects[0]
    write_outputs(project)
    write_receipt(project, "run-1")

    stage._validate_receipt(settings, project, "run-1")
    project.with_suffix(".b01").write_bytes(b"changed")
    with pytest.raises(ValueError, match="changed preparation artifact"):
        stage._validate_receipt(settings, project, "run-1")


def test_wine_worker_uses_ras_commander_preprocessing_api(tmp_path, monkeypatch):
    import types

    worker = load_container_script("windows_worker")
    calls = {}
    plan_path = tmp_path / "Model.p01"

    class FakeProject:
        pass

    class FakeResult:
        success = True
        plan_number = "01"
        geometry_number = "01"
        tmp_hdf_path = tmp_path / "Model.p01.tmp.hdf"
        b_file_path = tmp_path / "Model.b01"
        x_file_path = tmp_path / "Model.x01"
        elapsed_seconds = 1.5
        signal_source = "bco"
        full_result_copied = False
        timed_out = False
        error = None

        def __bool__(self):
            return True

    class FakeGeomPreprocessor:
        @staticmethod
        def clear_geompre_files(path, *, ras_object):
            calls["clear"] = (path, ras_object)

    class FakePlan:
        @staticmethod
        def get_plan_path(plan, *, ras_object):
            return plan_path

        @staticmethod
        def update_run_flags(path, *, ras_object, **flags):
            calls["flags"] = flags

    class FakePreprocess:
        @staticmethod
        def preprocess_plan(plan, **kwargs):
            calls["prepare"] = (plan, kwargs)
            return FakeResult()

    def fake_init(project, **kwargs):
        calls["init"] = kwargs

    module = types.ModuleType("ras_commander")
    module.GeomPreprocessor = FakeGeomPreprocessor
    module.RasPlan = FakePlan
    module.RasPreprocess = FakePreprocess
    module.RasPrj = FakeProject
    module.init_ras_project = fake_init
    monkeypatch.setitem(sys.modules, "ras_commander", module)
    monkeypatch.setattr(
        worker,
        "verify_ras_commander",
        lambda *args: {
            "ras_commander_wheel_sha256": "d" * 64,
            "ras_commander_distribution_version": "0.99.2",
        },
    )

    result = worker.run_worker(
        str(tmp_path / "Model.prj"),
        "01",
        "C:/HEC-RAS/6.5/Ras.exe",
        "C:/runtime/ras_commander.whl",
        "d" * 64,
        300,
        True,
    )

    assert result["success"] is True
    assert calls["init"]["accept_tcu"] is True
    assert calls["clear"][0] == plan_path
    assert calls["flags"] == {"geometry_preprocessor": True}
    assert calls["prepare"][1]["max_wait"] == 300
    assert calls["prepare"][1]["clear_existing"] is True
