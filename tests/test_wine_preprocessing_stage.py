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

PREPARE_IMAGE = "rascommander/hec-ras-wine-precompute_6.6@sha256:" + "a" * 64


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
        "str_prepare_version": "6.6",
        "str_wine_profile": "/runtime/hecras-6.6",
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
        "runtime": {"kind": "wine", "hec_ras_version": "6.6"},
        "arguments": {"num_cores": 2},
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


@pytest.mark.parametrize("value", [None, "1", "2", "8"])
def test_host_core_configuration_matches_container_and_solver(generated_models, value):
    config_path, _, _, projects = generated_models
    config = configparser.ConfigParser(interpolation=None)
    config.read(config_path)
    if value is None:
        config["03_run_hec_ras"].pop("int_prepare_cores_per_job")
    else:
        config["03_run_hec_ras"]["int_prepare_cores_per_job"] = value
    with config_path.open("w", encoding="utf-8") as stream:
        config.write(stream)
    settings = stage._load_settings(config_path)
    command = stage._container_command(settings, projects[0], 300, "core-test")
    expected = value or "2"
    assert command[command.index("--cpus") + 1] == expected
    assert command[command.index("--num-cores") + 1] == expected


@pytest.mark.parametrize("value", ["0", "9", "-1", "2.5", "True"])
def test_host_rejects_unsupported_core_configuration(generated_models, value):
    config_path, _, _, _ = generated_models
    config = configparser.ConfigParser(interpolation=None)
    config.read(config_path)
    config["03_run_hec_ras"]["int_prepare_cores_per_job"] = value
    with config_path.open("w", encoding="utf-8") as stream:
        config.write(stream)
    with pytest.raises(ValueError):
        stage._load_settings(config_path)


@pytest.mark.parametrize("num_cores", [None, 8, True, "2", 2.0])
def test_host_requires_receipt_to_confirm_requested_core_count(generated_models, num_cores):
    config, _, _, projects = generated_models
    settings = stage._load_settings(config)
    project = projects[0]
    write_outputs(project)
    write_receipt(project, "cores")
    receipt_path = project.parent / ".ras-commander/runs/cores/prepare.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["arguments"]["num_cores"] = num_cores
    receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(ValueError, match="wrong HEC-RAS core count"):
        stage._validate_receipt(settings, project, "cores")


@pytest.mark.parametrize("num_cores", [None, 1, 2, 8])
@pytest.mark.parametrize("applied", [True, False])
def test_wine_worker_uses_ras_commander_preprocessing_api(tmp_path, monkeypatch, num_cores, applied):
    import types

    worker = load_container_script("windows_worker")
    checks = types.ModuleType("model_checks")
    checks.preflight = lambda *args: ([], tmp_path / "Model.g01.hdf", {}, [])
    checks.normalize_inputs = lambda paths: []
    checks.validate_outputs = lambda *args: {"geometry": {"area": {}}, "temporary_plan": {"area": {}}}
    monkeypatch.setitem(sys.modules, "model_checks", checks)
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
        def set_num_cores(path, count, *, ras_object, refresh_dataframes):
            calls["cores"] = (path, count, ras_object, refresh_dataframes)

        @staticmethod
        def set_2d_flow_options(path, *, cores, include_default, ras_object):
            calls["mesh_cores"] = (path, cores, include_default, ras_object)

        @staticmethod
        def get_plan_value(path, key, *, ras_object):
            assert path == plan_path
            assert key == "UNET D2 Cores"
            return calls["cores"][1] if applied else None

        @staticmethod
        def update_run_flags(path, *, ras_object, **flags):
            calls["flags"] = flags

    class FakePreprocess:
        @staticmethod
        def preprocess_plan(plan, **kwargs):
            assert "cores" in calls
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

    arguments = (
        str(tmp_path / "Model.prj"),
        "01",
        "C:/HEC-RAS/6.6/Ras.exe",
        "C:/runtime/ras_commander.whl",
        "d" * 64,
        300,
        True,
    )
    options = {} if num_cores is None else {"num_cores": num_cores}
    if not applied:
        with pytest.raises(RuntimeError, match="no effective UNET D2 Cores setting"):
            worker.run_worker(*arguments, **options)
        assert "prepare" not in calls
        return
    result = worker.run_worker(*arguments, **options)

    assert result["success"] is True
    assert calls["init"]["accept_tcu"] is True
    assert calls["clear"][0] == plan_path
    assert calls["cores"] == (plan_path, num_cores or 2, calls["init"]["ras_object"], False)
    assert calls["mesh_cores"] == (plan_path, num_cores or 2, True, calls["init"]["ras_object"])
    assert result["num_cores"] == (num_cores or 2)
    assert calls["flags"] == {"geometry_preprocessor": True}
    assert calls["prepare"][1]["max_wait"] == 300
    assert calls["prepare"][1]["clear_existing"] is True
    assert calls["prepare"][1]["fix_line_endings"] is True


@pytest.mark.parametrize("script", ["prepare", "windows_worker"])
@pytest.mark.parametrize("value", ["0", "9", "-1", "2.0", "two"])
def test_core_cli_rejects_out_of_range_and_noninteger_values(script, value):
    module = load_container_script(script)
    args = ["prepare", "--project", "model.prj"] if script == "prepare" else [
        "--project", "model.prj", "--plan", "01", "--ras-executable", "Ras.exe",
        "--ras-commander-wheel", "runtime.whl", "--expected-ras-commander-wheel-sha256",
        "d" * 64, "--timeout", "300", "--result", "result.json",
    ]
    with pytest.raises(SystemExit) as exc:
        module.build_parser().parse_args(args + ["--num-cores", value])
    assert exc.value.code == 2


@pytest.mark.parametrize("num_cores", [0, 9, -1, 2.0, "2", True, None])
def test_core_validation_precedes_runtime_or_model_access(tmp_path, monkeypatch, num_cores):
    prepare = load_container_script("prepare")
    worker = load_container_script("windows_worker")

    def forbidden(*args, **kwargs):
        pytest.fail("invalid core count must fail before runtime access")

    monkeypatch.setattr(prepare, "load_runtime", forbidden)
    monkeypatch.setattr(worker, "verify_ras_commander", forbidden)
    with pytest.raises(prepare.JobError, match="integer from 1 to 8"):
        prepare.run_prepare("missing.prj", "01", 300, False, "invalid", tmp_path,
                            "runtime.json", "6.5", num_cores=num_cores)
    with pytest.raises(ValueError, match="integer from 1 to 8"):
        worker.run_worker("missing.prj", "01", "Ras.exe", "runtime.whl", "d" * 64,
                          300, False, num_cores=num_cores)
    assert list(tmp_path.iterdir()) == []


@pytest.mark.parametrize("num_cores", [None, 1, 8])
@pytest.mark.parametrize("returned_cores", [None, "requested"])
def test_controller_propagates_cores_and_records_receipt(tmp_path, monkeypatch, num_cores, returned_cores):
    from types import SimpleNamespace

    prepare = load_container_script("prepare")
    project = write_project(tmp_path, "core-test")
    prefix = tmp_path / "seed"
    prefix.mkdir()
    monkeypatch.setenv("RAS2FIM_SCRATCH_ROOT", str(tmp_path / "scratch"))
    runtime = {
        "identity": {"kind": "wine", "hec_ras_version": "6.5"},
        "prefix": prefix, "wheel_sha": "d" * 64,
        "windows_python": "python.exe", "ras_executable": "Ras.exe",
        "ras_commander_wheel": "runtime.whl",
    }
    monkeypatch.setattr(prepare, "load_runtime", lambda *args: runtime)
    commands = []
    expected = num_cores or 2

    def runner(command, environment, timeout):
        if command[0] == "winepath":
            return SimpleNamespace(returncode=0, stdout=command[-1], stderr="")
        commands.append(command)
        assert command[command.index("--num-cores") + 1] == str(expected)
        outputs = prepare.expected_outputs(project, "01", "01")
        for output in outputs:
            output.write_bytes(b"prepared artifact")
        payload = {
            "success": True, "plan": "01", "geometry": "01",
            "num_cores": expected if returned_cores == "requested" else None,
            "timed_out": False, "full_result_copied": False, "signal_source": "bco",
            "ras_commander_wheel_sha256": "d" * 64,
            "hdf_validation": {"geometry": {"area": {}}, "temporary_plan": {"area": {}}},
            **dict(zip(("tmp_hdf_path", "b_file_path", "x_file_path"), map(str, outputs))),
        }
        Path(command[command.index("--result") + 1]).write_text(json.dumps(payload))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    options = {} if num_cores is None else {"num_cores": num_cores}
    success, receipt_path = prepare.run_prepare(
        project, "01", 300, False, "cores", project.parent, "runtime.json", "6.5",
        runner=runner, **options,
    )
    receipt = json.loads(receipt_path.read_text())
    assert len(commands) == 1
    assert receipt["arguments"]["num_cores"] == expected
    assert success is (returned_cores == "requested")
    if not success:
        assert "different core count" in receipt["error"]["message"]


@pytest.mark.parametrize("profile", [None, "", "/controlled/wine-6.6"])
def test_runtime_mount_is_optional_for_bundled_images(generated_models, profile):
    config_path, _, _, projects = generated_models
    config = configparser.ConfigParser(interpolation=None)
    config.read(config_path)
    section = config["03_run_hec_ras"]
    if profile is None:
        section.pop("str_wine_profile")
    else:
        section["str_wine_profile"] = profile
    with config_path.open("w", encoding="utf-8") as stream:
        config.write(stream)
    settings = stage._load_settings(config_path)
    command = stage._container_command(settings, projects[0], 300, "bundled-test")
    mounts = [command[index + 1] for index, item in enumerate(command) if item == "--mount"]
    runtime_mounts = [mount for mount in mounts if "dst=/runtime/wine-seed" in mount]
    assert runtime_mounts == (["type=bind,src=" + profile + ",dst=/runtime/wine-seed,readonly"] if profile else [])
    assert "--read-only" in command


@pytest.mark.parametrize("relative,excluded", [
    ("prefix/drive_c/Program Files (x86)/HEC/HEC-RAS/6.5/Ras.exe", False),
    ("prefix/drive_c/Program Files (x86)/HEC/HEC-RAS/7.0.1/Ras.exe", True),
    ("prefix/drive_c/users/rasworker/AppData/Local/Temp/setup.exe", True),
    ("prefix/drive_c/users/rasworker/Documents/model.prj", True),
    ("prefix/drive_c/users/rasworker/AppData/Local/HEC/settings.xml", False),
    ("prefix/drive_c/users/rasworker/AppData/Local/pip/cache/download", True),
    ("prefix/drive_c/windows/system32/kernel32.dll", False),
    ("prefix/user.reg", False),
    ("prefix/drive_c/a0da2ed8e5bc9f2258/Setup.exe", True),
    ("prefix/drive_c/ProgramData/HEC/Installation Cache/setup.msi", True),
    ("prefix/drive_c/Python311/Lib/site-packages/pyogrio/tests/fixtures/test.prj", False),
])
def test_release_profile_keeps_runtime_and_excludes_private_or_temporary_files(monkeypatch, relative, excluded):
    monkeypatch.syspath_prepend(str(REPO_ROOT / "containers" / "hecras-prepare"))
    bundle = load_container_script("bundle_profile")
    assert bundle.excluded_path(relative, "6.5") is excluded
