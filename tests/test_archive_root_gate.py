from __future__ import annotations

import importlib.util
import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from archive_root_gate_common import (  # noqa: E402
    GateError,
    parse_archive,
    self_test_blake3,
    validate_artifact_manifest,
    write_artifact_manifest,
)


def load_script(module_name: str, file_name: str):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPTS / file_name)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


fixtures = load_script("archive_root_adversarial_for_test", "archive_root_adversarial.py")
panel = load_script("jshape_panel_for_test", "211_materialize_jshape_panel.py")
validator = load_script("archive_root_validator_for_test", "210_validate_archive_root.py")
driver = load_script("archive_root_driver_for_test", "209_benchmark_archive_root.py")


class ArchiveRootGateTests(unittest.TestCase):
    def test_collection_structural_counts_ignore_identity_encoding_bytes(self) -> None:
        left = {
            "segment_junctions": 10,
            "archive_routes": 20,
            "chunk_postings": 30,
            "raw_section_bytes": 40,
        }
        right = {**left, "raw_section_bytes": 41}
        self.assertEqual(
            validator.collection_structural_counts(left),
            validator.collection_structural_counts(right),
        )
        right["chunk_postings"] += 1
        self.assertNotEqual(
            validator.collection_structural_counts(left),
            validator.collection_structural_counts(right),
        )

    def test_replay_stderr_has_a_strict_mode_specific_grammar(self) -> None:
        cases = {
            "gene-stream": (
                "molecules 22156853 -> assigned 21197710 -> 10393179 UMIs in 3600217 entries | "
                "streaming | open+dict 0.1s, +anno 1.5s, +decode/replay 1.4s, total 3.8s\n"
            ),
            "gene-eager": (
                "molecules 22156853 -> assigned 21197710 -> 10393179 UMIs in 3600217 entries | "
                "eager | load 0.7s, +anno 1.5s, +replay 0.6s, total 3.4s\n"
            ),
            "velocity-stream": (
                "velocity: 22156853 molecules -> 14298196 UMIs in 4998654 (cell,gene) entries | "
                "load 0.5s, total 15.2s\n"
            ),
            "velocity-eager": (
                "velocity: 22156853 molecules -> 14298196 UMIs in 4998654 (cell,gene) entries | "
                "load 0.6s, total 14.9s\n"
            ),
        }
        for mode, text in cases.items():
            parsed = validator.parse_replay_stderr(text, mode)
            self.assertEqual(parsed["mode"], mode)
            self.assertGreater(parsed["counts"]["molecules"], 0)
        with self.assertRaises(GateError):
            validator.parse_replay_stderr(cases["gene-stream"].replace("22156853", "22156854"), "gene-stream")
        with self.assertRaises(GateError):
            validator.parse_replay_stderr(cases["gene-stream"] + "warning\n", "gene-stream")

    def test_replay_stderr_is_validated_in_exactness_and_resource_summary_passes(self) -> None:
        diagnostic = (
            "molecules 22156853 -> assigned 21197710 -> 10393179 UMIs in 3600217 entries | "
            "streaming | open+dict 0.1s, +anno 1.5s, +decode/replay 1.4s, total 3.8s\n"
        )
        timing = """Command being timed: "aie replay"
User time (seconds): 1.00
System time (seconds): 0.10
Elapsed (wall clock) time (h:mm:ss or m:ss): 0:01.00
Maximum resident set size (kbytes): 1024
File system inputs: 0
File system outputs: 0
Exit status: 0
"""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for label in ("v1", "v2"):
                directory = root / label
                directory.mkdir()
                (directory / "command.json").write_text(
                    json.dumps(
                        {
                            "schema": "gravlax.archive-root-command.v1",
                            "argv": ["aie", "replay"],
                            "cwd": str(root),
                            "environment": {"RAYON_NUM_THREADS": "24"},
                        },
                        sort_keys=True,
                    )
                    + "\n"
                )
                (directory / "stdout.txt").write_text("")
                (directory / "stderr.txt").write_text(diagnostic)
                (directory / "time.txt").write_text(timing)
            validator.command_files(root / "v1", replay_mode="gene-stream")
            summary = validator.paired_summary(
                [
                    {"arm": "v1", "label": "v1"},
                    {"arm": "v2", "label": "v2"},
                ],
                root,
                replay_mode="gene-stream",
            )
            self.assertEqual(summary["wall_seconds"]["v2_over_v1"], 1.0)
            self.assertEqual(summary["max_rss_kib"]["v2_over_v1"], 1.0)

    def test_artifact_manifest_uses_posix_lexical_path_order(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            first = run / "commands/v1-root-reverse/command.json"
            second = run / "commands/v1-root/time.txt"
            first.parent.mkdir(parents=True)
            second.parent.mkdir(parents=True)
            first.write_text("first\n")
            second.write_text("second\n")
            written = write_artifact_manifest(run)
            lines = (run / "artifact-manifest.tsv").read_text().splitlines()
            self.assertEqual(
                [line.split("\t", 1)[0] for line in lines[1:]],
                [
                    "commands/v1-root-reverse/command.json",
                    "commands/v1-root/time.txt",
                ],
            )
            self.assertEqual(validate_artifact_manifest(run), written)

    def test_driver_precreates_every_planned_output_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary) / "gate-a-run"
            expected = {
                ".",
                "archives",
                "archive-audits",
                "commands",
                "collections",
                "collections/build-benchmark",
                "replay",
                "replay/canonical",
                "replay/benchmark",
                "fixtures",
            }
            self.assertEqual(
                {
                    "." if path == run else path.relative_to(run).as_posix()
                    for path in driver.planned_output_directories(run)
                },
                expected,
            )
            driver.prepare_output_tree(run)
            self.assertEqual(
                {
                    "." if path == run else path.relative_to(run).as_posix()
                    for path in driver.planned_output_directories(run)
                    if path.is_dir()
                },
                expected,
            )
            with self.assertRaises(GateError):
                driver.prepare_output_tree(run)

    def test_dirty_or_draft_runs_are_never_promotable(self) -> None:
        hook = [Path("fixture-hook")]
        self.assertFalse(
            driver.is_promotable_run(
                "draft", code_dirty=False, scripts_dirty=False, fixture_hooks=hook
            )
        )
        self.assertFalse(
            driver.is_promotable_run(
                "final", code_dirty=True, scripts_dirty=False, fixture_hooks=hook
            )
        )
        self.assertFalse(
            driver.is_promotable_run(
                "final", code_dirty=False, scripts_dirty=True, fixture_hooks=hook
            )
        )
        self.assertFalse(
            driver.is_promotable_run(
                "final", code_dirty=False, scripts_dirty=False, fixture_hooks=[]
            )
        )
        self.assertTrue(
            driver.is_promotable_run(
                "final", code_dirty=False, scripts_dirty=False, fixture_hooks=hook
            )
        )

    def test_independent_blake3_vectors_and_streaming(self) -> None:
        self_test_blake3()

    def test_synthetic_canonical_layouts_roundtrip(self) -> None:
        cases = (
            [],
            [(b"x", 7, b"payload")],
            [(b"x" * 255, 7, b"payload")],
            [(b"one", 1, b"a"), (b"unknown.optional", 2, b"bc")],
        )
        for sections in cases:
            value = fixtures.make_archive(b"gravlax-aie-directory-root-v2\0", sections)
            with tempfile.NamedTemporaryFile(suffix=".aie") as handle:
                handle.write(value)
                handle.flush()
                parsed = parse_archive(
                    Path(handle.name), b"gravlax-aie-directory-root-v2\0", verify_payloads=True
                )
            self.assertEqual(parsed.version, 2)
            self.assertEqual(len(parsed.sections), len(sections))

    def test_root_mutation_fails_closed(self) -> None:
        value = bytearray(
            fixtures.make_archive(
                b"gravlax-aie-directory-root-v2\0", [(b"x", 7, b"payload")]
            )
        )
        value[-36] ^= 1
        with tempfile.NamedTemporaryFile(suffix=".aie") as handle:
            handle.write(value)
            handle.flush()
            with self.assertRaises(GateError):
                parse_archive(
                    Path(handle.name), b"gravlax-aie-directory-root-v2\0", verify_payloads=False
                )

    def test_panel_varint_and_selection_golden(self) -> None:
        value, cursor = panel.decode_varint(bytes.fromhex("ac02"), 0)
        self.assertEqual((value, cursor), (300, 2))
        self.assertEqual(
            panel.selection_digest("chr1", 14660, 185187),
            "4acc33d50a408e411d0cfd21a6dcb19ebe71ebea852b0ca873e5edf36e242790",
        )

    def test_scientific_projection_removes_only_execution_fields(self) -> None:
        value = {
            "totals": {"umis": 2},
            "planning": {"total_seconds": 1.0},
            "samples": [
                {
                    "sample": "A",
                    "archive": "/tmp/a.aie",
                    "actual_archive_bytes_read": 123,
                    "umis": 2,
                }
            ],
        }
        self.assertEqual(
            validator.scientific_projection(value),
            {"totals": {"umis": 2}, "samples": [{"sample": "A", "umis": 2}]},
        )

    def test_excluded_inventory_covers_archives_collections_and_matrices(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            project = Path(temporary)
            run = project / "runs/gate-a"
            archives = {}
            for sample in (*validator.SAMPLES, "D0"):
                path = run / "archives" / f"{sample}.v2.aie"
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(sample.encode())
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                archives[sample] = {
                    "sealed_bytes": path.stat().st_size,
                    "sealed_sha256": digest,
                }
            collection = run / "collections/root.aicollection"
            collection.parent.mkdir(parents=True)
            collection.write_bytes(b"collection")
            binary = project / "env/aie"
            binary.parent.mkdir(parents=True)
            binary.write_bytes(b"binary")
            matrix = run / "replay/canonical/v2-gene-stream"
            matrix.mkdir(parents=True)
            matrix_file = matrix / "matrix.mtx"
            matrix_file.write_bytes(b"matrix")
            matrix_sha = hashlib.sha256(matrix_file.read_bytes()).hexdigest()
            inventory = validator.excluded_artifact_inventory(
                project_root=project,
                run=run,
                run_manifest={"bytes": 123, "manifest_sha256": "0" * 64},
                binary=binary,
                binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
                archives=archives,
                replay_artifacts={
                    matrix: {
                        "files": [
                            {
                                "path": "matrix.mtx",
                                "bytes": matrix_file.stat().st_size,
                                "sha256": matrix_sha,
                            }
                        ],
                        "tree_sha256": "1" * 64,
                    }
                },
            )
            self.assertIn("runs/gate-a\tgate_a_run_tree\t123", inventory)
            self.assertIn("runs/gate-a/archives/A.v2.aie\trooted_aie", inventory)
            self.assertIn("runs/gate-a/collections/root.aicollection\tcollection", inventory)
            self.assertIn("runs/gate-a/replay/canonical/v2-gene-stream\tmatrix_tree", inventory)


if __name__ == "__main__":
    unittest.main()
