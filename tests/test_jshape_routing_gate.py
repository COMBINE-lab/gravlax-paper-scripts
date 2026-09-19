from __future__ import annotations

import copy
import importlib.util
import json
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def load_script(module_name: str, file_name: str):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPTS / file_name)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


common = load_script("jshape_routing_gate_common_for_test", "jshape_routing_gate_common.py")
driver = load_script("jshape_routing_driver_for_test", "212_benchmark_jshape_routing.py")
validator = load_script("jshape_routing_validator_for_test", "213_validate_jshape_routing.py")
adversarial = load_script("jshape_routing_adversarial_for_test", "214_jshape_routing_adversarial.py")


class JshapeRoutingGateTests(unittest.TestCase):
    def test_frozen_panel_and_metadata_are_self_consistent(self) -> None:
        panel = ROOT / "experiments/designs/archive-root-jshape-routing-panel.tsv"
        metadata = ROOT / "experiments/designs/archive-root-jshape-routing-panel.metadata.json"
        cases, value = common.load_frozen_panel(panel, metadata)
        self.assertEqual(len(cases), 96)
        self.assertEqual(len({case.locus for case in cases}), 96)
        self.assertEqual(
            {name: sum(case.stratum == name for case in cases) for name in common.STRATA},
            {name: 24 for name in common.STRATA},
        )
        self.assertEqual(value["panel"]["sha256"], common.PANEL_SHA256)
        self.assertEqual(
            {row["id"]: (row["bytes"], row["sha256"], row["archive_root"]) for row in value["archives"]},
            {
                sample: (
                    common.ROOTED_ARCHIVES[sample]["bytes"],
                    common.ROOTED_ARCHIVES[sample]["sha256"],
                    common.ROOTED_ARCHIVES[sample]["archive_root"],
                )
                for sample in common.SAMPLES
            },
        )

    def test_only_source_frozen_final_adapter_can_promote(self) -> None:
        adapter = common.validate_adapter(copy.deepcopy(common.DEFAULT_ADAPTER))
        self.assertEqual(adapter["status"], "FINAL")
        self.assertTrue(common.adapter_is_frozen_final(adapter))
        self.assertTrue(
            common.is_promotable_run(
                "final",
                adapter=adapter,
                code_dirty=False,
                scripts_dirty=False,
                fixture_hooks=[Path("fixture")],
                required_fixture_hooks=[Path("fixture")],
            )
        )
        self.assertFalse(
            common.is_promotable_run(
                "final",
                adapter=adapter,
                code_dirty=False,
                scripts_dirty=False,
                fixture_hooks=[Path("alternate")],
                required_fixture_hooks=[Path("fixture")],
            )
        )
        override = copy.deepcopy(adapter)
        override["query_route_extra"] = ["--not-frozen"]
        self.assertFalse(common.adapter_is_frozen_final(override))

    def test_frozen_schedule_has_eight_alternating_paired_blocks(self) -> None:
        build = common.build_schedule()
        self.assertEqual(len(build), 16)
        for block in range(1, 9):
            rows = [row for row in build if row["block"] == block]
            expected = ["fallback", "route"] if block % 2 else ["route", "fallback"]
            self.assertEqual([row["arm"] for row in rows], expected)
            self.assertEqual([row["position"] for row in rows], [1, 2])
        panel = common.build_schedule(f"panel-{index:03d}" for index in range(1, 97))
        self.assertEqual(len(panel), 96 * 16)
        self.assertEqual(len({row["label"] for row in panel}), 96 * 16)

    def test_command_adapter_isolated_and_shell_free(self) -> None:
        adapter = common.DEFAULT_ADAPTER
        archives = {sample: Path(f"/{sample}.aie") for sample in common.SAMPLES}
        fallback = common.build_command(Path("/aie"), archives, common.SAMPLES, Path("/f.aic"), "fallback", adapter)
        route = common.build_command(Path("/aie"), archives, common.SAMPLES, Path("/r.aic"), "route", adapter)
        extension = common.build_command(
            Path("/aie"), archives, common.SAMPLES[4:], Path("/e.aic"), "route", adapter,
            base=Path("/base.aic"),
        )
        self.assertNotIn("--shape-routes", fallback)
        self.assertIn("--shape-routes", route)
        self.assertEqual(extension.count("--base"), 1)
        self.assertEqual(extension[extension.index("--base") + 1], "/base.aic")
        case = common.QueryCase("panel-001", "junction", "chr1:10-20", "one_archive")
        self.assertEqual(
            common.query_command(Path("/aie"), Path("/r.aic"), case, "route", adapter),
            ["/aie", "collection", "junction", "/r.aic", "chr1:10-20", "--json", "--explain", "--top", "0"],
        )
        self.assertEqual(
            common.inspect_command(Path("/aie"), Path("/r.aic"), "route_full", adapter)[-1],
            "--verify-routes",
        )
        with self.assertRaises(common.GateError):
            common.render(["{unknown}"], {"aie": "/aie"})

    def test_scientific_projection_removes_only_execution_fields(self) -> None:
        value = {
            "schema": "x",
            "totals": {"umis": 7, "cells": 3},
            "samples": [
                {
                    "sample": "A",
                    "umis": 7,
                    "archive": "/a.aie",
                    "actual_archive_bytes_read": 99,
                    "chunks_decoded": 2,
                    "shape_route_used": True,
                    "shape_route_section": "s.0.10",
                    "shape_route_sidecar_bytes_read": 12,
                }
            ],
            "planning": {"total_logical_bytes_read": 101},
            "elapsed_seconds": 0.1,
        }
        self.assertEqual(
            common.scientific_projection(value),
            {"schema": "x", "totals": {"umis": 7, "cells": 3}, "samples": [{"sample": "A", "umis": 7}]},
        )

    def test_route_section_and_nearest_rank_contracts(self) -> None:
        self.assertTrue(common.canonical_route_section_name("s.0.0"))
        self.assertTrue(common.canonical_route_section_name("s.7.4294967295"))
        for bad in ("s.00.1", "s.1.-1", "j.0.1", "s.1.2.extra"):
            self.assertFalse(common.canonical_route_section_name(bad))
        self.assertEqual(common.nearest_rank(list(range(1, 97)), 0.95), 92)

    def test_output_tree_preflight_enumerates_every_child_parent(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary) / "run"
            common.prepare_output_tree(run)
            self.assertTrue(all(path.is_dir() for path in common.planned_output_directories(run)))
            self.assertIn(run / "collections/benchmark", common.planned_output_directories(run))
            self.assertIn(run / "commands/query-panel-benchmark", common.planned_output_directories(run))

    def test_query_accounting_requires_exact_sums_routes_and_phases(self) -> None:
        adapter = copy.deepcopy(common.DEFAULT_ADAPTER)
        adapter["schemas"]["junction"] = "test.junction.v1"
        case = common.QueryCase("panel-001", "junction", "chr1:10-20", "one_archive")
        value = {
            "schema": "test.junction.v1",
            "totals": {"umis": 1, "cells": 1},
            "planning": {
                "source_archive_identity_bytes_read": 10,
                "source_archive_execution_bytes_read": 90,
                "source_archive_bytes_read": 100,
                "collection_sidecar_bytes_read": 20,
                "shape_route_sidecar_payload_bytes_read": 5,
                "total_logical_bytes_read": 120,
                "route_blocks_loaded": 1,
                "routed_archives": 1,
                "fallback_archives": 0,
                "archives_opened": 1,
                "collection_load_seconds": 0.01,
                "identity_guard_seconds": 0.01,
                "route_planning_seconds": 0.01,
                "source_execution_seconds": 0.02,
                "total_seconds": 0.06,
            },
        }
        result = validator.validate_query_json(value, arm="route", case=case, adapter=adapter)
        self.assertEqual((result["source_bytes"], result["sidecar_bytes"], result["total_bytes"]), (100, 20, 120))
        bad = copy.deepcopy(value)
        bad["planning"]["total_logical_bytes_read"] = 121
        with self.assertRaises(common.GateError):
            validator.validate_query_json(bad, arm="route", case=case, adapter=adapter)

    def test_absent_source_payload_gate_distinguishes_identity_from_execution_io(self) -> None:
        accounting = {
            ("absent", arm, form): {
                "source_execution_bytes": 0,
                "source_bytes": 123,
            }
            for arm in ("fallback", "route")
            for form in ("root", "extension")
        }
        validator.require_absent_source_payload_zero(accounting)
        accounting[("absent", "route", "root")]["source_execution_bytes"] = 1
        with self.assertRaises(common.GateError):
            validator.require_absent_source_payload_zero(accounting)

    def test_canonical_build_reports_sorted_and_flattened_sample_ids(self) -> None:
        self.assertEqual(validator.canonical_build_reported_samples("base"), common.SAMPLES[:4])
        for form in ("root", "repeat", "reverse", "extension"):
            self.assertEqual(validator.canonical_build_reported_samples(form), common.SAMPLES)
        with self.assertRaises(common.GateError):
            validator.canonical_build_reported_samples("unknown")

    def test_build_accounting_rejects_forbidden_payloads(self) -> None:
        adapter = copy.deepcopy(common.DEFAULT_ADAPTER)
        adapter["schemas"]["build"] = "test.build.v1"
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "route.aicollection"
            output.write_bytes(b"12345")
            value = {
                "schema": "test.build.v1",
                "output": str(output),
                "collection_format_version": 4,
                "shape_routes_requested": True,
                "new_archives": 1,
                "segment_junctions": 3,
                "archive_routes": 2,
                "chunk_postings": 2,
                "raw_section_bytes": 20,
                "file_bytes": 5,
                "shape_routes": {
                    "archives": 1,
                    "sections": 1,
                    "exact_spans": 2,
                    "pairs": 3,
                    "raw_bytes": 10,
                    "compressed_bytes": 2,
                },
                "elapsed_seconds": 0.05,
                "source_io": {
                    "identity_content_bytes_read": 0,
                    "total_bytes_read": 100,
                    "shape_route_payload_bytes_read": 30,
                    "shape_route_source_bytes_read": 40,
                    "sections_read": [
                        "chroms", "index.chunks", "index.jpost", "index.junctions",
                        "meta", "rans.tables", "shapes",
                    ],
                    "archives": [
                        {
                            "id": "A",
                            "format_version": 2,
                            "identity_scheme": "aie-directory-root-v2",
                            "identity_content_bytes_read": 0,
                            "total_bytes_read": 100,
                            "shape_route_payload_bytes_read": 30,
                            "shape_route_source_bytes_read": 40,
                            "sections_read": [
                                "chroms", "index.chunks", "index.jpost", "index.junctions",
                                "meta", "rans.tables", "shapes",
                            ],
                        }
                    ],
                },
            }
            result = validator.validate_build_json(
                value,
                arm="route",
                output=output,
                adapter=adapter,
                expected_samples=("A",),
            )
            self.assertEqual(result["route_compressed_bytes"], 2)

            extension_output = Path(temporary) / "route-extension.aicollection"
            extension_output.write_bytes(b"12345")
            extension = copy.deepcopy(value)
            extension["output"] = str(extension_output)
            extension["source_io"]["archives"][0]["id"] = "B"
            extension["source_io"]["archives"].insert(
                0,
                {
                    "id": "A",
                    "format_version": 2,
                    "identity_scheme": "aie-directory-root-v2",
                    "identity_content_bytes_read": 0,
                    "total_bytes_read": 10,
                    "shape_route_payload_bytes_read": 0,
                    "shape_route_source_bytes_read": 0,
                    "sections_read": [],
                },
            )
            extension["source_io"]["total_bytes_read"] += 10
            validator.validate_build_json(
                extension,
                arm="route",
                output=extension_output,
                adapter=adapter,
                expected_samples=("A", "B"),
                expected_new_samples=("B",),
            )

            bad = copy.deepcopy(value)
            bad["source_io"]["sections_read"] = ["c0", *value["source_io"]["sections_read"]]
            with self.assertRaises(common.GateError):
                validator.validate_build_json(
                    bad,
                    arm="route",
                    output=output,
                    adapter=adapter,
                    expected_samples=("A",),
                )

    def test_route_inspect_requires_exact_bindings_and_bounded_sorted_descriptors(self) -> None:
        adapter = copy.deepcopy(common.DEFAULT_ADAPTER)
        adapter["schemas"]["inspect"] = "test.inspect.v1"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            collection = root / "route.aicollection"
            collection.write_bytes(b"x" * 100)
            source_directories = {}
            archives = []
            for ordinal, sample in enumerate(common.SAMPLES):
                source_path = root / f"{sample}.v2.aie"
                shapes = types.SimpleNamespace(name="shapes", committed_blake3="a" * 64)
                source_directories[sample] = types.SimpleNamespace(
                    path=source_path,
                    file_bytes=common.ROOTED_ARCHIVES[sample]["bytes"],
                    sections=(shapes,),
                )
                archives.append(
                    {
                        "id": sample,
                        "path": str(source_path),
                        "bytes": common.ROOTED_ARCHIVES[sample]["bytes"],
                        "archive_format_version": 2,
                        "native_identity": {
                            "scheme": "aie-directory-root-v2",
                            "blake3": common.ROOTED_ARCHIVES[sample]["archive_root"],
                        },
                        "encoded_sections_identity": {
                            "scheme": "aie-encoded-sections-v1",
                            "blake3": "b" * 64,
                        },
                        "chunks": 4,
                        "shape_routes": {
                            "codec_version": 1,
                            "archive_ordinal": ordinal,
                            "source_root": common.ROOTED_ARCHIVES[sample]["archive_root"],
                            "shapes_digest": "a" * 64,
                            "n_shapes": 10,
                            "sections": 1,
                            "exact_spans": 2,
                            "pairs": 3,
                            "compressed_bytes": 10,
                            "reconstruction_verified": True,
                            "descriptors": [
                                {
                                    "first_span": 10,
                                    "last_span": 20,
                                    "section_name": f"s.{ordinal}.10",
                                    "exact_spans": 2,
                                    "pairs": 3,
                                    "compressed_bytes": 10,
                                }
                            ],
                        },
                    }
                )
            value = {
                "schema": "test.inspect.v1",
                "format_version": 4,
                "path": str(collection),
                "file_bytes": 100,
                "layers": [
                    {
                        "path": str(collection),
                        "format_version": 4,
                        "root_digest": "c" * 64,
                        "archives": 8,
                        "junction_rows": 5,
                        "junction_segments": 2,
                        "shape_routes": {
                            "archives": 8,
                            "sections": 8,
                            "exact_spans": 16,
                            "pairs": 24,
                            "compressed_bytes": 80,
                            "reconstruction_verified": True,
                        },
                    }
                ],
                "reference": {"stamped": True, "algo": "x", "digest": "d" * 64},
                "chromosomes": ["chr1"],
                "archives": archives,
                "index": {
                    "segment_junction_rows": 5,
                    "global_junctions": None,
                    "archive_routes": 5,
                    "chunk_postings": 5,
                    "interval_chunks": 32,
                    "junction_bin_bp": 65536,
                    "junction_segments": 2,
                    "shape_route_archives": 8,
                    "shape_route_sections": 8,
                    "shape_route_exact_spans": 16,
                    "shape_route_pairs": 24,
                    "shape_route_compressed_bytes": 80,
                },
                "guard": {
                    "filesystem_identity": "size + nanosecond mtime + nanosecond ctime + device + inode",
                    "content_digest_recorded": True,
                    "content_digest_verified": False,
                    "shape_route_payloads_verified": True,
                    "shape_route_reconstruction_verified": True,
                },
                "io": {
                    "collection_sidecar_bytes_read": 90,
                    "source_identity_content_bytes_read": 0,
                    "source_identity_total_bytes_read": 80,
                    "route_verification_source_bytes_read": 100,
                },
                "elapsed_seconds": 0.2,
            }
            summary = validator.validate_collection_inspect(
                value,
                arm="route",
                full=True,
                adapter=adapter,
                collection=collection,
                layer_paths=[collection],
                source_directories=source_directories,
            )
            self.assertEqual(summary["archives"], 8)
            bad = copy.deepcopy(value)
            bad["archives"][0]["shape_routes"]["descriptors"][0]["exact_spans"] = 257
            with self.assertRaises(common.GateError):
                validator.validate_collection_inspect(
                    bad,
                    arm="route",
                    full=True,
                    adapter=adapter,
                    collection=collection,
                    layer_paths=[collection],
                    source_directories=source_directories,
                )

    def test_fixture_panel_matches_every_frozen_integrity_case(self) -> None:
        required_rejections = {
            "unknown_route_version",
            "incorrect_source_root",
            "incorrect_shapes_digest",
            "out_of_range_shape_id",
            "unsorted_pairs",
            "duplicate_pairs",
            "wrong_archive_ordinal",
            "wrong_span_bucket",
            "count_overflow",
            "checked_addition_overflow",
            "truncation",
            "checksum_mismatch",
            "trailing_bytes",
        }
        self.assertTrue(required_rejections.issubset(common.REQUIRED_FIXTURE_CASES))
        self.assertTrue(all(common.REQUIRED_FIXTURE_CASES[name] == "reject" for name in required_rejections))
        self.assertEqual(len(common.REQUIRED_FIXTURE_CASES), 26)
        self.assertEqual(
            set(common.REQUIRED_FIXTURE_CASES),
            set(common.EXTERNAL_FIXTURE_CASES)
            | {case for cases in common.RUST_FIXTURE_TESTS.values() for case in cases},
        )

    def test_adversarial_block_codec_and_manifest_locator_are_independent(self) -> None:
        block = {
            "version": 1,
            "archive": 0,
            "n_shapes": 4,
            "first": 10,
            "last": 20,
            "n_spans": 2,
            "spans": [
                {"span": 10, "pairs": [[0, 4], [0, 8], [2, 1]]},
                {"span": 20, "pairs": [[3, 2]]},
            ],
        }
        raw = adversarial.encode_block(block)
        self.assertEqual(adversarial.decode_block(raw), block)

        manifest = bytearray([0, 0, 0])

        def varint(value: int) -> None:
            adversarial.put_varint(manifest, value)

        def string(value: str) -> None:
            encoded = value.encode()
            varint(len(encoded))
            manifest.extend(encoded)

        varint(1)
        string("chr1")
        string("a" * 64)
        varint(1)
        string("A")
        string("/A.aie")
        for value in (100, 1, 2, 3, 4, 5, 6, 2):
            varint(value)
        for value in ("aie-directory-root-v2", "b" * 64, "c" * 64):
            string(value)
        varint(0)
        manifest.append(1)
        varint(1)
        manifest.extend(bytes.fromhex("b" * 64))
        manifest.extend(bytes.fromhex("d" * 64))
        varint(4)
        varint(1)
        varint(10)
        varint(20)
        string("s.0.10")
        for _ in range(4):
            varint(0)
        locations = adversarial.binding_locations(bytes(manifest))
        self.assertEqual(len(locations), 1)
        self.assertEqual(locations[0]["codec"], 1)
        self.assertEqual(manifest[locations[0]["source_root"]], 0xBB)

    def test_pass_outputs_roll_back_inventory_if_result_publish_fails(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result, inventory = root / "result.json", root / "inventory.tsv"
            real_replace = validator.os.replace

            def fail_second(source: Path, destination: Path) -> None:
                if Path(destination) == result:
                    raise OSError("synthetic result publish failure")
                real_replace(source, destination)

            with mock.patch.object(validator.os, "replace", side_effect=fail_second):
                with self.assertRaises(OSError):
                    validator.write_pass_outputs(
                        result_path=result,
                        result={"status": "PASS"},
                        inventory_path=inventory,
                        inventory="header\n",
                    )
            self.assertFalse(result.exists())
            self.assertFalse(inventory.exists())
            self.assertFalse(result.with_name(result.name + ".tmp").exists())
            self.assertFalse(inventory.with_name(inventory.name + ".tmp").exists())


if __name__ == "__main__":
    unittest.main()
