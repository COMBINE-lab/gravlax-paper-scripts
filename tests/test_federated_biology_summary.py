#!/usr/bin/env python3
"""Focused tests for the bounded-memory federation result reader."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "181_summarize_federated_biology.py"
SPEC = importlib.util.spec_from_file_location("federated_biology_summary", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class StreamingEventDocumentTests(unittest.TestCase):
    def test_visits_events_and_retains_top_level_metadata(self) -> None:
        document = {
            "coordinates": "0-based half-open",
            "events": [{"id": "event-a"}, {"id": "event-b", "value": [1, 2, 3]}],
            "planning": {"candidate_events": 3, "retained_events": 2},
            "min_row_informative": 20,
            "schema": "test.v1",
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "events.json"
            path.write_text(json.dumps(document, indent=2))
            events: list[dict] = []
            metadata = MODULE.visit_event_document(path, events.append)

        self.assertEqual(events, document["events"])
        self.assertNotIn("events", metadata)
        self.assertEqual(metadata["planning"], document["planning"])
        self.assertEqual(MODULE.metadata_row_threshold(metadata, path), 20)

    def test_rejects_truncated_document(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "truncated.json"
            path.write_text('{"events": [{"id": "event-a"}')
            with self.assertRaises(SystemExit):
                MODULE.visit_event_document(path, lambda event: None)


class CountSemanticsTests(unittest.TestCase):
    def test_legacy_count_semantics_are_preserved(self) -> None:
        summary = MODULE.ArmSummary(20028, 22953, None, [], {})
        self.assertEqual(MODULE.event_count_fields(summary), {"events_screened": 20028})

    def test_pushdown_count_semantics_are_explicit(self) -> None:
        summary = MODULE.ArmSummary(4703, 22953, 10, [], {})
        self.assertEqual(
            MODULE.event_count_fields(summary),
            {
                "catalogue_recurrent_events": 22953,
                "denominator_qualified_events": 4703,
                "legacy_summed_minimum_events": None,
                "legacy_summed_minimum_status": "not_evaluated_due_to_row_pushdown",
                "min_row_informative": 10,
            },
        )


if __name__ == "__main__":
    unittest.main()
