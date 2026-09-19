# Collection discovery benchmark

## find-events routed-novel-v32: 89.0 s wall, 11.87 GiB peak
```
{
 "candidate_entities": 2976274,
 "retained_entities": 205735,
 "archives_total": 8,
 "archives_opened": 8,
 "max_candidates": 10000000,
 "max_candidates_considered": 100000000,
 "max_routed_entries": 1000000000,
 "routed_target_associations": 22446538,
 "routed_chunk_postings": 229799,
 "max_exact_match_attempts": 2500000000,
 "exact_match_attempts": 324154860,
 "annotation_excluded_splice_candidates": 110,
 "planned_compressed_bytes": 588039413,
 "actual_archive_bytes_read": 1017509893,
 "source_archive_identity_bytes_read": 935628,
 "collection_sidecar_bytes_read": 9605128,
 "terminal_tail_available_archives": 0,
 "terminal_tail_unavailable_archives": 0,
 "terminal_tail_routed_chunks": 0
}
```
## find-events coord-novel-v32: 93.4 s wall, 11.73 GiB peak
```
{
 "candidate_entities": 2976274,
 "retained_entities": 205735,
 "archives_total": 8,
 "archives_opened": 8,
 "max_candidates": 10000000,
 "max_candidates_considered": 100000000,
 "max_routed_entries": 1000000000,
 "routed_target_associations": 22446538,
 "routed_chunk_postings": 229799,
 "max_exact_match_attempts": 2500000000,
 "exact_match_attempts": 324154860,
 "annotation_excluded_splice_candidates": 110,
 "planned_compressed_bytes": 588039413,
 "actual_archive_bytes_read": 1017509893,
 "source_archive_identity_bytes_read": 935628,
 "collection_sidecar_bytes_read": 8946933,
 "terminal_tail_available_archives": 0,
 "terminal_tail_unavailable_archives": 0,
 "terminal_tail_routed_chunks": 0
}
```
## find-events coord-novel-v32: unavailable ('list' object has no attribute 'keys')
## find-events routed-all: 91.5 s wall, 11.77 GiB peak
```
{
 "candidate_entities": 2976274,
 "retained_entities": 209570,
 "archives_total": 8,
 "archives_opened": 8,
 "max_candidates": 10000000,
 "max_candidates_considered": 100000000,
 "max_routed_entries": 1000000000,
 "routed_target_associations": 22446538,
 "routed_chunk_postings": 229799,
 "max_exact_match_attempts": 2500000000,
 "exact_match_attempts": 324154860,
 "annotation_excluded_splice_candidates": 110,
 "planned_compressed_bytes": 588039413,
 "actual_archive_bytes_read": 1017509893,
 "source_archive_identity_bytes_read": 935628,
 "collection_sidecar_bytes_read": 9605128,
 "terminal_tail_available_archives": 0,
 "terminal_tail_unavailable_archives": 0,
 "terminal_tail_routed_chunks": 0
}
```
## find-events routed-all: 1:31.48 wall, 12.34 GB peak; total_seconds 85.9; retained_entities 209,570 (same candidates and bytes read as the other arms)
## naive per-archive discovery: 22.9 s total wall over 200 runs; per archive: A 1.6s, B 1.5s, C 1.0s, D 5.5s, E 1.2s, F 1.4s, G 9.4s, H 1.2s