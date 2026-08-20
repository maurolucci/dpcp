# DPCP Instance Format

This directory contains DPCP instances split into four files per prefix:

- `<prefix>.dpcp.graph`
- `<prefix>.dpcp.partP`
- `<prefix>.dpcp.partQ`
- `<prefix>.dpcp.dict`

Example prefix:

- `r_N110_p0.25_n22_m22_i0`

Associated files:

- `r_N110_p0.25_n22_m22_i0.dpcp.graph`
- `r_N110_p0.25_n22_m22_i0.dpcp.partP`
- `r_N110_p0.25_n22_m22_i0.dpcp.partQ`
- `r_N110_p0.25_n22_m22_i0.dpcp.dict`

All vertex and partition identifiers are non-negative integers with 0-based indexing.

## 1) .dpcp.graph file

Describes the simple undirected graph.

### Header

First line:

`N:m`

- `N`: number of vertices (expected IDs: `0, 1, ..., N-1`)
- `m`: number of edges

### Body

Then exactly `m` lines follow, one per edge:

`u v`

- `u` and `v` are vertex IDs
- In this repository generators, each edge is written only once with `u < v`
- No self-loops are written (`u != v`)

## 2) .dpcp.partP file

Describes the first partition of vertices, with `n` parts.

### Header

First line:

`N:n`

- `N`: number of vertices (must match `.graph`)
- `n`: number of parts in `P`

### Body

Then `n` lines follow, one per part:

`pi k v1 v2 ... vk`

- `pi`: part ID in `P` (expected in `0..n-1`)
- `k`: number of vertices in that part
- `v1..vk`: IDs of vertices that belong to `P[pi]`

## 3) .dpcp.partQ file

Describes the second partition of vertices, with `m` parts.

### Header

First line:

`N:m`

- `N`: number of vertices (must match `.graph`)
- `m`: number of parts in `Q`

### Body

Then `m` lines follow, one per part:

`qj k v1 v2 ... vk`

- `qj`: part ID in `Q` (expected in `0..m-1`)
- `k`: number of vertices in that part
- `v1..vk`: IDs of vertices that belong to `Q[qj]`

## 4) .dpcp.dict file

Stores the mapping of each vertex to its double label `(pi, qj)`.

### Header

First line:

`N:n:m`

- `N`: number of vertices
- `n`: number of parts in `P`
- `m`: number of parts in `Q`

### Body

Then `N` lines follow:

`v pi qj`

- `v`: vertex ID
- `pi`: part of `P` to which `v` belongs
- `qj`: part of `Q` to which `v` belongs

This file is redundant with respect to `.partP` and `.partQ`, but it is useful for traceability and conversions.

## Expected Consistency Across Files

For a well-formed instance:

1. The headers containing `N` must match (`.graph`, `.partP`, `.partQ`, and also `.dict` if present).
2. In `.partP`, each vertex must appear exactly once across all parts.
3. In `.partQ`, each vertex must appear exactly once across all parts.
4. If `.dict` exists, for each vertex `v`, the pair `(pi, qj)` must match its membership in `.partP` and `.partQ`.

## Solver Usage Note

The main DPCP instance parser consumes `.graph`, `.partP`, and `.partQ`. The `.dict` file is used as an auxiliary file in generation and conversion scripts.