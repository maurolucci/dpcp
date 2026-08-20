# cfcol

Branch-and-Price and ILP implementations for DPCP (double partition coloring problem) and CFCP (conflict-free coloring problem).

## Related Work

This repository corresponds to the work:

Double Partition Coloring problem: a unified approach for conflict-free coloring of hypergraphs and other coloring problems in graphs

Authors:

- Mauro Lucci (UNR, CONICET) - mlucci@fceia.unr.edu.ar
- Graciela Nasini (UNR, CONICET) - nasini@fceia.unr.edu.ar
- Paola Tolomei (UNR, CONICET) - ptolomei@fceia.unr.edu.ar
- Luis Miguel Torres (EPN) - luis.torres@epn.edu.ec

## Requirements

- Linux
- g++ with C++20 support
- Boost Program Options
- IBM ILOG CPLEX and Concert
- Gurobi (for the Gurobi-based ILP solver)

Path assumptions are defined in [Makefile](Makefile):

- CPLEXDIR=/opt/ibm/ILOG/CPLEX_Studio2211/cplex
- CONCERTDIR=/opt/ibm/ILOG/CPLEX_Studio2211/concert
- GUROBI_HOME=/opt/gurobi1203/linux64

If your installation paths are different, update those variables before building.

## Build

From the project root:

	cd exactcolors
	make
	cd ..
	make

This builds the executable dpcp in the project root.

## Instance Prefix and Format

Input is provided as a prefix without the final extension.
For a prefix X, the solver reads:

- X.graph
- X.partP
- X.partQ

In this repository, most generated instances include .dpcp in the prefix itself.
For example, the prefix:

- instances/dpcp/er-2/r_N110_p0.25_n22_m22_i0.dpcp

corresponds to:

- instances/dpcp/er-2/r_N110_p0.25_n22_m22_i0.dpcp.graph
- instances/dpcp/er-2/r_N110_p0.25_n22_m22_i0.dpcp.partP
- instances/dpcp/er-2/r_N110_p0.25_n22_m22_i0.dpcp.partQ

Detailed format documentation is available in [instances/README.md](instances/README.md).

## Run

General usage:

	./dpcp -s SOLVER -f PREFIX1 [PREFIX2 ...] [options]

Solvers documented here:

- byp
- compact
- gurobi
- heur

## Real Examples

Single small DPCP instance:

	./dpcp -s byp -f instances/dpcp/small-tests/single_P_part_size1.dpcp

Random DPCP instance from er-2:

	./dpcp -s byp -f instances/dpcp/er-2/r_N110_p0.25_n22_m22_i0.dpcp -t 600

CFCP-derived instance set in one run:

	./dpcp -s compact -f instances/cfcp/open/myciel4.dpcp instances/cfcp/open/queen6_6.dpcp -n 2

Gurobi-based compact ILP with output directory:

	./dpcp -s gurobi -f instances/cfcp/open/huck.dpcp -o out/ -t 1200

## Main Parameters

- -s, --solver: solver name
- -f, --graph: list of input prefixes (without .graph/.partP/.partQ)
- -o, --out: output directory (optional)
- -t, --time: global time limit in seconds (default 900)
- -n, --repeat: number of repetitions (default 1)
- -v, --verbose: verbosity level (0, 1, 2; default 1)
- --ub: initial upper bound (default DBL_MAX)
- --relax: solve only the root node
- --relax-time: LP relaxation time limit in seconds when using --relax (default 0 means ignored)
- --preproc-off: disable preprocessing

Heuristic and pricing controls:

- --heur-initial: 0 none, 1 greedy 1-step, 2 greedy 2-step, 3 semi-greedy 2-step
- --heur-nodes: 0 none, 1 greedy 1-step, 2 greedy 2-step, 3 semi-greedy 2-step
- --heur-2step-variant: 1 DEG, 2 CDEG, 3 SCDEG
- --heur-semigreedy-alpha
- --heur-semigreedy-iter
- --heur-semigreedy-time
- --tree-search: 1 best-bound, 2 dfs
- --inherit-cols: 0..4 inheritance strategy
- --inherit-pool-max-cols
- --dummy-weight
- --pricing-method: 0..5
- --pricing-greedy-max-cols
- --pricing-max-cols-per-iter
- --pricing-greedy-alpha
- --pricing-exact-time
- --branching-variable: 1 FMS, 2 LNTT

For the complete and authoritative list:

	./dpcp --help

## Output

If -o is not provided, logs and stats are printed to stdout.

If -o DIR is provided, the following subdirectories are created:

- DIR/log
- DIR/debug
- DIR/stat
- DIR/sol
- DIR/col
- DIR/iter

Output files are named per run as:

- instance-stem-solver-run.log
- instance-stem-solver-run.debug
- instance-stem-solver-run.stat
- instance-stem-solver-run.sol
- instance-stem-solver-run.col
- instance-stem-solver-run.iter

If a stat file for the same instance-solver-run already exists, execution is skipped and stored stats are printed.

## Notes on Gurobi

- The solver option -s gurobi runs the compact ILP implementation backed by Gurobi.
- The project still links against CPLEX/Concert in the main build, so both toolchains must be correctly configured in [Makefile](Makefile).
- Gurobi runtime libraries must be visible to the linker and loader on your system.

## Project Structure

- [main.cpp](main.cpp): entry point and argument parsing
- [src](src): core algorithms
- [include](include): headers
- [exactcolors](exactcolors): coloring and MWIS routines
- [instances](instances): instance sets and generators
- [experiments](experiments): experiment scripts and analysis
- [Makefile](Makefile): build configuration

## License

See [LICENSE](LICENSE).
