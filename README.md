# cfcol

Branch-and-Price implementation for DPCP (double partition coloring problem) and CFCP (conflict-free coloring problem).

## Requirements

- Linux
- g++ with C++20 support
- IBM ILOG CPLEX and Concert
- Boost Program Options

The Makefile assumes CPLEX/Concert paths installed at:

- /opt/ibm/ILOG/CPLEX_Studio2211/cplex
- /opt/ibm/ILOG/CPLEX_Studio2211/concert

If your installation is different, edit [Makefile](Makefile).

## Build

From the project root:

	cd exactcolors
	make
	cd ..
	make

This generates the dpcp executable in the root directory.

## Input Format

Each instance is passed by file prefix. For a prefix X, the following files must exist:

- X.graph
- X.partP
- X.partQ

Example prefix:

- instances/test

## Run

General usage:

	./dpcp -s SOLVER -f PREFIX1 [PREFIX2 ...] [options]

Available solvers:

- byp
- compact
- heur
- feas-enum
- feas-ilp

Examples:

	./dpcp -s byp -f instances/test_small

	./dpcp -s byp -f instances/test_small instances/test_medium -t 600 -n 3

	./dpcp -s compact -f instances/test_small -o out/

## Main Parameters

- -s, --solver: solver type
- -f, --graph: list of instance prefixes
- -o, --out: optional output directory
- -t, --time: time limit in seconds (default 900)
- -v, --verbose: verbosity level 0..2
- -n, --repeat: number of repetitions
- --relax: solve root node only
- --ub: initial upper bound
- --preproc-off: disable preprocessing

Branch-and-Price parameters:

- --tree-search: 1 best-bound, 2 dfs
- --heur-root, --heur-nodes: node heuristic type
- --feas-root, --feas-nodes: feasibility check type
- --inherit-cols: column inheritance strategy
- --pricing-method: pricing strategy
- --pricing-greedy-max-cols: number of greedy pricing attempts
- --pricing-max-cols-per-iter: max columns added per pricing iteration
- --pricing-greedy-alpha: greedy pricing alpha
- --pricing-exact-time: max time for exact pricing
- --branching-variable: 1 FMS, 2 LNTT

For the full list:

	./dpcp --help

## Output

If -o is not provided, output is printed to stdout.

If -o DIR is provided, these subdirectories are created:

- DIR/log
- DIR/debug
- DIR/stat
- DIR/sol
- DIR/col
- DIR/iter

## Project Structure

- [main.cpp](main.cpp): program entry point and argument parsing
- [src](src): algorithm implementation
- [include](include): headers
- [exactcolors](exactcolors): external code for coloring/MWIS routines
- [instances](instances): test instances
- [experiments](experiments): scripts and experimental analysis
- [Makefile](Makefile): build configuration

## License

See [LICENSE](LICENSE).
