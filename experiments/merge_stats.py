import os
import sys
from typing import Any, Dict, List

import pandas as pd


def get_solver(name: str, subdir: str) -> str:
    parts = subdir.split(os.path.sep)[2:-1]
    return name + "_" + "_".join(parts)


def parse_int(value: str) -> int:
    return int(value)


def parse_float(value: str) -> float:
    return float(value)


def parse_first_line(l0: List[str], subdir: str) -> Dict[str, Any]:
    # Newest format (with RL initialization counters):
    # ... nNodesInt,nNodesFrac,nNodesGcp,nNodesTrivial,
    # nNodesInfeas,nNodesInfeasPrepro,nNodesInfeasCheck,nNodesInfeasAux,
    # gcpAvgTime,nsol,nsolHeur,nsolLR,nsolGCP,nsolTrivial,
    # ninitSol,ninitDummy,ninit
    if len(l0) >= 34:
        return {
            "instance": l0[0],
            "solver": get_solver(l0[1], subdir),
            "run": parse_int(l0[2]),
            "nvertices": parse_int(l0[3]),
            "nedges": parse_int(l0[4]),
            "nP": parse_int(l0[5]),
            "nQ": parse_int(l0[6]),
            "nvars": parse_int(l0[7]),
            "ncons": parse_int(l0[8]),
            "state": l0[9],
            "terminationReason": l0[10],
            "time": parse_float(l0[11]),
            "nodes": parse_int(l0[12]),
            "nodesLeft": parse_int(l0[13]),
            "lb": parse_float(l0[14]),
            "ub": parse_float(l0[15]),
            "gap": parse_float(l0[16]),
            "nNodesInt": parse_int(l0[17]),
            "nNodesFrac": parse_int(l0[18]),
            "nNodesGcp": parse_int(l0[19]),
            "nNodesTrivial": parse_int(l0[20]),
            "nNodesInfeas": parse_float(l0[21]),
            "nNodesInfeasPrepro": parse_int(l0[22]),
            "nNodesInfeasCheck": parse_int(l0[23]),
            "nNodesInfeasAux": parse_int(l0[24]),
            "gcpAvgTime": parse_float(l0[25]),
            "nsol": parse_int(l0[26]),
            "nsolHeur": parse_int(l0[27]),
            "nsolLR": parse_int(l0[28]),
            "nsolGCP": parse_int(l0[29]),
            "nsolTrivial": parse_int(l0[30]),
            "ninitSol": parse_int(l0[31]),
            "ninitDummy": parse_int(l0[32]),
            "ninit": parse_int(l0[33]),
        }

    # New format (stats.cpp prior to node-count reordering):
    # instance,solver,run,nvertices,nedges,nP,nQ,nvars,ncons,state,
    # terminationReason,time,nodes,nodesLeft,lb,ub,gap,
    # nNodesInfeas,nNodesInfeasPrepro,nNodesInfeasCheck,nNodesInfeasAux,
    # gcpAvgTime,nsol,nsolHeur,nsolLR,nsolGCP,nsolTrivial,nNodesTrivial,
    # nNodesInt,nNodesFrac,nNodesGcp
    if len(l0) >= 31:
        return {
            "instance": l0[0],
            "solver": get_solver(l0[1], subdir),
            "run": parse_int(l0[2]),
            "nvertices": parse_int(l0[3]),
            "nedges": parse_int(l0[4]),
            "nP": parse_int(l0[5]),
            "nQ": parse_int(l0[6]),
            "nvars": parse_int(l0[7]),
            "ncons": parse_int(l0[8]),
            "state": l0[9],
            "terminationReason": l0[10],
            "time": parse_float(l0[11]),
            "nodes": parse_int(l0[12]),
            "nodesLeft": parse_int(l0[13]),
            "lb": parse_float(l0[14]),
            "ub": parse_float(l0[15]),
            "gap": parse_float(l0[16]),
            "nNodesInfeas": parse_float(l0[17]),
            "nNodesInfeasPrepro": parse_int(l0[18]),
            "nNodesInfeasCheck": parse_int(l0[19]),
            "nNodesInfeasAux": parse_int(l0[20]),
            "gcpAvgTime": parse_float(l0[21]),
            "nsol": parse_int(l0[22]),
            "nsolHeur": parse_int(l0[23]),
            "nsolLR": parse_int(l0[24]),
            "nsolGCP": parse_int(l0[25]),
            "nsolTrivial": parse_int(l0[26]),
            "nNodesTrivial": parse_int(l0[27]),
            "nNodesInt": parse_int(l0[28]),
            "nNodesFrac": parse_int(l0[29]),
            "nNodesGcp": parse_int(l0[30]),
            "ninitSol": 0,
            "ninitDummy": 0,
            "ninit": 0,
        }

    # Intermediate format (without nNodesInt/nNodesFrac/nNodesGcp columns)
    if len(l0) >= 28:
        return {
            "instance": l0[0],
            "solver": get_solver(l0[1], subdir),
            "run": parse_int(l0[2]),
            "nvertices": parse_int(l0[3]),
            "nedges": parse_int(l0[4]),
            "nP": parse_int(l0[5]),
            "nQ": parse_int(l0[6]),
            "nvars": parse_int(l0[7]),
            "ncons": parse_int(l0[8]),
            "state": l0[9],
            "terminationReason": l0[10],
            "time": parse_float(l0[11]),
            "nodes": parse_int(l0[12]),
            "nodesLeft": parse_int(l0[13]),
            "lb": parse_float(l0[14]),
            "ub": parse_float(l0[15]),
            "gap": parse_float(l0[16]),
            "nNodesInfeas": parse_float(l0[17]),
            "nNodesInfeasPrepro": parse_int(l0[18]),
            "nNodesInfeasCheck": parse_int(l0[19]),
            "nNodesInfeasAux": parse_int(l0[20]),
            "gcpAvgTime": parse_float(l0[21]),
            "nsol": parse_int(l0[22]),
            "nsolHeur": parse_int(l0[23]),
            "nsolLR": parse_int(l0[24]),
            "nsolGCP": parse_int(l0[25]),
            "nsolTrivial": parse_int(l0[26]),
            "nNodesTrivial": parse_int(l0[27]),
            "nNodesInt": 0,
            "nNodesFrac": 0,
            "nNodesGcp": 0,
            "ninitSol": 0,
            "ninitDummy": 0,
            "ninit": 0,
        }

    # Old format (legacy merge_stats.py expectations):
    # instance,solver,run,nvertices,nedges,n,m,nvars,ncons,state,time,
    # nodes,nodesLeft,lb,ub,gap,ninfeas,ninfeasPrepro,ninfeasCheck,ninfeasAux,
    # nint,ngcp,gcpTime,nsol,nsolHeur,nsolLR,ntrivial (legacy names/order)
    if len(l0) >= 27:
        return {
            "instance": l0[0],
            "solver": get_solver(l0[1], subdir),
            "run": parse_int(l0[2]),
            "nvertices": parse_int(l0[3]),
            "nedges": parse_int(l0[4]),
            "nP": parse_int(l0[5]),
            "nQ": parse_int(l0[6]),
            "nvars": parse_int(l0[7]),
            "ncons": parse_int(l0[8]),
            "state": l0[9],
            "terminationReason": "",
            "time": parse_float(l0[10]),
            "nodes": parse_int(l0[11]),
            "nodesLeft": parse_int(l0[12]),
            "lb": parse_float(l0[13]),
            "ub": parse_float(l0[14]),
            "gap": parse_float(l0[15]),
            "nNodesInfeas": parse_float(l0[16]),
            "nNodesInfeasPrepro": parse_int(l0[17]),
            "nNodesInfeasCheck": parse_int(l0[18]),
            "nNodesInfeasAux": parse_int(l0[19]),
            "gcpAvgTime": parse_float(l0[22]),
            "nsol": parse_int(l0[23]),
            "nsolHeur": parse_int(l0[24]),
            "nsolLR": parse_int(l0[25]),
            "nsolGCP": 0,
            "nsolTrivial": parse_int(l0[26]),
            "nNodesTrivial": parse_int(l0[26]),
            "nNodesInt": parse_int(l0[20]),
            "nNodesFrac": 0,
            "nNodesGcp": parse_int(l0[21]),
            "ninitSol": 0,
            "ninitDummy": 0,
            "ninit": 0,
        }

    raise ValueError(f"Unexpected first-line stats format with {len(l0)} fields")


def parse_second_line(l1: List[str]) -> Dict[str, Any]:
    has_root_semigreedy_iters = len(l1) >= 23
    base_index = 1 if has_root_semigreedy_iters else 0

    return {
        "rootlb": parse_float(l1[0]),
        "rootub": parse_float(l1[1]),
        "rootHeurTime": parse_float(l1[2]),
        "rootSemigreedyIters": parse_int(l1[3]) if has_root_semigreedy_iters else 0,
        "rootFeasTime": parse_float(l1[3 + base_index]),
        "rootNCalls": parse_int(l1[4 + base_index]),
        "rootNCallsPool": parse_int(l1[5 + base_index]),
        "rootNCallsHeur": parse_int(l1[6 + base_index]),
        "rootNCallsMwis1": parse_int(l1[7 + base_index]),
        "rootNCallsMwis2": parse_int(l1[8 + base_index]),
        "rootNCallsExact": parse_int(l1[9 + base_index]),
        "rootNCols": parse_int(l1[10 + base_index]),
        "rootNColsPool": parse_int(l1[11 + base_index]),
        "rootNColsHeur": parse_int(l1[12 + base_index]),
        "rootNColsMwis1": parse_int(l1[13 + base_index]),
        "rootNColsMwis2": parse_int(l1[14 + base_index]),
        "rootNColsExact": parse_int(l1[15 + base_index]),
        "rootTime": parse_float(l1[16 + base_index]),
        "rootTimePool": parse_float(l1[17 + base_index]),
        "rootTimeHeur": parse_float(l1[18 + base_index]),
        "rootTimeMwis1": parse_float(l1[19 + base_index]),
        "rootTimeMwis2": parse_float(l1[20 + base_index]),
        "rootTimeExact": parse_float(l1[21 + base_index]),
    }


def parse_third_line(l2: List[str]) -> Dict[str, Any]:
    return {
        "otherNodesHeurTime": parse_float(l2[0]),
        "otherNodesFeasNCalls": parse_int(l2[1]),
        "otherNodesFeasTime": parse_float(l2[2]),
        "otherNodesNCalls": parse_int(l2[3]),
        "otherNodesNCallsPool": parse_int(l2[4]),
        "otherNodesNCallsHeur": parse_int(l2[5]),
        "otherNodesNCallsMwis1": parse_int(l2[6]),
        "otherNodesNCallsMwis2": parse_int(l2[7]),
        "otherNodesNCallsExact": parse_int(l2[8]),
        "otherNodesNCols": parse_int(l2[9]),
        "otherNodesNColsPool": parse_int(l2[10]),
        "otherNodesNColsHeur": parse_int(l2[11]),
        "otherNodesNColsMwis1": parse_int(l2[12]),
        "otherNodesNColsMwis2": parse_int(l2[13]),
        "otherNodesNColsExact": parse_int(l2[14]),
        "otherNodesTime": parse_float(l2[15]),
        "otherNodesTimePool": parse_float(l2[16]),
        "otherNodesTimeHeur": parse_float(l2[17]),
        "otherNodesTimeMwis1": parse_float(l2[18]),
        "otherNodesTimeMwis2": parse_float(l2[19]),
        "otherNodesTimeExact": parse_float(l2[20]),
    }


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: python merge_stats_v2.py <stats_path> <output_csv>")
        return 2

    stats_path = sys.argv[1]
    output_csv = sys.argv[2]

    rows: List[Dict[str, Any]] = []
    for subdir, _, files in os.walk(stats_path):
        for file in files:
            if os.path.splitext(file)[-1] != ".stat":
                continue

            stat_file = os.path.join(subdir, file)
            with open(stat_file, encoding="utf-8") as f:
                l0 = f.readline().strip().split(",")
                l1 = f.readline().strip().split(",")
                l2 = f.readline().strip().split(",")

            if len(l1) < 22 or len(l2) < 21:
                raise ValueError(f"Unexpected stats format in {stat_file}")

            row: Dict[str, Any] = {}
            row.update(parse_first_line(l0, subdir))
            row.update(parse_second_line(l1))
            row.update(parse_third_line(l2))
            rows.append(row)

    df = pd.DataFrame(rows)

    preferred_order = [
        "instance",
        "solver",
        "run",
        "nvertices",
        "nedges",
        "nP",
        "nQ",
        "nvars",
        "ncons",
        "state",
        "terminationReason",
        "time",
        "nodes",
        "nodesLeft",
        "lb",
        "ub",
        "gap",
        "nNodesInt",
        "nNodesFrac",
        "nNodesGcp",
        "nNodesTrivial",
        "nNodesInfeas",
        "nNodesInfeasPrepro",
        "nNodesInfeasCheck",
        "nNodesInfeasAux",
        "gcpAvgTime",
        "nsol",
        "nsolHeur",
        "nsolLR",
        "nsolGCP",
        "nsolTrivial",
        "ninitSol",
        "ninitDummy",
        "ninit",
        "rootlb",
        "rootub",
        "rootHeurTime",
        "rootSemigreedyIters",
        "rootFeasTime",
        "rootNCalls",
        "rootNCallsPool",
        "rootNCallsHeur",
        "rootNCallsMwis1",
        "rootNCallsMwis2",
        "rootNCallsExact",
        "rootNCols",
        "rootNColsPool",
        "rootNColsHeur",
        "rootNColsMwis1",
        "rootNColsMwis2",
        "rootNColsExact",
        "rootTime",
        "rootTimePool",
        "rootTimeHeur",
        "rootTimeMwis1",
        "rootTimeMwis2",
        "rootTimeExact",
        "otherNodesHeurTime",
        "otherNodesFeasNCalls",
        "otherNodesFeasTime",
        "otherNodesNCalls",
        "otherNodesNCallsPool",
        "otherNodesNCallsHeur",
        "otherNodesNCallsMwis1",
        "otherNodesNCallsMwis2",
        "otherNodesNCallsExact",
        "otherNodesNCols",
        "otherNodesNColsPool",
        "otherNodesNColsHeur",
        "otherNodesNColsMwis1",
        "otherNodesNColsMwis2",
        "otherNodesNColsExact",
        "otherNodesTime",
        "otherNodesTimePool",
        "otherNodesTimeHeur",
        "otherNodesTimeMwis1",
        "otherNodesTimeMwis2",
        "otherNodesTimeExact",
    ]

    existing = [c for c in preferred_order if c in df.columns]
    extra = [c for c in df.columns if c not in existing]
    df = df[existing + extra]

    df.to_csv(output_csv, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
