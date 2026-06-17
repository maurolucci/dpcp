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
    # instance,solver,run,nvertices,nedges,nP,nQ,nvars,ncons,state,
    # terminationReason,time,nodes,nodesLeft,lb,ub,gap,initialHeurValue,
    # initialHeurTime,initialSemigreedyIters,nNodesInt,nNodesFrac,nNodesGcp,
    # nNodesTrivial,nNodesInfeas,nNodesInfeasPrepro,nNodesInfeasCheck,
    # nNodesInfeasAux,gcpAvgTime,nsol,nsolHeur,nsolLR,nsolGCP,nsolTrivial,
    # ninitSol,ninitDummy,ninit
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
        "initialHeurValue": parse_float(l0[17]),
        "initialHeurTime": parse_float(l0[18]),
        "initialSemigreedyIters": parse_int(l0[19]),
        "nNodesInt": parse_int(l0[20]),
        "nNodesFrac": parse_int(l0[21]),
        "nNodesGcp": parse_int(l0[22]),
        "nNodesTrivial": parse_int(l0[23]),
        "nNodesInfeas": parse_float(l0[24]),
        "nNodesInfeasPrepro": parse_int(l0[25]),
        "nNodesInfeasCheck": parse_int(l0[26]),
        "nNodesInfeasAux": parse_int(l0[27]),
        "gcpAvgTime": parse_float(l0[28]),
        "nsol": parse_int(l0[29]),
        "nsolHeur": parse_int(l0[30]),
        "nsolLR": parse_int(l0[31]),
        "nsolGCP": parse_int(l0[32]),
        "nsolTrivial": parse_int(l0[33]),
        "ninitSol": parse_int(l0[34]),
        "ninitDummy": parse_int(l0[35]),
        "ninit": parse_int(l0[36]),
    }


def parse_second_line(l1: List[str]) -> Dict[str, Any]:
    # Current format (29 fields): rootNVertices,rootNEdges,rootNP,rootNQ,
    # rootlb,rootub,rootHeurTime,rootFeasTime,(empty),rootCgTime,
    # rootNCalls,rootNCallsPool,rootNCallsHeur,rootNCallsMwis1,rootNCallsMwis2,
    # rootNCallsExact,rootNCols,rootNColsPool,rootNColsHeur,rootNColsMwis1,
    # rootNColsMwis2,rootNColsExact,rootTime,rootTimePool,rootTimeHeur,
    # rootTimeMwis1,rootTimeMwis2,rootTimeExact
    return {
        "rootNVertices": parse_int(l1[0]),
        "rootNEdges": parse_int(l1[1]),
        "rootNP": parse_int(l1[2]),
        "rootNQ": parse_int(l1[3]),
        "rootlb": parse_float(l1[4]),
        "rootub": parse_float(l1[5]),
        "rootHeurTime": parse_float(l1[6]),
        "rootFeasTime": parse_float(l1[7]),
        "rootCgTime": parse_float(l1[8]), 
        "rootNCalls": parse_int(l1[9]),
        "rootNCallsPool": parse_int(l1[10]),
        "rootNCallsHeur": parse_int(l1[11]),
        "rootNCallsMwis1": parse_int(l1[12]),
        "rootNCallsMwis2": parse_int(l1[13]),
        "rootNCallsExact": parse_int(l1[14]),
        "rootNCols": parse_int(l1[15]),
        "rootNColsPool": parse_int(l1[16]),
        "rootNColsHeur": parse_int(l1[17]),
        "rootNColsMwis1": parse_int(l1[18]),
        "rootNColsMwis2": parse_int(l1[19]),
        "rootNColsExact": parse_int(l1[20]),
        "rootTime": parse_float(l1[21]),
        "rootTimePool": parse_float(l1[22]),
        "rootTimeHeur": parse_float(l1[23]),
        "rootTimeMwis1": parse_float(l1[24]),
        "rootTimeMwis2": parse_float(l1[25]),
        "rootTimeExact": parse_float(l1[26]),
    }


def parse_third_line(l2: List[str]) -> Dict[str, Any]:
    # Current format (21 fields): otherNodesHeurTime,otherNodesFeasNCalls,
    # otherNodesFeasTime,otherNodesNCalls,otherNodesNCallsPool,
    # otherNodesNCallsHeur,otherNodesNCallsMwis1,otherNodesNCallsMwis2,
    # otherNodesNCallsExact,otherNodesNCols,otherNodesNColsPool,
    # otherNodesNColsHeur,otherNodesNColsMwis1,otherNodesNColsMwis2,
    # otherNodesNColsExact,otherNodesTime,otherNodesTimePool,otherNodesTimeHeur,
    # otherNodesTimeMwis1,otherNodesTimeMwis2,otherNodesTimeExact
    return {
        "otherNodesHeurTime": parse_float(l2[0]),
        "otherNodesFeasNCalls": parse_int(l2[1]),
        "otherNodesFeasTime": parse_float(l2[2]),
        "otherNodesNCalls": parse_float(l2[3]),
        "otherNodesNCallsPool": parse_float(l2[4]),
        "otherNodesNCallsHeur": parse_float(l2[5]),
        "otherNodesNCallsMwis1": parse_float(l2[6]),
        "otherNodesNCallsMwis2": parse_float(l2[7]),
        "otherNodesNCallsExact": parse_float(l2[8]),
        "otherNodesNCols": parse_float(l2[9]),
        "otherNodesNColsPool": parse_float(l2[10]),
        "otherNodesNColsHeur": parse_float(l2[11]),
        "otherNodesNColsMwis1": parse_float(l2[12]),
        "otherNodesNColsMwis2": parse_float(l2[13]),
        "otherNodesNColsExact": parse_float(l2[14]),
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
        "initialHeurValue",
        "initialHeurTime",
        "initialSemigreedyIters",
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
        "rootNVertices",
        "rootNEdges",
        "rootNP",
        "rootNQ",
        "rootlb",
        "rootub",
        "rootHeurTime",
        "rootFeasTime",
        "rootCgTime",
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
