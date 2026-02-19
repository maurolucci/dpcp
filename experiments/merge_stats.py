import os
import sys
import pandas as pd


def get_solver(name, subdir):
    parts = subdir.split(os.path.sep)[2:-1]
    return name + "_" + "_".join(parts)


def main() -> int:
    statsPath = sys.argv[1]
    statsName = sys.argv[2]

    columns = [
        "instance",
        "solver",
        "run",
        "nvertices",
        "nedges",
        "nA",
        "nB",
        "nvars",
        "ncons",
        "state",
        "time",
        "nodes",
        "nodesLeft",
        "lb",
        "ub",
        "gap",
        "ninfeas",
        "ninfeasPrepro",
        "ninfeasCheck",
        "ninfeasAux",
        "nint",
        "ngcp",
        "gcpTime",
        "nsol",
        "nsolHeur",
        "nsolLR",
        "ntrivial",
        "bestTime",
        "bestIter",
        "rootlb",
        "rootub",
        "rootHeurTime",
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

    data = []
    for subdir, dirs, files in os.walk(statsPath):
        for file in files:
            if os.path.splitext(file)[-1] != ".stat":
                continue
            with open(os.path.join(subdir, file)) as f:

                # first line
                l0 = f.readline().strip().split(",")
                entry = [
                    l0[0],  # instance
                    get_solver(l0[1], subdir),  # solver
                    int(l0[2]),  # run
                    int(l0[3]),  # nvertices
                    int(l0[4]),  # nedges
                    int(l0[5]),  # nA
                    int(l0[6]),  # nB
                    int(l0[7]),  # nvars
                    int(l0[8]),  # ncons
                    l0[9],  # state
                    float(l0[10]),  # time
                    int(l0[11]),  # nodes
                    int(l0[12]),  # nodesLeft
                    float(l0[13]),  # lb
                    float(l0[14]),  # ub
                    float(l0[15]),  # gap
                    float(l0[16]),  # ninfeas
                    int(l0[17]),  # ninfeasPrepro
                    int(l0[18]),  # ninfeasCheck
                    int(l0[19]),  # ninfeasAux
                    int(l0[20]),  # nint
                    int(l0[21]),  # ngcp
                    float(l0[22]),  # gcpTime
                    int(l0[23]),  # nsol
                    int(l0[24]),  # nsolHeur
                    int(l0[25]),  # nsolLR
                    int(l0[26]),  # ntrivial
                    float(l0[27]),  # bestTime
                    int(l0[28]),
                ]  # bestIter

                # second line
                l1 = f.readline().strip().split(",")
                entry.extend(
                    [
                        float(l1[0]),  # rootlb
                        float(l1[1]),  # rootub
                        float(l1[2]),  # rootHeurTime
                        float(l1[3]),  # rootFeasTime
                        int(l1[4]),  # rootNCalls
                        int(l1[5]),  # rootNCallsPool
                        int(l1[6]),  # rootNCallsHeur
                        int(l1[7]),  # rootNCallsMwis1
                        int(l1[8]),  # rootNCallsMwis2
                        int(l1[9]),  # rootNCallsExact
                        int(l1[10]),  # rootNCols
                        int(l1[11]),  # rootNColsPool
                        int(l1[12]),  # rootNColsHeur
                        int(l1[13]),  # rootNColsMwis1
                        int(l1[14]),  # rootNColsMwis2
                        int(l1[15]),  # rootNColsExact
                        float(l1[16]),  # rootTime
                        float(l1[17]),  # rootTimePool
                        float(l1[18]),  # rootTimeHeur
                        float(l1[19]),  # rootTimeMwis1
                        float(l1[20]),  # rootTimeMwis2
                        float(l1[21]),  # rootTimeExact
                    ]
                )

                # third line
                l2 = f.readline().strip().split(",")
                entry.extend(
                    [
                        float(l2[0]),  # otherNodesHeurTime
                        int(l2[1]),  # otherNodesFeasNCalls
                        float(l2[2]),  # otherNodesFeasTime
                        int(l2[3]),  # otherNodesNCalls
                        int(l2[4]),  # otherNodesNCallsPool
                        int(l2[5]),  # otherNodesNCallsHeur
                        int(l2[6]),  # otherNodesNCallsMwis1
                        int(l2[7]),  # otherNodesNCallsMwis2
                        int(l2[8]),  # otherNodesNCallsExact
                        int(l2[9]),  # otherNodesNCols
                        int(l2[10]),  # otherNodesNColsPool
                        int(l2[11]),  # otherNodesNColsHeur
                        int(l2[12]),  # otherNodesNColsMwis1
                        int(l2[13]),  # otherNodesNColsMwis2
                        int(l2[14]),  # otherNodesNColsExact
                        float(l2[15]),  # otherNodesTime
                        float(l2[16]),  # otherNodesTimePool
                        float(l2[17]),  # otherNodesTimeHeur
                        float(l2[18]),  # otherNodesTimeMwis1
                        float(l2[19]),  # otherNodesTimeMwis2
                        float(l2[20]),  # otherNodesTimeExact
                    ]
                )

                data.append(entry)

    df = pd.DataFrame(data, columns=columns)
    df.to_csv(statsName, index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
