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
        "n",
        "m",
        "state",
        "terminationReason",
        "value",
        "totalTime",
        "totalIters",
        "bestTime",
        "bestIter",
    ]

    data = []
    for subdir, dirs, files in os.walk(statsPath):
        for file in files:
            if os.path.splitext(file)[-1] != ".stat":
                continue
            with open(os.path.join(subdir, file)) as f:

                # first line
                l0 = f.readline().strip().split(",")

                # Legacy format (13 fields):
                # instance,solver,run,nvertices,nedges,n,m,state,value,totalTime,totalIters,bestTime,bestIter
                # Current format (14 fields):
                # instance,solver,run,nvertices,nedges,n,m,state,terminationReason,value,totalTime,totalIters,bestTime,bestIter
                if len(l0) >= 14:
                    termination_reason = l0[8]
                    value = l0[9]
                    total_time = l0[10]
                    total_iters = l0[11]
                    best_time = l0[12]
                    best_iter = l0[13]
                elif len(l0) >= 13:
                    termination_reason = ""
                    value = l0[8]
                    total_time = l0[9]
                    total_iters = l0[10]
                    best_time = l0[11]
                    best_iter = l0[12]
                elif len(l0) >= 12:
                    termination_reason = ""
                    value = l0[8]
                    total_time = l0[9]
                    total_iters = l0[10]
                    best_time = l0[11]
                    best_iter = pd.NA
                else:
                    continue

                entry = [
                    l0[0],  # instance
                    get_solver(l0[1], subdir),  # solver
                    int(l0[2]),  # run
                    int(l0[3]),  # nvertices
                    int(l0[4]),  # nedges
                    int(l0[5]),  # n
                    int(l0[6]),  # m
                    l0[7],  # state
                    termination_reason,
                    value,
                    total_time,
                    total_iters,
                    best_time,
                    best_iter,
                ]

                data.append(entry)

    df = pd.DataFrame(data, columns=columns)
    df.to_csv(statsName, index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
