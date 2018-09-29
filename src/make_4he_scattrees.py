
import treeconvutils


runs = range(272, 316 + 1)


if __name__ == "__main__":
    treeconvutils.convert_runs_range(runs, max_processes=4)
