
import treeconvutils


runs = range(324, 350 + 1)

if __name__ == "__main__":
    treeconvutils.convert_runs_range(runs, max_processes=4)
