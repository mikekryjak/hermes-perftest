#!/usr/bin/env python3

import numpy as np
import pandas as pd
import argparse


def run_benchmark(path):

    df = pd.read_csv(path)

    print(df)

#------------------------------------------------------------
# PARSER
#------------------------------------------------------------

if __name__ == "__main__":
    
    # Define arguments
    parser = argparse.ArgumentParser(description = "Case monitor")
    parser.add_argument("path", type=str, help = "Path to results table")
    
    # Extract arguments and call function
    args = parser.parse_args()
    run_benchmark(args.path)
