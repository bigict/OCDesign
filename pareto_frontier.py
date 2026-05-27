"""
Pareto frontier select tool
Function: find the parete optimal solution according to the attributes
Input: csv file including the name column and all attribute columns
Output: the pareto frontier name
"""

import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Set
import argparse
import sys


def load_data(csv_file: str) -> pd.DataFrame:
    """
    load the date in the csv file
    """
    try:
        df = pd.read_csv(csv_file)
        print(f"successfully load date from: {csv_file}")
        print(f"data shape {df.shape}")
        print(f"attribute columns {list(df.columns)}")
        return df
    except Exception as e:
        print(f"failed to load data: {e}")
        sys.exit(1)


def parse_sort_directions(sort_str: str, columns: List[str]) -> Dict[str, int]:
    """
    parse the sort method
    format: "attribute1:max,attibute2:min,..."
    max: maximization is better
    min: minimization is better
    """
    directions = {}
    
    # defalut sort method:max
    if not sort_str:
        print("without specify the sort method, using defaule method:max")
        for col in columns:
            directions[col] = 1  # 1:max，-1:min
        return directions
    
    # parse the sort direction
    for item in sort_str.split(','):
        if ':' not in item:
            print(f"Warning: format error '{item}'，continue")
            continue
        
        col, direction = item.split(':')
        col = col.strip()
        direction = direction.strip().lower()
        
        if col not in columns:
            print(f"Warining: attribute '{col}' not in date，continue")
            continue
        
        if direction == 'max':
            directions[col] = 1
        elif direction == 'min':
            directions[col] = -1
        else:
            print(f"Warining: unknow sort direciton '{direction}'，using defalut:max")
            directions[col] = 1
    
    # set the default sort direction for unspecify attrubutes
    for col in columns:
        if col not in directions:
            directions[col] = 1  #  default:maximization is better
            print(f"attribute '{col}' using defaule direction: max")
    
    return directions


def normalize_data(df: pd.DataFrame, directions: Dict[str, int]) -> pd.DataFrame:
    """
    normalize data: change the data direction, if the sort method is min, the negtive the data  to make the data change to maximization is better
    """
    df_normalized = df.copy()
    
    for col, direction in directions.items():
        if col in df_normalized.columns:
            if direction == -1: 
                df_normalized[col] = -df_normalized[col]
    
    return df_normalized


def is_dominated(point_a: np.ndarray, point_b: np.ndarray) -> bool:
    """
    check whether a is dominated by a 
    b dominates a: b is not worse than a in all attribute, and better than a in at least 1 attibute
    """
    # check whether b is worse than a
    all_not_worse = np.all(point_b >= point_a)
    # check b is better than a in at least 1 attribute
    any_better = np.any(point_b > point_a)
    
    return all_not_worse and any_better


def find_pareto_frontier(df: pd.DataFrame, attribute_cols: List[str]) -> List[int]:
    """
    find the pareto frontier indices
    """
    n_points = len(df)
    pareto_indices = []
    
    # get the numpy array of the attribute
    data = df[attribute_cols].values
    
    for i in range(n_points):
        dominated = False
        
        for j in range(n_points):
            if i == j:
                continue
            
            # wheter the ponit i is dominate by another
            if is_dominated(data[i], data[j]):
                dominated = True
                break
        
        # if i is not dominated by another, i is the best point
        if not dominated:
            pareto_indices.append(i)
    
    return pareto_indices


def main():
    parser = argparse.ArgumentParser(description='Tool to select the parete frontier')
    parser.add_argument('csv_file', help='csv file path')
    parser.add_argument('--sort', '-s', help='sort method,eg: "property1:max,property2:min,..."')
    parser.add_argument('--name-col', '-n', default='Name', 
                       help='the name column, default is Name')
    parser.add_argument('--output', '-o', help='output the result to a csv file(not required)')
    
    args = parser.parse_args()
    
    # 1. load date from csv file
    df = load_data(args.csv_file)
    
    # 2. get the attribute columns that you specify through --sort or -s
    attribute_cols = [col for col in df.columns if col != args.name_col and col in args.sort]
    if args.name_col in df.columns:name_col = args.name_col
    else:name_col = df.columns[0]
    print(f"Warning:Name Column '{args.name_col}' does not exist，using the first column {name_col}")
    print(f"Attribute columns: {attribute_cols}")
    
    # 3. attribute sort method:max/min
    directions = parse_sort_directions(args.sort, attribute_cols)
    
    # 4. normalize the data: add negetive(-) or not 
    df_normalized = normalize_data(df, directions)
    
    # 5. find pareto frontier
    pareto_indices = find_pareto_frontier(df_normalized, attribute_cols)
    
    # 6. output the result 
    print("\n" + "="*50)
    print("Pareto frontier")
    print("="*50)
    
    if not pareto_indices:
        print("Cat find parete frontier")
        return
    
    print(f"Find {len(pareto_indices)} Pareto frontier:")
    print("-"*50)
    
    results = []
    for idx in pareto_indices:
        name = df.iloc[idx][name_col]
        attributes = {col: df.iloc[idx][col] for col in attribute_cols}
        results.append((name, attributes))
        
        print(f"Name: {name}")
        for col in attribute_cols:
            direction = "max" if directions[col] == 1 else "min"
            print(f"  {col}: {df.iloc[idx][col]} ({direction})")
        print("-"*30)
    
    # 7. save output result if specify the output file
    if args.output:
        # get the dateframe
        result_df = df.iloc[pareto_indices].copy()
        result_df['is_pareto'] = True
        
        # saveSV
        result_df.to_csv(args.output, index=False)
        print(f"\n All result are saved to: {args.output}")
    
    return results


if __name__ == "__main__":
    main()
