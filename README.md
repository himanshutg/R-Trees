# R-Trees in Spatial Searches

This repository contains a project implementation related to R-Trees, a data structure used for spatial searches and nearest neighbor queries. The project was completed as part of a team assignment in the course MA 2233 Data Structures and Applications Lab. The repository includes a C++ implementation of R-Trees and various methods for querying nearest neighbors in multi-dimensional spaces.

## Files in the Repository

1. **`code.cpp`**: 
   - This is the main implementation file where the R-Tree data structure is coded.
   - The code includes functionalities for creating an R-Tree, inserting points, splitting nodes, and querying nearest neighbors using different distance metrics (Euclidean, Manhattan, L1).
   - The code also evaluates the performance and effectiveness of the R-Tree structure when handling high-dimensional data.

2. **`R-Trees.pptx`**:
   - A presentation that provides a detailed explanation of R-Trees, including how they are constructed, how they function, and their application in spatial searches.
   - The presentation also covers theoretical concepts like Minimum Bounding Rectangles (MBR) and the pruning techniques used in nearest neighbor search.

3. **`MA 2233 Team Projects-R-Trees.pdf`**:
   - A problem statement document outlining the project requirements, challenges, and goals.
     
4. **`output.txt`**:
   - This file contains the output generated by running the code, including the results of nearest neighbor queries, the performance metrics, and any observations made during the execution.

## Project Overview

R-Trees are tree data structures used for indexing multi-dimensional information such as geographical coordinates. This project focuses on:

- **Creating 200 random points** within a 10-dimensional space.
- **Constructing an R-Tree** with different depths and testing its efficiency in finding the nearest neighbor to a given query point.
- **Evaluating different distance metrics** (Euclidean, Manhattan, L1) to determine which one provides the most accurate nearest neighbor results.
- **Analyzing the impact of dimensionality** by comparing results in 10 dimensions vs. 20 dimensions.
- **Investigating the computational complexity** associated with inserting, deleting nodes, and querying in the R-Tree.

## Real-Life Applications

- **Geographical searches**: Finding nearby restaurants, gas stations, or any point of interest based on a user’s location.
- **Military applications**: Detecting enemy submarines within a certain range.

## Challenges Faced

- **Implementation**: Deciding on the structure of nodes, managing leaf and non-leaf nodes, and implementing efficient splitting algorithms.
- **Efficiency**: Ensuring the R-Tree structure remains optimal in terms of space utilization and computational performance.

## Usage

To run the code, compile `code.cpp` using a C++ compiler and execute the binary. The results will be stored in `output.txt`.

## Credits

- **Shubham Vishwakarma**
- **Kaustubh Dandegaonkar**
- **Himanshu Jindal**

The project draws on key academic resources, including papers by Antonin Guttman and research by the University of Maryland.

---
