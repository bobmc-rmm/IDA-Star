# IDA-Star
Interative Deepening Astar Graph Search

This algorithm is relatively easy to implement since it does not require abstractions such as Lists and Queues. The graph data structues can be implemented as automatic instance variables of the recursive Search function. A segment of the search space is searched by DFS so that the required memory increases by one node per step. Many more nodes are created and released during recursive backtracking to the start node. Since all the nodes upto the current depth is visited, run time increases exponentially.

This is inspired by  *  https://algorithmsinsight.wordpress.com/graph-theory-2/ida-star-algorithm-in-general/
 * see https://en.wikipedia.org/wiki/Iterative_deepening_A*

The included example puzzles are not random. I created them by starting with the goal by moving the tiles to create puzzles of increasing difficulty.  I discovered failure when the depth parameter g=52. This was with the manhattan taxi heuristic.  Next I added Linear Conflicts by row. This reduced the number of nodes generated but was inconsistent. Finally, I added sorting of the successor nodes such that the nodes with lowest heuristic value were visited first. This enabled the solution of a g=53 puzzle.
