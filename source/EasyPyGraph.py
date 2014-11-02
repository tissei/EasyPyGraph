#!/usr/bin/env python
# -*- coding: ascii -*-

"""
PyGraph
~~~~~~~~~~~~~

[TODO -> Brief introduction to the PyGraph]

[TODO -> Brief introduction to the DiGraph Class]

WARNING! This is a prototype.

Copyright 2014 José Carlos S.A. Tissei and Lucas de Oliveira Teixeira

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   
"""

#imports

import heapq


class DiGraph:
    def __init__(self):
        self.adjList = {}

    # TODO -> Improve printing method!
    def __str__(self):
        return str( self.adjList )
    
    """
    Add an vertex to the adjList
    
    verify if the vertex dont exist already, if not adds him to the adjList
    
    @param vertex integer the vertex to be added
    @return       void
    """
    
    def addVertex(self, vertex):
        if vertex not in self.adjList:
            self.adjList[vertex] = {}

    """
    Add an edge between two vertexes
    
    adds the vertex to the adjList if they dont exist and add the weight
    between then in the adjList.
    
    @see addVertex method    
    @param vertex   integer the starting vertex
    @param adjacent integer the adjacent vertex
    @return         void
    """
            
    def addEdge(self, vertex, adjacent, weight = 1):
        self.addVertex(vertex)
        self.addVertex(adjacent)
        self.adjList[vertex][adjacent] = weight
        
    """
    Remove a vertex and all his edges
    
    remove the given vertex from the adjList and iterates trough all other
    vertices adjList to remove any edges connecting to him.
    
    @param vertex integer the starting vertex
    @return       void
    """

    def removeVertex(self, vertex):
        if vertex in self.adjList:
            del self.adjList[vertex]
            for adjacent in self.adjList:
                if vertex in self.adjList[adjacent]:
                    del self.adjList[adjacent][vertex]
    
    """
    Remove a edge between two vertexes
    
    delete the adjList index between the starting vertex and its adjacent
    
    @param vertex integer the starting vertex
    @param vertex integer the starting vertex adjacent
    @return       void
    """

    def removeEdge(self, vertex, adjacent):
        if vertex in self.adjList:
            if adjacent in self.adjList[vertex]:
                del self.adjList[vertex][adjacent]
                
    """            
    Verifies if a given vertex exists
    
    verify if the vertex is in the adjList
    
    @param integer vertex the vertex to be searched
    @return
    """
                
    def hasVertex(self, vertex):
        return vertex in self.adjList
    """
    Verifies if an edge exists between two vertex
    
    verify if the vertex is in the adjList and if it's adjacent is in the 
    vertex adjList.
    
    @param vertex   integer the starting vertex
    @param adjacent integer the adjacent vertex from the starting vertex
    @return         boolean 
    """

    def hasEdge(self, vertex, adjacent):
        return vertex in self.adjList and adjacent in self.adjList[vertex]
    
    """
    Returns the output degree of a vertex
    
    returns the lenght of the adjacent list from the given vertex
    
    @param vertex integer the vertex wich the degree will be calculated
    @return integer if the vertex exists and None if it dont
    """
    
    def outputDegree(self, vertex):
        if vertex in self.adjList:
            return len(self.adjList[vertex])
        return None
    
    """
    Returns the input degree of a given vertex
    
    iterate trough all vertexes in the adjList and verify if any of 
    them is adjacent to the given vertex.
    
    @param vertex integer the vertex wich the degree will be calculated
    @return integer if the vertex exists and None if it dont
    """
    
    def inputDegree(self, vertex):
        if vertex in self.adjList:
            degree = 0
            for adjacent in self.adjList:
                if vertex in self.adjList[adjacent]:
                    degree += 1
            return degree
        return None
    
    """
    Breadth-first search (BFS) is a strategy for searching in a graph when 
    search is limited to essentially two operations: (a) visit and inspect a 
    node of a graph; (b) gain access to visit the nodes that neighbor the 
    currently visited node.
    
    the BFS begins at a root node and inspects all the neighboring nodes. 
    Then for each of those neighbor nodes in turn, it inspects their neighbor 
    nodes which were unvisited, and so on.
    
    @param vertex root vertex
    @return list list with the search order
    """
    
    def breadthFirstSearch(self, start):
        visited = [start]
        queue = [start]
        while queue:
            t = queue.pop(0)
            if t == start and visited[0] != t: 
                return visited
            for adjacent in self.adjList[t]:
                if adjacent not in visited:
                    visited.append(adjacent)
                    queue.append(adjacent)
        return visited
    
    """
    Depth-first search (DFS) is an algorithm for traversing or searching 
    graphs. Starts with a root node and explores as far as possible along 
    each branch before backtracking.    
    
    iterates trough the starting vertex adjacents, select on of then and
    keep going trought it's adjacent until theres a vertex with no adjacents
    or all of its adjacent have already being visited, then comes back to the
    starting vertex, select another adjacent and reapeat the process until all
    vertexes are visited.
    
    @param vertex starting vertex
    @return list list with the search order
    """
    
    def depthFirstSearch(self, start):
        visited = []
        queue = [start]
        while lista:
            t = queue.pop()
            if t not in visited:
                visited.append(t)
                for adjacent in self.adjList[t]:
                    queue.append(adjacent)
        return visited

    """
    Verifies if a vertex can be reached from another vertex
    
    utilizes the depth search to verify if there's any way to reach the
    end vertex from the start vertex.
    
    @param vertex integer the starting vertex
    @param vertex integer the vertex to be reached
    @return       boolean 
    """
    
    def areConected(self, start, end):
         end in self.depthSearch(start): 
    
    """
    Verifies if a vertex can be reached back starting from himself
    
    utilizes the depth search to verify if there's any way to reach back
    the vertex from any possible conection from his adjacents.
    
    @param vertex integer the starting vertex
    @return       boolean 
    """
    
    def hasCicle(self, vertex):
        return vertex in self.depthSearch(vertex)
    
    """
    Verifies if there is no vertex that can be reached back starting from himself
    
    iterate trough adjList to test all vertexes and verifies if there is any way
    that they can be reached back starting from himselfs.
    
    @param vertex integer the starting vertex
    @return       boolean 
    """
    
    def isAcyclic(self):
        for vertex in self.adjList:
            for adjacent in self.adjList[vertex]:
                if vertex in self.adjList[adjacent]:
                    return False
        return True
    
    """
    Topological ordering of a graph is a linear ordering of its 
    vertices such that for every directed edge uv from vertex u to vertex v, 
    u comes before v in the ordering.
    
    search for all vertexes with input degree equal to 0 and put them in queue,
    iterate trough the queue breaking the conection between the vertexes and 
    its adjacents, if after this any adjacent has an input degree equal to 0
    the vertex is add to the queue, this method is repeated until the queue is
    empty.
    
    @return list list of ordered vertexes 
    """
    
    def topologicalSorting(self):
        adjListCopy = dict(self.adjList)
        outOnly = []
        count = 0
        for key in self.adjList:
            if self.inputDegree(key) == 0:
                outOnly.append(key)
        ordered = {}
        while len(outOnly) > 0:
            v = outOnly.pop(0)
            ordered[v] = count
            count += 1
            vList = adjListCopy[v].copy()
            for adjacent in vList:
                del adjListCopy[v][adjacent]
                if self.inputDegree(adjacent) == 0:
                       outOnly.append(adjacent)
        for key in adjListCopy:
            if len(adjListCopy[key]) > 0: return None
        else:
            return ordered
        
    """
    A bipartite graph is a graph whose vertices can be divided into two 
    disjoint sets U and V (that is, U and V are each independent sets) such 
    that every edge connects a vertex in U to one in V.
    
    The aproach choosen is the graph coloring were a graph is bipartite if it
    can be colored using only 2 colors.
    The algoritm iterates trough all vertexes and for each on of them test
    the attribution of two colors by checking if any of his adjacents is the
    same color, then iterates again trough all vertex checking if any of them
    are connected to the selected vertex has the same color, if none of those
    two conditions is true, then assing the color to the vertex and break the
    iteration.
    
    @return boolean
    """
        
    def isBipartite(self):
        colors = [0, 1]
        coloredGraph = {}
        for key in self.adjList: coloredGraph[key] = None
        for vertex in coloredGraph:
            if coloredGraph[vertex] == None:
                for color in colors:
                    checkAdj = True
                    for adjacent in self.adjList[vertex]:
                        if color == coloredGraph[adjacent]:
                            checkAdj = False
                            break
                    for isConected in self.adjList:
                        if vertex != isConected:
                            if vertex in self.adjList[isConected]:
                                if color == coloredGraph[isConected]:
                                    checkAdj = False
                                    break
                    if checkAdj:
                        coloredGraph[vertice] = color
                        break
        if None in coloredGraph: return False
        return True
        
    """
    Dijkstra's algorithm, is a graph search algorithm that solves the 
    single-source shortest path problem for a graph with non-negative edge 
    path costs, producing a shortest path tree.
    
    iterates trought all vertexes, put then on the ways list and verify if
    their weight is the lowest or if there is no weight saved, then put the
    adjacent vertexes is the ways list and repeat the process to each one,
    until all paths weights are calculated.
    

    @param  int  starting vertex
    @return list list with the shortest paths  and weights from the starting 
                 vertex to all vertexes.
    """
        
    def Dijkstra(self, v):
        dist = dict.fromkeys(self.adjList)
        ways = []
        heapq.heappush(ways, ( 0 , [v] ))
        while ways:
            ( weight , way ) = heapq.heappop(ways)
            last = way[ len(way) - 1 ]
            dist[last] = weight
            for neighbour in self.adjList[last]:
                if dist[neighbour] == None:
                    heapq.heappush(ways, ( weight + self.adjList[last][neighbour] , way + [neighbour] ) )
        return dist
    
        
    """
    The Bellman–Ford algorithm is an algorithm that computes shortest paths 
    from a single source vertex to all of the other vertices. It is slower 
    than Dijkstra's algorithm for the same problem, but more versatile, as 
    it is capable of handling graphs in which some of 
    the edge weights are negative numbers.
    
    iterates trought all vertexes, put then on the ways list and verify if
    their weight is the lowest or if there is no weight saved, then put the
    adjacent vertexes is the ways list and repeat the process to each one,
    until all paths weights are calculated.
    iterates trough all vertexes and it's adjcents and verify if there's any
    negative cicle.
    

    @param  int  starting vertex
    @return list list with the shortest paths  and weights from the starting 
                 vertex to all vertexes.
    """
            
    def BellmanFord(self, v):
        dist = dict.fromkeys(self.adjList)
        ways = []
        heapq.heappush(ways, ( 0 , [v] ))
        while ways:
            ( weight , way ) = heapq.heappop(ways)
            last = way[ len(way) - 1 ]
            dist[last] = weight
            for neighbour in self.adjList[last]:
                if dist[neighbour] == None:
                    heapq.heappush(ways, ( weight + self.adjList[last][neighbour] , way + [neighbour] ) )
        for vertex in self.adjList:
            for neighbour in self.adjList[vertex]:
                if dist[vertex] != None and dist[neighbour] != None:
                    if dist[vertex] > dist[neighbour] + self.adjList[vertex][neighbour]:
                        print 'negative cicle between ' + (str) (vertex) + ' and ' + (str) (neighbour) 
                        return None
        return dist

    """
    The Floyd–Warshall algorithm is a graph analysis algorithm for finding 
    shortest paths in a weighted graph with positive or negative edge weights 
    (but with no negative cycles). A single execution of the algorithm will 
    find the lengths (summed weights) of the shortest paths between all pairs 
    of vertices, though it does not return details of the paths themselves.
    
    Uses 3 iterations to make all possible comparisons between the vertexes
    and then compare 3 diferent combinations to get the weight, keep iterating
    until all the paths are found, the weights are put in a matrix and then
    the matrix is returned.
    
    @return [[]] matrix with the shortest paths from all vertexes to all vertexes.
    """
    
    def FloydWarshal(self):
        dist = self.adjList
        for k in dist:
            for i in dist:
                for j in dist:
                        if j not in dist[i]:
                            if k in dist[i] and j in dist[k]:
                                if i != j:
                                    dist[i][j] = dist[i][k] + dist[k][j]
                                else:
                                    dist[i][j] = 0
                        elif j in dist[i] and k in dist[i] and j in dist[k]:
                            if dist[i][j] > dist[i][k] + dist[k][j]:
                                dist[i][j] = dist[i][k] + dist[k][j] 
        return dist

digraph = DiGraph()



digraph.addVertex(10)



digraph.addEdge(10,2)
digraph.addEdge(2,5)
digraph.addEdge(3,2)
digraph.addEdge(5,7)
digraph.addEdge(2,10)
digraph.addEdge(2,4)
digraph.addEdge(4,7)
digraph.addEdge(4,8)

print digraph

print digraph.lengthSearch(2)

