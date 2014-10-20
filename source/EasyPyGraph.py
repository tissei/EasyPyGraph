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
    
    def addVertex(self, vertex):
        if vertex not in self.adjList:
            self.adjList[vertex] = {}

    def addEdge(self, vertex, adjacent, weight = 1):
        self.addVertex(vertex)
        self.addVertex(adjacent)
        self.adjList[vertex][adjacent] = weight

    def removeVertex(self, vertex):
        if vertex in self.adjList:
            del self.adjList[vertex]
            for adjacent in self.adjList:
                if vertex in self.adjList[adjacent]:
                    del self.adjList[adjacent][vertex]

    def removeEdge(self, vertex, adjacent):
        if vertex in self.adjList:
            if adjacent in self.adjList[vertex]:
                del self.adjList[vertex][adjacent]
                
    # Verify if a given vertex exists
    #
    # verify if the vertex is in the adjList
    #
    # @param integer vertex the vertex to be searched
    # @return
                
    def hasVertex(self, vertex):
        return vertex in self.adjList
    
    # Verify if an edge exists between two vertex
    #
    # verify if the vertex is in the adjList and if it's adjacent is in the 
    # vertex adjList.
    #
    # @param vertex   integer the starting vertex
    # @param adjacent integer the adjacent vertex from the starting vertex
    # @return         boolean True or False

    def hasEdge(self, vertex, adjacent):
        return vertex in self.adjList and adjacent in self.adjList[vertex]
    
    # Returns the output degree of a vertex
    #
    # returns the lenght of the adjacent list from the given vertex
    #
    # @param vertex integer the vertex wich the degree will be calculated
    # @return integer if the vertex exists and None if it dont
    
    def outputDegree(self, vertex):
        if vertex in self.adjList:
            return len(self.adjList[vertex])
        return None
    
    # Returns the input degree of a given vertex
    #
    # iterate trough all vertexes in the adjList and verify if any of 
    # them is adjacent to the given vertex.
    #
    # @param vertex integer the vertex wich the degree will be calculated
    # @return integer if the vertex exists and None if it dont
    
    def inputDegree(self, vertex):
        if vertex in self.adjList:
            degree = 0
            for adjacent in self.adjList:
                if vertex in self.adjList[adjacent]:
                    degree += 1
            return degree
        return None
    
    def lengthSearch(self, v):
        visitado = []
        lista = [v]
        while lista:
            corrente = lista.pop()
            visitado.append(corrente)
            for vizinho in self.adjList[corrente]:
                if vizinho not in visitado and vizinho not in lista:
                    lista.append(vizinho)
        return visitado
    
    def depthSearch(self, v):
        visitado = []
        lista = [v]
        while lista:
            corrente = lista.pop()
            if corrente not in visitado:
                visitado.append(corrente)
                for vizinho in self.adjList[corrente]:
                    lista.append(vizinho)
        return visitado

     def areConected(self, v, w):
         if w in self.depthSearch(v): return True
         return False

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
                        return 'negative cicle between ' + (str) (vertex) + ' and ' + (str) (neighbour) 
        return dist

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

