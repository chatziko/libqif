#include "Graph.h"
#include <iostream>
/*
This file belongs to the LIBQIF library.
A Quantitative Information Flow C++ Toolkit Library.
Copyright (C) 2013  Universidad Nacional de Río Cuarto(National University of Río Cuarto).
Author: Martinelli Fernán - fmartinelli89@gmail.com - Universidad Nacional de Río Cuarto (Argentina)
LIBQIF Version: 1.0
Date: 12th Nov 2013 
========================================================================
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

=========================================================================
*/
Graph::Graph(IntType vertex_num,StringType& edges)
{
    //traslate the string to a vector of integer 1pairs
    std::vector < std::pair<int,int> > new_edges;
    StringType copy=edges;
    StringType current;
    IntType pos,pairs=0,i=0;
    while(i<copy.size()){
        pos=copy.find(";");
        current= copy.substr(i,pos);
        i=pos+1;
        //vertex_num++;
        while(current[0]==' '){
            current = current.substr(1);
        }
        pos=current.find(" ");
        int left = std::atoi(current.substr(0,pos).c_str());
        int rigth = std::atoi(current.substr(pos).c_str());
        std::pair<int,int> new_pair = std::pair<int,int>(left,rigth);
        new_edges.push_back(new_pair);
        pairs++;
        if(i==0){
            copy="";
        }else{
            copy = copy.substr(i);
        }
        i=0;
    }

    Graph(vertex_num,new_edges);
}

Graph::Graph(IntType vertex_num, std::vector< std::pair<int, int> > & edges)
{
    //std::cerr << "HOME is not defined." << std::endl;
    if(vertex_num<0){
        throw 1;
    }
    V = vertex_num;
    adjacency = arma::mat(V,V);
    adjacency.zeros();
    //std::cerr << "HOME is not defined." << std::endl;
    /*
    for(int i=0; i<V; i++) {
        for(int j=0; j<V; j++) {
            adjacency(i,j) = 0;
        }
    }*/
    //std::cerr << "HOME is not defined." << std::endl;
    for(int i=0; i<edges.size(); i++) {
        int v1 = edges[i].first;
        int v2 = edges[i].second;
        
        adjacency(v1-1,v2-1) = 1;
        adjacency(v2-1,v1-1) = 1;
    }
    
    distances = arma::mat(V,V);
    // Just use Floyd-Warshall
    
    // Initialization (we use V+1 as infinity)
    for(int i=0; i<V; i++) {
        for(int j=0; j<V; j++) {
            if(i == j) distances(i,j) = 0;
            else if(adjacency(i,j) == 1) distances(i,j) = 1;
            else distances(i,j) = V+1;
        }
    }
    
    // Main algorithm
    for(int k=0; k<V; k++) {
        for(int i=0; i<V; i++) {
            for(int j=0; j<V; j++) {
                if(distances(i,j) > distances(i,k) + distances(k,j)) {
                    distances(i,j) = distances(i,k) + distances(k,j);
                }
            }
        }
    }
    
    // We set the infinity values to -1
    for(int i=0; i<V; i++) {
        for(int j=0; j<V; j++) {
            if(distances(i,j) > V) distances(i,j) = -1;
        }
    }
    
}

//Graph::~Graph()
//{
//    adjacency.~mat();
//    distances.~mat();
//}

IntType Graph::vertex_number()
{
    return V;
}

int Graph::get_distance(IntType v1, IntType v2)
{
    return distances(v1-1,v2-1);
}

bool Graph::is_an_edge(IntType v1, IntType v2)
{
    return adjacency(v1-1,v2-1)==1;
}