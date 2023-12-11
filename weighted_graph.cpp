//
//  Copyright (C) 2018 Satoshi Shimizu
//
//
//  This file is part of MECQ shown in the following paper:
//  - Satoshi Shimizu, Kazuami Yamaguchi, Sumio Masuda,
//    ``A Maximum Edge-Weight Clique Extraction Algorithm Based on Branch-and-Bound,''
//    https://arxiv.org/abs/1810.10258, 2018.
//
//  MECQ is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  MECQ is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with MECQ.  If not, see <http://www.gnu.org/licenses/>.
//
#include<cstdlib>
#include<cstdio>
#include"weighted_graph.h"

weighted_graph::weighted_graph(char *inFile)
{
    using namespace std;
    FILE *fp=fopen(inFile,"r");
    n=0;
    degree=NULL;
    edge_weight=NULL;
    nonadjacency_list=NULL;
    adjacency_list=NULL;
    edge_weight_list=NULL;
    vertex_weight=NULL;
    if(fp == NULL)
    {
        fprintf(stderr,"\"%s\" doesn't exist\n",inFile);
        exit(1);
    }
    else
    {
        char form;
        while(fscanf(fp," %c",&form) != EOF)
        {
            switch(form)
            {
                case 'c':
                    if( fscanf(fp,"%*[^\n]") != 0)
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    break;
                case 'd':
                    if( fscanf(fp,"%*[^\n]") != 0)
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    break;
                case 'v':
                    if( fscanf(fp,"%*[^\n]") != 0)
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    break;
                case 'x':
                    if( fscanf(fp,"%*[^\n]") != 0)
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    break;
                case 'p':
                    if( fscanf(fp,"%*s %d %d",&n,&m) != 2 )
                    {
                        fprintf(stderr,"input file error\n");
                        exit(1);
                    }
                    vertex_weight=new int[n+1];
                    degree=new int[n+1];
                    for(int i=0;i<n;i++)
                    {
                        vertex_weight[i]=0;
                        degree[i]=0;
                    }
                    vertex_weight[n]=0; //dummy vertex
                    degree[n]=0; //dummy vertex
                    edge_weight=new int*[n];
                    for(int i=0;i<n;i++)
                    {
                        edge_weight[i]=new int[n];
                        for(int j=0;j<n;j++)
                        {
                            edge_weight[i][j]=NOT_ADJACENT;
                        }
                    }
                    break;
                case 'e':
                    {
                        int v1,v2,w,ret;
                        ret = fscanf(fp, "%d %d %d",&v1,&v2,&w);
                        if( ret == 2 )
                        {
                            w=0;
                        }
                        else if( ret != 3 )
                        {
                            fprintf(stderr,"input file error\n");
                            exit(1);
                        }
                        edge_weight[v1-1][v2-1]=w;
                        edge_weight[v2-1][v1-1]=w;
                        degree[v1-1]++;
                        degree[v2-1]++;
                    }
                    break;
                case 'n':
                    {
                        int v,w;
                        if( fscanf(fp,"%d %d",&v,&w) != 2 )
                        {
                            fprintf(stderr,"input file error\n");
                            exit(1);
                        }
                        vertex_weight[v-1]=w;
                    }
                    break;
                default:
                    fprintf(stderr,"input file error\n");
                    exit(1);
                    break;
            }
        }
    }
    fclose(fp);

    adjacency_list=new int*[n];
    edge_weight_list=new int*[n];
    for(int i=0; i<n; i++)
    {
        int *ewi=edge_weight[i];
        int di=degree[i];
        adjacency_list[i]=new int[di];
        edge_weight_list[i]=new int[di];
        int *adji=adjacency_list[i];
        int *ewli=edge_weight_list[i];
        int k=0;
        for(int j=0; j<n; j++)
        {
            if(ewi[j]!=NOT_ADJACENT)
            {
                adji[k]=j;
                ewli[k]=ewi[j];
                k++;
            }
        }
    }
}

weighted_graph::~weighted_graph()
{
    delete[] degree;
    delete[] vertex_weight;
    for(int i=0; i<n; i++)
    {
        delete[] edge_weight[i];
    }
    delete[] edge_weight;
    if(adjacency_list != NULL)
    {
        for(int i=0; i<n; i++)
        {
            delete[] adjacency_list[i];
        }
        delete[] adjacency_list;
    }
    if(edge_weight_list != NULL)
    {
        for(int i=0; i<n; i++)
        {
            delete[] edge_weight_list[i];
        }
        delete[] edge_weight_list;
    }
    if(nonadjacency_list != NULL)
    {
        for(int i=0; i<n; i++)
        {
            delete[] nonadjacency_list[i];
        }
        delete[] nonadjacency_list;
    }
}

bool weighted_graph::is_clique(int* clique, int clique_size,int clique_weight)
{
    int weight_sum=0;
    for(int i=0; i<clique_size; i++)
    {
        int ci=clique[i];
        weight_sum += vertex_weight[ci];
        for(int j=i+1; j<clique_size; j++)
        {
            if(edge_weight[ci][clique[j]] == NOT_ADJACENT)
            {
                return false;
            }
            weight_sum += edge_weight[ci][clique[j]];
        }
    }
    if(weight_sum != clique_weight)
    {
        return false;
    }
    return true;
}

void weighted_graph::create_nonadjlist()
{
    nonadjacency_list=new int*[n];
    for(int i=0; i<n; i++)
    {
        int *ewi=edge_weight[i];
        int di=degree[i];
        nonadjacency_list[i]=new int[n-di];
        int *nadji=nonadjacency_list[i];
        int k=0;
        for(int j=0; j<n; j++)
        {
            if(ewi[j]==NOT_ADJACENT)
            {
                nadji[k++]=j;
            }
        }
    }
}

void weighted_graph::delete_nonadjlist()
{
    for(int i=0; i<n; i++)
    {
        delete[] nonadjacency_list[i];
    }
    delete[] nonadjacency_list;
    nonadjacency_list = NULL;
}
